#!/usr/bin/env python

# snpstat, version 1.0, 2014-03-14

from __future__ import division, print_function
import sys
import argparse
import numpy as np

class SNP(object):

    def __init__(self,ic,ia,informat):
        """
            ic = number of information columns before genotypes start
            ia = allele representation, 1 = 0/1/2, 2 = 11/13/33, 3 = 1 1/1 3/3 3
                 for 2 and 3, alleles are represented as 1/2/3/4
        """
        self.ped = {}
        self.pedlist = []
        self.mark = {}
        self.marklist = []
        self.sep = '\t'
        if informat in ['Plink','plink']:
            self.ic = 6
            self.ia = 3
            self.nc = 1
        elif informat in ['DMU','dmu']:
            self.ic = 1
            self.ia = 1
            self.nc = 0
        elif not informat:
            self.ic = ic
            self.ia = ia
            self.nc = 0
        else:
            sys.stderr.write('Unknown input format: "%s"\n' % informat)
            sys.exit(1)

    def tbasea(self,a,mark):
        # Transforms from 0/1/2 to 1/2/3/4
        m1 = self.mark[mark]['alleles'][0]
        m2 = self.mark[mark]['alleles'][1]
        if a == 1: return m1+m2
        if a == 0: return m1+m1
        if a == 2: return m2+m2
        return '00'

    def tbase012(self,a,mark):
        # Transform from 1/2/3/4 to 0/1/2
        try:
            m1 = self.mark[mark]['alleles'][0]
        except IndexError:
            if a[0] != a[1]: self.mark[mark]['alleles'] = [a[0],a[1]]
            elif a[0] != '0': self.mark[mark]['alleles'].append(a[0])
            m1 = a[0]
        try:
            m2 = self.mark[mark]['alleles'][1]
        except IndexError:
            if a[0] != '0' and a[0] != m1:
                self.mark[mark]['alleles'].append(a[0])
                m2 = a[0]
            elif a[1] != '0' and a[1] != m1:
                self.mark[mark]['alleles'].append(a[1])
                m2 = a[1]
            else:
                # The second allelel has not been encountered yet.
                m2 = 'X'
        if a[0] != a[1]: return '1'
        if a[0] == m1: return '0'
        if a[0] == m2: return '2'
        return np.nan

    def readGenos(self,genofile):
        """
            Reads the genotype file and converts the genotypes into a numpy array as 0/1/2
            The -a parameter will specify if the alleles are given as:
              0/1/2 (-a 1),
              11/13/33 (-a 2),
              1 1/1 3/3 3 (-a 3)
        """
        self.gen = np.zeros((len(self.ped),len(self.mark)))
        self.gen[:] = np.nan
        with open(genofile,'r') as fin:
            for line in fin:
                if line.startswith('#'): continue
                l = line.strip().split()
                if len(l) < 1: continue
                irow = self.ped[l[0]]['rank']
                for i,mark  in enumerate(self.marklist):
                    if mark not in self.mark: continue
                    icol = self.mark[self.marklist[i]]['rank']
                    if self.ia == 1:
                        a = l[i+self.ic]
                    elif self.ia == 2: 
                        a = self.tbase012(l[i+self.ic],mark)
                    elif self.ia == 3:
                        a = self.tbase012(l[i*2+self.ic]+l[i*2+1+self.ic],mark)
                    if a not in ['0','1','2']: a = np.nan
                    else: a = int(a)
                    self.gen[irow,icol] = a

    def readPedigree(self,pedfile,real=True):
        """
            Reads a pedigree from either a separate pedigree file or from the given genotype file (real=False)
            The first 3 columns have to be: sample_name, father and mother, regardless of input file
            If the number of information columns are either 1 or 2, it will assume that no parents are
            present or just the father, respectively.
        """
        with open(pedfile,'r') as fin:
            count = 0
            for line in fin:
                if line.startswith('#'): continue
                l = line.strip().split()
                name,father,mother = '0','0','0'
                if len(l) > 0: name = l[self.nc]
                if (real or self.ic > 1) and len(l) > 1: father = l[self.nc+1]
                if (real or self.ic > 2) and len(l) > 2: mother = l[self.nc+2]
                if name == '0': continue
                if name not in self.ped:
                    self.ped[name] = {'father':father,
                                      'mother':mother,
                                      'rank':count,
                                      'children':[]}
                    count += 1
                    self.pedlist.append(name)
                else:
                    sys.stderr.write('%s present more than once\n' % name)
        self.updatePed()

    def updatePed(self):
        # Assign children to parents and set sex of parents
        for n in self.ped:
            father,mother = self.ped[n]['father'],self.ped[n]['mother']
            if father in self.ped:
                #self.ped[father]['sex'] = '1'
                self.ped[father]['children'].append(n)
            if mother in self.ped:
                #self.ped[mother]['sex'] = '0'
                self.ped[mother]['children'].append(n)

    def readMarkers(self,markerfile):
        """
            Read a marker file, format is Plink map-file
            Will also accept the addition of the marker alleles as two extra columns
            and a marker file containing just one column of marker names.
        """
        with open(markerfile,'r') as fin:
            count = 0
            for line in fin:
                if line.startswith('#'): continue
                l = line.strip().split()
                if len(l) == 0: continue
                if len(l) == 6: chrom,name,distance,position,a1,a2 = l
                elif len(l) == 4:
                    chrom,name,distance,position = l # Plink
                    a1,a2 = [],[]
                elif len(l) == 1:
                    name = l[0]
                    chrom,pos,a1,a2 = '0',count,[],[]
                if name not in self.mark:
                    self.mark[name] = {'chrom':chrom,
                                       'pos':int(position),
                                       'alleles': a1+a2,
                                       'rank':count}
                    count += 1
                    self.marklist.append(name)

    def collectMarkers(self, ingeno):
        """ 
            Reads input file to search for the necessary marker information
            If the markers are not present as a comment line on top of the file,
            it will calculate the number of markers to be the same as the number of alleles
            after the information columns.
        """
        with open(ingeno,'r') as fin:
            for line in fin:
                if line.startswith('#'):
                    l = line.strip('#').strip().split()
                    for i,e in enumerate(l):
                        self.mark[e] = {'chrom':'0',
                                       'pos':i,
                                       'alleles': [],
                                       'rank':i}
                        self.marklist.append(e)
                    break
            else:
                l = line.strip().split()
                if self.ia == 3:
                    for i in xrange(0,len(l[self.ic:])//2):
                        self.mark[str(i)] = {'chrom':'0',
                                           'pos':i,
                                           'alleles': [],
                                           'rank':i}
                        self.marklist.append(str(i))
                else:
                    for i,e in enumerate(l[self.ic:]):
                        self.mark[str(i)] = {'chrom':'0',
                                           'pos':i,
                                           'alleles': [],
                                           'rank':i}
                        self.marklist.append(str(i))

#*****************************************************************************************************

    def findDiscords(self,anim,sire,dam):
        """ 
            Detects all mendelian discords within the given trio
            Will work even if one of the parents is missing
        """
        # Father
        try:
            res = self.gen[self.ped[sire]['rank'],:]*self.gen[self.ped[anim]['rank'],:]
            wrongP = res==-1
        except KeyError:
            wrongP = [False]*len(self.mark)
        # Mother
        try:
            res = self.gen[self.ped[dam]['rank'],:]*self.gen[self.ped[anim]['rank'],:]
            wrongM = res==-1
        except KeyError:
            wrongM = [False]*len(self.mark)
        # Trios
        try:
            res = self.gen[self.ped[sire]['rank'],:]*self.gen[self.ped[dam]['rank'],:] - (self.gen[self.ped[anim]['rank'],:]*self.gen[self.ped[anim]['rank'],:])
            wrongT = res==1
        except KeyError:
            wrongT = [False]*len(self.mark)
        return np.logical_or(np.logical_or(wrongP,wrongM),wrongT),wrongP,wrongM,wrongT

    def correctMendel(self):
        """ 
            Checks for mendelian discords and sets all conflicting alleles to missing
            It is possible the deletion could be handled more intelligently.
        """
        self.gen -= 1
        self.disc = np.zeros((len(self.ped),len(self.mark)))
        one = np.ones((1,len(self.mark)))
        for child in self.pedlist:
            father = self.ped[child]['father']
            mother = self.ped[child]['mother']
            index,indP,indM,indT = self.findDiscords(child,father,mother)
            self.disc[self.ped[child]['rank'],index] += 1
            if father in self.ped:
                self.disc[self.ped[father]['rank'],indP] += 1
            if mother in self.ped:
                self.disc[self.ped[mother]['rank'],indM] += 1
        self.gen[self.disc>0] = np.nan
        self.gen += 1

    def calcMAF(self,maffile):
        """ 
            Calculates statistics on allele distribution, MAF and error percentage for all markers
        """
        with open(maffile,'w') as fout:
            fout.write('#SNP\tMAF\ta0\ta1\ta2\tan\terr\tadjErr\n')
            for n in self.marklist:
                i = self.mark[n]['rank']
                d = self.gen[:,i]
                ind = np.isfinite(d)
                try:
                    MAF = sum(d[ind])/(2*(len(d[ind])))
                except ZeroDivisionError:
                    MAF = np.nan
                a0,a1,a2,an = len(d[d==0]),len(d[d==1]),len(d[d==2]),len(d[np.isnan(d)])
                if MAF>0.5: MAF = 1-MAF
                if len(d[ind]) > 0:
                    err = sum(self.disc[:,i])/(len(d[ind])+sum(self.disc[:,i]))
                else:
                    err = 0
                if MAF > 0:
                    relerr = err / (MAF+MAF+MAF*MAF)
                else:
                    relerr = err
                fout.write('%s\t%.5f\t%d\t%d\t%d\t%d\t%.3f\t%.3f\n' % (n,MAF,a0,a1,a2,an,err,relerr))

    def calcDist(self,maffile):
        """ 
            Calculates statistics on allele distribution and error percentage for the samples
        """
        with open(maffile,'a') as fout:
            fout.write('#Sample\tx\ta0\ta1\ta2\tan\terr\tadjErr\n')
            for n in self.pedlist:
                i = self.ped[n]['rank']
                d = self.gen[i,:]
                ind = np.isfinite(d)
                a0,a1,a2,an = len(d[d==0]),len(d[d==1]),len(d[d==2]),len(d[np.isnan(d)])
                if len(d[ind]) > 0:
                    err = sum(self.disc[i,:])/(len(d[ind])+sum(self.disc[i,:]))
                else:
                    err = 0
                MAF = 0
                relerr = err
                fout.write('%s\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\n' % (n,MAF,a0,a1,a2,an,err,relerr))

    def writeGeno(self,infile,outfile):
        """ 
            Writes corrected genotypes to the outputfile in the same format as the input files
            The information columns (specified by the -c parameter) are re-read from the input file
        """
            
        def trans1(a):
            if a == 0: return '0'
            if a == 1: return '1'
            if a == 2: return '2'
            return '-1'

        def trans2(a,m):
            if a == 0: return m[0]+m[0]
            if a == 1: return m[0]+m[1]
            if a == 2: return m[1]+m[1]
            return '00'

        def trans3(a,m):
            if a == 0: return m[0]+self.sep+m[0]
            if a == 1: return m[0]+self.sep+m[1]
            if a == 2: return m[1]+self.sep+m[1]
            return '0'+self.sep+'0'
        
        with open(infile,'r') as fin:
            fout = open(outfile,'w')
            mlist = [''.join(self.mark[m]['alleles']) for m in self.marklist]
            for e in zip(mlist,self.marklist):
                if len(e[0]) <2:
                    print(e)
            for line in fin:
                if line.startswith('#'):
                    fout.write(line)
                    continue
                l = line.strip().split()
                fout.write('%s' % self.sep.join(l[0:self.ic]))
                child = l[self.nc]
                if self.ia == 1:
                    fout.write('\t%s\n' % self.sep.join([trans1(g) for g in self.gen[self.ped[child]['rank'],:]]))
                elif self.ia == 2:
                    fout.write('\t%s\n' % self.sep.join([trans2(g,mark) for g,mark in zip(self.gen[self.ped[child]['rank'],:],mlist)]))
                elif self.ia == 3:
                    fout.write('\t%s\n' % self.sep.join([trans3(g,mark) for g,mark in zip(self.gen[self.ped[child]['rank'],:],mlist)]))
            fout.close()

def main():
    parser = argparse.ArgumentParser(description='Processes genotypes.')
    parser.add_argument('ingeno',help='Input genotypes file')
    parser.add_argument('outgeno',help='Output genotypes file')
    parser.add_argument('repfile',help='Output report file')
    parser.add_argument('-p','--pedigree',dest='pedigreefile',help='Pedigree file',default=None)
    parser.add_argument('-m','--markers',dest='markerfile',help='Marker file')
    parser.add_argument('-n','--informat',dest='informat',help='Format of input file (Plink/DMU)')
    parser.add_argument('-c',dest='infocol',type=int,help='Non-genotype columns', default=1)
    parser.add_argument('-a',dest='allele',type=int,help='Alleleformat, 1=0/1/2, 2=11/13/33, 3=1 1/1 3/3 3', default=1)
    parser.add_argument('-v','--verbose',action="store_true",help='Prints runtime info')
    args = parser.parse_args()
    gen = SNP(args.infocol,args.allele,args.informat)
    # Collect pedigree information
    if args.pedigreefile:
        gen.readPedigree(args.pedigreefile)
    else:
        gen.readPedigree(args.ingeno,False)
    # Collect marker information
    if args.markerfile:
        gen.readMarkers(args.markerfile)
    else:
        gen.collectMarkers(args.ingeno)
    gen.readGenos(args.ingeno)
    gen.correctMendel()
    gen.calcMAF(args.repfile)
    gen.calcDist(args.repfile)
    gen.writeGeno(args.ingeno, args.outgeno)
    

if __name__ == '__main__':
    main()
