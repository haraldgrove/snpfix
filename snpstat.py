#!/usr/bin/env python

# snpstat, version 1.0, 2014-03-14

from __future__ import division, print_function
import sys
import argparse
import numpy as np

class SNP(object):

    def __init__(self,ic,ia):
        """ ic = # columns before genotypes
            ia = allele representation, 1 = 0/1/2, 2 = 1/2/3/4
        """
        self.ped = {}
        self.pedlist = []
        self.mark = {}
        self.marklist = []
        self.ic = ic
        self.ia = ia

    def tbasea(self,a,mark):
        m1 = self.mark[mark]['alleles'][0]
        m2 = self.mark[mark]['alleles'][1]
        if a == 1: return m1+m2
        if a == 0: return m1+m1
        if a == 2: return m2+m2
        return '00'

    def tbase012(self,a,mark):
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
            else:
                m2 = 'X'
        if a[0] != a[1]: return '1'
        if a[0] == m1: return '0'
        if a[0] == m2: return '2'
        return np.nan

    def readGenos(self,genofile):
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
        """ First 3 columns: name,father,mother """
        with open(pedfile,'r') as fin:
            count = 0
            for line in fin:
                if line.startswith('#'): continue
                l = line.strip().split()
                name,father,mother = '0','0','0'
                if len(l) > 0: name = l[0]
                if (real or self.ic > 1) and len(l) > 1: father = l[1]
                if (real or self.ic > 2) and len(l) > 2: mother = l[2]
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
            Columns options:
              name,position,allele1,allele2,chromosome
              chromosome,rank,name,position
              name
        """
        with open(markerfile,'r') as fin:
            count = 0
            for line in fin:
                if line.startswith('#'): continue
                l = line.strip().split()
                if len(l) == 0: continue
                if len(l) == 6: chrom,name,distance,position,a1,a2 = l
                elif len(l) == 4: chrom,name,distance,position,a1,a2 = l,[],[] # Plink
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
        # Father
        try:
            res = self.gen[self.ped[sire]['rank'],:]*self.gen[self.ped[anim]['rank'],:]
            wrongP = res==-1
        except KeyError:
            wrongP = []
        # Mother
        try:
            res = self.gen[self.ped[dam]['rank'],:]*self.gen[self.ped[anim]['rank'],:]
            wrongM = res==-1
        except KeyError:
            wrongM = []
        # Trios
        try:
            res = self.gen[self.ped[sire]['rank'],:]*self.gen[self.ped[dam]['rank'],:] - (self.gen[self.ped[anim]['rank'],:]*self.gen[self.ped[anim]['rank'],:])
            wrongT = res==1
        except KeyError:
            wrongT = []
        return np.logical_or(np.logical_or(wrongP,wrongM),wrongT),wrongP,wrongM,wrongT

    def correctMendel(self):
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
            if a == 0: return m[0]+'\t'+m[0]
            if a == 1: return m[0]+'\t'+m[1]
            if a == 2: return m[1]+'\t'+m[1]
            return '0\t0'
        
        with open(infile,'r') as fin:
            fout = open(outfile,'w')
            mlist = [''.join(self.mark[m]['alleles']) for m in self.marklist]
            for line in fin:
                if line.startswith('#'):
                    fout.write(line)
                    continue
                l = line.strip().split()
                fout.write('%s' % '\t'.join(l[0:self.ic]))
                child = l[0]
                if self.ia == 1:
                    fout.write('\t%s\n' % '\t'.join([trans1(g) for g in self.gen[self.ped[child]['rank'],:]]))
                elif self.ia == 2:
                    fout.write('\t%s\n' % '\t'.join([trans2(g,mark) for g,mark in zip(self.gen[self.ped[child]['rank'],:],mlist)]))
                elif self.ia == 3:
                    fout.write('\t%s\n' % '\t'.join([trans3(g,mark) for g,mark in zip(self.gen[self.ped[child]['rank'],:],mlist)]))
            fout.close()

def main():
    parser = argparse.ArgumentParser(description='Processes genotypes.')
    parser.add_argument('ingeno',help='Input genotypes file')
    parser.add_argument('outgeno',help='Output genotypes file')
    parser.add_argument('repfile',help='Output report file')
    parser.add_argument('-p','--pedigree',dest='pedigreefile',help='Pedigree file')
    parser.add_argument('-m','--markers',dest='markerfile',help='Marker file')
    parser.add_argument('-c',dest='infocol',type=int,help='Non-genotype columns', default=1)
    parser.add_argument('-a',dest='allele',type=int,help='Alleleformat, 1=0/1/2, 2=11/13/33, 3=1 1/1 3/3 3(plink)', default=1)
    parser.add_argument('-v','--verbose',help='Prints runtime info')
    args = parser.parse_args()
    gen = SNP(args.infocol,args.allele)
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
