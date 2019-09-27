#!/usr/bin/env python

import sys
import argparse
import numpy as np
import gzip
import random
import time

PUNNETT = {'00': np.array([1.0,0.0,0.0]),
               '01': np.array([0.5,0.5,0.0]),
               '02': np.array([0.0,1.0,0.0]),
               '11': np.array([0.25,0.5,0.25]),
               '12': np.array([0.0,0.5,0.5]),
               '22': np.array([0.0,0.0,1.0])}

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
            if a[0] != a[1]:
                self.mark[mark]['alleles'] = [a[0],a[1]]
                m1 = a[0]
            elif a[0] != '0':
                self.mark[mark]['alleles'].append(a[0])
                m1 = a[0]
            else:
                m1 = 'X'
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
        if a[0] != '0':
            sys.stderr.write('ERROR: Marker %s has more than 2 alleles\n' % mark)
            sys.stderr.write('%s\t%s\t%s%s\n' % (a,mark,m1,m2))
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
        marklist = None
        with open(genofile,'r') as fin:
            for line in fin:
                if line.startswith('#'):
                    if not marklist: marklist = line.strip('#').strip().split()
                    continue
                l = line.strip().split()
                if len(l) < 1: continue
                try: irow = self.ped[l[self.nc]]['rank']
                except KeyError:
                    continue
                for i,mark  in enumerate(self.marklist):
                    if mark not in self.mark: continue
                    icol = self.mark[mark]['rank']
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

#**************************************************************************************

    def check3(self,obs):
        est = np.array([0.25,0.5,0.25])*sum(obs)
        res = sum(np.power(obs-est,2)/est)
        return res

    def check2(self,obs):
        est = np.array([0.5,0.5])*sum(obs)
        res = sum(np.power(obs-est,2)/est)
        return res

    def check1(self,obs):
        res = -1*np.log(np.power(0.5,obs))
        return res

    def calcStat(self,obs):
        res = {}
        indO = obs>0
        if sum(indO) == 3:
            val = self.check3(obs)
            res['11'] = val
            if val < 9.210: return res # 1:2:1 distribution at alpha == 0.01, v==2
            return {'nan':9999}
            if obs[0] < obs[2]:
                val = self.check2(obs[1:3])
                res['12'] = val
            else:
                val = self.check2(obs[0:2])
                res['01'] = val
            if val < 6.635: return res
            maxind = np.argmax(obs)
            val = self.check1(obs[maxind])
            if maxind == 0: res['00'] = val
            elif maxind == 1: res['02'] = val
            else: res['22'] = val
            return res
        if sum(indO) == 2:
            val = self.check2(obs[indO])
            if obs[0] < obs[2]: res['12'] = val
            else: res['01'] = val
            if val < 6.635: return res # 1:1 distribution at alpha == 0.01, v==1
            return {'nan':9999}
            val = self.check3(obs)
            res['11'] = val
            if val < 9.210: return res # 1:2:1 distribution at alpha == 0.01, v==2
            maxind = np.argmax(obs)
            val = self.check1(obs[maxind])
            if maxind == 0: res['00'] = val
            elif maxind == 1: res['02'] = val
            else: res['22'] = val
            return res
        if sum(indO) == 1:
            val = self.check1(obs[indO])
            if obs[0] > 0: res['00'] = val
            elif obs[1] > 0: res['02'] = val
            else: res['22'] = val
            return res # 1 distribution at probability p = 0.01
        return res

    def getfathers(self):
        self.fathers = {}
        self.mothers = {}
        for animal in self.pedlist:
            father,mother = self.ped[animal]['father'],self.ped[animal]['mother']
            if father == '0' or mother == '0':
                continue
            if father not in self.fathers:
                self.fathers[father] = {}
                self.fathers[father][mother] = [self.ped[animal]['rank']]
            elif mother not in self.fathers[father]:
                self.fathers[father][mother] = [self.ped[animal]['rank']]
            else:
                self.fathers[father][mother].append(self.ped[animal]['rank'])
            if mother not in self.mothers:
                self.mothers[mother] = {}
                self.mothers[mother][father] = [self.ped[animal]['rank']]
            elif father not in self.mothers[mother]:
                self.mothers[mother][father] = [self.ped[animal]['rank']]
            else:
                self.mothers[mother][father].append(self.ped[animal]['rank'])

    def inferAlleles(self,reportfile,logfile):
        #self.simulate(10,10)
        self.getfathers()
        fout = open(reportfile,'w')
        flog = open(logfile,'w')
        missing = {}
        flog.write('marker\tfather\tmother\tEstGeno\tstatistic\tObsGeno\ta0\ta1\ta2\n')
        for marker in self.marklist:
            fathrec = {}
            mothrec = {}
            genos = {}
            icol = self.mark[marker]['rank']
            for father in self.fathers:
                fathrec[father] = {'0':0,'1':0,'2':0,'rep':0}
                af = np.nan
                for mother in self.fathers[father]:
                    if mother not in mothrec: 
                        mothrec[mother] = {'0':0,'1':0,'2':0,'rep':0}
                    ind = self.fathers[father][mother]
                    gentyp = self.gen[ind,icol]
                    a0 = len(gentyp[gentyp==0])*1.0
                    a1 = len(gentyp[gentyp==1])*1.0
                    a2 = len(gentyp[gentyp==2])*1.0
                    stat = self.calcStat(np.array([a0,a1,a2]))
                    am = np.nan
                    # Prints all checked possibilities key=genotype, 
                    # stat[key] is the statistic for that genotype
                    for key in stat:
                        if stat[key] > 999:
                            continue
                        flog.write('%s\t%s\t%s\t%s\t%s\t%s%s\t%s\t%s\t%s\n' \
                                % (marker,father,mother,key,stat[key],af,am,a0,a1,a2))
                        genos[father,mother] = key
                        fathrec[father]['rep'] += 1
                        mothrec[mother]['rep'] += 1
                        if '0' in key:
                            fathrec[father]['0'] += 1
                            mothrec[mother]['0'] += 1
                        if '1' in key:
                            fathrec[father]['1'] += 1
                            mothrec[mother]['1'] += 1
                        if '2' in key:
                            fathrec[father]['2'] += 1
                            mothrec[mother]['2'] += 1
            for father in fathrec:
                a0 = fathrec[father]['0']
                a1 = fathrec[father]['1']
                a2 = fathrec[father]['2']
                r = fathrec[father]['rep']
                if r < 2: continue # only do assignment if the parent is present in 2 or more matings
                if a0 == r and a0>a1 and a0>a2:
                    a = '0'
                elif a1 == r and a1>a0 and a1>a2:
                    a = '1'
                elif a2 == r and a2>a0 and a2>a1:
                    a = '2'
                else:
                    a = 'nan'
                    continue
                fout.write('%s\t%s\t%s\t%s\n' % (marker,father,a,'M'))
                for mother in self.fathers[father]:
                    try:
                        g = genos[father,mother]
                    except KeyError:
                        continue
                    b = g.replace(a,'',1)
                    fout.write('%s\t%s\t%s\t%s\n' % (marker,mother,b,'F'))
            for mother in mothrec:
                a0 = mothrec[mother]['0']
                a1 = mothrec[mother]['1']
                a2 = mothrec[mother]['2']
                r = mothrec[mother]['rep']
                if r < 2: continue
                if a0 == r and a0>a1 and a0>a2:
                    a = '0'
                elif a1 == r and a1>a0 and a1>a2:
                    a = '1'
                elif a2 == r and a2>a0 and a2>a1:
                    a = '2'
                else:
                    a = 'nan'
                    continue
                fout.write('%s\t%s\t%s\t%s\n' % (marker,mother,a,'F'))
                for father in self.mothers[mother]:
                    try: g = genos[father,mother]
                    except KeyError: continue
                    b = g.replace(a,'',1)
                    fout.write('%s\t%s\t%s\t%s\n' % (marker,father,b,'M'))
        fout.close()
        flog.close()

def simulate(N,n):
    count1 = 0
    count2 = 0
    for i in xrange(N):
        res1 = [0,0,0]
        res2 = [0,0]
        for j in xrange(n):
            r = random.random()
            if r < 0.25: res1[0] += 1
            elif r < 0.75: res1[1] += 1
            else: res1[2] += 1
            if r < 0.5: res2[0] += 1
            else: res2[1] += 1
        if res1[0] < 5 or res1[2] < 5:
            count1 += 1
        if abs(res2[0]-res2[1]) > 14:
            count2 += 1 
    print('AaxAa:%.5f\nAAxAa:%.5f' % (count1/N,count2/N))
 
def main():
    parser = argparse.ArgumentParser(description='Processes genotypes.')
    parser.add_argument('ingeno',help='Input genotypes file')
    parser.add_argument('-o', '--outgeno',help='Output genotypes file', default='out.txt')
    parser.add_argument('-r','--repfile',help='Output report file', default='report.txt')
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
        print('No pedigreefile found, gathering pedigreedata from genotypefile.')
        gen.readPedigree(args.ingeno,False)
    # Collect marker information
    if args.markerfile:
        gen.readMarkers(args.markerfile)
    else:
        print('No markerfile found, gathering markerdata from genotypefile.')
        gen.collectMarkers(args.ingeno)
    gen.readGenos(args.ingeno)
    gen.inferAlleles(args.repfile,args.outgeno)
    

if __name__ == '__main__':
    t = time.time()
    main()
    sys.stdout.write('Time spent: %.3f\n' % (time.time()-t))
