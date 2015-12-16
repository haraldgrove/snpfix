#!/usr/bin/env python

# snpalleles
# Reports the distribution of genotypes for each marker

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

    def readGenos(self,genofile,outfile):
        """
            Reads the genotype file and converts the genotypes into a numpy array as 0/1/2
            The -a parameter will specify if the alleles are given as:
              0/1/2 (-a 1),
              11/13/33 (-a 2),
              1 1/1 3/3 3 (-a 3)
        """
        self.gen = np.zeros((len(self.ped),len(self.mark)))
        self.gen[:] = np.nan
        maf = {}
        with open(genofile,'r') as fin:
            for line in fin:
                if line.startswith('#'):
                    mlist = line.strip('#').strip().split()
                    continue
                l = line.strip().split()
                if len(l) < 1: continue
                for i,mark  in enumerate(mlist):
                    if mark not in self.mark: continue
                    if self.ia == 2: 
                        a = l[i+self.ic]
                    elif self.ia == 3:
                        a = l[i*2+self.ic]+l[i*2+1+self.ic]
                    if len(a) == 1: a = a+a
                    if a == '00': continue
                    if a[0] not in self.mark[mark]['alleles']: self.mark[mark]['alleles'][a[0]] = 0
                    if a[1] not in self.mark[mark]['alleles']: self.mark[mark]['alleles'][a[1]] = 0
                    self.mark[mark]['alleles'][a[0]] += 1
                    self.mark[mark]['alleles'][a[1]] += 1
        if not outfile: fout = sys.stdout
        else: fout = open(outfile,'w')
        for mark in self.marklist:
            fout.write('%s' % mark)
            for a in self.mark[mark]['alleles']:
                fout.write('\t%s\t%d' % (a,self.mark[mark]['alleles'][a]))
            fout.write('\n')
                    
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
                                       'alleles': {},
                                       'rank':count}
                    if len(a1) > 0 and len(a2) > 0:
                        self.mark[name]['alleles'][a1] = 0
                        self.mark[name]['alleles'][a2] = 0
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
                                       'alleles': {},
                                       'rank':i}
                        self.marklist.append(e)
                    break
            else:
                l = line.strip().split()
                if self.ia == 3:
                    for i in xrange(0,len(l[self.ic:])//2):
                        self.mark[str(i)] = {'chrom':'0',
                                           'pos':i,
                                           'alleles': {},
                                           'rank':i}
                        self.marklist.append(str(i))
                else:
                    for i,e in enumerate(l[self.ic:]):
                        self.mark[str(i)] = {'chrom':'0',
                                           'pos':i,
                                           'alleles': {},
                                           'rank':i}
                        self.marklist.append(str(i))

#*****************************************************************************************************

def main():
    parser = argparse.ArgumentParser(description='Processes genotypes.')
    parser.add_argument('ingeno',help='Input genotypes file')
    parser.add_argument('-r','--repfile',help='Output report file')
    parser.add_argument('-m','--markers',dest='markerfile',help='Marker file')
    parser.add_argument('-n','--informat',dest='informat',help='Format of input file (Plink/DMU)')
    parser.add_argument('-c',dest='infocol',type=int,help='Non-genotype columns', default=3)
    parser.add_argument('-a',dest='allele',type=int,help='Alleleformat, 1=0/1/2, 2=11/13/33, 3=1 1/1 3/3 3', default=3)
    parser.add_argument('-v','--verbose',action="store_true",help='Prints runtime info')
    args = parser.parse_args()
    gen = SNP(args.infocol,args.allele,args.informat)
    # Collect marker information
    if args.markerfile:
        gen.readMarkers(args.markerfile)
    else:
        gen.collectMarkers(args.ingeno)
    gen.readGenos(args.ingeno,args.repfile)
    

if __name__ == '__main__':
    main()
