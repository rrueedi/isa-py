#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 15:48:17 2018

@author: rico
"""

import isa
import pandas, numpy
from optparse import OptionParser

def main():    
   
    def fcb(option, opt, value, parser):
        setattr(parser.values, option.dest, value.split(','))
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    
    parser.add_option('-i','--inputfile',dest='inpfile',type='string')
    parser.add_option('-o','--outputfile',dest='outfile',type='string',default='itsi.csv')
    parser.add_option('--dsame',dest='dsame',type='float',default=0.80)
    parser.add_option('--dconv',dest='dconv',type='float',default=0.975)
    parser.add_option('--nseed',dest='nseed',type='int',default=100)
    parser.add_option('--seedsparsity',dest='seedsparsity',type='int',default=0)
    parser.add_option('--maxiter',dest='maxiter',type='int',default=50)
    parser.add_option('--sgc',dest='sgc',type='int',default=0)
    parser.add_option('--sgr',dest='sgr',type='int',default=1)
    parser.add_option('--thc',dest='thc',type='string',action='callback',callback=fcb,default=[1,2,3])
    parser.add_option('--thr',dest='thr',type='string',action='callback',callback=fcb,default=[1,2,3])
    parser.add_option('--norm',dest='norm',type='string',default='double')

    (options, args) = parser.parse_args()
    
    A = pandas.read_csv(options.inpfile,index_col=0,header=0)
    A = A.fillna(0)
    a = A.values
    sthr = [float(x) for x in options.thr]
    sthc = [float(x) for x in options.thc]
    rsSR, csSC, sROB, sTHR, sTHC = \
    isa.itersigal(a,\
              sgr=numpy.sign(options.sgr),\
              sgc=numpy.sign(options.sgc),\
              seedsparsity=options.seedsparsity,\
              nseed=options.nseed,\
              normalisation_method=options.norm,\
              dconverged=options.dconv,\
              dsame=options.dsame,\
              sthr=sthr,\
              sthc=sthc,\
              maxiter=options.maxiter)
    idx=['row_threshold','col_threshold','robustness']
    idx.extend(A.index)
    idx.extend(A.columns)
    nf = '{:0'+str(int(numpy.ceil(numpy.log10(0.5+len(sROB)))))+'d}'
    
    col = ['M'+nf.format(x) for x in range(len(sROB))]
    B=pandas.DataFrame(numpy.vstack([sTHR,sTHC,sROB,rsSR,csSC]),index=idx,columns=col)
    B.to_csv(options.outfile)
    
if __name__ == '__main__':
    main()
    