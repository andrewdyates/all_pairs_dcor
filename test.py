#!/usr/bin/python
"""
In [9]: %timeit COR = compute_all_ppc(M)
1000 loops, best of 3: 811 us per loop

In [10]: %timeit SCI = squareform((pdist(M, metric='correlation')-1)*-1)
10000 loops, best of 3: 148 us per loop

In [12]: %timeit DCOR =  compute_all_dcor(M)
100 loops, best of 3: 5.01 ms per loop

In [19]: %timeit DCORL = test.loop_dcor(M)
10 loops, best of 3: 56.4 ms per loop
"""
from __init__ import *
import matrix_io as mio
import script
import numpy as np
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
FNAME="nice.may3.Eg.expr.gold.celegans.csv"
import dcor
M = mio.load(FNAME)['M']

def loop_dcor(M):
  m = M.shape[0]
  D = np.zeros((m,m))
  for i,rowi in enumerate(M):
    for j,rowj in enumerate(M):
      D[i,j] = dcor.dcor(rowi,rowj)
  return D

def main():
  DCOR = script.main(fname=FNAME)
  print DCOR
  COR = compute_all_ppc_numpy(M)
  print COR
  SCI = compute_all_pcc_scipy(M)
  print np.all(np.abs(COR-SCI) < 0.0000000000001)
  DCOR = compute_all_dcor(M)
  DCORL = loop_dcor(M)
  print np.all(np.abs(DCOR-DCORL) < 0.0000000000001)


  
if __name__ == "__main__":
  main()

