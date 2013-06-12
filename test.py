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
import dcor

FNAME="nice.may3.Eg.expr.gold.celegans.csv"
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
  COR = compute_all_ppc_numpy(M)
  SCI = compute_all_pcc_scipy(M)
  print np.all(np.abs(COR-SCI) < 0.0000000000001)
  DCOR = compute_all_dcor(M)
  DCORL = loop_dcor(M)
  print np.all(np.abs(DCOR-DCORL) < 0.0000000000001)
  DCOR2 = compute_all_dcor_2(M)
  print np.all(np.abs(DCOR-DCOR2) < 0.0000000000001)
  DCOR3 = compute_all_dcor_3(M)
  print np.all(np.abs(DCOR-DCOR3) < 0.0000000000001)
  print np.all(np.abs(DCOR-DCOR3) < 0.1)
  
  print DCOR[5,3], DCOR[3,5]
  print DCOR3[5,3], DCOR3[3,5]

if __name__ == "__main__":
  main()

