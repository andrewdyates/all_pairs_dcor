#!/usr/bin/python
from __init__ import *
import matrix_io as mio
import script
import numpy as np
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
FNAME="nice.may3.Eg.expr.gold.celegans.csv"

def main():
  DCOR = script.main(fname=FNAME)
  print DCOR
  M = mio.load(FNAME)['M']
  COR = compute_all_ppc(M)
  print COR
  SCI = squareform((pdist(M, metric='correlation')-1)*-1)
  print SCI
  print COR-SCI


  
if __name__ == "__main__":
  main()

