#!/usr/bin/python
"""Compute all pairs euclidean distance correlation using matrix multiplications.
Uses n^2M memory.

Use:
  python all_pairs_dcor/script.py fname=nice.may3.Eg.expr.gold.celegans.csv
"""
import sys
import matrix_io as mio
from __init__ import *
import cPickle as pickle

def main(fname=None, pkl=True, algorithm="3", outtag="", **kwds):
  assert fname
  if isinstance(pkl, basestring) and pkl.lower() in ('f','false','none'): pkl = False
  print "Loading data from %s..." % fname
  D = mio.load(fname)
  print "Computing all pairs (euclidean) distance correlation from a (%d x %d) data matrix to a (%d x %d) result matrix..." % (D['M'].shape[0], D['M'].shape[1], D['M'].shape[0], D['M'].shape[0])
  print "Computing all pairs (euclidean) distance correlation..."
  
  if algorithm == "1":
    print "Using Algorithm 1: single dot product, n^2*m memory"
    DCOR = compute_all_dcor(D['M'], **kwds)
  elif algorithm == "2":
    print "Using Algorithm 2: multiple dot products, n*m memory"
    DCOR = compute_all_dcor_2(D['M'], **kwds)
  elif algorithm == "3":
    print "Using Algorithm 3: multiple dot products, n*m memory, n choose 2 savings"
    DCOR = compute_all_dcor_3(D['M'], **kwds)
  else:
    raise Exception, "Unknown algorithm %s" % algorithm

  if outtag and outtag[-1] != ".":
    outtag += "."
  fname_out = '%s.%sdcor.tab' % (fname, outtag)
  print "Saving %s..." % (fname_out)
  mio.save(DCOR, fname_out, fmt="%.4f", row_ids=D['row_ids'], col_ids=D['row_ids'])
  if pkl:
    fname_pkl_out = fname_out.rpartition('.')[0]+'.pkl'
    print "Saving %s..." % (fname_pkl_out)
    pickle.dump(DCOR, open(fname_pkl_out,"w"), protocol=-1)
  return DCOR

if __name__ == "__main__":
  args = dict([s.split('=') for s in sys.argv[1:]])
  print args
  main(**args)
