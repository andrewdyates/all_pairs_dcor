#!/usr/bin/python
import numpy as np

def compute_all_dcor(M):
  """Return matrix of all-row-pairs distance correlations."""
  Mat = np.matrix(M)
  m,n = Mat.shape #m: num rows, n: num cols
  DIST = np.ones((m,n**2))
  for i, row in enumerate(Mat):
    D = np.abs(row-row.T)
    u_row = np.mean(D,1).reshape(n,1)
    u_col = np.mean(D,0)
    u = np.mean(D)
    D = D - u_row - u_col + u
    DIST[i,] = D.ravel()
  DCOV = np.sqrt(np.dot(DIST,np.transpose(DIST)))
  DSTD = np.sqrt(np.diag(DCOV))
  DCOV = DCOV / DSTD / DSTD.reshape(m,1)
  return DCOV
  
