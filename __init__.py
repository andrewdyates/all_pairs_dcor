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
  DCOR = DCOV / DSTD / DSTD.reshape(m,1)
  return DCOR

def compute_all_ppc_numpy(M):
  """Return matrix of all-row-pairs correlation using numpy functions."""
  m,n = M.shape
  u_row = np.mean(M,1).reshape(m,1)
  Mat = M-u_row
  COV = np.dot(Mat,np.transpose(Mat))
  STD = np.sqrt(np.diag(COV))
  COR = COV / STD / STD.reshape(m,1)
  return COR
  
def compute_all_pcc_scipy(M):
  """More efficient compute_all_pcc with scipy dependency."""
  from scipy.spatial.distance import pdist
  from scipy.spatial.distance import squareform
  return np.fill_diagonal(squareform((pdist(M, metric='correlation')-1)*-1),1.0)
