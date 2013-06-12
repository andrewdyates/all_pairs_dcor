#!/usr/bin/python
import numpy as np

def compute_all_dcor_2(M):
  """Iteratively construct all-pairs distance matrices.
  Medium overhead, but much lower memory usage than compute_all_dcor.
  """
  if not isinstance(M, np.matrix):
    M = np.matrix(M)
  m,n = M.shape
  RowMeans = np.matrix(np.zeros((m,n)))
  for j in xrange(n):
    RowMeans[:,j] = np.mean(np.abs(M-M[:,j]),1)
  DMeans = np.mean(RowMeans,1)
  DCOV = np.zeros((m,m))
  for j in xrange(n):
    A = np.abs(M-M[:,j]) - RowMeans - RowMeans[:,j] + DMeans
    DCOV += np.dot(A,A.T)
  DCOV = np.sqrt(DCOV) / n
  DSTD = np.sqrt(np.diag(DCOV))
  DCOR = DCOV / DSTD / DSTD.reshape(m,1)
  return DCOR

def compute_all_dcor_3(M):
  """Iteratively construct all-pairs distance matrices with n choose 2 savings.
  Highest overhead, but algorithmically more efficient then compute_all_dcor2; same memory.
  """
  if not isinstance(M, np.matrix):
    M = np.matrix(M)
  m,n = M.shape
  RowMeans = np.matrix(np.zeros((m,n)))
  for j in xrange(n):
    RowMeans[:,j] = np.mean(np.abs(M-M[:,j]),1)
  DMeans = np.mean(RowMeans,1)
  DCOV = np.zeros((m,m))

  # upper triangle
  for j in xrange(n-1):
    A = np.abs(M[:,j+1:n]-M[:,j]) - RowMeans[:,j+1:n] - RowMeans[:,j] + DMeans
    DCOV += np.dot(A,A.T)
  DCOV *= 2
  # diagonal
  for j in xrange(n):
    A = DMeans - 2*RowMeans[:,j]
    DCOV += np.dot(A,A.T)
  DCOV = np.sqrt(DCOV) # / n**2 (term cancels out in dCor division)
  DSTD = np.sqrt(np.diag(DCOV))
  DCOR = DCOV / DSTD / DSTD.reshape(m,1)
  return DCOR
                      
    
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
  DCOV = np.sqrt(np.dot(DIST,np.transpose(DIST)))/n
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
  D = squareform((pdist(M, metric='correlation')-1)*-1)
  np.fill_diagonal(D, 1.0)
  return D
