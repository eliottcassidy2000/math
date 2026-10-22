from pathlib import Path
import ast,hashlib,numpy as np,sympy as s
p=str(Path(__file__).with_name('downloaded_G96_reduced.txt'));raw=open(p,'rb').read();G=[[int(x) for x in row.split(',')] for row in raw.decode().strip().strip('[]').split(';')];print('sha256',hashlib.sha256(raw).hexdigest(),flush=True);A=np.array(G);print('shape',A.shape,'symmetric',np.array_equal(A,A.T),'even diag',all(A.diagonal()%2==0),flush=True);print('eig extrema',np.linalg.eigvalsh(A)[[0,-1]],flush=True);print('det exact',s.polys.matrices.DomainMatrix.from_Matrix(s.Matrix(G)).det(),flush=True)
from fractions import Fraction as F
B=[[F(x) for x in r] for r in G]; piv=[]
for k in range(len(B)):
 d=B[k][k];piv.append(d)
 if d<=0:print('nonpositive pivot',k,d,flush=True);break
 for i in range(k+1,len(B)):
  for j in range(i,len(B)):
   B[j][i]=B[i][j]=B[i][j]-B[k][i]*B[k][j]/d
print('pivots',len(piv),'first',piv[:3],'last',piv[-3:],flush=True)
