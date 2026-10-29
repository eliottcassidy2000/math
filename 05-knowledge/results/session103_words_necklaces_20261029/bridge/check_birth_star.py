import sys,itertools
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'singularity'))
from native_kernel import observer
for m in [3,4,5]:
 fail=[];count=0
 edges=list(itertools.combinations(range(1,m),2))
 for bits in range(1<<len(edges)):
  R=[[0]*m for _ in range(m)];R[1][0]=1
  for j in range(2,m):R[0][j]=1
  for k,(i,j) in enumerate(edges):R[i][j]=(bits>>k)&1;R[j][i]=1-R[i][j]
  for q in [0,1]:
   n=m+3;A=[row+[0]*3 for row in R]+[[1]*m+[0]*3 for _ in range(3)]
   for i,j in [(m,m+1),(m+1,m+2),(m+2,m)]:A[i][j]=1
   A[q][m]=1;A[m][q]=0
   ports,M=observer(A);idx={p:i for i,p in enumerate(ports)}
   ok=all(row[idx[0,m]]-row[idx[m,0]]-row[idx[1,m]]+row[idx[m,1]]==0 for row in M)
   if not ok:fail.append((bits,q))
   count+=1
 if fail:raise ArithmeticError(fail)
 print(m,count,len(fail),flush=True)
