import sys,itertools
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'singularity'))
from native_kernel import observer
m=3;R=[[0,0,1],[1,0,0],[0,1,0]]
for l in [3,4,5]:
 fail=[];count=0;free=list(itertools.combinations(range(m,m+l-1),2))
 for bits in range(1<<len(free)):
  n=m+l;A=[row+[0]*l for row in R]+[[1]*m+[0]*l for _ in range(l)]
  for i in range(m+1,n-1):A[i][n-1]=1
  A[n-1][m]=1
  for k,(i,j) in enumerate(free):A[i][j]=(bits>>k)&1;A[j][i]=1-A[i][j]
  for q in [0,1]:
   A[q][m]=1;A[m][q]=0
   ports,M=observer(A);idx={p:i for i,p in enumerate(ports)}
   ok=all(row[idx[0,m]]-row[idx[m,0]]-row[idx[1,m]]+row[idx[m,1]]==0 for row in M)
   if not ok:fail.append((bits,q))
   count+=1;A[q][m]=0;A[m][q]=1
 if fail:raise ArithmeticError(fail)
 print(l,count,len(fail),flush=True)
