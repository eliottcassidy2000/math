from collections import Counter
from itertools import combinations
def paley(n):
 S={i*i%n for i in range(1,n)}
 return [[int(i!=j and (j-i)%n in S)for j in range(n)]for i in range(n)]
def double(A):
 n=len(A);B=[[0]*(2*n+1)for _ in range(2*n+1)]
 for i in range(n):
  B[n][i]=1;B[n+1+i][n]=1
  for j in range(n):
   B[i][j]=A[i][j];B[n+1+i][n+1+j]=A[j][i]
   B[i][n+1+j]=int(i==j or A[i][j]);B[n+1+j][i]=1-B[i][n+1+j]
 return B
def refine(A,marks=(),colors=None):
 n=len(A);c=list(colors) if colors is not None else [marks.index(i)+1 if i in marks else 0 for i in range(n)];history=[]
 while True:
  cells=[[v for v in range(n)if c[v]==k]for k in sorted(set(c))];history.append(sorted(map(len,cells)))
  records=[(c[v],tuple(sum(A[v][u]for u in C)for C in cells))for v in range(n)]
  ix={r:i for i,r in enumerate(sorted(set(records)))};d=[ix[r]for r in records]
  if len(set(d))==len(set(c)):return cells,history,c
  c=d
