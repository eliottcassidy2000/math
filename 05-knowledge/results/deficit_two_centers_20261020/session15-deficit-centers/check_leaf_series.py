#!/usr/bin/env python3
import json
from pathlib import Path
N=200
A=[0]*(N+1);A[0]=1
for n in range(1,N+1):A[n]=A[n-1]+(A[n-3] if n>=3 else 0)
a=lambda n:A[n] if n>=0 else 0
M=[a(n//2-2) if n%2==0 else 0 for n in range(N+1)]
L=[(a(n)-a(n//2))//2+M[n] if n%2==0 else (a(n)-a(n//2)-a(n//2-1))//2 for n in range(N+1)]
# Separate power-series inversion of the expanded rational denominator.
Q={0:1,1:-1,2:-1,5:1,6:-1,7:1,9:1};P={4:2,5:-1,7:-1}
R=[]
for n in range(N+1):R.append(P.get(n,0)-sum(c*R[n-k] for k,c in Q.items() if k and k<=n))
assert R==L
for r in range(5,N//2+1):assert a(r)-2*a(r-2)==a(r-7);assert L[2*r]==(a(2*r)-a(r-7))//2
for n in range(5,N+1):assert M[n]==(M[n-2] if n>=2 else 0)+(M[n-6] if n>=6 else 0)
rows=[dict(n=n,rigid=A[n],mixed_leaves=M[n],all_leaves=L[n]) for n in range(N+1)]
Path(__file__).with_name('leaf_series_coefficients.json').write_text(json.dumps(dict(checked_through=N,rows=rows),indent=2)+'\n')
print('PASS through',N);print('all leaves0..30:',L[:31]);print('mixed leaves0..30:',M[:31])
