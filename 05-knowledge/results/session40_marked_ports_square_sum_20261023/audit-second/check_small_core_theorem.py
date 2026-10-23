#!/usr/bin/env python3
import itertools,json
from pathlib import Path

def matrix(n,mask):
 A=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(itertools.combinations(range(n),2)):A[i][j]=(mask>>k)&1;A[j][i]=1-A[i][j]
 return A
def code(A,p):return tuple(A[p[i]][p[j]] for i in range(len(p)) for j in range(i+1,len(p)))
def require(b,msg):
 if not b:raise RuntimeError(msg)
regular=[]
for m in range(1024):
 A=matrix(5,m)
 if set(map(sum,A))=={2}:regular.append(A)
canon={min(code(A,p) for p in itertools.permutations(range(5))) for A in regular}
require(len(regular)==24 and len(canon)==1,'regular five uniqueness')
C=[[int((j-i)%3==1) for j in range(3)] for i in range(3)];blocks=[]
for m in range(512):
 B=[[(m>>(3*i+j))&1 for j in range(3)] for i in range(3)];rs=list(map(sum,B));cs=[sum(B[i][j] for i in range(3)) for j in range(3)]
 if len(set(rs+cs))!=1:continue
 q=rs[0]
 if q in(1,2):
  pi=[next(j for j in range(3) if B[i][j]==(q==1)) for i in range(3)]
  J=sum(2*C[pi[i]][pi[j]]-1 for i in range(3) for j in range(3) if C[i][j])
  converse=sum(2*(1-C[pi[i]][pi[j]])-1 for i in range(3) for j in range(3) if i!=j and not C[i][j])
  require(J in(-3,3) and J==converse,'shared chirality converse')
  gauge={tuple(B[(i+a)%3][(j+b)%3] for i in range(3) for j in range(3)) for a in range(3) for b in range(3)}
  require(len(gauge)==3,'rotation orbit')
 else:J=0
 blocks.append((q,J))
require(len(blocks)==14,'regular mixed blocks')
H=[]
for steps in({1,2,3},{1,2,4}):
 A=[[int((j-i)%7 in steps) for j in range(7)] for i in range(7)]
 h=sum(all(A[p[i]][p[i+1]] for i in range(6)) for p in itertools.permutations(range(7)))
 H.append(h);require(all(A[i][i+1] for i in range(6)),'common order-nine path')
 require(all(A[i][j]==A[(-j)%7][(-i)%7] for i in range(7) for j in range(7)),'self converse')
require(H==[175,189],'seven-cell hostile')
out={'status':'PASS','regular_five_presentations':len(regular),'regular_five_isotypes':len(canon),'regular_three_by_three_blocks':len(blocks),'seven_cell_H':H,'mixed_chirality_cases':12}
Path(__file__).with_name('small-core-audit-receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
