#!/usr/bin/env python3
"""Actual five-letter grammar: exact refinements and hostile isomorphism checks."""
from pathlib import Path
from collections import Counter
from itertools import permutations
import json
HERE=Path(__file__).resolve().parent
def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def mm(A,B):return [[sum(x*y for x,y in zip(r,c))for c in zip(*B)]for r in A]
def power(A,n):
 R=[[int(i==j)for j in range(len(A))]for i in range(len(A))]
 for _ in range(n):R=mm(R,A)
 return R
C=[[1,0,1],[1,1,1],[1,1,0]];K=[[1,0,1],[1,1,1],[1,1,1]]
def trace(M):return sum(M[i][i]for i in range(len(M)))
trans={'A':{'a':'A','b':'B'},'B':{'c':'A','d':'B','e':'B'}}
layer={'':'B'};rows=[]
# Coefficients of reciprocal D(t); independent scalar recurrence path.
g=[1,2,5]
for j in range(3,15):g.append(2*g[-1]+g[-2]-g[-3])
q=[1,3]
for j in range(2,15):q.append(3*q[-1]-q[-2])
for m in range(11):
 occ=Counter();first=Counter()
 for w in layer:
  hits=[i+1 for i in range(1,m)if w[i-1:i+1] in ('be','ee')]
  occ[len(hits)]+=1
  if hits:first[hits[0]]+=1
 require(occ[0]==sum(power(C,m)[1]),'safe transfer mismatch')
 require(len(layer)==q[m],'total grammar mismatch')
 for k in range(2,m+1):require(first[k]==g[k-2]*q[m-k],'first hit mismatch')
 # Evaluate the occurrence polynomial at several integers against transfer.
 for u in (0,1,2,3):
  M=[r[:]for r in C];M[2][2]=u
  require(sum(n*u**r for r,n in occ.items())==sum(power(M,m)[1]),'occurrence transfer mismatch')
 rows.append({'length':m,'total':len(layer),'occurrences':dict(sorted(occ.items())),'first_hit_end':dict(sorted(first.items()))})
 layer={w+a:t for w,state in layer.items()for a,t in trans[state].items()}
# Independent coefficients of the exact-r rational series.
from math import comb
def conv(a,b):return [sum(a[j]*b[k-j]for j in range(k+1)if j<len(a) and k-j<len(b))for k in range(11)]
for r in range(1,11):
 gp=[1]+[0]*10
 for _ in range(r+1):gp=conv(gp,g[:11])
 numerator=[0]*(r+1)+[(-1)**j*comb(2*r-2,j)for j in range(2*r-1)]
 coeff=conv(numerator,gp)
 for m in range(11):require(coeff[m]==rows[m]['occurrences'].get(r,0),'exact-r GF mismatch')
S=[[1,0,1],[1,2,0],[1,1,0]];T=[[1,1,1],[1,1,0],[1,0,0]]
require(mm(C,S)==mm(S,T),'intertwiner mismatch')
require(S[1]==[1,2,0] and [row[0]for row in S]==[1,1,1],'boundary intertwiner mismatch')
# Cyclic grammar words: not seed-B insertion words or tournament classes.
from math import gcd
cyclic_rows=[]
for m in range(1,9):
 cyclic_words=set()
 for start in ('A','B'):
  paths={'':start}
  for _ in range(m):paths={w+a:t for w,st in paths.items()for a,t in trans[st].items()}
  cyclic_words.update(w for w,st in paths.items()if st==start and all(w[i]+w[(i+1)%m] not in ('be','ee')for i in range(m)))
 require(len(cyclic_words)==trace(power(C,m)),'cyclic trace mismatch')
 necklaces={min(w[j:]+w[:j]for j in range(m))for w in cyclic_words}
 burnside=sum(sum(gcd(k,d)==1 for k in range(1,d+1))*trace(power(C,m//d))for d in range(1,m+1)if m%d==0)//m
 require(len(necklaces)==burnside,'cyclic Burnside mismatch')
 cyclic_rows.append({'length':m,'rooted':len(cyclic_words),'necklaces':len(necklaces)})
# Native code hostile records; fully independent permutation isomorphism test.
records=json.loads((HERE/'source_rotation_witness.json').read_text())['records']
def decode(bits):
 n=(1+__import__('math').isqrt(1+8*len(bits)))//2;A=[[0]*n for _ in range(n)];k=0
 for i in range(n):
  for j in range(i+1,n):A[i][j]=int(bits[k]);A[j][i]=1-int(bits[k]);k+=1
 return A
def iso(A,B):
 return any(all(A[i][j]==B[p[i]][p[j]]for i in range(len(A))for j in range(len(A)))for p in permutations(range(len(A)))if all(sum(A[i])==sum(B[p[i]])for i in range(len(A))))
hostiles=[]
for left,right,operation in [('ecb','cbe','rotation'),('de','ed','reversal')]:
 A,B=decode(records[left]['bits']),decode(records[right]['bits']);require(not iso(A,B),'unexpected isomorphism')
 hostiles.append({'left':left,'right':right,'operation':operation,'isomorphic':False,'order':len(A)})
out={'scope':'Actual classical grammar, not an all-order matrix nonsingularity proof','rows':rows,'cyclic_grammar':cyclic_rows,'hostiles':hostiles}
(HERE/'actual_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
