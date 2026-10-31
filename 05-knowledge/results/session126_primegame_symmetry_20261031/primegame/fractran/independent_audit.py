#!/usr/bin/env python3
"""Independent scratch audit: no imports from the lane's verifier."""
from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import prod
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent
P=(2,3,5,7,11,13,17,19,23,29)
F=tuple(map(Fraction, ('17/91','78/85','19/51','23/38','29/33','77/29','95/23','77/19','1/17','11/13','13/11','15/2','1/7','55/1')))
checks=0

def check(b,msg):
 global checks
 checks+=1
 if not b: raise ValueError(msg)

def nxt(n):
 for i,f in enumerate(F):
  y=n*f
  if y.denominator==1: return y.numerator,i
 raise ValueError('no next')

def val(e):return prod(p**a for p,a in zip(P,e))

def fac(n):
 e=[]
 for p in P:
  a=0
  while n%p==0:n//=p;a+=1
  e.append(a)
 check(n==1,'outside prime basis')
 return e

def controller(e):
 a,b,c,d=e[:4]
 s=next((p for p,x in zip(P[4:],e[4:]) if x),1)
 if s==1:
  if a:return [a-1,b+1,c+1,d]+[0]*6,11
  if d:return [a,b,c,d-1]+[0]*6,12
  return [a,b,c+1,d,1,0,0,0,0,0],13
 # Independent tabulated finite-state edges, guards in priority order.
 table={13:[(d,[0,0,0,-1],17,0),(True,[0,0,0,0],11,9)],17:[(c,[1,1,-1,0],13,1),(b,[0,-1,0,0],19,2),(True,[0,0,0,0],1,8)],19:[(a,[-1,0,0,0],23,3),(True,[0,0,0,1],11,7)],11:[(b,[0,-1,0,0],29,4),(True,[0,0,0,0],13,10)],29:[(True,[0,0,0,1],11,5)],23:[(True,[0,0,1,0],19,6)]}
 for enabled,diff,t,i in table[s]:
  if enabled:return [x+y for x,y in zip((a,b,c,d),diff)]+[int(p==t) for p in P[4:]],i

# Exhaustive one-step guard audit on all 7*3^4=567 legal control/counter states.
for a,b,c,d in product(range(3),repeat=4):
 for s in (1,11,13,17,19,23,29):
  e=[a,b,c,d]+[int(p==s) for p in P[4:]]
  y,i=nxt(val(e));g,j=controller(e)
  check((y,i)==(val(g),j),'controller mismatch')
  check(sum(fac(y)[4:])<=1,'one-hot not preserved')

# Check all stated trace rows using Fraction arithmetic, not denominator divisibility code.
pth=ROOT/'primegame_trace.tsv'
rows=pth.read_text().splitlines()
check(len(rows)==36983,'trace row count')
previous=None;outputs=[];maxbits=0
for line in rows[1:]:
 t,n,enc,j=line.split('\t');t,n=int(t),int(n); e=list(map(int,enc.split(',')))
 check(n==val(e),'trace value mismatch')
 if previous is not None:check(n==previous,'trace successor mismatch')
 if j:
  previous,i=nxt(n);check(i+1==int(j),'chosen instruction mismatch')
 else:previous=None
 if t and n&(n-1)==0:outputs.append([n.bit_length()-1,t])
 maxbits=max(maxbits,n.bit_length())
check(outputs==[[2,19],[3,69],[5,281],[7,710],[11,2375],[13,3893],[17,8102],[19,11361],[23,19268],[29,36981]],'output table mismatch')
check(maxbits==150,'max bit count')

# Integer divisibility-test macros via Fraction transitions; first return to canonical form.
macros=updates=0
for n in range(2,51):
 for d in range(1,n):
  state=5**n*7**d*13;word=[];hits=[];q,r=divmod(n,d)
  for t in range(1,100000):
   state,i=nxt(state);word.append(chr(65+i))
   if state&(state-1)==0:hits.append(state)
   e=fac(state)
   if e[0]==e[1]==0 and e[4:]==[0,1,0,0,0,0]:break
  expected=(5**n*7**(d-1)*13 if r else 5**(n+1)*7**n*13)
  check(state==expected,'macro endpoint')
  check(hits==([2**n] if d==1 else []),'macro pure-power outputs')
  prefix=('AB'*d+'J'+'EF'*d+'K')*q+'AB'*r+'A'
  suffix=('C'+'DG'*n+'H'+'EF'*(r-1)+'K' if r else 'I'+'L'*n+'M'*(d-1)+'N'+'EF'*n+'K')
  check(''.join(word)==prefix+suffix,'macro instruction word')
  check(t==(6*n+2*q+2 if r else 7*n+2*q+d+3),'macro step formula')
  macros+=1;updates+=t

# Determinant by fraction-free Bareiss elimination.
deltas=[[a-b for a,b in zip(fac(f.numerator),fac(f.denominator))] for f in F]
indices=[8,12,9,0,13,4,5,11,2,3]
A=[[deltas[j][i] for j in indices] for i in range(10)]
sign=1;old=1
for k in range(9):
 p=next(i for i in range(k,10) if A[i][k])
 if p!=k:A[p],A[k]=A[k],A[p];sign=-sign
 pivot=A[k][k]
 for i in range(k+1,10):
  for j in range(k+1,10):
   num=pivot*A[i][j]-A[i][k]*A[k][j]
   check(num%old==0,'Bareiss noninteger division');A[i][j]=num//old
 for i in range(k+1,10):A[i][k]=0
 old=pivot
det=sign*A[-1][-1];check(det==-1,'minor not unimodular')

# Exact known seed witnesses including count-shifts.
for seed,witness in [(7,[7,1,55,65,55]),(17,[17,1,55,65,55]),(55,[55,65,55]),(25,[25,1375,1625,1375]),(23,[23,95,385,455,85,78,66,58,154,182,34,2])]:
 for x,y in zip(witness,witness[1:]):check(nxt(x)[0]==y,'seed witness')
check(nxt(2)[0]==15,'seed15 join')
check(nxt(8)[0]==60,'seed60 join')
n=16;outs=[]
for t in range(1,2245):
 n,_=nxt(n)
 if n&(n-1)==0:outs.append([n.bit_length()-1,t])
check(outs==[[5,150],[7,579],[11,2244]],'seed16 outputs')

result={'status':'PASS','independent_method':'Fraction arithmetic for all trace and macro transitions; separate finite-control table; Bareiss determinant','trace_sha256':sha256(pth.read_bytes()).hexdigest(),'trace_updates':36981,'controller_states_checked':567,'macro_pairs_checked':macros,'macro_updates':updates,'unimodular_minor_determinant':det,'explicit_checks':checks,'defects':['Minor wording only: replace at most one control prime occurs to exponent one by all six control exponents are zero except possibly one equal to one.'],'mathematical_review':['FUS at 2^ell: tau=ell+1,omega=1,nonprime and nonsquarefree for ell>=2, hence defect ell-3.','Simultaneous prime relabeling is conjugacy, not fixed-rack coordinate symmetry.','CA construction is valid as variable-time single-head tape simulation, not a direct one-step local interpretation or a tested explicit rule.']}
(ROOT/'independent_audit_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
