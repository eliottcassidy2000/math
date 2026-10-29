#!/usr/bin/env python3
"""Exact weighted partition versus joint-observation algebra controls."""
from fractions import Fraction
from itertools import product
import random,json
CHECKS=0
def need(p,msg):
 global CHECKS
 CHECKS+=1
 if not p:raise AssertionError(msg)
def image(A,f):return tuple(sum(a*x for a,x in zip(row,f)) for row in A)
def addbasis(B,v):
 v=list(map(Fraction,v))
 for pivot,b in B:
  c=v[pivot]
  if c:v=[x-c*y for x,y in zip(v,b)]
 pivot=next((i for i,x in enumerate(v) if x),None)
 if pivot is None:return False
 c=v[pivot];v=tuple(x/c for x in v)
 B.append((pivot,v));B.sort(key=lambda z:z[0]);return True
def closure(As,outputs,joint=True):
 n=len(As[0]);B=[]
 for v in [(1,)*n,*outputs]:addbasis(B,v)
 while True:
  old=len(B);vs=[b for _,b in B]
  for A in As:
   for v in vs:addbasis(B,image(A,v))
  if joint:
   for v in vs:
    for w in vs:addbasis(B,tuple(x*y for x,y in zip(v,w)))
  if len(B)==old:return [v for _,v in B]
def cells_from(outputs,n):
 signatures=[tuple(f[i] for f in outputs) for i in range(n)]
 return tuple(tuple(j for j in range(n) if signatures[j]==signatures[i]) for i in range(n))
def refine(As,outputs):
 n=len(As[0]);cells=cells_from(outputs,n)
 while True:
  blocks=sorted(set(cells));sigs=[]
  for i in range(n):sigs.append((cells[i],tuple(sum(A[i][j] for j in C) for A in As for C in blocks)))
  nxt=tuple(tuple(j for j in range(n) if sigs[i]==sigs[j]) for i in range(n))
  if nxt==cells:return cells
  cells=nxt
def test(As,outputs):
 n=len(As[0]);B=closure(As,outputs);pi=refine(As,outputs)
 need(cells_from(B,n)==pi,'algebra versus refinement')
 need(len(B)==len(set(pi)),'dimension versus cells')
 for A in As:
  for C in set(pi):
   for block in set(pi):need(len({sum(A[i][j] for j in C) for i in block})==1,'equitable')
 return len(B),len(closure(As,outputs,False))
def main():
 cases=0
 for n in range(1,4):
  for bits in range(1<<(n*n)):
   A=[[(bits>>(i*n+j))&1 for j in range(n)] for i in range(n)]
   for f in range(1<<n):test([A],[tuple((f>>i)&1 for i in range(n))]);cases+=1
 n=4;pairs=[(i,j) for i in range(n) for j in range(i)]
 for bits in range(1<<len(pairs)):
  A=[[0]*n for _ in range(n)]
  for k,(i,j) in enumerate(pairs):A[i][j]=(bits>>k)&1
  for f in range(1<<n):test([A],[tuple((f>>i)&1 for i in range(n))]);cases+=1
 rng=random.Random(92001)
 for _ in range(120):
  n=rng.randrange(4,9);As=[[[rng.randrange(-2,3) for j in range(n)] for i in range(n)] for _ in range(2)]
  test(As,[tuple(rng.randrange(3) for i in range(n))]);cases+=1
 # Six-state shared-branch correlation.
 A=[[0]*6 for _ in range(6)]
 A[0][2]=A[0][5]=A[1][3]=A[1][4]=1
 level=(1,1,0,0,0,0);a=(0,0,0,0,1,1);b=(0,0,0,1,0,1)
 pa=refine([A],[level,a]);pb=refine([A],[level,b])
 common=tuple(tuple(j for j in range(6) if pa[i]==pa[j] and pb[i]==pb[j]) for i in range(6))
 need(len(set(pa))==len(set(pb))==3,'separate quotient sizes')
 need(len(set(common))==5,'joint compatible cells')
 d,l=test([A],[level,a,b]);need((d,l)==(6,4),'joint nonlinear cost')
 ab=image(A,tuple(x*y for x,y in zip(a,b)));need(ab[:2]==(1,0),'correlation witness')
 # Seven-state equal linear continuation, unequal second moment.
 A=[[0]*7 for _ in range(7)]
 A[0][2]=A[0][3]=1;A[1][4]=2;A[2][5]=2;A[3][6]=2;A[4][5]=A[4][6]=1
 f=(0,0,0,0,0,1,0);levels=[tuple(int(j in C) for j in range(7)) for C in [(0,1),(2,3,4),(5,6)]]
 B=closure([A],[f,*levels],False);need(all(v[0]==v[1] for v in B),'all linear responses agree')
 g=image(A,f);second=image(A,tuple(x*x for x in g));need(second[:2]==(4,2),'second moment separates')
 joint=refine([A],[f,*levels]);need(joint[0]!=joint[1],'stable partition separates')
 print(json.dumps({'status':'PASS','checks':CHECKS,'universe_cases':cases,'correlation_dimensions':[d,l],'variance_values':list(second[:2])},sort_keys=True))
if __name__=='__main__':main()
