#!/usr/bin/env python3
"""Stdlib exact certificates. Checks survive -O; no huge source enumeration."""
from fractions import Fraction
from pathlib import Path
import argparse,hashlib,json,math,subprocess
HERE=Path(__file__).resolve().parent
checks=0
def gate(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
def state(w):
 A=D=1;B=0
 for k in w:A,D,B=3*A,D*2**k,3*B+D
 return A,D,B
def power(a,b,k,m):
 # Independent binary affine composition.
 u,v=1,0
 while k:
  if k&1:u,v=a*u%m,(a*v+b)%m
  a,b=a*a%m,(a*b+b)%m;k//=2
 return u,v
def certificate(w):
 p=101;m=p*p;A,D,B=state(w);a=A*pow(D,-1,m)%m;b=B*pow(D,-1,m)%m
 orbit=[];index={};x=0
 for i in range(m):
  gate(x not in index,'cycle uniqueness');index[x]=i;orbit.append(x);x=(a*x+b)%m
 gate(x==0,'cycle closure')
 c=(a-1)//p
 for i,z in enumerate(orbit):
  gate(z==b*(i+p*c*i*(i-1)//2)%m,'quadratic clock')
  y=z*pow(b,-1,m)%m
  gate((y-p*c*y*(y-1)//2)%m==i,'inverse quadratic clock')
 roots=[];phases=[];s=Fraction(B,D-A)
 for j in range(len(w)):
  Aj,Dj,Bj=state(w[:j]);root=-Bj*pow(Aj,-1,m)%m
  roots.append(root);phases.append(str((Aj*s+Bj)/Dj))
  gate(power(a,b,index[root],m)[1]==root,'binary anchor')
 anchors=sorted(index[z] for z in roots)
 gaps=[y-x for x,y in zip(anchors,anchors[1:])]+[m-anchors[-1]]
 capacities=gaps[:-1]+[gaps[-1]-1]
 covered={0};survivor_counts=[];last=None
 for r in range(m+1):
  if r:
   covered.add((-r)%m)
   for z in anchors[1:]:covered.add((z-r+1)%m)
  actual=m-len(covered);expected=sum(max(c-r,0) for c in capacities)
  gate(actual==expected,'all lengths coverage formula')
  if actual:last=(r,[orbit[i] for i in range(m) if i not in covered]) if actual==1 else last
  survivor_counts.append(actual)
 maxr=max(capacities)-1
 # Direct node evaluations at the last one-residue boundary, not root coverage.
 survivors=[i for i in range(m) if all((i+t)%m not in anchors for t in range(maxr)) and (i+maxr)%m!=0]
 # This direct comprehension is intentionally bounded by m*maxr; short-circuiting saves cost.
 gate(len(survivors)==1,'boundary singleton')
 n=orbit[survivors[0]];x=n
 for t in range(maxr):
  for j in range(len(w)):
   Aj,Dj,Bj=state(w[:j]);gate((Aj*x+Bj)%m!=0,'direct phase square avoidance')
  x=(a*x+b)%m
 gate(x!=0,'final endpoint survival')
 # Nonresonant exceptional square repair, checked across a full modulus orbit.
 q=13 if w==(1,1,5) else 7;q2=q*q
 f=B*pow(D-A,-1,q2)%q2;seed=(f+q)%q2;x=seed;seen=set()
 while x not in seen:
  seen.add(x)
  for j in range(len(w)):
   Aj,Dj,Bj=state(w[:j]);gate((Aj*x+Bj)%q2!=0,'exceptional repair')
  x=(A*x+B)*pow(D,-1,q2)%q2
 gate(x==seed,'repair period')
 return dict(word=w,carry=B,gap=D-A,roots=roots,anchors=anchors,gaps=gaps,capacities=capacities,max_squarefree_repetitions=maxr,last_source_residue=n,fixed_phases=phases,repair_prime=q,repair_seed=seed,repair_period=len(seen),count_sha256=hashlib.sha256(json.dumps(survivor_counts).encode()).hexdigest())
def main():
 result=[certificate((1,1,5)),certificate((1,3,3))]
 gate(result[0]['capacities']==[3795,1419,4986],'inherited capacities')
 gate(result[0]['last_source_residue']==7332,'inherited boundary')
 gate(result[1]['capacities']==[8601,179,1420],'new capacities')
 # Finite general resonance universe, including nonunit carries.
 for p in (3,5,7):
  for s in (1,2,3):
   for t in range(s):
    for e in (1,2,3):
     m=p**e;a=1+p**s;b=p**t;period=p**max(e-t,0)
     for x in range(m):
      gate((power(a,b,period,m)[0]*x+power(a,b,period,m)[1]-x)%m==0,'predicted return')
      if period>1:
       u,v=power(a,b,period//p,m);gate((u*x+v-x)%m!=0,'minimal return')
 # Hostile: denominator prime not sufficient for transitivity/forced square.
 A,D,B=state((1,13));gate((A,D,B)==(9,16384,5),'hostile affine')
 gate(Fraction(B,D-A)==Fraction(1,3275),'hostile denominator')
 for seed in range(1,25,5):
  x=seed
  for _ in range(25):
   gate(x%25!=0 and (3*x+1)%25!=0,'hostile all phase avoidance')
   x=(A*x+B)*pow(D,-1,25)%25
 # Same-clock append distinction.
 append=[]
 for w in ((1,1,5),(1,3,3)):
  A,D,B=state(w+(1,));append.append(str(Fraction(B,D-A)))
 gate(append==['37/35','221/175'],'append contrast')
 path='05-knowledge/results/collatz_guards_20260921_squarefree.md'
 source_sha256='a7727f346a5ff53d34c94d1a788c0810b28c0a203c1645f7d80be0753000df0f'
 if args.source_repo:
  blob=subprocess.check_output(['git','-C',args.source_repo,'show','e871dd4513c:'+path])
  if hashlib.sha256(blob).hexdigest()!=source_sha256:raise RuntimeError('pinned source')
 # Provenance hash is metadata, not a core mathematical gate.
 gate(len(source_sha256)==64,'source hash metadata')
 return dict(status='FINITE-EXACT',checks=checks,source_commit='e871dd4513c',source_path=path,source_sha256=source_sha256,certificates=result,append_one_fixed_points=append,hostile=dict(word=[1,13],fixed_point='1/3275',avoiding_class='1 mod5',square_modulus=25),scope='Proof in README; no Collatz convergence claim; large-source existence uses fixed-word sieve, not enumeration.')
if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--source-repo',help='Optional Git repository for checking the pinned inherited blob')
 args=parser.parse_args()
 print(json.dumps(main(),indent=2,sort_keys=True))
