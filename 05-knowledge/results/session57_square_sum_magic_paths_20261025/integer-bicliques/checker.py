#!/usr/bin/env python3
"""Independent finite certificate checker; errors are reported with nonzero exit."""
import math,itertools,json,sys,traceback

def require(condition,message):
 if not condition:raise ValueError(message)

def common_factor(a,b,N):
 require(type(a) is int and type(b) is int and type(N) is int and 1<=a<b<=N,'invalid pair/universe')
 d=b-a;out=[]
 for p in range(1,math.isqrt(d)+1):
  if d%p:continue
  q=d//p
  if p>=q or (q-p)%2:continue
  u=(q-p)//2;x=u*u-a
  if 1<=x<=N and x not in (a,b):out.append(x)
 return sorted(out)

def roots_witness(L,R):
 require(len(set(L+R))==len(L)+len(R),'vertices collide')
 require(all(type(v) is int and v>0 for v in L+R),'invalid vertices')
 roots=[[math.isqrt(a+x) for x in R] for a in L]
 require(all(r*r==a+x for a,row in zip(L,roots) for x,r in zip(R,row)),'nonsquare cross-sum')
 return roots

def run():
 # All possible K2,4 with maximum vertex at most 274 have a unique sorted left pair.
 best=None;pair_count=0
 N=274
 for a in range(1,N):
  for b in range(a+1,N+1):
   C=common_factor(a,b,N);pair_count+=1
   brute=[x for x in range(1,N+1) if x not in (a,b) and math.isqrt(x+a)**2==x+a and math.isqrt(x+b)**2==x+b]
   require(C==brute,f'factor disagreement {a,b}')
   if len(C)>=4:
    item=(max(b,C[3]),a,b,C[:4])
    if best is None or item<best:best=item
 require(best==(274,15,255,[1,34,106,274]),f'wrong K24 minimum {best}')
 roots24=roots_witness([15,255],[1,34,106,274])
 # Independent graph definition by testing all pairs, not by enumerating square sums.
 N=673;A=[set() for _ in range(N+1)]
 for a in range(1,N):
  for b in range(a+1,N+1):
   if math.isqrt(a+b)**2==a+b:A[a].add(b);A[b].add(a)
 min33=None;eligible=0;tested=set()
 for a in range(1,N):
  for b in range(a+1,N+1):
   C=common_factor(a,b,N)
   require(set(C)==A[a]&A[b],f'adjacency mismatch {a,b}')
   if len(C)<3:continue
   eligible+=1
   for R in itertools.combinations(C,3):
    if R in tested:continue
    tested.add(R);L=sorted(A[R[0]]&A[R[1]]&A[R[2]])
    if len(L)>=3:
     item=(max(L[2],R[-1]),L[:3],list(R))
     if min33 is None or item<min33:min33=item
 require(min33==(673,[1,193,673],[3,168,483]),f'wrong K33 minimum {min33}')
 roots33=roots_witness(min33[1],min33[2])
 # Regression tests: parity obstruction, forbidden overlap, original positive control.
 require(all(common_factor(a,a+2,100)==[] for a in range(1,99)),'mod 4 obstruction failed')
 require(8 not in common_factor(8,248,300),'self-loop admitted')
 roots_witness([1,97],[3,24,99,528])
 for args in [(0,2,3),(2,1,3),(1,1,3),(1,5,4)]:
  try:common_factor(*args)
  except ValueError:pass
  else:raise ValueError('invalid input accepted')
 return {'status':'PASS','k24_minimum':274,'k24_roots':roots24,'k24_pairs_verified':pair_count,'k33_minimum':673,'k33_roots':roots33,'k33_seed_pairs':eligible,'k33_distinct_right_triples':len(tested)}
if __name__=='__main__':
 try:
  result=run();print(json.dumps(result,indent=2))
 except Exception as e:
  print(json.dumps({'status':'FAIL','type':type(e).__name__,'message':str(e)}),file=sys.stderr);sys.exit(1)
