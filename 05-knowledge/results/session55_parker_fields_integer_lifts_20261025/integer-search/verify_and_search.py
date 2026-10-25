#!/usr/bin/env python3
"""Exact certificates and complete bounded searches. Standard library only."""
import argparse, hashlib, json, math
from pathlib import Path
def require(condition, message="certificate invariant failed"):
 if not condition:raise RuntimeError(message)

def verify_payload(expected, observed):
 require(expected==observed,"certificate data mismatch")

def tamper_selftest(data):
 import copy
 outcomes={}
 for name in ("changed_content","missing_key","extra_key"):
  bad=copy.deepcopy(data)
  if name=="changed_content":bad["schema"]="tampered"
  elif name=="missing_key":bad.pop("schema")
  else:bad["unexpected_extra_key"]=True
  try:verify_payload(data,bad)
  except RuntimeError:outcomes[name]="rejected"
  else:raise RuntimeError("tamper test accepted "+name)
 return outcomes

LINES=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
def issquare(n): return n>=0 and math.isqrt(n)**2==n
def prime(p): return p>=2 and all(p%d for d in range(2,math.isqrt(p)+1))
def grid(c,x,y): return [c+x,c-x-y,c+y,c-x+y,c,c+x-y,c-y,c+x+y,c-x]
def assess(a,p):
 s=[sum(a[i] for i in line) for line in LINES]
 return dict(entries=a,line_sums=s,distinct=len(set(a))==9,nonnegative=min(a)>=0,integer_square_entries=[issquare(z) for z in a],modular_square_entries=[any(r*r%p==z%p for r in range(p)) for z in a],modular_magic=len({z%p for z in s})==1,integer_magic=len(set(s))==1,no_wrap=max(s)<p,modulus=p,canonical_representatives=all(0<=z<p for z in a),least_modular_roots=[next((r for r in range(p) if r*r%p==z%p),None) for z in a])
def enumerate_offsets(D,c):
 D=set(D); ordered=sorted(D)
 for x in ordered:
  for y in ordered:
   if y>=x:break
   if len({x,y,x-y,x+y})==4 and x-y in D and x+y in D:
    yield grid(c,x,y)
def weak_prime_search(bound=59):
 rows=[]
 for p in range(2,bound+1):
  if not prime(p):continue
  Q={r*r%p for r in range(p)}; witnesses=[]; tested=0
  for c in sorted(Q):
   if not 0<3*c<p:continue
   D={d for d in range(1,c+1) if c+d in Q and c-d in Q}
   tested+=len(D)*(len(D)-1)//2
   witnesses.extend(enumerate_offsets(D,c))
  require({tuple(z) for z in witnesses}==independent_weak_grids(p), 'certificate invariant failed')
  rows.append(dict(p=p,oriented_pair_tests=tested,canonical_grids=len(witnesses),first=witnesses[0] if witnesses else None,independent_corner_enumeration_agrees=True))
 return rows

def independent_weak_grids(p):
 # Different universe: all ordered choices of two corner entries from Q.
 Q={r*r%p for r in range(p)}; out=set()
 for c in Q:
  if not 0<3*c<p:continue
  for a in Q:
   for b in Q:
    z=grid(c,a-c,b-c)
    if min(z)<0 or max(z)>=p or len(set(z))!=9 or not all(v in Q for v in z):continue
    x,y=sorted((abs(a-c),abs(b-c)),reverse=True)
    out.add(tuple(grid(c,x,y)))
 return out

def divisor_count_products(B):
 # Gaussian representation formula, independent of Euclidean enumeration.
 spf=list(range(B+1))
 for p in range(2,math.isqrt(B)+1):
  if spf[p]!=p:continue
  for v in range(p*p,B+1,p):
   if spf[v]==v:spf[v]=p
 counts=[0]*(B+1)
 for e in range(1,B+1):
  v=e;prod=1
  while v>1:
   p=spf[v];k=0
   while v%p==0:v//=p;k+=1
   if p%4==1:prod*=2*k+1
  counts[e]=(prod-1)//2
 return counts

def sieve_offsets(B):
 """All positive d with e^2-d and e^2+d squares, 1<=e<=B."""
 banks=[None]*(B+1); primitive=0; insertions=0
 for m in range(2,math.isqrt(B)+1):
  for n in range(1,m):
   e0=m*m+n*n
   if e0>B:break
   if (m-n)%2==0 or math.gcd(m,n)!=1:continue
   primitive+=1;d0=4*m*n*(m*m-n*n)
   for k in range(1,B//e0+1):
    e=k*e0
    if banks[e] is None:banks[e]=set()
    banks[e].add(d0*k*k);insertions+=1
 return banks,primitive,insertions

def brute_offsets(e):
 c=e*e
 return {c-a*a for a in range(e) if issquare(2*c-a*a)}

def integer_search(B,independent_bound):
 banks,primitive,insertions=sieve_offsets(B)
 counts=divisor_count_products(B)
 histogram={}; tests=0; solutions=[]; digest=hashlib.sha256(); maxcount=0; maxcenters=[]
 for e in range(1,B+1):
  D=banks[e] or set(); n=len(D)
  require(n==counts[e], (e,n,counts[e]))
  if e<=independent_bound:require(D==brute_offsets(e), (e,D,brute_offsets(e)))
  for d in D:require(0<d<e*e and issquare(e*e-d) and issquare(e*e+d), 'certificate invariant failed')
  histogram[str(n)]=histogram.get(str(n),0)+1
  digest.update((str(e)+':'+','.join(map(str,sorted(D)))+'\n').encode())
  if n>maxcount:maxcount=n;maxcenters=[e]
  elif n==maxcount:maxcenters.append(e)
  # Four distinct magnitudes are necessary. No primitivity filter is used.
  if n<4:continue
  tests+=n*(n-1)//2
  for a in enumerate_offsets(D,e*e):
   require(all(issquare(z) for z in a) and len(set(a))==9, 'certificate invariant failed')
   require(len({sum(a[i] for i in l) for l in LINES})==1, 'certificate invariant failed')
   solutions.append(a)
 return dict(status='FINITE-EXACT',center_root_bound=B,inclusive=True,includes_nonprimitive=True,includes_zero_entries=True,primitive_pythagorean_parameters=primitive,sieve_insertions=insertions,offset_count_histogram=histogram,eligible_pair_tests=tests,maximum_offset_count=maxcount,first_ten_maximum_centers=maxcenters[:10],solutions=solutions,offset_banks_sha256=digest.hexdigest(),independent_gaussian_count_comparison_through=B,independent_brute_comparison_through=min(B,independent_bound),excluded_no_wrap_moduli_through=3*B*B)

def controls():
 require(list(enumerate_offsets({1,2,3,4},5))==[[8,1,6,3,5,7,4,9,2]], 'certificate invariant failed')
 require(not list(enumerate_offsets({1,2,3},5)), 'certificate invariant failed')
 weak=assess([28,3,17,5,16,27,15,29,4],59)
 lo=assess([8,1,6,3,5,7,4,9,2],311)
 cain=assess([r*r%29 for r in [9,11,1,6,0,14,12,16,8]],29)
 parker=assess([r*r for r in [29,1,47,41,37,1,23,41,29]],5003)
 constant=assess([4]*9,29)
 for a in [weak,lo]:require(a['distinct'] and a['modular_magic'] and a['no_wrap'] and all(a['modular_square_entries']) and not all(a['integer_square_entries']), 'certificate invariant failed')
 require(cain['distinct'] and cain['modular_magic'] and not cain['integer_magic'], 'certificate invariant failed')
 require(all(parker['integer_square_entries']) and not parker['integer_magic'] and not parker['distinct'], 'certificate invariant failed')
 require(constant['integer_magic'] and all(constant['integer_square_entries']) and constant['no_wrap'] and not constant['distinct'], 'certificate invariant failed')
 return dict(minimal_prime_no_wrap_residue_hostile=weak,lo_shu_hostile=lo,cain_mod29_hostile=cain,parker_integer_near_miss=parker,constant_positive_control_without_distinctness=constant)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--bound',type=int,default=100000);ap.add_argument('--independent-bound',type=int,default=1000);ap.add_argument('--output',type=Path,default=Path(__file__).with_name('certificates.json'));args=ap.parse_args()
 data=dict(schema='session55-parker-lifts-v1',controls=controls(),weak_prime_search=weak_prime_search(),integer_search=integer_search(args.bound,args.independent_bound))
 require(next(r['p'] for r in data['weak_prime_search'] if r['canonical_grids'])==59, 'certificate invariant failed')
 require(all(v=='rejected' for v in tamper_selftest(data).values()),'tamper self-test failed')
 if args.verify:verify_payload(data,json.loads(args.output.read_text()));print('VERIFIED existing certificate')
 else:args.output.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
 print(json.dumps(data['integer_search'],indent=2,sort_keys=True)); print('Checked' if args.verify else 'Wrote',args.output)
if __name__=='__main__':main()
