#!/usr/bin/env python3
"""Exact collision-sieve certificate; stdlib, no producer imports."""
import argparse, hashlib, json, math, time, statistics
from fractions import Fraction
from pathlib import Path
P=(3,5,7,11,13,17,19,23,31,43,47,67)
HERE=Path(__file__).resolve().parent

def need(x,msg='invariant failed'):
 if not x: raise RuntimeError(msg)
def coll(u,v,p):return u*v*(u*u-v*v)*(u*u-4*v*v)*(4*u*u-v*v)%p==0
def points(u,v,p):return {(s*u+t*v)%p for s in (-1,0,1) for t in (-1,0,1)}
def masks():return {p:tuple(sum(1<<v for v in range(p) if coll(u,v,p)) for u in range(p)) for p in P}
M=masks()
def banks(B):
 out=[set() for _ in range(B+1)]
 for m in range(2,math.isqrt(B)+1):
  for n in range(1,m):
   h=m*m+n*n
   if h>B:break
   if math.gcd(m,n)!=1 or (m-n)%2==0:continue
   d=4*m*n*(m*m-n*n)
   for k in range(1,B//h+1):out[k*h].add(k*k*d)
 return out

def census(B):
 D=banks(B);digest=hashlib.sha256();pairs=[];buckets=[];root_count=0
 for e,A in enumerate(D[1:],1):
  digest.update((str(e)+':'+','.join(map(str,sorted(A)))+'\n').encode())
  for d in A:
   for z in (e*e-d,e*e+d):need(math.isqrt(z)**2==z);root_count+=1
  if len(A)<4:continue
  a=sorted(A);buckets.append((e,a,A))
  pairs.extend((e,x,y) for i,x in enumerate(a) for y in a[:i])
 stages=[];remaining=pairs;standalone={}
 for p in P:
  standalone[str(p)]=sum(not coll(x,y,p) for e,x,y in pairs)
  after=[z for z in remaining if coll(z[1],z[2],p)]
  stages.append(dict(prime=p,before=len(remaining),rejected=len(remaining)-len(after),remaining=len(after)))
  remaining=after
 def full(xs):
  distinct=diff=both=0;solutions=[]
  for e,x,y in xs:
   if len({x,y,x-y,x+y})!=4:continue
   distinct+=1
   if x-y not in D[e]:continue
   diff+=1
   if x+y not in D[e]:continue
   both+=1;solutions.append((e,x,y))
  return dict(distinct_magnitudes=distinct,difference_gap_present=diff,both_gaps_present=both,solutions=solutions)
 baseline=full(pairs);filtered=full(remaining);need(baseline['solutions']==filtered['solutions'])
 # Independent direct square bank control, including zero small root.
 for e in range(1,min(B,1000)+1):
  c=e*e;direct={c-a*a for a in range(e) if math.isqrt(2*c-a*a)**2==2*c-a*a};need(D[e]==direct)
 def direct_run():return [(e,x,y) for e,x,y in pairs if len({x,y,x-y,x+y})==4 and x-y in D[e] and x+y in D[e]]
 def poly_run():return [(e,x,y) for e,x,y in pairs if all(coll(x,y,p) for p in P) and len({x,y,x-y,x+y})==4 and x-y in D[e] and x+y in D[e]]
 def mask_run():return [(e,x,y) for e,x,y in pairs if all((M[p][x%p]>>(y%p))&1 for p in P) and len({x,y,x-y,x+y})==4 and x-y in D[e] and x+y in D[e]]
 def fast_gap_run():return [(e,x,y) for e,x,y in pairs if x-y in D[e] and x!=2*y and x+y in D[e]]
 def reduced_mask_run():return [(e,x,y) for e,x,y in pairs if all((M[p][x%p]>>(y%p))&1 for p in reversed(P[4:])) and x!=2*y and x-y in D[e] and x+y in D[e]]
 def conditional_run():
  # Cheap exact gap lookup first; modular rejection only at its survivors.
  return [(e,x,y) for e,x,y in pairs if x-y in D[e] and x!=2*y and all((M[p][x%p]>>(y%p))&1 for p in P) and x+y in D[e]]
 def signature_run():
  out=[]
  for e,a,A in buckets:
   sig={p:[v%p for v in a] for p in P}
   # Each prime maps the y residue to a bitset of positions in this bank.
   index={p:{} for p in P}
   for p in P:
    for j,r in enumerate(sig[p]):index[p][r]=index[p].get(r,0)|(1<<j)
   for i,x in enumerate(a):
    bits=(1<<i)-1
    for p in P:
     allowed=M[p][sig[p][i]];keep=0
     for r,b in index[p].items():
      if (allowed>>r)&1:keep|=b
     bits &= keep
     if not bits:break
    while bits:
     b=bits&-bits;j=b.bit_length()-1;bits-=b;y=a[j]
     if x!=2*y and x-y in A and x+y in A:out.append((e,x,y))
  return out
 timings={}
 for name,fn in [('full_gap_baseline',direct_run),('difference_first_gap_baseline',fast_gap_run),('reduced_descending_mask_then_gaps',reduced_mask_run),('polynomial_then_gaps',poly_run),('mask_then_gaps',mask_run),('difference_then_mask_then_sum',conditional_run),('bank_signature_bitsets_then_gaps',signature_run)]:
  values=[]
  for _ in range(3):
   t=time.perf_counter();answer=fn();values.append(time.perf_counter()-t);need(answer==baseline['solutions'])
  timings[name]=dict(seconds=values,median_seconds=statistics.median(values))
 survivors_digest=hashlib.sha256(''.join(f'{e},{x},{y}\n' for e,x,y in remaining).encode()).hexdigest()
 # Include actual roots and Euclidean carries for a sample surviving candidate; no residue-only upgrade.
 sample=[]
 for e,x,y in remaining[:3]:
  c=e*e;z=[c+x,c-x-y,c+y,c-x+y,c,c+x-y,c-y,c+x+y,c-x]
  roots=[math.isqrt(v) if v>=0 and math.isqrt(v)**2==v else None for v in z]
  sample.append(dict(e=e,u=x,v=y,gcd_uv=math.gcd(x,y),entries=z,exact_roots=roots,residue_and_carry_67=[list(divmod(v,67)) for v in z]))
 result=dict(bound=B,pairs=len(pairs),eligible_centers=len(buckets),bank_sha256=digest.hexdigest(),verified_actual_roots=root_count,standalone_rejections=standalone,sequential=stages,remaining=len(remaining),survivor_sha256=survivors_digest,baseline_full_gap=baseline,filtered_full_gap=filtered,survivor_samples=sample,independent_brute_bank_bound=min(B,1000))
 csv='e,u,v,gcd_uv,root_c_minus_u,root_c_plus_u,root_c_minus_v,root_c_plus_v\n'+''.join(','.join(map(str,(e,x,y,math.gcd(x,y),math.isqrt(e*e-x),math.isqrt(e*e+x),math.isqrt(e*e-y),math.isqrt(e*e+y))))+'\n' for e,x,y in remaining)
 return result,timings,csv

def controls():
 rows=[];aff=Fraction(1);prim=Fraction(1);conditioned=[]
 for p in P:
  n=0
  for u in range(p):
   for v in range(p):need(coll(u,v,p)==(len(points(u,v,p))<9));n+=coll(u,v,p)
  lines=4 if p==3 else 6 if p==5 else 8
  need(n==1+lines*(p-1));aff*=Fraction(n,p*p);prim*=Fraction(lines,p+1)
  Q={r*r%p for r in range(p)}
  for c in (0,1):
   G=[d for d in range(p) if (c-d)%p in Q and (c+d)%p in Q]
   rejects=sum(not coll(u,v,p) for u in G for v in G)
   conditioned.append(dict(p=p,normalized_center=c,gaps=G,ordered_pairs=len(G)**2,collision_rejected=rejects))
   if p<=11:need(rejects==0)
  rows.append(dict(p=p,lines=lines,allowed_affine=n,total_affine=p*p,allowed_nonzero=n-1,total_nonzero=p*p-1))
 # Exhaustive small CRT product: no statistical independence assumption.
 mod=3*5*7*11;local=P[:4]
 count=sum(all(coll(u,v,p) for p in local) for u in range(mod) for v in range(mod))
 expected=math.prod(1+(4 if p==3 else 6 if p==5 else 8)*(p-1) for p in local);need(count==expected)
 # Positive collision all squares, hostile non-Parker field distinct offsets.
 need(all(coll(0,0,p) for p in P));need(not coll(1,3,11));need(coll(11,33,11))
 return dict(local=rows,opposite_gap_conditioned=conditioned,affine_CRT_density=str(aff),primitive_projective_CRT_density=str(prim),CRT_modulus=math.prod(P),small_CRT=dict(modulus=mod,allowed=count),content_hostile=dict(u=11,v=33,p=11,passes_before_division=True,passes_after_division=False))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--bound',type=int,default=100000);args=ap.parse_args()
 result,timings,csv=census(args.bound);data=dict(schema='session55-collision-sieve-v1',controls=controls(),search=result)
 path=HERE/'certificate.json'
 if args.verify:need(json.loads(path.read_text())==data,'certificate mismatch');need((HERE/'survivors.csv').read_text()==csv,'survivors mismatch');print('VERIFIED exact certificate')
 else:(HERE/'survivors.csv').write_text(csv);path.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n');(HERE/'benchmark.json').write_text(json.dumps(timings,indent=2,sort_keys=True)+'\n')
 print(json.dumps(result,indent=2));print(json.dumps(timings,indent=2))
if __name__=='__main__':main()
