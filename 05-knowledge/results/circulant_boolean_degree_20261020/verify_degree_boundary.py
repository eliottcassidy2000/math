#!/usr/bin/env python3
"""Self-contained exact controls for the Boolean-degree boundary.
All failures are explicit, so python -O leaves checking enabled.
"""
from collections import Counter
from itertools import combinations
import argparse,hashlib,json,math
from pathlib import Path
CHECKS=0

def check(ok,*context):
 global CHECKS
 CHECKS+=1
 if not ok:raise RuntimeError(context)

def shift(x,m,twisted=False):
 bits=[(x>>i)&1 for i in range(m)]
 out=[bits[-1]^int(twisted)]+bits[:-1]
 return sum(b<<i for i,b in enumerate(out))

def cycle(x,m,twisted=False):
 first=x;out=[]
 while True:
  out.append(x);x=shift(x,m,twisted)
  if x==first:return out
  check(len(out)<=2*m,'orbit failed to close',m,first)

def good(x,m):
 d=len(cycle(x,m));return sum((x>>i)&1 for i in range(d))%2==0

def mobius(table):
 coef=list(table);n=(len(coef)-1).bit_length()
 for i in range(n):
  for mask in range(1<<n):
   if mask&(1<<i):coef[mask]^=coef[mask^(1<<i)]
 return coef

def scalar(x,m,monomials):
 last=(x>>(m-1))&1
 z=[((x>>i)&1)^last for i in range(m-1)]
 return sum(math.prod(z[j] for j in mon) for mon in monomials)%2

def vector(x,m,monomials):
 # Direct inverse powers via 2m-1 forward twisted shifts, independently
 # of the search's inverse bit-shift implementation.
 answer=0
 for i in range(m):
  y=x
  for _ in range((2*m-i)%(2*m)):y=shift(y,m,True)
  answer|=scalar(y,m,monomials)<<i
 return answer

WITNESSES={
 4:[(1,),(0,1)],
 6:[(0,),(1,),(0,2),(0,3)],
 8:[(0,),(0,4),(1,5),(3,4),(4,5),(1,2,3),(3,5,6)],
}

def verify_witness(m,mons):
 values=[vector(x,m,mons) for x in range(1<<m)]
 target={x for x in range(1<<m) if good(x,m)}
 type_map={};fibers=Counter(values)
 for x,y in enumerate(values):
  check(good(y,m),'not target-valued',m,x,y)
  check(values[shift(x,m,True)]==shift(y,m),'equivariance',m,x)
  check(values[x^((1<<m)-1)]==y,'complement',m,x)
  a=min(cycle(x,m,True));b=min(cycle(y,m))
  check(a not in type_map or type_map[a]==b,'ill-defined type',m,x)
  type_map[a]=b
 check(set(values)==target,'labelled surjectivity',m)
 check(len(set(type_map.values()))==len(type_map),'type injection',m)
 check(set(type_map.values())=={min(cycle(y,m)) for y in target},'type surjectivity',m)
 full_anf=mobius([y&1 for y in values]);q_anf=mobius([values[z]&1 for z in range(1<<(m-1))])
 full_degree=max((s.bit_count() for s,a in enumerate(full_anf) if a),default=0)
 quotient_degree=max((s.bit_count() for s,a in enumerate(q_anf) if a),default=0)
 check(full_degree==quotient_degree==max(map(len,mons)),'wrong ANF degree',m)
 check({s for s,a in enumerate(q_anf) if a}=={sum(1<<i for i in mon) for mon in mons},'wrong ANF monomials',m)
 reps=sorted(type_map)
 return {'m':m,'prime':2*m+1 if all((2*m+1)%d for d in range(2,math.isqrt(2*m+1)+1)) else None,'degree':full_degree,'monomials':[list(mon) for mon in mons],'source_types':len(type_map),'target_labelled':len(target),'scalar_truth_weight':sum(y&1 for y in values),'quotient_Walsh':sum(1-2*(values[z]&1) for z in range(1<<(m-1))),'source_representatives':reps,'images':[values[x] for x in reps],'labelled_fiber_histogram':dict(sorted(Counter(fibers.values()).items()))}

def verify_walsh_divisibility():
 report=[]
 for r in range(1,5):
  for truth in range(1<<(1<<r)):
   values=[(truth>>x)&1 for x in range(1<<r)]
   coefficients=mobius(values)
   d=max((s.bit_count() for s,a in enumerate(coefficients) if a),default=0)
   if d:
    W=(1<<r)-2*truth.bit_count()
    check(W%(1<<((r+d-1)//d))==0,'Walsh divisibility',r,truth,d,W)
  report.append({'variables':r,'Boolean_functions':1<<(1<<r)})
 return report

def verify_dyadic_weight(m):
 seen=set();target_reps=[]
 for x in range(1<<m):
  if x not in seen:
   orb=cycle(x,m);seen.update(orb)
   if good(x,m):target_reps.append(min(orb))
 A=(1<<(m-1))//m
 check(len(target_reps)==A,'type count',m)
 check(sum(x.bit_count() for x in target_reps)==(m//2)*(A-1),'target weight sum',m)
 check(all(min(cycle(x^((1<<m)-1),m)) in target_reps for x in target_reps if x),'complement closure',m)
 if m<=8:
  source_reps={min(cycle(x,m,True)) for x in range(1<<m)}
  check(len(source_reps)==A,'source count',m)
  check(all(len(cycle(x,m,True))==2*m for x in source_reps),'source freedom',m)
 return {'m':m,'source_target_types':A,'sum_target_type_weights':sum(x.bit_count() for x in target_reps),'forced_quotient_Walsh':m,'lower_degree':(m-2+(m.bit_length()-1))//(m.bit_length()-1)}

def embed_source(x,e,q):
 return sum((((x>>(i%e))&1)^((i//e)&1))<<i for i in range(e*q))

def embed_target(x,e,q):return sum(((x>>(i%e))&1)<<i for i in range(e*q))

def verify_restrictions():
 rows=[]
 for e,q in [(4,3),(4,5),(8,3),(8,5)]:
  m=e*q
  source_types=set();embedded_types=set();small_target_types=set();large_target_types=set()
  for x in range(1<<e):
   y=embed_source(x,e,q)
   check(shift(y,m,True)==embed_source(shift(x,e,True),e,q),'source embedding equivariance',e,q,x)
   z=y
   for _ in range(2*e):z=shift(z,m,True)
   check(z==y,'source fixed subset',e,q,x)
   source_types.add(min(cycle(x,e,True)));embedded_types.add(min(cycle(y,m,True)))
   y=embed_target(x,e,q)
   check(shift(y,m)==embed_target(shift(x,e),e,q),'target embedding equivariance',e,q,x)
   check(good(x,e)==good(y,m),'primitive-root predicate',e,q,x)
   if good(x,e):small_target_types.add(min(cycle(x,e)));large_target_types.add(min(cycle(y,m)))
  check(len(source_types)==len(embedded_types)==len(small_target_types)==len(large_target_types),'small orbit equality',e,q)
  rows.append({'e':e,'q':q,'m':m,'small_types':len(source_types)})
 return rows

def prime(n):return n>=2 and all(n%d for d in range(2,math.isqrt(n)+1))
def examples():
 out=[]
 for p in [13,17,97,193,641,257,65537]:
  check(prime(p),'example not prime',p)
  m=(p-1)//2;e=m&-m;a=e.bit_length()-1
  out.append({'p':p,'m':m,'dyadic_part':e,'lower_degree':2 if e==2 and m>2 else (e-1+a-1)//a})
 return out

def main():
 global CHECKS
 ap=argparse.ArgumentParser();ap.add_argument('--output',default=None);args=ap.parse_args()
 result={'status':'FINITE-EXACT controls; all-order lower bound proved separately','witnesses':[verify_witness(m,mons) for m,mons in WITNESSES.items()],'walsh_exhaustive_controls':verify_walsh_divisibility(),'dyadic_weight_controls':[verify_dyadic_weight(m) for m in [4,8,16]],'restriction_controls':verify_restrictions(),'prime_examples':examples()}
 result['checks']=CHECKS;result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 s=json.dumps(result,indent=2,sort_keys=True)+'\n'
 if args.output:Path(args.output).write_text(s)
 print(s,end='')
if __name__=='__main__':main()
