"""Exact continuous-simplex certificate, new orders 8 and 9."""
from pathlib import Path
import itertools,json,math,hashlib
ROOT=Path(__file__).resolve().parent
def need(ok,label):
 if not ok:raise RuntimeError(label)
def mul(a,b,n):
 out={}
 for e,c in a.items():
  for f,d in b.items():
   g=tuple(x+y for x,y in zip(e,f))
   if max(g)<=n:out[g]=out.get(g,0)+c*d
 return {g:c for g,c in out.items() if c}
E={(1,0,0):1,(0,1,0):1,(0,0,1):1,(1,1,0):1,(1,0,1):1,(0,1,1):1,(1,1,1):3}
K={(2,2,0):1,(2,0,2):1,(0,2,2):1,(2,1,0):1,(2,0,1):1,(1,2,0):1,(1,0,2):1,(0,2,1):1,(0,1,2):1,(2,1,1):-1,(1,2,1):-1,(1,1,2):-1,(1,1,1):-6}
def stirling_profile(r,n):
 s=[1]+[0]*n
 for _ in range(r):s=[0]+[k*s[k]+s[k-1] for k in range(1,n+1)]
 return [math.factorial(k)*s[k] for k in range(n+1)]
def mixed(num,h,qs,n):
 return sum(c*math.comb(k+h,h)*qs[0][k+e[0]]*qs[1][k+e[1]]*qs[2][k+e[2]] for e,c in num.items() for k in range(n-max(e)+1))
records=[]
for n in range(8,10):
 orders=list(range(1 if n%2 else 2,n+1,2));qs=[stirling_profile(r,n) for r in orders];d=len(qs);triples=list(itertools.combinations_with_replacement(range(d),3));power={(0,0,0):1};L=math.factorial(n)//(2**(n//2))-1;coords=[]
 for h in range(1,3*n+3):
  num=mul(K,power,n)
  if h>=2:
   out={mon:0 for mon in triples}
   for t in triples:
    mult=6//math.prod(math.factorial(t.count(i)) for i in set(t));value=h*mult*mixed(num,h,[qs[i] for i in t],n)
    need(value%3==0,'integral homogeneous coefficient');value//=3
    lower=[i for i in t if i!=d-1];p=len(lower)
    for fill in itertools.product(range(d),repeat=3-p):
     mon=tuple(sorted(lower+list(fill)));out[mon]+=value*L**p
   need(all(c<=0 for c in out.values()),f'nonpositive simplex coefficients n={n},h={h}')
   coords.append({'h':h,'coefficients':[[list(mon),c] for mon,c in sorted(out.items())]})
  power=mul(power,E,n)
 # Above the cover order, every numerator monomial exceeds the retained box.
 need(not mul(K,power,n),'vanishing tail numerator')
 records.append({'n':n,'generator_orders':orders,'lower_weight_budget':L,'coordinates':coords})
# A positive coefficient at order ten records only a certificate boundary.
n=10;h=2;orders=[2,4,6,8,10];qs=[stirling_profile(r,n) for r in orders];d=5;L=math.factorial(n)//32-1
num=mul(K,E,n);out={mon:0 for mon in itertools.combinations_with_replacement(range(d),3)}
for t in out:
 mult=6//math.prod(math.factorial(t.count(i)) for i in set(t));value=h*mult*mixed(num,h,[qs[i] for i in t],n)
 need(value%3==0,'integral boundary coefficient');value//=3
 lower=[i for i in t if i!=d-1];p=len(lower)
 for fill in itertools.product(range(d),repeat=3-p):out[tuple(sorted(lower+list(fill)))]+=value*L**p
boundary={'n':10,'h':2,'simplex_monomial':[0,0,2],'positive_coefficient':out[(0,0,2)],'interpretation':'Failure of a sufficient coefficient-sign certificate only; no actual-profile counterexample.'}
need(boundary['positive_coefficient']==20098555126088950347840000,'explicit order-ten certificate boundary')
print(json.dumps({'status':'PASS','scope':'New order8 and9 coefficient certificates under H<=n!/2^floor(n/2), yielding the stated fixed-ambient actual-tournament corollary together with the earlier1..7 certificate. No claim past9.','next_order_certificate_boundary':boundary,'coordinate_certificates':sum(len(r['coordinates']) for r in records),'integer_coefficient_certificates':sum(len(c['coefficients']) for r in records for c in r['coordinates']),'records':records,'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'PROOF.md',ROOT/'PROVENANCE.md',ROOT/'verify.py']}},sort_keys=True,indent=2))
