"""Independent exact convolution and direct geometric-sum row audit."""
import importlib.util,itertools,json,random,hashlib
from pathlib import Path
D=Path(__file__).resolve().parent; P=D.parent/'oddperfect'
s=importlib.util.spec_from_file_location('candidate',P/'verify.py');a=importlib.util.module_from_spec(s);s.loader.exec_module(a)
def need(v,m):
 if not v:raise RuntimeError(m)
def direct(S,q):
 out=[]
 for i,p in enumerate(S):
  rr=[]
  # The inherited primitive-divisor theorem bounds n by max(r-1).
  for n in range(2,max(S)):
   if not(n%4==2 if p==q else n%2==1):continue
   x=sum(p**k for k in range(n));two=0
   while x%2==0:x//=2;two+=1
   v=[]
   for r in S:
    k=0
    while x%r==0:x//=r;k+=1
    v.append(k)
   if x==1 and two==int(p==q):v[i]-=n-1;rr.append(tuple(v))
  out.append(rr)
 return out
cases=0
ps=[3,5,7,11,13,17,19,23]
for k in [2,3,4]:
 for S in itertools.combinations(ps,k):
  for q in S:
   if q%4==1:need(direct(S,q)==a.rows(S,q),'row catalogue');cases+=1
rng=random.Random(484);tests=0
for i in range(1000):
 d=rng.randrange(1,5);R=[[tuple(rng.randrange(-5,6)for _ in range(d))for _ in range(rng.randrange(1,5))]for _ in range(rng.randrange(1,5))]
 exact={(0,)*d:1}
 for rr in R:
  nxt={}
  for x,v in exact.items():
   for row in rr:
    y=tuple(b+c for b,c in zip(x,row));nxt[y]=nxt.get(y,0)+v
  exact=nxt
 m=a.moduli(R,d)
 need(all(all(-m[j]<x[j]<m[j]for j in range(d))for x in exact),'no-wrap range')
 need(a.count_mod(R,m)==exact.get((0,)*d,0),'integer convolution');tests+=1
need(a.count_mod([[(0,0)],[(0,0)]],[1,1])==1,'modulus one')
need(a.count_mod([[(-1,),(1,)],[(1,),(-1,)]],[1])==4,'small modulus accepts extra choices')
need(a.count_mod([[(-1,),(1,)],[(1,),(-1,)]],[3])==2,'exact multiplicity')
r={'status':'PASS','independent_direct_row_cases':cases,'independent_integer_convolutions':tests,'modulus_one_and_multiplicity_controls':True,'proof_sha256':hashlib.sha256((P/'PROOF.md').read_bytes()).hexdigest(),'scope':'Fixed prescribed support only; no global OPN theorem.'}
print(json.dumps(r,indent=2,sort_keys=True))
