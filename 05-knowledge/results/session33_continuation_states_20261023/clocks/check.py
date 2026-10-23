from itertools import product
import json
from pathlib import Path
checks=0
def ck(x):
 global checks
 checks+=1
 if not x: raise RuntimeError('failed')
def rank(c,S): return sum(v for A,v in enumerate(c) if A&S)
def service(c,d):
 # Independent unit-owner augmenting matching, not the Hall expression.
 owners=[A for A,v in enumerate(c) for _ in range(v)]
 match=[None]*len(owners)
 def aug(j,seen):
  for i,A in enumerate(owners):
   if A>>j&1 and i not in seen:
    seen.add(i)
    if match[i] is None or aug(match[i],seen):
     match[i]=j; return True
  return False
 return sum(aug(j,set()) for j,v in enumerate(d) for _ in range(v))
n=3; full=(1<<n)-1
signatures=set(); cases=0
for bits in product(range(2),repeat=full):
 c=(0,)+bits
 r=[rank(c,S) for S in range(full+1)]
 signatures.add(tuple(r))
 f=[r[full]-r[full^T] for T in range(full+1)]
 for A in range(1,full+1):
  inv=sum((-1)**((A^T).bit_count())*f[T] for T in range(full+1) if T&A==T)
  ck(inv==c[A])
 for d in product(range(4),repeat=n):
  defect=max(sum(d[j] for j in range(n) if S>>j&1)-r[S] for S in range(full+1))
  served=service(c,d)
  ck(sum(d)-served==defect)
  # Neighborhood-closure compression preserves maximum deficiency.
  active=[A for A in range(1,full+1) if c[A]]
  def closure(S):
   neigh={A for A in active if A&S}
   return sum(1<<j for j in range(n) if all(A in neigh for A in active if A>>j&1))
  closed=[S for S in range(full+1) if closure(S)==S]
  ck(defect==max(sum(d[j] for j in range(n) if S>>j&1)-r[S] for S in closed))
  cases+=1
ck(len(signatures)==128)
# Equal singleton rank and total capacity, different pair response.
a=[0]*8;b=[0]*8
a[3]=a[4]=1;b[5]=b[2]=1
ck(all(rank(a,S)==rank(b,S) for S in (1,2,4,7)))
ck(service(a,(1,1,0))==1 and service(b,(1,1,0))==2)
# Artificial source splitting with identical eligibility is exactly neutral.
# Temporal traces differ, cut-wise profiles do not: persistent owner vs handoff.
traces1=[(1,1)]; traces2=[(1,0),(0,1)]
ck(all(sum(x[t] for x in traces1)==sum(x[t] for x in traces2) for t in (0,1)))
u=[min(S.bit_count(),2) for S in range(8)]
f=[u[7]-u[7^T] for T in range(8)]
coeff=[sum((-1)**((A^T).bit_count())*f[T] for T in range(8) if T&A==T) for A in range(8)]
ck(coeff[7]==-1 and all(coeff[A]==1 for A in (3,5,6)))
result={'status':'PASS','checks':checks,'mask_profiles':128,'demand_cases':cases,'universe':'3 obligation labels; nonempty mask capacities in {0,1}; each demand in {0,1,2,3}','hostile':{'bank_A':{'12':1,'3':1},'bank_B':{'13':1,'2':1},'demand':[1,1,0],'served':[1,2]},'claims_not_tested':['physical LRC realizability','native annulus extraction','uniform 153N','universal Collatz entry']}
print(json.dumps(result,indent=2))
Path(__file__).with_name('receipt.json').write_text(json.dumps(result,indent=2)+'\n')
