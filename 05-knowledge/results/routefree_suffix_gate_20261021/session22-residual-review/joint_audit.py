"""Independent exact joint-gate extension, separate from pinned core audit."""
import json,itertools,random
from pathlib import Path
from collections import Counter
from audit import A,q,b,require

def add(v,w):return tuple(v[k]+w[k] for k in (0,1))
def signed(vs):
 return [tuple(sum(e*v[k] for e,v in zip(es,vs)) for k in (0,1)) for es in itertools.product((-1,1),repeat=len(vs))]
def corners(vs):
 c=Counter(vs);total=tuple(sum(v[k] for v in vs) for k in (0,1))
 return [tuple(2*sum(c[A[j%12]]*A[j%12][k] for j in range(i,i+6))-total[k] for k in (0,1)) for i in range(12)]
def joint(old,new):return max(q(n)+b(o,n) for o in corners(old) for n in corners(new))
def run():
 rng=random.Random(14422021);tests=0
 for no in range(6):
  for nn in range(6):
   for _ in range(4):
    old=[rng.choice(A+((0,0),)) for j in range(no)];new=[rng.choice(A+((0,0),)) for j in range(nn)]
    G=joint(old,new)
    require(G==max(q(n)+b(o,n) for o in signed(old) for n in signed(new)))
    Do=max(map(q,signed(old)));Dn=max(map(q,signed(new)))
    require(G>=Dn and (G-Dn)**2<=4*Do*Dn)
    tests+=1
 old=[(1,0),(1,2)];new=[(-1,1),(0,1)];G=joint(old,new);mx=-10**9
 for choices in itertools.product(range(7),repeat=4):
  oo=[(0,0)]*3;nn=[(0,0)]*3
  for i,(z,v) in enumerate(zip(choices,old+new)):
   if z:
    u=(z-1)//2;sg=1 if z%2 else -1;t=oo if i<2 else nn
    t[u]=add(t[u],tuple(sg*x for x in v))
  value=sum(q(n)+b(o,n) for o,n in zip(oo,nn))
  require(value<=G);mx=max(mx,value)
 return {'brute_individual_sign_comparisons':tests,'partial_assignment_cases':2401,'joint_G':G,'actual_quadratic_mixed_max':mx,'scope':'Algebraic independent test of exact joint extension; no new native census'}
if __name__=='__main__':
 blob=json.dumps(run(),sort_keys=True,indent=2)+'\n'
 Path(__file__).with_name('joint-receipt.json').write_text(blob);print(blob)
