from math import comb
import json
from sympy import Matrix,ZZ
from sympy.matrices.normalforms import smith_normal_form

def vp(a,p):
 a=abs(int(a));v=0
 if not a:return 10**9
 while a%p==0:v+=1;a//=p
 return v
checks=[]
for p in [2,3,5,7]:
 for e in [1,2]:
  for m in range(1,7):
   for hs in [(0,1),(0,1,2),(0,1,p),(0,2,3)]:
    n=len(hs);xs=[p**e*h for h in hs]
    M=Matrix([[comb(c,r)*x**(c-r) if c>=r else 0 for c in range(n*m)] for x in xs for r in range(m)])
    S=smith_normal_form(M,domain=ZZ)
    a=sorted(vp(S[i,i],p) for i in range(n*m))
    pos=[v for v in a if v]
    first=e+vp(m,p);mult=2 if p==2 and e==1 and m%2==0 else 1
    if a.count(0)!=m or min(pos)!=first or pos.count(first)!=mult:raise Exception((p,e,m,hs,a,first,mult))
    checks.append(dict(p=p,e=e,m=m,h=hs,partition=a))
print(json.dumps({'full_smith_cases':len(checks),'all_pass':True,'examples':[x for x in checks if x['p']==2 and x['e']==1 and x['m']==4]},indent=2))
