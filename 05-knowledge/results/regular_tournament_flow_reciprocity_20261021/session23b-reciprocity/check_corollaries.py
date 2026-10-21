import sympy as S
from pathlib import Path
import hashlib,json
q=S.symbols('q');checks=0
def need(v,label):
 global checks;checks+=1
 if not v:raise RuntimeError(label)
arcs=sorted((i,j) for i in range(7) for j in range(7) if (j-i)%7 in {1,2,4})
expected={2:[[1,0,0,0,-1,0,0],[0,0,1,1,-1,-1,0]],3:[[1,0,-2,1]]}
obsall={}
for mult,sign in [(2,1),(3,-1)]:
 seen=set();obs=[]
 for e in arcs:
  if e in seen:continue
  orb=[];f=e
  while f not in orb:
   orb.append(f);u,v=f;f=((mult*u)%7,(mult*v)%7) if sign==1 else ((mult*v)%7,(mult*u)%7)
  need(set(orb)<=set(arcs),'actual arcs');seen.update(orb);obs.append(orb)
 B=S.zeros(7,len(obs))
 for c,orb in enumerate(obs):
  for u,v in orb:B[u,c]-=1;B[v,c]+=1
 rr=B.rref()[0];rows=[list(rr.row(i)) for i in range(7) if any(rr.row(i))]
 need(rows==expected[mult],'independent incidence constraints')
 obsall[mult]={'orbits':obs,'rank':B.rank(),'dimension':len(obs)-B.rank()}
U3=(1+q**3+q**6+q**9+q**12)/((1-q**3)*(1-q**6)**2*(1-q**9)**2)
V6=(1+q**15)/((1-q**6)*(1-q**12)*(1-q**18))
U7=1/(1-q**7)**3
for F in (U3,V6,U7):need(S.cancel(F.subs(q,1/q)+q**21*F)==0,'negative reciprocity shift21')
x,y=S.symbols('x y')
# Derivative summation of (a+f+1) over x^a y^f.
G=1/((1-x)*(1-y));weighted=S.diff(G,x)*x+S.diff(G,y)*y+G
need(S.cancel(weighted-(1-x*y)/((1-x)**2*(1-y)**2))==0,'geometric derivative sum')
need(S.cancel(weighted.subs({x:q**9,y:q**6})/(1-q**3)**2-U3)==0,'U3 derivation')
need(S.cancel((1+q**15)/((1-q**18)*(1-q**12)*(1-q**6))-V6)==0,'two parity classes')
for r in range(1,31):
 d=r*(2*r-1);need((d+r)//2==r*r and (d-r*r)%2==0,'anti involution dimension')
 need(d-r*r==r*(r-1),'cumulative exponent')
R=(1+2*q**3+3*q**4+q**5+3*q**6+q**7+3*q**8+2*q**9+q**12)/((1-q**3)**3*(1-q**4)**2*(1-q**5))
f=(R+4/(1-q**5)**2)/5
k=(1+q**4+q**5+q**6+q**7+q**8+q**12)/((1-q**3)*(1-q**5)*(1-q**6)*(1-q**8))
a=S.limit((1-q)**6*f,q,1);b=S.limit((1-q)**4*k,q,1)
need(S.cancel((b/S.factorial(3))/(a/S.factorial(5)))==S.Rational(2100,17),'pointwise constant')
need(S.cancel((b/S.factorial(4))/(a/S.factorial(6)))==S.Rational(3150,17),'cumulative constant')
root=Path(__file__).resolve().parent
print(json.dumps({'status':'PASS','checks':checks,'independent_arc_orbits':obsall,'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'REGULAR_FLOW_RECIPROCITY.md',root/'COROLLARIES.md',Path(__file__)]}},indent=2))
