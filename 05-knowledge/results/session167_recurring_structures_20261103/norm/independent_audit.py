#!/usr/bin/env python3
"""Independent symbolic, orbit, descent, and brute Hamilton controls."""
import itertools, json
import sympy as S
s,u,d,x=S.symbols('s u d x'); f=x**3-x*x-x-1
N=d**3+2*d*d*s+d*d*u+2*d*s*s-2*d*s*u-d*u*u+2*s**3-2*s*u*u+u**3
z=(x*x-x)*s+x*u+d
C=lambda v:(v[1],sum(v),v[0])
def need(a,msg):
 if not a: raise RuntimeError(msg)
need(S.resultant(f,z,x).expand()==N,'resultant norm identity')
need(S.expand(N.subs({s:u,u:s+u+d,d:s}, simultaneous=True)-N)==0,'symbolic invariance')
need(S.rem(z.subs({s:u,u:s+u+d,d:s},simultaneous=True)-x*z,f,x)==0,'literal intertwiner')
norm=S.lambdify((s,u,d),N,'math')
roots=[]; equal_moves=0
for L in range(1,51):
 for a in range(L+1):
  for b in range(L-a+1):
   v=(a,b,L-a-b)
   need(norm(*v)>0,'positive norm on full test simplex')
   w=v
   while w[1]>=w[0]+w[2]:
    inv=(w[2],w[0],w[1]-w[0]-w[2])
    need(sum(inv)<=sum(w),'nonincreasing inverse sum')
    if sum(inv)==sum(w):
     equal_moves+=1
     need(w[0]==w[2]==0 and inv[1]<inv[0]+inv[2], 'equal sum reaches root')
    need(norm(*w)==norm(*inv),'descent invariant')
    w=inv
   need(norm(*w)*1792>sum(w)**3,'strict cubic cone bound sampled')
   if norm(*v)==2: need(w==(1,0,0),'all test norm-two points reduce to unique root')
   if L<=15 and b<a+v[2] and norm(*v)==2: roots.append(v)
need(roots==[(1,0,0)],'independent finite root enumeration')
# Literal enumeration of Hamilton permutations, not subset dynamic programming.
state=(0,1,1)
for n in range(2,9):
 def arc(a,b): return b==a+1 or a>=b+2
 paths=[p for p in itertools.permutations(range(n)) if all(arc(a,b) for a,b in zip(p,p[1:]))]
 starts=sum(p[0]==n-1 for p in paths)
 smaller=sum(p[0]==n-2 and all(arc(a,b) for a,b in zip(p,p[1:])) for p in itertools.permutations(range(n-1)))
 need((starts,len(paths)-starts,smaller)==state,'brute Hamilton endpoint state')
 state=C(state)
# Direct orbit enumeration and scalar-period checks, independently of matrix powers.
period_results={}
for p,t in [(3,13),(5,31),(47,46)]:
 unseen={(a,b,c) for a in range(p) for b in range(p) for c in range(p) if norm(a,b,c)%p==2}
 lengths=[]
 while unseen:
  root=next(iter(unseen)); v=root; length=0
  while True:
   unseen.remove(v); length+=1; v=tuple(x%p for x in C(v))
   if v==root: break
  lengths.append(length)
 need(set(lengths)=={t},'all modular shell orbit lengths')
 for k in (1,2,3):
  m=p**k; initial=(1,2,0); v=initial; counts=[]
  while True:
   counts.append((v[0]+v[1])%m); v=tuple(x%m for x in C(v))
   if v==initial: break
  period=t*p**(k-1)
  need(len(counts)==period,'direct state period')
  for q in S.factorint(period):
   candidate=period//q
   need(any(counts[i]!=counts[(i+candidate)%period] for i in range(period)), 'direct scalar period minimality')
 period_results[p]={'shell_orbits':len(lengths),'base_period':t,'direct_prime_power_exponents':[1,2,3]}
# Explicit polynomial factorizations and higher-degree identities.
split_difference=S.Poly(S.expand((x-5)*(x-17)*(x-26)-f),x)
need(all(int(c)%47==0 for c in split_difference.all_coeffs()),'47 split factorization')
for k,pp,ee,inv,v1,v2 in [
 (4,1+x+x*x,1+x,2-x+2*x*x-x**3,[1,1,2,2],[1,2,2,1]),
 (5,1+x,1+x*x+x**3,1+3*x+4*x*x+x**3-2*x**4,[1,1,1,1,2],[1,1,1,2,1])]:
 fk=x**k-sum(x**j for j in range(k))
 poly=lambda vs:sum(c*x**i for i,c in enumerate(vs))
 need(S.rem(x**(k-2 if k==4 else 4)*pp-poly(v1),fk,x)==0,'first factor identity')
 need(S.rem(ee*pp-poly(v2),fk,x)==0,'second factor identity')
 need(S.resultant(fk,pp,x)==(3 if k==4 else 2),'p norm')
 need(S.resultant(fk,ee,x)==1,'epsilon norm')
 need(S.rem(ee*inv-1,fk,x)==0,'inverse identity')
q=x**4-x**3-x**2-1
need(q.subs(x,S.Rational(19,10))>0,'quintic upper inequality start')
print(json.dumps({'status':'PASS','symbolic_resultant_and_invariance':True,
 'simplex_total_sum_limit':50,'equal_sum_events_checked':equal_moves,
 'complete_norm_two_roots':roots,'brute_Hamilton_n_max':8,
 'direct_modular_orbits':period_results,'higher_degree_factor_and_inverse_certificates':True},indent=2))
