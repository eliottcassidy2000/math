from pathlib import Path
import json, importlib.util, itertools
D=Path(__file__).parent; F=D/'frozen'
C=[tuple(map(int,s.split())) for s in (F/'CORE_STATES.txt').read_text().splitlines()]
spec=importlib.util.spec_from_file_location('inverse',F/'inverse_alpha.py'); inv=importlib.util.module_from_spec(spec);spec.loader.exec_module(inv)
def compile(c,p,r,b):
 h,q,A,G,T,S,E=C[c]
 # Prepend prefix vertices right to left, then append suffix left to right.
 for i in range(p.bit_length()):
  e=(p>>i)&1
  q,A,G,T,S,E=2*q,A+e*q,2*G+e*h,T+e*(G+S),e*h,E
 for i in reversed(range(b)):
  e=(r>>i)&1
  q,A,G,T,S,E=2*q,2*A+e*h,G+e*q,T+e*(A+E),S,e*h
 return 2*A+h-T,A,T
# Enumerate every grammar candidate with alpha <= target_cap-h, without alpha stripping.
cap=1000; represented={}; candidates=0
for c,(h,q,A,*_) in enumerate(C):
 for p in range((cap-h-A)//q+1):
  if p==0 and A==0:continue
  a0=q*p+A;b=0
  while (1<<b)*a0<=cap-h:
   bound=min((1<<b)-1,(cap-h-(1<<b)*a0)//h)
   for r in range(bound+1):
    n,alpha,tau=compile(c,p,r,b);candidates+=1
    if n<=cap:represented.setdefault(n,(c,p,r,b))
   b+=1
for n in range(1,cap+1):
 z=inv.find(n)
 assert (z['status']=='FOUND')==(n in represented),(n,z,represented.get(n))
 if z['status']=='FOUND':assert compile(z['core'],z['P'],z['R'],z['b'])==(n,z['alpha'],z['tau'])
# Literal recorded graph validation independent of producer graph function.
cores=[(8,1,3,6),(2,4,9,19,7),(2,8,3,21,7)]
large=json.loads((F/'LARGE_GRAPHS.json').read_text())
for z in large:
 c=z['core'];p=z['P'];r=z['R'];b=z['b'];a=p.bit_length();ci=0 if c<16 else 1+(c-16)//32;cut=c if c<16 else (c-16)%32;core=cores[ci];k=len(core);n=a+k+b+2;out=z['outgoing_masks'];assert len(out)==n
 selected=[bool(p>>(a-i-1)&1) for i in range(a)]+[bool(cut>>i&1) for i in range(k)]+[bool(r>>(b-i-1)&1) for i in range(b)]
 for i in range(n):
  assert out[i]>=0 and out[i]>>n==0 and not(out[i]>>i&1)
  for j in range(i+1,n):
   edge=bool(core[i-a]>>(j-a)&1) if a<=i<j<a+k else (selected[i] if j==n-1 and i<n-2 else False)
   # For x edges edge means x->i, unlike ordinary increasing edges.
   expected=(not edge) if j==n-1 and i<n-2 else (edge if a<=i<j<a+k else True)
   assert bool(out[i]>>j&1)==expected and bool(out[j]>>i&1)!=expected
 cyc=z['Hamilton_cycle'];assert sorted(cyc)==list(range(n));assert all(out[u]>>v&1 for u,v in zip(cyc,cyc[1:]+cyc[:1]))
 assert compile(c,p,r,b)[0]==z['N']
# Exhaust elementary f bound/equality and suffix constant/interval inequalities.
for p in range(1,100000):
 f=inv.f(p);assert 0<=f<=p-1;assert (f==p-1)==(p&(p+1)==0)
for c,(h,q,A,G,T,S,E) in enumerate(C):
 assert T<=A and S<=h and G<=q-h
 for b in range(8):
  for r in range(1<<b):
   const=compile(c,0,r,b)[0];rev=int(format(r,f'0{b}b')[::-1],2) if b else 0;K=q*((1<<(b+1))-rev)-G
   assert const>=h and K>=q*(1<<b)+h
print(json.dumps(dict(all_passed=True,inverse_targets=cap,direct_candidates=candidates,large_graphs=len(large),f_controls=99999),indent=2))
