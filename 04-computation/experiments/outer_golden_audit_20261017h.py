#!/usr/bin/env python3
"""Separate reconstruction of compressed witnesses and finite clock controls."""
import argparse,json,math
from pathlib import Path
from fractions import Fraction as Q
from collatz_pair_mass_20261013b import step,parents

def need(ok,m):
 if not ok:raise RuntimeError(m)
def T(n):return (3*n+1)//2 if n%2 else n//2

def graph_audit():
 folder=Path('04-computation/experiments/outer_golden_20261017h')
 old=json.loads(Path('04-computation/experiments/root_positive_transport_20261017g/golden_shared.json').read_text())
 result=[]
 for name in ('golden','rotate','conjugate','signature'):
  r=json.loads((folder/(name+'.json')).read_text());n=r['vertices'];edges=r['edges'];checks=0
  need(len(edges)==len(set(map(tuple,edges))) and all(0<=a<b<n for a,b in edges),'simple graph')
  need(r['patterns']==old['patterns'],'all canonical port patterns retained')
  recipe=r.get('witness_recipe')
  for i,pat in enumerate(old['patterns']):
   c0=list(map(int,old['witnesses'][i]))
   if recipe is None:col=list(map(int,r['witnesses'][i]))
   elif recipe['kind']=='cyclic_colour_shift':col=c0+[(c+recipe['shift'])%6 for c in c0]
   elif recipe['kind']=='fold_to_inherited_graph':col=[c0[j] for j in recipe['fold']]
   elif recipe['kind']=='marked_cartesian_cycle':
    col=c0+[-1]*(n-len(c0))
    for p,ids in enumerate(r['attachments']):
     for j,a in enumerate(ids):col[a]=(c0[p]+recipe['cycle_colours'][j])%6
   else:raise RuntimeError('unknown witness representation')
   need(len(col)==n and col[:8]==pat and all(0<=c<6 for c in col),'literal source and colour range')
   for a,b in edges:need(col[a]!=col[b],'edge verified in reconstructed witness')
   checks+=len(edges)
  need(all(r['four_colour'][a]!=r['four_colour'][b] for a,b in edges),'four colour witness')
  if name=='signature':
   nb=[set() for _ in range(n)]
   for a,b in edges:nb[a].add(b);nb[b].add(a)
   live=set(range(n));peeled=[]
   while True:
    batch=[a for a in sorted(live) if a>=len(c0) and len(nb[a]&live)<=5]
    if not batch:break
    for a in batch:live.remove(a);peeled.append(a)
   need(len(live)==len(c0),'new exterior has empty relative six-core')
  result.append(dict(name=name,vertices=n,edges=len(edges),boundary_edge_checks=checks,patterns=570))
 return result

def triple_controls():
 count=0
 for N in range(1,258):
  for d in range(11):
   lo=Q(2*N+1,2**d);hi=Q(4*N+1,2**d);q=3**d
   hs=[h for h in range(math.floor(lo)+1,math.floor(hi)+1) if (q*h-1)%2 and (q*h-1)%8!=5]
   covered=set();triple_sum=0
   g=lambda h:((q*h-1)%7)%3-1
   for h0 in hs:
    if h0%3 or not (lo<h0-8 and h0+8<=hi):continue
    triple=(h0-8,h0,h0+8)
    need(all(h in hs for h in triple) and not(covered&set(triple)),'disjoint complete triples')
    covered.update(triple);triple_sum+=2*g(h0)-g(h0-8)-g(h0+8)
   residual=[h for h in hs if h not in covered]
   weight=lambda h:2 if h%3==0 else -1
   need(sum(abs(weight(h)) for h in residual)<=18,'uniform cropped-triple debt')
   need(abs(sum(weight(h)**2 for h in hs)-Q(3*N,2**(d+1)))<=18,'outer squared energy with cropped endpoints')
   need(sum(weight(h)*g(h) for h in hs)==triple_sum+sum(weight(h)*g(h) for h in residual),'exact signed decomposition for arbitrary observer')
   count+=1
 # Actual conditional translation cylinders after three steps.
 K=8;back={(0,0):1}
 for _ in range(K):back=step(back,parents)
 pairs=0
 for d in range(8):
  for r in (1,3,7):
   p=next(h for h in range(24) if h%3==0 and (3**d*h-1)%8==r)
   for direction in (-1,1):
    hits=0;e=3**(d+(3 if r==7 else 2))
    for z in range(2**K):
     h=24*(2**40+z)+p;a=3**d*h-1;b=a+direction*8*3**d;oa=ob=0
     for unused in range(3+K):oa+=a%2;ob+=b%2;a=T(a);b=T(b)
     hits+=a==b and oa==ob
    need(hits==back.get((0,e),0),'conditioned cylinder mass equals named translation law')
    pairs+=2**K
 return dict(cropped_interval_cases=count,conditioned_literal_pairs=pairs)

def golden_tower():
 out=[]
 for k in range(1,6):
  mod=3**k;length=8*3**(k-1);seen=set();norms=set();cycles=0
  for a in range(mod):
   for b in range(mod):
    if (a%3,b%3)==(0,0) or (a,b) in seen:continue
    x,y=a,b;cyc=[]
    while (x,y) not in seen:
     seen.add((x,y));cyc.append((x,y));x,y=(x+y)%mod,x
    values={(x*x-x*y-y*y)%mod for x,y in cyc}
    need(len(cyc)==length and len(values)==2,'full primitive orbit length and two opposite norms')
    label=tuple(sorted(values));need(label not in norms,'distinct cycles have distinct norm-pair label');norms.add(label);cycles+=1
  x,y=1,0;powers=[]
  for _ in range(length):powers.append((x,y));x,y=(x+y)%mod,x
  scalars=[];a=1
  for _ in range(2*3**(k-1)):scalars.append(a);a=2*a%mod
  joint={(a*x%mod,a*y%mod) for a in scalars for x,y in powers}
  need(joint==seen and cycles==3**(k-1),'complementary generators cover all primitive pairs')
  out.append(dict(k=k,modulus=mod,primitive_pairs=len(seen),Q_orbits=cycles,Q_order=length,scalar_order=len(scalars)))
 return out

def energy_scope_controls():
 # An abstract non-ROOT component is a type-level hostile, not a claimed
 # positive-integer Collatz counterexample. ROOT projection cannot be unitary.
 c={"root":3,"unrooted":5};f={"root":1,"unrooted":0}
 total=sum(v*v for v in c.values())
 projected=sum((c[u]*f[u])**2 for u in c)
 need((total,projected)==(34,9),'ROOT projection retains an unknown energy complement')
 mod=9;A=3;a=1;j=1
 primitive=lambda r:((3*r+1)*pow(2,-a,mod)-(4**j-1)//3)*pow(4,-j,mod)%mod
 ordinary=lambda r:(3*r+1)*pow(2,-A,mod)%mod
 need((primitive(0),ordinary(0))==(1,8),'primitive affine phase is not ordinary Syracuse phase')
 # Explicit physical clock versus the canonical stripped exponent.
 # 11 -> 17 -> 13 = B(3), stripping the terminal braid costs two in the
 # sum of homogeneous exponents, but is not an extra elapsed clock step.
 x=11;clock=0;exponents=[];braid=0;seed=11
 for _ in range(2):
  val=3*seed+1;a=(val&-val).bit_length()-1;t=val>>a;j=0
  while t%8==5:t=(t-1)//4;j+=1
  exponents.append(a+2*j)
  for unused in range(a+2*braid):x=T(x);clock+=1
  seed=t;braid=j
 need((seed,braid,x,clock,sum(exponents))==(3,1,13,3,5),'terminal braid corrects the integer clock')
 # Independent scalar source compiler and literal signed push; compare
 # squared vector difference, avoiding floating square-root comparisons.
 def row(N,rounded=False):
  out={}
  for seed in range(1,4*N+1,2):
   if seed%8==5:continue
   m=seed;j=0;co=0
   while (4**j*seed if rounded else m)<=4*N:
    y=4**j*seed if rounded else m;r=(seed+j)%3
    co+=([0,-3,3] if y<=2*N else [-1,-1,2])[r]
    j+=1;m=4*m+1
   if co:out[seed]=co
  return out
 def push(c):
  out={}
  for seed,co in c.items():
   a=2 if seed%8==1 else 1;t=(3*seed+1)//2**a
   while t%8==5:t=(t-1)//4
   out[t]=out.get(t,0)+(-1)**a*co
  return {t:co for t,co in out.items() if co}
 cases=0
 for N in list(range(1,65))+[127,128,129,1023,1024,1025]:
  exact=row(N);rounded=row(N,True)
  bound=4*(((4*N).bit_length()-1)//2)
  error=exact.get(1,0)-rounded.get(1,0)
  for depth in range(41):
   keys=exact.keys()|rounded.keys()
   diff2=sum((exact.get(t,0)-rounded.get(t,0))**2 for t in keys)
   need(diff2<=bound**2,'source-rounding energy radius uniform in exact transport depth')
   if N&(N-1)==0:
    E=sum(v*v for v in exact.values());H=sum(v*v for v in rounded.values())
    need(E==H+2*error*rounded.get(1,0)+error*error,'dyadic ROOT-atom energy correction')
   exact=push(exact);rounded=push(rounded);cases+=1
 return dict(rounding_energy_profiles=cases,abstract_energy=total,ROOT_projection=projected,primitive_phase_mod9=1,ordinary_phase_mod9=8,source=11,terminal_integer=x,terminal_seed=seed,terminal_braid=braid,actual_clock=clock,homogeneous_clock=sum(exponents))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args()
 r=dict(status='FINITE-EXACT separate formula and witness replay; self-audit',graphs=graph_audit(),triples=triple_controls(),golden_tower=golden_tower(),energy_scope=energy_scope_controls())
 Path(a.output).write_text(json.dumps(r,sort_keys=True,separators=(',',':'))+'\n');print(json.dumps(r,sort_keys=True))
