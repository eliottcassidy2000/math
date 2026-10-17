#!/usr/bin/env python3
"""Exact port conjugates, trace connectors, and prime-19 component charges."""
import argparse,hashlib,json
from collections import defaultdict,Counter
from fractions import Fraction as F
from pathlib import Path
import sympy as sp
import cubic_hub_coupling_20261017a as H
import charge_ports_eleven_20261011 as G
checks=0
def need(ok,label):
 global checks
 checks+=1
 if not ok:raise RuntimeError(label)

def charge(v):
 z=H.mul(H.j,tuple(v[8:16]));return tuple(int(v[k]-v[k+16]+8*z[k])%19 for k in range(8))

def controls():
 hub=H.hub();d=hub[4];gamma=H.scale(H.add(H.scale(H.I,5),H.j),F(1,6))
 need(H.norm(d)==H.I and H.norm(gamma)==H.I,'half and j connectors are unit directions')
 need(H.sub(H.add(d,H.conj(d)),H.I)==H.scale(H.I,F(1,2)),'explicit half connector')
 need(H.sub(H.scale(gamma,6),H.scale(H.I,5))==H.j,'explicit j connector')
 need(H.sub(H.scale(H.I,2),H.add(gamma,H.conj(gamma)))==H.scale(H.I,F(1,3)),'explicit thirds connector')
 w=sp.symbols('w');f=w**3+2*w*w-4*w-4;traces=[]
 for r in (0,1):
  q=(w+r)**2+1;N=int(sp.resultant(f,q,w));p=sp.rem(2*N-4*N*sp.invert(q,f),f,w);traces.append(p)
 need(sp.expand(5*(traces[1]+64)-6*(traces[0]+26))==64*w,'dyadic trace connector recovers w without dividing by 11')
 alphas=[]
 for r in (1,2,3,4):
  e=next(e for e in range(55) if e%5==r and e%11==1)
  def autom(x):
   z=G.O
   for i,c in enumerate(x):
    if c:z=G.add(z,G.sc(G.powers[(i*e)%55],c))
   return z
  alphas.append([autom(p) for p in G.ports])
 need(all(G.mul(G.add(p[2],G.sc(G.I,-1)),G.conj(G.add(p[2],G.sc(G.I,-1))))==G.I for p in alphas),'four port-2 conjugates on unit circle centred at 1')
 need(len({p[2] for p in alphas})==4 and len({p[4] for p in alphas})==4,'four distinct circle conjugates')
 need(all(G.mul(p[4],G.conj(p[4]))==G.I for p in alphas),'port 4 conjugates centred at 0')
 need(all(G.add(p[3],G.conj(p[3]))==G.I for p in alphas) and len({p[3] for p in alphas})==4,'four distinct collinear port-3 conjugates')
 need({p[5] for p in alphas}=={G.ports[5],G.ports[7]},'rhombus conjugates exchanged')
 doubles=[]
 for sign in (-1,1):
  t=H.scale(H.add(H.scale(H.I,7),H.scale(H.basis(5),sign*3)),F(1,8))
  need(H.norm(H.sub(t,H.I))==H.I and H.norm(H.sub(t,H.scale(H.I,F(3,4))))==H.I,'two exact double-port points')
  need(H.norm(t)!=H.I,'double-port point does not see port 0')
  doubles.append([str(c) for c in t])
 zero=(0,)*8;u3=(H.I,H.conj(H.a),H.scale(H.a,-1));scaled=tuple(int(76*c) for block in u3 for c in block)
 need(charge(scaled)==zero,'mixed cubic unit has zero fixed-prime charge')
 pa=tuple(int(76*c) for c in H.a)+(0,)*16
 need(charge(pa)!=zero,'pure coefficient a cannot be unit-reachable')
 return dict(trace_polynomials=[str(t) for t in traces],double_port_points=doubles,
             nonzero_coefficient_charge=charge(pa),cubic_unit_charge=charge(scaled))

def golden_generator():
 u,j,t,v=sp.symbols('u j t v');N=2*j*u+t*(1-u);D=2*j+t*(1-u);value=N/D
 bar=value.subs({u:1/u,j:-j},simultaneous=True)
 need(sp.cancel(value*bar-1)==0,'golden-cubic Cayley generator has exact norm one')
 defect=sp.expand(N*D.subs(t,1-t)-N.subs(t,1-t)*D)
 need(sp.expand(defect-2*j*(2*t-1)*(u-1)**2)==0,'golden conjugate difference is nonzero in the specified field')
 inverse=(2*j*v+t*(v-1))/(2*j+t*(v-1))
 need(sp.cancel(inverse.subs(v,value)-u)==0,'mixed generator retains the full cubic coordinate')
 x=sp.symbols('s');f=(x-t)**3+2*(x-t)**2-4*(x-t)-4
 rem=sp.rem(sp.expand(f),t*t-t-1,t);aa=sp.expand(rem).coeff(t,0);bb=sp.expand(rem).coeff(t,1)
 h=sp.resultant(f,t*t-t-1,t)
 need(sp.gcd(aa,bb)==1,'sum coordinate recovers golden coordinate without exceptional common root')
 need(sp.expand(h-(x**6+x**5-14*x**4-11*x**3+30*x*x+13*x-19))==0,'exact primitive sextic with norm minus nineteen')
 need(sp.rem((4*t-3)*(4*(1-t)-3),t*t-t-1,t)==-19,'golden port difference has norm minus nineteen')
 return dict(formula=str(value),golden_conjugate_cross_numerator=str(sp.factor(defect)),inverse=str(inverse),
             primitive_sum_polynomial=str(h),golden_decoder=str(-aa/bb),
             status='PROVED norm-one primitive generator of M(sqrt5) over B; new cross-edge forcing OPEN')

def golden_probe():
 phi=G.add(G.I,G.add(G.a,G.conj(G.a)));psi=G.add(G.I,G.sc(phi,-1));points=G.ports+[phi,psi]
 E,surv,den=G.graph(points);new=[e for e in E if e[1]>=8]
 need(new==[(2,8),(4,9)],'golden quadratic escape has exactly two pendant port edges')
 beta=(0,1,2,0,2,1,2,1)
 need(all(beta[i]!=beta[j] for i,j in E if j<8),'declared hostile boundary is proper')
 for r in (9,10):
  need(r*r%19==5,'both real golden residue branches')
  a=(1+r)*10%19;b=(1-a)%19
  need(15 in (a,b),'one golden contact aliases port 6 modulo nineteen')
  need(beta[6]==beta[2]==beta[4],'both nineteen branches obstruct this residue-only boundary lift')
 return dict(coordinates=[[str(c) for c in phi],[str(c) for c in psi]],edges=E,
             boundary_hostile=list(beta),golden_residue_pairs=[[5,15],[15,5]],
             conclusion='The two pendants preserve all 570 boundary patterns, but no residue-only nineteen colouring of the expanded carrier realizes the declared pattern.')

def residue_boundary():
 c=(4,3,2,4,3,2,1,0,4,1,0,4,3,2,1,0,2,1,0)
 points=[(x,y) for x in range(19) for y in range(19)];index={p:i for i,p in enumerate(points)}
 directions=[p for p in points if (p[0]*p[0]+p[1]*p[1])%19==1]
 E=[(index[(x,y)],index[((x+a)%19,(y+b)%19)]) for x,y in points for a,b in directions]
 base=[c[(x+y)%19] for x,y in points];ports=[index[(0,0)],index[(1,0)],index[(15,0)]]
 need(len(directions)==20 and all(base[i]!=base[j] for i,j in E),'inherited nineteen-residue five-colouring on every directed edge')
 maps=[]
 for which in ('distinct','0=6','1=6'):
  col=base.copy()
  if which!='distinct':
   col[ports[2]]=5;col[ports[0 if which=='0=6' else 1]]=5
  need(all(col[i]!=col[j] for i,j in E),'six-colour repaired residue map')
  p=[col[i] for i in ports]
  need(p[0]!=p[1] and ((p[0]==p[2])==(which=='0=6')) and ((p[1]==p[2])==(which=='1=6')),'exact three-port relation pattern')
  maps.append(dict(pattern=which,colouring=col,port_colours=p))
 return dict(points=361,directed_edges=len(E),source_residues=ports,colour_maps=maps,
             conclusion='Every proper six-colour pattern on ports 0,1,6 lifts to the entire cubic carrier; external old ports have no new incidences.')

def hensel_colouring(P,E,sig):
 # One genuine fixed-prime branch, to precision 19^2. Charge is its
 # integrality test; the next digit gives the actual residue of x/76.
 modulus=361;r=27;s=159;w=99
 need((r*r-7)%modulus==0 and (s*s-11)%modulus==0 and (w**3+2*w*w-4*w-4)%modulus==0,'exact Hensel branches to nineteen squared')
 def mul(a,b):return ((a[0]*b[0]-a[1]*b[1])%modulus,(a[0]*b[1]+a[1]*b[0])%modulus)
 def inv(a):
  n=pow((a[0]*a[0]+a[1]*a[1])%modulus,-1,modulus);return (a[0]*n%modulus,-a[1]*n%modulus)
 u=mul((w+3,-s),inv((w+3,s)));powers=[(1,0),u,mul(u,u)];ev=[]
 for uk in powers:
  for j in range(8):
   c=(r if j&1 else 1)*(s if j&2 else 1);ev.append(mul(uk,(0,c) if j&4 else (c,0)))
 nums={i:tuple(sum(v*e[k] for v,e in zip(p,ev))%modulus for k in (0,1)) for i,p in P.items()};refs={};res={}
 for i in sorted(P):
  z=nums[i];ref=refs.setdefault(sig[i],z);delta=tuple((a-b)%modulus for a,b in zip(z,ref))
  need(all(a%19==0 for a in delta),'equal full charges give integral differences at chosen branch')
  res[i]=tuple((a//19)*5%19 for a in delta)
 need([res[i] for i in (0,1,6)]==[(0,0),(1,0),(15,0)],'actual marked source residues after exact division by seventy-six')
 need(all(((res[i][0]-res[j][0])**2+(res[i][1]-res[j][1])**2)%19==1 for i,j in E if i in P and j in P),'every exact unit edge maps to a nonzero residue unit')
 c=(4,3,2,4,3,2,1,0,4,1,0,4,3,2,1,0,2,1,0);col={i:c[sum(z)%19] for i,z in res.items()}
 need(all(col[i]!=col[j] for i,j in E if i in P and j in P),'independent five-colouring of the whole bank by Hensel reduction')
 pe=[(i,j) for i,j in E if j<8];patterns=[]
 def visit(t):
  n=len(t)
  if n==8:patterns.append(tuple(t));return
  for c0 in range(min(6,max(t,default=-1)+2)):
   if all(t[i]!=c0 for i,j in pe if j==n):visit(t+[c0])
 visit([]);zero=(0,)*8;root=[i for i in P if sig[i]==zero];other=[i for i in P if sig[i]!=zero]
 for beta in patterns:
  merged=None
  if beta[0]==beta[6]:merged={(0,0),(15,0)}
  elif beta[1]==beta[6]:merged={(1,0),(15,0)}
  def template(i):return 5 if merged is not None and res[i] in merged else col[i]
  rename={template(i):beta[i] for i in (0,1,6)}
  for a,b in zip((a for a in range(6) if a not in rename),(b for b in range(6) if b not in rename.values())):rename[a]=b
  values={i:rename[template(i)] for i in root};values.update({i:col[i] for i in other});values.update(enumerate(beta))
  need(all(values[i]!=values[j] for i,j in E),'every boundary explicitly extends by the independent Hensel colour maps')
 return dict(sqrt7=r,sqrt11=s,w=w,u=u,modulus=modulus,source_residues=[res[i] for i in (0,1,6)],
             residue_coordinates={str(i):list(z) for i,z in res.items()},carrier_colours={str(i):a for i,a in col.items()},
             boundary_patterns_replayed=len(patterns))

def run(path):
 raw=Path(path).read_bytes();g=json.loads(raw);P={int(i):v for i,v in g['coordinates'].items()};E=g['edges']
 need(g['denominator']==76,'declared denominator bank')
 sig={i:charge(v) for i,v in P.items()};zero=(0,)*8
 need(all(sig[i]==sig[j] for i,j in E if i in P and j in P),'every exact induced unit edge preserves charge')
 adj=defaultdict(list)
 for i,j in E:adj[i].append(j);adj[j].append(i)
 seen=set();sizes=[];root=set()
 for i in range(g['vertices']):
  if i in seen:continue
  stack=[i];seen.add(i);c=set()
  while stack:
   a=stack.pop();c.add(a)
   for b in adj[a]:
    if b not in seen:seen.add(b);stack.append(b)
  sizes.append(len(c))
  if 0 in c:root=c
 need(all(sig[i]==zero for i in root if i in sig),'all source-connected carrier vertices have zero charge')
 missing=sorted(i for i in P if sig[i]==zero and i not in root)
 return dict(source_sha256=hashlib.sha256(raw).hexdigest(),hensel_colouring=hensel_colouring(P,E,sig),charge_fibres=len(set(sig.values())),zero_charge_vertices=sum(s==zero for s in sig.values()),
             finite_components=len(sizes),component_size_histogram=sorted(Counter(sizes).items()),source_component_size=len(root),
             source_component_carrier_vertices=sum(i in P for i in root),zero_charge_outside_source_component=missing)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--graph',required=True);a=p.parse_args();c=controls()
 print(json.dumps(dict(status='PROVED algebra and FINITE-EXACT graph projection; paths guaranteed by companion theorem, not individually constructed',controls=c,residue_boundary=residue_boundary(),golden_probe=golden_probe(),golden_generator=golden_generator(),graph=run(a.graph),checks=checks),sort_keys=True,separators=(',',':')))
