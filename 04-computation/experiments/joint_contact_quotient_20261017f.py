#!/usr/bin/env python3
"""Exact joint-contact direction classification and a four-state fibre quotient.

The quotient is a sufficient colouring certificate, never a forcing test.
The source is the full A + u M + v1 M + v2 M + v3 M grid.
"""
import json,argparse,hashlib
from pathlib import Path
from itertools import product,combinations
from fractions import Fraction as Q
import sympy as sp
from pysat.solvers import Solver
from escaped_ports_20261017d import make,need
from golden_boundary_20261017c import PORT_EDGES,patterns
from contact_catalog_verify_20261017e import add,mul,norm,sub,scale


ONE=(Q(1),Q(0),Q(0),Q(0))
def bar(v):return (v[0],-v[1],-v[2],v[3])
def fm(a,b):
 v=0
 while b:
  if b&1:v^=a
  b>>=1;a<<=1
  if a&4:a^=7
 return v
def f4(v,eta,rho):
 a,b,c,d=v
 coeff=[a-5*b-c+5*d,2*c-10*d,6*b-6*d,12*d]
 out=0
 for q,r in zip(coeff,[1,eta,rho,fm(eta,rho)]):
  need(q.denominator%2==1,'coefficient integral in the 1,eta,rho,eta*rho basis at2')
  if q.numerator&1:out^=r
 return out


def algebra(catalog,selected):
 V=[tuple(Q(c) for c in v) for v in catalog['differences']]
 ds=[V[catalog['rows'][i]['d']] for i in selected];es=[V[catalog['rows'][i]['e']] for i in selected]
 determinants=[]
 for i,j in combinations(range(3),2):
  d,e,D,E=ds[i],es[i],ds[j],es[j]
  a=[mul(e,bar(D)),add(mul(d,bar(D)),mul(e,bar(E))),mul(d,bar(E))]
  b=[mul(bar(d),E),add(mul(bar(d),D),mul(bar(e),E)),mul(bar(e),D)]
  minors=[(k,l,sub(mul(a[k],b[l]),mul(a[l],b[k]))) for k,l in combinations(range(3),2)]
  good=next((k,l,v) for k,l,v in minors if any(v))
  determinants.append(dict(pair=[i,j],rows=list(good[:2]),minor=[str(x) for x in good[2]]))
 # Three independent nonsquare-character rows prove all seven nontrivial
 # squareclass products nonsquare, even after adjoining zeta5.
 witnesses=[];span={0}
 for p in sp.primerange(31,2000):
  if p%15!=1:continue
  js=sp.sqrt_mod(-11,p,all_roots=True)
  if not js:continue
  zs=[z for z in range(2,p) if pow(z,15,p)==1 and pow(z,5,p)!=1 and pow(z,3,p)!=1]
  ws=[w for w in range(p) if (w**3+2*w*w-4*w-4)%p==0]
  for z,j,w in product(zs,js,ws):
   h=(1+2*pow(z,5,p))%p
   def ev(v,conjugated=False):
    if conjugated:v=bar(v)
    return sum(c.numerator*pow(c.denominator,-1,p)*r for c,r in zip(v,[1,j,h,h*j]))%p
   den=(w+3+j)%p;num=(w+3-j)%p
   if not den or not num:continue
   u=num*pow(den,-1,p)%p;ui=pow(u,-1,p);ns=[];Dvals=[]
   for d,e in zip(ds,es):
    N=(ev(d)+u*ev(e))*(ev(d,True)+ui*ev(e,True))%p;ns.append(N)
    if not N:break
    Dvals.append((1-4*pow(N,-1,p))%p)
   if len(Dvals)!=3 or 0 in Dvals:continue
   mask=sum((pow(D,(p-1)//2,p)==p-1)<<i for i,D in enumerate(Dvals))
   if mask not in span:
    witnesses.append(dict(prime=int(p),zeta15=int(z),sqrt_minus11=int(j),w=int(w),N=ns,D=Dvals,mask=mask))
    span|={v^mask for v in list(span)}
   if len(span)==8:return ds,es,dict(squareclass_witnesses=witnesses,character_span=sorted(span),cross_minors=determinants,ambient_degree=384,direction_families=8)
 raise RuntimeError('independent quadratic squareclasses not yet witnessed')


def main(output):
 old=Path('04-computation/experiments/marked_contact_current_20261017e')
 catalog=json.loads((old/'catalog.json').read_text());selected=[1,162,315]
 ds,es,alg=algebra(catalog,selected)
 f,z,j,alpha,ports,units,eta,cj,escape,reduce=make('sqrt3')
 h=f.sub(f.scale(eta,2),f.one);hj=f.mul(h,j)
 embed=lambda v:sum_field(f,[f.scale(x,c) for x,c in zip([f.one,j,h,hj],v)])
 d=[embed(v) for v in ds];e=[embed(v) for v in es]
 rho=f.scale(f.add(f.scalar(5),j),Q(1,6));D=f.add(f.one,eta)
 M=[f.zero,f.one,eta,D,rho,f.mul(rho,eta),f.mul(rho,D)]
 A=list(ports)
 for a,b in product(ports,M):
  c=f.add(a,b)
  if c not in A:A.append(c)
 IA={v:i for i,v in enumerate(A)};IM={v:i for i,v in enumerate(M)}
 unit=lambda a:f.mul(a,f.apply(a,cj))==f.one
 EA=[(i,j) for i,j in combinations(range(len(A)),2) if unit(f.sub(A[j],A[i]))]
 EM=[(i,j) for i,j in combinations(range(7),2) if unit(f.sub(M[j],M[i]))]
 U=[tuple(Q(c) for c in v) for v in catalog['units']];UF=[embed(v) for v in U]
 need(len(A)==54 and len(EA)==146 and len(EM)==11,'inherited exact seed and ports')
 contacts=[];source_counts=[]
 for i in range(3):
  rows=[];count=0
  for qi,q in enumerate(UF):
   da,db=f.mul(d[i],q),f.mul(e[i],q)
   aa=[(a,IA[f.add(v,da)]) for a,v in enumerate(A) if f.add(v,da) in IA]
   bb=[(b,IM[f.add(v,db)]) for b,v in enumerate(M) if f.add(v,db) in IM]
   cc=[(c,IM[f.add(v,q)]) for c,v in enumerate(M) if f.add(v,q) in IM]
   count+=len(aa)*len(bb)*len(cc)
   if bb and cc:
    rows += [(a,b,qi) for a,b in aa]
  need(count%2==0,'two directed copies per contact edge');contacts.append(rows);source_counts.append(count//2*49)
 # An additive F4-valued homomorphism on each spindle coordinate.
 # We search a tiny quotient, preserving every contact edge and all ports.
 attempts=[];chosen=None
 pats=patterns(8,PORT_EDGES)
 for er,rr,r,t1,t2,t3 in product((2,3),(2,3),(1,2,3),(1,2,3),(1,2,3),(1,2,3)):
  ts=[t1,t2,t3];phi=lambda v:f4(v,er,rr)
  if any(phi(ds[i])^fm(r,phi(es[i]))^ts[i]==0 for i in range(3)):continue
  edges=set()
  for a in range(len(A)):
   for x,y in combinations(range(4),2):edges.add(tuple(sorted((4*a+x,4*a+y))))
  for a,b in EA:
   for x in range(4):edges.add((4*a+x,4*b+x))
  for i,rows in enumerate(contacts):
   for a,b,qi in rows:
    delta=fm(fm(r,phi(es[i]))^ts[i],phi(U[qi]))
    for x in range(4):edges.add(tuple(sorted((4*a+x,4*b+(x^delta)))))
  edges=sorted(edges);n=4*len(A)
  # Try a four-colouring first as a cheap hostile to forcing.
  with solver(n,edges,4) as sol:
   ok=sol.solve()
   colour=extract(sol,n,4) if ok else None
  attempts.append(dict(eta=er,rho=rr,r=r,t=ts,edges=len(edges),four_colourable=ok))
  if not ok:continue
  need(all(colour[a]!=colour[b] for a,b in edges),'literal four-colour certificate')
  # Put the eight literal ports first for the independent receipt checker.
  order=[4*a for a in range(8)]+[a for a in range(n) if a not in {4*b for b in range(8)}]
  rev={v:i for i,v in enumerate(order)};ed=sorted(tuple(sorted((rev[a],rev[b]))) for a,b in edges)
  if [q for q in ed if q[1]<8]!=sorted(tuple(sorted(q)) for q in PORT_EDGES):
   attempts[-1]['extra_port_edge']=True;continue
  witnesses=[];blocked=[]
  with solver(n,ed,6) as sol:
   for pat in pats:
    if not sol.solve(assumptions=[6*a+c+1 for a,c in enumerate(pat)]):blocked.append(pat);continue
    col=extract(sol,n,6);need(col[:8]==list(pat) and all(col[a]!=col[b] for a,b in ed),'quotient boundary witness')
    witnesses.append(''.join(map(str,col)))
  attempts[-1]['blocked_patterns']=len(blocked)
  if blocked:continue
  coeffs=[tuple(Q(c) for c in v) for v in [[0,0,0,0],[1,0,0,0],['1/2',0,'1/2',0],['3/2',0,'1/2',0],['5/6','1/6',0,0],['5/12','1/12','5/12','1/12'],['5/4','1/4','5/12','1/12']]]
  pm=[phi(v) for v in coeffs]
  source_hash=hashlib.sha256();quotient_vertices=set()
  for a,b,c1,c2,c3 in product(range(len(A)),range(7),range(7),range(7),range(7)):
   x=fm(r,pm[b])^fm(t1,pm[c1])^fm(t2,pm[c2])^fm(t3,pm[c3]);v=rev[4*a+x]
   quotient_vertices.add(v);source_hash.update(v.to_bytes(2,'little'))
  # Independent replay accepts exact peeling order, which may be empty.
  adj=[set() for _ in range(n)]
  for a,b in ed:adj[a].add(b);adj[b].add(a)
  live=set(range(n));peel=[]
  while True:
   batch=[a for a in sorted(live) if a>=8 and len(adj[a]&live)<=5]
   if not batch:break
   for a in batch:live.remove(a);peel.append(a)
  chosen=dict(eta=er,rho=rr,r=r,t=ts,spindle_reduction=pm,source_projection_sha256=source_hash.hexdigest(),
   vertices=n,edges=ed,patterns=pats,witnesses=witnesses,four_colour=[colour[a] for a in order],
   peeling_order=peel,relative_six_core=len(live)-8,quotient_order=order)
  break
 need(chosen is not None,'no complete positive quotient found in the declared bank')
 result=dict(status='PROVED direction classification + FINITE-EXACT source-preserving positive quotient; forcing OPEN',
  algebra=alg,selected=selected,source_vertices=len(A)*7**4,source_axis_edges=len(EA)*7**4+4*len(A)*len(EM)*7**3,
  source_mixed_edges=source_counts,source_total_edges=len(EA)*7**4+4*len(A)*len(EM)*7**3+sum(source_counts),
  contact_rows=contacts,attempts=attempts,**chosen)
 Path(output).write_text(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ('edges','patterns','witnesses','four_colour','contact_rows','quotient_order')},sort_keys=True))


def sum_field(f,values):
 r=f.zero
 for v in values:r=f.add(r,v)
 return r
def solver(n,edges,k):
 s=Solver(name='g3')
 for a in range(n):
  s.add_clause([k*a+c+1 for c in range(k)])
  for c,d in combinations(range(k),2):s.add_clause([-k*a-c-1,-k*a-d-1])
 for a,b in edges:
  for c in range(k):s.add_clause([-k*a-c-1,-k*b-c-1])
 return s
def extract(s,n,k):
 m=set(v for v in s.get_model() if v>0);return [next(c for c in range(k) if k*a+c+1 in m) for a in range(n)]

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args();main(a.output)
