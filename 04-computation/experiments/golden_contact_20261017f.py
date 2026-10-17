#!/usr/bin/env python3
"""A golden contact outside the biquadratic contact carrier; exact full graph.

Complete graph formula is checked independently by all-pair reduction and
exact norm tests in F[s]/(s^2-D). Original ports stay first and literal.
"""
import json, argparse, hashlib
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
from collections import Counter
import numpy as np
import sympy as sp
from pysat.solvers import Solver
from escaped_ports_20261017d import make, need
from golden_boundary_20261017c import PORT_EDGES, patterns


def real_bounds():
 h0,h1=Q(47,10),Q(24,5);w0,w1=h0-3,h1-3
 P=lambda w:w**3+2*w*w-4*w-4
 need(P(w0)<0<P(w1),'isolated real cubic branch')
 need(Q(33,10)**2<11<Q(10,3)**2 and Q(11,5)**2<5<Q(9,4)**2,'rational radical bounds')
 need((h0*h0-11)/(h0*h0+11)>Q(1,3) and (h1*h1-11)/(h1*h1+11)<Q(9,25),'real part of u')
 need(2*h0*Q(33,10)/(h1*h1+11)>Q(9,10) and 2*h1*Q(10,3)/(h0*h0+11)<1,'imaginary part of u')
 lower=Q(13,5)+1+2*(Q(13,10)*Q(1,3)-1)
 upper=Q(21,8)+1+2*(Q(21,16)*Q(9,25)-Q(9,10)**2)
 need(0<lower<upper<4,'genuine complex contact at the intended real branch')
 return dict(w_interval=[str(w0),str(w1)],N_interval=[str(lower),str(upper)])


def build(kind):
 real=real_bounds()
 f,z,j,alpha,ports,units,eta,conj,escape,unused=make('sqrt3')
 u=units[0];baseD=f.add(f.one,eta);d=f.add(f.one,alpha);shift_e=f.one;Z=f.add(d,f.mul(u,shift_e))
 bar=lambda a:f.apply(a,conj)
 N=f.mul(Z,bar(Z));D=f.sub(f.one,f.scale(f.inv(N),4))
 need(bar(D)==D,'D is real')
 rho=f.scale(f.add(f.scalar(5),j),Q(1,6))
 M=[f.zero,f.one,eta,baseD,rho,f.mul(rho,eta),f.mul(rho,baseD)]
 # Find a degree-one place with nonsquare D; denominators are checked.
 def reduce(a,conjugated=False):
  if conjugated:a=bar(a)
  out=0
  for coeff,ex in zip(a,f.basis):
   need(coeff.denominator%p!=0,'good-prime integral coefficients')
   val=coeff.numerator*pow(coeff.denominator,-1,p)
   for x,k in zip(roots,ex):val=val*pow(x,k,p)%p
   out=(out+val)%p
  return out
 witness=None
 for p in sp.primerange(31,1000):
  if p%15!=1:continue
  zs=[x for x in range(2,p) if pow(x,15,p)==1 and pow(x,3,p)!=1 and pow(x,5,p)!=1]
  js=sp.sqrt_mod(-11,p,all_roots=True)
  ws=[x for x in range(p) if (x**3+2*x*x-4*x-4)%p==0]
  for roots in product(zs,js,ws):
   if any(c.denominator%p==0 for c in D):continue
   dn=reduce(D);nn=reduce(N)
   if dn and pow(dn,(p-1)//2,p)==p-1:witness=(int(p),list(map(int,roots)),int(dn),int(nn));break
  if witness:break
 need(witness is not None,'genuine contact quadratic extension')
 p,roots,dn,nn=witness
 for poly,r in zip(f.polys,roots):need(sum(c*pow(r,k,p) for k,c in enumerate(poly))%p==0,'split-place defining relation')
 # x=A+B*s, s^2=D, conjugation sends s to -s.
 def norm(x):
  a,b=x;aa,bb=bar(a),bar(b)
  return f.sub(f.mul(a,aa),f.mul(D,f.mul(b,bb))),f.sub(f.mul(b,aa),f.mul(a,bb))
 v=(f.scale(Z,Q(-1,2)),f.scale(Z,Q(1,2)))
 need(norm(v)==(f.one,f.zero),'unit contact generator')
 need(norm((f.add(Z,v[0]),v[1]))==(f.one,f.zero),'designed three-coordinate unit edge')
 A=list(ports);seen=set(A)
 for port in ports:
  for m in M:
   t=f.add(port,m)
   if t not in seen:seen.add(t);A.append(t)
 addresses=[(i,0,0) for i in range(8)]
 addresses += [a for a in product(range(len(A)),range(7),range(7)) if a not in set(addresses)]
 index={a:i for i,a in enumerate(addresses)}
 coords=[]
 for a,b,c in addresses:
  coords.append((f.add(f.add(A[a],f.mul(u,M[b])),f.mul(v[0],M[c])),f.mul(v[1],M[c])))
 need(len(set(coords))==len(coords),'injective coordinate addresses')
 def unit_edges(pts):
  out=[]
  for a in range(len(pts)):
   for b in range(a+1,len(pts)):
    delta=f.sub(pts[b],pts[a])
    if f.mul(delta,bar(delta))==f.one:out.append((a,b))
  return out
 EA,EM=unit_edges(A),unit_edges(M)
 need(len(EM)==11,'Moser spindle positive control')
 edges=set()
 for a,b in EA:
  for c,e in product(range(7),repeat=2):edges.add(tuple(sorted((index[a,c,e],index[b,c,e]))))
 for a,b in EM:
  for c,e in product(range(len(A)),range(7)):
   edges.add(tuple(sorted((index[c,a,e],index[c,b,e]))))
   edges.add(tuple(sorted((index[c,e,a],index[c,e,b]))))
 axis_count=len(edges)
 unitsM=sorted({f.sub(M[b],M[a]) for a,b in EM}|{f.sub(M[a],M[b]) for a,b in EM})
 IA={v:i for i,v in enumerate(A)};IM={v:i for i,v in enumerate(M)}
 mixed=[]
 for i,(a,b,c) in enumerate(addresses):
  for delta in unitsM:
   aa=IA.get(f.add(A[a],f.mul(d,delta)));bb=IM.get(f.add(M[b],f.mul(shift_e,delta)));cc=IM.get(f.add(M[c],delta))
   if None not in (aa,bb,cc):
    k=index[aa,bb,cc]
    if i<k:mixed.append((i,k));edges.add((i,k))
 need(len(edges)==axis_count+len(mixed),'mixed edges differ from all axis edges')
 # Independent all-pair sieve evaluates both embedding and conjugate embedding.
 # At this split base place they are NOT Frobenius conjugates; keep both.
 r=np.array([[reduce(a),reduce(b),reduce(a,True),-reduce(b,True)%p] for a,b in coords],dtype=np.int64)
 actual=[];cache={};candidates=0
 for a in range(len(coords)):
  ds=(r[a+1:]-r[a])%p
  good=np.flatnonzero(((ds[:,0]*ds[:,2]+dn*ds[:,1]*ds[:,3])%p==1)&((ds[:,0]*ds[:,3]+ds[:,1]*ds[:,2])%p==0))
  for off in good:
   b=a+1+int(off);candidates+=1
   delta=(f.sub(coords[b][0],coords[a][0]),f.sub(coords[b][1],coords[a][1]))
   if delta not in cache:cache[delta]=norm(delta)==(f.one,f.zero)
   if cache[delta]:actual.append((a,b))
 need(set(actual)==edges,'all-pair exact norm graph equals the structural graph formula')
 edges=sorted(edges)
 need([e for e in edges if e[1]<8]==sorted(tuple(sorted(e)) for e in PORT_EDGES),'original eight ports unchanged')
 neighbours=[set() for _ in coords]
 for a,b in edges:neighbours[a].add(b);neighbours[b].add(a)
 live=set(range(len(coords)));order=[]
 while True:
  batch=[a for a in sorted(live) if a>=8 and len(neighbours[a]&live)<=5]
  if not batch:break
  for a in batch:live.remove(a);order.append(a)
 pats=patterns(8,PORT_EDGES);witnesses=[];blocked=[]
 with Solver(name='g3') as sol:
  for a in live:
   sol.add_clause([6*a+c+1 for c in range(6)])
   for c in range(6):
    for b in range(c):sol.add_clause([-6*a-c-1,-6*a-b-1])
  for a,b in edges:
   if a in live and b in live:
    for c in range(6):sol.add_clause([-6*a-c-1,-6*b-c-1])
  for pat in pats:
   if not sol.solve(assumptions=[6*a+c+1 for a,c in enumerate(pat)]):blocked.append(pat);continue
   model=set(v for v in sol.get_model() if v>0);col=[-1]*len(coords)
   for a in live:col[a]=next(c for c in range(6) if 6*a+c+1 in model)
   for a in reversed(order):col[a]=next(c for c in range(6) if c not in {col[b] for b in neighbours[a]})
   need(col[:8]==list(pat) and all(col[a]!=col[b] for a,b in edges),'all literal ports and edge-colour constraints')
   witnesses.append(''.join(map(str,col)))
 need(not blocked,'potential forcing detected; retain a separate certified UNSAT proof before promotion')
 with Solver(name='g3') as four:
  for a in range(len(coords)):four.add_clause([4*a+c+1 for c in range(4)])
  for a,b in edges:
   for c in range(4):four.add_clause([-4*a-c-1,-4*b-c-1])
  need(four.solve(),'finite golden graph four-colour positive control')
  model={v for v in four.get_model() if v>0};col4=[next(c for c in range(4) if 4*a+c+1 in model) for a in range(len(coords))]
  need(all(col4[a]!=col4[b] for a,b in edges),'explicit golden four-colouring')
 origin={index[IA[m],b,c] for m in M for b,c in product(range(7),repeat=2)}
 return dict(status='FINITE-EXACT induced graph; all boundary extensions witnessed',vertices=len(coords),edges=edges,
  kind=kind,field_dimension=96,real_embedding_bounds=real,four_colour=col4,base_degree=48,irreducibility_place=dict(p=p,zeta15=roots[0],j=roots[1],w=roots[2],N=nn,D=dn),
  first_axis_vertices=len(A),first_axis_edges=len(EA),spindle_edges=len(EM),axis_edges=axis_count,mixed_edges=mixed,
  origin_vertices=len(origin),origin_edges=sum(a in origin and b in origin for a,b in edges),origin_mixed=sum(a in origin and b in origin for a,b in mixed),
  all_pair_sieve_candidates=candidates,exact_distinct_differences=len(cache),relative_six_core=len(live)-8,peeling_order=order,
  source_addresses=addresses,first_axis_coordinates=[[str(c) for c in a] for a in A],spindle_coordinates=[[str(c) for c in a] for a in M],
  quadratic_relation=[str(c) for c in D],contact_sum=[str(c) for c in Z],patterns=pats,witnesses=witnesses)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--kind',choices=['golden'],default='golden');a=p.parse_args()
 d=build(a.kind);Path(a.output).write_text(json.dumps(d,sort_keys=True,separators=(',',':'))+'\n')
 print(json.dumps({k:v for k,v in d.items() if k not in ('edges','mixed_edges','source_addresses','first_axis_coordinates','spindle_coordinates','quadratic_relation','contact_sum','patterns','witnesses','peeling_order','four_colour')},sort_keys=True))
