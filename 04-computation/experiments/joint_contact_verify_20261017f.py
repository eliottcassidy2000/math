#!/usr/bin/env python3
"""Replay every classified source edge, local characters, and the F4 field test.

Source completeness uses the proved direction classification in the companion
note; this replay independently rebuilds all instances of those directions.
No SAT solver or numerical tolerance is used.
"""
import argparse,json,hashlib
import sympy as sp
from fractions import Fraction as Q
from itertools import product,combinations
from pathlib import Path
from escaped_ports_20261017d import make,need
from contact_catalog_verify_20261017e import add,mul,sub


def bar(x):return x[0],-x[1],-x[2],x[3]
def ff(a,b):
 c=0
 for k in range(2):
  if b>>k&1:c^=a<<k
 return c^7 if c&4 else c

def residue(x):
 a,b,c,d=x
 basis=[a-5*b-c+5*d,2*c-10*d,6*b-6*d,12*d]
 out=0
 for coeff,v in zip(basis,[1,2,2,3]):
  need(coeff.denominator%2==1,'2-integral biquadratic coefficients')
  if coeff.numerator&1:out^=v
 return out


def run(path):
 R=json.loads(Path(path).read_text())
 cat=json.loads(Path('04-computation/experiments/marked_contact_current_20261017e/catalog.json').read_text())
 V=[tuple(map(Q,x)) for x in cat['differences']]
 ds=[V[cat['rows'][i]['d']] for i in R['selected']];es=[V[cat['rows'][i]['e']] for i in R['selected']]
 chars=[]
 for q in R['algebra']['squareclass_witnesses']:
  p,z,j,w=[q[k] for k in ['prime','zeta15','sqrt_minus11','w']]
  need(pow(z,15,p)==1 and pow(z,3,p)!=1 and pow(z,5,p)!=1 and (j*j+11)%p==0 and (w**3+2*w*w-4*w-4)%p==0,'place relations')
  h=(1+2*pow(z,5,p))%p;need((h*h+3)%p==0,'sqrt-3 branch')
  def ev(v):return sum(c.numerator*pow(c.denominator,-1,p)*r for c,r in zip(v,[1,j,h,h*j]))%p
  u=(w+3-j)*pow(w+3+j,-1,p)%p;vals=[];ns=[]
  for d,e in zip(ds,es):
   n=(ev(d)+u*ev(e))*(ev(bar(d))+pow(u,-1,p)*ev(bar(e)))%p
   ns.append(n);vals.append((1-4*pow(n,-1,p))%p)
  need(ns==q['N'] and vals==q['D'],'independent local norm and discriminant')
  mask=sum((pow(v,(p-1)//2,p)==p-1)<<i for i,v in enumerate(vals))
  need(mask==q['mask'] and all(vals),'unramified quadratic characters');chars.append(mask)
 need(all(any((a&b).bit_count()%2 for b in chars) for a in range(1,8)),'every nontrivial squareclass product detected')
 minors=[]
 for i,k in combinations(range(3),2):
  d,e,D,E=ds[i],es[i],ds[k],es[k]
  p=[mul(e,bar(D)),add(mul(d,bar(D)),mul(e,bar(E))),mul(d,bar(E))]
  q=[mul(bar(d),E),add(mul(bar(d),D),mul(bar(e),E)),mul(bar(e),D)]
  candidates=[sub(mul(p[a],q[b]),mul(p[b],q[a])) for a,b in combinations(range(3),2)]
  need(any(any(c) for c in candidates),'cross-character rank two');minors.append([[str(c) for c in v] for v in candidates])
 # Local field has u=1, eta=rho=omega; the three Hensel polynomials split simply.
 zbar=[residue(add(d,e)) for d,e in zip(ds,es)];rootsets=[]
 need(zbar==[2,1,2],'three contact sums are units at the fixed place above2')
 for z in zbar:
  roots=[v for v in range(4) if ff(v,v)^ff(z,v)^ff(z,z)==0]
  need(len(roots)==2 and 0 not in roots,'two simple nonzero Hensel roots')
  need(all(ff(v,ff(v,v))==1 for v in roots),'conjugate residue is inverse residue');rootsets.append(roots)
 # A rejected equal-norm detour, retained with its exact obstruction.
 H,Y=sp.symbols('H Y');P=H**3-7*H**2+11*H-1
 eq=sp.Poly(sp.resultant(P,33*Y-(22-H**2),H),Y)
 need(eq.all_coeffs()==[35937,-42471,12243,67],'equal-norm contact polynomial')
 normD=Q(-67,35937);need(normD.numerator*pow(normD.denominator,-1,8)%8==5,'real cubic norm is nonsquare at2')
 # Reconstruct source sets and all source directions independently.
 f,z,j,alpha,ports,units,eta,cj,escape,red=make('sqrt3');h=f.sub(f.scale(eta,2),f.one)
 def emb(x):
  v=f.zero
  for c,t in zip(x,[f.one,j,h,f.mul(h,j)]):v=f.add(v,f.scale(t,c))
  return v
 rho=f.scale(f.add(f.scalar(5),j),Q(1,6));D=f.add(f.one,eta)
 M=[f.zero,f.one,eta,D,rho,f.mul(rho,eta),f.mul(rho,D)]
 A=list(ports)
 for a,b in product(ports,M):
  v=f.add(a,b)
  if v not in A:A.append(v)
 IA={v:i for i,v in enumerate(A)};IM={v:i for i,v in enumerate(M)}
 unit=lambda x:f.mul(x,f.apply(x,cj))==f.one
 EA=[(a,b) for a,b in combinations(range(len(A)),2) if unit(f.sub(A[b],A[a]))]
 EM=[(a,b) for a,b in combinations(range(7),2) if unit(f.sub(M[b],M[a]))]
 U=sorted({f.sub(M[b],M[a]) for a,b in EM}|{f.sub(M[a],M[b]) for a,b in EM})
 pm=R['spindle_reduction'];need(pm==[0,1,2,3,2,3,1],'literal spindle reductions')
 rev={v:i for i,v in enumerate(R['quotient_order'])};edges={tuple(e) for e in R['edges']}
 def project(a):
  x=ff(R['r'],pm[a[1]])
  for t,c in zip(R['t'],a[2:]):x^=ff(t,pm[c])
  return rev[4*a[0]+x]
 digest=hashlib.sha256()
 for a in product(range(len(A)),range(7),range(7),range(7),range(7)):digest.update(project(a).to_bytes(2,'little'))
 need(digest.hexdigest()==R['source_projection_sha256'],'every source address retained')
 count=0
 def check(a,b):
  nonlocal count
  x,y=project(a),project(b);need(x!=y and tuple(sorted((x,y))) in edges,'source edge survives quotient');count+=1
 for a,b in EA:
  for tail in product(range(7),repeat=4):check((a,)+tail,(b,)+tail)
 for slot in range(1,5):
  for a,b in EM:
   for tail in product(range(len(A)),range(7),range(7),range(7)):
    src=list(tail);dst=list(tail);src.insert(slot,a);dst.insert(slot,b);check(src,dst)
 need(count==R['source_axis_edges'],'independently counted axis edges');mixed=[]
 for slot,(d,e) in enumerate(zip(ds,es),2):
  before=count;d,e=emb(d),emb(e)
  for q in U:
   aa=[(a,IA[f.add(v,f.mul(d,q))]) for a,v in enumerate(A) if f.add(v,f.mul(d,q)) in IA]
   bb=[(b,IM[f.add(v,f.mul(e,q))]) for b,v in enumerate(M) if f.add(v,f.mul(e,q)) in IM]
   cc=[(c,IM[f.add(v,q)]) for c,v in enumerate(M) if f.add(v,q) in IM]
   for (a,Aa),(b,Bb),(c,Cc),rest in product(aa,bb,cc,product(range(7),repeat=2)):
    src=[a,b];dst=[Aa,Bb];r=iter(rest)
    for s in range(2,5):
     v=c if s==slot else next(r);src.append(v);dst.append(Cc if s==slot else v)
    if tuple(src)<tuple(dst):check(src,dst)
  mixed.append(count-before)
 need(mixed==R['source_mixed_edges'] and count==R['source_total_edges'],'every classified source edge replayed')
 for colour in [R['four_colour']]+[[int(c) for c in word] for word in R['witnesses']]:
  need(len(colour)==R['vertices'] and all(colour[a]!=colour[b] for a,b in edges),'literal quotient colouring')
 need(len(R['witnesses'])==570 and all(list(map(int,c[:8]))==p for p,c in zip(R['patterns'],R['witnesses'])),'570 literal boundary extensions')
 return dict(status='FINITE-EXACT independent replay',source_vertices=R['source_vertices'],source_edges_replayed=count,source_mixed_edges=mixed,
  squareclass_masks=chars,local_contact_residues=zbar,local_root_sets=rootsets,quotient_vertices=R['vertices'],quotient_edges=len(edges),boundary_extensions=570,
  source_projection_sha256=digest.hexdigest(),equal_norm_polynomial=list(map(int,eq.all_coeffs())),equal_norm_norm=str(normD))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('input');a=p.parse_args();print(json.dumps(run(a.input),sort_keys=True))
