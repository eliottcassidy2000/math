#!/usr/bin/env python3
"""A signature-directed quadratic escape and an anchored anti-invariant unit walk.
All old ports and all 5400 inherited vertices are retained. Only new pairs
need fresh induced-edge tests; old-old distances cannot change on extension.
"""
import json,argparse
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
from pysat.solvers import Solver
from escaped_ports_20261017d import make,need,P
from golden_boundary_20261017c import PORT_EDGES,patterns
from shared_golden_contact_20261017g import real_bounds


def signature():
 real=real_bounds();lower=Q(13,5)+1+2*(Q(13,10)*Q(1,3)+Q(9,10)**2)
 need(lower==Q(913,150) and lower>4,'conjugate embedding discriminant positive')
 places=[]
 for z,j,w in ((7,19,25),(10,19,25)):
  p=31
  need(pow(z,15,p)==1 and pow(z,3,p)!=1 and pow(z,5,p)!=1,'primitive fifteenth root')
  need((j*j+11)%p==0 and (w**3+2*w*w-4*w-4)%p==0,'base defining relations')
  alpha=pow(z,3,p);u=(w+3-j)*pow(w+3+j,-1,p)%p
  n=(1+alpha+u)*(1+pow(alpha,-1,p)+pow(u,-1,p))%p
  d=(1-4*pow(n,-1,p))%p;r=(4-n)%p
  places.append(dict(p=p,z=z,j=j,w=w,N=n,D=d,R=r,characters=[pow(d,15,p),pow(r,15,p)]))
 need([x['characters'] for x in places]==[[30,30],[30,1]],'independent quadratic characters of D and R')
 return dict(real_bounds=real,other_N_lower=str(lower),places=places,extension_degree=192)


def run():
 cert=signature();old=json.load(open('04-computation/experiments/root_positive_transport_20261017g/golden_shared.json'))
 f,z,j,alpha,ports,units,eta,conj,escape,reduce=make('sqrt3');bar=lambda a:f.apply(a,conj)
 Z=tuple(map(Q,old['contact_sum']));D=tuple(map(Q,old['quadratic_relation']))
 r=f.scale(f.mul(Z,f.inv(bar(Z))),-1);ri=bar(r)
 need(f.mul(r,ri)==f.one,'new unit square r has norm one')
 flip=f.substitution([z,f.scale(j,-1),f.gen(2)])
 need(f.apply(Z,flip)==f.add(f.add(f.one,alpha),bar(units[0])),'real embedding reverses u only')
 A=[tuple(map(Q,a)) for a in old['first_axis_coordinates']];C=[tuple(map(Q,a)) for a in old['spindle_coordinates']]
 v=(f.scale(Z,Q(-1,2)),f.scale(Z,Q(1,2)))
 uc=[f.mul(units[0],c) for c in C];vc=[(f.mul(v[0],c),f.mul(v[1],c)) for c in C]
 # Coordinates are a0+a1*s+b*e. e^2=r; conjugation e -> e/r.
 base=[(f.add(f.add(A[a],uc[b]),vc[c][0]),vc[c][1],f.zero) for a,b,c in old['source_addresses']]
 a=ports[5]
 path=[f.zero,f.one,f.sub(f.one,a),f.scalar(Q(1,4)),f.scalar(Q(5,4)),f.sub(f.scalar(Q(5,4)),a),f.scalar(Q(1,2))]
 L=path+[f.scale(x,-1) for x in path[1:]]
 need(len(set(L))==13,'anchored walk has thirteen distinct vertices')
 coords=list(base);index={x:i for i,x in enumerate(coords)};attachments=[]
 for port in ports:
  ids=[]
  for b in L:
   x=(port,f.zero,b)
   if x not in index:index[x]=len(coords);coords.append(x)
   ids.append(index[x])
  attachments.append(ids)
 def norm(x):
  a,c,b=x;aa,cc,bb=bar(a),bar(c),bar(b)
  return (f.add(f.sub(f.mul(a,aa),f.mul(D,f.mul(c,cc))),f.mul(b,bb)),
          f.sub(f.mul(c,aa),f.mul(a,cc)),
          f.add(f.mul(b,aa),f.mul(f.mul(a,bb),ri)),
          f.add(f.scale(f.mul(b,cc),-1),f.mul(f.mul(c,bb),ri)))
 unit=(f.one,f.zero,f.zero,f.zero)
 def mulmod(a,b):return np.column_stack(((a[:,0]*b[:,0]-a[:,1]*b[:,1])%P,(a[:,0]*b[:,1]+a[:,1]*b[:,0])%P))
 def barmod(a):return a*np.array([1,-1],dtype=np.int64)%P
 rr=np.array([[*reduce(a),*reduce(c),*reduce(b)] for a,c,b in coords],dtype=np.int64)
 dn,im=reduce(D);need(im==0,'real discriminant reduction');ir=np.array(reduce(ri),dtype=np.int64)
 n=len(base);edges=set(map(tuple,old['edges']));cache={};candidates=0
 for i in range(n,len(coords)):
  delta=(rr[:i]-rr[i])%P;a=delta[:,:2];c=delta[:,2:4];b=delta[:,4:]
  c0=(a[:,0]**2+a[:,1]**2-dn*(c[:,0]**2+c[:,1]**2)+b[:,0]**2+b[:,1]**2)%P
  c1=(c[:,1]*a[:,0]-c[:,0]*a[:,1])%P
  e0=(mulmod(b,barmod(a))+mulmod(mulmod(a,barmod(b)),np.broadcast_to(ir,a.shape)))%P
  e1=(-mulmod(b,barmod(c))+mulmod(mulmod(c,barmod(b)),np.broadcast_to(ir,a.shape)))%P
  good=np.flatnonzero((c0==1)&(c1==0)&np.all(e0==0,axis=1)&np.all(e1==0,axis=1))
  for jj in good:
   j0=int(jj);candidates+=1;x=tuple(f.sub(coords[i][k],coords[j0][k]) for k in range(3))
   if x not in cache:cache[x]=norm(x)==unit
   if cache[x]:edges.add((j0,i))
 edges=sorted(edges);es=set(edges);need([e for e in edges if e[1]<8]==sorted(tuple(sorted(e)) for e in PORT_EDGES),'old literal port relation')
 cycle=[0,1,2,3,4,5,6,12,11,10,9,8,7,0]
 for ids in attachments:
  need(all(tuple(sorted((ids[a],ids[b]))) in es for a,b in zip(cycle,cycle[1:])),'two six-step paths plus conjugate unit edge')
 # sigma fixes each port and flips e; the edge joining +/-e/2 cannot fold.
 need(all(tuple(sorted((ids[6],ids[12]))) in es for ids in attachments),'eight fixed-port anti-invariant edges')
 # Exact Cartesian relation: old graph union (port graph square C13).
 predicted=set(map(tuple,old['edges']))
 for ids in attachments:
  predicted.update(tuple(sorted((ids[a],ids[b]))) for a,b in zip(cycle,cycle[1:]))
 for a,b in PORT_EDGES:
  predicted.update(tuple(sorted((attachments[a][j],attachments[b][j]))) for j in range(13))
 need(predicted==es,'all new edges form precisely the marked Cartesian odd-cycle context')
 colours=[None]*13
 for k,a in enumerate(cycle[:-1]):colours[a]=k%2 if k<12 else 2
 need(colours[0]==0 and all(colours[a]!=colours[b] for a,b in zip(cycle,cycle[1:])),'fixed-root three-colouring of the odd cycle')
 def extend(col0,mod):
  col=list(map(int,col0))+[-1]*(len(coords)-n)
  for p,ids in enumerate(attachments):
   for j,a in enumerate(ids):
    c=(int(col0[p])+colours[j])%mod
    if col[a]>=0:need(col[a]==c,'shared old ports fixed by product extension')
    col[a]=c
  need(all(0<=c<mod for c in col) and all(col[a]!=col[b] for a,b in edges),'complete constructive product colouring')
  return col
 pats=patterns(8,PORT_EDGES)
 for pat,col0 in zip(pats,old['witnesses']):need(extend(col0,6)[:8]==list(pat),'all inherited boundary patterns')
 four=extend(old['four_colour'],4)
 recipe=dict(kind='marked_cartesian_cycle',cycle_colours=colours,modulus=6)
 return dict(status='FINITE-EXACT degree-192 signature escape, complete induced extension and positive witnesses',
  signature=cert,vertices=len(coords),old_vertices=n,edges=edges,new_edges=len(edges)-len(old['edges']),attachments=attachments,
  sieve_prime=P,sieve_candidates=candidates,distinct_exact_norms=len(cache),patterns=pats,witness_recipe=recipe,four_colour=four,
  unit_square=[str(c) for c in r],path_coefficients=[[str(c) for c in b] for b in L])

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args();r=run()
 Path(a.output).write_text(json.dumps(r,sort_keys=True,separators=(',',':'))+'\n')
 print(json.dumps({k:v for k,v in r.items() if k not in ('edges','attachments','patterns','witnesses','four_colour','unit_square','path_coefficients','witness_recipe')},sort_keys=True))
