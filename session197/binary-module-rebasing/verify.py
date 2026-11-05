from pathlib import Path
from itertools import product
import hashlib,json
from core import inventories,require,tournament
from rebase import *

P=Path(__file__).resolve().parent
def main():
 grids=[]
 for t,j,m in product((1,2),(1,2),(1,2,3)):
  T,p=construction(t,j,m);z=inventories(T)
  require(z['H']==hp(t,j,m),'positive tensor grid')
  grids.append(dict(t=t,j=j,m=m,rows=T,H=z['H']))
 responses=[]
 for t,j in product((1,2,3),(1,2)):
  T,p=construction(t,j,1);z=inventories(T);r=p['root'];W=append(T,{r});w=inventories(W)
  require((z['H'],z['e'][r],z['R'][r])==(hp(t,j,1),ending(t,j),cover(t,j)),'response grid')
  require(w['H']==zcount(t,j)==z['H']-z['e'][r]+z['R'][r],'Z grid')
  responses.append(dict(t=t,j=j,rows=T,profile=[z['H'],z['e'][r],z['R'][r]],Z_rows=W,Z_H=w['H']))
 controls=0;maps=0
 for t,j,m in product(range(1,7),range(1,7),range(1,6)):
  T,p=construction(t,j,m);require(tournament(T),'tournament');cycle=strong_cycle(t,j,m,p)
  require(sorted(cycle)==list(range(len(T)))and all(T[a]>>b&1 for a,b in zip(cycle,cycle[1:]+cycle[:1])),'Hamilton cycle')
  for u in p['U']:require(cover_check(T,u,{p['b']}),'U rank1 cover')
  for z in p['Z']:require(cover_check(T,z,{p['a']}),'Z rank1 cover')
  for r in p['R']:require(cover_check(T,r,{p['a'],p['b']}),'R rank2 cover')
  if m==1:
   W=append(T,{p['root']});new=len(T)
   for u in p['U']:require(cover_check(W,u,{p['b'],new}),'augmented U rank2 cover')
   for z in p['Z']:require(cover_check(W,z,{p['a']}),'augmented Z rank1 cover')
  if j>=2:move_right(t,j,m);maps+=1
  if t>=2:move_left(t,j,m);maps+=1
  controls+=1
 incomparability=[]
 for j in range(2,33):
  h=hp(1,j,1);ea,Ra=ending(1,j),cover(1,j);eb,Rb=ending(2,j-1),cover(2,j-1)
  A=h-ea+Ra;B=h-eb+Rb
  require(A<B<h+2*Ra and B<3*Ra and h<=Ra and h<=Rb,'one-gadget separation inequalities')
  incomparability.append(dict(j=j,H=h,A=A,B=B))
 T,p=construction(1,2,1);U,q=construction(2,1,1)
 require(sorted(x.bit_count()for x in T)!=sorted(x.bit_count()for x in U),'nonisomorphic hostile')
 out=dict(status='PASS',count_tensor_grid=grids,response_tensor_grid=responses,rank_and_cycle_controls=controls,physical_replacement_maps=maps,one_gadget_separation=incomparability,example_degree_sequences=[sorted(x.bit_count()for x in T),sorted(x.bit_count()for x in U)])
 encoded=json.dumps(out,indent=2,sort_keys=True)+'\n';(P/'VALIDATION.json').write_text(encoded);print(encoded,end='')
if __name__=='__main__':main()
