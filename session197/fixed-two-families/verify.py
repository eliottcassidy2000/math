from core import inventories,cycles,literal_h,root_profile,glue,birth,require
from families import g,scalar,word,graph,prefix_graph,inverse
from itertools import product,combinations
from pathlib import Path
import json

def profile4(T,r):
 p=inventories(T);d=p['allH'][((1<<len(T))-1)^(1<<r)]
 return p['H'],p['e'][r],p['R'][r],d,p['s'][r]
def converse_map(T,U):
 n=len(T);f=[1,0]+list(range(n-1,1,-1))
 return len(U)==n and all(((T[a]>>b)&1)==((U[f[b]]>>f[a])&1)for a in range(n)for b in range(n)if a!=b)
def labelled(n):
 E=list(combinations(range(n),2))
 for b in range(1<<len(E)):
  T=[0]*n
  for k,(a,c)in enumerate(E):T[a if b>>k&1 else c]|=1<<(c if b>>k&1 else a)
  yield T

def main():
 seed_checks=0;data=[]
 for u in range(2,33):
  D,p=g(u),u%2;U=prefix_graph(u);I=inventories(U);h=D+u-p
  require((I['H'],I['e'][0],I['R'][0],I['e'][1],I['R'][1],I['s'][0],I['s'][1],cycles(U))==(h,D-p,D+u,1,2*u-1,1,1,1),'prefix identities')
  T,_=birth(U,[0,1],2,2)
  require(profile4(T,0)==(D+3*u-p-2,u-1,2*D+4*u-p-1,D+u-p,2)and cycles(T)==1,'N seed')
  T,_=birth(U,[0,1],0,2);T,_=birth(T,[0,1],3,2)
  require(profile4(T,0)==(3*D+6*u-2-2*p,D+u-1,3*D+8*u-1-p,D+2*u,2)and cycles(T)==2,'M seed')
  seed_checks+=3
 for u in range(1,33):
  D,p=g(u),u%2;T=graph('P',u,1)
  require(profile4(T,1)==(5*D+9*u-p,D+3*u,5*D+10*u-p,D+4*u,2)and cycles(T)==2,'P seed');seed_checks+=1
 cell=0
 for n in range(1,4):
  for T in labelled(n):
   for a in range(n):
    h,e,R=root_profile(T,a);V,_=glue(T,{a},[0,5,1])
    require(root_profile(V,n+1)==(h+4*R,R,4*h-e+6*R),'literal transitive D cell');cell+=1
 total=literal=0
 for kind in ('N','M','P'):
  lo=1 if kind=='P' else 2
  for u in range(lo,9):
   for v in range(lo,9):
    T=graph(kind,u,v);h=scalar(kind,u,v);require(inventories(T)['H']==h and cycles(T)==2,'family scalar/cycles')
    require(converse_map(T,graph(kind,v,u)),'actual converse graph map')
    # Removing the two permanent ports leaves exactly the birth order.
    require(all(T[a]>>b&1 for a in range(2,len(T))for b in range(a+1,len(T))),'feedback two')
    if len(T)<=8:require(literal_h(T)==h,'literal permutation H');literal+=1
    if kind=='N':require(h>=3*u*v+1,'N product bound')
    elif kind=='M':require(h>=6*u*v+u+v+1,'M product bound')
    else:require(h>=6*u*v+3*u+3*v+1,'P product bound')
    total+=1
  # Oversized rectangle is independent of the claimed inverse region.
  raw={}
  for u in range(lo,102):
   for v in range(u,102):
    h=scalar(kind,u,v)
    if h<=101:raw.setdefault(h,[]).append((u,v))
  for h in range(1,102):require(inverse(kind,h)==raw.get(h,[]),'finite inverse test')
 for kind,u,v,h in [('N',2,2,13),('M',2,2,29),('M',4,2,55),('P',1,24,231),('P',8,3,247)]:
  T=graph(kind,u,v);require(inventories(T)['H']==h and cycles(T)==2,'displayed witness')
  data.append(dict(family=kind,u=u,v=v,H=h,C=2,rows=T,word=word(kind,u,v)))
 result=dict(status='PASS',prefix_gate_profile_checks=seed_checks,independent_D_cell_graph_checks=cell,physical_family_graphs=total,literal_permutation_checks=literal,converse_maps=total,inverse_rectangle_limit=101,examples=data)
 text=json.dumps(result,sort_keys=True,indent=2)+'\n';(Path(__file__).parent/'VALIDATION.json').write_text(text);print(text,end='')
if __name__=='__main__':main()
