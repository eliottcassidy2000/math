#!/usr/bin/env python3
"""Literal replay of rigid-guard hub theorem against complete n<=8 atlases.
Historical producer checker with path/output-only packaging adaptations.
Requires the authenticated full atlas, including its n8 subdirectory.
"""
from pathlib import Path
import json
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("--atlas-dir", type=Path, required=True)
args=parser.parse_args()
ROOT=Path(__file__).resolve().parent
OLD=args.atlas_dir

def schema(n):
 c=[(i,j) for j in range(n-1,1,-1) for i in range(j-1)]
 return c,[c.index((n-1-j,n-1-i)) for i,j in c]
def decode(t,n,c):
 a=[0]*n
 for i in range(n-1):a[i]|=1<<(i+1)
 for k,(i,j) in enumerate(c):
  x,y=(j,i) if t>>k&1 else(i,j);a[x]|=1<<y
 return a
def encode(a,p,c):return sum((1-((a[p[i]]>>p[j])&1))<<k for k,(i,j) in enumerate(c))
def paths(a):
 n=len(a);out=[]
 def rec(p,used):
  if len(p)==n:out.append(tuple(p));return
  poss=(a[p[-1]] if p else(1<<n)-1)&~used
  while poss:
   b=poss&-poss;poss^=b;rec(p+(b.bit_length()-1,),used|b)
 rec((),0);return out
def converse(a):return [((1<<len(a))-1)^r^(1<<i) for i,r in enumerate(a)]
def flip(a,edges):
 b=a[:]
 for x,y in edges:b[x]^=1<<y;b[y]^=1<<x
 return b
def partner(a,p):return flip(converse(a),zip(p,p[1:]))
def c3(a):
 n=len(a);return n*(n-1)*(n-2)//6-sum(r.bit_count()*(r.bit_count()-1)//2 for r in a)
def components(a):
 n=len(a);r=[x|(1<<i) for i,x in enumerate(a)]
 for k in range(n):
  for i in range(n):
   if r[i]>>k&1:r[i]|=r[k]
 unused=set(range(n));out=[]
 while unused:
  v=min(unused);s={w for w in unused if(r[v]>>w&1)and(r[w]>>v&1)};out.append(s);unused-=s
 out.sort(key=lambda s:-sum((a[next(iter(s))]>>w)&1 for z in out if z!=s for w in z))
 return out
def require(b,msg):
 if not b:raise RuntimeError(msg)

summaries=[];certs=[]
for n in range(3,9):
 d=json.loads((OLD/('n8' if n==8 else '')/f'census_n{n}.json').read_text());c,sig=schema(n);mask=(1<<len(c))-1
 def refl(t):return sum(((t>>k)&1)<<j for k,j in enumerate(sig))
 def edge(t):return min(t,mask^t,refl(t),mask^refl(t))
 st={'n':n,'eligible_nodes':0,'eligible_incident_edges':0,'literal_four_walks':0,'minimum_hub_triangle_change':None,'strong_nonself_nodes':0,'minimum_descending_edges_nonself':None,'self_converse_single_descent_nodes':[]}
 down=[0]*len(d['nodes'])
 for u,v,m in d['paired_black_edges']:
  if d['nodes'][u]['c3']>d['nodes'][v]['c3']:down[u]+=m
  if d['nodes'][v]['c3']>d['nodes'][u]['c3']:down[v]+=m
 for row in d['nodes']:
  a=decode(row['representative'],n,c);comps=components(a);degree=row['reflection_pairs'];u=row['node']
  if len(comps)==1:
   if not row['SC']:
    st['strong_nonself_nodes']+=1;require(down[u]>=2,('double descent',n,u));st['minimum_descending_edges_nonself']=down[u] if st['minimum_descending_edges_nonself'] is None else min(st['minimum_descending_edges_nonself'],down[u])
   elif down[u]==1:st['self_converse_single_descent_nodes'].append({'node':u,'representative':row['representative'],'degree':degree,'c3':row['c3']})
  if degree<2 or len(comps)<2:continue
  k=0
  while k<len(comps) and len(comps[k]) in (1,3):k+=1
  ell=len(comps)
  while ell>0 and len(comps[ell-1]) in (1,3):ell-=1
  if k==0 or ell==len(comps) or k==len(comps):continue
  left=sum(map(len,comps[:k]));right=sum(map(len,comps[ell:]))
  if max(left,right)<2:continue
  if left<2:a=converse(a)
  groups={}
  for p in paths(a):
   t=encode(a,p,c)
   if t==refl(t):continue
   groups.setdefault((p[0],p[1],p[-1]),{}).setdefault(edge(t),p)
  triple,g=max(groups.items(),key=lambda kv:len(kv[1]));require(len(g)==degree,('normalization',n,u,len(g),degree))
  st['eligible_nodes']+=1;st['eligible_incident_edges']+=degree
  aa,b,z=triple;ds=[r.bit_count() for r in a];delta=2*ds[aa]-ds[b]-ds[z]-3
  require(ds[aa]+ds[b]>n and delta>0,('score guard',n,u,ds,triple))
  st['minimum_hub_triangle_change']=delta if st['minimum_hub_triangle_change'] is None else min(st['minimum_hub_triangle_change'],delta)
  checks={}
  for e,p in g.items():
   s=partner(a,p);q=p[1:]+p[:1]
   require(all(s[x]>>y&1 for x,y in zip(q,q[1:])),('rotation path',n,u,e))
   tq=encode(s,q,c);require(tq!=refl(tq),('rotation black',n,u,e))
   v=partner(s,q);require(v==flip(a,[(aa,b),(aa,z)]),('common literal hub',n,u,e));require(c3(v)-c3(a)==delta,('triangle change',n,u,e))
   require(d['pi'][encode(a,p,c)]==u,('source fibre',n,u,e));require(d['pi'][mask^encode(a,p,c)]!=u,('strong partner',n,u,e))
   checks[e]=(p,q,edge(tq),d['pi'][mask^tq])
  require(len({r[3] for r in checks.values()})==1,('same merged hub',n,u))
  hub=next(iter(checks.values()))[3];require(hub!=u,('nonfolding hub',n,u))
  keys=list(g)
  for e in keys:
   f=next(x for x in keys if x!=e);p,q,r,_=checks[e];pp,qq,rr,_=checks[f]
   require(e not in (f,r,rr),('target avoidance',n,u,e,f,r,rr));st['literal_four_walks']+=1
   certs.append({'n':n,'source_node':u,'source_representative':row['representative'],'source_conversed':left<2,'target_edge':e,'other_edge':f,'first_path':p,'other_path':pp,'rotation_edges':[r,rr],'hub_node':hub,'hub_c3_change':delta})
 summaries.append(st);print(st,flush=True)
require({'status':'all literal checks passed','orders':summaries}==json.loads((ROOT/'controls'/'boundary_hub_checks.json').read_text()),'Bundled full-atlas summary mismatch')
require(json.loads(json.dumps(certs))==json.loads((ROOT/'controls'/'boundary_hub_certificates.json').read_text()),'Bundled certificate mismatch')
print('Matched bundled full-atlas replay result.',flush=True)
