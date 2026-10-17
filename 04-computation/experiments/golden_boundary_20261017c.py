#!/usr/bin/env python3
"""Exact six-colour lifts; finite residue graphs are upper-bound carriers only."""
import argparse, hashlib, json
from itertools import combinations
from pathlib import Path

P=19
TABLE=(4,3,2,4,3,2,1,0,4,1,0,4,3,2,1,0,2,1,0)
PORT_EDGES=((0,1),(1,2),(2,3),(3,4),(4,0),(0,5),(5,6),(6,7),(7,0))
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
def patterns(n,edges):
 out=[]
 def walk(a):
  if len(a)==n:out.append(tuple(a));return
  for c in range(min(5,max(a,default=-1)+1)+1):
   b=a+[c]
   if all(b[i]!=b[j] for i,j in edges if max(i,j)<len(b)):walk(b)
 walk([]);return out
def q(z):return (z[0]*z[0]+z[1]*z[1])%P
def sub(a,b):return ((a[0]-b[0])%P,(a[1]-b[1])%P)
def residue_graph():
 points=[(x,y) for x in range(P) for y in range(P)]
 directions=[z for z in points if q(z)==1]
 edges=set()
 for a,(x,y) in enumerate(points):
  for dx,dy in directions:
   b=((x+dx)%P)*P+(y+dy)%P
   edges.add(tuple(sorted((a,b))))
 need(len(directions)==20 and len(edges)==3610,'complete residue unit graph')
 need(all(TABLE[(points[a][0]+points[a][1])%P]!=TABLE[(points[b][0]+points[b][1])%P] for a,b in edges),'all edges replay inherited five-colouring')
 return points,sorted(edges)
def template(vals,pattern):
 # Recolour literal marked points assigned the fresh colour, not whole fibres.
 for fresh in range(6):
  forward={};backward={}
  for v,c in zip(vals,pattern):
   if c==fresh:continue
   if v in forward and forward[v]!=c or c in backward and backward[c]!=v:break
   forward[v]=c;backward[c]=v
  else:
   unused=[c for c in range(6) if c!=fresh and c not in backward]
   for v in range(5):
    if v not in forward:forward[v]=unused.pop(0)
   return dict(fresh=fresh,permutation=[forward[v] for v in range(5)])
 return None
def projection_params():
 return [(a,b,t) for a in range(P) for b in range(P) if (a*a+b*b)%P==2 for t in range(P)]
def base_colours(points,param):
 a,b,t=param;return [TABLE[(a*x+b*y+t)%P] for x,y in points]
def small_lift(points,edges,pat,param,model):
 vals=base_colours(points,param);perm=model['permutation'];fresh=model['fresh']
 need(len(set(perm+[fresh]))==6,'five-colour permutation and unused sixth')
 need(all(c==fresh or perm[v]==c for v,c in zip(vals,pat)),'every marked colour is obtained by the declared single-point repair')
 need(all(pat[a]!=pat[b] for a,b in edges),'recoloured literal points form an independent colour class')
def branches():
 for ar,ai in [(2,4),(2,15),(7,3),(7,16)]:
  phi=(1+2*ar)%P;root5=(2*phi-1)%P
  for sign in (-1,1):
   jy=7*sign*root5*pow(8,-1,P)%P
   pts=[(0,0),(1,0),((1+ar)%P,ai),((1+ar+ar*ar-ai*ai)%P,(ai+2*ar*ai)%P),((-ar)%P,ai),(17,-jy%P),(15,0),(17,jy)]
   need(all(q(sub(pts[a],pts[b]))==1 for a,b in PORT_EDGES),'every literal old edge reduces to a residue edge')
   yield dict(alpha=[ar,ai],sqrt11=7*sign,phi=phi,points=pts)
def circle_colourings():
 D={1,2,4,6,9,10,13,15,17,18};by=[[] for _ in range(P)]
 for k in range(1,5):
  for A in combinations(range(P),k):
   if all((a-b)%P not in D for a,b in combinations(A,2)):
    mask=sum(1<<a for a in A)
    for a in A:by[a].append(mask)
 need(not any(all((a-b)%P not in D for a,b in combinations(A,2)) for A in combinations(range(P),5)),'no independent five-set, so every independent class has size at most four')
 sols=[]
 def rec(mask,blocks):
  if not mask:sols.append(blocks);return
  left=5-len(blocks)
  if left<=0 or mask.bit_count()>4*left:return
  a=(mask&-mask).bit_length()-1
  for z in by[a]:
   if z&mask==z:rec(mask^z,blocks+[z])
 rec((1<<P)-1,[])
 def partition(vals):return sorted(sum(1<<i for i,c in enumerate(vals) if c==v) for v in set(vals))
 known={tuple(partition([TABLE[(i+t)%P] for i in range(P)])) for t in range(P)}
 need({tuple(sorted(s)) for s in sols}==known,'all five-colour partitions of this circulant are the 19 translates')
 return len(sols)
def run():
 points,edges=residue_graph();params=projection_params()
 # A norm-two row is (1,1) after a norm-one rotation; test all directions independently.
 directions=[z for z in points if q(z)==1]
 for a,b,t in params:
  need(all((a*x+b*y)%P in {1,2,4,6,9,10,13,15,17,18} for x,y in directions),'each projection has exactly the inherited forbidden differences')
 internal=[(0,0),(1,0),(17,15),(15,0),(17,4)]
 internal_edges=((0,1),(0,2),(2,3),(3,4),(0,4))
 small=patterns(5,internal_edges);small_models=[]
 for pat in small:
  for param in params:
   m=template(base_colours(internal,param),pat)
   if m is not None:
    small_lift(internal,internal_edges,pat,param,m)
    small_models.append(dict(pattern=pat,projection=param,**m));break
  else:raise RuntimeError('internal pattern not extended')
 need(len(small)==12,'all 12 canonical internal patterns')
 old=patterns(8,PORT_EDGES);need(len(old)==570,'all 570 original six-colour patterns')
 # Replay the internal construction after relabelling for every original pattern.
 mixed=[]
 for pat in old:
  pp=[pat[i] for i in (0,1,5,6,7)]
  found=None
  for param in params:
   m=template(base_colours(internal,param),pp)
   if m is not None:found=(param,m);break
  need(found is not None,'every old boundary has a literal internal repair')
  param,m=found;small_lift(internal,internal_edges,pp,param,m)
  gold=base_colours([(5,0),(15,0)],param)
  cols=[m['permutation'][x] for x in gold]
  for j,target in enumerate((pat[2],pat[4])):
   if cols[j]==target:cols[j]=m['fresh']
   need(cols[j]!=target,'two golden neighbours repaired separately')
  mixed.append(dict(pattern=pat,projection=param,**m,golden_colours=cols))
 models={};bs=list(branches())
 for b in bs:
  for param in params:
   vals=base_colours(b['points'],param)
   for pat in old:
    if pat in models:continue
    m=template(vals,pat)
    if m is not None:
     small_lift(b['points'],PORT_EDGES,pat,param,m)
     models[pat]=dict(pattern=pat,branch=bs.index(b),projection=param,**m)
 missing=[p for p in old if p not in models]
 return dict(status='FINITE-EXACT universal colouring templates; field and incidence hypotheses proved in companion note',
             residue_vertices=len(points),residue_edges=len(edges),circle_partitions=circle_colourings(),branches=bs,
             internal_templates=small_models,mixed_templates=mixed,cyclotomic_templates=list(models.values()),
             remaining_template_patterns=missing,checks=checks)
if __name__=='__main__':
 print(json.dumps(run(),sort_keys=True,separators=(',',':')))
