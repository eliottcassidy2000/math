#!/usr/bin/env python3
"""Exact fold/inflate groupoid checks; standard library only."""
from itertools import combinations,permutations
from functools import lru_cache
import json

def req(c,s):
 if not c:raise RuntimeError(s)
@lru_cache(None)
def pairs(n):return tuple(combinations(range(n),2))
@lru_cache(None)
def perms(n):return tuple(permutations(range(n)))
def matrix(n,k):
 a=[[0]*n for _ in range(n)]
 for j,(u,v) in enumerate(pairs(n)):a[u][v]=(k>>j)&1;a[v][u]=1-a[u][v]
 return a
def encode(a,p=None):
 if p is None:p=tuple(range(len(a)))
 return sum(a[p[u]][p[v]]<<j for j,(u,v) in enumerate(pairs(len(a))))
def canon(a):return min(encode(a,p) for p in perms(len(a)))
def anti(a,r):return all(r[r[i]]==i for i in range(len(a))) and all(a[i][j]==a[r[j]][r[i]] for i,j in pairs(len(a)))
def put(a,u,v,b):a[u][v]=b;a[v][u]=1-b

def inflate(m,x,y,color=None):
 a=matrix(m,x);b=matrix(m,y);n=2*m+int(color is not None);t=[[0]*n for _ in range(n)]
 for i in range(m):put(t,2*i,2*i+1,1)
 for i,j in pairs(m):
  put(t,2*i,2*j,a[i][j]);put(t,2*i+1,2*j+1,1-a[i][j]);g=a[i][j]^b[i][j]
  put(t,2*i,2*j+1,g);put(t,2*j,2*i+1,g)
 if color is not None:
  for i in range(m):put(t,2*m,2*i,(color>>i)&1);put(t,2*i+1,2*m,(color>>i)&1)
 rho=tuple([i^1 for i in range(2*m)]+([2*m] if color is not None else []))
 req(anti(t,rho),'inflation not anti-involutive');return t,rho

def fold(t,rho):
 req(anti(t,rho),'bad fold input');fixed=[i for i in range(len(t)) if rho[i]==i];req(len(fixed)<=1,'too many centers')
 uv=[]
 for i in range(len(t)):
  if i<rho[i]:uv.append((i,rho[i]) if t[i][rho[i]] else(rho[i],i))
 m=len(uv);x=y=0
 for k,(i,j) in enumerate(pairs(m)):
  u,v=uv[i];w,z=uv[j];aa=t[u][w];g=t[u][z]
  req(g==t[w][v],'cross bit not symmetric');x|=aa<<k;y|=(aa^g)<<k
 c=sum(t[fixed[0]][u]<<i for i,(u,v) in enumerate(uv)) if fixed else None
 return m,x,y,c

def pair_transform(m,x,y,c,p):
 xx=encode(matrix(m,x),p);yy=encode(matrix(m,y),p)
 cc=sum(((c>>p[i])&1)<<i for i in range(m)) if c is not None else None
 return xx,yy,cc

def pair_canon(m,x,y,c):return min(pair_transform(m,x,y,c,p) for p in perms(m))
def pair_aut(m,x,y,c):return sum(pair_transform(m,x,y,c,p)==(x,y,c) for p in perms(m))

def labelled_orbits(n):
 covered={};representatives=[]
 for k in range(1<<len(pairs(n))):
  if k in covered:continue
  a=matrix(n,k);orb={encode(a,p) for p in perms(n)};root=min(orb);representatives.append(root)
  covered.update((h,root) for h in orb)
 return covered,representatives

def main():
 small=[];checks=0
 for n in range(7):
  orbit,roots=labelled_orbits(n);full=(1<<len(pairs(n)))-1;rows=[]
  for k in roots:
   if orbit[full^k]!=k:continue
   t=matrix(n,k);autos=[p for p in perms(n) if encode(t,p)==k];antis=[p for p in perms(n) if anti(t,p)]
   req(antis,'SC class lacks anti-involution');keys=set()
   for r in antis:
    m,x,y,c=fold(t,r);keys.add(pair_canon(m,x,y,c));tt,rr=inflate(m,x,y,c)
    req(canon(tt)==k,'fold/inflate changed class');central=sum(all(p[r[i]]==r[p[i]] for i in range(n)) for p in autos)
    req(central==pair_aut(m,x,y,c),'marked automorphism mismatch');req(len(autos)==len(antis)*central,'conjugacy orbit size')
    checks+=1
   req(len(keys)==1,'fold class depends on anti-involution');rows.append({'class':k,'automorphisms':len(autos),'anti_involutions':len(antis),'fold_class':next(iter(keys)),'fold_automorphisms':pair_aut(*fold(t,antis[0]))})
  req(len(rows)==[1,1,1,2,2,8,12][n],'SC small census');small.append({'n':n,'ordinary_classes':len(roots),'self_converse':len(rows),'rows':rows})
 pair_counts=[]
 for m in range(5):
  bits=len(pairs(m));lim=1<<bits
  for colored in (False,True):
   seen=set();reps=[]
   for x in range(lim):
    for y in range(lim):
     for c in range(1<<m if colored else 1):
      cc=c if colored else None;key=(x,y,cc)
      if key in seen:continue
      orb={pair_transform(m,x,y,cc,p) for p in perms(m)};seen.update(orb);reps.append(key)
      t,r=inflate(m,x,y,cc);req(fold(t,r)==(m,x,y,cc),'literal inverse failure');checks+=1
   pair_counts.append({'m':m,'colored':colored,'labelled_pairs':len(seen),'isomorphism_classes':len(reps)})
 req([x['isomorphism_classes'] for x in pair_counts]==[1,1,1,2,2,8,12,88,176,2752],'pair count bridge')
 c3=next(x for x in small[3]['rows'] if x['automorphisms']==3)
 req(c3['fold_automorphisms']==1,'C3 automorphism hostile absent')
 if '--test-failing-gate' in __import__('sys').argv:req(False,'intentional failure control')
 print(json.dumps({'status':'PASS','fold_controls':checks,'small_tournament_census':small,'ordered_pair_counts':pair_counts,'C3_stabilizer_hostile':c3},sort_keys=True,indent=2))
if __name__=='__main__':main()
