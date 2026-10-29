#!/usr/bin/env python3
"""Exact signed modular cycle indices; stdlib only, no imported census."""
from itertools import permutations, product
from functools import lru_cache
from fractions import Fraction
from collections import Counter
import json
from pathlib import Path

def need(condition, message):
 if not condition: raise RuntimeError(message)

def matrix(n,m):
 a=[[0]*n for _ in range(n)]; k=0
 for i in range(n):
  for j in range(i+1,n): a[i][j]=(m>>k)&1; a[j][i]=1-a[i][j]; k+=1
 return tuple(map(tuple,a))
def sub(a,b): return tuple(tuple(a[i][j] for j in b) for i in b)
def dual(a): return tuple(tuple(a[j][i] for j in range(len(a))) for i in range(len(a)))
def cyc(p):
 seen=set(); out=[]
 for i in range(len(p)):
  if i not in seen:
   j=i; c=0
   while j not in seen: seen.add(j); c+=1; j=p[j]
   out.append(c)
 return tuple(sorted(out))
def cycles(p):
 seen=set(); out=[]
 for i in range(len(p)):
  if i not in seen:
   j=i; c=[]
   while j not in seen: seen.add(j); c.append(j); j=p[j]
   out.append(c)
 return out
@lru_cache(None)
def maps(a,b,e=1):
 if len(a)!=len(b): return ()
 n=len(a)
 # Degree filter is independent of modular decomposition.
 da=list(map(sum,a)); db=list(map(sum,b)); ans=[]
 for p in permutations(range(n)):
  if any(db[p[i]] != (da[i] if e==1 else n-1-da[i]) for i in range(n)): continue
  if all(b[p[i]][p[j]]==(a[i][j] if e==1 else 1-a[i][j]) for i in range(n) for j in range(i+1,n)): ans.append(p)
 return tuple(ans)
def brute(a,e):
 d=Counter(cyc(p) for p in maps(a,a,e)); z=len(maps(a,a,1))
 return {k:Fraction(v,z) for k,v in d.items()}
def children(a):
 n=len(a); full=(1<<n)-1; mods=[]
 for s in range(1,full+1):
  inside=[i for i in range(n) if s>>i&1]
  if all(len({a[j][i] for i in inside})==1 for j in range(n) if not s>>j&1): mods.append(s)
 strong=[s for s in mods if all(not(s&t) or s&t==s or s&t==t for t in mods)]
 proper=[s for s in strong if s!=full]
 maximal=[s for s in proper if not any(s!=t and s&t==s for t in proper)]
 return [tuple(i for i in range(n) if s>>i&1) for s in maximal]
def polymul(a,b):
 c={}
 for x,v in a.items():
  for y,w in b.items():
   k=tuple(sorted(x+y)); c[k]=c.get(k,0)+v*w
 return c
@lru_cache(None)
def recursive(a,e):
 if len(a)==1:return {(1,):Fraction(1)}
 bs=children(a); ts=[sub(a,b) for b in bs]; q=sub(a,[b[0] for b in bs]); k=len(bs)
 hp=[p for p in maps(q,q,1) if all(maps(ts[i],ts[p[i]],1) for i in range(k))]
 he=[p for p in maps(q,q,e) if all(maps(ts[i],ts[p[i]],e) for i in range(k))]
 out={}
 for p in he:
  term={():Fraction(1)}
  for c in cycles(p):
   ell=len(c); poly=recursive(ts[c[0]],e**ell)
   term=polymul(term,{tuple(ell*d for d in x):v for x,v in poly.items()})
  for x,v in term.items(): out[x]=out.get(x,0)+v/len(hp)
 return out

def image_mask(s,p,complement=False):
 t=sum(1<<p[i] for i in range(len(p)) if s>>i&1)
 return t^((1<<len(p))-1) if complement else t

def attachment(a):
 n=len(a); g=maps(a,a,1); anti=maps(a,a,-1)
 orbits={min(image_mask(s,p) for p in g) for s in range(1<<n)}
 total=sum(v*2**len(c) for c,v in recursive(a,1).items())
 need(total==len(orbits),'ordinary attachment orbit mismatch')
 if not anti: fixed=0
 else:
  h=anti[0]; fixed=sum(min(image_mask(image_mask(s,h,True),p) for p in g)==s for s in orbits)
 predicted=sum(v*2**len(c) for c,v in recursive(a,-1).items() if all(d%2==0 for d in c))
 need(fixed==predicted,'signed attachment orbit mismatch')
 need(n%2==0 or fixed==0,'odd-base marked converse obstruction failed')
 return int(total),int(fixed)
def compose(q,ts):
 addr=[(i,j) for i,t in enumerate(ts) for j in range(len(t))]
 return tuple(tuple(ts[i][u][v] if i==j else q[i][j] for j,v in addr) for i,u in addr)
def encode(poly): return {','.join(map(str,k)):str(v) for k,v in sorted(poly.items())}
def main():
 checks=0; sc=0; summaries=[]; reps_by_n={}; descent=[]
 for n in range(1,6):
  count=0; reps=[]
  for m in range(1<<(n*(n-1)//2)):
   a=matrix(n,m)
   for e in (1,-1): need(recursive(a,e)==brute(a,e),('signed recursion mismatch',n,m,e))
   attachment(a); count+=1; sc+=bool(maps(a,a,-1))
   if not any(maps(a,b,1) for b in reps): reps.append(a)
   if n%2 and maps(a,a,-1):
    fixed_vertices={i for f in maps(a,a,-1) for i in range(n) if f[i]==i}
    root=min(fixed_vertices)
    need({g[root] for g in maps(a,a,1)}==fixed_vertices,'anti-fixed vertex orbit mismatch')
  reps_by_n[n]=reps
  summaries.append({'n':n,'labelled_tournaments':count,'ordinary_isotypes':len(reps),'self_converse_isotypes':sum(bool(maps(a,a,-1)) for a in reps)}); checks+=count
  if n>=3 and n%2:
   predicted=sum(attachment(a)[1] for a in reps_by_n[n-1])
   actual=sum(bool(maps(a,a,-1)) for a in reps)
   need(predicted==actual,'odd unpointed descent mismatch')
   descent.append({'order':n,'sum_even_base_R':predicted,'direct_self_converse_isotypes':actual})
 c3=matrix(3,5); one=matrix(1,0); edge=matrix(2,1)
 examples={
  'C3':c3,
  'transitive4':matrix(4,63),
  'strong4_C3_edge_single_single':compose(c3,[edge,one,one]),
  'A_C3_to_sink':compose(edge,[c3,one]),
  'six_join_dual':compose(edge,[c3,c3]),
  'C3_edge_edge_edge':compose(c3,[edge,edge,edge]),
  'C3_to_transitive3':compose(edge,[c3,matrix(3,7)]),
  'C3_C3_single_single':compose(c3,[c3,one,one]),
  'eight_dual_blocks':compose(edge,[compose(edge,[c3,one]),compose(edge,[one,c3])]),
 }
 result={}
 for name,a in examples.items():
  for e in (1,-1): need(recursive(a,e)==brute(a,e),('composite recursion mismatch',name,e))
  x,y=attachment(a)
  result[name]={'n':len(a),'aut':len(maps(a,a,1)),'anti':len(maps(a,a,-1)),'Z_plus':encode(recursive(a,1)),'Z_minus':encode(recursive(a,-1)),'attachment_orbits':x,'self_converse_marked_attachment_orbits':y}
 need(recursive(examples['transitive4'],1)==recursive(examples['strong4_C3_edge_single_single'],1) and recursive(examples['transitive4'],-1)==recursive(examples['strong4_C3_edge_single_single'],-1),'same-index hostile failed')
 need(not maps(examples['transitive4'],examples['strong4_C3_edge_single_single'],1),'same-index hostile unexpectedly isomorphic')
 out={'status':'FINITE-EXACT','universe':summaries,'labelled_cases':checks,'self_converse_labelled_cases':sc,'tests':['brute permutations versus strong-module signed recursion','direct binary attachment orbits versus ordinary cycle index','direct converse-fixed attachment orbits versus even-cycle signed index'],'examples':result,'odd_unpointed_descent':descent}
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
