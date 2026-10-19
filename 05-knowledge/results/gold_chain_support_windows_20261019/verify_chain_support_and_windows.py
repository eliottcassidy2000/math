#!/usr/bin/env python3
"""Exact, bounded controls for chain-support/rank-window theorems.
The all-size claims are proved in the note, not extrapolated from these tests.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
from collections import defaultdict,deque
from math import comb
import json,hashlib
HERE=Path(__file__).resolve().parent
checks=0
def need(ok,msg):
 global checks
 if not ok:raise RuntimeError(msg)
 checks+=1

def close(p):
 p=list(p)
 for k in range(len(p)):
  for v in range(len(p)):
   if p[v]>>k&1:p[v]|=p[k]
 return None if any(p[v]>>v&1 for v in range(len(p)))else tuple(p)

def add(p,*rels):
 q=list(p)
 for a,b in rels:q[b]|=1<<a
 return close(q)

@lru_cache(None)
def count(p):
 if p is None:return 0
 full=(1<<len(p))-1
 @lru_cache(None)
 def f(s):
  if s==full:return 1
  return sum(f(s|1<<v)for v in range(len(p))if not s>>v&1 and not p[v]&~s)
 return f(0)

def words(p):
 out=[]
 def f(s,w):
  if len(w)==len(p):out.append(tuple(w));return
  for v in range(len(p)):
   if not s>>v&1 and not p[v]&~s:f(s|1<<v,w+[v])
 f(0,[]);return out

def upsets(p):return tuple(sum(1<<v for v in range(len(p))if p[v]>>u&1)for u in range(len(p)))
def degrees(p):return [len(p)-1-(a|b).bit_count()for a,b in zip(p,upsets(p))]
def chains(m,n):return tuple((1<<i)-1 for i in range(m))+tuple(((1<<j)-1)<<m for j in range(n))
def ordinal(p,q):return tuple(p)+tuple(((1<<len(p))-1)|(z<<len(p))for z in q)
def chain(n):return tuple((1<<i)-1 for i in range(n))
def root(t):
 p=chains(t,t);rels=[]
 for i in range(t-2):rels.extend([(i,t+i+2),(t+i,i+2)])
 rels.extend([(0,t),(t-1,2*t-1)])
 return add(p,*rels)
def before(w,a,b):return w.index(a)<w.index(b)
def distance(a,b):return sum(x!=y for x,y in zip(a,b))
def induced(p,mask):
 vs=[v for v in range(len(p))if mask>>v&1]
 return tuple(sum(1<<i for i,a in enumerate(vs)if p[v]>>a&1)for v in vs),vs

def maxmatching(adj,nr):
 left=[-1]*len(adj);right=[-1]*nr
 for s in range(len(adj)):
  queue=deque([s]);parents={};seen={s};end=None
  while queue and end is None:
   a=queue.popleft()
   for b in adj[a]:
    if b in parents:continue
    parents[b]=a
    if right[b]<0:end=b;break
    nxt=right[b]
    if nxt not in seen:seen.add(nxt);queue.append(nxt)
  if end is not None:
   b=end
   while b>=0:
    a=parents[b];old=left[a];left[a]=b;right[b]=a;b=old
 return sum(b>=0 for b in left)

def chain_tests():
 certpath=HERE/'chain_support_sharp_certificates.json';cert=json.loads(certpath.read_text());out=[]
 need([z['m']for z in cert]==list(range(2,8)),'exact chain certificate scope')
 for row in cert:
  m=row['m'];p=chains(m,m);W=words(p);A=[w for w in W if w[0]==0];B=[w for w in W if w[0]==m]
  C=comb(2*m,m);need(len(W)==C==count(p),'chain count independent checks')
  need(W==sorted(W),'canonical chain word indexing')
  pairs=row['matching'];need(len(pairs)==len(A)==C//2,'perfect chain matching size')
  need({i for i,j in pairs}==set(range(len(A)))and {j for i,j in pairs}==set(range(len(B))),'perfect chain matching capacities')
  need(all(distance(A[i],B[j])<=m+2 for i,j in pairs),'sharp-support chain edges')
  U=tuple(range(2*m));V=tuple(range(m-1))+(m,m-1)+tuple(range(m+1,2*m));Z=(m,)+tuple(range(m))+tuple(range(m+1,2*m))
  need([w for w in B if distance(U,w)<=m+1]==[Z],'U unique neighbor')
  need([w for w in B if distance(V,w)<=m+1]==[Z],'V unique neighbor')
  for w in A:
   r=w.index(m);z=(m,)+w[:r]+w[r+1:]
   need(z in set(B)and distance(w,z)==r+1,'run rotation legal and sharp constructed cost')
  out.append({'m':m,'words':C,'perfect_matching':len(pairs),'support':m+2,'lower_bound':m+2})
 for m in range(1,5):
  for n in range(1,5):
   W=words(chains(m,n));A=[w for w in W if w[0]==0];B=[w for w in W if w[0]==m]
   adj=[[j for j,b in enumerate(B)if distance(a,b)<=3]for a in A]
   need(maxmatching(adj,len(B))==comb(m+n-2,m-1),'exact support-three chain formula')
 return out,hashlib.sha256(certpath.read_bytes()).hexdigest()

def root_tests():
 out=[]
 for t in range(2,15):
  p=root(t);h=3**(t-2);g=(3*h+1)//2
  need(count(p)==g and count(add(p,(t,1)))==h,'root counts by ideal DP')
  need(max(degrees(p))<=3,'root degree bound')
  if t>=4:need(3*g<5*h and 3*h<2*g,'root interior heavy range')
  out.append({'t':t,'g':g,'h':h,'max_degree':max(degrees(p))})
 cases=[]
 for m in range(2,8):
  C=comb(2*m,m);t=2;h=1
  while h<=C//2:t+=1;h*=3
  g=(3*h+1)//2;P=ordinal(root(t),chains(m,m));R=add(P,(t,1));second=(2*t,2*t+m)
  E,b=count(P),count(R);c=count(add(R,second))
  need((E,b,c)==(g*C,h*C,h*C//2),'full ordinal counts by independent ideal DP')
  need(max(degrees(P))==max(3,m),'full construction degree')
  cap=h*(C//2-1);threshold=2*b-E
  need(threshold-cap==h-C//2>0,'strict support barrier consequence')
  need(E-b-c==C//2>0,'strict heavy Gold margin')
  # Both small-branch cells are equal; independently check by two comparisons.
  light=add(P,(1,t));csmall=count(add(light,second))
  need(2*csmall==E-b and E-(E-b)-csmall>0,'strict smaller-branch fixed-second margin')
  cases.append({'m':m,'t':t,'vertices':len(P),'degree':max(degrees(P)),'E':E,'b':b,'capacity_required':threshold,'support_m_plus_1_cap':cap,'shortfall':threshold-cap,'actual_and_support_m_plus_2_margin':C//2,'exact_minimum_support':m+2})
 return out,cases

def window_test(p,x,y):
 U=upsets(p);n=len(p);W=words(p);l=min(p[x].bit_count(),p[y].bit_count())+1;r=n-min(U[x].bit_count(),U[y].bit_count());w=r-l+1
 inc=lambda v:((1<<n)-1)^((1<<v)|p[v]|U[v])
 need(w<=(inc(x)|inc(y)).bit_count(),'sharp union-neighborhood window bound')
 fibers=defaultdict(list)
 for z in W:
  need(l<=z.index(x)+1<=r and l<=z.index(y)+1<=r,'canonical rank window contains pair')
  fibers[(z[:l-1],z[r:])].append(z)
 raw=0;signed=[]
 for (pre,post),ws in fibers.items():
  mask=sum(1<<v for v in ws[0][l-1:r]);q,vs=induced(p,mask);inside=words(q)
  need(len(ws)==len(inside),'fiber is full induced extension set')
  need({z[l-1:r]for z in ws}=={tuple(vs[i]for i in z)for z in inside},'fiber bijection full words')
  a=sum(before(z,x,y)for z in ws);bb=len(ws)-a;raw+=min(a,bb);signed.append(a-bb)
 ideals=[s for s in range(1<<n)if all(not s>>v&1 or not p[v]&~s for v in range(n))]
 Ilist=[s for s in ideals if s.bit_count()==l-1];Jlist=[s for s in ideals if s.bit_count()==r];ideal_sum=0;total=0
 for I in Ilist:
  for J in Jlist:
   if I&~J:continue
   S=J^I;O=((1<<n)-1)^J;q,vs=induced(p,S)
   need(x in vs and y in vs,'all ideal windows contain pair')
   a=count(add(q,(vs.index(x),vs.index(y))));bb=count(q)-a
   mult=count(induced(p,I)[0])*count(induced(p,O)[0]);ideal_sum+=mult*min(a,bb);total+=mult*count(q)
 need(ideal_sum==raw and total==len(W),'independent nested-ideal capacity formula')
 a=sum(before(z,x,y)for z in W);bb=len(W)-a
 need(2*raw==len(W)-sum(abs(z)for z in signed),'window signed-defect identity')
 if all(z>=0 for z in signed)or all(z<=0 for z in signed):need(raw==min(a,bb),'common-sign saturation')
 return {'n':n,'pair':[x,y],'window':[l,r],'fibers':len(fibers),'capacity':raw,'global_smaller_side':min(a,bb),'mixed_bias':min(signed)<0<max(signed)}

def main():
 cert,sha=chain_tests();roots,cases=root_tests();window_records=[]
 for p in[chains(2,2),chains(3,3),root(4),ordinal(chain(2),ordinal(chains(2,2),chain(2)))]:
  for x,y in combinations(range(len(p)),2):
   if (p[x]>>y|p[y]>>x)&1:continue
   window_records.append(window_test(p,x,y))
 need(any(z['n']==4 and z['pair']==[0,2]and z['capacity']==2 and z['global_smaller_side']==3 and z['mixed_bias']for z in window_records),'window sign-cancellation hostile')
 return {'status':'PASS','scope':'General claims proved in readable note; exact sharp chain repairs only m=2,...,7; finite root and rank-window controls','sharp_chain_certificates':cert,'certificate_sha256':sha,'root_count_controls':roots,'strict_Gold_support_threshold_cases':cases,'window_controls':window_records,'checks':checks}

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
