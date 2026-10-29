#!/usr/bin/env python3
"""Exact labelled tournaments n<=5; selected SCC chains; no external libraries."""
import itertools as it, json
from collections import Counter
CHECKS=0
def check(x):
 global CHECKS
 CHECKS+=1
 if not x: raise AssertionError(CHECKS)
def tour(n,mask):
 a=[[0]*n for _ in range(n)]
 for k,(u,v) in enumerate(it.combinations(range(n),2)): a[u][v]=(mask>>k)&1; a[v][u]=1-a[u][v]
 return a
def paths(a):
 n=len(a)
 def rec(p,left):
  if not left: yield tuple(p)
  for v in left:
   if not p or a[p[-1]][v]: yield from rec(p+[v],left-{v})
 return list(rec([],set(range(n))))
def blue(a,p):
 return all(a[p[i]][p[j]]==a[p[-1-j]][p[-1-i]] for i in range(len(p)) for j in range(i+1,len(p)))
def flip(a,edges):
 b=[r[:] for r in a]
 for u,v in edges: b[u][v]^=1; b[v][u]^=1
 return b
def pedges(p): return {tuple(sorted((u,v))) for u,v in zip(p,p[1:])}
def native(a,p):return flip(a,set(it.combinations(range(len(a)),2))-pedges(p))
def switch(a,S):return flip(a,{(u,v) for u,v in it.combinations(range(len(a)),2) if (u in S)!=(v in S)})
def legal(a,p):return all(a[u][v] for u,v in zip(p,p[1:]))
def scc(a):
 n=len(a);r=[row[:] for row in a]
 for i in range(n):r[i][i]=1
 for k in range(n):
  for i in range(n):
   for j in range(n):r[i][j]|=r[i][k]&r[k][j]
 left=set(range(n)); out=[]
 while left:
  v=min(left);b={u for u in left if r[u][v] and r[v][u]};out.append(b);left-=b
 return sorted(out,key=lambda b:-sum(a[next(iter(b))][next(iter(c))] for c in out if c!=b))
def sub(a,B):return [[a[u][v] for v in sorted(B)] for u in sorted(B)]
def word(a,p):return tuple(a[p[i]][p[j]] for i in range(len(a)) for j in range(i+1,len(a)))
def fibres(a):return {word(a,p) for p in paths(a)}
def revword(w,n):
 a=tour(n,0)
 for b,(u,v) in zip(w,it.combinations(range(n),2)):a[u][v]=b;a[v][u]=1-b
 return tuple(a[n-1-j][n-1-i] for i in range(n) for j in range(i+1,n))
def recurrence(a):
 bs=[sub(a,B) for B in scc(a)]; fs=[fibres(b) for b in bs]
 b=1
 for i in range(len(bs)//2):
  j=len(bs)-1-i
  if len(bs[i])!=len(bs[j]) or {revword(w,len(bs[i])) for w in fs[i]}!=fs[j]:return 0
  b*=len(fs[i])
 if len(bs)%2:
  i=len(bs)//2;b*=sum(w==revword(w,len(bs[i])) for w in fs[i])
 return b
census=[]; two_step_switch=None
for n in range(1,6):
 c=Counter()
 for m in range(1<<(n*(n-1)//2)):
  a=tour(n,m);ps=paths(a);fs={word(a,p) for p in ps};bf={word(a,p) for p in ps if blue(a,p)}
  check(len(bf)==recurrence(a));c['tournaments']+=1;c['paths']+=len(ps)
  for p in ps:
   b=native(a,p);check(legal(b,p));check(native(b,p)==a);check(blue(a,p)==blue(b,p))
   # Converse carries the reversed ordered path.
   ac=flip(a,set(it.combinations(range(n),2)))
   check(native(ac,p[::-1])==flip(b,set(it.combinations(range(n),2))))
   for bits in range(1<<(n-1)):
    S={v for v in range(n) if bits>>v&1};sw=switch(a,S)
    check(legal(sw,p)==(not S));check(native(sw,p)==switch(b,S))
    if n>=3:check(sw!=b)
   if n<=4:
    for q in paths(b):
     d=native(b,q);check(d==flip(a,pedges(p)^pedges(q)))
     if d!=a and two_step_switch is None:
      for bits in range(1<<(n-1)):
       S={v for v in range(n) if bits>>v&1}
       if d==switch(a,S):two_step_switch={'n':n,'mask':m,'P':p,'Q':q,'S':sorted(S)}
 census.append({'n':n,**c})
def cycle(n):return [[0 if u==v else int((v-u)%n <= (n-1)//2) for v in range(n)] for u in range(n)]
def join(bs):
 n=sum(map(len,bs));a=[[0]*n for _ in range(n)];off=0
 for b in bs:
  for u in range(len(b)):
   for v in range(len(b)):a[off+u][off+v]=b[u][v]
   for v in range(off+len(b),n):a[off+u][v]=1
  off+=len(b)
 return a
controls=[]
for label,bs in [('C3>C3',[cycle(3)]*2),('C3>1>C3',[cycle(3),[[0]],cycle(3)]),('T5>T5',[cycle(5)]*2),('C3>T5>C3',[cycle(3),cycle(5),cycle(3)])]:
 a=join(bs);ps=paths(a);f={word(a,p) for p in ps};b={word(a,p) for p in ps if blue(a,p)};check(len(b)==recurrence(a))
 # Every nonsingleton old block is destroyed by its incident boundary port.
 for p in ps:
  k=native(a,p);off=0
  for child in bs:
   B=set(range(off,off+len(child)));off+=len(child)
   module=all(len({k[v][u] for u in B})<=1 for v in set(range(len(a)))-B)
   check(module==(len(B)==1))
 controls.append({'name':label,'n':len(a),'H':len(ps),'fibre':len(f),'blue':len(b),'black':len(f)-len(b)})
# Sharp order-seven two-native-step switching witness; all 512 within-part fillings.
P=(0,4,1,5,2,6,3); Q=(1,6,0,5,3,4,2); S={0,1,2,3}
cut={(u,v) for u,v in it.combinations(range(7),2) if (u in S)!=(v in S)}
check(not pedges(P)&pedges(Q));check(pedges(P)|pedges(Q)==cut)
inside=sorted(set(it.combinations(range(7),2))-cut);colour_counts=Counter();witness=None
for bits in range(1<<len(inside)):
 a=tour(7,0)
 for k,(u,v) in enumerate(inside):a[u][v]=(bits>>k)&1;a[v][u]=1-a[u][v]
 for u,v in zip(P,P[1:]):a[u][v]=1;a[v][u]=0
 for u,v in zip(Q,Q[1:]):a[u][v]=0;a[v][u]=1
 b=native(a,P);c=native(b,Q)
 check(legal(a,P));check(legal(b,Q));check(c==switch(a,S));check(c!=a)
 colors=('blue' if blue(a,P) else 'black')+'-'+('blue' if blue(b,Q) else 'black')
 colour_counts[colors]+=1
 if witness is None and colors=='black-black':witness={'adjacency':a,'P':P,'Q':Q,'S':sorted(S)}
# Complete extremal skeleton classification with P fixed in canonical order.
canonical=tuple(range(7)); cuts={frozenset((u,v) for u,v in it.combinations(range(7),2) if (u in SS)!=(v in SS)) for bits in range(1,1<<6) for SS in [{v for v in range(7) if bits>>v&1}]}
skeletons=[q for q in it.permutations(range(7)) if frozenset(pedges(canonical)^pedges(q)) in cuts]
expected=(2,5,0,3,6,1,4)
check(set(skeletons)=={expected,expected[::-1]})
sharp={'n':7,'canonical_Q_skeletons':skeletons,'within_part_fillings':1<<len(inside),'colour_counts':dict(colour_counts),'black_black_witness':witness}
print(json.dumps({'status':'PASS','checks':CHECKS,'universe':'all labelled tournaments n=1..5; every directed Hamilton path; switches modulo complementary subsets; all legal two-step words n<=4; selected chains','sharp_two_step_switch':sharp,'census':census,'two_step_switch_witness':two_step_switch,'controls':controls},indent=2))
