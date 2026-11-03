#!/usr/bin/env python3
"""Independent projective puzzle audit. No producer imports. Explicit checks.
Lines use dual linear equations, words act on a board of labeled tokens,
and the subgroup is enumerated by generator closure rather than normal forms.
"""
from itertools import product,combinations
from collections import deque
from fractions import Fraction
import argparse,json
p=argparse.ArgumentParser();p.add_argument('--mutate',choices=['single-swap','inverse-h','reciprocal-adjacency']);args=p.parse_args()
checks=0
def need(x,s):
 global checks
 checks+=1
 if not x:raise RuntimeError(s)
P=[(1,a,b) for a in range(3) for b in range(3)]+[(0,1,a) for a in range(3)]+[(0,0,1)]
need(len(P)==13 and len(set(P))==13,'13 normalized points')
classes={frozenset((tuple(v),tuple(2*x%3 for x in v))) for v in product(range(3),repeat=3) if any(v)}
need(len(classes)==13,'13 scalar classes')
need(all(sum(v in cl for v in P)==1 for cl in classes),'one normalized representative per projective point')
# Dual homogeneous equations: each normalized covector defines a line.
lines={frozenset(i for i,v in enumerate(P) if sum(a*b for a,b in zip(n,v))%3==0) for n in P}
need(len(lines)==13 and all(len(L)==4 for L in lines),'13 four-point projective lines')
for a,b in combinations(range(13),2):need(sum({a,b}<=L for L in lines)==1,'unique line through two points')
for L,M in combinations(lines,2):need(len(L&M)==1,'two projective lines meet once')
I=tuple(range(13))
def compose(a,b):return tuple(a[b[i]] for i in range(13))
def inv(s):return tuple(s.index(i) for i in range(13))
def power(s,k):
 z=I
 for _ in range(k):z=compose(s,z)
 return z

def run(path):
 board=list(range(13));steps=[]
 need(path[0]==path[-1]==0,'closed based word')
 for a,b in zip(path,path[1:]):
  need(board[a]==0 and a!=b,'hole at declared start of legal move')
  L=next(L for L in lines if {a,b}<=L);c,d=sorted(L-{a,b})
  board[a],board[b]=board[b],board[a]
  if not(args.mutate=='single-swap' and (a,b)==(0,1)):board[c],board[d]=board[d],board[c]
  need(board[b]==0,'hole at declared endpoint')
  steps.append({'from':a,'to':b,'line':sorted(L),'other_swap':[c,d]})
 # board[pos]=original token; invert to get position action on original labels.
 return tuple(board.index(token) for token in range(13)),steps

g,gs=run([0,1,3,2,4,0]);h,hs=run([0,1,7,12,9,4,0])
def cycperm(*cycles):
 s=list(I)
 for c in cycles:
  for a,b in zip(c,c[1:]+c[:1]):s[a]=b
 return tuple(s)
need(g==cycperm((1,4,2,12,3),(6,11,10,7,8)),'g pinned cycle action')
need(h==cycperm((1,4,12,2),(6,10,11,8)),'h pinned cycle action')
need(power(g,5)==I and g!=I,'g has order five')
need(power(h,4)==I and power(h,2)!=I,'h has order four')
if args.mutate=='inverse-h':h=inv(h)
need(compose(compose(h,g),inv(h))==power(g,2),'global conjugation exponent two')
V=(3,1,4,2,12)
need(all(g[V[i]]==V[(i+1)%5] for i in range(5)),'orbit translation')
need(all(h[V[i]]==V[2*i%5] for i in range(5)),'orbit multiplier')
# Breadth-first generated subgroup, independent of asserted semidirect forms.
G={I};queue=deque([I])
while queue:
 s=queue.popleft()
 for generator in (g,h,inv(g),inv(h)):
  t=compose(generator,s)
  if t not in G:G.add(t);queue.append(t)
 need(len(G)<=100,'bounded subgroup enumeration guard')
need(len(G)==20,'generated subgroup order twenty')
restrictions={tuple(V.index(s[v]) for v in V) for s in G}
need(len(restrictions)==20,'five-point action faithful')
need(restrictions=={tuple((a*i+b)%5 for i in range(5)) for a in range(1,5) for b in range(5)},'restriction exactly AGL1(5)')
S={1,2} if args.mutate!='reciprocal-adjacency' else {1,4}
E=frozenset((V[i],V[(i+s)%5]) for i in range(5) for s in S)
need(all(((a,b) in E)+((b,a) in E)==1 for a,b in combinations(V,2)),'exactly one directed edge per unordered pair')
need(all(sum(a==u for a,b in E)==2 for u in V),'regular outdegrees')
def act(s,E):return frozenset((s[a],s[b]) for a,b in E)
family={act(s,E) for s in G}
need(len(family)==4,'four-object orbit')
need(act(g,E)==E and act(h,E)!=E,'translation preserves, multiplier changes tournament')
need(act(power(h,2),E)==frozenset((b,a) for a,b in E),'normalizer square gives converse')
need(sum(act(s,E)==E for s in G)==5,'single-object stabilizer is C5')
need(sum(all(act(s,F)==F for F in family) for s in G)==5,'family kernel is C5')
A=[[int((u,v) in E) for v in V] for u in V]
P1=[[int(g[u]==v) for v in V] for u in V]
def mm(A,B):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]
P2=mm(P1,P1)
need(A==[[P1[i][j]+P2[i][j] for j in range(5)] for i in range(5)],'adjacency equals row P+P squared')
M=[[int(i==j) for j in range(5)] for i in range(5)];tr=[];c=[Fraction(1)]
for k in range(1,6):M=mm(M,A);tr.append(sum(M[i][i] for i in range(5)))
for k in range(1,6):c.append(-sum(c[k-i]*tr[i-1] for i in range(1,k+1))/k)
need(c==[1,0,0,-5,-5,-2],'regular polynomial via Newton traces')
L=frozenset([0,1,2,12]);image=frozenset(g[x] for x in L)
need(L in lines and image==frozenset([0,3,4,12]) and image not in lines,'direct noncollineation witness')
need(26*24*18//2==5616 and 5616%5!=0,'PGL order obstruction')
print(json.dumps({'status':'PASS','checks':checks,'points':P,'g_images':g,'h_images':h,'g_steps':gs,'h_steps':hs,
 'orbit':V,'group_order':len(G),'family_size':len(family),'adjacency':A,'charpoly_descending':list(map(int,c)),
 'noncollineation':{'line':sorted(L),'image':sorted(image)},
 'scope':'Explicit projective-plane puzzle witnesses and their chosen five-orbit only; no global M12 action on fixed orbit, projective collineation, Baby Monster embedding or Hecke map.'},indent=2,sort_keys=True))
