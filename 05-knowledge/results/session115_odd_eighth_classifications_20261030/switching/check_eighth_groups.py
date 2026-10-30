#!/usr/bin/env python3
from itertools import permutations,product
from fractions import Fraction as F
from collections import Counter
import json
checks=0
def gate(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def comp(a,b):return tuple(a[x] for x in b)
def inv(a):return tuple(a.index(i) for i in range(len(a)))
def order(a):
 seen=set();ans=1
 from math import lcm
 for i in range(len(a)):
  if i in seen:continue
  j=i;n=0
  while j not in seen:seen.add(j);n+=1;j=a[j]
  ans=lcm(ans,n)
 return ans
def closure(gens,n):
 I=tuple(range(n)); S={I};todo=[I]
 while todo:
  a=todo.pop()
  for b in gens:
   c=comp(a,b)
   if c not in S:S.add(c);todo.append(c)
 return frozenset(S)
def parity(a):return sum(a[i]>a[j] for i in range(len(a)) for j in range(i+1,len(a)))%2
I8=tuple(range(8))
GL=[]
for cols in product(range(1,8),repeat=3):
 a=tuple((cols[0] if x&1 else 0)^(cols[1] if x&2 else 0)^(cols[2] if x&4 else 0) for x in range(8))
 if len(set(a))==8:GL.append(a)
gate(len(GL)==168,'GL size')
# Exhaust all odd subgroups by adjoining each group element; every prefix subgroup of an odd subgroup stays odd.
odd={frozenset([I8])};todo=list(odd)
while todo:
 H=todo.pop()
 for g in GL:
  if g in H:continue
  J=closure(list(H)+[g],8)
  if len(J)%2 and J not in odd:odd.add(J);todo.append(J)
gate(Counter(map(len,odd))=={1:1,3:28,7:8,21:8},'odd subgroup census')
rows=[]
for H in sorted(odd,key=lambda h:(len(h),sorted(h))):
 G=[tuple(h[x]^v for x in range(8)) for h in H for v in range(8)]
 direct=F(sum(order(g)%2==0 for g in G),len(G))
 fix=Counter(sum(h[x]==x for x in range(8)) for h in H)
 formula=1-sum((F(1,k)*count for k,count in fix.items()),F(0))/len(H)
 gate(direct==formula,('affine formula',len(H)))
 gate(direct=={1:F(7,8),3:F(5,8),7:F(1,8),21:F(3,8)}[len(H)],'eighth spectrum')
 weights=[];unseen=set(range(1,8));orbit_sizes=[]
 while unseen:
  v=min(unseen);orb={h[v] for h in H};unseen-=orb;orbit_sizes.append(len(orb))
  central=[g for g in G if all(g[x^v]==(g[x]^v) for x in range(8))]
  weights.append(F(sum(order(g)%2 for g in central),len(central)))
 expected={1:[F(1,8)]*7,3:[F(3,8),F(1,8),F(1,8)],7:[F(1,8)],21:[F(3,8)]}[len(H)]
 gate(sorted(weights)==sorted(expected),('affine decoration spectrum',len(H)))
 gate(sum(weights,F(0))==direct,('affine weighted identity',len(H)))
 rows.append({'H_order':len(H),'even_mass':str(direct),'fixed_vector_histogram':dict(sorted(fix.items())),'translation_orbit_sizes':sorted(orbit_sizes),'decoration_weights':sorted(map(str,weights))})
# Independent finite checks of arbitrary-Sylow-D8 theorem.
def find_D8(G):
 I=tuple(range(len(G[0])))
 for r in G:
  if order(r)!=4:continue
  ri=inv(r)
  for s in G:
   if order(s)==2 and comp(comp(s,r),s)==ri:
    P=closure([r,s],len(r))
    if len(P)==8:return P
 raise RuntimeError('No D8')
def audit_group(name,G):
 P=find_D8(G); gate(len(G)%8==0 and len(G)%16!=0,('Sylow',name))
 Os={g:order(g) for g in G}; remaining={g for g in P if Os[g] in (2,4)};decor=[]
 while remaining:
  x=min(remaining);cx={g for g in G if comp(g,x)==comp(x,g)}
  cl={comp(comp(g,x),inv(g)) for g in G}; remaining-=cl
  rho=F(sum(Os[g]%2 for g in cx),len(cx))
  decor.append((Os[x],str(rho),len(cx)))
 c=sum(o==2 for o,w,z in decor); direct=F(sum(o%2==0 for o in Os.values()),len(G))
 gate(direct==F(2*c+1,8),(name,'D8 formula',direct,c))
 gate(sum((F(w) for _,w,_ in decor),F(0))==direct,(name,'weighted identity'))
 gate(sorted(w for o,w,z in decor if o==4)==['1/4'],(name,'order4'))
 gate(sorted(w for o,w,z in decor if o==2)==sorted(['1/8']+['1/4']*(c-1)),(name,'involutions'))
 return {'name':name,'order':len(G),'involution_classes':c,'even_mass':str(direct),'decorations':sorted(decor)}
S4=list(permutations(range(4)));D8=list(closure([(1,2,3,0),(0,3,2,1)],4))
GL7=[tuple(g[x+1]-1 for x in range(7)) for g in GL]
examples=[('D8',D8),('S4',S4),('GL3(2)',GL7)]
for n in (6,7):examples.append(('A'+str(n),[a for a in permutations(range(n)) if parity(a)==0]))
for name,G in examples[:3]:
 n=len(G[0]);GG=[tuple(g)+tuple(n+(j+k)%3 for j in range(3)) for g in G for k in range(3)]
 examples.append((name+' x C3',GG))
controls=[audit_group(name,G) for name,G in examples]
print(json.dumps({'status':'PASS','checks':checks,'odd_subgroup_census':dict(sorted(Counter(map(len,odd)).items())),'affine_representatives':[next(row for row in rows if row['H_order']==h) for h in(1,3,7,21)],'D8_controls':controls},indent=2))
