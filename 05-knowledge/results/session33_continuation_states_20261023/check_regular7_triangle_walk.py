#!/usr/bin/env python3
import itertools,json
from collections import Counter
from fractions import Fraction
V=range(7);E=list(itertools.combinations(V,2));index={e:i for i,e in enumerate(E)}
def arc(m,a,b):return bool(m>>index[tuple(sorted((a,b)))]&1)==(a<b)
def hp(m):
 d=[[0]*7 for _ in range(128)]
 for a in V:d[1<<a][a]=1
 for s in range(128):
  for a in V:
   if d[s][a]:
    for b in V:
     if not(s>>b&1) and arc(m,a,b):d[s|1<<b][b]+=d[s][a]
 return sum(d[-1])
def flips(m):
 for a,b,c in itertools.combinations(V,3):
  if arc(m,a,b)==arc(m,b,c)==arc(m,c,a):
   yield m^sum(1<<index[e] for e in itertools.combinations((a,b,c),2))
def relabel(m,p):
 out=0
 for a,b in E:
  x,y=(p[a],p[b]) if arc(m,a,b) else (p[b],p[a])
  if x<y:out|=1<<index[x,y]
 return out
def circ(S):return sum(1<<i for i,(a,b) in enumerate(E) if (b-a)%7 in S)
p=circ({1,2,4});b=next(flips(p));r=circ({1,2,3});reps=[p,b,r];H=[189,171,175]
checks=0
def check(x,msg):
 global checks
 checks+=1
 if not x:raise ValueError(msg)
orbits=[{relabel(m,q) for q in itertools.permutations(V)} for m in reps]
check(list(map(len,orbits))==[240,1680,720],'orbit sizes')
check(sum(map(len,orbits))==len(set.union(*orbits)),'disjoint')
# Independent score-pruned exhaustive generation, not an isomorphism-library census.
regular=set();deg=[0]*7;left=[6]*7
def generate(i,m):
 if i==21:
  check(deg==[3]*7,'terminal regularity');regular.add(m);return
 a,b=E[i];left[a]-=1;left[b]-=1
 for winner,bit in [(a,1),(b,0)]:
  deg[winner]+=1
  if all(deg[v]<=3<=deg[v]+left[v] for v in (a,b)):generate(i+1,m|(bit<<i))
  deg[winner]-=1
 left[a]+=1;left[b]+=1
generate(0,0)
check(regular==set.union(*orbits),'complete census')
classes={m:i for i,O in enumerate(orbits) for m in O}
Q=[[0,14,0],[2,9,3],[0,7,7]]
for i,O in enumerate(orbits):
 check(hp(reps[i])==H[i],'DP HP')
 brute=sum(all(arc(reps[i],a,b) for a,b in zip(q,q[1:])) for q in itertools.permutations(V))
 check(brute==H[i],'word HP')
 for m in O:
  cnt=Counter(classes[x] for x in flips(m))
  check([cnt[j] for j in range(3)]==Q[i],'equitable quotient')
# Characteristic polynomial verified by multiplication: Q^3-16Q^2+14Q+196I=0.
def mul(A,B):return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
Q2=mul(Q,Q);Q3=mul(Q2,Q)
check(all(Q3[i][j]-16*Q2[i][j]+14*Q[i][j]+196*(i==j)==0 for i in range(3) for j in range(3)),'characteristic polynomial')
pi=[Fraction(1,11),Fraction(7,11),Fraction(3,11)]
check(all(sum(pi[i]*Q[i][j] for i in range(3))==14*pi[j] for j in range(3)),'stationary')
t=[0,Fraction(10),Fraction(12)]
for i in (1,2):check(t[i]==1+sum(Fraction(Q[i][j],14)*t[j] for j in range(3)),'hitting time')
# Explicit missed-chart move: baseline one-line flip has an unrestricted move to H175.
chart=set()
lines=[tuple((i+j)%7 for j in (0,1,3)) for i in V]
lm=[sum(1<<index[e] for e in itertools.combinations(sorted(L),2)) for L in lines]
for mask in range(128):
 x=p
 for i in V:
  if mask>>i&1:x^=lm[i]
 chart.add(x)
hostile=next((m,x) for m in chart for x in flips(m) if classes[x]==2)
check(hostile[1] not in chart,'unrestricted move exits Fano chart')
print(json.dumps({'status':'PASS','checks':checks,'universe':'all2640 labelled regular tournaments on7 vertices; every cyclic-triangle reversal','representatives':reps,'Hamiltonian_counts':H,'orbit_sizes':list(map(len,orbits)),'quotient_neighbor_counts':Q,'stationary':[str(x) for x in pi],'stationary_H':str(sum(x*y for x,y in zip(pi,H))),'mean_Paley_hitting':[str(x) for x in t],'mean_Paley_return':11,'charpoly':'(x-14)(x^2-2x-14)','chart_exit_masks':hostile},sort_keys=True,indent=2))
