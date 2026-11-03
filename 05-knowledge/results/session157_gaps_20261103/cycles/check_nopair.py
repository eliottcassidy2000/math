from itertools import combinations,permutations
from collections import defaultdict
E=list(combinations(range(6),2)); fixed={(i,i+1):True for i in range(5)};fixed[0,5]=False;free=[e for e in E if e not in fixed]
def orbit_key(o):
 def a(u,v): return o[min(u,v),max(u,v)]==(u<v)
 ss=[]
 for shift in range(6):
  for dual in (0,1):
   f=lambda x:(shift+(x if dual==0 else -x))%6
   ss.append(''.join('1' if a(f(i),f(j)) != bool(dual) else '0' for i,j in free))
 return min(ss)
def counts(o):
 def a(u,v):return o[min(u,v),max(u,v)]==(u<v)
 cycles=[]
 for k in (3,5):
  for vs in combinations(range(6),k):
   for p in permutations(vs[1:]):
    z=(vs[0],)+p
    if all(a(z[i],z[(i+1)%k]) for i in range(k)):cycles.append(set(z))
 t=[x for x in cycles if len(x)==3]
 return len(t),len(cycles),sum(not (x&y) for x,y in combinations(t,2))
groups=defaultdict(list)
for m in range(512):
 o=dict(fixed);o.update((e,bool(m>>i&1)) for i,e in enumerate(free))
 groups[orbit_key(o)].append((m,counts(o)))
expected = {
'000000001':(4,7,0,12), '000000101':(5,9,0,12),
'000000111':(4,8,0,6), '000001101':(6,12,0,12),
'000100001':(6,12,0,12), '000100111':(6,11,0,12),
'000110011':(7,16,0,12),
}
actual={}
for key, ms in sorted(groups.items()):
    if len({x[1] for x in ms}) != 1: raise RuntimeError('symmetry mismatch')
    t,C,D=ms[0][1]
    if D==0: actual[key]=(t,C,D,len(ms))
if len(groups)!=52 or actual!=expected: raise RuntimeError('classification mismatch')
if any(C==10 for t,C,D,size in actual.values()): raise RuntimeError('counterexample')
print('PASS: 512 orientations, 52 orbits; 7 D=0 orbits match the table; C != 10.')
