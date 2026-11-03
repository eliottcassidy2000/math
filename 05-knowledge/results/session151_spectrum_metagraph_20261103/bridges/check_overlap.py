def require(x):
 if not x: raise AssertionError("bridge check failed")
from itertools import combinations,permutations
from collections import Counter
for n in range(3,7):
 p=tuple(range(n)); fixed={(i,i+1) for i in range(n-1)}; pairs=list(combinations(range(n),2)); chords=[e for e in pairs if e not in fixed]
 stats=Counter()
 for bits in range(1<<len(chords)):
  edges=set(fixed)|{e for i,e in enumerate(chords) if bits>>i&1}
  def arc(E,a,b): return ((min(a,b),max(a,b)) in E)==(a<b)
  K=set(fixed)|set(chords).difference(edges)
  hs=[Counter(),Counter()]; ps=[set(),set()]
  for q in permutations(p):
   r=sum(tuple(sorted((a,b))) in fixed for a,b in zip(q,q[1:]))
   for j,E in enumerate((edges,K)):
    if all(arc(E,a,b) for a,b in zip(q,q[1:])): hs[j][r]+=1; ps[j].add(q)
  require(ps[0]&ps[1]=={p})
  require({q[::-1] for q in ps[0] if not any(tuple(sorted((a,b))) in fixed for a,b in zip(q,q[1:]))}=={q for q in ps[1] if not any(tuple(sorted((a,b))) in fixed for a,b in zip(q,q[1:]))})
  closed=arc(edges,n-1,0); require(hs[0][n-2]==(n-1 if closed else 0))
  require(hs[1][n-2]==(0 if closed else n-1))
  d=sum(hs[0].values())-sum(hs[1].values())
  require(d==((n-1) if closed else -(n-1))+sum(hs[0][r]-hs[1][r] for r in range(1,n-2)))
  stats[d]+=1
 print(f'n={n}: {1<<len(chords)} marked tilings checked; delta range [{min(stats)},{max(stats)}]; distinct delta={len(stats)}')
# Explicit closed-side monotonicity counterexample, independent of the profile summary.
n=5
P=tuple(range(n)); F={(i,i+1) for i in range(n-1)}
E=set(combinations(range(n),2))-{(0,4)}
K=F|(set(combinations(range(n),2))-F-E)
def profile(E):
 def arc(a,b): return ((min(a,b),max(a,b)) in E)==(a<b)
 result=[0]*n
 for q in permutations(P):
  if all(arc(a,b) for a,b in zip(q,q[1:])):
   result[sum(tuple(sorted((a,b))) in F for a,b in zip(q,q[1:]))]+=1
 return result
assert profile(E)==[0,4,0,4,1] and profile(K)==[0,4,6,0,1]
print('n=5 explicit closed/open counterexample: profiles [0,4,0,4,1] and [0,4,6,0,1], H=9 and H=11')
