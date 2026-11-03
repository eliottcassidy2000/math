from itertools import combinations
n=7
forward={(0,1),(0,3),(0,4),(0,5),(1,2),(2,3),(2,6),(3,4),(3,6),(4,5),(4,6),(5,6)}
z=[[0]*n for _ in range(n)]
for i,j in combinations(range(n),2):
 z[i][j]=int((i,j) in forward);z[j][i]=1-z[i][j]
def check(c,msg):
 if not c: raise ValueError(msg)
def scorekey(a):
 s=tuple(sorted(map(sum,a))); t=tuple(sorted(n-1-x for x in s)); return min(s,t)
def switch(a,p):
 edges={frozenset((p[i],p[i+1])) for i in range(n-1)}
 check(all(a[p[i]][p[i+1]] for i in range(n-1)),'not a path')
 b=[r[:] for r in a]
 for i,j in combinations(range(n),2):
  if frozenset((i,j)) not in edges: b[i][j],b[j][i]=b[j][i],b[i][j]
 return b
def blue(a,p):
 return all(a[p[i]][p[j]]==a[p[n-1-j]][p[n-1-i]] for i,j in combinations(range(n),2))
w=z
cycle={frozenset((g,(g+1)%n)) for g in range(n)}
z=[r[:] for r in w]
for i,j in combinations(range(n),2):
 if frozenset((i,j)) not in cycle: z[i][j],z[j][i]=z[j][i],z[i][j]
leaves=[]; ps=[]
for g in range(n):
 h=(g+1)%n
 check(w[g][h]==1,'cycle missing')
 a=[r[:] for r in w];a[g][h]=0;a[h][g]=1
 p=[(g+1+i)%n for i in range(n)]
 check(switch(a,p)==z,'bad hub spoke')
 check(switch(z,p)==a,'bad involution')
 check(not blue(a,p),'unexpected blue')
 leaves.append(a);ps.append(p)
keys=[scorekey(z)]+[scorekey(a) for a in leaves]
check(len(set(keys))==n+1,'merger scores not distinct')
rules=[({0,1,2},3),({2,3,4},5),({5},6),({6},0)]
for g in range(n):
 f=next(dst for domain,dst in rules if g in domain)
 u=switch(switch(leaves[g],ps[g]),ps[f]);check(u==leaves[f],'invalid macro')
 diffs=[(i,j) for i,j in combinations(range(n),2) if leaves[g][i][j]!=u[i][j]]
 check(set(diffs)=={tuple(sorted((g,(g+1)%n))),tuple(sorted((f,(f+1)%n)))},'bad macro support')
 print(f'{g} -> {f}:  {scorekey(leaves[g])},  black-black')
print('center',scorekey(z))
print('all checks passed')
def H(a):
 dp=[[0]*n for _ in range(1<<n)]
 for i in range(n): dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   if dp[mask][i]:
    for j in range(n):
     if not (mask>>j)&1 and a[i][j]: dp[mask|(1<<j)][j]+=dp[mask][i]
 return sum(dp[-1])
print('H center, leaves:',[H(a) for a in [z]+leaves])
from itertools import permutations
def H_brute(a):
 return sum(all(a[p[i]][p[i+1]] for i in range(n-1)) for p in permutations(range(n)))
counts=[H(a) for a in [z]+leaves]
check(counts==[H_brute(a) for a in [z]+leaves],'DP and permutation counts disagree')
check(len(set(counts))==n+1,'path counts not distinct')
print('DP counts verified by permutation enumeration')
