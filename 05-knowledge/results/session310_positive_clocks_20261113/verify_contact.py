from itertools import permutations,product
from collections import Counter

def require(x):
    if not x: raise ValueError('exact contact control failed')
P=list(permutations(range(6)))
distribution=Counter()
for k in range(512):
    def arc(i,j):
        if i//3==j//3: return (j-i)%3==1
        if i<j: return k>>(3*i+j-3)&1
        return 1-arc(j,i)
    h=sum(all(arc(i,j) for i,j in zip(p,p[1:])) for p in P)
    require((h==9)==(k in (0,511)))
    distribution[h]+=1
require(distribution=={9:2,17:18,23:36,27:12,29:90,31:36,33:72,37:144,41:36,43:36,45:30})
f=[0,3,6,6]
h=0
for runs in product(range(1,4),repeat=3):
    for a in range(3):
        if tuple(sum((a+t)%3==i for t in range(sum(runs))) for i in range(3))==runs:
            h+=f[runs[0]]*f[runs[1]]*f[runs[2]]
require(h==3159)
print('PASS: 512 mixed-pair orientations and the cyclic triple visit count')
