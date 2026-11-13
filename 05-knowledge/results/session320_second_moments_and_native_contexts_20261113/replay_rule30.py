from collections import Counter, deque
from itertools import combinations
from fractions import Fraction

I=(0,1,2,3)
G=((0,2,3,3),(2,0,1,1))
def compose(a,b):return tuple(a[j] for j in b)
def trace(a):return sum(i==j for i,j in enumerate(a))
def M(w):
    t=I
    for c in w:t=compose(t,G[int(c)])
    return t
def req(x):
    if not x:raise RuntimeError('replay mismatch')
S={I}; q=deque([I])
while q:
    a=q.popleft()
    for b in G:
        c=compose(a,b)
        if c not in S:S.add(c);q.append(c)
S=sorted(S); L=S.copy();L.remove(I)
req(len(S)==41 and Counter(len(set(a)) for a in S)==Counter({1:4,2:30,3:6,4:1}))
D={(a,b):trace(compose(a,b)) for a in S for b in S}
def partition(base,probes):
    by={}
    for a in base:by.setdefault(tuple(D[a,b] for b in probes),[]).append(a)
    return set(frozenset(v) for v in by.values())

allp=partition(S,S)
req(len(allp)==23)
req(partition(S,[M(w) for w in ('','0','1','01')])==allp)
req(all(len(partition(S,P))<23 for P in combinations(S,3)))
req(all(len(partition(L,P))<22 for P in combinations(S,3)))
req(all(D[M('01'),M(w)]==D[M('10'),M(w)] for w in ('','0','1')))
req(D[M('01'),M('01')]!=D[M('10'),M('01')])
future=partition(L,S);req(len(future)==22)
levels=[{I}]
for n in range(11):
    levels.append({compose(a,b) for a in levels[-1] for b in G})
req([len(x) for x in levels]==[1,2,4,8,16,24,30,34,36,36,36,36])
req(levels[8]==levels[11])
for m in range(3,11):req(partition(L,levels[m])==future)
for m in range(3):req(all(len({D[a,b] for a in S})>1 for b in levels[m]))
from itertools import product
req([w for w in map(''.join,product('01',repeat=3)) if len({D[a,M(w)] for a in S})==1]==['010'])
print('23 full classes, 22 nonempty; four minimal probes; each depth m>=3 separates')

# Transfer T on test vectors and Cayley-Hamilton certificate.
def T(f):return {a:sum(f[compose(a,b)] for b in G) for a in S}
diff={a:D[M('01'),a]-D[M('10'),a] for a in S}
req(set(diff.values())=={-1,0,1})
for sign in (-1,1):
    f={a:int(diff[a]==sign) for a in S};u=[f]
    for n in range(4):u.append(T(u[-1]))
    req([v[I] for v in u]==[0,0,1,2,3])
    residual={a:u[4][a]-2*u[3][a]+u[2][a]-2*u[1][a]+2*u[0][a] for a in S}
    for _ in S:
        req(residual[I]==0)
        residual=T(residual)
print('canonical pair sign GF z^2/[(1-z)(1-z-2z^3)] verified')

# General signed-balance transfer at memory length two.
for a in S:
    for m in (2,3,4,5):
        def total(f,m):
            for _ in range(m):f=T(f)
            return f[I]
        req(total({b:D[a,b] for b in S},m)==1<<m)
print('same-depth signed balance controls passed')
