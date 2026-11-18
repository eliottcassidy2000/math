from itertools import product
from math import isqrt

def h(adj):
    n = len(adj)
    dp = [[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i] = 1
    for s in range(1,1<<n):
        for j, k in enumerate(dp[s]):
            if not k: continue
            t = adj[j] & ~s
            while t:
                b = t & -t
                dp[s|b][b.bit_length()-1] += k
                t ^= b
    return sum(dp[-1])

def apex(w):
    n = len(w)
    adj = [(1<<n)-(1<<(i+1)) for i in range(n)]+[0]
    for i,b in enumerate(w):
        if b: adj[-1] |= 1<<i
        else: adj[i] |= 1<<n
    return adj

def prediction(w, i):
    if w[i] != 0 or any(w[:i]): return False
    tail = w[i+1:]
    b = next((j for j,x in enumerate(tail) if x),len(tail))
    s = tail[b:]
    if b < 1 or s.count(0)>1: return False
    return 1+(1<<b)

total=0
for n in range(1,9):
    values={w:h(apex(w)) for w in product(range(2),repeat=n)}
    for w,a in values.items():
        for i,x in enumerate(w):
            if x: continue
            k=values[w[:i]+(1,)+w[i+1:]]
            ratio=k//a if k>a and k%a==0 else False
            assert ratio == prediction(w,i),(w,i,a,k,ratio)
            total+=1
print('direct-DP rising spoke tests:',total)
def from_code(n,code):
    adj=[0]*n
    e=0
    for i in range(n):
        for j in range(i+1,n):
            a,b=(j,i) if code>>e&1 else (i,j)
            adj[a]|=1<<b
            e+=1
    return adj

for n,code,p in [(4,4,5),(7,179,35),(7,185,49)]:
    assert h(from_code(n,code))==p
    print('literal count:',n,code,p)
for j in range(1,30):
    a=1<<j
    for e in [-1,1]:
        q=a*a+(a+e)**2
        assert (isqrt(q)**2==q)==(j==2 and e==-1)
print('Pell check: only adjacent 3,4 among powers of two tested')

def lift(code,a,b):
    original=from_code(7,code)
    # Ground set 0, A, B, 4, 5, 6. Vertex 3 is compressed into B.
    nodes=[0]+[1]*a+[2]*b+[4,5,6]
    return [sum(1<<j for j,v in enumerate(nodes)
                if (w==v and i<j) or (original[w]>>v&1))
            for i,w in enumerate(nodes)]

checks=0
for a,b in product(range(1,5),repeat=2):
    u,v=1<<a,1<<b
    base=(u-1)*(v-1)+2
    values=[base, (2*u-2)*v-u+3, 3*base, (6*u-6)*v-3*u+11]
    for c,val in zip((6,14,524294,524302),values):
        assert h(lift(c,a,b))==val
        checks+=1
    assert values[3]-7*base==(4-v)*(u-1)-6
    assert values[3]-values[1]-values[2]+base==2*(u-1)*v+2
print('two-switch lift:',checks,'literal DP checks')
