"""Standalone finite checks for the explicit pair and physical reader in 23_chan.md."""
from itertools import combinations, permutations
from math import comb, gcd

rows = [
    (4584474270165959682, 77, 2024, 21789),
    (9169081464381770758, 15, 960, 14175),
    (9160127196319778056, 35, 1352, 16947),
    (9169072685987728394, 33, 1440, 18225),
]

def subsets(s):
    a = s
    while True:
        yield a
        if not a:
            break
        a = (a - 1) & s

def profile(code):
    n = 8
    A = [[(code >> (n*i+j)) & 1 for j in range(n)] for i in range(n)]
    for i in range(n):
        if not (A[i][i] == 0): raise AssertionError("slot proof check failed")
        for j in range(i):
            if not (A[i][j] + A[j][i] == 1): raise AssertionError("slot proof check failed")
    dp = [[0]*n for _ in range(1 << n)]
    h = [1] + [0]*((1 << n)-1)
    for s in range(1, 1 << n):
        for j in range(n):
            if (s >> j) & 1:
                r = s ^ (1 << j)
                dp[s][j] = sum(dp[r][i]*A[i][j] for i in range(n)) if r else 1
        h[s] = sum(dp[s])
    q2 = [sum(h[a]*h[s^a] for a in subsets(s)) for s in range(1 << n)]
    q3 = sum(q2[a]*h[255^a] for a in subsets(255))
    return h[-1], q2[-1], q3

def reader(code):
    n, v = 6, 5
    A = [[0]*n for _ in range(n)]
    for bit, (i,j) in enumerate(combinations(range(n),2)):
        A[i][j] = (code >> bit) & 1
        A[j][i] = 1-A[i][j]
    for i in range(n):
        reach = {i}
        for _ in range(n):
            reach |= {y for x in reach for y in range(n) if A[x][y]}
        if not (len(reach) == n): raise AssertionError("slot proof check failed")
    F = [0]*(n+1)
    for w in permutations(range(v)):
        m, o = 0, A[v][w[0]]+A[w[-1]][v]
        for x,y in zip(w,w[1:]):
            can_visit = A[x][v] and A[v][y]
            if not A[x][y]:
                if not can_visit:
                    break
                m += 1
            else:
                o += can_visit
        else:
            for k in range(o+1):
                F[m+k] += comb(o,k)
    return F

p = [profile(code) for code,*_ in rows]
if not (p == [r[1:] for r in rows]): raise AssertionError("slot proof check failed")
a = tuple(x*y for x,y in zip(p[0],p[1]))
b = tuple(x*y for x,y in zip(p[2],p[3]))
if not (a == (1155,1943040,308859075)): raise AssertionError("slot proof check failed")
if not (b == (1155,1946880,308859075)): raise AssertionError("slot proof check failed")
if not ((rows[0][0]^rows[2][0]).bit_count() == 14): raise AssertionError("slot proof check failed")
if not ((rows[1][0]^rows[3][0]).bit_count() == 8): raise AssertionError("slot proof check failed")
readers = (29130,29002,29000)
if not (tuple(x^y for x,y in zip(readers,readers[1:])) == (128,2)): raise AssertionError("slot proof check failed")
for code,c in zip(readers,(-3,-1,1)):
    r = reader(code)
    if not (r[1:] == [12+18+c,12+9,4,0,0,0]): raise AssertionError("slot proof check failed")
    if not (r[0] == 12+c): raise AssertionError("slot proof check failed")
# The nonconstant reader coefficients in ordered nonempty cover basis are
# (12+18+c,12+9,4) by Q3=3H+3P2+P3 and Q2=2H+P2.
for s in range(20):
    box = [[4*t[2]*3**s+9*t[1]*2**s+c*t[0] for c in (-3,-1,1)] for t in (a,b)]
    if not (all(u[1]-u[0] == u[2]-u[1] == 2310 for u in box)): raise AssertionError("slot proof check failed")
    if not (all(y-x == 34560*2**s for x,y in zip(*box))): raise AssertionError("slot proof check failed")
    if not (gcd(2310,34560*2**s) == 30): raise AssertionError("slot proof check failed")
for m in range(2,11):
    g = [a[1]**(m-k)*b[1]**k for k in range(m+1)]
    for k in range(m-1):
        if not (g[k+2]-2*g[k+1]+g[k] == 3840**2*a[1]**(m-k-2)*b[1]**k): raise AssertionError("slot proof check failed")
print("PASS: four profiles, eleven internal edits, three strong literal readers, both strides")
