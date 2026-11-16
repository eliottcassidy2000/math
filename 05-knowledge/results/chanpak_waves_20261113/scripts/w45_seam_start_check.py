import random, itertools, sys
sys.path.insert(0, '.')
from w45_hk import ham_split, sub
rng = random.Random(45)
bad = 0
for trial in range(200):
    n = rng.randint(2, 7)
    out = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < .5: out[i] |= 1 << j
            else: out[j] |= 1 << i
    A = 0
    H, st, en = ham_split(n, out)
    N = n + 3; x, y, z = n, n + 1, n + 2
    o = out[:] + [0, 0, 0]
    o[x] |= 1 << y; o[y] |= 1 << z; o[z] |= 1 << x
    for t in (x, y, z):
        for i in range(n): o[t] |= 1 << i
    o[x] &= ~(1 << A); o[A] |= 1 << x
    Hs = ham_split(N, o)[0]
    ob = out[:] + [0]; b = n
    ob[A] |= 1 << b
    for i in range(1, n): ob[b] |= 1 << i
    Hp = ham_split(n + 1, ob)[0]
    if Hs != 3 * Hp + 2 * st[A]: bad += 1
print("seam identity H = 3H(T+) + 2*H_start(A):", 200 - bad, "of 200 hold")
# one-apex tape closed form
def onetape(u):
    m = len(u); A = m; out = [0] * (m + 1)
    for i in range(m):
        for j in range(i + 1, m): out[i] |= 1 << j
        if u[i]: out[i] |= 1 << A
        else: out[A] |= 1 << i
    return ham_split(m + 1, out)[0]
def cform(u):
    m = len(u); c = 0
    for b in range(m):
        for a in range(b + 1, m):
            if u[b] == 0 and u[a] == 1:
                d = a - b; c += 1 if d == 1 else 2 ** (d - 2)
    return c
badc = 0
for m in range(1, 11):
    for u in itertools.product((0, 1), repeat=m):
        if onetape(u) != 1 + 2 * cform(u): badc += 1
print("one-apex closed form c(u) = sum over (0 before 1) of w(dist), w=1,1,2,4,..: failures", badc)
