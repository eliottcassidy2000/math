from math import comb, isqrt

def need(ok):
    if not ok:
        raise ValueError("certificate failed")

def trim(a):
    while a and not a[-1]:
        a.pop()
    return a

def mul(a, b, p):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] = (c[i+j] + x*y) % p
    return trim(c)

def rem(a, h, p):
    a = trim(a[:])
    while len(a) >= len(h):
        i, x = len(a) - len(h), a[-1]
        for j, y in enumerate(h):
            a[i+j] = (a[i+j] - x*y) % p
        trim(a)
    return a

def power(a, n, h, p):
    b = [1]
    while n:
        if n & 1:
            b = rem(mul(a, b, p), h, p)
        a = rem(mul(a, a, p), h, p)
        n //= 2
    return b

def evolve(p, lam, n):
    row = {0: 1}
    for _ in range(n):
        nxt = {}
        for x, v in row.items():
            for dx, a in ((0, 1+lam), (1, 1), (-1, lam)):
                nxt[x+dx] = (nxt.get(x+dx, 0) + a*v) % p
        row = nxt
    return row.get(0, 0)

def direct(p, lam, n):
    return sum((comb(n,j)%p)**2 * pow(lam,j,p) for j in range(n+1)) % p

def digit(p, lam, n):
    a = 1
    while n:
        n, r = divmod(n, p)
        a = a * direct(p, lam, r) % p
    return a

got = []
for p in range(5, 128):
    if any(p % d == 0 for d in range(2, isqrt(p)+1)):
        continue
    m = (p-1)//2
    h = [comb(m,j)**2 % p for j in range(m+1)]
    A = [256*x%p for x in [1,-3,6,-7,6,-3,1]]
    B = [0, 0, 1, p-2, 1]
    if power(A,p,h,p) == rem(mul(A,power(B,p-1,h,p),p),h,p):
        got.append(p)
    for lam in (0,1,2,p-1):
        need(evolve(p,lam,m) == direct(p,lam,m))
        need(evolve(p,lam,2*p+3) == digit(p,lam,2*p+3))
need(got == [5,7,11,13,17,19,23,29,31,41,47,59,71])
print("p=5..127:", *got)
