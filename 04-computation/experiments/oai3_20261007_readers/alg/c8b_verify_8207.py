from fractions import Fraction as Fr
def f(word, x=Fr(-1)):
    for c in word:
        x = (3*x + 1) / 2**c
    return x
for a in range(1, 6):
    u, v = (2, 2, 10, a), (6, 3, 2, 1, a + 2)
    assert f(u) == f(v), a
    print("a=%d: f_(2,2,10,%d)(-1) = f_(6,3,2,1,%d)(-1) = %s" % (a, a, a + 2, f(u)))
# genuine Collatz switch (THM-4555 (iv), D = 1): odd n with word 1^r (2,2,10,a,...) merges with m = (n+1)/2 - 1
def U(n):
    n = 3*n + 1
    e = 0
    while n % 2 == 0:
        n //= 2; e += 1
    return n, e
def word(n, k):
    w = []
    for _ in range(k):
        n, e = U(n); w.append(e)
    return w
found = 0
for n in range(3, 2*10**7, 2):
    w = word(n, 6)
    # r leading ones then (2,2,10)
    r = 0
    while r < len(w) and w[r] == 1: r += 1
    if r >= 1 and w[r:r+3] == [2, 2, 10]:
        m = (n + 1)//2 - 1
        j = r + 4
        x, y = n, m
        for _ in range(j):
            x, _ = U(x); y, _ = U(y)
        print("n = %d (word %s...), m = (n+1)/2-1 = %d: U^%d(n) = %d, U^%d(m) = %d, equal: %s"
              % (n, w, m, j, x, j, y, x == y))
        found += 1
        if found >= 3: break
