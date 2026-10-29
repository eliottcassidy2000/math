"""Attachment lemma check.  For a word w (k letters, sum H, D = 2^H - 3^k), marked point x = C(w)/D, and a letter v,
the extended word wv has C(wv) = 3 C(w) + 2^H and D' = 2^(H+v) - 3^(k+1).  Claim: for every prime power p^j | D'
with p not dividing 6 (2^v - 3) D, p^j | C(wv) iff x = 1/(2^v - 3) (mod p^j), the fixed point of the appended
branch x -> (3x+1)/2^v.  Locking corollary: if 2^(v_i) = 2^v (mod p) for every letter and ord_p(3 2^-v) does not
divide k, then every point of the cycle is = 1/(2^v - 3) (mod p)  [p not dividing 6 (2^v-3) D].
Also: other multipliers a (C(wv) = a C(w) + 2^H, fixed point 1/(2^v - a))."""
import random
from math import gcd
from sympy import factorint, n_order
def C_of(w, a):
    s, out, kk = 0, 0, len(w)
    for j, v in enumerate(w):
        out += a ** (kk - 1 - j) * 2 ** s
        s += v
    return out
random.seed(7)
tests = fails = 0
lock_tests = lock_fails = 0
for trial in range(4000):
    a = random.choice((3, 5, 7))
    k = random.randint(1, 6)
    w = tuple(random.randint(1, 7) for _ in range(k))
    v = random.randint(1, 7)
    H = sum(w)
    D = 2 ** H - a ** k
    Dn = 2 ** (H + v) - a ** (k + 1)
    if D == 0 or Dn == 0 or 2 ** v == a:
        continue
    C = C_of(w, a); Cn = C_of(w + (v,), a)
    assert Cn == a * C + 2 ** H
    for p, e in factorint(abs(Dn)).items():
        if p in (2, 3) or (a % p == 0) or D % p == 0 or (2 ** v - a) % p == 0:
            continue
        for j in range(1, e + 1):
            pj = p ** j
            lhs = Cn % pj == 0
            x = C * pow(D, -1, pj) % pj
            fixed = pow((2 ** v - a) % pj, -1, pj)
            rhs = x == fixed
            tests += 1
            if lhs != rhs:
                fails += 1
    # locking corollary at small primes p
    for p in (5, 7, 11, 13, 17, 19, 23):
        if a % p == 0 or D % p == 0:
            continue
        r = pow(2, w[0], p)
        if not all(pow(2, vi, p) == r for vi in w):
            continue
        alpha = a * pow(r, -1, p) % p
        if alpha == 1 or (2 ** w[0] - a) % p == 0:
            continue
        if k % n_order(alpha, p) == 0:
            continue
        fixed = pow((r - a) % p, -1, p)
        pts = [C_of(w[i:] + w[:i], a) * pow(D, -1, p) % p for i in range(k)]
        lock_tests += 1
        if any(x != fixed for x in pts):
            lock_fails += 1
print('attachment lemma: %d prime-power tests, %d failures' % (tests, fails))
print('locking corollary: %d tests, %d failures' % (lock_tests, lock_fails))
# the 7x+169 cycle through 79 and the integer Collatz orbit of 79
def U(n, a, b):
    t = a * n + b
    while t % 2 == 0:
        t //= 2
    return t
n, cyc = 79, []
while n not in cyc:
    cyc.append(n); n = U(n, 7, 169)
print('7x+169 orbit of 79:', cyc, 'returns to', n)
n, orb = 79, [79]
while n != 1:
    n = U(n, 3, 1); orb.append(n)
print('3x+1 (odd-step) orbit of 79:', orb)
