"""Locking theorem (THM-4849 (D)), both directions, and the interval-locking primes of the rotational tournaments.

Let w = (v_1..v_k) be a cycle word, H = sum v_i, v >= 1 an appended letter, D' = 2^(H+v) - 3^(k+1), and p^e | D'
with p not dividing 6 (2^v - 3).  Claim: every point x_i = C(rot^i w)/D(w) of the cycle is = xi_v = 1/(2^v - 3)
(mod p^e) iff v_i = v (mod ord_(p^e)(2)) for every letter.  (Then p does not divide D(w), automatically.)
So the attachment response at p^e is marking-independent and maximal exactly for these "locked" words.

Second part: the interval word (1^(k-1), k+2) of R_q (connection set {1..k}, q = 2k+1) is locked under the appended
letter v = 1 exactly at the primes of gcd(2^m - 1, 3^m - 1), m = (q+1)/2; checked for every odd prime q < 100.
"""
import itertools
from math import gcd
from sympy import factorint, n_order, primerange


def C_of(w):
    s, out, k = 0, 0, len(w)
    for j, v in enumerate(w):
        out += 3 ** (k - 1 - j) * 2 ** s
        s += v
    return out


def necklaces(k, maxletter):
    seen = set()
    for w in itertools.product(range(1, maxletter + 1), repeat=k):
        rots = [w[i:] + w[:i] for i in range(k)]
        m = min(rots)
        if m not in seen:
            seen.add(m)
            yield m


tests = agree = locked = 0
for k in range(1, 6):
    for w in necklaces(k, 7 if k <= 4 else 5):
        H = sum(w)
        D = 2 ** H - 3 ** k
        if D == 0:
            continue
        rots = [w[i:] + w[:i] for i in range(k)]
        for v in range(1, 9):
            Dn = 2 ** (H + v) - 3 ** (k + 1)
            for p, e in factorint(abs(Dn)).items():
                if p in (2, 3) or (2 ** v - 3) % p == 0:
                    continue
                for j in range(1, e + 1):
                    pj = p ** j
                    o = n_order(2, pj)
                    letters = all((vi - v) % o == 0 for vi in w)
                    if letters:
                        assert D % p != 0           # automatic
                    if D % p == 0:
                        continue                    # points not p-integral; outside the statement
                    xi = pow((2 ** v - 3) % pj, -1, pj)
                    pts = all(C_of(r) * pow(D, -1, pj) % pj == xi for r in rots)
                    tests += 1
                    agree += (pts == letters)
                    locked += pts
print('locking theorem: %d (word, letter, prime power) tests, %d agree, %d locked' % (tests, agree, locked))

print('\ninterval words of R_q under appended v = 1: locked primes vs gcd(2^m - 1, 3^m - 1), m = (q+1)/2')
for q in primerange(3, 100):
    k = (q - 1) // 2
    m = (q + 1) // 2
    w = (1,) * (k - 1) + (k + 2,)
    H = sum(w)
    D = 2 ** H - 3 ** k
    Dn = 2 ** (H + 1) - 3 ** (k + 1)
    rots = [w[i:] + w[:i] for i in range(k)]
    lock = []
    for p in factorint(abs(Dn)):
        if p in (2, 3) or D % p == 0:
            continue
        if all(C_of(r) * pow(D, -1, p) % p == (-1) % p for r in rots):    # xi_1 = 1/(2-3) = -1
            lock.append(p)
    g = gcd(2 ** m - 1, 3 ** m - 1)
    pred = sorted(factorint(g)) if g > 1 else []
    print('  q=%2d  m=%2d  locked %s  gcd primes %s  %s' % (q, m, sorted(lock), pred, 'ok' if sorted(lock) == pred else 'MISMATCH'))
