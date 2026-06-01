"""Lonely Runner Conjecture: p-adic tree of denominators investigation.
Pure Python exact integer arithmetic."""
from math import gcd
from fractions import Fraction

def phi(q):
    # Euler totient
    result = q
    p = 2
    n = q
    while p*p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result

def is_safe(r, q, n):
    """residue r in [0,q), safe iff min(r,q-r) >= q/n, i.e. n*min(r,q-r) >= q."""
    mr = r if r <= q-r else q-r
    return n*mr >= q

def lonely_count_at_q(speeds, q):
    """Count a in (Z/q)* such that all runners safe at t=a/q. n = len(speeds)+1.
    Returns (count, phi(q))."""
    n = len(speeds)+1
    cnt = 0
    tot = 0
    for a in range(1, q):
        if gcd(a, q) != 1:
            continue
        tot += 1
        ok = True
        for v in speeds:
            r = (v*a) % q
            if not is_safe(r, q, n):
                ok = False
                break
        if ok:
            cnt += 1
    return cnt, tot

def lonely_count_all_a(speeds, q):
    """Count ALL a in [1,q-1] (not just units) such that lonely. Used for smallest-q search."""
    n = len(speeds)+1
    cnt = 0
    for a in range(1, q):
        ok = True
        for v in speeds:
            r = (v*a) % q
            if not is_safe(r, q, n):
                ok = False
                break
        if ok:
            cnt += 1
    return cnt

def f_pk(speeds, p, k):
    """lonely fraction f(p^k) = |O ∩ B|/phi(p^k) over units a."""
    q = p**k
    cnt, tot = lonely_count_at_q(speeds, q)
    return Fraction(cnt, tot) if tot else Fraction(0), cnt, tot

def V_target(m):
    """V = (1 - 2/n)^m, n=m+1."""
    n = m+1
    return Fraction(n-2, n)**m

if __name__ == "__main__":
    print("phi(8)=", phi(8), "phi(81)=", phi(81))
