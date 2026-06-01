"""Lonely Runner Conjecture: Galois-orbit / equidistribution investigation.
Pure Python, exact integer arithmetic, no internet.

Setup: distinct positive integer speeds v_1..v_m, n=m+1, threshold 1/n.
At rational time t=a/q (q prime, a in (Z/q)*), runner i sits at residue
v_i*a mod q. SAFE iff min(r,q-r) >= q/n, i.e. r in safe box
B_q = [ceil(q/n), q-ceil(q/n)]. LONELY at a/q iff every runner is safe.

Galois orbit O_q = {(v_i a mod q)_i : a in (Z/q)*}, size q-1.
LRC at denominator q  <=>  O_q meets B_q^m.
"""
from math import gcd
from fractions import Fraction


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    i = 3
    while i * i <= n:
        if n % i == 0:
            return False
        i += 2
    return True


def primes_up_to(N):
    return [p for p in range(2, N + 1) if is_prime(p)]


def ceil_div(a, b):
    return -((-a) // b)


def safe_box(q, n):
    """Return (lo, hi): residue r is safe iff lo <= r <= hi."""
    lo = ceil_div(q, n)          # ceil(q/n)
    hi = q - lo
    return lo, hi


def lonely_count(speeds, q, n):
    """Number of a in (Z/q)* for which all runners are safe (q prime)."""
    lo, hi = safe_box(q, n)
    cnt = 0
    for a in range(1, q):  # all nonzero residues = (Z/q)* since q prime
        ok = True
        for v in speeds:
            r = (v * a) % q
            if r < lo or r > hi:
                ok = False
                break
        if ok:
            cnt += 1
    return cnt


def is_galois_lonely(speeds, q, n):
    return lonely_count(speeds, q, n) > 0


def fraction_lonely(speeds, q, n):
    return Fraction(lonely_count(speeds, q, n), q - 1)


def box_volume(q, n, m):
    """Exact volume of B_q^m as fraction: ((hi-lo+1)/q)^m, the per-runner
    measure of being in the (closed) safe box among all residues mod q."""
    lo, hi = safe_box(q, n)
    width = hi - lo + 1
    return Fraction(width, q) ** m


def limiting_volume(n, m):
    """(1 - 2/n)^m as a Fraction."""
    return Fraction(n - 2, n) ** m


def smallest_galois_lonely_prime(speeds, n, pmax=2000):
    for q in primes_up_to(pmax):
        if q <= n:  # box empty or trivial when q too small relative to n
            # safe box: lo=ceil(q/n), need lo <= hi i.e. box nonempty
            lo, hi = safe_box(q, n)
            if hi < lo:
                continue
        if is_galois_lonely(speeds, q, n):
            return q
    return None
