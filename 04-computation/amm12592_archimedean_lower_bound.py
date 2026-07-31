#!/usr/bin/env python3
"""Certified lower bounds on rho(m) for AMM 12592, by expansion at u = -1.

In deviation form the within-block system (THM-3008) is homogeneous:

    sum_k E_k(u) (1+u)^{L_k}  =  +- u^{m-1},        L_k = m-1-k-a_k,
    |[u^i] E_k| <= C(a_k, i),   [u^i]E_k = C(a_k,i)  (mod 2).

TWO REDUCTIONS.

(A) MOD 2 -- no obstruction, ever.  The forced half-integer slots are the i
    with C(a_k,i) odd, i.e. i a binary submask of a_k, and
    sum_{i subset a} u^i = (1+u)^a over F_2.  So the forced pattern reduces to
        sum_k (1+u)^{a_k+L_k} = sum_k (1+u)^{m-1-k} = [(1+u)^m - 1]/u,
    which is u^{m-1} over F_2 because m is a power of two.  Every coefficient
    below the top vanishes identically, for EVERY depth profile: inside a
    shell there is no 2-adic obstruction at all, and the one surviving parity
    is exactly the free label of z = 1^m.  Parity is a red herring here; the
    obstruction is archimedean.

(B) AT u = -1 -- a genuine size obstruction.  Put u = -1+v, so
    (1+u)^{L_k} = v^{L_k} and the right side becomes +-(-1+v)^{m-1}.
    Comparing coefficients of v^d and using
        |[v^r] E_k(-1+v)| <= sum_i C(a_k,i) C(i,r) = C(a_k,r) 2^{a_k-r},
    a necessary condition for the profile (a_k) is

        C(m-1, d)  <=  sum_{k : 0 <= d-L_k <= a_k} C(a_k, d-L_k) 2^{a_k-d+L_k}
                                                        for every d.   (ARCH)

    (ARCH) is monotone in each a_k, costs O(m^2), and already beats the
    classical bound: it refutes ratio 3/2 at m = 8, ratio 14/9 at m = 16, and
    ratio 64/41 at m = 32, whereas THM-2160 S6.2 only gives 3/2 for all m.
"""
import sys
from fractions import Fraction
from math import comb

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(1000000)


def arch_ok(m, a):
    """(ARCH): returns (ok, first failing d)."""
    L = [m - 1 - k - a[k] for k in range(m)]
    for d in range(m):
        lhs = comb(m - 1, d)
        rhs = 0
        for k in range(m):
            r = d - L[k]
            if 0 <= r <= a[k]:
                rhs += comb(a[k], r) << (a[k] - r)
            if rhs >= lhs:
                break
        if lhs > rhs:
            return False, d
    return True, None


def profile(m, C):
    a = []
    for k in range(m):
        n = m + k
        cap = Fraction(C) * n - n - 1
        if cap < 0:
            return None
        a.append(min(m - 1 - k, int(cap)))
    return a


def lower_bound(m):
    """Largest candidate ratio refuted by (ARCH); rho(m) is strictly above it."""
    cands = sorted({Fraction(m + k + 1 + aa, m + k)
                    for k in range(m) for aa in range(m - k)})
    lo, hi, best = 0, len(cands) - 1, None
    while lo <= hi:                       # (ARCH) is monotone in C
        mid = (lo + hi) // 2
        a = profile(m, cands[mid])
        ok = False if a is None else arch_ok(m, a)[0]
        if ok:
            hi = mid - 1
        else:
            best = cands[mid]
            lo = mid + 1
    return best


if __name__ == "__main__":
    ms = [int(x) for x in sys.argv[1:]] or [2 ** r for r in range(2, 13)]
    print("   m    certified rho(m) >        (2 - bound)")
    for m in ms:
        b = lower_bound(m)
        if b is None:
            print(f"{m:5d}    (nothing refuted)")
        else:
            print(f"{m:5d}    {str(b):>18s} = {float(b):.6f}   {float(2-b):.6f}")
        sys.stdout.flush()
