#!/usr/bin/env python3
"""
monad-explorer-2026-06-06-S701 (part b)
Full single-runner deletion census of AP_n: for each a in {1,...,n-1}, compute
M(AP_n \ {a}).  Identifies EXACTLY which single deletions raise the gap above 1/n
and by how much.  Tests the S701 claim that the self-conjugate runner n/2 is the
UNIQUE maximal "guard": its removal alone doubles M to 2/n.

Also verifies the elementary proof's two halves:
  (LB) at t = 1/(n/2), every runner of AP_n\{n/2} is at distance >= 2/n;
  (UB) the even subset {2,4,...,n-2} = 2*AP_{n/2} is a sub-config with M = 2/n,
       so M(AP_n\{n/2}) <= 2/n.  Hence equality.
"""
from fractions import Fraction
from math import gcd


def frac(x):
    r = x % 1
    return min(r, 1 - r)


def candidate_times(V):
    c = set()
    for v in V:
        v = abs(v)
        for k in range(2 * v):
            c.add(Fraction(2 * k + 1, 2 * v) % 1)
    for i in range(len(V)):
        for j in range(len(V)):
            for s in (1, -1):
                d = V[i] + s * V[j]
                if d:
                    d = abs(d)
                    for k in range(d + 1):
                        c.add(Fraction(k, d) % 1)
    c.discard(Fraction(0))
    return c


def M_exact(V):
    V = list(V)
    best = Fraction(0)
    bt = None
    for t in candidate_times(V):
        mn = min(frac(v * t) for v in V)
        if mn > best:
            best, bt = mn, t
    return best, bt


if __name__ == "__main__":
    print("=" * 76)
    print("S701b: single-runner deletion census of AP_n (which deletions raise M?)")
    print("=" * 76)
    for n in range(4, 17, 2):
        AP = list(range(1, n))
        base = Fraction(1, n)
        rows = []
        for a in AP:
            R = [v for v in AP if v != a]
            M, t = M_exact(R)
            rows.append((a, M, t))
        # the ones that strictly raise M
        raised = [(a, M, t) for (a, M, t) in rows if M > base]
        print(f"\n--- n={n}  (1/n={base}, 2/n={Fraction(2,n)}, n/2={n//2}) ---")
        for (a, M, t) in rows:
            tag = ""
            if a == n // 2:
                tag = "  <== self-conjugate n/2"
            flag = f"  RAISED to {M}={float(M):.4f} @t={t}" if M > base else "  (stays 1/n)"
            print(f"  delete a={a:>2}: M={M}{flag}{tag}")
        maxM = max(M for (_, M, _) in rows)
        argmax_a = [a for (a, M, _) in rows if M == maxM]
        print(f"  >> max gap after a single deletion: {maxM}={float(maxM):.5f}"
              f" at a={argmax_a}  (2/n={Fraction(2,n)})")
        print(f"  >> n/2 is the UNIQUE maximal guard? "
              f"{argmax_a == [n//2] and maxM == Fraction(2,n)}")

    print()
    print("=" * 76)
    print("Elementary proof check: M(AP_n\\{n/2}) = 2/n via even sub-AP skeleton")
    print("=" * 76)
    for n in range(4, 21, 2):
        m = n // 2
        R = [v for v in range(1, n) if v != m]
        evens = [2 * i for i in range(1, m)]        # {2,4,...,n-2} = 2*AP_m
        # LB witness t=1/m: min distance over R
        t = Fraction(1, m)
        lb = min(frac(v * t) for v in R)
        # M of the even skeleton (should equal M(AP_m) = 1/m = 2/n)
        Me, _ = M_exact(evens)
        Map_m, _ = M_exact(list(range(1, m)))
        MR, tR = M_exact(R)
        ok = (lb == Fraction(2, n)) and (Me == Fraction(2, n)) and (MR == Fraction(2, n)) \
             and (Map_m == Fraction(1, m))
        print(f"  n={n:>2} m={m:>2}: minDist(R,1/m)={lb}={float(lb):.4f}  "
              f"M(evens=2*AP_m)={Me}  M(AP_m)={Map_m}  M(R)={MR}@{tR}  2/n={Fraction(2,n)}  "
              f"{'OK' if ok else 'FAIL'}")
