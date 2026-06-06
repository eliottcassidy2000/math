#!/usr/bin/env python3
"""
monad-explorer-2026-06-06-S702
THE COMPLETE single-runner deletion profile of AP_n, unifying both halves:

    M(AP_n \\ {a}) = | 1/a                        if  n/2 <= a <= n-1     (S701)
                    | 2 / D*(n,a)                 if  1 <= a < n/2        (S702)
      where  D*(n,a) = min { D >= n+a : gcd(D,a) = 1 }.

Verifies:
 (1) the unified closed form against exact M (n=4..40), 0 failures;
 (2) the optimum is always a "c=2 two-binder" optimum on the lower half:
     M * (b1+b2) = 2 with b1+b2 = D* (binders straddle 0 at residues +-2);
 (3) corollaries: argmax_a M = n/2 (the guard, even n; M=2/n) and
     argmin_a M = n-1 (fastest runner, M=1/(n-1)) -> the deletion-criticality
     spread is [1/(n-1), 2/n].
"""
from fractions import Fraction
from math import gcd


def M_exact_int(V):
    mx = max(V)
    bn, bd, bt = 0, 1, None
    for D in range(2, 2 * mx + 1):
        for p in range(1, D):
            c = D
            for v in V:
                r = (v * p) % D
                d = r if r <= D - r else D - r
                if d < c:
                    c = d
                    if c == 0:
                        break
            if c == 0:
                continue
            if c * bd > bn * D:
                bn, bd, bt = c, D, (p, D)
    return Fraction(bn, bd), bt, bn  # gap, (p,D), numerator c at that D


def Dstar(n, a):
    D = n + a
    while gcd(D, a) != 1:
        D += 1
    return D


def closed_form(n, a):
    if 2 * a >= n:
        return Fraction(1, a)
    return Fraction(2, Dstar(n, a))


if __name__ == "__main__":
    print("=" * 78)
    print("S702: UNIFIED single-deletion profile of AP_n  (both halves)")
    print("=" * 78)
    fails = c2_fails = 0
    total = 0
    for n in range(4, 41):
        for a in range(1, n):
            R = [v for v in range(1, n) if v != a]
            M, (p, D), c = M_exact_int(R)
            cf = closed_form(n, a)
            total += 1
            if M != cf:
                fails += 1
                print(f"  FAIL n={n} a={a}: M={M} closed-form={cf}")
            # c=2 check on lower half (two-binder straddle, sum=D*)
            if 2 * a < n:
                binders = [v for v in R
                           if min((v * p) % D, D - (v * p) % D) == c]
                bsum_ok = any(x + y == D for x in binders for y in binders if x < y)
                if not (c == 2 and bsum_ok and D == Dstar(n, a)):
                    c2_fails += 1
                    print(f"  c2 FAIL n={n} a={a}: c={c} D={D} D*={Dstar(n,a)} "
                          f"binders={binders}")
    print(f"\n  unified closed form: {total-fails}/{total} "
          f"({'ALL PASS' if fails==0 else f'{fails} FAILS'})  [n=4..40]")
    print(f"  lower-half c=2 two-binder (b1+b2=D*): "
          f"{'ALL PASS' if c2_fails==0 else f'{c2_fails} FAILS'}")

    print("\n  Deletion-criticality spread (argmax = guard, argmin = load-bearing):")
    for n in range(6, 25):
        prof = {a: closed_form(n, a) for a in range(1, n)}
        amax = max(prof, key=lambda a: prof[a])
        amin = min(prof, key=lambda a: prof[a])
        print(f"   n={n:>2}: max M={prof[amax]} at a={amax}"
              f"{' (=n/2 guard)' if amax==n//2 else ''}; "
              f"min M={prof[amin]} at a={amin}"
              f"{' (=n-1 fastest)' if amin==n-1 else ''}")
