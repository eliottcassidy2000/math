#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S35 -- PRECISE, assumption-challenging study of the r=2
double-13-lift covering (the (G) residual, branch 2).

r=2 shape (i,j): the AP {1..12} with speed i -> i+13a, speed j -> j+13b.
kps-S47 claims: every non-AP member clears at some q<=25 (Q0=25), height-uniform,
hard shapes lift speed 6/12.  CHALLENGES to test:
 (A) Is (i+13a) mod q independent of a?  NO (verified) -- so clearing q DEPENDS on
     (a,b); "height-uniform" means the covering SET {q<=25} is fixed, but which q
     clears varies.  So uniformity = a covering-system claim over (a mod q, b mod q).
 (B) Does Q0=25 hold over a WIDE (a,b) range (not just kps's [0,25])?  Look for an
     ESCAPE (a,b) with no clearing q<=25.
 (C) Which shapes are genuinely hard (max clearing q), and is speed 6/12 the reason?
 (D) The covering STRUCTURE: for a hard shape, the map (a mod q, b mod q) -> clears?
     Is it a finite decidable check (small lcm) or does it need an Erdos argument?

Clearing at (q,c): mu=ceil(2q/25); all v*c mod q in [mu, q-mu].
"""
import itertools
from math import gcd
from fractions import Fraction as F


def mu_of(q):
    return -(-2 * q // 25)  # ceil(2q/25)


def clears_at(vlist, q):
    """does SOME unit c clear all speeds at modulus q (v*c mod q off band)?  return c or None."""
    mu = mu_of(q)
    if 2 * mu > q:
        return None
    for c in range(1, q):
        if gcd(c, q) != 1:
            continue
        if all(mu <= (v * c) % q <= q - mu for v in vlist):
            return c
    return None


def min_clear_q(vlist, qmax=60):
    for q in range(6, qmax + 1):
        if clears_at(vlist, q) is not None:
            return q
    return None


def shape_family(i, j, a, b):
    base = [v for v in range(1, 13) if v not in (i, j)]
    fam = sorted(base + [i + 13 * a, j + 13 * b])
    if len(set(fam)) != 12 or any(v <= 0 for v in fam):
        return None
    return tuple(fam)


def main():
    print("=" * 88)
    print("(B)+(C): min clearing q over (a,b), per shape -- does Q0=25 hold on a WIDE range?")
    print("=" * 88)
    RNG = range(0, 40)  # wide: a,b to 39 (heights to ~1+13*39=508) -- beyond kps's [0,25]
    worst = {}
    escapes = []
    for i, j in itertools.combinations(range(1, 13), 2):
        mq = 0
        esc = 0
        for a in RNG:
            for b in RNG:
                if a == 0 and b == 0:
                    continue  # the AP
                fam = shape_family(i, j, a, b)
                if fam is None:
                    continue
                q = min_clear_q(fam, qmax=60)
                if q is None:
                    esc += 1
                    escapes.append((i, j, a, b, fam))
                else:
                    mq = max(mq, q)
        worst[(i, j)] = (mq, esc)
    # report the hardest shapes
    ranked = sorted(worst.items(), key=lambda kv: -kv[1][0])
    print(f"  {'shape':>8} {'max clearing q':>15} {'#escapes(q>60)':>15}")
    for (i, j), (mq, esc) in ranked[:14]:
        flag = "  <-- HARD" if mq >= 20 else ""
        print(f"  {str((i,j)):>8} {mq:>15} {esc:>15}{flag}")
    print()
    allmax = max(mq for mq, esc in worst.values())
    totesc = sum(esc for mq, esc in worst.values())
    print(f"  MAX clearing q over ALL shapes and (a,b) in [0,40)^2: {allmax}")
    print(f"  ESCAPES (no clear at q<=60): {totesc}  {'-> Q0<=60 holds' if totesc==0 else '-> INVESTIGATE'}")
    print(f"  => kps Q0=25 {'CONFIRMED wider' if allmax<=25 else 'CHALLENGED: max q='+str(allmax)+' > 25'}")

    print()
    print("=" * 88)
    print("(D): covering STRUCTURE for the hardest shape -- clearing q as fn of (a mod q, b mod q)")
    print("=" * 88)
    ih, jh = ranked[0][0]
    print(f"  hardest shape (i,j)=({ih},{jh}); for each (a,b) in [0,26)^2, the MIN clearing q:")
    # tally which q's are used and whether a small covering set suffices
    used_q = {}
    for a in range(0, 26):
        for b in range(0, 26):
            if a == 0 and b == 0:
                continue
            fam = shape_family(ih, jh, a, b)
            if fam is None:
                continue
            q = min_clear_q(fam)
            used_q[q] = used_q.get(q, 0) + 1
    print(f"  min-clearing-q distribution: {dict(sorted(used_q.items()))}")
    print(f"  distinct moduli needed: {sorted(k for k in used_q if k)}")
    # is clearing periodic in a with period q? test: does (a) -> clears at q depend only on a mod q?
    print()
    print("  PERIODICITY test (is clearing at q a function of a mod q, b mod q?):")
    for q in [17, 19, 23]:
        ok = True
        for a in range(0, q):
            for b in range(0, q):
                f1 = shape_family(ih, jh, a, b)
                f2 = shape_family(ih, jh, a + q, b + q)
                if f1 and f2:
                    c1 = clears_at(f1, q)
                    c2 = clears_at(f2, q)
                    if (c1 is None) != (c2 is None):
                        ok = False
        print(f"    q={q}: clearing at q depends only on (a mod q, b mod q): {ok}")
    print()
    print("  => if periodic, the covering is a FINITE decidable check over (a,b) mod lcm(moduli);")
    print("     lcm of the needed moduli sets the Lean-decidability (small => decide-able).")


if __name__ == "__main__":
    main()
