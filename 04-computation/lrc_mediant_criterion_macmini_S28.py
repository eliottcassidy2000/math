#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S28 -- HYP-4572: THE MEDIANT-ATTAINER ACHIEVABILITY CRITERION.

opus-S118 cracked the construction I kept missing (4 sessions): the crux mediant
3/(3N+2) attainer is the BORDERED-AP family

        F(N) = {1,2,...,N} \ {N-1}  u  {3(N-1)}
             = {1,...,N-2, N, 3(N-1)}       (N speeds)

- far element  3(N-1)
- binding pair {5, 3(N-1)} sums to q = 3N+2 (antipodal residues +-3 at maximizer)
- verified: N=7 -> {1,2,3,4,5,7,18}; N=13 -> {1..11,13,36}.

opus's correction: NOT primality (N=5, q=17 prime, FAILS); lead N=1 mod 6.

GOAL:
 (A) build F(N) for N=5..25, compute exact M, test M == 3/(3N+2) (achievable);
 (B) find the EXACT criterion (which N) + confirm/refute N=1 mod 6;
 (C) DERIVE it: maximizer t=a/q with 5a==3 mod q; the base residues {v*a mod q}
     must avoid {0,+-1,+-2}; a = 3*inv5 mod q -- report the residue picture and
     the FIRST base element that collides (the obstruction) at failing N;
 (D) pin the N=12 obstruction (q=38=2*19): show WHICH base element lands in
     {0,+-1,+-2} mod 38 and tie it to the 2*19 factorization.
"""
import itertools
from fractions import Fraction as F
from math import gcd


def exact_M(W):
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    for v in W:
        dens.add(2 * v)
    best = F(0)
    seen = set()
    for s in dens:
        if s == 0:
            continue
        for j in range(1, s):
            t = F(j, s)
            if t in seen:
                continue
            seen.add(t)
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
    return best


def primitive(W):
    g = 0
    for v in W:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in W)) if g > 1 else tuple(sorted(W))


def family(N):
    """opus-S118 bordered-AP mediant attainer: {1..N}\{N-1} u {3(N-1)}."""
    base = [v for v in range(1, N + 1) if v != N - 1]
    return primitive(tuple(sorted(set(base) | {3 * (N - 1)})))


def dist_to_int(x):
    """||x|| = distance to nearest integer, x a Fraction."""
    r = x - int(x)
    if r < 0:
        r += 1
    return min(r, 1 - r)


def residue_picture(N):
    """At q=3N+2, a = 3*inv(5) mod q (from binding 5a==3). Return the residues
    v*a mod q for each base speed v, flag those within 2 of 0 (mod q)."""
    q = 3 * N + 2
    if gcd(5, q) != 1:
        return q, None, None  # 5 not invertible (q divisible by 5)
    inv5 = pow(5, -1, q)
    a = (3 * inv5) % q
    W = family(N)
    picture = []
    bad = []
    for v in W:
        r = (v * a) % q
        d = min(r, q - r)  # distance to 0 mod q
        picture.append((v, r, d))
        if d < 3:  # lands in forbidden {0,+-1,+-2}
            bad.append((v, r, d))
    return q, a, (picture, bad)


def main():
    print("=" * 84)
    print("MEDIANT ATTAINER  F(N) = {1..N}\\{N-1} u {3(N-1)}  -- achievability of 3/(3N+2)")
    print("=" * 84)
    print(f"  {'N':>3} {'q=3N+2':>7} {'q factor':>12} {'F(N) (primitive)':>34}")
    print(f"      {'M(F(N))':>10} {'3/q':>8} {'ACHIEVE?':>9} {'N mod 6':>8}")
    achievable = []
    for N in range(5, 26):
        q = 3 * N + 2
        W = family(N)
        # factor q
        fac = []
        m = q
        d = 2
        while d * d <= m:
            while m % d == 0:
                fac.append(d)
                m //= d
            d += 1
        if m > 1:
            fac.append(m)
        facstr = "*".join(map(str, fac)) + (" (PRIME)" if len(fac) == 1 else "")
        M = exact_M(W)
        target = F(3, q)
        ach = (M == target)
        if ach:
            achievable.append(N)
        Wshort = str(list(W)) if len(W) <= 16 else str(list(W)[:15]) + "..]"
        print(f"  {N:>3} {q:>7} {facstr:>12} {Wshort:>34}")
        print(f"      {str(M):>10} {str(target):>8} {('YES' if ach else 'no'):>9} {N % 6:>8}")
    print()
    print(f"  ACHIEVABLE N: {achievable}")
    print(f"  N mod 6 of achievable: {sorted(set(N % 6 for N in achievable))}")
    # check N=1 mod 6
    n1mod6 = [N for N in range(5, 26) if N % 6 == 1]
    print(f"  N==1 mod 6 in range: {n1mod6}")
    print(f"  criterion 'achievable <=> N==1 mod 6' holds: "
          f"{sorted(achievable) == sorted(n1mod6)}")

    print()
    print("=" * 84)
    print("RESIDUE PICTURE at maximizer a=3*inv5 mod q -- WHY each N (fails/works)")
    print("=" * 84)
    for N in [5, 7, 9, 11, 12, 13, 15, 19]:
        q, a, res = residue_picture(N)
        if res is None:
            print(f"  N={N:>2} q={q}: 5 | q, binding pair degenerate (skip)")
            continue
        picture, bad = res
        status = "ACHIEVE (all base residues clear {0,+-1,+-2})" if not bad \
            else f"FAIL: {len(bad)} base speed(s) land in forbidden zone"
        print(f"  N={N:>2} q={q:>3} a={a:>3}: {status}")
        if bad:
            for v, r, d in bad:
                print(f"        speed {v:>3} -> residue {r:>3} mod {q} (dist {d} < 3)  "
                      f"[collision: v*a == {r}]")

    print()
    print("=" * 84)
    print("THE N=12 CRUX  (q=38=2*19)  -- the obstruction in detail")
    print("=" * 84)
    q, a, res = residue_picture(12)
    if res:
        picture, bad = res
        print(f"  q=38=2*19, a=3*inv5 mod 38 = {a}")
        print(f"  full residue picture (speed -> v*a mod 38, dist to 0):")
        line = "    "
        for v, r, d in picture:
            tag = "  <<BAD" if d < 3 else ""
            line = f"    {v:>3} -> {r:>3} (d={d}){tag}"
            print(line)
        print(f"  => {len(bad)} collision(s): the mediant 3/38 is NOT attained by F(12).")
        # tie to 2*19: look at residues mod 2 and mod 19
        print()
        print("  CRT view (q=38=2*19): residue r <-> (r mod 2, r mod 19).")
        print("  forbidden zone {0,+-1,+-2} mod 38 = {0,1,2,36,37}.")
        for v, r, d in bad:
            print(f"     bad speed {v}: r={r} = (mod2 {r%2}, mod19 {r%19}); "
                  f"the collision is at the {'2' if r%19 not in (0,1,2,17,18) else '19'}-part")


if __name__ == "__main__":
    main()
