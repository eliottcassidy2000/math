#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S25 -- HYP-4542: THE RECURSIVE-ACROSS-n PATTERN toward the
HEIGHT UPPER BOUND (opus-S114's reductive target).

The second gap (1/n, 2/(2n-1)) is n-specific (nonempty n=7:5/33, n=8:3/23;
empty n=13).  Does the ACHIEVABLE gap fraction's numerator c CAP and the min
gap-member HEIGHT stay bounded (=> q<=13c bounded => finite check via residue
bridge)?  Search each n=6..13 for gap members (2D GAPs + lifts + random),
validate vs the known members, track: count, min M, max c, height range,
Stern-Brocot depth.
"""
import itertools, random
from fractions import Fraction as F
from math import gcd

def float_M(W):
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    for v in W:
        dens.add(2 * v)
    best = 0.0
    for s in dens:
        if s == 0:
            continue
        for j in range(1, s):
            t = j / s
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
    return best

def exact_M(W):
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    for v in W:
        dens.add(2 * v)
    best = F(0); seen = set()
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

def candidates(n, seed=0):
    """generate likely gap-member families for LRC(n): n-1 speeds."""
    k = n - 1
    rng = random.Random(seed)
    out = set()
    AP = list(range(1, k + 1))
    q1, q2 = n, 2 * n - 1
    # 2D generalized APs {a + i*d1 + j*d2}
    for d1 in range(1, 2 * n):
        for d2 in range(1, n):
            for L1 in range(2, k):
                L2 = -(-k // L1)  # ceil
                elts = set()
                for i in range(L1):
                    for j in range(L2):
                        elts.add(1 + i * d1 + j * d2)
                        if len(elts) >= k:
                            break
                    if len(elts) >= k:
                        break
                W = primitive(tuple(sorted(elts)))
                if len(W) == k:
                    out.add(W)
    # lifts of the AP {1..k}
    for _ in range(15000):
        W = list(AP)
        for _i in range(rng.randint(1, 3)):
            W[rng.randrange(k)] += rng.choice([q1, q2, 2*q1, q1+q2]) * rng.randint(1, 2)
        out.add(primitive(tuple(sorted(set(W)))))
    # random low-height
    for _ in range(15000):
        W = primitive(tuple(sorted(rng.sample(range(1, 3 * n + 6), k))))
        if len(W) == k:
            out.add(W)
    return [W for W in out if len(W) == k]

KNOWN = {7: (1,5,6,11,16,17), 8: (1,2,3,4,5,7,18)}

def main():
    print("=" * 84)
    print("RECURSIVE-ACROSS-n: achievable gap fractions & height, per n")
    print("=" * 84)
    print(f"  {'n':>3} {'gap':>18} {'#members':>9} {'min M (=2nd val)':>17} "
          f"{'max c':>6} {'height range':>14} {'known?':>7}")
    for n in range(6, 14):
        k = n - 1
        lo, hi = F(1, n), F(2, 2 * n - 1)
        flo, fhi = float(lo), float(hi)
        fams = candidates(n, seed=n)
        # ensure known member is in the pool (validation)
        if n in KNOWN:
            fams.append(KNOWN[n])
        members = []
        for W in fams:
            fm = float_M(W)
            if flo - 1e-6 < fm < fhi + 1e-6:
                M = exact_M(W)
                if lo < M < hi:
                    members.append((M, W))
        members = sorted(set(members))
        if members:
            minM = members[0][0]
            maxc = max(M.numerator for M, W in members)
            heights = [max(W) for M, W in members]
            hr = f"{min(heights)}..{max(heights)}"
            known_ok = "YES" if (n not in KNOWN or any(W == KNOWN[n] for M, W in members)) else "MISS"
        else:
            minM = None; maxc = 0; hr = "--"; known_ok = "n/a" if n not in KNOWN else "MISS!"
        print(f"  {n:>3} {f'(1/{n},2/{2*n-1})':>18} {len(members):>9} "
              f"{(str(minM) if minM else 'EMPTY'):>17} {maxc:>6} {hr:>14} {known_ok:>7}")
    print()
    print("  READING: max c = largest achievable numerator (Stern-Brocot depth).")
    print("  If max c CAPS (bounded) and DECREASES to 0 (empty) with n, and height")
    print("  stays ~linear, then a gap member's value denom q<=13c is bounded =>")
    print("  the residue-bridge census is FINITE = the height upper bound.")
    print("  (validation: n=7,8 must show known members found, not MISS.)")

if __name__ == "__main__":
    main()
