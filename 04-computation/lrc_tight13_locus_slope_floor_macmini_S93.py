#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S93 -- HYP-3840 part 2: EXHAUSTIVE-ish tight-13-set locus mod 14
and the exact slope floor at N=14.

Families enumerated (residue classification HYP-3750: permutation OR duplication+drop):
  (A) single-lift permutation type: {1..13}\{x} u {x+14j}, j=1..12  (residues = perm)
  (B) dup+drop: {1..13}\{v} u {s+14j}, s != v, j=1..12              (drop v, dup s)
  (C) double dup+drop: {1..13}\{v1,v2} u {s1+14j1, s2+14j2}         (j<=4, small sweep)
  (D) dilations: c*S for tight S found, gcd(c,14)=1 -- slope INVARIANT (checked)
Also: q=9 (n=8) full dup+drop sweep to explain klein-S48 census zero via the
unit-residue lemma + the second-order obstruction.

Element cap: 13 + 14*12 = 181 (covers deep-well scale n(n-1)=182 territory).
"""
from fractions import Fraction as F
from math import gcd
import itertools, sys

def dist(x):
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def M_exact(S):
    Sl = sorted(set(S)); dens = set()
    for i, v in enumerate(Sl):
        for w in Sl[i:]:
            dens.add(v + w)
            if w > v: dens.add(w - v)
    best = F(0)
    for den in dens:
        for m in range(1, den):
            t = F(m, den)
            mn = min(dist(v * t) for v in Sl)
            if mn > best:
                best = mn
    return best

def M_at_least(S, thr):
    """Quick check M(S) >= thr: search witness on binding grid only above thr."""
    Sl = sorted(set(S)); dens = set()
    for i, v in enumerate(Sl):
        for w in Sl[i:]:
            dens.add(v + w)
            if w > v: dens.add(w - v)
    for den in dens:
        for m in range(1, den):
            t = F(m, den)
            if min(dist(v * t) for v in Sl) >= thr:
                return True
    return False

def danger_intervals(S, r):
    iv = []
    for v in S:
        half = r / v
        for k in range(v):
            c = F(k, v); a, b = c - half, c + half
            if a < 0: iv.append((a + 1, F(1), v)); iv.append((F(0), b, v))
            elif b > 1: iv.append((a, F(1), v)); iv.append((F(0), b - 1, v))
            else: iv.append((a, b, v))
    return iv

def lonely_measure(S, r):
    iv = sorted(danger_intervals(S, r)); tot = F(0); cur_a = cur_b = None
    for a, b, v in iv:
        if cur_b is None: cur_a, cur_b = a, b
        elif a <= cur_b: cur_b = max(cur_b, b)
        else: tot += cur_b - cur_a; cur_a, cur_b = a, b
    if cur_b is not None: tot += cur_b - cur_a
    return 1 - tot

def breakpoints_below(S, rmax):
    pts = set(); Sl = sorted(set(S))
    for i, v in enumerate(Sl):
        for w in Sl[i:]:
            for den in ([v + w] + ([w - v] if w > v else [])):
                if den <= 0: continue
                d = 1
                while F(d, den) < rmax:
                    pts.add(F(d, den)); d += 1
    return sorted(pts)

def critical_slope(S, N):
    rmax = F(1, N)
    grid = breakpoints_below(S, rmax)
    last = grid[-1] if grid else F(0)
    r1 = last + (rmax - last) * F(1, 3); r2 = last + (rmax - last) * F(2, 3)
    c1 = lonely_measure(S, r1) / (1 - N * r1)
    c2 = lonely_measure(S, r2) / (1 - N * r2)
    assert c1 == c2, f"not linear in last cell for {S}"
    return c1

def cap(q):
    return F(2, q) * sum(F(1, a) for a in range(1, q) if gcd(a, q) == 1)

def sweep(q, jmax, do_double=False, cap_elem=None):
    n = q - 1
    base = set(range(1, q))
    tight = {}
    # (A) permutation single-lift
    for x in range(1, q):
        for j in range(1, jmax + 1):
            S = tuple(sorted((base - {x}) | {x + q * j}))
            if len(S) != n: continue
            if M_exact(S) == F(1, q):
                tight[S] = ("perm-lift", x, x, j)
    # (B) dup+drop
    for v in range(1, q):
        for s in range(1, q):
            if s == v: continue
            for j in range(1, jmax + 1):
                S = tuple(sorted((base - {v}) | {s + q * j}))
                if len(S) != n: continue
                if M_exact(S) == F(1, q):
                    tight[S] = ("dup+drop", v, s, j)
    # (C) double dup+drop (small sweep)
    if do_double:
        for v1, v2 in itertools.combinations(range(1, q), 2):
            for s1 in range(1, q):
                for s2 in range(s1, q):
                    for j1 in range(1, 5):
                        for j2 in range(j1 if s1 == s2 else 1, 5):
                            e1, e2 = s1 + q * j1, s2 + q * j2
                            if e1 == e2: continue
                            S = tuple(sorted((base - {v1, v2}) | {e1, e2}))
                            if len(S) != n: continue
                            # cheap pre-filter: witness at 1/q must survive
                            if min(dist(v * F(1, q)) for v in S) < F(1, q) and \
                               all(min(dist(v * F(a, q)) for v in S) < F(1, q)
                                   for a in range(1, q) if gcd(a, q) == 1):
                                continue
                            if M_exact(S) == F(1, q):
                                tight[S] = ("2xdup+drop", (v1, v2), (s1, s2), (j1, j2))
    return tight

def main():
    print("=" * 78)
    print("TIGHT LOCUS + SLOPE FLOOR, q=14 (13 runners)")
    print("=" * 78)
    q = 14
    tight = sweep(q, jmax=12, do_double=True)
    c_ap = cap(q)
    print(f"C_AP(14) = {c_ap} = {float(c_ap):.6f}")
    print(f"non-AP tight sets found: {len(tight)}")
    slopes = [(critical_slope(list(S), q), S, info) for S, info in sorted(tight.items())]
    slopes.sort()
    for c, S, info in slopes:
        rel = "= C_AP" if c == c_ap else ("< C_AP  <<< BEATER" if c < c_ap else "> C_AP")
        print(f"  slope {c} = {float(c):.6f}  {rel}   {info}   S={list(S)}")
    floor = min([c_ap] + [c for c, _, _ in slopes])
    print(f"\n=> slope floor at N=14 over enumerated tight locus: {floor} = {float(floor):.6f}")

    print()
    print("=" * 78)
    print("WHY NO UNIT-DUP TIGHT SET AT q=14 BUT ONE AT q=8: the second constraint")
    print("=" * 78)
    # For dup+drop tight: drop v (non-unit forced by unit-residue lemma). What forces
    # which s work?  Report the covering failure certificate for a few candidates.
    q = 14
    base = set(range(1, q))
    for (v, s, j) in [(2, 5, 1), (2, 3, 1), (4, 9, 1), (6, 11, 1), (12, 10, 1), (12, 11, 1)]:
        S = sorted((base - {v}) | {s + q * j})
        M = M_exact(S)
        status = "TIGHT" if M == F(1, q) else f"M={M}={float(M):.5f} {'>' if M > F(1,q) else '<'} 1/14"
        print(f"  drop {v}, dup {s} via {s+q*j}: {status}")

    print()
    print("=" * 78)
    print("q=9 (n=8): explain klein-S48 census ZERO at n=8")
    print("=" * 78)
    q = 9
    tight9 = sweep(q, jmax=10, do_double=False)
    print(f"q=9 tight non-AP sets: {len(tight9)}")
    for S, info in tight9.items():
        c = critical_slope(list(S), q)
        print(f"  {info}  slope={c}={float(c):.6f}  S={list(S)}")
    # non-units mod 9: 3,6.  Try all drops in {3,6} explicitly and report failure mode.
    base = set(range(1, q))
    for v in (3, 6):
        for s in range(1, q):
            if s == v: continue
            S = sorted((base - {v}) | {s + 9})
            M = M_exact(S)
            if M >= F(1, 9):
                print(f"  drop {v}, dup {s} via {s+9}: M={M} {'TIGHT' if M==F(1,9) else '>1/9 (not tight, LRC-safe)'}")

    print()
    print("=" * 78)
    print("q=12 and q=16 quick sweeps (even-modulus family: does a beater exist?)")
    print("=" * 78)
    for q in (12, 16):
        t = sweep(q, jmax=6, do_double=False)
        c_ap_q = cap(q)
        print(f"q={q}: C_AP={c_ap_q}={float(c_ap_q):.6f}; non-AP tight found={len(t)}")
        for S, info in sorted(t.items()):
            c = critical_slope(list(S), q)
            rel = "=" if c == c_ap_q else ("<  BEATER" if c < c_ap_q else ">")
            print(f"   slope {float(c):.6f} {rel} C_AP  {info}  S={list(S)}")

if __name__ == "__main__":
    main()
