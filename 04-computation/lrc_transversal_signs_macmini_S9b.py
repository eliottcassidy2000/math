#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S9b -- HYP-4352 (cont.): THE MOD-25 TRANSVERSAL SIGN-CHOICE
enumeration -- sharpening the sibling-S7 transversal AP-rigidity.

SYNTHESIS: a gap member is (S9) non-clearing mod 13 AND 25 + (THM-369) contains
a multiple of each m<=12.  The sieve needs a multiple of 5 (non-unit mod 25);
the transversal needs to hit all 10 unit pairs {a,25-a}.  With 12 elements this
forces ~10 units (one per pair, a SIGN CHOICE) + ~2 non-units (multiples of 5).
So a gap member's mod-25 profile is determined UP TO the 2^10 = 1024 sign
choices of which element of each pair.  The AP {1..12} picks the SMALL element
of each pair + {5,10}.

This enumerates all 1024 sign choices (with the two mult-of-5 slots = 5,10) at
CANONICAL lift (smallest positive residues), computes exact M, and asks: is the
AP's sign-choice the ONLY one avoiding a clearing (>=2/25) or does the whole
transversal family sit above the window (so no gap member)?
"""
import itertools
from fractions import Fraction as F
from math import gcd

def dist_int(x, q):
    r = x % q
    return min(r, q - r)

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

def nonclearing_mod(W, q):
    thr = F(2, 25) * q
    residues = [v % q for v in W]
    for a in range(1, q):
        if all(dist_int(r * a, q) >= thr for r in residues):
            return False
    return True

RHO = F(2, 25); FLOOR = F(1, 13)

# the 10 unit pairs mod 25
UNITS = [u for u in range(1, 25) if gcd(u, 25) == 1]
PAIRS = [(u, 25 - u) for u in UNITS if u < 25 - u]   # 10 pairs

def main():
    print("=" * 78)
    print("MOD-25 TRANSVERSAL SIGN-CHOICE ENUMERATION (1024 choices)")
    print("canonical lift: pick a or 25-a per pair; add mult-of-5 slots {5,10}")
    print("=" * 78)
    print(f"  10 unit pairs: {PAIRS}")
    ap_signs = tuple(0 for _ in PAIRS)   # 0 = pick small element (the AP)
    results = {"gap": [], "floor": [], "clear": [], "below": []}
    nonclear_both = []      # non-clearing mod BOTH 13 and 25
    for signs in itertools.product([0, 1], repeat=10):
        chosen = [PAIRS[i][signs[i]] for i in range(10)]   # 10 unit residues
        W = tuple(sorted(chosen + [5, 10]))                 # + the mult-of-5 slots
        if len(set(W)) != 12:
            continue
        M = exact_M(W)
        if FLOOR < M < RHO:
            bucket = "gap"
        elif M == FLOOR:
            bucket = "floor"
        elif M >= RHO:
            bucket = "clear"
        else:
            bucket = "below"
        results[bucket].append((signs, W, M))
        if nonclearing_mod(W, 13) and nonclearing_mod(W, 25):
            nonclear_both.append((signs, W, M))
    total = sum(len(v) for v in results.values())
    print(f"  sign-choices with 12 distinct elements: {total}")
    print(f"  by M-bucket: gap={len(results['gap'])}, floor(1/13)={len(results['floor'])}, "
          f"clear(>=2/25)={len(results['clear'])}, below={len(results['below'])}")
    print()
    # the AP's sign-choice
    ap_entry = [(s, W, M) for s, W, M in results['floor'] + results['clear'] + results['gap'] + results['below']
                if s == ap_signs]
    if ap_entry:
        s, W, M = ap_entry[0]
        print(f"  AP sign-choice (all small): W={list(W)}, M={M} "
              f"({'= 1/13 the tight AP' if M==FLOOR else '?'})")
    print(f"  non-clearing mod BOTH 13 and 25: {len(nonclear_both)}")
    for s, W, M in nonclear_both[:12]:
        tag = " <-- AP" if s == ap_signs else (" IN-GAP!" if FLOOR < M < RHO else "")
        print(f"    signs={s} M={M} W={list(W)}{tag}")
    print()
    if results['gap']:
        print(f"  *** {len(results['gap'])} IN-GAP sign-choices (canonical lift):")
        for s, W, M in results['gap'][:10]:
            print(f"    signs={s} M={M} W={list(W)}")
        print("  NOTE: canonical lift only -- these need the full lift check, but")
        print("  their EXISTENCE at canonical lift flags candidate gap transversals.")
    else:
        print("  ZERO in-gap sign-choices at canonical lift: every mod-25 transversal")
        print("  (with mult-of-5 slots) is the AP (1/13) or clears (>=2/25) -- the")
        print("  transversal AP-rigidity holds at canonical lift.")
    # how many sign-choices are BELOW (M<1/13)?  those violate LRC13 => can't be
    # canonical-lift realizable as a covering family (interesting to see)
    print(f"\n  M-value distribution among sign-choices (canonical lift):")
    dist = {}
    for b in results:
        for s, W, M in results[b]:
            dist[M] = dist.get(M, 0) + 1
    for M in sorted(dist)[:20]:
        print(f"    M={M} = {float(M):.5f}  x{dist[M]}  {'(AP floor)' if M==FLOOR else ('(<gap)' if M<FLOOR else ('(GAP)' if M<RHO else ''))}")

if __name__ == "__main__":
    main()
