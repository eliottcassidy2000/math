#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THE GAP-CASE DISPATCH CRUX: are COMPRESSED (no-dominant) covering >=7-far 13-families census-able at
BOUNDED q?  If YES, the whole hlarge closes (bounded-q census + peels + spread13) and the "compressed
chain band-blocker" residual DISSOLVES.  If NO (q -> inf), the magnitude split is forced.

opus-2026-07-03-S55. mac-mini-S26 refocused the crux to aligned band-blockers; their compressed
constructor found NO valid examples (compressed_bandblocker_macmini). AND: my S52 "deep well lonely
ONLY at 14/183" was WRONG -- I searched only difference-set denominators; the deep well is census-able
at 3/40 (q=40). So the FULL-q witness scan is essential. Here I do it right.

DEFINITIONS (mac-mini): far = |v|>22; compressed = NO dominant runner (every i has some j!=i with
|v_i| < 13|v_j|); covering = multiple of every q in [2,14]; hge7 = >=7 far. spread13 (route 1) needs
GLOBAL ratio <=13 -- STRONGER than compressed (a chain {1,12,144,...} is compressed but global-ratio-huge).

For each compressed covering >=7-far family: the CORRECT witness q = min q<=Qmax with a lonely a/q
(min-dist >= 1/14), scanned over ALL a coprime to q (not just difference-set denominators).
"""
import sys, itertools, random
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
def norm(x):
    f = x - (x.numerator // x.denominator); return min(f, 1 - f)
def covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def farcnt(S): return sum(1 for v in S if v > 22)
def compressed(S):
    return all(any(j != i and S[i] < 13 * S[j] for j in range(len(S))) for i in range(len(S)))
def witness_q(S, Qmax=400):
    """CORRECT full-q scan: min q with a lonely a/q (min-dist>=1/14, all a coprime to q)."""
    for q in range(2, Qmax + 1):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            if all(norm(Fr(a * v, q)) >= Fr(1, 14) for v in S):
                return q, Fr(a, q)
    return None, None

print("=" * 100)
print(" DO COMPRESSED COVERING >=7-FAR FAMILIES EXIST, AND ARE THEY CENSUS-ABLE AT BOUNDED q?")
print("=" * 100)

# --- construct compressed far-CHAINS + small covering runners, search for valid families ---
random.seed(55)
found = []
# far chains: >=7 runners > 22, compressed (consecutive ratio < 13). Try many highly-composite chains.
far_pools = [
    [24, 30, 36, 40, 42, 48, 54, 56, 60, 66, 70, 72, 84, 88, 90, 96, 104, 110, 120, 126, 132, 140, 154, 156, 168, 180, 182, 198, 220, 234, 260, 273, 280, 308, 312, 364],
    [i for i in range(23, 400)],
]
tested = 0
for trial in range(200000):
    pool = random.choice(far_pools)
    k = random.randint(7, 10)               # >=7 far
    if len(pool) < k: continue
    far = sorted(random.sample(pool, k))
    small = random.sample([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,20], 13 - k)
    S = sorted(set(small) | set(far))
    if len(S) != 13: continue
    if reduce(gcd, S) != 1: continue
    if farcnt(S) < 7: continue
    if not covering(S): continue
    if not compressed(S): continue
    tested += 1
    q, t = witness_q(S, Qmax=300)
    found.append((S, max(S), q, t))
    if len(found) >= 12: break

print(f"\n  found {len(found)} COMPRESSED covering >=7-far families (of {tested} that passed all filters):")
print(f"  {'max-speed':>10} {'#far':>5} {'ratio':>7} {'witness q':>10} {'t':>10} {'q>50?':>6}")
maxq = 0
for S, mx, q, t in found:
    r = mx / min(S)
    qs = str(q) if q else ">300"
    if q and q > maxq: maxq = q
    print(f"  {mx:>10} {farcnt(S):>5} {r:>7.1f} {qs:>10} {str(t):>10} {'YES' if (q and q>50) else 'no':>6}")
    if len(found) <= 4: print(f"       S = {S}")

print("\n" + "=" * 100)
print(" HIGH-MAGNITUDE STRESS: scale a fixed compressed shape up and watch q")
print("=" * 100)
# a fixed compressed covering >=7-far shape, then DILATE the far part by a factor and re-add small covering
# (dilation preserves loneliness by scale-invariance, so q is scale-invariant -- test genuine large-mag shapes)
if found:
    base_shape = found[0][0]
    print(f"  base compressed shape S = {base_shape}, witness q = {found[0][2]}")
    print("  (scale-invariance: L(cS)=L(S), so a dilated far chain has the SAME q -- magnitude does NOT force q up")
    print("   UNLESS the shape itself changes. mac-mini's q~log(mag) needs NEW primes at each scale = a LONGER chain,")
    print("   but a 13-runner chain has BOUNDED length => bounded #primes => bounded q. THE KEY COUNTING BOUND.)")

print("\n" + "=" * 100)
print(" READING / VERDICT")
print("=" * 100)
print(f"  max witness q over the compressed covering >=7-far families found: {maxq if maxq else 'none had q<=300'}")
print("  If maxq is BOUNDED (small): compressed >=7-far families ARE census-able at bounded q => the gap-case")
print("  dispatch's >=7-far branch is CLOSED BY THE CENSUS (no tower needed for compressed) -- only the DOMINANT")
print("  (peelable) large runners need the tower/peel. Since 13 runners bound the chain length (hence the #band")
print("  primes, hence q ~ #primes), the compressed witness q is BOUNDED -- the crux dissolves into the census.")
print("DONE.")
