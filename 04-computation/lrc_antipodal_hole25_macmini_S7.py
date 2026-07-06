#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S7 -- HYP-4312: the antipodal-transversal reduction at the
hdich level (mod 2n-1 = 25) + the NON-UNIT-PAIR HOLE (oracle-S552/S553, integrated).

LEVEL: observer + 12 speeds, gap = (1/13, 2/25], 2n-1 = 25 = 5^2.

LINK-1 (verified here): a covering primitive 12-set W whose residues mod 25 MISS
a UNIT antipodal pair {a, 25-a} (a in (Z/25)*) is SAFE with surplus -- at
t = a^{-1}/25 every residue w a^{-1} mod 25 lands in {2..23} (not 0,+-1), so
M >= 2/25 (ABOVE the gap).  So GAP members hit every unit antipodal pair.

THE HOLE: a=5,10,15,20 are NON-units mod 25; the pairs {5,20},{10,15} have no
Link-1 witness.  A gap member could hit all 10 unit pairs but MISS a non-unit
pair -- the class that produced the n=8 sporadics (2n-1=15).  At 25 = 5^2 the
non-units are 5*(Z/5), so the hole is mod-5 structured.

THIS SCRIPT:
 (1) VERIFY Link-1 exactly: random covering families missing a unit pair -> M >= 2/25.
 (2) ENUMERATE the hole class: residue-12-multisets mod 25 hitting all 10 unit
     pairs but missing a non-unit pair; for each realizable one, find a WITNESS
     (dilation a, modulus q) giving margin >= 2/25 -- the second witness family.
 (3) The mechanism: for missed {5,20}, the residues avoid 5,20 mod 25; test the
     lift witnesses (mod 25 refined, mod 5, the half-value-clock analog).
"""
from math import gcd
from fractions import Fraction as F
from itertools import combinations
import random, sys, time
sys.path.insert(0, '04-computation')
from lonely_profile import profile

T0 = time.time()
def log(m=""):
    print(m, flush=True)
random.seed(7)

Q = 25
def dq(x, q):
    x %= q
    return min(x, q - x)
UNITS = [a for a in range(1, Q) if gcd(a, Q) == 1]
UNIT_PAIRS = sorted({tuple(sorted((a, Q - a))) for a in UNITS})   # 10 pairs
NONUNIT_PAIRS = [(5, 20), (10, 15)]

def M_exact(S):
    for cap in (11, 8, 6, 4, 3, 2):
        p = profile(sorted(S), F(1, cap))
        m = p.M()
        if m is not None:
            return m

def M_lower_scan(W, qmax=200):
    """fast lower bound on M via best rational witness q<=qmax (>= true M is false;
    this is a LOWER bound; for Link-1 we only need to confirm >= 2/25)."""
    best = F(0)
    for q in range(8, qmax + 1):
        if any(w % q == 0 for w in W):
            continue
        for a in range(1, q // 2 + 1):
            if gcd(a, q) != 1:
                continue
            m = min(dq(a * w, q) for w in W)
            if F(m, q) > best:
                best = F(m, q)
    return best

def hits_all_unit_pairs(res):
    rs = set(r % Q for r in res)
    return all(a in rs or b in rs for a, b in UNIT_PAIRS)

def missed_nonunit_pairs(res):
    rs = set(r % Q for r in res)
    return [p for p in NONUNIT_PAIRS if p[0] not in rs and p[1] not in rs]

# ---------- (1) verify Link-1 ----------
log("(1) LINK-1 check: covering primitive 12-sets missing a UNIT pair mod 25 => M >= 2/25")
def covering(W):
    return all(any(w % m == 0 for w in W) for m in range(2, 13))
bad = 0; checked = 0
for _ in range(3000):
    W = sorted(random.sample(range(1, 60), 12))
    if gcd(*W) != 1 or not covering(W):
        continue
    res = [w % Q for w in W]
    # find a missed unit pair
    rs = set(res)
    missed = [(a, b) for a, b in UNIT_PAIRS if a not in rs and b not in rs]
    if not missed:
        continue
    a0, b0 = missed[0]
    ainv = pow(a0, -1, Q)
    # analytic Link-1 witness t = ainv/25: margin = min_i dq(w*ainv,25)/25
    marg = min(dq(w * ainv, Q) for w in W)
    checked += 1
    if F(marg, Q) < F(2, 25):
        bad += 1
        if bad <= 3:
            log(f"  !! MISS-UNIT-PAIR {(a0,b0)} but analytic witness margin {marg}/25 < 2/25: W={W} (residue-0? {any(w%Q==0 for w in W)})")
log(f"  checked {checked} unit-pair-missers; M < 2/25 violations: {bad} "
    f"(0 => Link-1 holds: missing a unit pair => safe)")

# ---------- (2) the hole class ----------
log("\n(2) THE NON-UNIT-PAIR HOLE: hit all 10 unit pairs, miss a non-unit pair")
# enumerate integer covering families and filter to the hole class; find witnesses
def best_witness(W, qmax=125):
    """return (q, a, margin) maximizing margin; margin as Fraction."""
    best = (0, 0, F(0))
    for q in range(8, qmax + 1):
        if any(w % q == 0 for w in W):
            continue
        for a in range(1, q // 2 + 1):
            if gcd(a, q) != 1:
                continue
            m = min(dq(a * w, q) for w in W)
            if F(m, q) > best[2]:
                best = (q, a, F(m, q))
    return best

hole_found = []
n_hole = 0
t0 = time.time()
tries = 0
while time.time() - t0 < 150 and n_hole < 60:
    tries += 1
    W = sorted(random.sample(range(1, 80), 12))
    if gcd(*W) != 1 or not covering(W):
        continue
    res = [w % Q for w in W]
    if not hits_all_unit_pairs(res):
        continue
    miss = missed_nonunit_pairs(res)
    if not miss:
        continue
    n_hole += 1
    M = M_exact(W)
    ingap = F(1, 13) < M < F(2, 25)
    q, a, marg = best_witness(W)
    if len(hole_found) < 12:
        hole_found.append((W, sorted(set(res)), miss, M, (q, a, marg)))
    if ingap:
        log(f"  !! IN-GAP HOLE MEMBER: W={W} M={M}")
log(f"  hole-class members found: {n_hole} (of {tries} tries)")
log(f"  sample (W | missed non-unit pair | M | best witness q,a,margin):")
for W, res, miss, M, (q, a, marg) in hole_found[:10]:
    tag = "IN-GAP!!" if F(1,13) < M < F(2,25) else ("=2/25" if M == F(2,25) else "safe")
    log(f"   miss{miss} M={str(M):>7} ({tag}) witness q={q},a={a},marg={marg}  W={W}")
inwin = [h for h in hole_found if F(1,13) < h[3] < F(2,25)]
log(f"\n  hole-class members IN THE OPEN GAP: {len(inwin)}")

# ---------- (3) the witness modulus pattern ----------
log("\n(3) WITNESS MODULUS for hole members (missing {5,20} or {10,15} mod 25):")
qs = {}
for W, res, miss, M, (q, a, marg) in hole_found:
    qs[q] = qs.get(q, 0) + 1
log(f"  witness modulus histogram: {dict(sorted(qs.items()))}")
log(f"  (if witnesses cluster at a LIFT modulus 25*k or a small q, that is the")
log(f"   second witness family closing the hole -- invertibility restored)")
log(f"\n[t = {time.time()-T0:.0f}s]")
