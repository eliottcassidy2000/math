# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont44: THE THREE-GAP DISJUNCTION for the AP corner of the divisor-complete residual.
#
# CONTEXT (fleet): the residual (divisor-complete => M>1/14) splits into
#   [near-AP corner longest-AP>=k-5=8: CLOSED by LEM-012 klein-S196 gap-split]
#   [spread bulk longest-AP<=7 (99%, S237): opus-S238 shows needs full irreducible window disjunction].
# opus-S236: for an AP the residues {(a+jd)p mod q} are THEMSELVES an AP mod q (step dp) => clearing =
#   "an AP on Z/q avoiding the danger arc {0,+-1}" = a THREE-GAP / Steinhaus statement (the structured,
#   provable regime). THIS pins the three-gap disjunction EXPLICITLY in the modulus-q (speed-family) framing.
#
# THE CONSECUTIVE-REDUCTION LEMMA (the extreme three-gap case, elementary + Lean-ready):
#   For q with gcd(d,q)=1, the multiplier p = d^{-1} mod q sends the 13-term AP v_i=a+i*d to the residues
#   {a*d^{-1} + i mod q : i=0..12} = 13 CONSECUTIVE residues, a block of length 13. For q in [16,28] the
#   danger arc is exactly {0,+-1}={q-1,0,1} (ceil(q/14)=2, safe band [2,q-2]). A 13-consecutive block
#   {r,...,r+12} (r=a*d^{-1} mod q) avoids {q-1,0,1}  <=>  it does not wrap AND r in [2, q-14].
#   Symmetric multiplier p=-d^{-1} gives the reflected block, clearing <=> r in [14, q-2].
#   => AP clears at q via alpha=+-1  <=>  r=a*d^{-1} mod q in [2,q-14] U [14,q-2].
# THE DISJUNCTION: clear at SOME q in a bounded window [16,Q]. This script (1) validates the lemma exactly,
# (2) tests the alpha=+-1 disjunction over divisor-complete APs, (3) reports the minimal window Q needed,
# (4) sanity-checks that the FULL three-gap (all multipliers p) is never worse than alpha=+-1.
from math import gcd
from functools import reduce
from itertools import product

def safe_band(q):            # danger = residues r with NOT( ceil(q/14) <= r <= q-ceil(q/14) )
    lo = -(-q // 14)         # ceil(q/14)
    return lo, q - lo        # safe = [lo, q-lo]
def clears_full(v, q):       # true clearing: some multiplier p makes ALL residues safe
    lo, hi = safe_band(q)
    return any(all(lo <= (vi * p) % q <= hi for vi in v) for p in range(1, q))
def clears_consec(a, d, q):  # alpha=+-1 (consecutive) reduction certificate
    if gcd(d, q) != 1: return None
    dinv = pow(d, -1, q)
    r = (a * dinv) % q
    lo, hi = safe_band(q)            # for q in [16,28]: lo=2, hi=q-2
    # block {r,...,r+12} avoids danger  <=>  r in [lo, q-1-12-... ]; compute directly & exactly:
    def block_ok(rr):
        res = [(rr + i) % q for i in range(13)]
        return all(lo <= x <= hi for x in res)
    if block_ok(r): return ('+', r)
    if block_ok((-r) % q): return ('-', (-r) % q)  # p = -d^{-1}
    return None
def is_DC(v):  return all(any(x % D == 0 for x in v) for D in range(2, 15))
def prim(v):   return reduce(gcd, v) == 1

def main():
    # (1) EXACT validation of the block-fit formula r in [2,q-14] (alpha=+, q in [16,28], danger {0,+-1})
    print("(1) block-fit lemma check (q in [16,28], danger={0,+-1}):  block {r..r+12} safe  <=>  r in [2,q-14]")
    ok = True
    for q in range(16, 29):
        lo, hi = safe_band(q)
        if (lo, hi) != (2, q - 2): ok = False
        for r in range(q):
            res = [(r + i) % q for i in range(13)]
            block_ok = all(lo <= x <= hi for x in res)
            formula  = (2 <= r <= q - 14)
            if block_ok != formula: ok = False; print("   MISMATCH q,r=", q, r)
    print(f"    formula r in [2,q-14] matches exact block-fit for ALL q in [16,28], all r: {ok}")

    # (2) the alpha=+-1 disjunction over divisor-complete 13-term APs
    aps = []
    for d in range(1, 60):
        for a in range(1, 200):
            v = [a + i * d for i in range(13)]
            if prim(v) and is_DC(v):
                aps.append((a, d, v))
    print(f"\n(2) primitive divisor-complete 13-term APs (d<=59,a<=199): {len(aps)}")
    worstQ = 0; fails = 0; ex = None
    for a, d, v in aps:
        got = None
        for q in range(16, 40):
            c = clears_consec(a, d, q)
            if c is not None: got = (q, c); break
        if got is None:
            fails += 1
            if ex is None: ex = (a, d, v)
        else:
            if got[0] > worstQ: worstQ = got[0]
    print(f"    alpha=+-1 (consecutive) disjunction clears via window [16,{worstQ}];  FAILURES = {fails}")
    if ex: print("    example alpha=+-1 failure:", ex[:2], "(would need full three-gap, alpha!=+-1)")

    # (3) does FULL clearing ever need a q that alpha=+-1 misses in the window?  (three-gap >= consecutive)
    gap = 0
    for a, d, v in aps:
        qc = next((q for q in range(15, 40) if q % 14 and clears_full(v, q)), None)
        qs = next((q for q in range(16, 40) if clears_consec(a, d, q) is not None), None)
        if qc is not None and qs is not None and qs > qc + 0: gap = max(gap, qs - qc)
    print(f"    max( first-consec-q  -  first-full-clear-q ) over corner = {gap}  (0 => consecutive is as early as full)")

    # (4) the tight showcase: {2..14} -- 13 consecutive fit EXACTLY in [2,14] mod 16
    a, d = 2, 1
    print(f"\n(4) showcase {{2..14}}: consecutive cert at q=16 =>", clears_consec(2, 1, 16),
          "(residues {2..14} = safe band [2,14] mod 16, exact fit)")

if __name__ == "__main__":
    main()
