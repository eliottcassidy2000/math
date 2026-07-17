#!/usr/bin/env python3
"""Recon for HYP-7265 (death-star-S52): the two-circle deep certificate.

Canonical family v = (1..13).  Claim: bandCount(p) >= 6  <=>
  (I)  exists s in {0,1}: 14*6*|p - s*q| < q      [integer circle]
  (II) 14*6*|2p - q| < q                           [half circle]
Check exactness at every q in a range (all p), boundary widths, and the
exact deep-count formula that follows.
"""
def bandcount(q, p):
    c = 0
    for v in range(1, 14):
        r = (v * p) % q
        if 14 * r < q or 14 * r > 13 * q: c += 1
    return c

mismatch = []
qtested = 0
for q in range(15, 1200):
    qtested += 1
    for p in range(1, q):
        deep = bandcount(q, p) >= 6
        circI = (14 * 6 * p < q) or (14 * 6 * (q - p) < q)
        circII = 14 * 6 * abs(2 * p - q) < q
        cert = circI or circII
        if deep != cert:
            mismatch.append((q, p, bandcount(q, p), circI, circII))
            if len(mismatch) <= 8:
                print("MISMATCH", q, p, "bc =", bandcount(q, p),
                      "circI", circI, "circII", circII)
print(f"two-circle certificate over q in [15,1200): {qtested} moduli,"
      f" mismatches = {len(mismatch)}", "PASS" if not mismatch else "FAIL")

# exact deep-count formula if certificate holds:
# N_deep = #circI + #circII - #overlap
# circI: p <= floor((q-1)/84) or p >= q - floor((q-1)/84): 2*floor((q-1)/84)
# circII: |2p - q| < q/84: p in ((q - q/84)/2, (q + q/84)/2): count parity-dep
def direct_deep(q):
    return sum(1 for p in range(1, q) if bandcount(q, p) >= 6)
def formula(q):
    b = (q - 1) // 84
    nI = 2 * b
    nII = sum(1 for p in range(1, q) if 84 * abs(2 * p - q) < q)
    # overlap: both circles: |2p-q| < q/84 and p < q/84-ish: needs q tiny
    nboth = sum(1 for p in range(1, q)
                if 84 * abs(2 * p - q) < q and (84 * p < q or 84 * (q - p) < q))
    return nI + nII - nboth
bad = [q for q in range(15, 1200) if direct_deep(q) != formula(q)]
print(f"exact deep-count identity: mismatched q = {bad[:6]}",
      "PASS" if not bad else "FAIL")
# the half-circle count closed form: #{p : |2p-q| < q/84} = distinguish parity
import math
def nII_closed(q):
    # |2p - q|: odd q -> odd values 1,3,..; even q -> even values 0,2,..
    # 84|2p-q| < q <=> |2p-q| <= ceil(q/84) - 1
    B = (q - 1) // 84  # max allowed |2p-q| is <= B' where 84*x < q <=> x <= (q-1)//84
    if q % 2 == 1:
        # 2p - q odd: values +-1, +-3, ... <= B
        return 2 * ((B + 1) // 2)
    else:
        # values even: 0, +-2, ... <= B
        return 1 + 2 * (B // 2)
bad2 = [q for q in range(15, 1200)
        if sum(1 for p in range(1, q) if 84 * abs(2 * p - q) < q) != nII_closed(q)]
print(f"half-circle closed form: mismatched q = {bad2[:6]}",
      "PASS" if not bad2 else "FAIL")
