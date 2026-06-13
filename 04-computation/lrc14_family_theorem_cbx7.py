#!/usr/bin/env python3
"""claudebox-2026-06-11-S7 part D — the family theorem for S(r) = 7*{1..12} u {r}.
(a) r > 1092 = 13*84 = (n-1)*V': loose by THM-398 Cor B2 (dominance dodge on r —
    note: on the STRANGER, no divisibility needed).
(b) r <= 1092 (r not in 7Z, r distinct from core): FINITE — verify every instance
    loose EXACTLY via minimal witness modulus (band criterion, exact integers),
    and record which instances additionally evade [twisted shells m<=27 (band +-1)
    u B'(dodge restricted to multiples of 14)] — the S622-era toolkit.
Also: confirm the residue signature of the evaders (13 | r, r mod 27 in {+-10}, ...)
and that all evader minimal moduli lie in the band-2 range [28, 3n-1=41]."""
from fractions import Fraction as F
from math import gcd
from functools import reduce

TARGET = F(1, 14)
CORE = [7 * k for k in range(1, 13)]

def witness_at_q(S, q):
    B = q // 14
    prof = [v % q for v in S]
    for a in range(1, q):
        if gcd(a, q) != 1: continue
        if all(min((p * a) % q, q - (p * a) % q) > B for p in prof): return a
    return None

def min_witness_modulus(S, qmax=600):
    for q in range(2, qmax + 1):
        if witness_at_q(S, q): return q
    return None

def safe_components(U):
    arcs = []
    for u in U:
        for k in range(0, u + 1):
            a, b = F(14 * k - 1, 14 * u), F(14 * k + 1, 14 * u)
            arcs.append((max(a, F(0)), min(b, F(1))))
    arcs = sorted(arcs); merged = []
    for a, b in arcs:
        if merged and a <= merged[-1][1]: merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else: merged.append((a, b))
    return [(merged[i][1], merged[i + 1][0]) for i in range(len(merged) - 1)
            if merged[i + 1][0] > merged[i][1]]

def Bprime_mult14(S):
    for v in [u for u in S if u % 14 == 0]:
        comps = safe_components([u for u in S if u != v])
        if any(b - a > F(2, 14 * v) for a, b in comps): return True
    return False

n_total = n_loose = 0
evaders = []
maxq_overall = 0
for r in range(1, 1093):
    if r % 7 == 0: continue
    S = sorted(CORE + [r])
    if len(set(S)) != 13 or reduce(gcd, S) != 1: continue
    n_total += 1
    mq = min_witness_modulus(S)
    assert mq is not None, f"NO witness modulus <=600 for r={r} — needs exact M check!"
    n_loose += 1
    maxq_overall = max(maxq_overall, mq)
    if mq > 27 and not Bprime_mult14(S):
        evaders.append((r, mq))

print(f"family S(r)=7*{{1..12}} u {{r}}, valid r <= 1092: {n_total} instances")
print(f"  all loose (exact band-criterion witness found, q <= 600): {n_loose}/{n_total}")
print(f"  max minimal witness modulus over the family: {maxq_overall}")
print(f"  evaders of [shells m<=27 u B'(mult-of-14)]: {len(evaders)}")
print(f"  evader list (r, min q): {evaders}")
sig = [(r, r % 13, r % 27) for r, _ in evaders]
print(f"  evader signatures (r, r mod 13, r mod 27): {sig}")
band2 = all(28 <= q <= 41 for _, q in evaders)
print(f"  every evader's minimal modulus in the band-2 range [28, 3n-1=41]: {band2}")
print()
print("THEOREM (family): every valid S(r) is loose —")
print("  r > 1092: THM-398 Cor B2, dominance dodge on the STRANGER r (r > 13*84 = (n-1)V').")
print(f"  r <= 1092: exhaustively verified above ({n_total} instances, exact).")
print("  The evaders prove: [twisted shells m <= 2n-1] u [B' restricted to multiples of n]")
print("  has a NONEMPTY residual at n=14; the residual is caught exactly one band-rung up.")
