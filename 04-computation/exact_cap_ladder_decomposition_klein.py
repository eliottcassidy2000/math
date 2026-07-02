#!/usr/bin/env python3
"""
klein-2026-07-01-S91 -- HYP-3837: THE EXACT CAP LADDER + the clean-branch decomposition
                        + HYP-3838: the d=3 clean-branch criterion (data-driven)

PART A (HYP-3837). cap_j := min over j-subsets P of {1..13} of Lambda_P(1/14)
(the small-part good-set measure; THM-576 identifies cap = pairwise-avoidance for j<=3
and exact constants at j=4,5). Compute the WHOLE ladder j=1..13 exactly (all 2^13-1
subsets), verify THM-576, extend, and give margins + minimizers.

PART B. Decomposition test (mac-mini S96 sec-1 / hp0cap): expand cap_8 and cap_9's
minimizing P via inclusion-exclusion with EXACT d-fold overlap values; classify each
term clean-branch (origin nest 2r/v_max) vs wrapped; exhibit the caps as clean + slice
arithmetic.

PART C (HYP-3838). All triples u<v<w<=20 at r=1/14: exact |D_u^D_v^D_w|; classify
clean (=2r/w); test criterion candidates; list the exact boundary.
"""
from fractions import Fraction as F
import itertools

R = F(1, 14)

def danger(v, r=R):
    rv = r / v
    out = []
    for a in range(v + 1):
        lo, hi = F(a, v) - rv, F(a, v) + rv
        lo, hi = max(lo, F(0)), min(hi, F(1))
        if hi > lo:
            out.append((lo, hi))
    return out

def merge(ivs):
    ivs = sorted(ivs)
    out = []
    for lo, hi in ivs:
        if out and lo <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], hi))
        else:
            out.append((lo, hi))
    return out

def measure(ivs):
    return sum(hi - lo for lo, hi in ivs)

def inter(A, B):
    out = []
    for lo1, hi1 in A:
        for lo2, hi2 in B:
            lo, hi = max(lo1, lo2), min(hi1, hi2)
            if hi > lo:
                out.append((lo, hi))
    return out

DANGERS = {v: merge(danger(v)) for v in range(1, 21)}

def lam(P):
    return 1 - measure(merge([iv for v in P for iv in DANGERS[v]]))

# ---------------- A. the ladder ----------------
print("=" * 96)
print("A. THE EXACT CAP LADDER: cap_j = min over j-subsets P of {1..13} of Lambda_P(1/14)")
print("=" * 96)
KNOWN = {1: F(78, 91), 2: F(66, 91), 3: F(55, 91), 4: F(1979, 4004), 5: F(2243, 5880)}
ladder = {}
for j in range(1, 14):
    best, arg, second = None, None, None
    for P in itertools.combinations(range(1, 14), j):
        L = lam(P)
        if best is None or L < best:
            second = best
            best, arg = L, P
        elif second is None or L < second:
            if L != best:
                second = L
    ladder[j] = (best, arg, second)
    known = KNOWN.get(j)
    chk = "" if known is None else f"  THM-576: {known}  MATCH={best == known}"
    marg = "" if second is None else f"  margin-to-2nd={float(second - best):.5f}"
    print(f"  j={j:2d} (k={13-j:2d})  cap = {str(best):>16} = {float(best):.6f}  argmin={arg}{chk}{marg}")

# ---------------- B. decomposition of cap_8, cap_9 ----------------
print("\n" + "=" * 96)
print("B. CLEAN-BRANCH DECOMPOSITION of the j=4, j=5 caps (binding rows k=9, k=8)")
print("=" * 96)
def dfold(P):
    """Exact |intersection of dangers over P|."""
    cur = DANGERS[P[0]]
    for v in P[1:]:
        cur = merge(inter(cur, DANGERS[v]))
    return measure(cur)

for j in (4, 5):
    cap, P, _ = ladder[j]
    print(f"\n  j={j}: minimizer P = {P}, cap = {cap}")
    total = F(1)
    for d in range(1, j + 1):
        sgn = -1 if d % 2 == 1 else 1
        terms = []
        for Q in itertools.combinations(P, d):
            m = dfold(Q)
            clean = (m == 2 * R / max(Q))
            terms.append((Q, m, clean))
            total += sgn * m
        ncl = sum(1 for _, _, c in terms if c)
        print(f"    d={d}: {len(terms)} terms, clean(=2r/vmax): {ncl}/{len(terms)}")
        if d >= 2:
            for Q, m, c in terms:
                tag = "CLEAN" if c else f"WRAPPED (excess {m - 2*R/max(Q)})"
                print(f"       |^{Q}| = {str(m):>10}  {tag}")
    print(f"    inclusion-exclusion total = {total}  == cap? {total == cap}")

# ---------------- C. d=3 criterion hunt ----------------
print("\n" + "=" * 96)
print("C. THE d=3 CLEAN-BRANCH CRITERION (all u<v<w <= 20, r = 1/14)")
print("=" * 96)
rows = []
for u, v, w in itertools.combinations(range(1, 21), 3):
    m = dfold((u, v, w))
    clean = (m == 2 * R / w)
    rows.append((u, v, w, m, clean))
nclean = sum(1 for r_ in rows if r_[4])
print(f"  triples: {len(rows)}, clean: {nclean}, wrapped: {len(rows) - nclean}")

# criterion candidates
def cand_pairwise(u, v, w): return v + w <= 14
def cand_sum_le(u, v, w): return v + w <= 14 + u
def cand_uvw(u, v, w): return u + v + w <= 21
for name, cand in [("pairwise v+w<=14 (sufficient?)", cand_pairwise),
                   ("v+w <= 14+u", cand_sum_le), ("u+v+w <= 21", cand_uvw)]:
    fp = sum(1 for (u, v, w, m, c) in rows if cand(u, v, w) and not c)
    fn = sum(1 for (u, v, w, m, c) in rows if not cand(u, v, w) and c)
    print(f"  candidate [{name}]: predicts-clean-but-wrapped {fp}, predicts-wrapped-but-clean {fn}")

print("\n  boundary sample -- wrapped triples with smallest v+w-u:")
wrapped = sorted([(v + w - u, (u, v, w), m) for (u, v, w, m, c) in rows if not c])[:8]
for key, t, m in wrapped:
    print(f"    {t}: v+w-u = {key}, |^| = {m} (2r/w = {2*R/t[2]})")
print("  clean triples with LARGEST v+w-u:")
cleanb = sorted([(v + w - u, (u, v, w), m) for (u, v, w, m, c) in rows if c], reverse=True)[:8]
for key, t, m in cleanb:
    print(f"    {t}: v+w-u = {key}, |^| = {m}")

print("\nDONE.")

# ---------------- C'. gcd-reduction law test ----------------
print("\n" + "=" * 96)
print("C'. GCD-REDUCTION LAW: |^ D_Q| = 2r/max(Q/g), g = gcd(Q)?  (all 1140 triples <= 20)")
print("=" * 96)
from math import gcd as _g
ok, bad = 0, []
for u, v, w in itertools.combinations(range(1, 21), 3):
    m = dfold((u, v, w))
    g = _g(_g(u, v), w)
    pred = 2 * R / (w // g)
    if m == pred: ok += 1
    else: bad.append(((u, v, w), m, pred, g))
print(f"  gcd-law holds: {ok}/1140;  violations: {len(bad)}")
for t, m, p, g in bad[:10]:
    print(f"    {t} (g={g}): |^| = {m} vs predicted {p}  (excess {m - p})")
if bad:
    prim_bad = [b for b in bad if b[3] == 1]
    print(f"  violations with gcd=1 (primitive): {len(prim_bad)}")
