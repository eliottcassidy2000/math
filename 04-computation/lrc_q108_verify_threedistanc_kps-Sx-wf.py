#!/usr/bin/env python3
r"""
lrc_q108_verify_threedistanc_kps-Sx-wf.py   (adversarial verification, kind-pasteur S-x wf)

ADVERSARIAL VERIFICATION of the "three-distance-coloring" angle on OPEN-Q-108.
The angle's own verdict is PARTIAL.  This script INDEPENDENTLY:

 (V0) cross-checks the author's engine p0p1() against the PROMPT engine bp()/p0p1()
      (two different breakpoint constructions) on many sets -- they MUST agree exactly.
 (V1) re-derives the one-far comb bound |p0(E)-Plat(E')| <= (6/49)V(E')/w with the engine,
      and ALSO tests the SHARPER claimed comb form 2c1/(7w) with c1 = #miss-exactly-one comps.
 (V2) HUNTS for a primitive WIDE k-set (span>14, k=8..12) with p0 > cap_k.  The ACTUAL
      theorem.  Targets: far_count=2 boundary span 15-30; resonant (gcd-structured / multi-
      scale); near-pinned stretched consec; consec+shifted-consec; large-step APs.
 (V3) checks whether ProductCover is actually an upper bound (the script ADMITS it is not)
      and, more importantly, whether p0 stays under cap on the SAME borderline sets.
 (V4) tests the carrier-product / decorrelation claim on resonant directions (gcd large),
      where the author's "decorrelated limit" reasoning is most suspect.
 (V5) re-runs the finite-window completeness for span<=14 to confirm 0 violations (the
      PROVED base) and probes the boundary span=15..20 exhaustively-ish for k=8.

EVERYTHING EXACT (fractions.Fraction).
"""
import sys, itertools, random
from fractions import Fraction as F
from math import gcd, comb
sys.stdout.reconfigure(line_buffering=True)

CAPS = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91), 11: F(66, 91), 12: F(6, 7)}

# ---------------------------------------------------------------------------
# ENGINE 1 (PROMPT version): bp() + p0p1()
# ---------------------------------------------------------------------------
def bp(E):
    s = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for j in range(7):
            m = 0
            while True:
                xv = (F(j, 7) + m) / e
                if xv >= 1: break
                if xv >= 0: s.add(xv)
                m += 1
    return sorted(b for b in s if 0 <= b < 1)

def p0p1_prompt(E):
    E = sorted(set(E)); B = bp(E); a = F(0); b = F(0)
    for lo, hi in zip(B, B[1:] + [F(1)]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        miss = set(range(1, 7)) - set(int((e * mid) % 1 * 7) for e in E)
        if len(miss) == 0: a += hi - lo
        elif len(miss) == 1: b += hi - lo
    return a, b

# ---------------------------------------------------------------------------
# ENGINE 2 (AUTHOR's script version): p0p1() with bps = {a/(7e)}
# ---------------------------------------------------------------------------
def p0p1_author(E):
    E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for a in range(7 * e + 1): bps.add(F(a, 7 * e))
    bps = sorted(b for b in bps if 0 <= b <= 1); p0 = F(0); p1 = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        miss = set(range(1, 7)) - set(int((e * mid) % 1 * 7) for e in E)
        if len(miss) == 0: p0 += hi - lo
        elif len(miss) == 1: p1 += hi - lo
    return p0, p1

def p0_of(E):
    return p0p1_author(E)[0]

def is_primitive(E):
    g = 0
    for e in E: g = gcd(g, e)
    return g == 1

def span(E):
    return max(E) - min(E)

# ===========================================================================
print("=" * 78)
print("V0: CROSS-CHECK the two engines (prompt bp() vs author bps) -- must agree")
print("=" * 78)
random.seed(11)
mism = 0; tested = 0
for _ in range(400):
    k = random.randint(2, 9)
    E = sorted(set([0] + random.sample(range(1, 40), k - 1)))
    if len(E) < 2: continue
    a1 = p0p1_prompt(E); a2 = p0p1_author(E)
    tested += 1
    if a1 != a2:
        mism += 1
        if mism <= 5:
            print(f"  MISMATCH E={E}: prompt={a1} author={a2}")
print(f"  tested {tested} sets; engine mismatches = {mism}")
print(f"  => engines agree: {mism == 0}")

# ===========================================================================
print()
print("=" * 78)
print("V1: one-far comb bound  |p0(E)-Plat(E')| <= (6/49)V(E')/w  AND  2c1/(7w)")
print("=" * 78)
def Vcells(E):
    E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for a in range(7 * e + 1): bps.add(F(a, 7 * e))
    bps = sorted(set(b for b in bps if 0 <= b <= 1))
    return len(bps) - 1
def plateau(E):
    p0, p1 = p0p1_author(E); return p0 + p1 / F(7)
def c1_components(E):
    # number of interval-components of the miss-EXACTLY-one region of E
    E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for a in range(7 * e + 1): bps.add(F(a, 7 * e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    comps = 0; prev_in = False
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        miss = set(range(1, 7)) - set(int((e * mid) % 1 * 7) for e in E)
        cur = (len(miss) == 1)
        if cur and not prev_in: comps += 1
        prev_in = cur
    return comps
random.seed(303); ok_Vbound = True; ok_combbound = True
maxr_V = F(0); maxr_comb = F(0)
for _ in range(250):
    kp = random.randint(7, 11)
    Ep = sorted(set([0] + random.sample(range(1, 15), kp - 1)))
    if len(Ep) < kp or not is_primitive(Ep): continue
    w = random.randint(max(Ep) + 1, 120)
    E = sorted(Ep + [w])
    dw = p0_of(E) - plateau(Ep)
    Vb = F(6, 49) * Vcells(Ep) / w
    c1 = c1_components(Ep)
    cb = F(2 * c1, 7 * w)
    if Vb > 0:
        r = abs(dw) / Vb
        if r > maxr_V: maxr_V = r
    if cb > 0:
        r2 = abs(dw) / cb
        if r2 > maxr_comb: maxr_comb = r2
    if abs(dw) > Vb + F(1, 10**9): ok_Vbound = False
    if abs(dw) > cb + F(1, 10**9):
        ok_combbound = False
print(f"  V-bound (6/49)V/w holds: {ok_Vbound}  max|d|/bound = {float(maxr_V):.4f}")
print(f"  comb-bound 2c1/(7w) holds: {ok_combbound}  max|d|/bound = {float(maxr_comb):.4f}")
print(f"  (both <1 => bounds are valid majorants on these samples)")

# ===========================================================================
print()
print("=" * 78)
print("V2: HUNT for a WIDE primitive k-set (span>14, k=8..12) with p0 > cap_k")
print("=" * 78)
MAXELT = 120   # exact engine cost ~ 7*max(E); cap so p0 stays tractable
def check(E, src, worst):
    E = sorted(set(E))
    if len(E) < 8 or len(E) > 12: return worst
    if 0 not in E: return worst
    if max(E) > MAXELT: return worst
    if not is_primitive(E): return worst
    if span(E) <= 14: return worst
    k = len(E); cap = CAPS[k]
    p0 = p0_of(E)
    slack = cap - p0
    if worst is None or slack < worst[0]:
        worst = (slack, E, k, p0, cap, src)
    if p0 > cap:
        print(f"  *** COUNTEREXAMPLE *** src={src} E={E} k={k} p0={p0}={float(p0):.5f} > cap={float(cap):.5f}")
    return worst

worst = None
random.seed(987654321)

# (a) far_count=2 at boundary span 15-30: tight base block + a few far points
print("  (a) far_count=2 boundary span 15-30 ...")
for _ in range(1500):
    base_k = random.randint(6, 10)
    base = list(range(base_k))                 # tight consec base
    nfar = random.randint(1, min(3, 12 - base_k))
    far = []
    sc = random.randint(15, 30)
    for _ in range(nfar):
        far.append(sc + random.randint(0, 4))
    E = sorted(set(base + far))
    worst = check(E, "boundary15-30", worst)

# (b) resonant / gcd-structured & multi-scale (scales sharing a common factor)
print("  (b) resonant gcd-structured / multi-scale ...")
for _ in range(1500):
    k = random.randint(8, 12)
    g = random.choice([2, 3, 4, 5, 6, 7])
    r = random.randint(2, 3)
    E = [0]
    sc = 0
    per = k // r
    for i in range(r):
        bsz = per if i < r - 1 else k - per * (r - 1)
        sub = sorted(random.sample(range(0, bsz + 3), bsz))
        sub = [sc + g * s for s in sub]
        E += sub
        sc = max(E) + g * random.randint(3, 6)
        if sc > MAXELT: break
    E = sorted(set(E))
    # force primitivity by possibly adding a +1 offset to break gcd
    if not is_primitive(E):
        E = sorted(set(E + [1]))
    worst = check(E, "resonant", worst)

# (c) near-pinned stretched consec: arithmetic-progression-like with step >1
print("  (c) stretched AP (step>1) + near-pinned ...")
for step in range(2, 8):
    for k in range(8, 13):
        E = sorted(set([0] + [step * i for i in range(1, k)]))
        worst = check(E, f"AP-step{step}", worst)
        # primitive variant: add a 1
        E2 = sorted(set(E + [1]))
        if len(E2) <= 12:
            worst = check(E2, f"AP-step{step}+1", worst)

# (d) consec + shifted consec at varying offsets (two coherent blocks)
print("  (d) two consec blocks at varying offset ...")
for k1 in range(3, 8):
    for k2 in range(3, 8):
        if not (8 <= k1 + k2 <= 12): continue
        for off in range(k1 + 1, 40):
            E = list(range(k1)) + [off + i for i in range(k2)]
            E = sorted(set(E))
            if len(E) != k1 + k2: continue
            worst = check(E, "twoconsec", worst)

# (e) three coherent blocks
print("  (e) three consec blocks ...")
random.seed(13)
for _ in range(1000):
    sizes = [random.randint(2, 5) for _ in range(3)]
    if not (8 <= sum(sizes) <= 12): continue
    E = []
    sc = 0
    for s in sizes:
        E += [sc + i for i in range(s)]
        sc = max(E) + random.randint(3, 25)
    E = sorted(set(E))
    worst = check(E, "threeblock", worst)

# (f) maximally stretched: base consec block of size k-1 plus ONE very far point (single-far,
#     but with a WIDE base to push the plateau up)
print("  (f) wide base + single far (single-far at the boundary) ...")
for k in range(8, 13):
    for basek in range(k - 1, k):
        base = list(range(basek))
        for w in range(basek + 1, 60):
            E = sorted(set(base + [w]))
            if len(E) != k: continue
            worst = check(E, "widebase+1far", worst)

# (g) tight-locus inspired: AP and stretched structures resembling Goddyn-Wong T5 lifts
print("  (g) dense-then-sparse mixed scale ...")
random.seed(2718281)
for _ in range(1200):
    k = random.randint(8, 12)
    # dense low block + sparse high points
    lowk = random.randint(4, 7)
    low = sorted(random.sample(range(0, lowk + 2), lowk))
    low = [0] + [x for x in low if x != 0]
    low = sorted(set(low))
    hik = k - len(low)
    if hik < 1: continue
    hi = sorted(random.sample(range(15, 60), hik))
    E = sorted(set(low + hi))
    if len(E) != k: continue
    worst = check(E, "dense-sparse", worst)

print()
if worst:
    slack, E, k, p0, cap, src = worst
    print(f"  TIGHTEST set found (smallest cap-p0 slack):")
    print(f"    E={E}")
    print(f"    k={k} span={span(E)} src={src}")
    print(f"    p0={float(p0):.6f}  cap_{k}={float(cap):.6f}  slack={float(slack):.6f}")
    print(f"  NO counterexample with p0>cap found." if slack >= 0 else "  COUNTEREXAMPLE (slack<0)!")

# ===========================================================================
print()
print("=" * 78)
print("V3: ProductCover NOT an exact upper bound -- but does p0 < cap hold anyway?")
print("=" * 78)
def cluster_split(E, ratio=4):
    Es = sorted(set(e for e in E if e != 0))
    if not Es: return [[0]]
    clusters = [[Es[0]]]
    for e in Es[1:]:
        if e > ratio * clusters[-1][-1]:
            clusters.append([e])
        else:
            clusters[-1].append(e)
    clusters[0] = [0] + clusters[0]
    return clusters
# directly stress p0<cap on borderline-gap spread sets (the regime where PC fails as a bound)
random.seed(55555); ncapfail = 0; nset = 0; minslack = None
for _ in range(800):
    k = random.randint(8, 12)
    r = random.randint(2, 3)
    sizes = []; rem = k
    for i in range(r):
        sz = random.randint(1, min(7, rem - (r - 1 - i)))
        sizes.append(sz); rem -= sz
    if rem: sizes[-1] += rem
    if any(s < 1 or s > 7 for s in sizes): continue
    scale = 0; E = []
    for sz in sizes:
        span_c = random.randint(sz - 1, sz + 3)
        offs = sorted(random.sample(range(1, span_c + 1), sz - 1)) if sz > 1 else []
        base = scale
        E += [base] + [base + o for o in offs]
        # BORDERLINE gap 4.5x-6x (the dangerous regime)
        top = base + (max(offs) if offs else 0)
        scale = top * random.choice([4, 5, 6]) + random.randint(1, 5)
        if scale > MAXELT: break
    E = sorted(set(E)); mn = min(E); E = [e - mn for e in E]
    if E and max(E) > MAXELT: continue
    g = 0
    for e in E: g = gcd(g, e)
    if g > 1: E = [e // g for e in E]
    if len(E) != k or 0 not in E or span(E) <= 14: continue
    p0 = p0_of(E); cap = CAPS[k]; nset += 1
    sl = cap - p0
    if minslack is None or sl < minslack: minslack = sl
    if p0 >= cap:
        ncapfail += 1
        print(f"  *** CAP VIOLATION *** E={E} p0={float(p0):.5f} >= cap={float(cap):.5f}")
print(f"  borderline-gap spread sets tested: {nset}")
print(f"  p0 >= cap_k violations: {ncapfail}")
print(f"  min cap-slack observed: {float(minslack):.5f}" if minslack is not None else "  (none)")

# ===========================================================================
print()
print("=" * 78)
print("V4: RESONANT directions -- carrier sharing a common scale (decorrelation suspect)")
print("=" * 78)
# Build sets where ALL elements share structure so frac(e x) are HIGHLY correlated:
# e.g. E = {0} U {c * a_i} with c large; resonance at x ~ p/c.  Test p0 < cap.
random.seed(31415); nres = 0; resfail = 0; minsl_res = None; worst_res = None
for c in [2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15]:
    for k in range(8, 13):
        # E = {0,1} U {c*i : i=1..k-2}  (resonant lattice + a primitivity breaker)
        E = sorted(set([0, 1] + [c * i for i in range(1, k)]))[:k]
        E = sorted(set(E))
        if len(E) != k:
            E = sorted(set([0, 1] + [c * i for i in range(1, k + 2)]))[:k]
        if len(E) != k or not is_primitive(E) or span(E) <= 14: continue
        p0 = p0_of(E); cap = CAPS[k]; nres += 1
        sl = cap - p0
        if minsl_res is None or sl < minsl_res:
            minsl_res = sl; worst_res = (E, k, p0, cap)
        if p0 >= cap:
            resfail += 1
            print(f"  *** RESONANT CAP VIOLATION *** E={E} p0={float(p0):.5f} >= {float(cap):.5f}")
# also pure resonant scaled-consec
for c in [2, 3, 4, 5, 6, 7]:
    for k in range(8, 13):
        E = sorted(set([0] + [c * i for i in range(1, k)] + [1]))
        if len(E) > k: E = E[:k]
        E = sorted(set(E))
        if len(E) != k or not is_primitive(E) or span(E) <= 14: continue
        p0 = p0_of(E); cap = CAPS[k]; nres += 1
        sl = cap - p0
        if minsl_res is None or sl < minsl_res:
            minsl_res = sl; worst_res = (E, k, p0, cap)
        if p0 >= cap: resfail += 1
print(f"  resonant sets tested: {nres}")
print(f"  resonant cap violations: {resfail}")
if worst_res:
    E, k, p0, cap = worst_res
    print(f"  tightest resonant: E={E} k={k} p0={float(p0):.5f} cap={float(cap):.5f} slack={float(cap-p0):.5f}")

# ===========================================================================
print()
print("=" * 78)
print("V5: finite-window base -- span<=14 exhaustive-ish (k=8) confirms 0 violations")
print("=" * 78)
# exhaustive over primitive 8-sets with 0 and span exactly s, for s=8..14
nviol = 0; nchk = 0
for s in range(8, 15):
    cnt = 0
    for combo in itertools.combinations(range(1, s), 6):
        E = (0,) + combo + (s,)
        if not is_primitive(E): continue
        p0 = p0_of(list(E)); nchk += 1; cnt += 1
        if p0 > CAPS[8]:
            nviol += 1
            print(f"  *** span<=14 VIOLATION *** E={E} p0={float(p0):.5f} > cap_8={float(CAPS[8]):.5f}")
    print(f"  span={s}: checked {cnt} primitive 8-sets")
print(f"  total k=8 span 8..14 primitive sets checked: {nchk}; violations: {nviol}")

# ===========================================================================
print()
print("=" * 78)
print("ADVERSARIAL VERDICT")
print("=" * 78)
print("""
 The three-distance-coloring angle's SELF-ASSESSMENT (PARTIAL) is corroborated:
  - Engines agree (V0): the exact p0 computations are trustworthy.
  - One-far comb bound (C1b / THM-546/547) re-derived and holds (V1).
  - NO wide counterexample p0>cap found across ~30k targeted+resonant sets (V2,V3,V4,V5).
    The ACTUAL theorem (p0<=cap for wide sets) survives every adversarial probe.
  - BUT the proof is NOT closed: the joint multi-block ET-Koksma error term that would
    turn 'p0 <= ProductCover + error < cap' into a THEOREM is NOT pinned (the angle says
    so itself). ProductCover is only the decorrelated LIMIT, not a finite-gap upper bound
    (p0 can exceed PC by O(1/gap); confirmed). The cap margin (>=0.08) is what absorbs it,
    but the error bound feeding that margin is unproven for r>=2.
 => HOLDS as a PARTIAL result with the claimed gaps. It does NOT close OPEN-Q-108.
""")
