#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADVERSARIAL VERIFICATION of the LRC(14) "finite-check-to-threshold" claim
(file: lrc14_ck_finite-check-to-threshold_kps-Sx-wf.py), kind-pasteur Sx-wf.

The claimed RESULT is PARTIAL: the finite check is sound CONDITIONAL on the
imported bound  C(k) := sup_{E',w} w|Delta_w(E')| <= 2.71*(k-2).

This verifier does four independent things:

  (1) RE-DERIVE the exact constants (cap_k, Q(k-1), margin_k, threshold T_k) and
      the dovetail identity p0(E) = Phi(E') + Delta_w, EXACTLY.

  (2) HUNT for (E', w) with  w|Delta_w| > 2.71*(k-2)  -- the load-bearing bound.
      Aggressive search: multi-scale RESONANT w = lcm of cluster scales, w =
      0 mod many e's; large #scales; large sigma; geometric / arithmetic mixed cores.
      One violation makes the imported C(k) bound FALSE (with witness).

  (3) STRESS the per-cluster "slope <= 2.71" certification: test_D in the sibling
      uses only 4 hand-picked geometric designs.  We scan a far broader family of
      multi-cluster cores and measure worst w|Delta_w|/(r-1) to see if 2.71 holds.

  (4) Verify the FINITE-CHECK machinery itself:
      - fast p0 == reference p0 (exact);
      - upward-closure monotonicity (the prune's rigorous basis);
      - the upward-closure completion bound is a VALID upper bound;
      - re-run the C=3 EXHAUSTIVE finite check independently and confirm
        0 violations + consec argmax.

EXACT Fraction arithmetic throughout.  Floats only for printing / search ranking.
"""
import sys, itertools, random
from fractions import Fraction as F
from functools import reduce
from math import gcd, comb

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


# ======================================================================
# EXACT machinery (independent re-implementation)
# ======================================================================
def _missed_dist(E):
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for a in range(0, 7 * e + 1):
            bps.add(F(a, 7 * e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [F(0)] * 7
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i + 1]
        if hi == lo:
            continue
        mid = (lo + hi) / 2
        hit = set()
        for e in E:
            v = e * mid
            v = v - (v.numerator // v.denominator)
            hit.add((v.numerator * 7) // v.denominator)
        t = sum(1 for j in range(1, 7) if j not in hit)
        p[t] += hi - lo
    return p

def p0_ref(E):
    return _missed_dist(E)[0]

def Phi(E):
    p = _missed_dist(E)
    return p[0] + F(1, 7) * p[1]

def _lcm(a, b):
    return a // gcd(a, b) * b

_ALL6 = 0b1111110

def p0_fast(E):
    Enz = sorted(set(x for x in E if x != 0))
    if not Enz:
        return F(0)
    L = reduce(_lcm, Enz)
    D = 7 * L
    den2 = 2 * L
    bps = {0, D}
    for e in Enz:
        step = D // (7 * e)
        x = 0
        for _ in range(7 * e + 1):
            bps.add(x)
            x += step
    bps = sorted(bps)
    num = 0
    for i in range(len(bps) - 1):
        lo = bps[i]; hi = bps[i + 1]
        if hi <= lo:
            continue
        midnum = lo + hi
        mask = 0
        for e in Enz:
            mask |= 1 << ((e * midnum // den2) % 7)
        if mask & _ALL6 == _ALL6:
            num += hi - lo
    return F(num, D)

def primitive(E):
    nz = [e for e in E if e != 0]
    return len(nz) > 0 and reduce(gcd, nz) == 1


# ----------------------------------------------------------------------
# w*Delta_w engine (the EXACT tool from the prompt, verbatim semantics)
# ----------------------------------------------------------------------
def G0(y):
    y = y - int(y)
    if y < 0:
        y = y + 1
    return y * F(6, 7) if y < F(1, 7) else F(6, 49) - (y - F(1, 7)) * F(1, 7)

def _cells_1miss(Ep):
    Ep = sorted(set(Ep)); bps = {F(0), F(1)}
    for e in Ep:
        if e == 0:
            continue
        for a in range(0, 7 * e + 1):
            bps.add(F(a, 7 * e))
    bps = sorted(b for b in bps if 0 <= b <= 1); out = []
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i + 1]
        if hi == lo:
            continue
        mid = (lo + hi) / 2; hit = set()
        for e in Ep:
            v = e * mid; v = v - (v.numerator // v.denominator)
            hit.add((v.numerator * 7) // v.denominator)
        m = [j for j in range(1, 7) if j not in hit]
        if len(m) == 1:
            out.append((lo, hi, m[0]))
    return out

def wDelta_signed(Ep, w):
    cells = _cells_1miss(Ep)
    return sum((G0(w * b - F(s, 7)) - G0(w * a - F(s, 7)) for (a, b, s) in cells), F(0))

def wDelta(Ep, w):
    return abs(wDelta_signed(Ep, w))   # EXACT Fraction


# ----------------------------------------------------------------------
# scale clusters (same definition as the sibling: ratio>3 splits)
# ----------------------------------------------------------------------
def scale_clusters(Ep, ratio=3):
    es = sorted(e for e in Ep if e > 0)
    if not es:
        return []
    clusters = [[es[0]]]
    for a, b in zip(es, es[1:]):
        if b > a * ratio:
            clusters.append([b])
        else:
            clusters[-1].append(b)
    return clusters

def num_scales(Ep, ratio=3):
    return len(scale_clusters(Ep, ratio))


# ======================================================================
# (1) EXACT CONSTANTS + DOVETAIL IDENTITY
# ======================================================================
CAP = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91)}
C_CLUSTER = F(271, 100)

def part1_constants():
    print("=" * 78)
    print("(1) EXACT CONSTANTS  (re-derived independently)")
    print("=" * 78)
    info = {}
    for k in (8, 9, 10):
        cap = CAP[k]
        Q = Phi(list(range(k - 1)))
        margin = cap - Q
        Cload = C_CLUSTER * (k - 2)
        T = Cload / margin
        Tc = T.numerator // T.denominator + (1 if T.numerator % T.denominator else 0)
        p0c = p0_ref(list(range(k)))
        info[k] = (cap, Q, margin, Tc)
        print(f"  k={k}: cap={cap}={float(cap):.5f}  Q(k-1)={Q}={float(Q):.5f}")
        print(f"        margin={margin}={float(margin):.5f}  C(k)<=2.71*(k-2)={float(Cload):.2f}")
        print(f"        T_k=C/margin={float(T):.2f}  ceil={Tc}   p0(consec_k)={float(p0c):.5f}")
        print(f"        p0(consec_k) < cap? {p0c < cap}   (slack {float(cap-p0c):.5f})")
    # expected
    exp = {8: (F(2243,5880), F(289,1470), F(1087,5880), 88),
           9: (F(1979,4004), F(621,1715), F(129643,980980), 144),
           10:(F(55,91),     F(1229,2744),F(5583,35672), 139)}
    allok = all(info[k] == exp[k] for k in (8,9,10))
    print(f"\n  constants match claimed values EXACTLY: {allok}")
    return info, allok

def part1b_dovetail(ntests=200, seed=7):
    print("\n" + "=" * 78)
    print("(1b) DOVETAIL identity  p0(E) = Phi(E') + Delta_w   (Delta_w = wDelta/w, signed)")
    print("=" * 78)
    rng = random.Random(seed)
    nbad = 0; nfastbad = 0
    for _ in range(ntests):
        k1 = rng.choice([6, 7, 8, 9])
        Ep = [0] + sorted(rng.sample(range(1, 22), k1 - 1))
        w = rng.randint(max(Ep) + 1, max(Ep) + 30)
        if w in Ep:
            continue
        E = sorted(Ep + [w])
        lhs = p0_ref(E) - Phi(Ep)          # = Delta_w  exactly (signed)
        rhs = wDelta_signed(Ep, w) / w      # = Delta_w  from engine (signed)
        if lhs != rhs:
            nbad += 1
            if nbad <= 3:
                print(f"    MISMATCH Ep={Ep} w={w}: lhs={lhs} rhs={rhs}")
        if p0_fast(E) != p0_ref(E):
            nfastbad += 1
    print(f"  dovetail identity exact: {ntests-nbad}/{ntests}  ({'OK' if nbad==0 else 'FAIL'})")
    print(f"  fast p0 == ref p0      : {ntests-nfastbad}/{ntests}  ({'OK' if nfastbad==0 else 'FAIL'})")
    return nbad == 0 and nfastbad == 0


# ======================================================================
# (2) HUNT for w|Delta_w| > 2.71*(k-2)
# ======================================================================
def resonant_ws(Ep, extra=60):
    """Candidate resonant w: lcm of cluster scales, lcm of all, multiples of small e's,
       and a dense band.  Returns sorted set of w to test."""
    es = sorted(e for e in Ep if e > 0)
    cands = set()
    # lcm of all
    L = reduce(_lcm, es)
    # multiples and divisors-of-products
    cl = scale_clusters(Ep)
    cluster_lcms = [reduce(_lcm, c) for c in cl]
    base = set([L] + cluster_lcms)
    # products of subsets of cluster lcms (resonance with several clusters)
    for r in range(1, len(cluster_lcms) + 1):
        for sub in itertools.combinations(cluster_lcms, r):
            base.add(reduce(_lcm, sub))
    for b in list(base):
        for mult in range(1, 6):
            cands.add(b * mult)
            for d in (-2, -1, 1, 2):
                if b * mult + d > 1:
                    cands.add(b * mult + d)
    # multiples of each e (resonant with a single runner)
    for e in es:
        for mult in range(1, 30):
            cands.add(e * mult)
            cands.add(e * mult + 1)
            cands.add(e * mult - 1)
    # dense low band
    for w in range(2, max(Ep) * 4 + extra):
        cands.add(w)
    return sorted(w for w in cands if w >= 2)

def hunt_violation(label, cores_iter, threshold_fn, wbudget_note=""):
    """cores_iter yields (E', k_of_full_set).  threshold_fn(k) = 2.71*(k-2).
       For each core, scan resonant w; report worst w|Delta| and any violation."""
    worst = F(0); worst_at = None; nviol = 0; viol_list = []
    ncores = 0
    for Ep, k in cores_iter:
        ncores += 1
        thr = threshold_fn(k)
        for w in resonant_ws(Ep):
            v = wDelta(Ep, w)
            if v > worst:
                worst = v; worst_at = (tuple(Ep), w, k)
            if v > thr:
                nviol += 1
                if len(viol_list) < 12:
                    viol_list.append((v, tuple(Ep), w, k, thr))
    print(f"  [{label}] scanned {ncores} cores {wbudget_note}")
    if worst_at:
        Ep, w, k = worst_at
        print(f"    worst w|Delta_w| = {float(worst):.4f} at E'={Ep} w={w}  "
              f"(k={k}, thr=2.71*(k-2)={float(threshold_fn(k)):.2f}, "
              f"#clusters={num_scales(Ep)})")
    print(f"    violations of 2.71*(k-2): {nviol}")
    for v, Ep, w, k, thr in sorted(viol_list, reverse=True)[:6]:
        print(f"      *** w|Delta|={float(v):.4f} > {float(thr):.2f}  E'={Ep} w={w} k={k}")
    return worst, worst_at, nviol, viol_list

def part2_hunt():
    print("\n" + "=" * 78)
    print("(2) HUNT: is  C(k) <= 2.71*(k-2)  ?   (look for violations)")
    print("=" * 78)
    thr = lambda k: float(C_CLUSTER) * (k - 2)
    allviol = 0

    # --- 2a: multi-scale RESONANT geometric cores (the worst case per the prompt) ---
    print("\n  2a) Multi-scale resonant cores: r clusters of size sz at scales S^i.")
    def geo_cores():
        for r in range(2, 6):
            for S in (5, 7, 8, 10, 12, 16, 20):
                for sz in (2, 3):
                    core = set()
                    for i in range(r):
                        b = 0 if i == 0 else S ** i
                        core |= set(range(b, b + sz))
                    core = sorted(core)
                    # k = full set size = |E'| + 1
                    kfull = len(core) + 1
                    if kfull > 11 or not primitive(core):
                        continue
                    yield core, kfull
    w, wa, nv, _ = hunt_violation("geo", geo_cores(), thr,
                                  "(resonant w incl. lcm of clusters)")
    allviol += nv

    # --- 2b: arithmetic-progression clusters at multiple scales (collapse-friendly) ---
    print("\n  2b) AP-block clusters at scales chosen so lcm collapses many breakpoints.")
    def ap_cores():
        for r in range(2, 5):
            for base in (6, 12, 24, 30, 36, 60):
                for sz in (2, 3):
                    core = set([0])
                    for i in range(1, r):
                        c = base * i
                        core |= set(range(c, c + sz))
                    core |= set(range(1, sz))
                    core = sorted(core)
                    kfull = len(core) + 1
                    if kfull > 11 or not primitive(core) or len(core) < 4:
                        continue
                    yield core, kfull
    w, wa, nv, _ = hunt_violation("ap", ap_cores(), thr)
    allviol += nv

    # --- 2c: HIGHLY composite resonant w; cores made of divisors of a smooth number ---
    print("\n  2c) Smooth-number cores (divisors of 2^a*3^b*5^c) -- maximal resonance.")
    def smooth_cores():
        smooths = [60, 120, 180, 210, 420, 840, 2520]
        for N in smooths:
            divs = sorted(d for d in range(1, N + 1) if N % d == 0)
            # pick small + scattered divisors to span scales
            for k in (8, 9, 10):
                sub = [d for d in divs if d <= N]
                # take a spread: smallest few + a few large
                if len(sub) < k - 1:
                    continue
                small = sub[:max(2, (k - 1) // 2)]
                large = sub[-(k - 1 - len(small)):]
                core = sorted(set([0]) | set(small) | set(large))
                if len(core) != k or not primitive(core):
                    # adjust
                    core = sorted(set([0]) | set(sub[:k - 1]))
                    if len(core) != k or not primitive(core):
                        continue
                yield core, k
    w, wa, nv, _ = hunt_violation("smooth", smooth_cores(), thr)
    allviol += nv

    # --- 2d: random multi-scale cores with random resonant w ---
    print("\n  2d) Random multi-scale cores (|E'|=7..9), resonant + random w.")
    rng = random.Random(2024)
    def rand_cores():
        for _ in range(4000):
            r = rng.randint(2, 4)
            scales = sorted(rng.sample([3,4,5,6,7,8,10,12,15,20,30,40,60,80,120], r))
            core = set([0, 1])
            for sc in scales:
                w0 = rng.randint(1, 3)
                core |= set(range(sc, sc + rng.choice([1,2,3])))
            core = sorted(core)
            kfull = len(core) + 1
            if kfull < 8 or kfull > 11 or not primitive(core):
                continue
            yield core, min(kfull, 10)
    w, wa, nv, _ = hunt_violation("rand", rand_cores(), thr)
    allviol += nv

    # --- 2e: LARGE sigma stress (big runners; check the sigma-bound is what matters) ---
    print("\n  2e) Large-sigma cores (big runners) -- sigma-bound says |Delta|<=(6/7)sigma/w.")
    def bigsigma_cores():
        for top in (50, 80, 120, 200):
            for r in range(2, 4):
                core = set([0, 1, 2])
                step = top // r
                for i in range(1, r + 1):
                    c = step * i
                    core |= set(range(c, c + 2))
                core = sorted(core)
                kfull = len(core) + 1
                if kfull > 11 or not primitive(core) or len(core) < 7:
                    continue
                yield core, min(kfull, 10)
    w, wa, nv, _ = hunt_violation("bigsigma", bigsigma_cores(), thr)
    allviol += nv

    print(f"\n  => TOTAL violations of C(k)<=2.71*(k-2) across all hunts: {allviol}")
    return allviol


# ======================================================================
# (3) STRESS the per-cluster slope <= 2.71
# ======================================================================
def part3_slope_stress():
    print("\n" + "=" * 78)
    print("(3) STRESS the slope:  worst w|Delta_w| / (r-1)  over MANY multi-cluster cores")
    print("=" * 78)
    print("  sibling test_D used only 4 hand-picked designs.  We scan broadly and report")
    print("  the worst ratio per #clusters r.  If this exceeds 2.71, the slope is wrong.")
    best_by_r = {}   # r -> (ratio, core, w, raw)
    rng = random.Random(99)
    designs = []
    # geometric, varied
    for r in range(2, 6):
        for S in (4,5,6,7,8,9,10,12,15,20,25,30):
            for sz in (2, 3):
                core = set()
                for i in range(r):
                    b = 0 if i == 0 else S ** i
                    core |= set(range(b, b + sz))
                core = sorted(core)
                if len(core) <= 11 and primitive(core):
                    designs.append(core)
    # linear-scale (not geometric): clusters at c, 2c, 3c...
    for r in range(2, 6):
        for c in (6, 8, 10, 12, 16, 20, 30):
            for sz in (2, 3):
                core = set([0]) | set(range(1, sz))
                for i in range(1, r):
                    core |= set(range(c*i, c*i + sz))
                core = sorted(core)
                if 4 <= len(core) <= 11 and primitive(core):
                    designs.append(core)
    seen = set()
    for core in designs:
        key = tuple(core)
        if key in seen:
            continue
        seen.add(key)
        r = num_scales(core)
        if r < 2:
            continue
        # scan w generously incl. resonances
        worst = F(0); ww = 0
        for w in resonant_ws(core):
            v = wDelta(core, w)
            if v > worst:
                worst = v; ww = w
        ratio = float(worst) / (r - 1)
        cur = best_by_r.get(r)
        if cur is None or ratio > cur[0]:
            best_by_r[r] = (ratio, tuple(core), ww, float(worst))
    print(f"  scanned {len(seen)} distinct multi-cluster cores")
    maxratio = 0.0
    for r in sorted(best_by_r):
        ratio, core, w, raw = best_by_r[r]
        maxratio = max(maxratio, ratio)
        flag = "  <-- EXCEEDS 2.71!" if ratio > 2.71 else ""
        print(f"    r={r}: worst w|Delta|={raw:.4f} at w={w}  ratio raw/(r-1)={ratio:.3f}"
              f"  core={core}{flag}")
    print(f"\n  => max observed slope (raw/(r-1)) = {maxratio:.3f}  "
          f"({'within 2.71' if maxratio <= 2.71 else 'EXCEEDS 2.71 -- bound questionable'})")
    return maxratio


# ======================================================================
# (4) FINITE-CHECK machinery verification
# ======================================================================
def part4_upward_closure(ntests=400, seed=11):
    print("\n" + "=" * 78)
    print("(4a) UPWARD-CLOSURE: E1 subset E2 => p0(E1) <= p0(E2)   (prune's basis)")
    print("=" * 78)
    rng = random.Random(seed)
    nbad = 0
    for _ in range(ntests):
        base = [0] + sorted(rng.sample(range(1, 40), rng.randint(4, 8)))
        add = rng.sample(range(1, 60), rng.randint(1, 4))
        E2 = sorted(set(base) | set(add))
        if p0_fast(E2) < p0_fast(base):
            nbad += 1
            if nbad <= 3:
                print(f"    VIOLATION base={base} add={add}: p0(super)<p0(sub)")
    print(f"  monotone (super >= sub): {ntests-nbad}/{ntests}  ({'OK' if nbad==0 else 'FAIL'})")
    return nbad == 0

def conservative_bnb(k, Tc, cap):
    """Independent re-implementation of the pruned exhaustive B&B (no node budget)."""
    universe = list(range(1, Tc))
    nU = len(universe)
    best = [F(0), None]
    n_leaf = [0]; n_pruned = [0]; n_cert = [0]; n_viol = [0]; viol = []
    prefix = [0]
    def dfs(start):
        need = k - len(prefix)
        if need == 0:
            E = tuple(prefix)
            if not primitive(E):
                return
            n_leaf[0] += 1
            v = p0_fast(E)
            if v > best[0]:
                best[0] = v; best[1] = E
            if v > cap:
                n_viol[0] += 1
                if len(viol) < 10:
                    viol.append((v, E))
            return
        rem = nU - start
        if rem < need:
            return
        if len(prefix) >= 3 and rem <= 28:
            ub = p0_fast(prefix + universe[start:])
            if ub <= cap:
                n_pruned[0] += 1
                n_cert[0] += comb(rem, need)
                return
        last = nU - need + 1
        for i in range(start, last):
            prefix.append(universe[i])
            dfs(i + 1)
            prefix.pop()
    dfs(0)
    return n_leaf[0], n_pruned[0], n_cert[0], best[0], best[1], n_viol[0], viol

def part4_finitecheck(info):
    print("\n" + "=" * 78)
    print("(4b) RE-RUN exhaustive finite check at C=3 threshold (independent)")
    print("=" * 78)
    allok = True
    for k in (8, 9, 10):
        cap, Q, margin, _Tc = info[k]
        T = F(3) / margin
        Tc = T.numerator // T.denominator + (1 if T.numerator % T.denominator else 0)
        nleaf, npr, ncert, mx, arg, nv, viol = conservative_bnb(k, Tc, cap)
        consec = tuple(range(k))
        is_consec = (arg == consec)
        ok = (nv == 0 and is_consec)
        allok = allok and ok
        print(f"  k={k}: C=3 T_k={Tc}  leaves={nleaf:,} cert-subtrees={npr:,} "
              f"(cover {ncert:,})")
        print(f"        max p0={float(mx):.5f} at {'consec' if is_consec else arg}  "
              f"viol={nv}  cap={float(cap):.5f}  {'PASS' if ok else 'FAIL'}")
        for v, E in sorted(viol, reverse=True)[:3]:
            print(f"          *** VIOLATION p0={float(v):.5f} at {E}")
    return allok


# ======================================================================
# (5) THE CRITICAL GAP CHECK: does the wide bulk really stay < cap?
#     The claim says span-monotone decay + Koksma handles wide sets at T_k=88/144/139.
#     We can't enumerate C(143,8), but we CAN look for high-p0 wide sets adversarially.
# ======================================================================
def part5_wide_hunt():
    print("\n" + "=" * 78)
    print("(5) WIDE-SET HUNT: is there a wide primitive set with p0 close to / above cap?")
    print("=" * 78)
    print("  Part D rests on 'wide => small p0' (sampled, not proved).  Adversarially")
    print("  search wide sets (large span, large max) for high p0.")
    rng = random.Random(31337)
    for k in (8, 9, 10):
        cap = CAP[k]
        worst = F(0); wa = None
        # random wide sets across the full peel band
        Tband = {8: 88, 9: 144, 10: 139}[k]
        for _ in range(60000):
            top = rng.randint(2 * k, Tband - 1)
            mid = rng.sample(range(1, top), k - 2)
            E = tuple(sorted(set([0] + mid + [top])))
            if len(E) != k or not primitive(E):
                continue
            v = p0_fast(E)
            if v > worst:
                worst = v; wa = E
        # also: structured wide sets (AP with large diff, near-consec shifted)
        for d in range(1, 18):
            for off in range(0, 5):
                E = tuple([0] + [off + d * i for i in range(1, k)])
                E = tuple(sorted(set(E)))
                if len(E) != k or not primitive(E):
                    continue
                v = p0_fast(E)
                if v > worst:
                    worst = v; wa = E
        print(f"  k={k}: cap={float(cap):.5f}  worst wide p0 found = {float(worst):.5f} "
              f"at {wa}  {'(<cap OK)' if worst <= cap else '*** >= cap!'}")
    return True


def main():
    print("ADVERSARIAL VERIFICATION: LRC(14) finite-check-to-threshold")
    print("EXACT Fraction arithmetic.\n")
    info, c_ok = part1_constants()
    dt_ok = part1b_dovetail()
    nviol = part2_hunt()
    slope = part3_slope_stress()
    mono_ok = part4_upward_closure()
    fc_ok = part4_finitecheck(info)
    part5_wide_hunt()

    print("\n" + "=" * 78)
    print("VERDICT")
    print("=" * 78)
    print(f"  (1)  exact constants reproduce            : {c_ok}")
    print(f"  (1b) dovetail identity exact + fast==ref  : {dt_ok}")
    print(f"  (2)  C(k)<=2.71*(k-2) violations found    : {nviol}  "
          f"({'NONE -> bound survives hunt' if nviol==0 else 'VIOLATED'})")
    print(f"  (3)  worst per-cluster slope observed     : {slope:.3f}  "
          f"({'<=2.71' if slope<=2.71 else '>2.71 EXCEEDS'})")
    print(f"  (4a) upward-closure monotonicity          : {mono_ok}")
    print(f"  (4b) C=3 exhaustive finite check (0 viol) : {fc_ok}")
    print()
    print("  NOTE: this verifier cannot CLOSE the imported gap (the non-resonant")
    print("  single-cluster Koksma constant 2.71 is numeric, not closed-form).  It can")
    print("  only fail to falsify it.  The finite check (4b) and constants (1) are")
    print("  rigorous and reproduce.  The load-bearing C(k) bound remains PARTIAL.")


if __name__ == "__main__":
    main()
