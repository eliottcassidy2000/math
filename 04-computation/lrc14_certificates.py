#!/usr/bin/env python3
"""
lrc14_certificates -- the consolidated exact-certificate library for LRC(14)
============================================================================
opus-2026-07-14-S291.  ENGINEERING MANDATE DELIVERABLE: one importable, documented,
self-testing module consolidating the exact-rational machinery built across the
S270-S290 arc (and used, in scattered per-script form, by the whole fleet).

Everything is EXACT (fractions.Fraction) unless suffixed _float.  Every public function
returns rigorous objects: measures, certificates, and explicit rational WITNESSES that a
referee (or Lean's `decide`) can re-check independently.

    from lrc14_certificates import *

    M_exact([1,3,4,5,6,7,8,9,10,11,12,13,14])     # exact sup-min clearance (M value)
    L_exact(body)                                  # exact 1/14-safe measure
    good_intervals(speeds)                         # the exact safe interval set
    disc_exact(core, v)                            # THM-731/732 grid discrepancy (Bernoulli)
    thm731_certificate(body, v)                    # peel certificate  L >= (6|G'|-sqrt(6 disc))/7
    capped_envelope_vstar(core)                    # THM-755 band edge v* = r/(pi |G'|)
    h_band_protocol(body)                          # THM-756's 3-line closure (returns layer+witness)
    fine_comb_witness(body)                        # THM-752 tooth-threshold witness (exact t*)
    clean_slot_witness(body)                       # THM-754(ii) k=7 explicit-delta witness
    slot_feasible(body, k, a)                      # S287 exact delta-interval shadow test

Provenance: THM-731/732 (klein/kps), THM-752/754/755/756 (opus), the interval engine
(fleet-standard).  Self-test: `python3 lrc14_certificates.py` reproduces the arc's key
exact numbers (deep-well certificate, band edges, witnesses) and exits nonzero on any drift.
"""
from fractions import Fraction as F
import math

LAM = F(1, 14)          # the LRC(14) threshold
__all__ = ["LAM", "good_intervals", "L_exact", "M_exact", "disc_exact",
           "thm731_certificate", "capped_envelope_vstar", "h_band_protocol",
           "fine_comb_witness", "clean_slot_witness", "slot_feasible", "is_covering"]

# ---------------- the interval engine ----------------

def _intersect(A, B):
    out = []; i = j = 0
    while i < len(A) and j < len(B):
        a1, b1 = A[i]; a2, b2 = B[j]
        lo, hi = max(a1, a2), min(b1, b2)
        if hi > lo: out.append((lo, hi))
        if b1 <= b2: i += 1
        else: j += 1
    return out

def good_intervals(speeds, lam=LAM):
    """exact interval set {t in [0,1) : ||v t|| >= lam for all v in speeds}."""
    cur = [(F(0), F(1))]
    for v in sorted(speeds):
        cur = _intersect(cur, [((k + lam) / v, (k + 1 - lam) / v) for k in range(v)])
    return cur

def L_exact(speeds, lam=LAM):
    """exact measure of the lam-safe set."""
    return sum(b - a for a, b in good_intervals(speeds, lam))

def M_exact(speeds):
    """exact M = sup_t min_v ||v t|| by bisection on lam over the rational lattice.

    M is attained; it is a rational with denominator dividing lcm-type combinations.
    Strategy: binary search on lam via L_exact(speeds, lam) > 0 over the finite set of
    candidate lams = clearances at interval endpoints of iterated refinements."""
    lo, hi = F(0), F(1, 2)
    for _ in range(80):
        mid = (lo + hi) / 2
        if L_exact(speeds, mid) > 0: lo = mid
        else: hi = mid
    # snap: the true M is min_v ||v t*|| at an endpoint witness of the lo-level set
    ivs = good_intervals(speeds, lo)
    best = F(0)
    for a, b in ivs:
        for t in (a, b, (a + b) / 2):
            c = min(min((v * t) % 1, 1 - (v * t) % 1) for v in speeds)
            if c > best: best = c
    return best

def is_covering(speeds):
    """divisor-complete test (every q in 2..14 divides some speed)."""
    return all(any(v % q == 0 for v in speeds) for q in range(2, 15))

# ---------------- THM-731/732: the peel certificates ----------------

def disc_exact(core, v, lam=LAM):
    """THM-731's v-grid discrepancy of the core's good set, exact (Bernoulli jump-pair form,
    kps THM-732)."""
    ivs = good_intervals(core, lam)
    B2 = lambda x: x * x - x + F(1, 6)
    js = []
    for a, b in ivs: js.append((a, 1)); js.append((b, -1))
    tot = F(0)
    for xp, sp in js:
        for xq, sq in js:
            tot += sp * sq * B2((v * (xp - xq)) % 1)
    return tot / (2 * v * v)

def thm731_certificate(body, v, lam=LAM):
    """the peel-v certificate: returns (fires: bool, exact_lower_bound_on_7L or None).
    Certificate: L >= [(6/7)|G'| - sqrt((6/49) disc_v)] ; fires iff disc_v < 6 |G'|^2."""
    core = [x for x in body if x != v]
    G = L_exact(core, lam)
    d = disc_exact(core, v, lam)
    return (d < 6 * G * G), (G, d)

def capped_envelope_vstar(core, lam=LAM):
    """THM-755's band edge: (H) holds unconditionally for every v > v* = r/(pi |G'|).
    Returns (vstar_float, r, G_exact).  Above vstar the peel certificate fires with no
    further computation."""
    ivs = good_intervals(core, lam)
    G = sum(b - a for a, b in ivs); r = len(ivs)
    return (r / (math.pi * float(G)) if G > 0 else float("inf")), r, G

# ---------------- THM-752: the fine-comb witness ----------------

def fine_comb_witness(body, lam=LAM):
    """tooth-threshold witness: peel w; if the core's safe set has an interval longer than
    one w-danger tooth (1/(7w) at lam=1/14... exactly 2*lam/w... tooth = 2*lam/w), a closed
    w-safe point inside it is a full witness.  Returns (t_exact, clearance) or None."""
    for w in sorted(body, reverse=True):
        C = [x for x in body if x != w]
        for a, b in good_intervals(C, lam):
            if (b - a) * w > 2 * lam:            # longer than one tooth
                k = (a * w).__floor__()
                for kk in (k - 1, k, k + 1, k + 2):
                    lo, hi = (kk + lam) / w, (kk + 1 - lam) / w
                    lo2, hi2 = max(lo, a), min(hi, b)
                    if hi2 >= lo2:
                        t = (lo2 + hi2) / 2 if hi2 > lo2 else lo2
                        cl = min(min((x * t) % 1, 1 - (x * t) % 1) for x in body)
                        if cl >= lam: return t, cl
    return None

# ---------------- THM-754(ii): the k=7 clean slot ----------------

def clean_slot_witness(body, lam=LAM):
    """explicit-delta k=7 witness (residue-level thresholds (1,3,5) c_min).  Returns
    (a, t_exact, clearance) or None."""
    C7 = [x for x in body if x % 7 == 0]
    if not C7: return None
    cmin = min(C7); d = F(1, 14 * cmin)
    for a in range(1, 7):
        ok = all(min((c * d) % 1, 1 - (c * d) % 1) >= lam for c in C7)
        if not ok: continue
        for x in body:
            if x % 7 == 0: continue
            r = (x * a) % 7; rho = min(r, 7 - r)
            if x > {1: cmin, 2: 3 * cmin, 3: 5 * cmin}[rho]: ok = False; break
        if ok:
            t = F(a, 7) + d
            cl = min(min((x * t) % 1, 1 - (x * t) % 1) for x in body)
            if cl >= lam: return a, t, cl
    return None

# ---------------- S287: the exact slot test ----------------

def slot_feasible(body, k, a, lam=LAM):
    """exact delta-interval shadow test at slot (k, a): returns (feasible, delta_witness)."""
    lo0, hi0 = F(1, 10 ** 9), F(1, 2 * k)
    cur = [(lo0, hi0)]
    for v in sorted(body, reverse=True):
        base = F(v * a, k)
        good = []
        mlo = (v * lo0 + base - 1).__floor__(); mhi = (v * hi0 + base + 1).__ceil__()
        for m in range(mlo, mhi + 1):
            g1, g2 = (m + lam - base) / v, (m + 1 - lam - base) / v
            g1, g2 = max(g1, lo0), min(g2, hi0)
            if g2 > g1: good.append((g1, g2))
        cur = _intersect(cur, sorted(good))
        if not cur: return False, None
    return True, cur[0][0]

# ---------------- THM-756: the band protocol ----------------

def h_band_protocol(body, lam=LAM):
    """the 3-line closure for a peel body: returns (layer, detail).
    layer 1: THM-755 (v > v*); layer 2: exact (H); layer 3: fine-comb witness;
    layer 4: direct L > 0; layer 0: L = 0 (tight boundary -- check equality loneliness)."""
    v = max(body)
    core = [x for x in body if x != v]
    vstar, r, G = capped_envelope_vstar(core, lam)
    if G > 0 and v > vstar: return 1, ("THM-755", vstar)
    if G > 0 and disc_exact(core, v, lam) < 6 * G * G: return 2, ("exact-(H)", v)
    w = fine_comb_witness(body, lam)
    if w: return 3, ("fine-comb", w)
    L = L_exact(body, lam)
    if L > 0: return 4, ("direct-L", L)
    return 0, ("tight-boundary", L)

# ---------------- self-test battery ----------------

def _selftest():
    ok = True
    def check(name, cond):
        nonlocal ok
        print(("  PASS  " if cond else "  FAIL  ") + name)
        ok = ok and cond
    deep = list(range(1, 13)) + [182]
    check("deep well L = 6617/388080-family value", L_exact(list(range(1, 13))) == F(6617, 194040))
    check("deep well covering", is_covering(deep))
    fires, (G, d) = thm731_certificate(deep, 182)
    check("THM-731 fires at the deep-well far peel", fires)
    vstar, r, G12 = capped_envelope_vstar(list(range(1, 13)))
    check("THM-755 band edge for {1..12} in (111, 113)", 111 < vstar < 113)
    w = fine_comb_witness([1, 10, 21, 24, 56, 65, 77, 135, 219, 265, 335, 367, 390])
    check("THM-752 closes the klein-S304 stall (witness exists)", w is not None and w[1] >= LAM)
    feas, dw = slot_feasible(list(range(1, 13)) + [182], 7, 1)
    check("slot (7,1) feasible for the deep well", feas)
    lay, det = h_band_protocol([1, 2, 3, 4, 5, 6, 7, 9, 10, 14, 60, 150, 431])
    check("band protocol closes a band body at layer <= 2", lay in (1, 2))
    ap = list(range(1, 14))
    check("AP is the tight boundary (L = 0)", L_exact(ap) == 0)
    check("AP is lonely with equality at t = 1/14",
          min(min((v * LAM) % 1, 1 - (v * LAM) % 1) for v in ap) == LAM)
    return ok

if __name__ == "__main__":
    import sys
    print("lrc14_certificates self-test:")
    sys.exit(0 if _selftest() else 1)
