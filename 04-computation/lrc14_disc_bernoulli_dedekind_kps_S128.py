#!/usr/bin/env python3
"""
lrc14_disc_bernoulli_dedekind_kps_S128.py
=========================================
kind-pasteur-2026-07-13-S128.  THM-732 + HYP-6495.

THM-732(i) — EXACT IDENTITY. For G a finite union of intervals on R/Z with signed edge set E
(left endpoint sigma=+1, right endpoint sigma=-1), and any positive integer v:

    disc_v(G) := sum_{m != 0} |ghat(mv)|^2
              =  (1/(2 v^2)) * sum_{e,e' in E} sigma_e sigma_e' B2({v(e - e')}),

where B2(t) = t^2 - t + 1/6.  Equivalently (Poisson, as in THM-731):
    disc_v = (1/v) sum_{j<v} A(j/v) - |G|^2,   A(tau) = |G cap (G - tau)|.

Proof: ghat(k) = (1/(2 pi i k)) sum_e sigma_e e(-k e); expand |ghat(mv)|^2; use
sum_{m!=0} e(m theta)/m^2 = 2 pi^2 B2({theta}) (abs. convergent). QED.

THM-732(ii) — RATIONALITY. LRC good-set edges are rationals (14k+-1)/(14w) => disc_v in Q,
EXACTLY computable; THM-731's certificate condition  disc_v < 6 |G'_{~v}|^2  is decidable in
exact integer arithmetic => L>0 becomes a finite rational-arithmetic PROOF per family
(through the PROVED chain: peeling identity + |eps_v|^2 <= (6/49) disc_v).

THM-732(iii) — FAR-ELEMENT TAIL. |B2({x})| <= 1/6 => disc_v <= (2r)^2/(12 v^2) = r^2/(3 v^2).
Certificate holds whenever v > v0 := r*sqrt(2)/(6 |G'|), with r, |G'| those of the FIXED
leave-one-out set (independent of v!).  => entire far-element RAYS {base, v} close:
tail (v > v0) + finitely many exact checks (covering v <= v0) + THM-366 (non-covering v).

This script:
 (1) [MISTAKE-136 rule] verifies identity (i) EXACTLY — Fraction arithmetic, both sides,
     random interval unions AND the four covering families' leave-one-out good sets;
 (2) exact rational certificates for the 4 families (deep well / near-AP residue / {2..14} /
     variant), incl. exact peeling-identity and (*) consistency asserts;
 (3) the ray closures with explicit v0 and the finite exact-check lists;
 (4) HYP-6495 numerics: (w,w') block decomposition, exposure fractions, collapse diagnostics.

Everything below runs in EXACT rational arithmetic (fractions.Fraction); floats only for display.
"""
from fractions import Fraction as F
from math import sqrt
import random, sys

ONE = F(1)
SIXTH = F(1, 6)

# ----------------------------------------------------------------------------------------
# exact good-set machinery
# ----------------------------------------------------------------------------------------

def bad_pieces(w):
    """Arcs D_w = {||w t|| < 1/14} as linear pieces on [0,1], each tagged with generator w."""
    r = F(1, 14 * w)
    out = []
    for k in range(w):
        c = F(k, w)
        lo, hi = c - r, c + r
        if lo < 0:
            out.append((F(0), hi, w))
            out.append((lo + 1, ONE, w))
        elif hi > 1:
            out.append((lo, ONE, w))
            out.append((F(0), hi - 1, w))
        else:
            out.append((lo, hi, w))
    return out

def good_intervals(speeds):
    """Maximal good intervals of {t : ||w t|| >= 1/14 for all w in speeds}, circularly.
    Returns (goods, r, meas): goods = list of (a, b, wa, wb) with a<b (b may exceed 1 when the
    interval wraps 0), wa/wb = generating speeds of the left/right edge. Edges mod 1: a%1, b%1."""
    pieces = []
    for w in speeds:
        pieces += bad_pieces(w)
    pieces.sort(key=lambda p: (p[0], p[1]))
    comps = []  # [lo, hi, w_lo, w_hi]
    for lo, hi, w in pieces:
        if comps and lo <= comps[-1][1]:
            if hi > comps[-1][1]:
                comps[-1][1] = hi
                comps[-1][3] = w
        else:
            comps.append([lo, hi, w, w])
    assert comps
    n = len(comps)
    goods = []
    for i in range(n):
        a, wa = comps[i][1], comps[i][3]
        j = (i + 1) % n
        b, wb = comps[j][0] + (ONE if j == 0 else 0), comps[j][2]
        if a < b:
            goods.append((a, b, wa, wb))
    meas = sum(b - a for a, b, _, _ in goods)
    return goods, len(goods), meas

def edges_of(goods):
    """Signed edges mod 1: (position, sigma, generator_speed)."""
    E = []
    for a, b, wa, wb in goods:
        E.append((a if a < 1 else a - 1, +1, wa))
        E.append((b if b < 1 else b - 1, -1, wb))
    return E

def normalize_union(intervals):
    """Split wrap pieces to land in [0,1], return sorted disjoint (lo,hi)."""
    out = []
    for a, b in intervals:
        a %= 1
        bb = a + (b - a) if False else None
    # simpler: caller passes (a,b) with possible b>1
    out = []
    for a, b in intervals:
        if b <= 1:
            out.append((a, b))
        else:
            out.append((a, ONE))
            out.append((F(0), b - 1))
    out.sort()
    return out

def union_measure_intersection(U1, U2):
    """|U1 cap U2| for sorted disjoint interval lists on [0,1]."""
    tot = F(0)
    i = j = 0
    while i < len(U1) and j < len(U2):
        a = max(U1[i][0], U2[j][0])
        b = min(U1[i][1], U2[j][1])
        if a < b:
            tot += b - a
        if U1[i][1] < U2[j][1]:
            i += 1
        else:
            j += 1
    return tot

def autocorr_at(goods, tau):
    """A(tau) = |G cap (G - tau)| exactly."""
    U1 = normalize_union([(a, b) for a, b, _, _ in goods])
    shifted = []
    for a, b, _, _ in goods:
        aa, bb = a - tau, b - tau
        aa %= 1
        bb = aa + (b - a)
        shifted.append((aa, bb))
    U2 = normalize_union(shifted)
    return union_measure_intersection(U1, U2)

# ----------------------------------------------------------------------------------------
# THM-732(i): the two exact forms of disc_v
# ----------------------------------------------------------------------------------------

def B2(x):
    """B2({x}) for Fraction x."""
    t = x - (x.numerator // x.denominator)
    return t * t - t + SIXTH

def disc_pair(goods, v, want_blocks=False):
    """(1/(2v^2)) sum_{e,e'} s s' B2({v(e-e')}) — THM-732(i) RHS, exact."""
    E = edges_of(goods)
    S = F(0)
    blocks = {}
    for e, s, we in E:
        for f, u, wf in E:
            term = s * u * B2(v * (e - f))
            S += term
            if want_blocks:
                key = (we, wf)
                blocks[key] = blocks.get(key, F(0)) + term
    d = S / (2 * v * v)
    return (d, blocks) if want_blocks else d

def disc_def(goods, meas, v):
    """(1/v) sum_j A(j/v) - |G|^2 — the THM-731 definition, exact."""
    tot = F(0)
    for j in range(v):
        tot += autocorr_at(goods, F(j, v))
    return tot / v - meas * meas

# ----------------------------------------------------------------------------------------
# family-level exact certificate (THM-731 chain + THM-732)
# ----------------------------------------------------------------------------------------

def eps_exact(goods, meas, v):
    """eps_v = |D_v cap G'| - (1/7)|G'| exactly."""
    U1 = normalize_union([(a, b) for a, b, _, _ in goods])
    Dv = normalize_union([(a, b) for a, b, _ in bad_pieces(v)])
    inter = union_measure_intersection(U1, Dv)
    return inter - meas / 7

def certify(base, v, verbose=True, check_def=False, want_blocks=False):
    """Full exact certificate for family base+[v], peeling v. Returns dict of exact quantities."""
    goods, r, meas = good_intervals(base)
    dp = disc_pair(goods, v, want_blocks=want_blocks)
    if want_blocks:
        dp, blocks = dp
    else:
        blocks = None
    assert dp >= 0, "positivity violated -> formula bug"
    if check_def:
        dd = disc_def(goods, meas, v)
        assert dd == dp, "THM-732(i) FAILS: pair form != definition (%s vs %s)" % (dd, dp)
    eps = eps_exact(goods, meas, v)
    # (*) exact:
    star_ok = eps * eps <= F(6, 49) * dp
    # peeling identity exact (needs full-family good measure):
    _, _, Lfull = good_intervals(list(base) + [v])
    peel_ok = (Lfull == F(6, 7) * meas - eps)
    # certificate: 6|G'|^2 - disc_v > 0  <=>  L_cert > 0
    margin = 6 * meas * meas - dp
    Lcert = float(F(6, 7) * meas) - sqrt(float(F(6, 49) * dp))
    return dict(goods=goods, r=r, meas=meas, disc=dp, eps=eps, Lfull=Lfull,
                star_ok=star_ok, peel_ok=peel_ok, margin=margin, Lcert=Lcert, blocks=blocks)

# ----------------------------------------------------------------------------------------
# (1) identity verification — MISTAKE-136 rule: verify the identity ITSELF, exactly
# ----------------------------------------------------------------------------------------

print("=" * 100)
print("THM-732(i) IDENTITY VERIFICATION (exact Fraction arithmetic, both sides)")
print("=" * 100)

random.seed(20260713)
ok_all = True
for trial in range(8):
    # random disjoint rational intervals on the circle
    npts = random.randint(2, 5) * 2
    pts = sorted(random.sample([F(a, 840) for a in range(840)], npts))
    goods = []
    for i in range(0, npts, 2):
        goods.append((pts[i], pts[i + 1], 0, 0))
    meas = sum(b - a for a, b, _, _ in goods)
    for v in (5, 7, 12):
        lhs = disc_def(goods, meas, v)
        rhs = disc_pair(goods, v)
        match = (lhs == rhs)
        ok_all &= match
        if not match:
            print("  TOY FAIL trial=%d v=%d: def=%s pair=%s" % (trial, v, lhs, rhs))
print("  toy unions (8 random unions x v in {5,7,12}): EXACT equality on all -> %s" % ok_all)
assert ok_all

FAMS = [
    ("deep well {1..12,182}",        list(range(1, 13)),            182),
    ("near-AP residue {1..11,13,84}", list(range(1, 12)) + [13],     84),
    ("small-far {2..14}",             list(range(2, 14)),            14),
    ("variant {1,3..13,182}",         [1] + list(range(3, 14)),      182),
]

for name, base, v in FAMS:
    goods, r, meas = good_intervals(base)
    lhs = disc_def(goods, meas, v)
    rhs = disc_pair(goods, v)
    print("  %-34s peel v=%-4d r=%-3d  def==pair: %s   disc_v = %s ~ %.6e"
          % (name, v, r, lhs == rhs, rhs, float(rhs)))
    assert lhs == rhs, "identity fails on %s" % name
print("  => THM-732(i) VERIFIED exactly (definition == Bernoulli edge-pair sum, as rationals).")

# ----------------------------------------------------------------------------------------
# (2) exact rational certificates for the four covering families
# ----------------------------------------------------------------------------------------

print()
print("=" * 100)
print("EXACT RATIONAL CERTIFICATES (THM-731 chain, all quantities in Q; floats display-only)")
print("  need: margin := 6|G'_{~v}|^2 - disc_v > 0   (<=> L_cert=(6/7)|G'|-sqrt((6/49)disc) > 0)")
print("=" * 100)

for name, base, v in FAMS:
    res = certify(base, v, want_blocks=True)
    print("\n%s   (peel v=%d)" % (name, v))
    print("   |G'_{~v}| = %s ~ %.6f    r = %d intervals   L_full = %s ~ %.6f"
          % (res['meas'], float(res['meas']), res['r'], res['Lfull'], float(res['Lfull'])))
    print("   disc_v    = %s ~ %.6e" % (res['disc'], float(res['disc'])))
    print("   eps_v     = %s ~ %.6e" % (res['eps'], float(res['eps'])))
    print("   peeling identity L=(6/7)|G'|-eps EXACT: %s ;  (*) |eps|^2<=(6/49)disc EXACT: %s"
          % (res['peel_ok'], res['star_ok']))
    print("   CERT margin 6|G'|^2 - disc = %s ~ %.3e  ->  %s   (L_cert ~ %.6f, true L ~ %.6f)"
          % (res['margin'], float(res['margin']),
             "POSITIVE: L>0 PROVED (exact arithmetic)" if res['margin'] > 0 else "not certified",
             res['Lcert'], float(res['Lfull'])))
    assert res['peel_ok'] and res['star_ok']

# ----------------------------------------------------------------------------------------
# (3) THM-732(iii): far-element tail + ray closures
# ----------------------------------------------------------------------------------------

print()
print("=" * 100)
print("RAY CLOSURES: {base, v} for ALL v.  tail: v > v0 = r*sqrt(2)/(6|G'|)  (disc<=r^2/(3v^2)<6|G'|^2)")
print("  below v0: covering v = exact checks; non-covering v = THM-366 (t=1/q witness).")
print("=" * 100)

RAYS = [
    ("{1..12, v}      (deep-well ray)",  list(range(1, 13)),        182),  # covering iff 182|v
    ("{1..11,13, v}   (residue ray)",    list(range(1, 12)) + [13],  84),  # covering iff 84|v
    ("{2..13, v}",                       list(range(2, 14)),         14),  # covering iff 14|v
    ("{1,3..13, v}",                     [1] + list(range(3, 14)),   14),  # covering iff 14|v
]

for name, base, step in RAYS:
    goods, r, meas = good_intervals(base)
    v0sq = F(r * r, 18) / (meas * meas)
    v0 = sqrt(float(v0sq))
    print("\n%s   r=%d, |G'|=%.6f  =>  v0 ~ %.1f" % (name, r, float(meas), v0))
    # exact checks: covering v = step*j <= v0
    j = 1
    all_ok = True
    checked = []
    while step * j * step * j <= v0sq:
        v = step * j
        if v in base:
            j += 1
            continue
        dp = disc_pair(goods, v)
        m = 6 * meas * meas - dp
        checked.append((v, m > 0, float(m)))
        all_ok &= (m > 0)
        j += 1
    # also certify the first ray point beyond v0 as a sanity illustration of the tail
    vtail = step * j
    dptail = disc_pair(goods, vtail)
    crude = F(r * r, 3) / (vtail * vtail)
    print("   exact checks (covering v<=v0): %s" %
          (", ".join("v=%d:%s(margin %.2e)" % (v, "OK" if ok else "FAIL", m) for v, ok, m in checked)
           if checked else "none needed"))
    print("   tail sanity at first v>v0 (v=%d): disc=%.3e <= crude r^2/(3v^2)=%.3e < 6|G'|^2=%.3e : %s"
          % (vtail, float(dptail), float(crude), float(6 * meas * meas),
             dptail <= crude < 6 * meas * meas))
    assert dptail <= crude < 6 * meas * meas
    status = "CLOSED (all v)" if all_ok else "NOT closed"
    print("   RAY %s: %s  [tail v>v0 by THM-732(iii); covering v<=v0 exact; non-covering v THM-366]"
          % (name.split()[0], status))

# ----------------------------------------------------------------------------------------
# (4) HYP-6495: block decomposition + exposure + collapse diagnostics
# ----------------------------------------------------------------------------------------

print()
print("=" * 100)
print("HYP-6495 NUMERICS: (w,w') blocks of the pair sum S = 2 v^2 disc_v ; exposure; collapse")
print("=" * 100)

for name, base, v in FAMS[:2]:
    res = certify(base, v, want_blocks=True)
    goods, r, meas = res['goods'], res['r'], res['meas']
    E = edges_of(goods)
    print("\n%s  (peel v=%d): 2r=%d edges" % (name, v, len(E)))
    # exposure fraction per speed
    from collections import Counter
    cnt = Counter(w for _, _, w in E)
    exposure = ", ".join("w=%d:%d/%d" % (w, cnt.get(w, 0), 2 * w) for w in sorted(base))
    print("   exposed edges per speed (of 2w endpoints): %s" % exposure)
    # top blocks
    blocks = res['blocks']
    tot = sum(blocks.values())
    diag = sum(val for (wa, wb), val in blocks.items() if wa == wb)
    print("   S = 2 v^2 disc = %.4f ; diagonal-blocks (w=w') sum = %.4f ; off = %.4f"
          % (float(tot), float(diag), float(tot - diag)))
    top = sorted(blocks.items(), key=lambda kv: -abs(float(kv[1])))[:8]
    print("   top blocks by |mass|: " + ", ".join("(%d,%d):%.3f" % (wa, wb, float(val))
                                                  for (wa, wb), val in top))
    # collapse diagnostic: crude vs actual
    crude = F((2 * r) ** 2, 6)
    print("   crude bound on S = (2r)^2/6 = %.1f ; actual S = %.4f ; collapse factor %.1fx"
          % (float(crude), float(tot), float(crude) / float(tot) if tot else float('inf')))
    print("   disc*3v^2/r^2 = %.4f (1.0 = crude saturated) ; disc*v^2*6/(2r) = %.4f (~1 = pure-diagonal scale)"
          % (float(res['disc'] * 3 * v * v / (r * r)), float(res['disc'] * v * v * 6 / (2 * r))))

print("\ndone.")
