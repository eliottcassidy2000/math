#!/usr/bin/env python3
"""
lrc_gap_filters_kps_S3.py -- HYP-4105 (kind-pasteur-2026-07-05-S3)

Verification for the S3 bricks (all integer arithmetic, fast):

(P1) CERTIFICATE COMPLETENESS constants: margin beta' at real t* => atom cert
     (round(s t*), s, ceil(beta s)) for s >= B/(2(beta'-beta)).  Random check.
(P2) GAP FILTER logic on the only known sub-2/25 sets (dilated APs):
     covering q<=12 in every unit direction + mod-13/14 near-unit pinning.
(P3) MERGE-DENOMINATOR EXCLUSION (THM-592): exact M is attained on the grid
     d/(v+w); gap window (1/13, 2/25) forces d >= 3, v+w >= 38 => w_max >= 19.
     Verify grid attainment vs fine scan; print the (d, v+w) gap-value table.
(P4) FILTERED GAP CENSUS [1, H]: all primitive 12-subsets passing
     {covering 2..12} + {mod-13 pinning shape} + {two elements summing >= 38}:
     exact M for each survivor; confirm NONE lands in (1/13, 2/25).
"""

from math import gcd, ceil
from fractions import Fraction
from itertools import combinations
import random

# ---------------------------------------------------------------- exact M (int)
def exact_M(ws):
    """M(S) via the THM-592 merge grid d/(v+w) (+ peaks d/(2v) via i=j).
    Integer arithmetic; returns a Fraction."""
    L = sorted(set(ws))
    best_n, best_d = 0, 1  # best margin as fraction best_n/best_d
    qs = set()
    n = len(L)
    for i in range(n):
        for j in range(i, n):
            qs.add(L[i] + L[j])
    for q in qs:
        # margin at t = d/q for all d: precompute residues
        for d in range(1, q // 2 + 1):
            # NOTE: do NOT skip gcd(d,q)>1 -- the binding time m/(v+w) need not be
            # in lowest terms (e.g. {6,10,18,54}: t = 1/4 = 4/16, pair 6+10)
            m = q  # min over runners of min(r, q-r)
            for w in ws:
                r = (w * d) % q
                r = r if r <= q - r else q - r
                if r < m:
                    m = r
                    if m * best_d <= best_n * q:  # cannot beat current best
                        break
            if m * best_d > best_n * q:
                best_n, best_d = m, q
    return Fraction(best_n, best_d)

def margin_at_frac(ws, p, q):
    """min margin at t = p/q, as Fraction."""
    m = q
    for w in ws:
        r = (w * p) % q
        r = min(r, q - r)
        if r < m:
            m = r
    return Fraction(m, q)

# ---------------------------------------------------------------- P1
def verify_completeness(trials=2000, seed=7):
    rng = random.Random(seed)
    fails = 0
    for _ in range(trials):
        kk = rng.randint(3, 12)
        ws = sorted(rng.sample(range(1, 300), kk))
        B = max(ws)
        t_star = Fraction(rng.randint(1, 999), 1000)
        bp = margin_at_frac(ws, t_star.numerator, t_star.denominator)
        if bp == 0:
            continue
        beta = bp * Fraction(4, 5)
        delta = bp - beta
        s = int(Fraction(B, 2 * delta)) + 1
        k = round(s * t_star)
        mu = ceil(beta * s)
        ok = all(mu <= (w * k) % s <= s - mu for w in ws)
        if not ok:
            fails += 1
            print(f"  P1 FAIL: ws={ws} t*={t_star} bp={bp} s={s}")
    print(f"[P1] certificate completeness: {trials} random instances, {fails} failures")
    dd = Fraction(14, 169) - Fraction(2, 25)
    print(f"[P1b] rigidity-slack instantiation: delta = 14/169 - 2/25 = {dd}; "
          f"s >= B/(2 delta) = {Fraction(1, 2) / dd} B ~= {float(Fraction(1,2)/dd):.1f} B")
    return fails == 0

# ---------------------------------------------------------------- P2
def verify_gap_filters():
    """The only sets with M < 2/25 are (by the verified rigidity) the dilated APs.
    Verify the filter conclusions hold on them (sanity of the Lean lemmas)."""
    fails = 0
    for c in [1, 2, 3, 5, 7, 13, 20]:
        ws = [c * j for j in range(1, 13)]
        M = exact_M(ws)
        assert M == Fraction(1, 13), (c, M)
        # covering filter: q <= 12, all units a
        for q in range(2, 13):
            for a in range(1, q):
                if gcd(a, q) > 1:
                    continue
                if not any((w * a) % q == 0 for w in ws):
                    fails += 1
                    print(f"  P2 FAIL covering: c={c} q={q} a={a}")
        # 13/14 near-unit pinning
        for q in (13, 14):
            for a in range(1, q):
                if gcd(a, q) > 1:
                    continue
                if not any((w * a) % q in (0, 1, q - 1) for w in ws):
                    fails += 1
                    print(f"  P2 FAIL pinning: c={c} q={q} a={a}")
    print(f"[P2] gap filters on dilated APs c*{{1..12}} (the entire sub-2/25 locus): "
          f"{fails} violations")
    return fails == 0

# ---------------------------------------------------------------- P3
def verify_merge_form(trials=400, seed=13):
    rng = random.Random(seed)
    bad = 0
    for _ in range(trials):
        kk = rng.randint(4, 12)
        ws = sorted(rng.sample(range(1, 60), kk))
        M = exact_M(ws)
        # independent check: fine scan over all p/q, q <= 130, cannot beat the grid
        fine = Fraction(0)
        for q in range(2, 131):
            for p in range(1, q):
                if gcd(p, q) > 1:
                    continue
                m = margin_at_frac(ws, p, q)
                if m > fine:
                    fine = m
        if fine > M:
            bad += 1
            print(f"  P3 FAIL: fine beats grid on {ws}: {fine} > {M}")
    print(f"[P3] merge-grid attainment (grid >= fine scan q<=130): {trials} sets, "
          f"{bad} failures")
    print("[P3b] gap-window merge values d/(v+w) in (1/13, 2/25), v+w <= 60 "
          "(NOTE: none with v+w < 38; d >= 3 always):")
    lo, hi = Fraction(1, 13), Fraction(2, 25)
    rows = []
    for ssum in range(2, 61):
        for d in range(1, ssum // 2 + 1):
            val = Fraction(d, ssum)
            if lo < val < hi:
                rows.append((d, ssum, val))
    for d, ssum, val in rows:
        print(f"      d={d}, v+w={ssum}: {val} = {float(val):.6f}")
    dmin = min(r[0] for r in rows)
    smin = min(r[1] for r in rows)
    print(f"      => minimal depth d = {dmin}, minimal sum v+w = {smin} "
          f"(gap violator's binding pair sums >= {smin} => w_max >= {(smin + 1) // 2})")
    return bad == 0

# ---------------------------------------------------------------- P4
def filtered_census(H=24):
    total = 0
    pass_cov = 0
    passed = 0
    gap_hits = []
    lo, hi = Fraction(1, 13), Fraction(2, 25)
    Ms = {}
    for ws in combinations(range(1, H + 1), 12):
        g = 0
        for w in ws:
            g = gcd(g, w)
            if g == 1:
                break
        if g != 1:
            continue
        total += 1
        # covering (a=1): every q in 2..12 divides some element
        cov = True
        for q in range(2, 13):
            if not any(w % q == 0 for w in ws):
                cov = False
                break
        if not cov:
            continue
        pass_cov += 1
        # mod-13 pinning shape, all units
        ok13 = all(any((w * a) % 13 in (0, 1, 12) for w in ws) for a in range(1, 13))
        if not ok13:
            continue
        # merge exclusion: two largest must sum >= 38
        srt = sorted(ws)
        if srt[-1] + srt[-2] < 38:
            continue
        passed += 1
        M = exact_M(ws)
        Ms[ws] = M
        if lo < M < hi:
            gap_hits.append((ws, M))
            print(f"  P4 GAP HIT: {ws} M={M}")
    print(f"[P4] census [1,{H}]: {total} primitive 12-subsets; {pass_cov} covering; "
          f"{passed} pass ALL gap filters; {len(gap_hits)} in the gap (MUST be 0)")
    if Ms:
        mn = min(Ms.values())
        args = [w for w, m in Ms.items() if m == mn]
        print(f"      min M among filter-passers: {mn} at {args[:3]}")
    return len(gap_hits) == 0

if __name__ == "__main__":
    print("=" * 74)
    print("HYP-4105: cert completeness + gap filters + merge exclusion")
    print("=" * 74)
    ok1 = verify_completeness()
    ok2 = verify_gap_filters()
    ok3 = verify_merge_form()
    ok4 = filtered_census()
    print("=" * 74)
    print(f"SUMMARY: P1 {'OK' if ok1 else 'FAIL'}; P2 {'OK' if ok2 else 'FAIL'}; "
          f"P3 {'OK' if ok3 else 'FAIL'}; P4 {'OK' if ok4 else 'FAIL'}")
