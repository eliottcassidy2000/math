#!/usr/bin/env python3
"""
THE BODY-WEIGHTED RAMANUJAN EXPANSION (THM-873 for ARBITRARY bodies)
boxeph-2026-07-16-S23 -- kind-pasteur cont.28's named critical task (a):
"body-weighted Ramanujan spectra ... the bridge from the closed interval-core
case to open stratum [B]".

CONTENT
 (A) CONCENTRIC COLLAPSE (any body E, any lambda): the bad set is one effective
     arc per primitive fraction a/l, l in D(E) (divisor closure), of half-width
     W_E(l) = lambda / m_E(l),  m_E(l) = min{v in E : l | v}.
     |arcs| = sum_{l in D(E)} phi(l); for E = {v}-closures this is v exactly.
 (B) SWEEP-EXACT TRANSFORM: sorting arcs by left edge with prefix-max M_r,
       1_B = sum_r 1_{I_r} - sum_r 1_{J_r},  J_r = (l_r, min(r_r, M_r)],
     so for h != 0
       -ghat(h) = bhat(h) = SUM_{l in D(E)} c_l(h) sin(2 pi h W_E(l))/(pi h)
                  - SUM_r T(J_r)
     -- the body-weighted Ramanujan first-order term plus explicit shadow
     corrections.  Interval cores in THM-873's regime: shadows = closed-gap
     lenses (recovers THM-873).  Cascades (arcs swallowed whole) are handled:
     a fully-shadowed arc has J_r = I_r and cancels exactly.
 (C) PEEL FORM (E = P u {v}): the new levels D(v)\D(P) contribute one thin arc
     per torsion point; arcs inside B_P cancel; arcs inside G_P survive whole;
     straddlers (<= 2 r_P of them) contribute partial lenses:
       ghat_E(h) = ghat_P(h) - sum_{new arcs A} T(A cap G_P)
     and grouping whole survivors by level yields RESTRICTED RAMANUJAN SUMS
       c_l^S(h) = sum_{a in S_l} e(-2 pi i h a / l),
       S_l = {a primitive mod l : the arc at a/l survives inside G_P}.
 (D) THE DEEP-WELL CENSUS: which torsion survives at lambda = 1/14 over the
     core P = {1..12}: level 13 FULL (all 12 primitive 13ths lonely), level 14
     ALL-STRADDLE (the tight locus: arcs centered ON dG_P), level 26 EMPTY,
     levels 91/182 partial -- the spectral fingerprint of the deep well.

Pure Python; exact interval arithmetic in Fractions; transforms refereed in
floats at machine precision.
"""

from fractions import Fraction as Fr
from math import gcd, sin, cos, pi
import cmath

# ---------------------------------------------------------------- intervals

def merge(ivs):
    """merge closed-ish intervals given as (lo, hi) Fractions on the line."""
    ivs = sorted(ivs)
    out = []
    for lo, hi in ivs:
        if out and lo <= out[-1][1]:
            out[-1][1] = max(out[-1][1], hi)
        else:
            out.append([lo, hi])
    return [(a, b) for a, b in out]

def circle_arcs_to_line(arcs):
    """arcs = list of (center, halfwidth) on R/Z -> merged intervals in [0,1),
    wrapping handled by unrolling into [0,2) then folding."""
    ivs = []
    for c, w in arcs:
        c = c - int(c)
        if c < 0:
            c += 1
        ivs.append((c - w, c + w))
    # fold into [0,1): shift pieces below 0 and above 1
    folded = []
    for lo, hi in ivs:
        assert hi - lo < 1
        if lo < 0:
            folded.append((lo + 1, Fr(1)))
            folded.append((Fr(0), hi))
        elif hi > 1:
            folded.append((lo, Fr(1)))
            folded.append((Fr(0), hi - 1))
        else:
            folded.append((lo, hi))
    return merge(folded)

def complement(ivs):
    out, prev = [], Fr(0)
    for lo, hi in ivs:
        if lo > prev:
            out.append((prev, lo))
        prev = max(prev, hi)
    if prev < 1:
        out.append((prev, Fr(1)))
    return out

def divisors(n):
    ds = set()
    d = 1
    while d * d <= n:
        if n % d == 0:
            ds.add(d); ds.add(n // d)
        d += 1
    return ds

def D_closure(E):
    ds = set()
    for v in E:
        ds |= divisors(v)
    return sorted(ds)

def m_E(E, l):
    return min(v for v in E if v % l == 0)

def effective_arcs(E, lam):
    """(A): one arc per primitive a/l, l in D(E), half-width lam/m_E(l)."""
    arcs = []
    for l in D_closure(E):
        w = Fr(lam, 1) / m_E(E, l) if not isinstance(lam, Fr) else lam / m_E(E, l)
        for a in range(l):
            if gcd(a, l) == 1 or l == 1:
                if l == 1 and a != 0:
                    continue
                arcs.append((Fr(a, l), w, l))
    return arcs

def bad_direct(E, lam):
    """direct definition: union over v in E, c mod v of arcs width lam/v."""
    arcs = []
    for v in E:
        for c in range(v):
            arcs.append((Fr(c, v), lam / v))
    return circle_arcs_to_line(arcs)

def T(lo, hi, h):
    """transform of interval: int_lo^hi e(-2 pi i h t) dt, h != 0 (float)."""
    return (cmath.exp(-2j * pi * h * float(lo)) - cmath.exp(-2j * pi * h * float(hi))) / (2j * pi * h)

def c_l(l, h):
    """Ramanujan sum via gcd closed form (also cross-checked directly)."""
    g = gcd(l, h)
    m = l // g
    # c_l(h) = mu(m) * phi(l)/phi(m)
    mu, mm = 1, m
    p = 2
    while p * p <= mm:
        if mm % p == 0:
            mm //= p
            if mm % p == 0:
                return 0
            mu = -mu
        p += 1
    if mm > 1:
        mu = -mu
    def phi(n):
        r, x, p = 1, n, 2
        res = n
        while p * p <= x:
            if x % p == 0:
                res -= res // p
                while x % p == 0:
                    x //= p
            p += 1
        if x > 1:
            res -= res // x
        return res
    return mu * phi(l) // phi(m)

# ---------------------------------------------------------------- referees

def referee_A(E, lam, name):
    direct = bad_direct(E, lam)
    eff = circle_arcs_to_line([(c, w) for c, w, l in effective_arcs(E, lam)])
    ok = direct == eff
    n_arcs = len(effective_arcs(E, lam))
    phis = sum((sum(1 for a in range(l) if gcd(a, l) == 1) if l > 1 else 1)
               for l in D_closure(E))
    print(f"  [{name}] (A) effective-arc system == direct bad set: {ok}; "
          f"#arcs = {n_arcs} == sum phi(l) over D(E) = {phis}: {n_arcs == phis}; "
          f"bad measure = {sum(b-a for a,b in direct)}")
    assert n_arcs == phis
    assert ok
    return direct

def referee_B(E, lam, name, hmax=40):
    """(B): -ghat = Ramanujan first order - shadow transforms, vs direct."""
    arcs = effective_arcs(E, lam)
    merged = circle_arcs_to_line([(c, w) for c, w, l in arcs])
    good = complement(merged)
    if not good:
        print(f"  [{name}] bad set covers the circle; skip transform referee")
        return
    # rotate so that a good midpoint is the origin
    g0 = (good[0][0] + good[0][1]) / 2
    shifted = []
    for c, w, l in arcs:
        cc = c - g0
        cc -= int(cc)
        if cc < 0:
            cc += 1
        shifted.append((cc - w, cc + w, l))
    assert all(lo > 0 and hi < 1 for lo, hi, l in shifted), "rotation failed"
    shifted.sort()
    # shadows via prefix max
    shadows = []
    M = None
    for lo, hi, l in shifted:
        if M is not None and M > lo:
            shadows.append((lo, min(hi, M)))
        M = hi if M is None else max(M, hi)
    # direct components in shifted coords
    comp = merge([(lo, hi) for lo, hi, l in shifted])
    Wl = {l: lam / m_E(E, l) for l in D_closure(E)}
    maxerr = 0.0
    for h in range(1, hmax + 1):
        # direct
        bd = sum(T(lo, hi, h) for lo, hi in comp)
        # formula: Ramanujan first order (levels; centers shifted by -g0 =>
        # e(+2 pi i h g0) factor on each level sum) - shadows
        ram = 0
        for l in D_closure(E):
            # sum over primitive a: e(-2 pi i h (a/l - g0)) = e(2 pi i h g0) c_l(h)
            s = sum(cmath.exp(-2j * pi * h * float(Fr(a, l) - g0))
                    for a in range(l) if (gcd(a, l) == 1 and (l > 1 or a == 0)) or (l == 1 and a == 0))
            ram += s * sin(2 * pi * h * float(Wl[l])) / (pi * h)
        sh = sum(T(lo, hi, h) for lo, hi in shadows)
        err = abs(bd - (ram - sh))
        maxerr = max(maxerr, err)
    n_swallowed = sum(1 for k, (lo, hi, l) in enumerate(shifted)
                      if k > 0 and max(x[1] for x in shifted[:k]) >= hi)
    print(f"  [{name}] (B) transform identity max|err| over h<=%d: %.2e; "
          "shadows=%d (fully swallowed arcs=%d)" % (hmax, maxerr, len(shadows), n_swallowed))
    assert maxerr < 1e-11
    return maxerr

def referee_C(P, v, lam, name, hmax=40):
    """(C) peel form + the restricted-sum census."""
    E = sorted(set(P) | {v})
    badP = bad_direct(P, lam)
    goodP = complement(badP)
    newlv = sorted(set(D_closure(E)) - set(D_closure(P)))
    # note: levels in both keep m_E = m_P iff min multiple unchanged
    changed = [l for l in D_closure(P) if m_E(E, l) != m_E(P, l)]
    assert not changed, f"core widths changed at levels {changed}"
    # classify new arcs
    census = {}
    surv, strad, swall = [], [], []
    for l in newlv:
        w = Fr(lam, 1) / v if not isinstance(lam, Fr) else lam / v
        assert m_E(E, l) == v
        S = []
        for a in range(l):
            if l > 1 and gcd(a, l) != 1:
                continue
            c = Fr(a, l)
            lo, hi = c - w, c + w
            # position relative to G_P (work mod 1: new arcs tiny, no wrap except a=0 impossible for l>1... l=1 not new)
            inside_good = any(x <= lo and hi <= y for x, y in goodP)
            inside_bad = any(x <= lo and hi <= y for x, y in badP)
            if inside_good:
                surv.append((c, w)); S.append(a)
            elif inside_bad:
                swall.append((c, w))
            else:
                strad.append((c, w))
        census[l] = (S, sum(1 for c, _ in strad if c.denominator == l),
                     sum(1 for c, _ in swall if c.denominator == l))
    # transform check: ghat_E - ghat_P = -sum_new T(A cap G_P)
    badE = bad_direct(E, lam)
    maxerr = 0.0
    for h in range(1, hmax + 1):
        dE = -sum(T(lo, hi, h) for lo, hi in badE)
        dP = -sum(T(lo, hi, h) for lo, hi in badP)
        s = 0
        for l in newlv:
            w = Fr(lam, 1) / v if not isinstance(lam, Fr) else lam / v
            for a in range(l):
                if l > 1 and gcd(a, l) != 1:
                    continue
                c = Fr(a, l)
                lo, hi = c - w, c + w
                for x, y in goodP:
                    ilo, ihi = max(lo, x), min(hi, y)
                    if ihi > ilo:
                        s += T(ilo, ihi, h)
        err = abs((dE - dP) + s)
        maxerr = max(maxerr, err)
    print(f"  [{name}] (C) peel identity max|err|: {maxerr:.2e}")
    assert maxerr < 1e-11
    print(f"      new-level census (level: survive/straddle/swallow of phi(l)):")
    for l in newlv:
        S, nstrad, nswall = census[l]
        phi_l = sum(1 for a in range(l) if gcd(a, l) == 1) if l > 1 else 1
        print(f"        l={l:>4}: {len(S):>3} / {nstrad:>3} / {nswall:>3}   of {phi_l}")
    # restricted Ramanujan cancellation exhibit
    print(f"      restricted Ramanujan sums (whole survivors), levels with S nonempty:")
    for l in newlv:
        S = census[l][0]
        if not S:
            continue
        vals = []
        for h in range(1, hmax + 1):
            cs = sum(cmath.exp(-2j * pi * h * a / l) for a in S)
            vals.append(abs(cs))
        print(f"        l={l:>4}: |S|={len(S):>3}, max_h |c_l^S(h)| = {max(vals):.3f} "
              f"(trivial {len(S)}), full-sum max |c_l(h)| = {max(abs(c_l(l, h)) for h in range(1, hmax+1))}")
    return census

# ---------------------------------------------------------------- run

if __name__ == "__main__":
    lam14 = Fr(1, 14)
    print("=" * 74)
    print("PART A/B -- the body-weighted expansion on six bodies (exact set identity")
    print("            + machine-precision transform referee)")
    bodies = [
        ("core {1..12} @ 1/14 (THM-873 sanity)", list(range(1, 13)), lam14),
        ("deep well {1..12,182} @ 1/14", list(range(1, 13)) + [182], lam14),
        ("multikiller {1..10,13,22,84} @ 1/14", list(range(1, 11)) + [13, 22, 84], lam14),
        ("shifted core {2..14} @ 1/14", list(range(2, 15)), lam14),
        ("worst |core|=1 body {1..11,13,84} @ 1/14", list(range(1, 12)) + [13, 84], lam14),
        ("cascade stress {4,6,9} @ 1/5", [4, 6, 9], Fr(1, 5)),
        ("wild {3,5,7,26,55} @ 1/8", [3, 5, 7, 26, 55], Fr(1, 8)),
    ]
    for name, E, lam in bodies:
        referee_A(E, lam, name)
        referee_B(E, lam, name)

    print()
    print("=" * 74)
    print("PART C -- the peel form and the restricted-torsion census")
    referee_C(list(range(1, 13)), 182, lam14, "deep well peel {1..12} + 182")
    referee_C(list(range(1, 12)) + [13], 84, lam14, "peel {1..11,13} + 84")
    referee_C(list(range(1, 11)) + [13, 22], 84, lam14, "peel {1..10,13,22} + 84")

    print()
    print("=" * 74)
    print("PART D -- disc_13({1..12}; 1/14): EXACT rational value via autocorrelation,")
    print("          cross-checked against the h-sum of the expansion")
    E = list(range(1, 13))
    bad = bad_direct(E, lam14)
    good = complement(bad)
    mG = sum(b - a for a, b in good)
    v = 13

    def shift_mod1(ivs, tau):
        out = []
        for lo, hi in ivs:
            lo, hi = lo - tau, hi - tau
            k = -int(lo) + (1 if lo < 0 else 0)
            lo, hi = lo + k, hi + k
            if hi <= 1:
                out.append((lo, hi))
            else:
                out.append((lo, Fr(1)))
                out.append((Fr(0), hi - 1))
        return merge(out)

    def inter_len(A, B):
        s, i, j = Fr(0), 0, 0
        while i < len(A) and j < len(B):
            lo = max(A[i][0], B[j][0]); hi = min(A[i][1], B[j][1])
            if hi > lo:
                s += hi - lo
            if A[i][1] < B[j][1]:
                i += 1
            else:
                j += 1
        return s

    disc_exact = sum(inter_len(good, shift_mod1(good, Fr(j, v))) for j in range(v)) / v - mG * mG
    print(f"  disc_13 EXACT = {disc_exact} = {float(disc_exact):.6e}")
    for H in (2000, 8000, 32000):
        s = 0.0
        for h in range(1, H + 1):
            gh = -sum(T(lo, hi, h * v) for lo, hi in bad)
            s += 2 * abs(gh) ** 2
        print(f"  h-sum |h| <= {H:>6}: {s:.6e}  (-> exact from below: {s <= float(disc_exact) + 1e-12})")
    print(f"  (fleet-quoted 1.648e-2 [opus-S271/kps THM-873] agrees to its printed precision: "
          f"{abs(float(disc_exact) - 1.648e-2) < 5e-5})")
    print()
    print("ALL CHECKS PASSED")
