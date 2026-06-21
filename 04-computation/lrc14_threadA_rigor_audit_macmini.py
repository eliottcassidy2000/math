"""
THREAD A adversarial rigor audit of the LRC(14) ledger (mac-mini, 2026-06-21).

INDEPENDENT re-derivation (NOT reusing the existing lrc14/lrc_q108 helpers) of:
  (1) caps cap_k = min_{|P|=13-k} meas(G_P), G_P={x: ||p x||>=1/14 forall p in P}.
  (2) meas(S7(E)) via exact rational breakpoints; consecutive argmax check (L3).
  (3) the 2D cell discrepancy D_{p,q} = sum_ij |mu(i,j)-1/49|, and the claimed
      bounds D<=14/p, D<=12/(7q), sup(D*p), sup(D*q).
  (4) the plateau P2(B) and the resonance correction R(p/q)=p0_inf-P2, |R|<=D.
  (5) the apex D=0 law: D_{p,q}=0 <=> 7|pq.

All exact (Fraction).  Adversarial: look for a finite check that DISAGREES with
the canon numbers.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations


# ---------- (1) caps -------------------------------------------------------
def meas_GP(P, half=F(1, 14)):
    """meas{x in [0,1): ||p x|| >= half for all p in P}.
    ||t|| = dist to nearest integer.  ||p x||>=h  <=>  frac(px) in [h,1-h].
    Breakpoints from each p: x where px hits an integer +- h.
    """
    if not P:
        return F(1)
    bps = set([F(0), F(1)])
    for p in P:
        # frac(px) = h  or  = 1-h  =>  px = m+h or m+1-h
        for m in range(0, p + 1):
            for off in (half, 1 - half):
                x = F(m + off, 1) / p
                if 0 <= x <= 1:
                    bps.add(x)
    bps = sorted(bps)
    total = F(0)
    for a, b in zip(bps, bps[1:]):
        mid = (a + b) / 2
        ok = True
        for p in P:
            fr = (p * mid) - int(p * mid)
            if fr < 0:
                fr += 1
            if not (half <= fr <= 1 - half):
                ok = False
                break
        if ok:
            total += b - a
    return total


def caps():
    """cap_k = min over |P|=13-k of meas(G_P), P subset {1..13}."""
    out = {}
    universe = list(range(1, 14))
    for k in range(8, 14):
        psz = 13 - k
        best = None
        for P in combinations(universe, psz):
            m = meas_GP(list(P))
            if best is None or m < best:
                best = m
        out[k] = best
    return out


# ---------- (2) meas(S7(E)) -----------------------------------------------
def sector(fr):
    # fr in [0,1); returns 0..6
    s = int(fr * 7)
    return 6 if s >= 7 else s


def meas_S7(E):
    """meas{x: {sector(frac(e x)): e in E} == {0..6}}."""
    bps = set([F(0), F(1)])
    for e in E:
        if e == 0:
            continue
        # sector wall: frac(e x) = j/7  => e x = m + j/7
        for m in range(0, e + 1):
            for j in range(0, 7):
                x = F(7 * m + j, 7 * e)
                if 0 <= x <= 1:
                    bps.add(x)
    bps = sorted(bps)
    total = F(0)
    for a, b in zip(bps, bps[1:]):
        mid = (a + b) / 2
        hit = set()
        for e in E:
            fr = (e * mid) - int(e * mid)
            if fr < 0:
                fr += 1
            hit.add(sector(fr))
        if len(hit) == 7:
            total += b - a
    return total


def L3_consec_argmax(k, maxE):
    """Among primitive E with 0 in E, |E|=k, max(E)<=maxE, check consecutive
    {0..k-1} is the argmax of meas(S7).  Return (consec_val, any_beater)."""
    consec = list(range(k))
    cv = meas_S7(consec)
    beater = None
    # iterate over E = {0} u (k-1)-subset of {1..maxE}
    for rest in combinations(range(1, maxE + 1), k - 1):
        E = (0,) + rest
        g = E[1]
        for x in E[2:]:
            g = gcd(g, x)
        if g != 1:
            continue
        v = meas_S7(list(E))
        if v > cv:
            beater = (E, v)
            break
    return cv, beater


# ---------- (3) discrepancy D_{p,q} ---------------------------------------
def cell_law(p, q):
    """mu(i,j) = meas{v in [0,1): (sector(q v), sector(p v)) = (i,j)}.
    breakpoints: q v in m+j/7 and p v in m+j/7.
    Returns 7x7 list of Fractions summing to 1."""
    bps = set([F(0), F(1)])
    for c in (q, p):
        for m in range(0, c + 1):
            for j in range(0, 7):
                x = F(7 * m + j, 7 * c)
                if 0 <= x <= 1:
                    bps.add(x)
    bps = sorted(bps)
    mu = [[F(0)] * 7 for _ in range(7)]
    for a, b in zip(bps, bps[1:]):
        mid = (a + b) / 2
        fq = (q * mid) - int(q * mid)
        if fq < 0:
            fq += 1
        fp = (p * mid) - int(p * mid)
        if fp < 0:
            fp += 1
        mu[sector(fq)][sector(fp)] += b - a
    return mu


def D_pq(p, q):
    mu = cell_law(p, q)
    s = F(0)
    for i in range(7):
        for j in range(7):
            s += abs(mu[i][j] - F(1, 49))
    return s


# ---------- (4) plateau & resonance ---------------------------------------
def g_B(B):
    """g_B(i,j) = x-measure{ {sector(b x): b in B} u {i,j} == {0..6} }.
    Returns 7x7 of Fractions."""
    bps = set([F(0), F(1)])
    for b in B:
        if b == 0:
            continue
        for m in range(0, b + 1):
            for j in range(0, 7):
                x = F(7 * m + j, 7 * b)
                if 0 <= x <= 1:
                    bps.add(x)
    bps = sorted(bps)
    # for each cell, accumulate length where base-hit set u {i,j} = full
    g = [[F(0)] * 7 for _ in range(7)]
    for a, b2 in zip(bps, bps[1:]):
        mid = (a + b2) / 2
        hit = set()
        for b in B:
            fr = (b * mid) - int(b * mid)
            if fr < 0:
                fr += 1
            hit.add(sector(fr))
        length = b2 - a
        missing = set(range(7)) - hit
        for i in range(7):
            for j in range(7):
                if missing <= {i, j}:  # base + {i,j} covers all
                    g[i][j] += length
    return g


def P2(B):
    g = g_B(B)
    s = F(0)
    for i in range(7):
        for j in range(7):
            s += g[i][j]
    return s / 49


def p0_inf(B, p, q):
    mu = cell_law(p, q)
    g = g_B(B)
    s = F(0)
    for i in range(7):
        for j in range(7):
            s += mu[i][j] * g[i][j]
    return s


# ---------- run ------------------------------------------------------------
if __name__ == "__main__":
    print("=== (1) caps cap_k = min_{|P|=13-k} meas(G_P) ===")
    canon_caps = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91),
                  11: F(66, 91), 12: F(6, 7), 13: F(1)}
    cp = caps()
    for k in range(8, 14):
        match = "OK" if cp[k] == canon_caps[k] else "MISMATCH"
        print(f"  cap_{k} = {cp[k]} = {float(cp[k]):.5f}  canon {canon_caps[k]}  [{match}]")

    print("\n=== (2/L3) consec argmax of meas(S7), small finite check ===")
    for k, maxE in [(8, 13), (9, 12), (10, 12)]:
        cv, beater = L3_consec_argmax(k, maxE)
        print(f"  k={k} maxE={maxE}: meas(S7(consec))={cv}={float(cv):.5f}  "
              f"cap_{k}={float(cp[k]):.5f}  consec<=cap:{cv<=cp[k]}  "
              f"beater:{beater}")

    print("\n=== (3) discrepancy D_{p,q}, bounds 14/p and 12/(7q) ===")
    ratios = [(2, 1), (3, 2), (4, 3), (5, 3), (5, 4), (7, 4), (7, 6), (8, 5),
              (9, 5), (15, 7), (28, 25), (7, 5), (6, 5), (8, 7)]
    viol_14p = 0
    sup_Dp = F(0)
    sup_Dq = F(0)
    for (p, q) in ratios:
        if gcd(p, q) != 1:
            continue
        D = D_pq(p, q)
        b14 = F(14, p)
        b12 = F(12, 7 * q)
        ok14 = D <= b14
        sup_Dp = max(sup_Dp, D * p)
        sup_Dq = max(sup_Dq, D * q)
        if not ok14:
            viol_14p += 1
        print(f"  p/q={p}/{q}: D={D}={float(D):.5f}  14/p={float(b14):.4f}[{ok14}]  "
              f"D*p={float(D*p):.4f}  D*q={float(D*q):.4f}  7|pq:{(p%7==0 or q%7==0)}  D=0:{D==0}")
    print(f"  violations of D<=14/p: {viol_14p}")
    print(f"  sup(D*p) over sample = {sup_Dp}={float(sup_Dp):.5f} (canon claim 20/7={float(F(20,7)):.5f})")
    print(f"  sup(D*q) over sample = {sup_Dq}={float(sup_Dq):.5f} (canon claim 12/7={float(F(12,7)):.5f})")

    print("\n=== (3b) exhaustive D<=14/p over all coprime 1<p/q<=43/20, q<=40 ===")
    viol = 0
    sup_dp = F(0)
    sup_dq = F(0)
    tested = 0
    for q in range(1, 41):
        for p in range(q + 1, q * 43 // 20 + 1):
            if gcd(p, q) != 1:
                continue
            if F(p, q) > F(43, 20):
                continue
            D = D_pq(p, q)
            tested += 1
            if D > F(14, p):
                viol += 1
                print(f"    VIOLATION p/q={p}/{q} D={D} > 14/p={F(14,p)}")
            sup_dp = max(sup_dp, D * p)
            sup_dq = max(sup_dq, D * q)
            # apex law check
            apex = (p % 7 == 0 or q % 7 == 0)
            if (D == 0) != apex:
                print(f"    APEX LAW VIOLATION p/q={p}/{q} D=0?{D==0} 7|pq?{apex}")
    print(f"  tested {tested} coprime ratios q<=40 in (1,43/20]; D<=14/p violations: {viol}")
    print(f"  sup(D*p)={sup_dp}={float(sup_dp):.5f}; sup(D*q)={sup_dq}={float(sup_dq):.5f}")

    print("\n=== (4) plateau P2(B) and resonance R(p/q)=p0_inf-P2, |R|<=D ===")
    bases = {8: [0, 2, 4, 6, 8, 10], 9: [0, 2, 4, 6, 8, 10, 12],
             10: [0, 2, 4, 6, 8, 10, 12, 14]}
    for k, B in bases.items():
        plat = P2(B)
        print(f"  k={k} B={B}: P2={plat}={float(plat):.5f}  cap_{k}={float(cp[k]):.5f}  "
              f"margin cap-P2={float(cp[k]-plat):.5f}")
        worstR = F(0)
        for (p, q) in [(2, 1), (3, 2), (4, 3), (5, 3), (5, 4), (7, 4), (7, 6)]:
            pinf = p0_inf(B, p, q)
            R = pinf - plat
            D = D_pq(p, q)
            ok = abs(R) <= D
            worstR = max(worstR, abs(R))
            below_cap = pinf <= cp[k]
            print(f"     p/q={p}/{q}: p0_inf={float(pinf):.5f} R={float(R):.5f} "
                  f"|R|<=D({float(D):.4f}):{ok} p0_inf<=cap:{below_cap}")
        print(f"     worst |R| = {float(worstR):.5f}")
