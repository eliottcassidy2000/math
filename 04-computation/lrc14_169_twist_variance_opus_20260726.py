#!/usr/bin/env python3
"""THM-2334 Section 7 cheapest decisive step: the exact 169-twist variance test.

Object (THM-2334, eqs (28)-(31), (41)-(42), specialized per THM-2349 to a
first-depth-one row at the k=2 delayed clock, R = 13^k = 169 = 13^(lambda_j+1)):

    H(ell) = delta_hat(m) * e_13(m ell_c)
             * hat{ prod_i I_i(w_i t + ell_i/13) J_i(R w_i t + R ell_i/13) }(X)
             * conj( hat{ prod_i I_i(w_i t + ell_i/13) }(Y) ),

for the 169 characters ell in G^ = L_13^perp / <w> (THM-2309 (15a)-(15b),
THM-2334 (35)-(37)).  Because 13 | R, the word translations R ell_i / 13 are
integers (THM-2334 (44)), so the word factor is twist-neutral and

    H(ell) = C0 * gamma(ell),
    C0     = delta_hat(m) / (4 pi^2 X Y)   (real, nonzero),
    gamma(ell) = e_13(m ell_8) * A_X(ell) * conj(B_Y(ell)),

where A_X(ell) = (2 pi i X) hat{E_Q^ell}(X) and B_Y(ell) = (2 pi i Y)
hat{E^ell}(Y) are CYCLOTOMIC INTEGERS (finite +/- sums of roots of unity of
order dividing NN below), because E^ell and E_Q^ell = E^ell * (1_Q o T^2) are
finite unions of rational intervals.

Decision object (THM-2334 (42)-(43)):

    sum_{q != 0} |A(q)|^2 = (1/169) sum_ell |H(ell) - Hbar|^2 > 0
        <=>  the bank {H(ell)} is NOT constant
        <=>  some ell has gamma(ell) != gamma(0).

Exact decision, no floats: gamma(ell) - gamma(0) lies in Z[zeta_NN].  We map
Z[zeta_NN] -> F_p by a ring homomorphism zeta_NN |-> h, where p = c*NN + 1 is
prime and h has exact multiplicative order NN mod p.  A NONZERO image proves
the difference is a nonzero algebraic number, hence H(ell) != H(0), hence the
variance is strictly positive.  (A zero image proves nothing; we then consult
a second independent homomorphism and report the outcome without a claim.)
Floats appear only in descriptive magnitudes, never in the decision.

Canonical row/word stratum (THM-2349 ordering (21), scalars = the canon's
typed row THM-2309 eq (25), live strict valuation profile (1,3,5)):

    guard H = 1 (odd 13-unit), units q = (14, 27, 40, 53, 66),
    blockers c_1 = 13 (owner j, depth 1), c_2 = 13^3, c_3 = 2*13^5 (deepest),
    w = (1, 14, 27, 40, 53, 66, 13, 2197, 742586).

    C_H = {x : ||H x|| > 1/7},   D_v = {x : ||v x|| < 1/14},
    A_0 = C_H \\ U_i D_{q_i},
    E_1 = A_0 n D_{c_1} \\ (D_{c_2} u D_{c_3})            (THM-2349 (22)),
    R_1 words (THM-2305 (3)):  Q_{1,{a}} = A_0 n D_{c_2} \\ (D_{c_1} u D_{c_3}),
                               Q_{1,{b}} = A_0 n D_{c_3} \\ (D_{c_1} u D_{c_2}),
                               Q_{1,{a,b}} = A_0 n D_{c_2} n D_{c_3} \\ D_{c_1},
    stratum E_{1,sigma,2} = E_1 n T^{-2} Q_{1,sigma},  T x = 13 x  (THM-2349 (31)-(32)).

THM-2309 eq (25) is typed but NOT asserted to be a scalar cover; no cover
hypothesis is needed anywhere below (THM-2309 Sec. 1).  The result is a
statement about the canonical typed row's marked current, not about a
hypothetical covering row.

Script: 04-computation/lrc14_169_twist_variance_opus_20260726.py
Output: 05-knowledge/results/lrc14_169_twist_variance_opus_20260726.out
"""

import cmath
import math
import sys
import time
from bisect import bisect_right
from fractions import Fraction
from math import gcd

T_START = time.time()


def log(msg=""):
    print(msg, flush=True)


def elapsed():
    return f"[{time.time() - T_START:8.1f}s]"


# --------------------------------------------------------------------------
# Row data and scales
# --------------------------------------------------------------------------
W = (1, 14, 27, 40, 53, 66, 13, 13**3, 2 * 13**5)
GUARD, OWNER, TA, TB = 0, 6, 7, 8
UNIT_IDX = (1, 2, 3, 4, 5)
K_CLOCK = 2
RDIL = 13**K_CLOCK  # 169


def nu13(n):
    v = 0
    while n % 13 == 0:
        n //= 13
        v += 1
    return v


assert W[GUARD] % 2 == 1 and W[GUARD] % 13 != 0, "guard must be odd 13-unit"
assert all(W[i] % 13 != 0 for i in UNIT_IDX), "units must be 13-units"
assert len({W[i] for i in UNIT_IDX}) == 5, "units distinct"
PROFILE = (nu13(W[OWNER]), nu13(W[TA]), nu13(W[TB]))
assert PROFILE == (1, 3, 5), PROFILE
assert 5 <= PROFILE[2] <= 19 and 1 <= PROFILE[1] < PROFILE[2], "THM-2349 universe"
assert RDIL % 13 == 0, "THM-2334 (44): word translations trivial"

LCM_W = 1
for v in W:
    LCM_W = LCM_W * v // gcd(LCM_W, v)
T_DEN = 182 * LCM_W                    # common denominator of all t-breakpoints
NN = RDIL * T_DEN                      # cyclotomic order for all exponentials
assert all(T_DEN % (182 * v) == 0 for v in W)
assert T_DEN % 91 == 0 and NN % 13 == 0

NN_PRIMES = []
_rem = NN
for q in (2, 3, 5, 7, 11, 13, 53):
    if _rem % q == 0:
        NN_PRIMES.append(q)
        while _rem % q == 0:
            _rem //= q
assert _rem == 1, f"unexpected factor left in NN: {_rem}"

TWO_PI = 2.0 * math.pi


def cn(e):
    """Numeric value of zeta_NN^e (e already reduced mod NN). Descriptive only."""
    return cmath.exp(2j * math.pi * (e / NN))


def fmtc(z):
    return f"({z.real:.6e}{z.imag:+.6e}j)"


# --------------------------------------------------------------------------
# Modular embeddings  Z[zeta_NN] -> F_p,  zeta_NN -> h  (exact order NN)
# --------------------------------------------------------------------------
def is_prime(n):
    if n < 2:
        return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % sp == 0:
            return n == sp
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def find_embedding(skip_c=()):
    c = 1
    while True:
        if c not in skip_c:
            p = c * NN + 1
            if is_prime(p):
                for g in range(2, 200):
                    h = pow(g, (p - 1) // NN, p)
                    if h == 1:
                        continue
                    if all(pow(h, NN // q, p) != 1 for q in NN_PRIMES):
                        return c, p, h
        c += 1


# --------------------------------------------------------------------------
# Exact interval machinery (integer endpoints on the T_DEN scale)
# --------------------------------------------------------------------------
def in_comb(i, elli):
    """Sorted disjoint intervals of {t : ||w_i t + elli/13|| < 1/14}."""
    wsp = W[i]
    U = T_DEN // (182 * wsp)
    lo = (-13 - 14 * elli) % 182
    out = []
    for n in range(wsp):
        s = (lo + 182 * n) * U
        t = s + 26 * U
        if t <= T_DEN:
            out.append((s, t))
        else:
            out.append((s, T_DEN))
            out.append((0, t - T_DEN))
    out.sort()
    return out


def subtract_comb(iv, wsp, PD, lo, hi):
    """Subtract the periodic windows ((lo+PD*n)/(PD*wsp),(hi+PD*n)/(PD*wsp))."""
    U = T_DEN // (PD * wsp)
    step = PD * U
    lenW = (hi - lo) * U
    base = (lo % PD) * U
    out = []
    ap = out.append
    for A, B in iv:
        p0 = A - lenW + 1
        k0 = -((base - p0) // step)          # ceil((p0-base)/step)
        p = base + k0 * step
        cur = A
        while p < B:
            we = p + lenW
            if we > cur:
                if p > cur:
                    ap((cur, p))
                cur = we
                if cur >= B:
                    break
            p += step
        if cur < B:
            ap((cur, B))
    return out


def build_set(pattern, ell):
    """pattern: dict index -> 'in' | 'out' | 'gout'; ell: 9-tuple of ints."""
    ins = [i for i, mmode in pattern.items() if mmode == "in"]
    start = min(ins, key=lambda i: W[i])
    iv = in_comb(start, ell[start])
    for i, mmode in pattern.items():
        if mmode == "gout":
            iv = subtract_comb(iv, W[i], 91, -13 - 7 * ell[i], 13 - 7 * ell[i])
    rest = sorted((W[i], i) for i, mmode in pattern.items()
                  if i != start and mmode in ("in", "out"))
    for _, i in rest:
        e = ell[i]
        if pattern[i] == "out":
            iv = subtract_comb(iv, W[i], 182, -13 - 14 * e, 13 - 14 * e)
        else:  # intersect with 'in' comb == subtract its complement windows
            iv = subtract_comb(iv, W[i], 182, 13 - 14 * e, 169 - 14 * e)
    return iv


def check_intervals(iv):
    last = -1
    for a, b in iv:
        assert 0 <= a < b <= T_DEN and a >= last, "interval list corrupt"
        last = b
    return sum(b - a for a, b in iv)


PAT_E = {0: "gout", 1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
         6: "in", 7: "out", 8: "out"}
PAT_QA = {0: "gout", 1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
          6: "out", 7: "in", 8: "out"}
PAT_QB = {0: "gout", 1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
          6: "out", 7: "out", 8: "in"}
PAT_QAB = {0: "gout", 1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
           6: "out", 7: "in", 8: "in"}
ZELL = (0,) * 9


# --------------------------------------------------------------------------
# X-side sweep:  A_X = (2 pi i X) hat{ E^ell * (1_Q o T^2) }(X)
#   pieces of E^ell n T^{-2}Q enumerated exactly; endpoints on the NN scale.
# Returns (per-embedding sums, complex sum, exact overlap length on NN scale).
# --------------------------------------------------------------------------
def x_sweep(E_iv, Q_iv, Q_starts, X, mods, tabs, want_values=True):
    nQ = len(Q_iv)
    nm = len(mods)
    S_mod = [0] * nm
    S_c = 0.0j
    ov = 0
    if want_values:
        wfp = [pow(mi[1], (-X * T_DEN) % NN, mi[0]) for mi in mods]
        wfm = [pow(mi[1], (X * T_DEN) % NN, mi[0]) for mi in mods]
        wfp_c = cn((-X * T_DEN) % NN)
        wfm_c = cn((X * T_DEN) % NN)
    for A, B in E_iv:
        LA = 169 * A
        sA = LA % T_DEN
        span = 169 * (B - A)
        assert span < T_DEN
        sE = sA + span
        idx = bisect_right(Q_starts, sA) - 1
        off = 0
        if idx < 0:
            idx = nQ - 1
            off = -T_DEN
        if want_values:
            be = (-X * (LA - sA)) % NN
            b_mod = [pow(mi[1], be, mi[0]) for mi in mods]
            b_c = cn(be)
            acc = [0] * nm
            acc_c = 0.0j
            wf = [1] * nm
            wf_c = 1.0 + 0.0j
            if off:
                wf = list(wfm)
                wf_c = wfm_c
        while True:
            qa0, qb0 = Q_iv[idx]
            qa = qa0 + off
            qb = qb0 + off
            if qa >= sE:
                break
            if qb > sA:
                lo = sA if qa < sA else qa
                hi = sE if qb > sE else qb
                if hi > lo:
                    ov += hi - lo
                    if want_values:
                        for j in range(nm):
                            p, h = mods[j]
                            if lo == qa:
                                vlo = tabs[j][0][idx] * wf[j] % p
                            else:
                                vlo = pow(h, (-X * lo) % NN, p)
                            if hi == qb:
                                vhi = tabs[j][1][idx] * wf[j] % p
                            else:
                                vhi = pow(h, (-X * hi) % NN, p)
                            acc[j] = (acc[j] + vlo - vhi) % p
                        vlo_c = tabs[-1][0][idx] * wf_c if lo == qa else cn((-X * lo) % NN)
                        vhi_c = tabs[-1][1][idx] * wf_c if hi == qb else cn((-X * hi) % NN)
                        acc_c += vlo_c - vhi_c
            idx += 1
            if idx == nQ:
                idx = 0
                off += T_DEN
                if want_values:
                    for j in range(nm):
                        wf[j] = wf[j] * wfp[j] % mods[j][0]
                    wf_c *= wfp_c
        if want_values:
            for j in range(nm):
                S_mod[j] = (S_mod[j] + b_mod[j] * acc[j]) % mods[j][0]
            S_c += b_c * acc_c
    return S_mod, S_c, ov


def make_tabs(Q_iv, X, mods):
    tabs = []
    for p, h in mods:
        tlo = [pow(h, (-X * a) % NN, p) for a, _ in Q_iv]
        thi = [pow(h, (-X * b) % NN, p) for _, b in Q_iv]
        tabs.append((tlo, thi))
    tlo_c = [cn((-X * a) % NN) for a, _ in Q_iv]
    thi_c = [cn((-X * b) % NN) for _, b in Q_iv]
    tabs.append((tlo_c, thi_c))
    return tabs


# --------------------------------------------------------------------------
# Y-side:  sum over endpoints of e(-F t)  ->  (2 pi i F) hat{set}(F)
# --------------------------------------------------------------------------
def endpoint_sum(iv, F, mods):
    nm = len(mods)
    S_mod = [0] * nm
    S_c = 0.0j
    for a, b in iv:
        ea = (-F * 169 * a) % NN
        eb = (-F * 169 * b) % NN
        for j in range(nm):
            p, h = mods[j]
            S_mod[j] = (S_mod[j] + pow(h, ea, p) - pow(h, eb, p)) % p
        S_c += cn(ea) - cn(eb)
    return S_mod, S_c


# --------------------------------------------------------------------------
# GF(13) linear algebra for G^ = L_13^perp / <w>
# --------------------------------------------------------------------------
def gf13_nullspace(rows):
    m = [list(r) for r in rows]
    ncol = 9
    piv = []
    r = 0
    for ccol in range(ncol):
        sel = None
        for rr in range(r, len(m)):
            if m[rr][ccol] % 13:
                sel = rr
                break
        if sel is None:
            continue
        m[r], m[sel] = m[sel], m[r]
        inv = pow(m[r][ccol], 11, 13)
        m[r] = [(x * inv) % 13 for x in m[r]]
        for rr in range(len(m)):
            if rr != r and m[rr][ccol] % 13:
                f = m[rr][ccol]
                m[rr] = [(m[rr][k] - f * m[r][k]) % 13 for k in range(ncol)]
        piv.append(ccol)
        r += 1
    free = [c for c in range(ncol) if c not in piv]
    basis = []
    for fc in free:
        v = [0] * ncol
        v[fc] = 1
        for ri, pc in enumerate(piv):
            v[pc] = (-m[ri][fc]) % 13
        basis.append(tuple(v))
    return basis, len(piv)


def gf13_rank(vecs):
    m = [list(v) for v in vecs]
    rank = 0
    for ccol in range(9):
        sel = None
        for rr in range(rank, len(m)):
            if m[rr][ccol] % 13:
                sel = rr
                break
        if sel is None:
            continue
        m[rank], m[sel] = m[sel], m[rank]
        inv = pow(m[rank][ccol], 11, 13)
        m[rank] = [(x * inv) % 13 for x in m[rank]]
        for rr in range(len(m)):
            if rr != rank and m[rr][ccol] % 13:
                f = m[rr][ccol]
                m[rr] = [(m[rr][k] - f * m[rank][k]) % 13 for k in range(9)]
        rank += 1
    return rank


# ==========================================================================
# MAIN
# ==========================================================================
def main():
    log("THM-2334 Section 7 decisive step: exact 169-twist variance test")
    log("script=04-computation/lrc14_169_twist_variance_opus_20260726.py")
    log("machine=opus  date=2026-07-26")
    log("")
    log("[1] canonical row/word stratum (THM-2349 / THM-2309 eq (25))")
    log(f"    w=(H,q1..q5,c1,c2,c3)={W}")
    log(f"    valuation profile nu13(c1,c2,c3)={PROFILE}  (strict, in THM-2349's 165-row universe)")
    log(f"    owner j=1 (c_1={W[OWNER]}), targets a=c_2={W[TA]}, b=c_3={W[TB]} (deepest comb leg)")
    log(f"    delayed clock k={K_CLOCK}, word dilation R=13^k={RDIL}=13^(lambda_j+1) (THM-2334 (3))")
    log(f"    T_DEN={T_DEN}")
    log(f"    NN=169*T_DEN={NN}  (cyclotomic order; primes {NN_PRIMES})")
    log("")

    # -- embeddings ---------------------------------------------------------
    c1c, P1, H1 = find_embedding()
    c2c, P2, H2 = find_embedding(skip_c=(c1c,))
    for (pp, hh) in ((P1, H1), (P2, H2)):
        assert is_prime(pp) and pow(hh, NN, pp) == 1
        assert all(pow(hh, NN // q, pp) != 1 for q in NN_PRIMES)
    log("[2] exact ring homomorphisms Z[zeta_NN] -> F_p  (zeta_NN -> h, order NN)")
    log(f"    p1={P1}  (={c1c}*NN+1)   h1={H1}")
    log(f"    p2={P2}  (={c2c}*NN+1)   h2={H2}")
    log("    order checks: PASS")
    log("")

    # -- sets ---------------------------------------------------------------
    log("[3] exact interval sets (integer endpoints over T_DEN)")
    for i in range(9):
        assert sum(b - a for a, b in in_comb(i, 0)) == T_DEN // 7
    E0 = build_set(PAT_E, ZELL)
    lenE0 = check_intervals(E0)
    mu_E = Fraction(lenE0, T_DEN)
    log(f"    E_1: intervals={len(E0)}  measure={mu_E}  ~{float(mu_E):.6f}")
    assert mu_E > 0, "E_1 must have positive measure"

    QSETS = {}
    for name, pat in (("{a}", PAT_QA), ("{b}", PAT_QB), ("{a,b}", PAT_QAB)):
        qiv = build_set(pat, ZELL)
        qlen = check_intervals(qiv)
        QSETS[name] = qiv
        log(f"    Q_1,{name}: intervals={len(qiv)}  measure={Fraction(qlen, T_DEN)}"
            f"  ~{qlen / T_DEN:.6f}")
    log(f"    {elapsed()} sets built")
    log("")

    # -- word stratum choice ------------------------------------------------
    log("[4] delayed word stratum E_1 n T^-2 Q_1,sigma (THM-2349 (31)-(32))")
    sigma = None
    Q_iv = None
    for name in ("{a}", "{b}", "{a,b}"):
        qiv = QSETS[name]
        qstarts = [a for a, _ in qiv]
        _, _, ov = x_sweep(E0, qiv, qstarts, 1, [], [], want_values=False)
        mu = Fraction(ov, NN)
        log(f"    measure(E_1 n T^-2 Q_1,{name}) = {mu}  ~{float(mu):.8f}")
        if sigma is None and mu > 0:
            sigma, Q_iv, mu_strat = name, qiv, mu
    assert sigma is not None, "no positive word stratum at k=2"
    Q_starts = [a for a, _ in Q_iv]
    log(f"    chosen canonical word sigma={sigma} (first positive in THM-2305 order)")
    log(f"    stratum measure (exact rational) = {mu_strat} > 0")
    log(f"    {elapsed()}")
    log("")

    # -- frequency triangle (X, Y=X+m*c3, m) --------------------------------
    log("[5] marked triangle X, Y=X+m*c_3, gcd(m,91)=1 (THM-2349 (5)-(7), THM-2327 (32)-(33))")
    mods2 = [(P1, H1), (P2, H2)]
    XCANDS = [13 * (1 + 13 * j) for j in (0, 1, -1, 2, -2, 3, -3, 4, -4)]
    MCANDS = [m for m in range(1, 25) if gcd(m, 91) == 1]
    chosen = None
    for m in MCANDS[:6]:
        for X in XCANDS:
            Y = X + m * W[TB]
            eX_mod, eX_c = endpoint_sum(E0, X, mods2)
            eY_mod, eY_c = endpoint_sum(E0, Y, mods2)
            tabs = make_tabs(Q_iv, X, mods2)
            fX_mod, fX_c, _ = x_sweep(E0, Q_iv, Q_starts, X, mods2, tabs)
            # require the p1 embedding (used for the whole bank) to certify all
            # three triangle nonvanishings directly
            ok = bool(fX_mod[0] and eX_mod[0] and eY_mod[0])
            if ok:
                chosen = (X, m, Y, fX_mod, fX_c, eX_mod, eX_c, eY_mod, eY_c)
                break
        if chosen:
            break
    assert chosen is not None, "no certified nonzero triangle found in search window"
    X, m, Y = chosen[0], chosen[1], chosen[2]
    fX_mod, fX_c, eX_mod, eX_c, eY_mod, eY_c = chosen[3:]
    kappa = (X // 13) % 13
    assert nu13(X) == 1 and nu13(Y) == 1 and (Y // 13 - X // 13) % 13 == 0
    assert gcd(m, 91) == 1
    log(f"    X={X}  m={m}  Y={Y}  kappa={kappa}  nu13(X)=nu13(Y)=1  gcd(m,91)=1")
    log(f"    certified nonzero (image under p1,p2 / numeric |hat|):")
    log(f"      (2 pi i X)*f_hat(X)  images={fX_mod}  |f_hat(X)| ~{abs(fX_c / (TWO_PI * X)):.3e}")
    log(f"      (2 pi i X)*e_hat(X)  images={eX_mod}  |e_hat(X)| ~{abs(eX_c / (TWO_PI * X)):.3e}")
    log(f"      (2 pi i Y)*e_hat(Y)  images={eY_mod}  |e_hat(Y)| ~{abs(eY_c / (TWO_PI * Y)):.3e}")
    dm_num = math.sin(math.pi * m / 7.0) / (math.pi * m)
    log(f"    deep leg delta_hat(m)=sin(pi m/7)/(pi m) ~{dm_num:.6f} != 0 (7 does not divide m)")
    log(f"    {elapsed()}")
    log("")

    # -- the 169 twists G^ --------------------------------------------------
    log("[6] target twist group G^ = L_13^perp/<w> (THM-2309 star+graft, u0=q5, grafts q1->a, q2->b)")
    u0 = 5
    Ppiv = [OWNER, 0, 1, 2, 3, 4]
    rows = []
    for k in Ppiv:
        r = [0] * 9
        r[u0] = (r[u0] + W[k]) % 13
        r[k] = (r[k] - W[u0]) % 13
        rows.append(r)
    ga = rows[Ppiv.index(1)]
    ga[u0] = (ga[u0] + W[TA]) % 13
    ga[TA] = (ga[TA] - W[u0]) % 13
    gb = rows[Ppiv.index(2)]
    gb[u0] = (gb[u0] + W[TB]) % 13
    gb[TB] = (gb[TB] - W[u0]) % 13
    wmod = tuple(v % 13 for v in W)
    for r in rows:
        assert sum(a * b for a, b in zip(r, wmod)) % 13 == 0
    assert gf13_rank(rows) == 6
    null_basis, rk = gf13_nullspace(rows)
    assert len(null_basis) == 3
    assert gf13_rank(list(null_basis) + [wmod]) == 3  # w in L^perp
    vsel = []
    cur = [wmod]
    for nb in null_basis:
        if gf13_rank(cur + [nb]) > gf13_rank(cur):
            cur.append(nb)
            vsel.append(nb)
    assert len(vsel) == 2
    v1, v2 = vsel
    log(f"    L rank=6; L^perp basis contains w; coset generators")
    log(f"    v1={v1}")
    log(f"    v2={v2}")
    reps = {}
    for al in range(13):
        for be in range(13):
            reps[(al, be)] = tuple((al * v1[i] + be * v2[i]) % 13 for i in range(9))
    assert len(set(reps.values())) == 169
    log("    169 coset representatives enumerated")
    log("")

    # -- twist bank ---------------------------------------------------------
    log("[7] twist bank H(ell): per-twist exact images gamma(ell) mod p1 + numeric H(ell)")
    log(f"    gamma(ell) = e_13(m*ell_8) * (2 pi i X)hat(E_Q^ell)(X) * conj((2 pi i Y)hat(E^ell)(Y))")
    log(f"    H(ell) = C0*gamma(ell),  C0 = delta_hat({m})/(4 pi^2 X Y)  (real, nonzero)")
    mods1 = [(P1, H1)]
    tabs1 = make_tabs(Q_iv, X, mods1)
    e13_ph_1 = [pow(H1, (m * l8 * (NN // 13)) % NN, P1) for l8 in range(13)]
    gam_p1 = {}
    H_num = {}
    ecount = 0
    for al in range(13):
        for be in range(13):
            ell = reps[(al, be)]
            E_iv = build_set(PAT_E, ell)
            ecount += len(E_iv)
            AX_mod, AX_c, _ = x_sweep(E_iv, Q_iv, Q_starts, X, mods1, tabs1)
            BYc_mod, BYc_c = endpoint_sum(E_iv, -Y, mods1)   # conj side: e(+Y t)
            g = e13_ph_1[ell[8] % 13] * AX_mod[0] % P1 * BYc_mod[0] % P1
            gam_p1[(al, be)] = g
            ph_c = cmath.exp(2j * math.pi * m * ell[8] / 13.0)
            H_num[(al, be)] = (dm_num * ph_c * (AX_c / (2j * math.pi * X))
                              * (BYc_c / (-2j * math.pi * Y)))
        log(f"    {elapsed()} alpha={al} done (avg |E^ell| intervals so far: {ecount // (13 * (al + 1))})")

    # gauge (representative-independence) check, THM-2334 (6)/(29)
    log("")
    log("[8] gauge check: H(ell + w) = H(ell) exactly (THM-2334 representative independence)")
    for base_ab in ((0, 0), (1, 0), (0, 1), (3, 7)):
        ell = reps[base_ab]
        ellw = tuple((ell[i] + wmod[i]) % 13 for i in range(9))
        E_iv = build_set(PAT_E, ellw)
        AX_mod, AX_c, _ = x_sweep(E_iv, Q_iv, Q_starts, X, mods1, tabs1)
        BYc_mod, BYc_c = endpoint_sum(E_iv, -Y, mods1)
        g_w = e13_ph_1[ellw[8] % 13] * AX_mod[0] % P1 * BYc_mod[0] % P1
        same = (g_w - gam_p1[base_ab]) % P1 == 0
        ph_c = cmath.exp(2j * math.pi * m * ellw[8] / 13.0)
        h_w = (dm_num * ph_c * (AX_c / (2j * math.pi * X))
               * (BYc_c / (-2j * math.pi * Y)))
        num_err = abs(h_w - H_num[base_ab]) / max(1e-300, abs(H_num[base_ab]))
        log(f"    ell(alpha,beta)={base_ab}: mod-p1 equal={same}  numeric rel err={num_err:.2e}")
        assert same, "gauge invariance failed -- convention or implementation error"
    log("    gauge checks: PASS")
    log("")

    # -- decision -----------------------------------------------------------
    log("[9] DECISION: exact nonconstancy test of the 169-twist bank (THM-2334 (42)-(43))")
    g0 = gam_p1[(0, 0)]
    witnesses = [(al, be) for (al, be), g in gam_p1.items()
                 if (g - g0) % P1 != 0]
    log(f"    twists with gamma(ell) != gamma(0) certified mod p1: {len(witnesses)} / 168 nonzero twists")
    if witnesses:
        wal, wbe = witnesses[0]
        log(f"    first witness: (alpha,beta)=({wal},{wbe})  ell={reps[(wal, wbe)]}")
        log(f"      gamma_p1(witness)={gam_p1[(wal, wbe)]}   gamma_p1(0)={g0}")
        # independent confirmation under second embedding
        mods2b = [(P2, H2)]
        tabs2b = make_tabs(Q_iv, X, mods2b)
        e13_ph_2 = [pow(H2, (m * l8 * (NN // 13)) % NN, P2) for l8 in range(13)]
        conf = {}
        for ab in ((0, 0), (wal, wbe)):
            E_iv = build_set(PAT_E, reps[ab])
            AX_mod, _, _ = x_sweep(E_iv, Q_iv, Q_starts, X, mods2b, tabs2b)
            BYc_mod, _ = endpoint_sum(E_iv, -Y, mods2b)
            conf[ab] = e13_ph_2[reps[ab][8] % 13] * AX_mod[0] % P2 * BYc_mod[0] % P2
        agree = (conf[(wal, wbe)] - conf[(0, 0)]) % P2 != 0
        log(f"      second embedding p2 also distinguishes witness: {agree}")
        log("")
        log("    VERDICT: NONCONSTANT twist bank -- CERTIFIED (exact ring-homomorphism")
        log("    certificate: gamma(ell)-gamma(0) has nonzero image in F_p1, hence is a")
        log("    nonzero element of Z[zeta_NN], hence H(ell) != H(0) as complex numbers).")
        log("    By THM-2334 (42)-(43), for this row/word stratum and triangle the")
        log("    nonzero-target mass is strictly positive:")
        log("        sum_{q != 0} |A(q)|^2 = (1/169) sum_ell |H(ell)-Hbar|^2 > 0.")
    else:
        log("    all 169 gamma(ell) images coincide mod p1 -- no positivity certificate;")
        log("    running full bank under second embedding p2 for corroboration...")
        mods2b = [(P2, H2)]
        tabs2b = make_tabs(Q_iv, X, mods2b)
        e13_ph_2 = [pow(H2, (m * l8 * (NN // 13)) % NN, P2) for l8 in range(13)]
        gam_p2 = {}
        for al in range(13):
            for be in range(13):
                E_iv = build_set(PAT_E, reps[(al, be)])
                AX_mod, _, _ = x_sweep(E_iv, Q_iv, Q_starts, X, mods2b, tabs2b)
                BYc_mod, _ = endpoint_sum(E_iv, -Y, mods2b)
                gam_p2[(al, be)] = (e13_ph_2[reps[(al, be)][8] % 13]
                                    * AX_mod[0] % P2 * BYc_mod[0] % P2)
        w2 = [ab for ab, g in gam_p2.items() if (g - gam_p2[(0, 0)]) % P2 != 0]
        log(f"    twists differing under p2: {len(w2)}")
        log("    VERDICT: bank constant under two independent exact homomorphisms;")
        log("    positivity NOT certified (consistent with an all-cancelling boundary,")
        log("    cf. THM-2333/THM-2334 Sec. 10 hostiles). Constancy itself not proved.")
    log("")

    # -- exact inverse DFT: ALL 169 target aggregates (klein MSG-2171/2172) --
    log("[9b] EXACT inverse DFT of the gamma bank (klein strengthening request):")
    log("     Anum(q) := sum_ell gamma(ell) z13^-(ell[TA] qa + ell[TB] qb) mod p1,")
    log("     z13 = H1^(NN/13) (exact order-13 root in F_p1);  Anum(q) = 169 A(q)/C0.")
    log("     A nonzero F_p1 image proves the cyclotomic target-aggregate numerator")
    log("     is a nonzero algebraic number, hence A(q) != 0 as a complex number.")
    z13 = pow(H1, NN // 13, P1)
    z13p = [pow(z13, k, P1) for k in range(13)]
    # pairing well-definedness: ell -> (ell[TA], ell[TB]) mod 13 must biject
    pair_map = {(al, be): (reps[(al, be)][TA] % 13, reps[(al, be)][TB] % 13)
                for al in range(13) for be in range(13)}
    assert len(set(pair_map.values())) == 169, "TA/TB pairing degenerate on coset reps"
    log("     pairing ell -> (ell[TA], ell[TB]) is a bijection onto F_13^2: PASS")
    Aq_exact = {}
    for qa in range(13):
        for qb in range(13):
            s = 0
            for ab, g in gam_p1.items():
                ta, tb = pair_map[ab]
                s = (s + g * z13p[(-(ta * qa + tb * qb)) % 13]) % P1
            Aq_exact[(qa, qb)] = s
    nz = [q for q, s in Aq_exact.items() if s % P1 != 0]
    nz_targets = [q for q in nz if q != (0, 0)]
    log(f"     nonzero A(q) images mod p1: {len(nz)} / 169  "
        f"(nonzero targets q != 0: {len(nz_targets)} / 168)")
    # control A: sum_q Anum(q) = 169 gamma(0)  (forward DFT at ell = 0)
    ctrlA = (sum(Aq_exact.values()) - 169 * gam_p1[(0, 0)]) % P1 == 0
    log(f"     control A (sum_q A = H(0)):        {'PASS' if ctrlA else 'FAIL'}")
    assert ctrlA
    # control B: sign-reversal (positive-sign pairing = A(-q)): same census
    nz_rev = 0
    for qa in range(13):
        for qb in range(13):
            s = 0
            for ab, g in gam_p1.items():
                ta, tb = pair_map[ab]
                s = (s + g * z13p[(ta * qa + tb * qb) % 13]) % P1
            if s % P1 != 0:
                nz_rev += 1
    log(f"     control B (sign-reversal census):  {'PASS' if nz_rev == len(nz) else 'FAIL'}"
        f"  ({nz_rev} nonzero)")
    assert nz_rev == len(nz)
    # control C: forward reconstruction 169 gamma(ell) = sum_q Anum(q) z13^{+<ell,q>}
    for ab in ((0, 0), (2, 5)):
        ta, tb = pair_map[ab]
        s = 0
        for (qa, qb), v in Aq_exact.items():
            s = (s + v * z13p[(ta * qa + tb * qb) % 13]) % P1
        ok = (s - 169 * gam_p1[ab]) % P1 == 0
        log(f"     control C (forward at ell{ab}):  {'PASS' if ok else 'FAIL'}")
        assert ok
    if len(nz_targets) == 168:
        log("     VERDICT: FULL TARGET-PLANE SUPPORT -- every one of the 168 nonzero")
        log("     target aggregates A(q) is a NONZERO complex number (exact certificate).")
        log("     This is strictly stronger than nonconstancy of the bank.  Scope is")
        log("     unchanged: typed row (not a scalar cover), unrestricted A(q) (NOT the")
        log("     all-91-unit projector B(q) of THM-2334 (49)); no row is removed.")
    log("")

    # -- descriptive numerics ----------------------------------------------
    log("[10] descriptive numerics (floats; NOT part of the decision)")
    Hbar = sum(H_num.values()) / 169.0
    var = sum(abs(v - Hbar) ** 2 for v in H_num.values()) / 169.0
    log(f"    H(0) ~ {fmtc(H_num[(0, 0)])}")
    log(f"    Hbar ~ {fmtc(Hbar)}")
    log(f"    variance (1/169) sum |H-Hbar|^2 ~ {var:.6e}")
    Aq = {}
    for qa in range(13):
        for qb in range(13):
            s = 0.0j
            for al in range(13):
                for be in range(13):
                    ell = reps[(al, be)]
                    ph = cmath.exp(-2j * math.pi * (ell[TA] * qa + ell[TB] * qb) / 13.0)
                    s += ph * H_num[(al, be)]
            Aq[(qa, qb)] = s / 169.0
    sumA = sum(Aq.values())
    log(f"    sum_q A(q) ~ {fmtc(sumA)}   (THM-2334 (40): should equal H(0); "
        f"err={abs(sumA - H_num[(0, 0)]):.2e})")
    off_mass = sum(abs(v) ** 2 for q, v in Aq.items() if q != (0, 0))
    log(f"    sum_(q!=0) |A(q)|^2 ~ {off_mass:.6e}   (Parseval vs variance err="
        f"{abs(off_mass - var):.2e})")
    top = sorted(((abs(v), q) for q, v in Aq.items() if q != (0, 0)), reverse=True)[:5]
    log(f"    A(0,0) ~ {fmtc(Aq[(0, 0)])}")
    for mag, q in top:
        log(f"    |A{q}| ~ {mag:.6e}")
    log("")

    log("[11] scope and mismatch notes")
    log("    * THM-2334 (31)-(32) 'products' enter through their Fourier coefficients")
    log("      at the marked frequencies X and Y (eqs (28)-(29),(41)); the variance (42)")
    log("      is therefore NOT a rational number: it equals")
    log("      [delta_hat(m)/(4 pi^2 X Y)]^2 * (cyclotomic real >= 0).  The task-summary")
    log("      phrase 'exact positive rational variance' cannot be met literally; the")
    log("      file-faithful exact decision is the certified nonconstancy above.")
    log("    * The word translations R*ell_i/13 are integers (THM-2334 (44)), so the")
    log("      transported word is target-neutral, as implemented.")
    log("    * The row is the canon's typed row THM-2309 (25) (guard=1 odd, profile")
    log("      (1,3,5)); it is typed but not asserted to be a scalar cover.  The")
    log("      computation decides the variance object for THIS canonical typed row's")
    log("      marked current; it does not by itself close any covering row or LRC(14).")
    log("")
    log(f"{elapsed()} all checks passed")


if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    main()
