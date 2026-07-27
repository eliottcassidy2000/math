#!/usr/bin/env python3
"""STAGE-2 THETA-SLAVED CONTRACTION at the newly opened r = 5 window
(canonical typed row, sigma = {b}, packet clock K = 2).

TASK (2026-07-28).  THM-2575/THM-2577 proved that on the canonical typed row
THM-2309 (25)

    w = (H,q1..q5,c1,c2,c3) = (1,14,27,40,53,66,13,2197,742586),

owner j = 1, word sigma = {b}, EVERY positive-return packet clock has
first-collision index r = 5, so the collision depth is d = 13^5 and the
THM-2471 (44)-(47) deep-root sidecar sits in its unique-affine-invariant
window:  h = c_3/d = 2, and the ONLY lawful (u,t)-coupling of the deep probe
to the collision stalk is through

    theta = t - 2u   (mod 13)          (THM-2471 (46)-(47), current branch).

This script runs the Stage-2 theta-slaved toothpick contraction there:
build ONE joint exact Boolean-fibre-product table on ONE common base,
gate it against the THM-2471 colour laws (J(0) = I_5 anchor), centre it
(lawful common-base ANOVA), contract it in the THM-2512 (12)-(14) toothpick
shape with the F_13 output index v SLAVED through theta, and decide the
resulting exact cyclotomic number S in Q(zeta_91).

OBJECT (all conventions logged; files override summaries):
    e = 1_{E_1},   f = 1_{Q_{1,{b}}} P^2 1_{E_1}      (THM-2471 (23), K = 2),
    d = 13^5,      U = P_d f,   V = P_d e             (THM-2471 (25), r = 5),
    A(y,u) = U((y+u)/13),  F(y,v) = V((y+v)/13)       (THM-2471 (4)-(5)).
    Stalk sheets (THM-2471 (28)):  w_u = (y+u)/13,  X_{u,a} = (w_u+a)/d,
    a in Z/dZ (current branch, carrying 1_Q and the packet ancestry).

RESPONSE-SIDE REFINEMENT (THM-2449 (14)-(16)-shaped Boolean factors on
LINKED NODES OF ONE FINITE ANCESTRY STALK -- inside one y-integral, never
as separately integrated controls; MISTAKE-281/293 discipline):
  * owner-clock factor d_{1,ell}(x) = d(13 x - ell/7) at the chart point
    w_u:   13 w_u = y + u  ==>  d(y - ell/7).  ROOT-BLIND at the stalk
    (u drops out as an integer): the owner clock reads the base exactly.
    cell_ell := {y : ||y - ell/7|| < 1/14},  7 cells tiling T a.e.
  * deep probe Delta_t(x) = d(c_3 x - t/13) on the current deep sheets:
    c_3 X_{u,a} = 2(w_u + a) == 2(y+u)/13  (mod 1)  -- SHEET-INDEPENDENT
    (THM-2471 (45) with d_0 = 1, h = 2).  Hence
       Delta_t(X_{u,a}) = 1_{||(2y - theta)/13 mod 1|| < 1/14},
       theta = (t - 2u) mod 13,
    the sidecar's unique affine invariant, realized.  As base sets:
       DEEP_0 = [0, 13/28),  DEEP_1 = (1/28, 27/28),  DEEP_2 = (15/28, 1),
       DEEP_theta = EMPTY for theta = 3..12   (2y/13 covers only [0,2/13)),
    and the exact root-count law at the stalk (THM-2403 (26) pushed through
    the sidecar):  sum_theta DEEP_theta = 2 - d(2y mod 1)  pointwise a.e.
    Per fixed root u, exactly the deep labels t == 2u, 2u+1, 2u+2 (mod 13)
    can fire: the slaving transports the deep support with the root.
  * word factor Q^eps_ell(R x) at response exponent R = 13^k
    (THM-2449 (15); the THM-2409 skew-diagonal is mandatory).  CLOCK
    ARITHMETIC for sigma={b}:
    the THM-2449 (25)-(31) eventual clock law uses classes
    k == K13 (mod ord_{D0}(13)) beyond an unspecified k_0, with
    T_DEN = 13^K13 D0 and K13 = 6.  The class arithmetic is word-independent,
    but does NOT prove 6 >= k_0.  Here k = 6 is instead evaluated directly
    and exactly.  It is the smallest sheet-return exponent in the class,
    eps = 13^6 mod 7 = +1.  At k = 6 the word factor is SHEET-EXACT on the
    current branch:  R X_{u,a} = 13(w_u + a) = y + u + 13a == y (mod 1),
    for every sheet a and root u; k = 6 is the unique clock with this
    property (k = 5 gives w_u, k = 7 gives 13y).  Hence
       Q^{+1}_ell(R X_{u,a}) = QW_ell(y) := T_b(y) * (1 - d(13y - ell/7)),
    T_b = the sigma={b} terminal word with its source-safe factor deleted
    (THM-2449 (15)):  T_b = A_0 n D_{c3} \\ D_{c2}.

THE ONE TABLE (single integral over the common base y, per entry):
    N(u, q, ell, theta)
      = int_T A(y,u) F(y,q) cell_ell(y) DEEP_theta(y) QW_ell(y) dy,
    u = current root, q = source root (s = u - q mod 13 the colour shift),
    ell in F_7, theta in F_13 (the slaved deep invariant; columns 3..12
    are exactly zero because their windows are empty).  Everything is
    paired at the same base point INSIDE the integral BEFORE any
    marginalization or DFT (MISTAKE-281).  Auxiliary members of the same
    integral family (for gates): P(u,q) (no refiner), cell-only, deep-only,
    d(2y)-refined, and word+cell (GW) / word+cell+d2 (GCD2) tables.

GATES:
  (g1) root-marginal colours: C(s) = sum_q P(q+s,q), J(k) = (1/169) *
       sum_s C(s) zeta^{-ks};  J(0) must equal the THM-2575/2577 anchor
       I_5 = 48602521488933856/337437093630814766589;  sum_s C(s) = 169 I_5;
       J(k) != 0 for all k != 0;  sum_k J(k) = 0; Hermitian symmetry.
  (g2) zero-sum laws: UV = 0 a.e. => every diagonal pair profile is
       IDENTICALLY zero (u = q), so C(0) = 0 and every refined diagonal
       entry vanishes.
  (g3) partition/law gates per pair and on the aggregate: sum_ell cell = 1
       a.e.;  sum_theta DEEP_theta = 2 - d(2y);  the refined versions
       sum_theta N = 2*GW - GCD2 entrywise.
  (g4) route-2 aggregate: the whole aggregate table is recomputed from the
       folded service profiles a(y) = sum_u A(y,u), b(y) = sum_q F(y,q)
       (different pairing geometry) and must agree entrywise.
  (g5) anchors: mu(E_1), |E_1|, mu(Q_{1,{b}}), |Q_{1,{b}}|, int f = rho^(2)({b})
       = 35505957232/16132831966251 (census), int U = int f, int V = mu(E).

CENTRING + CONTRACTION (THM-2512 (1)-(5), (12)-(14)):
    Aresp[ell][theta] = sum_{u,q} N(u,q,ell,theta)   (response table,
      = int_T a b cell_ell DEEP_theta QW_ell dy  by Fubini INSIDE the one
      integral; verified by route (g4));
    I_(ell,theta) = Aresp - row - col + grand   (all centring terms are
      marginals of the ONE table);  d(v,ell) := I_(ell,v)  (transpose);
    fibrewise version: the same centring on each (u, s) ancestry fibre,
    with sum_fibres = aggregate checked exactly (linearity).
    TOOTHPICK (output index v = theta, i.e. v IS the slaved deep invariant):
      R_{tau,a0,c}(v)   = sum_ell d(v - tau*rep(a0*ell + c), ell),
      Theta_{tau,a0,c}(alpha) = sum_v R(v) zeta13^{-alpha v},
      Psi_{tau,a0}(alpha,beta) = sum_c Theta_c(alpha) zeta7^{-beta c},
    rep: F_7 -> {0..6} c Z -> F_13 (THM-2508 convention).
    DECISIVE NUMBER:  S := Psi_{1,1}(1,1) at the fixed primitive quadruple
    (tau,a0,alpha,beta) = (1,1,1,1) -- chosen from the chart conventions
    alone (first primitive quadruple in the fixed order), independent of
    every ancestry point.  S in Q(zeta_91); decided by exact reduction of
    its integer numerator polynomial mod Phi_91 (zeta13 = x^7, zeta7 = x^13,
    x = zeta_91; the exponent map (i,j) -> 7i+13j mod 91 is the CRT
    bijection), cross-checked at the five split primes
    547, 911, 1093, 2003, 2549 (all == 1 mod 91), and cross-checked against
    the THM-2512 (15) factorization Psi = K_{alpha tau, beta} *
    dtilde(alpha, -beta a0).  By THM-2512 (19) (rational table, primitive
    quadruple):  S != 0  iff  the interaction is nonzero  iff  ALL 5184
    primitive cut coefficients fire.
CONTROLS:
    HOSTILE H1 (must vanish): beta = 0 -- Psi_{1,1}(1,0) = 0 exactly by the
      row-zero law (THM-2512 (4)).
    HOSTILE H2 (must vanish): u-independent slaving -- replace the slaved
      column family by its theta-average (the deep label carries no
      root-coupled information); the centred interaction is then 0
      identically and the same pipeline must output S = 0 exactly.
    POSITIVE: J(0) = I_5 anchor; sum_s C(s) = 169 I_5; int f = census rho;
      the (15) factorization identity; route-2 equality.

VERIFICATION MATRIX:
  (a) toy brute-force self-test of the ENTIRE pipeline on a hostile little
      model (grid 1260, base p = 3, packet clock 9, depth 27, deep speed
      54 = 2*27 so h' = 2 and theta' = t - 2u mod 3, cells mod 5, word set):
      definitional midpoint evaluation of every table entry -- INCLUDING
      literal stalk-sheet evaluation of the deep/word/cell factors at
      X = (w_u + a)/27 with sheet-independence checked -- against the
      production engine functions;
  (b) production slaving spot-check: exact Fraction evaluation of
      Delta_t(X_{u,a}) and window memberships at sampled stalk points vs
      the DEEP/CELL/QW interval sets (boundary hits detected and excluded);
  (c) the gates (g1)-(g5) above;
  (d) split primes + factorization identity + THM-2512 (19) consistency.

SCOPE (MISTAKE-281/283/286/287 discipline, hard):
  * TYPED row (THM-2309 (25)); NOT an asserted scalar cover;
  * the theta = t - 2u map is THM-2471's deep-root-sidecar candidate
    realization at the r = 5 window; it is NOT "THM-2512's bridge"
    (MISTAKE-286): THM-2512 Sec 5 is generic and contains no r, theta, or
    clock; a nonzero S here is a positive instance of THE SPECIFIC
    candidate, not a bridge theorem;
  * no physical current, no row exclusion, no LRC(14)/JC(2) progress is
    claimed; every pairing sits inside one integral on one common base
    BEFORE marginalization/DFT (MISTAKE-281); no square factor is
    discarded anywhere (MISTAKE-287: all decisions are exact integer
    zero-tests, no localization step occurs);
  * all decisions are integer/Fraction/Z[zeta]-exact; floats only in ~.

Script: 04-computation/lrc14_stage2_theta_contraction_opus_20260728.py
Output: 05-knowledge/results/lrc14_stage2_theta_contraction_opus_20260728.out
"""

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


def require(cond, msg):
    if not cond:
        raise SystemExit(f"REQUIRE FAILED: {msg}")


# --------------------------------------------------------------------------
# Row data and grid (identical to the audited census/bridge scripts)
# --------------------------------------------------------------------------
W = (1, 14, 27, 40, 53, 66, 13, 13**3, 2 * 13**5)
GUARD, OWNER, TA, TB = 0, 6, 7, 8
UNIT_IDX = (1, 2, 3, 4, 5)
C1, C2, C3 = W[OWNER], W[TA], W[TB]


def nu13(n):
    v = 0
    while n % 13 == 0:
        n //= 13
        v += 1
    return v


require((nu13(C1), nu13(C2), nu13(C3)) == (1, 3, 5), "valuation profile")
LAMBDA = nu13(C1)
KCLK = LAMBDA + 1
RPKT = 13**KCLK
require(RPKT == 169, "packet ancestry clock R = 13^(lambda+1) = 169")

LCM_W = 1
for v in W:
    LCM_W = LCM_W * v // gcd(LCM_W, v)
T_DEN = 182 * LCM_W
require(T_DEN == 297836897838480, "T_DEN anchor")
require(nu13(T_DEN) == 6, "13-adic valuation of T_DEN must be 6")
require(T_DEN % 28 == 0 and T_DEN % 14 == 0, "grid resolves /14 and /28")

RCOLL = 5                            # r = 5 (THM-2575/2577, sigma={b}, K=2)
DCOLL = 13**RCOLL                    # d = c_1 * 13^(r-1) = 13^5
HDEEP = C3 // DCOLL
require(gcd(C3, DCOLL) == DCOLL and HDEEP == 2 and HDEEP % 13 != 0,
        "r=5 window arithmetic: d_0 = 1, h = C/d = 2")
I5 = Fraction(48602521488933856, 337437093630814766589)  # THM-2575/2577
RHO_B = Fraction(35505957232, 16132831966251)            # census rho^(2)({b})

K13 = 6
D0 = T_DEN // 13**K13
require(D0 % 13 != 0, "D0 must be a 13-unit")
x0, ORD = 13 % D0, 1
while x0 != 1:
    x0 = x0 * 13 % D0
    ORD += 1
    require(ORD < 10**7, "order search overflow")
KRESP = 6                            # smallest exact sheet-return exponent
EPS = pow(13, KRESP, 7)
require(EPS == 1, "eps = 13^6 mod 7 = +1")


# --------------------------------------------------------------------------
# Interval machinery (verbatim from the audited census script)
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


def check_intervals(iv, top):
    last = -1
    for a, b in iv:
        require(0 <= a < b <= top and a >= last, "interval list corrupt")
        last = b
    return sum(b - a for a, b in iv)


def intersect_lists(x, y):
    """Intersection of two sorted disjoint half-open interval lists."""
    out = []
    i = j = 0
    nx, ny = len(x), len(y)
    while i < nx and j < ny:
        a = x[i][0] if x[i][0] > y[j][0] else y[j][0]
        b = x[i][1] if x[i][1] < y[j][1] else y[j][1]
        if a < b:
            out.append((a, b))
        if x[i][1] < y[j][1]:
            i += 1
        else:
            j += 1
    return out


def merge_touch(iv):
    """Sorted merge of an interval collection, gluing touching endpoints."""
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1]:
            if b > out[-1][1]:
                out[-1] = (out[-1][0], b)
        else:
            out.append((a, b))
    return out


def measure(iv):
    return sum(b - a for a, b in iv)


PAT_E = {0: "gout", 1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
         6: "in", 7: "out", 8: "out"}
PAT_QB = {0: "gout", 1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
          6: "out", 7: "out", 8: "in"}
PAT_TB = {0: "gout", 1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
          7: "out", 8: "in"}          # word with source-safe factor DELETED
ZELL = (0,) * 9

# --------------------------------------------------------------------------
# Exact profile machinery (verbatim from the audited base-only bridge script)
# --------------------------------------------------------------------------
def weighted_fold(pieces, mult, grid):
    """Numerator profile of mult * P_mult h for step h given as integer-
    weighted pieces (a, b, w) on [0, grid).  Mass conservation asserted."""
    const = 0
    events = []
    ap = events.append
    tot_in = 0
    for a, b, w in pieces:
        if w == 0:
            continue
        tot_in += w * (b - a)
        L = mult * (b - a)
        qq, rem = divmod(L, grid)
        const += qq * w
        if rem:
            st = (mult * a) % grid
            en = st + rem
            if en <= grid:
                ap((st, w))
                ap((en, -w))
            else:
                ap((st, w))
                ap((grid, -w))
                ap((0, w))
                ap((en - grid, -w))
    events.sort()
    starts = [0]
    vals = [const]
    cur = const
    for pos, dv in events:
        if pos == starts[-1]:
            cur += dv
            vals[-1] = cur
        else:
            cur += dv
            starts.append(pos)
            vals.append(cur)
    require(all(v >= 0 for v in vals), "negative weighted fold value")
    require(profile_mass(starts, vals, grid) == mult * tot_in,
            "weighted fold mass not conserved")
    return starts, vals


def profile_mass(starts, vals, grid):
    tot = 0
    n = len(starts)
    for i in range(n):
        nxt = starts[i + 1] if i + 1 < n else grid
        tot += vals[i] * (nxt - starts[i])
    return tot


def extract_window(starts, vals, u, p, grid):
    """Profile of y -> h((y+u)/p) on [0,grid): chart-window pullback
    (THM-2471 (4)-(5)); window u of h is [u*grid/p, (u+1)*grid/p)."""
    T_p = grid // p
    lo = u * T_p
    hi = lo + T_p
    n = len(starts)
    i = bisect_right(starts, lo) - 1
    ns = []
    nv = []
    while i < n and starts[i] < hi:
        a = starts[i] if starts[i] > lo else lo
        v = vals[i]
        if ns and nv[-1] == v:
            pass
        else:
            ns.append(p * a - u * grid)
            nv.append(v)
        i += 1
    require(ns and ns[0] == 0, "window profile must start at 0")
    return ns, nv


def product_cum(s1, v1, s2, v2, grid):
    """Merged product profile of two profiles + cumulative integral arrays."""
    n1 = len(s1)
    n2 = len(s2)
    i = j = 0
    ps = [0]
    pv = [v1[0] * v2[0]]
    while True:
        b1 = s1[i + 1] if i + 1 < n1 else grid
        b2 = s2[j + 1] if j + 1 < n2 else grid
        b = b1 if b1 < b2 else b2
        if b >= grid:
            break
        if b == b1:
            i += 1
        if b == b2:
            j += 1
        val = v1[i] * v2[j]
        if val != pv[-1]:
            ps.append(b)
            pv.append(val)
    pc = [0] * len(ps)
    for t in range(1, len(ps)):
        pc[t] = pc[t - 1] + pv[t - 1] * (ps[t] - ps[t - 1])
    ptot = pc[-1] + pv[-1] * (grid - ps[-1])
    return ps, pv, pc, ptot


def pquery(ps, pv, pc, P):
    i = bisect_right(ps, P) - 1
    return pc[i] + pv[i] * (P - ps[i])


def set_integral(prof, iv):
    """int of the (ps,pv,pc,ptot) profile over the interval list iv."""
    ps, pv, pc, _ = prof
    acc = 0
    for a, b in iv:
        acc += pquery(ps, pv, pc, b) - pquery(ps, pv, pc, a)
    return acc


def member(iv_starts, iv_ends, pos):
    """Half-open membership pos in U [a,b) for sorted disjoint intervals."""
    i = bisect_right(iv_starts, pos) - 1
    return i >= 0 and pos < iv_ends[i]


# --------------------------------------------------------------------------
# Exact Z[zeta_13] layer (verbatim bridge conventions)
# --------------------------------------------------------------------------
def canon(b):
    t = b[12]
    return tuple(x - t for x in b[:12])


def is_zero_b(b):
    return all(x == 0 for x in canon(b))


def badd(a, c):
    return [x + y for x, y in zip(a, c)]


def bconj(a):
    return [a[(13 - j) % 13] for j in range(13)]


def factor_str(n):
    out = []
    for pp in (2, 3, 5, 7, 11, 13, 53):
        e = 0
        while n % pp == 0:
            n //= pp
            e += 1
        if e:
            out.append(f"{pp}^{e}" if e > 1 else f"{pp}")
    if n > 1:
        out.append(f"[cofactor {n}]")
    return " * ".join(out) if out else "1"


# --------------------------------------------------------------------------
# Response-side refiner sets on the base (derivations in the docstring)
# --------------------------------------------------------------------------
def build_cells():
    """cell_ell = {y : ||y - ell/7|| < 1/14}, half-open realization."""
    cells = []
    T7 = T_DEN // 7
    T14 = T_DEN // 14
    for ell in range(7):
        a = (ell * T7 - T14) % T_DEN
        b = a + 2 * T14
        if b <= T_DEN:
            cells.append([(a, b)])
        else:
            cells.append([(0, b - T_DEN), (a, T_DEN)])
    return cells


def build_deep():
    """DEEP_theta = {y : ||(2y - theta)/13 mod 1|| < 1/14}, theta in F_13,
    half-open.  Nonempty exactly for theta in {0,1,2} (docstring)."""
    T28 = T_DEN // 28
    deep = [[] for _ in range(13)]
    deep[0] = [(0, 13 * T28)]
    deep[1] = [(T28, 27 * T28)]
    deep[2] = [(15 * T28, 28 * T28)]
    return deep


def build_d2():
    """{y : ||2y mod 1|| < 1/14}, half-open."""
    T28 = T_DEN // 28
    return [(0, T28), (13 * T28, 15 * T28), (27 * T28, 28 * T28)]


def check_deep_law(deep, d2):
    """Exact pointwise law sum_theta DEEP_theta = 2 - 1_{d2} on [0,T)."""
    bps = {0, T_DEN}
    for th in range(13):
        for a, b in deep[th]:
            bps.add(a)
            bps.add(b)
    for a, b in d2:
        bps.add(a)
        bps.add(b)
    bps = sorted(bps)
    d2s = [a for a, _ in d2]
    d2e = [b for _, b in d2]
    dst = [[a for a, _ in deep[th]] for th in range(13)]
    den = [[b for _, b in deep[th]] for th in range(13)]
    for i in range(len(bps) - 1):
        mid = (bps[i] + bps[i + 1]) // 2
        cnt = sum(1 for th in range(13) if member(dst[th], den[th], mid))
        want = 2 - (1 if member(d2s, d2e, mid) else 0)
        require(cnt == want, f"deep root-count law fails on piece {i}")


# --------------------------------------------------------------------------
# Toy brute-force self-test of the ENTIRE pipeline (route (a))
# --------------------------------------------------------------------------
def toy_test():
    """Hostile little model: base sets on grid 90 scaled to the engine grid
    GT = 1260; base p = 3, packet clock 9, depth 27, deep speed C2 = 2*27
    = 54 (h = 2, slaved theta = (t - 2u) mod 3), cells mod 5 of width 1/5,
    word set Wp.  Definitional midpoint evaluation (integer exact), with
    LITERAL stalk-sheet evaluation of every response factor at
    X = (w_u + a)/27 (sheet-independence checked), against the production
    engine functions."""
    GT = 1260
    sc = GT // 90
    Ep = [(0 * sc, 7 * sc), (17 * sc, 21 * sc), (34 * sc, 44 * sc),
          (60 * sc, 61 * sc), (85 * sc, 90 * sc)]
    Qp = [(4 * sc, 14 * sc), (40 * sc, 63 * sc), (77 * sc, 88 * sc)]
    Wp = [(2 * sc, 30 * sc), (50 * sc, 58 * sc), (66 * sc, 90 * sc)]
    p, Rp, dp, qm = 3, 9, 27, 5
    # engine sets: cells ||y-ell/5||<1/10 -> [252ell-126, 252ell+126) mod GT
    cells = []
    for ell in range(qm):
        a = (ell * 252 - 126) % GT
        b = a + 252
        cells.append([(a, b)] if b <= GT else [(0, b - GT), (a, GT)])
    # deep windows ||(2y-t)/3 mod 1|| < 1/14 -> theta in {0,1,2} only:
    deepw = [[] for _ in range(3)]
    deepw[0] = [(0, 135)]            # y in [0, 3/28)
    deepw[1] = [(495, 765)]          # y in (11/28, 17/28)
    deepw[2] = [(1125, 1260)]        # y in (25/28, 1)
    # ---- engine route (production functions on the toy grid)
    n_st, n_v = weighted_fold([(a, b, 1) for a, b in Ep], Rp, GT)
    fp = []
    np_ = len(n_st)
    for qa, qb in Qp:
        i = bisect_right(n_st, qa) - 1
        while True:
            pa = n_st[i]
            pb = n_st[i + 1] if i + 1 < np_ else GT
            a = qa if qa > pa else pa
            b = qb if qb < pb else pb
            if a < b and n_v[i]:
                fp.append((a, b, n_v[i]))
            if pb >= qb:
                break
            i += 1
    Ust, Uv = weighted_fold(fp, dp, GT)
    Vst, Vv = weighted_fold([(a, b, 1) for a, b in Ep], dp, GT)
    UWt = [extract_window(Ust, Uv, u, p, GT) for u in range(p)]
    VWt = [extract_window(Vst, Vv, u, p, GT) for u in range(p)]
    DEN = Rp * dp * dp * GT
    eng = [[[[0] * 3 for _ in range(qm)] for _ in range(p)] for _ in range(p)]
    for u in range(p):
        for q in range(p):
            prof = product_cum(UWt[u][0], UWt[u][1],
                               VWt[q][0], VWt[q][1], GT)
            for ell in range(qm):
                base = intersect_lists(Wp, cells[ell])
                for th in range(3):
                    G = intersect_lists(base, deepw[th])
                    eng[u][q][ell][th] = set_integral(prof, G)
    # ---- definitional route (integer arithmetic, midpoints y = m/2520)
    def ind_scaled(iv, num, den_to_gt):
        # membership of num/(GT*den_to_gt) in iv (on the GT grid)
        xx = num % (GT * den_to_gt)
        for a, b in iv:
            if a * den_to_gt <= xx < b * den_to_gt:
                return 1
        return 0

    fmemo = {}

    def fprime(n):                   # count over 9 at n/204120
        r = fmemo.get(n)
        if r is None:
            if not ind_scaled(Qp, n, 162):
                r = 0
            else:
                r = sum(ind_scaled(Ep, n + 204120 * c, 1458)
                        for c in range(9))
            fmemo[n] = r
        return r

    Umemo = {}

    def Uprime(m, u):                # count over 243 at (y+u)/3
        key = m + 2520 * u
        r = Umemo.get(key)
        if r is None:
            r = sum(fprime(m + 2520 * (u + 3 * a)) for a in range(dp))
            Umemo[key] = r
        return r

    Vmemo = {}

    def Vprime(m, u):                # count over 27 at (y+u)/3
        key = m + 2520 * u
        r = Vmemo.get(key)
        if r is None:
            r = sum(ind_scaled(Ep, m + 2520 * (u + 3 * a), 162)
                    for a in range(dp))
            Vmemo[key] = r
        return r

    ddef = [[[[0] * 3 for _ in range(qm)] for _ in range(p)]
            for _ in range(p)]
    dst = [[a for a, _ in deepw[t]] for t in range(3)]
    den_ = [[b for _, b in deepw[t]] for t in range(3)]
    for i in range(GT):
        m = 2 * i + 1                # y = m/2520, midpoint of grid cell i
        Uvals = [Uprime(m, u) for u in range(p)]
        Vvals = [Vprime(m, u) for u in range(p)]
        wordy = ind_scaled(Wp, m, 2)
        cellv = []
        for ell in range(qm):
            zc = (m - 504 * ell) % 2520      # y - ell/5 in m-units (/2520)
            cellv.append(1 if min(zc, 2520 - zc) < 252 else 0)
        # literal stalk factors, subsampled sheet exhaustion:
        arng = range(dp) if i % 40 == 0 else (0, 5, 26)
        for u in range(p):
            for t in range(3):
                th = (t - 2 * u) % 3
                vals = set()
                for a in arng:
                    # z = 54*X - t/3, X = (m/2520 + u + 3a)/81:
                    #   54*X - t/3 = (9n - 11340 t)/34020, n = m+2520(u+3a)
                    n = m + 2520 * (u + 3 * a)
                    zm = (9 * n - 11340 * t) % 34020
                    vals.add(1 if (zm < 2430 or zm > 31590) else 0)
                    # word at RX = 81*X = n/2520 -- must equal 1_W(y):
                    require(ind_scaled(Wp, n, 2) == wordy,
                            "toy: word factor not base-exact at stalk")
                require(len(vals) == 1,
                        "toy: deep factor sheet-DEPENDENT (impossible)")
                dv = vals.pop()
                require(dv == (1 if member(dst[th], den_[th], i) else 0),
                        "toy: slaving formula mismatch vs base window")
                if dv:
                    for q in range(p):
                        uv = Uvals[u] * Vvals[q] * wordy
                        if uv:
                            for ell in range(qm):
                                if cellv[ell]:
                                    ddef[u][q][ell][th] += uv
    ok_nonzero = 0
    for u in range(p):
        for q in range(p):
            for ell in range(qm):
                for th in range(3):
                    require(Fraction(ddef[u][q][ell][th], 243 * 27 * GT)
                            == Fraction(eng[u][q][ell][th], DEN),
                            f"toy mismatch at {(u, q, ell, th)}")
                    if eng[u][q][ell][th]:
                        ok_nonzero += 1
    require(ok_nonzero > 0, "toy table degenerate (test has no content)")
    log(f"    toy pipeline (GT=1260, p=3, Rpkt=9, depth=27, deep speed 54,")
    log(f"    cells mod 5, word set): all {p * p * qm * 3} table entries ==")
    log(f"    definitional stalk-level brute force ({ok_nonzero} nonzero):"
        f" PASS")
    log("    sheet-independence of the deep factor, base-exactness of the")
    log("    word factor at R*X, root-blindness of the cell factor, and the")
    log("    slaving theta = (t - 2u) mod 3 all verified LITERALLY on the")
    log("    toy stalk sheets: PASS")


# --------------------------------------------------------------------------
# Production slaving spot-check (route (b)): exact Fraction evaluation of
# the THM-2449 factors at sampled stalk sheets vs the base window sets
# --------------------------------------------------------------------------
def slaving_spot_check(cells, deep, Tb_iv, QW):
    rng = 987654321
    Tl = T_DEN
    csta = [[a for a, _ in c] for c in cells]
    cend = [[b for _, b in c] for c in cells]
    dsta = [[a for a, _ in deep[t]] for t in range(13)]
    dend = [[b for _, b in deep[t]] for t in range(13)]
    tsta = [a for a, _ in Tb_iv]
    tend = [b for _, b in Tb_iv]
    qsta = [[a for a, _ in QW[l]] for l in range(7)]
    qend = [[b for _, b in QW[l]] for l in range(7)]
    npts, nskip, nchk = 0, 0, 0
    ASAMP = (0, 1, 12345, DCOLL - 1)
    for _ in range(160):
        rng = (1103515245 * rng + 12345) % (1 << 62)
        ny = rng % Tl
        y = Fraction(ny, Tl)
        npts += 1
        for u in range(13):
            for t in range(13):
                th = (t - 2 * u) % 13
                # direct window test ||c3*X - t/13|| < 1/14 at sheets:
                vals = set()
                bad = False
                for a in ASAMP:
                    X = ((y + u) / 13 + a) / DCOLL
                    z = C3 * X - Fraction(t, 13)
                    z -= z.numerator // z.denominator  # frac part in [0,1)
                    dist = min(z, 1 - z)
                    if dist == Fraction(1, 14):
                        bad = True
                        break
                    vals.add(1 if dist < Fraction(1, 14) else 0)
                if bad:
                    nskip += 1
                    continue
                require(len(vals) == 1,
                        "production deep factor sheet-DEPENDENT (impossible)")
                got = vals.pop()
                want = 1 if member(dsta[th], dend[th], ny) else 0
                require(got == want,
                        f"slaving mismatch at u={u}, t={t}")
                nchk += 1
        # owner-clock factor at chart point (root-blind) + word factor:
        for l in range(7):
            z = y - Fraction(l, 7)
            z -= z.numerator // z.denominator
            dist = min(z, 1 - z)
            if dist != Fraction(1, 14):
                got = 1 if dist < Fraction(1, 14) else 0
                want = 1 if member(csta[l], cend[l], ny) else 0
                require(got == want, f"cell window mismatch at ell={l}")
            z = 13 * y - Fraction(l, 7)
            z -= z.numerator // z.denominator
            dist = min(z, 1 - z)
            if dist != Fraction(1, 14):
                gw = (1 if member(tsta, tend, ny) else 0) * \
                     (0 if dist < Fraction(1, 14) else 1)
                want = 1 if member(qsta[l], qend[l], ny) else 0
                require(gw == want, f"word skew window mismatch at ell={l}")
    log(f"    sampled {npts} exact base points x 13 roots x 13 deep labels")
    log(f"    x 4 sheets a in {{0,1,12345,13^5-1}}: {nchk} slaving identities")
    log(f"    Delta_t(X_(u,a)) == DEEP_((t-2u) mod 13)(y) verified exactly")
    log(f"    ({nskip} boundary hits skipped); cell and word skew windows")
    log("    verified against direct ||.|| Fraction tests: PASS")


# --------------------------------------------------------------------------
# Phi_91 / Q(zeta_91) layer for the decisive S
# --------------------------------------------------------------------------
def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    out[i + j] += ai * bj
    return out


def poly_divexact(a, b):
    """Exact division of integer polynomials (remainder must be 0)."""
    a = a[:]
    db = len(b) - 1
    require(b[-1] == 1, "divisor must be monic here")
    out = [0] * (len(a) - db)
    for i in range(len(a) - 1, db - 1, -1):
        c = a[i]
        if c:
            out[i - db] = c
            for j in range(db + 1):
                a[i - db + j] -= c * b[j]
    require(all(x == 0 for x in a), "polynomial division not exact")
    return out


def make_phi91():
    """Phi_91(x) = (x^91-1)(x-1) / ((x^13-1)(x^7-1)), degree 72, exact."""
    num = poly_mul([-1] + [0] * 90 + [1], [-1, 1])
    den = poly_mul([-1] + [0] * 12 + [1], [-1] + [0] * 6 + [1])
    # divide num by den (den monic):
    phi = poly_divexact(num, den)
    require(len(phi) == 73 and phi[-1] == 1, "Phi_91 degree/monic check")
    require(sum(abs(c) for c in phi) > 0, "Phi_91 nonzero")
    return phi


PHI91 = make_phi91()


def bucket_to_poly91(B):
    """13x7 bucket matrix (coeff of zeta13^i zeta7^j) -> integer polynomial
    in x = zeta_91 via zeta13 = x^7, zeta7 = x^13 (CRT exponent bijection
    (i,j) -> 7i + 13j mod 91)."""
    p = [0] * 91
    for i in range(13):
        for j in range(7):
            if B[i][j]:
                p[(7 * i + 13 * j) % 91] += B[i][j]
    return p


def poly_mod_phi91(p):
    """Reduce an integer polynomial (degree < 91) mod Phi_91; length 72."""
    a = p[:] + [0] * max(0, 73 - len(p))
    for i in range(len(a) - 1, 71, -1):
        c = a[i]
        if c:
            a[i] = 0
            for j in range(73):
                a[i - 72 + j] -= c * PHI91[j]
    return a[:72]


def bucket_mul_91(B1, B2):
    """Multiply two 13x7 bucket matrices in Z[zeta13, zeta7] (cyclic conv)."""
    out = [[0] * 7 for _ in range(13)]
    for i1 in range(13):
        row = B1[i1]
        for j1 in range(7):
            c1 = row[j1]
            if c1:
                for i2 in range(13):
                    r2 = B2[i2]
                    ii = (i1 + i2) % 13
                    orow = out[ii]
                    for j2 in range(7):
                        if r2[j2]:
                            orow[(j1 + j2) % 7] += c1 * r2[j2]
    return out


def psi_bucket(NUMT, tau, a0, alpha, beta):
    """THM-2512 (12)-(14) toothpick + double DFT on the integer defect
    dnum(v, ell) = NUMT[ell][v] (the transposed interaction numerators).
    Returns the 13x7 bucket matrix of Psi_{tau,a0}(alpha,beta)."""
    B = [[0] * 7 for _ in range(13)]
    for c in range(7):
        for v in range(13):
            acc = 0
            for l in range(7):
                acc += NUMT[l][(v - tau * ((a0 * l + c) % 7)) % 13]
            if acc:
                B[(-alpha * v) % 13][(-beta * c) % 7] += acc
    return B


def dtilde_bucket(NUMT, alpha, gamma):
    """dtilde(alpha, gamma) = sum_{s,ell} dnum(s,ell) z13^(-alpha s)
    z7^(-gamma ell)."""
    B = [[0] * 7 for _ in range(13)]
    for l in range(7):
        for s in range(13):
            c = NUMT[l][s]
            if c:
                B[(-alpha * s) % 13][(-gamma * l) % 7] += c
    return B


def K_bucket(u, beta):
    """K_{u,beta} = sum_{j=0}^{6} (z13^{-u} z7^{-beta})^j."""
    B = [[0] * 7 for _ in range(13)]
    for j in range(7):
        B[(-u * j) % 13][(-beta * j) % 7] += 1
    return B


def interaction_nums(Anum):
    """Anum: 7x13 integer numerators over a common den.  Returns the 7x13
    integer numerators of 91*I_(ell,th) over the same den (THM-2512 (3))."""
    rows = [sum(r) for r in Anum]
    cols = [sum(Anum[l][t] for l in range(7)) for t in range(13)]
    tot = sum(rows)
    return [[91 * Anum[l][t] - 7 * rows[l] - 13 * cols[t] + tot
             for t in range(13)] for l in range(7)]


# ==========================================================================
def main():
    log("STAGE-2 THETA-SLAVED CONTRACTION at the r = 5 window"
        " (canonical typed row, sigma={b}, K=2)")
    log("script=04-computation/lrc14_stage2_theta_contraction_opus"
        "_20260728.py")
    log("machine=opus  script-date=2026-07-28  run-date=2026-07-27")
    log("")
    log("[1] object and logged conventions")
    log(f"    row w={W}   profile nu13(c1,c2,c3)=(1,3,5)")
    log(f"    owner j=1 (c_1={C1}); word sigma={{b}}; packet clock K=2"
        f" (R_pkt=169)")
    log(f"    e = 1_(E_1);  f = 1_(Q_1,{{b}}) P^2 1_(E_1)   (THM-2471 (23))")
    log(f"    r = 5 (THM-2575/2577, all-clock law), d = 13^5 = {DCOLL}:")
    log(f"    U = P_d f, V = P_d e;  A(y,u)=U((y+u)/13), F(y,q)=V((y+q)/13)")
    log(f"    sidecar: delta=gcd(c_3,d)={DCOLL}, d_0=1, h=c_3/d={HDEEP},"
        f" 13 does not divide h")
    log("    => unique affine invariant theta = t - 2u (THM-2471 (46)-(47)),")
    log("       u = CURRENT-branch collision root (the branch of (45))")
    log("    stalk sheets (THM-2471 (28)): w_u=(y+u)/13, X_(u,a)=(w_u+a)/d")
    log("    THM-2449-shaped Boolean factors on LINKED NODES OF ONE FINITE")
    log("    ANCESTRY STALK (not one circle point; MISTAKE-293): exact")
    log("    identities re-verified in [2],[3]")
    log("      owner clock d(13x - ell/7) at x=w_u: 13w_u = y+u == y (mod 1)")
    log("        => cell_ell(y) = 1_(||y - ell/7|| < 1/14)   [root-blind]")
    log("      deep probe d(c_3 x - t/13) at x=X_(u,a):")
    log("        c_3 X_(u,a) = 2(w_u+a) == 2(y+u)/13 (mod 1) [sheet-indep]")
    log("        => Delta_t(X_(u,a)) = DEEP_theta(y), theta=(t-2u) mod 13,")
    log("           DEEP_0=[0,13/28), DEEP_1=(1/28,27/28), DEEP_2=(15/28,1),")
    log("           DEEP_theta = EMPTY for theta in 3..12;")
    log("           law: sum_theta DEEP_theta = 2 - d(2y mod 1)  a.e.")
    log("           (per root u exactly t in {2u,2u+1,2u+2} mod 13 can fire)")
    log("      word factor Q^eps_ell(R x) at x=X_(u,a), response clock"
        " R=13^k:")
    log(f"        CLOCK ARITHMETIC (THM-2449 (25)-(31) uses this class only")
    log(f"        eventually; no assertion that k=6 exceeds its k_0):")
    log(f"        T_DEN = 13^6 * D0, D0 = {D0}, ord_D0(13) = {ORD};")
    log(f"        class k == 6 (mod {ORD}); smallest exact sheet-return")
    log(f"        exponent k = {KRESP}, evaluated directly; eps=+{EPS}")
    log("        at k=6: R X_(u,a) = 13(w_u+a) = y + u + 13a == y (mod 1)")
    log("        (sheet-exact AND root-blind; unique such clock), so")
    log("        Q^(+1)_ell(R X) = QW_ell(y) = T_b(y)(1 - d(13y - ell/7)),")
    log("        T_b = A_0 n D_c3 \\ D_c2 (source-safe factor deleted,")
    log("        THM-2449 (15); THM-2409 skew-diagonal mandatory)")
    log("    THE ONE TABLE (one integral over the common base y per entry,")
    log("    MISTAKE-281: paired BEFORE any marginalization/DFT):")
    log("      N(u,q,ell,theta) = int_T A(y,u) F(y,q) cell_ell(y)")
    log("                          DEEP_theta(y) QW_ell(y) dy")
    log("      u = current root, q = source root, s = u - q (colour shift)")
    log("    danger d(z)=1_(||z||<1/14); half-open [A,B) realizations on the")
    log(f"    T_DEN = 182*lcm(w) = {T_DEN} grid (boundary differences are")
    log("    null sets; step-function integrals unaffected); no floats in")
    log("    any decision")
    log("    DFT signs: J(k)-side zeta13^(-ks) (THM-2471 (8)); toothpick")
    log("    R_(tau,a0,c)(v) = sum_ell d(v - tau*rep(a0*ell+c), ell),")
    log("    Psi_(tau,a0)(alpha,beta) = sum_(c,v) R(v) z13^(-alpha v)")
    log("    z7^(-beta c), rep: F_7 -> {0..6} (THM-2508/2512 (12)-(14));")
    log("    OUTPUT INDEX v := theta (the slaved deep invariant)")
    log("")

    log("[2] machinery self-test (hostile toy model, stalk-level"
        " definitional)")
    toy_test()
    log("")

    log("[3] exact interval sets, anchors, and refiner laws")
    E = build_set(PAT_E, ZELL)
    lenE = check_intervals(E, T_DEN)
    require(Fraction(lenE, T_DEN) == Fraction(1882176, 28589561),
            "mu(E_1) anchor")
    require(len(E) == 57072, "E_1 interval count anchor")
    log(f"    E_1: intervals={len(E)}  measure={Fraction(lenE, T_DEN)}:"
        f" PASS")
    QB = build_set(PAT_QB, ZELL)
    lenQB = check_intervals(QB, T_DEN)
    require(Fraction(lenQB, T_DEN) == Fraction(4839079319, 190921088358),
            "mu(Q_1,{b}) anchor (census)")
    require(len(QB) == 131762, "Q_1,{b} interval count anchor (census)")
    log(f"    Q_1,{{b}}: intervals={len(QB)}"
        f"  measure={Fraction(lenQB, T_DEN)}: PASS (census anchors)")
    TB = build_set(PAT_TB, ZELL)
    lenTB = check_intervals(TB, T_DEN)
    log(f"    T_b (word, source-safe deleted): intervals={len(TB)}"
        f"  measure={Fraction(lenTB, T_DEN)}")
    require(lenTB > lenQB, "T_b must strictly contain Q_b in measure")
    QB2 = subtract_comb(TB, C1, 182, -13, 13)
    require(merge_touch(QB2) == merge_touch(QB),
            "T_b minus D_c1 comb != Q_1,{b} (word build inconsistent)")
    log("    gate: T_b minus the D_c1 comb == Q_1,{b} interval-exact: PASS")
    cells7 = build_cells()
    allc = []
    for l in range(7):
        require(measure(cells7[l]) == T_DEN // 7, "cell measure 1/7")
        allc += cells7[l]
    require(merge_touch(allc) == [(0, T_DEN)], "cells do not tile T")
    for l1 in range(7):
        for l2 in range(l1 + 1, 7):
            require(intersect_lists(cells7[l1], cells7[l2]) == [],
                    "cells overlap")
    log("    owner-clock cells (at the stalk): 7 half-open cells of measure")
    log("    1/7 tile the circle exactly (a.e.): PASS")
    deep13 = build_deep()
    d2set = build_d2()
    check_deep_law(deep13, d2set)
    require(sum(measure(deep13[t]) for t in range(13)) ==
            2 * T_DEN - measure(d2set), "deep law measure")
    log("    slaved deep windows: DEEP_theta empty for theta in 3..12;")
    log("    exact law sum_theta DEEP_theta = 2 - d(2y mod 1) verified")
    log("    piecewise on every breakpoint interval: PASS")
    QW = [subtract_comb(TB, C1, 182, 26 * l - 13, 26 * l + 13)
          for l in range(7)]
    for l in range(7):
        check_intervals(QW[l], T_DEN)
    log(f"    word factors QW_ell = T_b * (1 - d(13y - ell/7)):"
        f" intervals={[len(QW[l]) for l in range(7)]}")
    GW_set = [intersect_lists(QW[l], cells7[l]) for l in range(7)]
    G_set = [[intersect_lists(GW_set[l], deep13[t]) for t in range(3)]
             for l in range(7)]
    GCD2_set = [intersect_lists(GW_set[l], d2set) for l in range(7)]
    log(f"    refiner products GW=QW n cell, G=GW n DEEP, GCD2=GW n d2"
        f" built;")
    log(f"    G interval counts per ell: "
        f"{[[len(G_set[l][t]) for t in range(3)] for l in range(7)]}")
    log("")
    log("    production slaving spot-check (route (b)):")
    slaving_spot_check(cells7, deep13, TB, QW)
    log("")

    log("[4] exact profiles f, U, V and their anchors")
    n_st, n_v = weighted_fold([(a, b, 1) for a, b in E], RPKT, T_DEN)
    fp = []
    np_ = len(n_st)
    for qa, qb in QB:
        i = bisect_right(n_st, qa) - 1
        while True:
            pa = n_st[i]
            pb = n_st[i + 1] if i + 1 < np_ else T_DEN
            a = qa if qa > pa else pa
            b = qb if qb < pb else pb
            if a < b and n_v[i]:
                fp.append((a, b, n_v[i]))
            if pb >= qb:
                break
            i += 1
    int_f_num = sum(w * (b - a) for a, b, w in fp)     # over 169*T_DEN
    require(Fraction(int_f_num, RPKT * T_DEN) == RHO_B,
            "int f != census rho^(2)({b}) anchor")
    log(f"    {elapsed()} f = 1_Q P^2 1_E: {len(fp)} weighted pieces;")
    log(f"    int f = {Fraction(int_f_num, RPKT * T_DEN)}"
        f" == census rho^(2)({{b}}): PASS")
    Ust, Uv = weighted_fold(fp, DCOLL, T_DEN)          # num/den 169*13^5
    Vst, Vv = weighted_fold([(a, b, 1) for a, b in E], DCOLL, T_DEN)
    require(profile_mass(Ust, Uv, T_DEN) == DCOLL * int_f_num,
            "int U != int f")
    require(profile_mass(Vst, Vv, T_DEN) == DCOLL * lenE, "int V != mu(E)")
    log(f"    {elapsed()} U = P_(13^5) f: {len(Ust)} pieces"
        f" (num/den {RPKT * DCOLL});")
    log(f"    V = P_(13^5) e: {len(Vst)} pieces (num/den {DCOLL});")
    log("    int U == int f and int V == mu(E_1) (P_m preserves integrals):"
        " PASS")
    pu = product_cum(Ust, Uv, Vst, Vv, T_DEN)
    require(pu[3] == 0 and all(v == 0 for v in pu[1]),
            "U*V not identically zero (contradicts first collision r=5)")
    log("    U*V == 0 IDENTICALLY as a profile (first-collision"
        " disjointness,")
    log("    THM-2471 (26)): PASS")
    UWl = [extract_window(Ust, Uv, u, 13, T_DEN) for u in range(13)]
    VWl = [extract_window(Vst, Vv, u, 13, T_DEN) for u in range(13)]
    log(f"    {elapsed()} 13 chart windows extracted for U and V")
    log("")

    DENC = RPKT * DCOLL * DCOLL * T_DEN
    log("[5] THE ONE TABLE: 169 pair profiles x 46 refiner integrals"
        " (exact)")
    log(f"    entry denominators: DENC = 169 * 13^10 * T_DEN = {DENC}")
    log(f"    DENC factored: {factor_str(DENC)}")
    P = [[0] * 13 for _ in range(13)]
    NG = [[None] * 13 for _ in range(13)]     # NG[u][q][l][th], th<3
    NGW = [[None] * 13 for _ in range(13)]
    NGC = [[None] * 13 for _ in range(13)]
    aggC = [0] * 7
    aggDeep = [0] * 3
    aggD2 = 0
    zero73 = [[0] * 3 for _ in range(7)]
    for u in range(13):
        for q in range(13):
            prof = product_cum(UWl[u][0], UWl[u][1],
                               VWl[q][0], VWl[q][1], T_DEN)
            ptot = prof[3]
            P[u][q] = ptot
            if u == q:
                require(ptot == 0 and all(v == 0 for v in prof[1]),
                        f"diagonal pair (u=q={u}) not identically zero")
                NG[u][q] = zero73
                NGW[u][q] = [0] * 7
                NGC[u][q] = [0] * 7
                continue
            cellv = [set_integral(prof, cells7[l]) for l in range(7)]
            require(sum(cellv) == ptot,
                    f"cell partition gate fails at pair ({u},{q})")
            d2v = set_integral(prof, d2set)
            deepv = [set_integral(prof, deep13[t]) for t in range(3)]
            require(sum(deepv) == 2 * ptot - d2v,
                    f"deep root-count gate fails at pair ({u},{q})")
            gw = [set_integral(prof, GW_set[l]) for l in range(7)]
            gc = [set_integral(prof, GCD2_set[l]) for l in range(7)]
            gg = [[set_integral(prof, G_set[l][t]) for t in range(3)]
                  for l in range(7)]
            for l in range(7):
                require(sum(gg[l]) == 2 * gw[l] - gc[l],
                        f"refined deep gate fails at ({u},{q},ell={l})")
            NG[u][q] = gg
            NGW[u][q] = gw
            NGC[u][q] = gc
            for l in range(7):
                aggC[l] += cellv[l]
            for t in range(3):
                aggDeep[t] += deepv[t]
            aggD2 += d2v
        log(f"    {elapsed()} current root u={u:2d}: 13 pairs done"
            " (gates PASS)")
    totP = sum(P[u][q] for u in range(13) for q in range(13))
    require(Fraction(totP, DENC) == 169 * I5,
            "sum_(u,q) P != 169*I_5 (THM-2471 (7) service anchor)")
    log("    per-pair gates (cell partition, deep root-count law, refined")
    log("    deep law) all PASS at all 169 pairs; diagonal pairs vanish")
    log("    IDENTICALLY (zero-sum law (g2)): PASS")
    log(f"    total service sum_(u,q) int A F = 169*I_5 ="
        f" {169 * I5}: PASS")
    log("")

    log("[6] route-2 independent aggregate (folded service profiles)")
    Upieces = [(Ust[i], Ust[i + 1] if i + 1 < len(Ust) else T_DEN, Uv[i])
               for i in range(len(Ust))]
    Vpieces = [(Vst[i], Vst[i + 1] if i + 1 < len(Vst) else T_DEN, Vv[i])
               for i in range(len(Vst))]
    aSt, aV = weighted_fold(Upieces, 13, T_DEN)   # a(y) = sum_u A(y,u)
    bSt, bV = weighted_fold(Vpieces, 13, T_DEN)   # b(y) = sum_q F(y,q)
    prof_ab = product_cum(aSt, aV, bSt, bV, T_DEN)
    require(prof_ab[3] == totP, "route-2: int a*b != sum of pair totals")
    for l in range(7):
        require(set_integral(prof_ab, cells7[l]) == aggC[l],
                f"route-2 cell mismatch at ell={l}")
    require(set_integral(prof_ab, d2set) == aggD2, "route-2 d2 mismatch")
    for t in range(3):
        require(set_integral(prof_ab, deep13[t]) == aggDeep[t],
                f"route-2 deep mismatch at theta={t}")
    for l in range(7):
        require(set_integral(prof_ab, GW_set[l]) ==
                sum(NGW[u][q][l] for u in range(13) for q in range(13)),
                f"route-2 GW mismatch at ell={l}")
        require(set_integral(prof_ab, GCD2_set[l]) ==
                sum(NGC[u][q][l] for u in range(13) for q in range(13)),
                f"route-2 GCD2 mismatch at ell={l}")
        for t in range(3):
            require(set_integral(prof_ab, G_set[l][t]) ==
                    sum(NG[u][q][l][t] for u in range(13)
                        for q in range(13)),
                    f"route-2 G mismatch at (ell={l},theta={t})")
    log(f"    {elapsed()} a = fold_13 U ({len(aSt)} pieces), b = fold_13 V"
        f" ({len(bSt)} pieces);")
    log("    EVERY aggregate refiner integral recomputed from the a*b")
    log("    profile (different pairing geometry) == pairwise sums: PASS")
    log("")

    log("[7] gate (g1): root-marginal colours of the r=5 stalk")
    Cs = [sum(P[(q + s) % 13][q] for q in range(13)) for s in range(13)]
    require(Cs[0] == 0, "C(0) != 0")
    require(all(c >= 0 for c in Cs), "negative C(s)")
    require(Fraction(sum(Cs), DENC) == 169 * I5,
            "sum_s C(s) != 169*I_5 (THM-2471 (14))")
    DEN_J = 169 * DENC
    JI = []
    for k in range(13):
        bkt = [0] * 13
        for s in range(13):
            bkt[(-k * s) % 13] += Cs[s]
        JI.append(bkt)
    require(canon(JI[0]) == canon([sum(Cs)] + [0] * 12), "J(0) bucket shape")
    require(Fraction(sum(Cs), DEN_J) == I5,
            "J(0) != I_5 ANCHOR (THM-2575/2577 triple value)")
    log("    J(0) = I_5 = 48602521488933856/337437093630814766589: PASS"
        "   (ANCHOR GATE)")
    accJ = [0] * 13
    for k in range(13):
        accJ = badd(accJ, JI[k])
    require(is_zero_b(accJ), "sum_k J(k) != 0")
    nzJ = sum(1 for k in range(1, 13) if not is_zero_b(JI[k]))
    require(nzJ == 12, "some nonzero colour J(k) vanishes")
    for k in range(13):
        require(canon(JI[(13 - k) % 13]) == canon(bconj(JI[k])),
                "Hermitian symmetry broken")
    log("    C(0) = 0; sum_s C(s) = 169*I_5; J(k) != 0 for all 12 k != 0;")
    log("    sum_k J(k) = 0; Hermitian J(13-k) = conj J(k)  (THM-2471")
    log("    (9)-(10) at r=5): ALL PASS")
    log(f"    C(s) exact numerators over DENC: {Cs}")
    log(f"    J(k) canonical coords over DEN_J = 169*DENC = {DEN_J}:")
    for k in range(13):
        log(f"      k={k:2d}: {list(canon(JI[k]))}")
    log("")
    return (E, lenE, QB, TB, cells7, deep13, d2set, QW, GW_set, G_set,
            GCD2_set, P, NG, NGW, NGC, Cs, JI, DENC, totP)


def peval(poly, w, mod):
    acc = 0
    for c in reversed(poly):
        acc = (acc * w + c) % mod
    return acc


def stage2(state):
    (E, lenE, QB, TB, cells7, deep13, d2set, QW, GW_set, G_set,
     GCD2_set, P, NG, NGW, NGC, Cs, JI, DENC, totP) = state

    log("[8] the response array, fibrewise centring, and the aggregate"
        " interaction")
    Anum = [[0] * 13 for _ in range(7)]
    for u in range(13):
        for q in range(13):
            gg = NG[u][q]
            for l in range(7):
                for t in range(3):
                    Anum[l][t] += gg[l][t]
    log("    Aresp[ell][theta] = sum_(u,q) N(u,q,ell,theta)  (Fubini inside")
    log("    the ONE integral; equals int a b cell_ell DEEP_theta QW_ell dy")
    log("    by route (g4)); columns theta = 3..12 are EXACT zeros (empty")
    log("    windows).  Exact numerators over DENC:")
    for l in range(7):
        log(f"      ell={l}: theta=0..2: {Anum[l][:3]}  | theta=3..12: all 0")
    require(all(Anum[l][t] == 0 for l in range(7) for t in range(3, 13)),
            "theta>=3 columns must be zero")
    require(any(Anum[l][t] for l in range(7) for t in range(3)),
            "response array degenerate (all zero)")
    # fibrewise centring (per (u,s) ancestry fibre; linearity check)
    NUM = interaction_nums(Anum)
    NUMsum = [[0] * 13 for _ in range(7)]
    fib_nz = 0
    fib_tot = 0
    per_s_nz = [0] * 13
    for u in range(13):
        for s in range(13):
            q = (u - s) % 13
            Af = [[NG[u][q][l][t] if t < 3 else 0 for t in range(13)]
                  for l in range(7)]
            NUMf = interaction_nums(Af)
            fib_tot += 1
            if any(NUMf[l][t] for l in range(7) for t in range(13)):
                fib_nz += 1
                per_s_nz[s] += 1
            for l in range(7):
                for t in range(13):
                    NUMsum[l][t] += NUMf[l][t]
    require(NUMsum == NUM,
            "sum of fibrewise interactions != aggregate interaction")
    log(f"    fibrewise ANOVA (per (u,s) ancestry fibre): {fib_nz} of"
        f" {fib_tot} fibres")
    log(f"    carry nonzero interaction; per shift s: {per_s_nz}")
    log("    (s = 0 fibres vanish identically -- zero-sum law (g2))")
    log("    linearity: sum of fibrewise interactions == centring of the")
    log("    aggregate (all centring terms marginals of the ONE table):"
        " PASS")
    DEN_I = 91 * DENC
    log(f"    aggregate interaction numerators 91*I_(ell,theta) over"
        f" DEN_I = 91*DENC:")
    for l in range(7):
        log(f"      ell={l}: {NUM[l]}")
    inz = any(NUM[l][t] for l in range(7) for t in range(13))
    log(f"    aggregate interaction nonzero: {inz}")
    log("")

    log("[9] THE DECISIVE THETA-SLAVED TOOTHPICK CONTRACTION")
    NUMT = NUM                      # dnum(v, ell) = NUM[ell][v] (transpose
    #                                 read happens inside psi_bucket)
    TAU, A0, ALPHA, BETA = 1, 1, 1, 1
    log(f"    fixed primitive quadruple (tau,a0,alpha,beta) ="
        f" ({TAU},{A0},{ALPHA},{BETA})")
    log("    (first primitive quadruple in the fixed chart order; chosen")
    log("    from the residue-chart conventions alone -- independent of")
    log("    every ancestry point, row datum, and table value)")
    log("    output index v = theta (the slaved deep invariant t - 2u):")
    log("    S := Psi_(1,1)(1,1) = sum_(c,v,ell) d(v - rep(ell+c), ell)")
    log("         z13^(-v) z7^(-c),   d(v,ell) = I_(ell,v)")
    Bpsi = psi_bucket(NUMT, TAU, A0, ALPHA, BETA)
    Praw = bucket_to_poly91(Bpsi)
    Sred = poly_mod_phi91(Praw)
    SNZ = any(c != 0 for c in Sred)
    DEN_S = 91 * DENC
    log(f"    S in Q(zeta_91), numerator reduced mod Phi_91 (deg 72 basis")
    log(f"    1, x, ..., x^71; x = zeta_91, zeta13 = x^7, zeta7 = x^13),")
    log(f"    over DEN_S = 91*DENC = {DEN_S}")
    log(f"    DEN_S factored: {factor_str(DEN_S)}")
    log("    reduced numerator coefficients:")
    for i in range(0, 72, 8):
        log(f"      x^{i:2d}..x^{min(i + 7, 71):2d}: {Sred[i:i + 8]}")
    # numeric display only (never a decision):
    import math as _m
    re_ = sum(c * _m.cos(2 * _m.pi * i / 91) for i, c in enumerate(Sred))
    im_ = sum(c * _m.sin(2 * _m.pi * i / 91) for i, c in enumerate(Sred))
    log(f"    ~ S = ({re_ / DEN_S:.6e}) + ({im_ / DEN_S:.6e}) i"
        "   [display only]")
    log(f"    DECISION: S {'!=' if SNZ else '=='} 0   (exact reduction"
        " mod Phi_91)")
    # THM-2512 (15) factorization identity:
    KB = K_bucket((ALPHA * TAU) % 13, BETA)
    DB = dtilde_bucket(NUMT, ALPHA, (-BETA * A0) % 7)
    PRB = bucket_mul_91(KB, DB)
    require(poly_mod_phi91(bucket_to_poly91(PRB)) == Sred,
            "THM-2512 (15) factorization identity fails")
    log("    THM-2512 (15) identity Psi = K_(alpha tau,beta) *"
        " dtilde(alpha,-beta a0):")
    log("    verified exactly in Z[x]/Phi_91: PASS  (positive control)")
    inz = any(NUM[l][t] for l in range(7) for t in range(13))
    require(SNZ == inz,
            "THM-2512 (19) consistency broken: S-nonzero != I-nonzero")
    log(f"    THM-2512 (19) consistency (rational table, primitive"
        f" quadruple):")
    log(f"    S != 0  <=>  interaction != 0:  both {SNZ}: PASS")
    for (t2, a2, al2, be2) in ((2, 3, 5, 4), (12, 6, 12, 6), (7, 2, 3, 1)):
        B2 = psi_bucket(NUMT, t2, a2, al2, be2)
        S2 = poly_mod_phi91(bucket_to_poly91(B2))
        require(any(S2) == inz,
                f"quadruple ({t2},{a2},{al2},{be2}) breaks the Galois"
                " all-or-all law")
        K2 = K_bucket((al2 * t2) % 13, be2)
        D2b = dtilde_bucket(NUMT, al2, (-be2 * a2) % 7)
        require(poly_mod_phi91(bucket_to_poly91(bucket_mul_91(K2, D2b)))
                == S2, "factorization fails at spot quadruple")
        log(f"    spot quadruple ({t2},{a2},{al2},{be2}):"
            f" nonzero={bool(any(S2))}, factorization PASS")
    if SNZ:
        log("    => by THM-2512 (17)/(19) (rational table, coprime Galois")
        log("       transitivity), ALL 5184 primitive cut coefficients of")
        log("       this theta-slaved bundle are nonzero.")
    log("")

    log("[10] split-prime cross-checks (547, 911, 1093, 2003, 2549)")
    nz_primes = 0
    for PP in (547, 911, 1093, 2003, 2549):
        require((PP - 1) % 91 == 0, f"{PP} not split in Q(zeta_91)")
        w = None
        for g in range(2, 100):
            cand = pow(g, (PP - 1) // 91, PP)
            if pow(cand, 7, PP) != 1 and pow(cand, 13, PP) != 1:
                w = cand
                break
        require(w is not None and pow(w, 91, PP) == 1,
                f"no order-91 element found mod {PP}")
        require(peval(PHI91, w, PP) == 0, f"Phi_91(w) != 0 mod {PP}")
        vr = peval(Praw, w, PP)
        vd = peval(Sred, w, PP)
        require(vr == vd, f"raw/reduced evaluation mismatch mod {PP}")
        if vr != 0:
            nz_primes += 1
        log(f"    p={PP}: w(order 91)={w},  S mod p = {vr}"
            f"  ({'NONZERO' if vr else 'zero'})")
    if SNZ:
        require(nz_primes > 0,
                "S != 0 but vanished at all five split primes (implausible;"
                " abort)")
        log(f"    {nz_primes}/5 split primes certify S != 0 independently:"
            " PASS")
    else:
        require(nz_primes == 0, "S == 0 but a split prime disagrees")
        log("    all five split primes consistent with S == 0: PASS")
    log("")

    log("[11] hostile controls (both must vanish exactly)")
    B0 = psi_bucket(NUMT, TAU, A0, ALPHA, 0)
    S0 = poly_mod_phi91(bucket_to_poly91(B0))
    require(all(c == 0 for c in S0), "H1: beta=0 coefficient must vanish")
    log("    H1 (beta = 0): Psi_(1,1)(1,0) = 0 exactly -- the row-zero law")
    log("    THM-2512 (4) kills the untwisted cut phase: PASS (vanishes)")
    Anum13 = [[sum(Anum[l]) for _t in range(13)] for l in range(7)]
    NUMflat = interaction_nums(Anum13)
    require(all(NUMflat[l][t] == 0 for l in range(7) for t in range(13)),
            "H2: flat-slaved interaction must vanish identically")
    Bf = psi_bucket(NUMflat, TAU, A0, ALPHA, BETA)
    Sf = poly_mod_phi91(bucket_to_poly91(Bf))
    require(all(c == 0 for c in Sf), "H2: flat-slaved S must vanish")
    log("    H2 (u-independent slaving): replacing the slaved column family")
    log("    by its theta-average (deep label decoupled from the root, the")
    log("    r<=4 base-only shape) kills the interaction AND the contraction")
    log("    identically: PASS (vanishes)")
    log("    => the nonzero signal (if any) is carried by the theta = t - 2u")
    log("       coupling itself, not by the ell-structure alone.")
    log("")

    log("[12] positive controls (marginals reproducing known exact values)")
    log("    * J(0) = I_5 (THM-2575/2577 anchor): PASS  ([7])")
    log("    * sum_s C(s) = 169*I_5 (THM-2471 (14)): PASS  ([5],[7])")
    log("    * int f = rho^(2)({b}) census anchor: PASS  ([4])")
    log("    * mu(E_1), |E_1|, mu(Q_1,{b}), |Q_1,{b}| anchors: PASS  ([3])")
    log("    * THM-2512 (15) factorization identity at 4 quadruples: PASS")
    log("    * route-2 aggregate equality (46 integrals): PASS  ([6])")
    log("")

    log("[13] scope and verdict (MISTAKE-281/283/286/287/293 typing)")
    log("    * TYPED row THM-2309 (25); NOT an asserted scalar cover.")
    log("    * Every pairing in this run sits INSIDE one integral over one")
    log("      common base point before any marginalization or DFT")
    log("      (MISTAKE-281); fibre labels (u,s,ell,theta) are carried, and")
    log("      all centring terms are marginals of the one table.")
    log("    * The owner cell lives at w_u, the deep/word factors at X_(u,a),")
    log("      and the source at Y_(q,e').  These are linked nodes of one")
    log("      finite Boolean stalk, NOT one circle point and NOT THM-2449's")
    log("      one-variable H^R table (MISTAKE-293).  k=6 is direct-exact;")
    log("      no claim that it exceeds THM-2449's eventual k_0.")
    log("    * The theta = t - 2u map is THM-2471 (44)-(47)'s deep-root")
    log("      sidecar invariant at the r = 5 window -- THM-2471's candidate")
    log("      realization, NOT 'THM-2512's bridge' (MISTAKE-286): THM-2512")
    log("      Section 5 is generic and contains no r, theta, or clock.")
    log("    * No physical current, no row exclusion, no LRC(14)/JC(2)")
    log("      progress is claimed.  The live ledger remains 165.")
    if SNZ:
        log("    VERDICT: S != 0 (exact).  This is the FIRST realized")
        log("    theta-slaved toothpick contraction pairing that is nonzero")
        log("    on a common ancestry base at the r = 5 window: a positive")
        log("    instance of THM-2471's deep-root-sidecar candidate")
        log("    realization on the canonical typed row at (j=1, sigma={b},")
        log("    K=2, response clock k=6).  All 5184 primitive cut")
        log("    coefficients of the slaved bundle fire (THM-2512 (17)/(19)")
        log("    applied to this rational array).  It is NOT a THM-2512")
        log("    bridge theorem and NOT a physical current.")
    else:
        log("    VERDICT: S == 0 (exact).  The theta-slaved contraction")
        log("    cancels exactly at this window: a sharp negative for the")
        log("    theta ansatz on this row/word/clock, redirecting the")
        log("    program to other realizations (other rows, clocks, or")
        log("    non-affine couplings).")
    log("")
    log(f"{elapsed()} all checks passed")


if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    state = main()
    stage2(state)
