#!/usr/bin/env python3
"""THE BASE-ONLY BRIDGE PAIRING on the canonical typed row (THM-2560 redirect (A)).

TASK.  THM-2560 (A) proved that on the canonical typed row THM-2309 (25)

    w = (1,14,27,40,53,66,13,2197,742586),  owner j = 1, word sigma = {a},

the first owner-normalized collision is at depth r = 3 (d = 13^3 = 2197) and
the deep phase is BASE-ONLY (h = C/d = 338, 13 | h): the lawful coupling to
the deep leg runs through the base coordinate y alone.  The redirect the
files license is a base-only bridge: pair the THM-2471 (8) collision-colour
data of the r = 3 stalk against a lawful y-side cut probe carrying an F_7
label.  Per HYP-9032 P2 the natural probe is the OWNER-CLOCK partition: the
seven translates d(c_1 y - ell/7) of THM-2449 (14) at j = 1 (c_1 = 13),
which partition the circle a.e.

OBJECT.  With the SAME packet as the audited realization companion
(04-computation/lrc14_ancestry_realization_stage_opus_20260727.py):

    e = 1_{E_1},   f = 1_{Q_{1,{a}}} P^2 1_{E_1},   d = 13^3 = 2197,
    U = P_d f,     V = P_d e            (THM-2471 (25), r = 3),
    A(y,u) = U((y+u)/13),  F(y,u) = V((y+u)/13),  u in F_13   (THM-2471 (5)),
    alpha_k(y) = (1/13) sum_u A(y,u) zeta_13^(-ku),
    phi_k(y)   = (1/13) sum_u F(y,u) zeta_13^(-ku)            (THM-2471 (8)),

define, for k in F_13 and ell in F_7,

    W(k,ell) = int_T alpha_k(y) conj(phi_k(y)) 1_{cell_ell}(y) dy,

    cell_ell = { y : d(13 y - ell/7) = 1 }
             = { y : frac(13 y) in ((2 ell - 1)/14, (2 ell + 1)/14) mod 1 }.

DISINTEGRATION CONVENTION (logged; THM-2471's exact definitions).
  * zeta_13 = a fixed abstract primitive 13th root of unity; all cyclotomic
    arithmetic is exact in Z[zeta_13] via the 12-dim basis 1,zeta,..,zeta^11
    with zeta^12 = -(1 + zeta + ... + zeta^11)  (reduction mod Phi_13).
  * transform sign: alpha_k carries zeta^(-ku); the pairing conjugates phi_k;
    both amplitudes carry the 1/13 normalization of THM-2471 (8).
  * pointwise colour density:  j_k(y) := alpha_k(y) conj(phi_k(y))
      = (1/169) sum_{s in F_13} c_s(y) zeta_13^(-ks),
      c_s(y) := sum_u A(y,u+s) F(y,u)      (the THM-2471 (13) integrand),
    so J(k) = int_T j_k(y) dy = (1/169) sum_s C(s) zeta^(-ks), C(s) = int c_s,
    and J(0) = C-total/169 = I_r  (THM-2471 (7),(14): sum_s C(s) = 169 I_r).
  * cellwise:  C_ell(s) := int_{cell_ell} c_s(y) dy   and
      W(k,ell) = (1/169) sum_s C_ell(s) zeta_13^(-ks).
    EVERYTHING is one integral over the common base y (the physical circle):
    the ell-restriction happens INSIDE the y-integral, before the finite
    root-DFT and before any marginalization (MISTAKE-281/283 discipline:
    pairing on one common base BEFORE marginalization).
  * cells realized half-open [ (2ell-1)/14, (2ell+1)/14 ) pulled back; the
    probe's open windows differ on a null set; step-function integrals agree.
  * grid: all sets are unions of half-open [A,B) with integer endpoints on
    the T_DEN = 182*lcm(w) grid; profiles U,V are integer-valued numerators
    over 13^5 resp. 13^3; every C_ell(s) is an integer over
    DENC = 169 * 2197^2 * T_DEN = 13^8 * T_DEN.  No floats in any decision.

QUESTIONS ANSWERED (exactly):
  (i)  both marginals of W: sum_{ell} W(k,ell) = J(k)  (GATE: must reproduce
       the THM-2471 colours; J(0) must equal the triple-verified
       I_3 = 9926558757352/109707098520974955 of THM-2560), and
       sum_{k} W(k,ell)  (= C_ell(0)/13, which is 0 iff the pointwise
       disjointness A(y,u)F(y,u) = 0 survives cellwise -- it must, since
       UV = 0 a.e. at the first collision).
  (ii) is W genuinely (k,ell)-mixed?  Exact tests: ell-independence,
       k-independence, and full rank-1 (product) test via all 1638 2x2
       minors over the field Q(zeta_13).
  If mixed: test the three HYP-9032 transplant inputs on the centred array
       Wc(k,ell) = W(k,ell) - J(k)/7:
       (T1) row-zero over F_7 (holds by construction) AND column-zero over
            F_13 (holds because sum_k W(.,ell) = 0 and sum_k J(k) = 0 --
            the exact analogue of THM-2512 (4)+(5));
       (T2) h-independence dichotomy: k-independent => 0 (structural, via
            column-zero, same argument as THM-2512 (5)); plus the instance;
       (T3) denominator floor: exact universal denominator of the entries'
            zeta-basis coordinates, factored, and min nonzero |coordinate|.

VERIFICATION MATRIX (independent routes, all exact):
  (a) toy brute-force self-test of the ENTIRE pipeline (grid 90, p = 3,
      packet clock 9, depth 27, clock cells mod 5): definitional midpoint
      Fraction evaluation of C_ell(s) == both engine routes;
  (b) anchors: mu(E_1), mu(Q_{1,{a}}), interval counts, int f = mu(F) =
      21376087/17907461390, int U = int f, int V = mu(E_1) (P_m preserves
      integrals; fold mass conservation asserted inside every fold);
  (c) cell partitions: both the 7x13-interval owner-clock partition of T and
      its 7x169-interval pullback partition are exact tilings;
  (d) ROUTE 1 (window route): extract A(.,u), F(.,u) as exact profiles from
      the 13 chart windows, per-pair product profiles, per-cell integrals;
      per-pair gate: the 7 cell integrals must resum to the full-circle one;
  (e) ROUTE 2 (rotation route, different pairing geometry): C_ell(s) =
      13 int_T U(x) V(x - s/13) 1_{frac(169 x) in win_ell}(x) dx  (proved by
      the substitution x = (y+u')/13 chartwise; then frac(13 y) = frac(169 x));
      independent full-table equality on all 91 entries + 13 globals;
  (f) gates: C_ell(0) = 0 for every ell; all C_ell(s) >= 0; row sums;
      sum_s C(s) = 169 I_3;  J(0) = I_3;  sum_k J(k) = 0;  J(k) != 0 for all
      k != 0 (THM-2471 (9)-(10));  Hermitian symmetry W(13-k,ell) =
      conj(W(k,ell)) exactly.

SCOPE (MISTAKE-281/283 discipline, hard):
  * the row is TYPED (THM-2309 (25)); NOT an asserted scalar cover;
  * no physical current is claimed, no row exclusion, no LRC(14)/JC(2)
    progress beyond what THM-2471/2560/2449/2512 literally license here;
  * W is an exact pairing of the r = 3 collision colours against the lawful
    owner-clock probe on one common base; it is NOT (yet) the THM-2512 d_A
    host, and no transplant theorem is claimed -- only the three inputs
    (T1)-(T3) are tested on this array;
  * all decisions are integer/Fraction/Z[zeta_13]-exact; floats only in ~.

Script: 04-computation/lrc14_base_only_bridge_opus_20260728.py
Output: 05-knowledge/results/lrc14_base_only_bridge_opus_20260728.out
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
# Row data and grid (identical to the audited 2026-07-27 scripts)
# --------------------------------------------------------------------------
W = (1, 14, 27, 40, 53, 66, 13, 13**3, 2 * 13**5)
GUARD, OWNER, TA, TB = 0, 6, 7, 8
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

DCOLL = 13**3                      # d = c_1 * 13^(r-1), r = 3 (THM-2560 (A))
I3 = Fraction(9926558757352, 109707098520974955)   # THM-2560 triple-verified


# --------------------------------------------------------------------------
# Interval machinery (verbatim from the audited scripts)
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


PAT_E = {0: "gout", 1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
         6: "in", 7: "out", 8: "out"}
PAT_QA = {0: "gout", 1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
          6: "out", 7: "in", 8: "out"}
ZELL = (0,) * 9


# --------------------------------------------------------------------------
# Generic exact profile machinery (grid-parameterized; shared by toy + main)
# A "profile" is (starts, vals): starts[0] = 0, piecewise constant integer
# values on [starts[i], starts[i+1]) covering [0, grid).
# --------------------------------------------------------------------------
def weighted_fold(pieces, mult, grid):
    """Profile of N(y) = sum_pieces w * #{preimages of y under x->mult*x mod 1
    lying in [a/grid, b/grid)} = sum_{b'<mult} h((y+b')/mult) numerators,
    i.e. the exact numerator profile of mult * P_mult h for h = step function
    with integer numerator w on each piece.  Mass conservation asserted."""
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
    """Profile of y -> h((y+u)/p) on [0,grid): the chart-window pullback
    THM-2471 (4)-(5).  Window u of h is [u*grid/p, (u+1)*grid/p);
    y-numerator = p*X - u*grid."""
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


def rotate_profile(starts, vals, sh, grid):
    """Profile of x -> h(x - sh/grid): piece [a,b) moves to [a+sh, b+sh)."""
    n = len(starts)
    pieces = []
    for i in range(n):
        a = starts[i]
        b = starts[i + 1] if i + 1 < n else grid
        a2 = (a + sh) % grid
        b2 = a2 + (b - a)
        if b2 <= grid:
            pieces.append((a2, b2, vals[i]))
        else:
            pieces.append((a2, grid, vals[i]))
            pieces.append((0, b2 - grid, vals[i]))
    pieces.sort()
    ns = []
    nv = []
    for a, b, v in pieces:
        ns.append(a)
        nv.append(v)
    require(ns and ns[0] == 0, "rotated profile must start at 0")
    return ns, nv


def product_cum(s1, v1, s2, v2, grid):
    """Merged product profile of two profiles + cumulative integral arrays.
    Returns (ps, pv, pc, ptot): value pv[i] on [ps[i], ps[i+1]),
    pc[i] = integral numerator up to ps[i], ptot = full-circle numerator."""
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


def clock_cells(p, q, grid, depthmap):
    """q interval lists: cell_ell = {y : frac(depthmap*y) in
    [(2ell-1)/(2q), (2ell+1)/(2q)) mod 1}, half-open, on the grid."""
    require(grid % (2 * q * depthmap) == 0, "grid does not resolve the cells")
    step = grid // depthmap
    u2q = grid // (2 * q * depthmap)
    cells = []
    for ell in range(q):
        pieces = []
        for n in range(depthmap):
            a = ((2 * ell - 1) * u2q + n * step) % grid
            b = a + 2 * u2q
            if b <= grid:
                pieces.append((a, b))
            else:
                pieces.append((a, grid))
                pieces.append((0, b - grid))
        pieces.sort()
        require(sum(bb - aa for aa, bb in pieces) == grid // q,
                "cell measure must be 1/q")
        cells.append(pieces)
    return cells


def check_partition(lists, grid):
    allp = sorted(iv for L in lists for iv in L)
    cur = 0
    for a, b in allp:
        require(a == cur and b > a, "cells do not tile the circle")
        cur = b
    require(cur == grid, "cell tiling incomplete")


# --------------------------------------------------------------------------
# The bridge engine (both routes + internal gates); used by toy AND main
# --------------------------------------------------------------------------
def bridge_tables(E_iv, Q_iv, p, qmod, Rpack, depth, grid, verbose=False):
    """Exact p x qmod table  Cnum[s][ell] = DEN * C_ell(s),
    DEN = Rpack * depth^2 * grid, where
      f = 1_Q * P_Rpack 1_E,  U = P_depth f,  V = P_depth 1_E,
      A(y,u) = U((y+u)/p),    F(y,u) = V((y+u)/p),
      C_ell(s) = int_{cell_ell} sum_u A(y,u+s) F(y,u) dy,
      cell_ell = {y : frac(p*y) in [(2ell-1)/(2qmod),(2ell+1)/(2qmod)) mod 1}.
    Route 1 (windows+cells) fills the table; per-(u0,s) pair the 7 cell
    integrals are asserted to resum to the full-circle integral.
    Route 2 (rotation+pullback cells) recomputes every entry independently:
      C_ell(s) = p * int U(x) V(x - s/p) 1_{frac(p^2 x) in win_ell} dx,
    and the two tables are asserted equal.  Returns dict."""
    # f = 1_Q * P_Rpack 1_E  (integer numerators over Rpack)
    n_st, n_v = weighted_fold([(a, b, 1) for a, b in E_iv], Rpack, grid)
    fp = []
    np_ = len(n_st)
    for qa, qb in Q_iv:
        i = bisect_right(n_st, qa) - 1
        while True:
            pa = n_st[i]
            pb = n_st[i + 1] if i + 1 < np_ else grid
            a = qa if qa > pa else pa
            b = qb if qb < pb else pb
            if a < b and n_v[i]:
                fp.append((a, b, n_v[i]))
            if pb >= qb:
                break
            i += 1
    int_f_num = sum(w * (b - a) for a, b, w in fp)     # over Rpack*grid
    Ust, Uv = weighted_fold(fp, depth, grid)           # over Rpack*depth
    Vst, Vv = weighted_fold([(a, b, 1) for a, b in E_iv], depth, grid)
    # P_m preserves integrals:
    require(profile_mass(Ust, Uv, grid) == depth * int_f_num,
            "int U != int f")
    lenE = sum(b - a for a, b in E_iv)
    require(profile_mass(Vst, Vv, grid) == depth * lenE, "int V != mu(E)")
    if verbose:
        log(f"    f: {len(fp)} weighted pieces, int f = "
            f"{Fraction(int_f_num, Rpack * grid)}")
        log(f"    U = P_{depth} f: {len(Ust)} pieces (num/den "
            f"{Rpack * depth});  V = P_{depth} 1_E: {len(Vst)} pieces "
            f"(num/den {depth})")
    cells = clock_cells(p, qmod, grid, p)
    check_partition(cells, grid)
    pulls = clock_cells(p, qmod, grid, p * p)
    check_partition(pulls, grid)
    if verbose:
        log(f"    owner-clock cells: {qmod} cells x {p} intervals, exact"
            f" tiling: PASS;  pullback cells: {qmod} x {p * p}, exact"
            f" tiling: PASS")
    # ---- ROUTE 1: window route
    Aprof = [extract_window(Ust, Uv, u, p, grid) for u in range(p)]
    Fprof = [extract_window(Vst, Vv, u, p, grid) for u in range(p)]
    Cnum = [[0] * qmod for _ in range(p)]
    Cglob = [0] * p
    for u0 in range(p):
        s2, v2 = Fprof[u0]
        for s in range(p):
            u1 = (u0 + s) % p
            s1, v1 = Aprof[u1]
            ps, pv, pc, ptot = product_cum(s1, v1, s2, v2, grid)
            Cglob[s] += ptot
            tot_cells = 0
            for ell in range(qmod):
                acc = 0
                for a, b in cells[ell]:
                    acc += pquery(ps, pv, pc, b) - pquery(ps, pv, pc, a)
                Cnum[s][ell] += acc
                tot_cells += acc
            require(tot_cells == ptot,
                    "cell integrals do not resum to the full circle")
    if verbose:
        log(f"    {elapsed()} route 1 (windows+cells) done; all"
            f" {p * p} per-pair cell-partition gates: PASS")
    # ---- ROUTE 2: rotation route
    T_p = grid // p
    C2 = [[0] * qmod for _ in range(p)]
    for s in range(p):
        rs, rv = rotate_profile(Vst, Vv, s * T_p, grid)   # V(x - s/p)
        ps, pv, pc, ptot = product_cum(Ust, Uv, rs, rv, grid)
        require(p * ptot == Cglob[s], "route-2 full-circle mismatch")
        for ell in range(qmod):
            acc = 0
            for a, b in pulls[ell]:
                acc += pquery(ps, pv, pc, b) - pquery(ps, pv, pc, a)
            C2[s][ell] = p * acc
    require(C2 == Cnum, "route-2 table != route-1 table")
    if verbose:
        log(f"    {elapsed()} route 2 (rotation+pullback) done: full"
            f" {p}x{qmod} table and all {p} globals AGREE with route 1")
    require(all(x >= 0 for row in Cnum for x in row),
            "negative C entry (impossible for nonneg step functions)")
    for s in range(p):
        require(sum(Cnum[s]) == Cglob[s], "row sum != global")
    return {"C": Cnum, "Cglob": Cglob, "DEN": Rpack * depth * depth * grid,
            "int_f_num": int_f_num, "nU": len(Ust), "nV": len(Vst),
            "nf": len(fp)}


# --------------------------------------------------------------------------
# Toy brute-force self-test of the ENTIRE pipeline (route (a))
# --------------------------------------------------------------------------
def toy_test():
    """Hostile little model: grid 90, base p = 3, packet clock R' = 9,
    collision depth d' = 27, owner-clock cells mod 5, sets with wraps and
    unequal pieces.  The engine table (both routes) must equal the
    definitional midpoint-evaluation Fraction table."""
    grid, p, qmod, Rp, depth = 90, 3, 5, 9, 27
    Ep = [(0, 3), (17, 21), (34, 44), (60, 61), (85, 90)]
    Qp = [(4, 14), (40, 63), (77, 88)]
    res = bridge_tables(Ep, Qp, p, qmod, Rp, depth, grid)
    DEN = res["DEN"]

    def mk_ind(iv):
        return lambda x: 1 if any(Fraction(a, grid) <= x < Fraction(b, grid)
                                  for a, b in iv) else 0

    indE = mk_ind(Ep)
    indQ = mk_ind(Qp)
    fmemo = {}

    def f_at(x):
        r = fmemo.get(x)
        if r is None:
            if not indQ(x):
                r = Fraction(0)
            else:
                r = Fraction(sum(indE(((x + c) / Rp) % 1)
                                 for c in range(Rp)), Rp)
            fmemo[x] = r
        return r

    def U_at(x):
        return sum(f_at(((x + b) / depth) % 1)
                   for b in range(depth)) / depth

    def V_at(x):
        return Fraction(sum(indE(((x + b) / depth) % 1)
                            for b in range(depth)), depth)

    Cb = [[Fraction(0)] * qmod for _ in range(p)]
    for t in range(grid):
        y = Fraction(2 * t + 1, 2 * grid)
        tt = (p * y) % 1
        ell = int(tt * qmod + Fraction(1, 2)) % qmod   # nearest ell/qmod
        Au = [U_at(((y + u) / p) % 1) for u in range(p)]
        Fu = [V_at(((y + u) / p) % 1) for u in range(p)]
        for s in range(p):
            for u in range(p):
                Cb[s][ell] += Au[(u + s) % p] * Fu[u]
    ok_nonzero = 0
    for s in range(p):
        for ell in range(qmod):
            require(Cb[s][ell] / grid == Fraction(res["C"][s][ell], DEN),
                    f"toy mismatch at (s,ell)=({s},{ell})")
            if res["C"][s][ell]:
                ok_nonzero += 1
    require(ok_nonzero > 0, "toy table degenerate (test has no content)")
    log(f"    toy pipeline (grid 90, p=3, R'=9, d'=27, cells mod 5):"
        f" all {p}x{qmod} entries == definitional brute force"
        f" ({ok_nonzero} nonzero): PASS")


# --------------------------------------------------------------------------
# Exact Z[zeta_13] layer: 13-bucket integer vectors, canonical 12-dim form
# --------------------------------------------------------------------------
def canon(b):
    t = b[12]
    return tuple(x - t for x in b[:12])


def is_zero_b(b):
    return all(x == 0 for x in canon(b))


def badd(a, c):
    return [x + y for x, y in zip(a, c)]


def bsub(a, c):
    return [x - y for x, y in zip(a, c)]


def bscale(a, m):
    return [m * x for x in a]


def bmul(a, c):
    r = [0] * 13
    for i, ai in enumerate(a):
        if ai:
            for j, cj in enumerate(c):
                if cj:
                    r[(i + j) % 13] += ai * cj
    return r


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


# ==========================================================================
def main():
    log("THE BASE-ONLY BRIDGE PAIRING on the canonical typed row"
        " (THM-2560 redirect (A))")
    log("script=04-computation/lrc14_base_only_bridge_opus_20260728.py")
    log("machine=opus  date=2026-07-28")
    log("")
    log("[1] object and logged conventions")
    log(f"    row w={W}  owner j=1 (c_1={C1}), word sigma={{a}},"
        f" packet clock R=13^2=169")
    log(f"    r = 3 (THM-2560 (A)), d = 13^3 = {DCOLL}:  U = P_d f,"
        f"  V = P_d e")
    log("    e = 1_(E_1);  f = 1_(Q_1,{a}) P^2 1_(E_1)   (THM-2471 (23))")
    log("    A(y,u) = U((y+u)/13), F(y,u) = V((y+u)/13)  (THM-2471 (4)-(5))")
    log("    alpha_k, phi_k with 1/13 normalization and zeta^(-ku) sign"
        " (THM-2471 (8))")
    log("    colour density  j_k(y) = alpha_k(y) conj(phi_k(y))"
        " = (1/169) sum_s c_s(y) zeta^(-ks),")
    log("      c_s(y) = sum_u A(y,u+s) F(y,u)      (THM-2471 (13) integrand)")
    log("    owner-clock probe (THM-2449 (14), j=1): cell_ell ="
        " {y : d(13y - ell/7)=1}")
    log("      realized half-open: frac(13y) in [(2ell-1)/14,(2ell+1)/14)"
        " mod 1  (a.e. equal)")
    log("    W(k,ell) = int_T j_k(y) 1_(cell_ell)(y) dy"
        "  -- ONE integral over the common base y;")
    log("    the ell-restriction sits INSIDE the y-integral, BEFORE the"
        " root-DFT and BEFORE")
    log("    any marginalization (MISTAKE-281/283 common-base discipline).")
    log("    disintegration: C_ell(s) = int_(cell_ell) c_s;"
        "  W(k,ell) = (1/169) sum_s C_ell(s) zeta^(-ks);")
    log("    J(k) = sum_ell W(k,ell) = (1/169) sum_s C(s) zeta^(-ks)"
        "   (THM-2471 (8),(13),(15))")
    log("    zeta = abstract primitive 13th root; basis 1..zeta^11,"
        " zeta^12 = -(1+...+zeta^11)")
    log(f"    grid: half-open [A,B) on T_DEN = 182*lcm(w) = {T_DEN};"
        " no floats in decisions")
    log("")

    log("[2] machinery self-test (hostile toy model, definitional brute"
        " force)")
    toy_test()
    log("")

    log("[3] exact interval sets and anchors")
    E = build_set(PAT_E, ZELL)
    lenE = check_intervals(E, T_DEN)
    require(Fraction(lenE, T_DEN) == Fraction(1882176, 28589561),
            "mu(E_1) anchor")
    require(len(E) == 57072, "E_1 interval count anchor")
    log(f"    E_1: intervals={len(E)}  measure={Fraction(lenE, T_DEN)}: PASS")
    Q = build_set(PAT_QA, ZELL)
    lenQ = check_intervals(Q, T_DEN)
    require(Fraction(lenQ, T_DEN) == Fraction(143103830843, 5727632650740),
            "mu(Q_1,{a}) anchor")
    require(len(Q) == 22478, "Q_1,{a} interval count anchor")
    log(f"    Q_1,{{a}}: intervals={len(Q)}"
        f"  measure={Fraction(lenQ, T_DEN)}: PASS")
    log("")

    log("[4] profiles, cells, and the exact 13 x 7 table (both routes)")
    res = bridge_tables(E, Q, 13, 7, 169, DCOLL, T_DEN, verbose=True)
    Cnum = res["C"]
    Cglob = res["Cglob"]
    DENC = res["DEN"]
    require(DENC == 169 * DCOLL * DCOLL * T_DEN, "DENC shape")
    require(Fraction(res["int_f_num"], 169 * T_DEN)
            == Fraction(21376087, 17907461390),
            "int f != stored mu(F) anchor")
    log("    int f == stored measure(E_1 n T^-2 Q_1,{a}) ="
        " 21376087/17907461390: PASS")
    log("    int U == int f and int V == mu(E_1)  (P_m preserves"
        " integrals): PASS")
    log("")

    log("[5] gates on the raw table  (DENC = 169 * 2197^2 * T_DEN ="
        f" {DENC})")
    log(f"    DENC factored: {factor_str(DENC)}")
    # C_ell(0) = 0 cellwise (first-collision disjointness survives cells)
    require(all(x == 0 for x in Cnum[0]),
            "C_ell(0) != 0: UV = 0 a.e. violated cellwise")
    log("    C_ell(0) = 0 for every ell  (UV = 0 a.e. at the first"
        " collision, cellwise): PASS")
    require(Cglob[0] == 0, "C(0) != 0")
    tot = sum(Cglob)
    require(Fraction(tot, DENC) == 169 * I3,
            "sum_s C(s) != 169 * I_3 (THM-2471 (14) / THM-2560 anchor)")
    log("    sum_s C(s) = 169*I_3 = 9926558757352/649154429118195"
        "  (I_3 of THM-2560): PASS")
    log("    all 91 entries nonnegative rationals: PASS")
    log("")
    log("    The exact table DENC*C_ell(s)  (rows s=0..12, columns"
        " ell=0..6):")
    for s in range(13):
        log(f"      s={s:2d}: {Cnum[s]}")
    log("    row sums (= DENC*C(s)):")
    for s in range(13):
        log(f"      s={s:2d}: {Cglob[s]}")
    log("")

    log("[6] the colour table W(k,ell) in Z[zeta_13]/(169*DENC) and its"
        " marginals")
    # integer bucket vectors over DEN_W = 169*DENC
    DEN_W = 169 * DENC
    WI = [[None] * 7 for _ in range(13)]
    for k in range(13):
        for ell in range(7):
            b = [0] * 13
            for s in range(13):
                b[(-k * s) % 13] += Cnum[s][ell]
            WI[k][ell] = b
    JI = []
    for k in range(13):
        b = [0] * 13
        for s in range(13):
            b[(-k * s) % 13] += Cglob[s]
        JI.append(b)
    # GATE: ell-sum (the k-indexed marginal) reproduces J(k)
    for k in range(13):
        acc = [0] * 13
        for ell in range(7):
            acc = badd(acc, WI[k][ell])
        require(canon(acc) == canon(JI[k]),
                f"sum_ell W(k,ell) != J(k) at k={k}")
    log("    GATE  sum_ell W(k,ell) == J(k) for every k in F_13: PASS")
    # row-zero along k
    for ell in range(7):
        acc = [0] * 13
        for k in range(13):
            acc = badd(acc, WI[k][ell])
        require(is_zero_b(acc), f"sum_k W(k,ell) != 0 at ell={ell}")
    log("    LAW   sum_k  W(k,ell) == 0 exactly for every ell in F_7"
        "  (= C_ell(0)/13): PASS")
    # J gates
    require(canon(JI[0]) == canon([tot] + [0] * 12), "J(0) bucket shape")
    require(Fraction(tot, DEN_W) == I3, "J(0) != I_3")
    log("    J(0) = I_3 = 9926558757352/109707098520974955: PASS"
        "   (GATE, triple-verified value)")
    accJ = [0] * 13
    for k in range(13):
        accJ = badd(accJ, JI[k])
    require(is_zero_b(accJ), "sum_k J(k) != 0")
    log("    sum_k J(k) = 0, hence sum_(k!=0) J(k) = -I_3"
        "  (THM-2471 (10)): PASS")
    nzJ = sum(1 for k in range(1, 13) if not is_zero_b(JI[k]))
    require(nzJ == 12, "some nonzero colour J(k) vanishes")
    log("    J(k) != 0 for all 12 k != 0  (THM-2471 (9)): PASS")
    for k in range(13):
        for ell in range(7):
            require(canon(WI[(13 - k) % 13][ell])
                    == canon(bconj(WI[k][ell])),
                    "Hermitian symmetry broken")
    log("    Hermitian symmetry W(13-k,ell) = conj(W(k,ell)) exact: PASS")
    log("")
    log(f"    J(k) exactly (canonical coords n_0..n_11 over DEN_J = 169*DENC"
        f" = {DEN_W};")
    log("     J(k) = (n_0 + n_1 zeta + ... + n_11 zeta^11)/DEN_J):")
    for k in range(13):
        log(f"      k={k:2d}: {list(canon(JI[k]))}")
    log("")
    log("    W(k,ell) exactly (canonical coords over the same DEN_J);"
        " k=0 column is rational:")
    for ell in range(7):
        log(f"      W(0,ell={ell}) = {Fraction(canon(WI[0][ell])[0], DEN_W)}")
    log("    full mixed entries:")
    for k in range(1, 13):
        for ell in range(7):
            log(f"      k={k:2d} ell={ell}: {list(canon(WI[k][ell]))}")
    log("")

    log("[7] mixedness decision (all tests exact over Q(zeta_13))")
    ell_dep = any(canon(WI[k][ell]) != canon(WI[k][0])
                  for k in range(13) for ell in range(1, 7))
    k_dep = any(canon(WI[k][ell]) != canon(WI[0][ell])
                for ell in range(7) for k in range(1, 13))
    log(f"    ell-dependence (exists k: W(k,.) nonconstant): "
        f"{'YES' if ell_dep else 'NO'}")
    log(f"    k-dependence   (exists ell: W(.,ell) nonconstant): "
        f"{'YES' if k_dep else 'NO'}")
    nz_minor = 0
    minors_total = 0
    first_nz = None
    for k1 in range(13):
        for k2 in range(k1 + 1, 13):
            for e1 in range(7):
                for e2 in range(e1 + 1, 7):
                    minors_total += 1
                    m = bsub(bmul(WI[k1][e1], WI[k2][e2]),
                             bmul(WI[k1][e2], WI[k2][e1]))
                    if not is_zero_b(m):
                        nz_minor += 1
                        if first_nz is None:
                            first_nz = (k1, k2, e1, e2)
    require(minors_total == 1638, "minor count")
    log(f"    2x2 minors over Q(zeta_13): {nz_minor} nonzero of"
        f" {minors_total}")
    if nz_minor:
        log(f"    first nonzero minor at (k1,k2;ell1,ell2) = {first_nz}"
            "  => rank >= 2: W is NOT a")
        log("    product a(k)b(ell) in any field factorization, and not"
            " ell-independent.")
        verdict_mixed = True
    else:
        verdict_mixed = False
        log("    all minors vanish => rank <= 1: W is a product/degenerate;"
            " see shape below.")
    if verdict_mixed:
        log("    VERDICT: W is GENUINELY (k,ell)-MIXED.")
    else:
        if not ell_dep:
            log("    VERDICT: DEGENERATE -- W(k,ell) = J(k)/7 exactly"
                " (ell-independent).")
        else:
            log("    VERDICT: DEGENERATE -- rank-1 product, ell-dependent.")
    log("")

    log("[8] the three transplant inputs (HYP-9032 Sec. 2) on the centred"
        " array")
    log("    Wc(k,ell) := W(k,ell) - J(k)/7   (the ell-centring; integer"
        " buckets over 7*169*DENC)")
    DEN_C = 7 * DEN_W
    WcI = [[bsub(bscale(WI[k][ell], 7), JI[k]) for ell in range(7)]
           for k in range(13)]
    # (T1)
    for k in range(13):
        acc = [0] * 13
        for ell in range(7):
            acc = badd(acc, WcI[k][ell])
        require(is_zero_b(acc), "T1 row-zero fails (impossible)")
    for ell in range(7):
        acc = [0] * 13
        for k in range(13):
            acc = badd(acc, WcI[k][ell])
        require(is_zero_b(acc), "T1 column-zero fails")
    log("    (T1) row-zero over F_7:    sum_ell Wc(k,ell) = 0 for every k:"
        "  HOLDS (by centring)")
    log("         column-zero over F_13: sum_k Wc(k,ell) = 0 for every ell:"
        "  HOLDS -- because")
    log("         sum_k W(k,ell) = 0 (computed law, [6]) and sum_k J(k) = 0;")
    log("         this is the exact analogue of THM-2512 (4)+(5) for this"
        " array.")
    log("    NOTE  the RAW array W fails F_7-row-zero: sum_ell W(k,ell) ="
        " J(k) != 0 for k != 0;")
    log("          the law that holds RAW is along k (13-primary"
        " invisibility, HYP-9032 obstacle 3).")
    nz_c = sum(1 for k in range(13) for ell in range(7)
               if not is_zero_b(WcI[k][ell]))
    log(f"    Wc nonzero entries: {nz_c} of 91"
        f"  (Wc = 0 would mean W(k,ell) = J(k)/7 exactly)")
    # (T2)
    kc_dep = any(canon(WcI[k][ell]) != canon(WcI[0][ell])
                 for ell in range(7) for k in range(1, 13))
    log("    (T2) h-independence dichotomy (h = k, the F_13 index):"
        " STRUCTURAL --")
    log("         if Wc(k,ell) = b(ell) independent of k, column-zero"
        " gives 13 b(ell) = 0,")
    log("         hence Wc = 0  (same argument as THM-2512 (5)).  Instance:"
        f" Wc is"
        f" {'k-DEPENDENT' if kc_dep else 'k-independent'}"
        f" and {'NONZERO' if nz_c else 'ZERO'};")
    log("         the dichotomy is satisfied"
        + (" non-vacuously." if (kc_dep and nz_c) else " (degenerate"
           " branch)."))
    # (T3)
    from math import lcm as _lcm
    D_all = 1
    minabs = None
    for k in range(13):
        for ell in range(7):
            for cco in canon(WcI[k][ell]):
                if cco:
                    fr = Fraction(cco, DEN_C)
                    D_all = _lcm(D_all, fr.denominator)
                    a = abs(fr)
                    if minabs is None or a < minabs:
                        minabs = a
    log("    (T3) denominator floor (zeta-basis coordinate convention,"
        " logged):")
    log(f"         universal denominator D of all Wc coordinates:"
        f" D = {D_all}")
    log(f"         D factored: {factor_str(D_all)}")
    log(f"         min nonzero |coordinate| = {minabs}"
        f"  (>= 1/D by construction)")
    D_raw = 1
    minabs_raw = None
    for k in range(13):
        for ell in range(7):
            for cco in canon(WI[k][ell]):
                if cco:
                    fr = Fraction(cco, DEN_W)
                    D_raw = _lcm(D_raw, fr.denominator)
                    a = abs(fr)
                    if minabs_raw is None or a < minabs_raw:
                        minabs_raw = a
    log(f"         (raw W: universal denominator {D_raw} ="
        f" {factor_str(D_raw)}; min nonzero |coord| = {minabs_raw})")
    log("")

    log("[9] per-cell colour saturation (rationality argument,"
        " disintegrated)")
    dead = []
    for ell in range(7):
        nzk = sum(1 for k in range(1, 13) if not is_zero_b(WI[k][ell]))
        allzero = all(Cnum[s][ell] == 0 for s in range(13))
        log(f"    cell ell={ell}: nonzero mixed colours {nzk}/12;"
            f"  W(0,ell) {'>' if canon(WI[0][ell])[0] else '='} 0")
        if allzero:
            dead.append(ell)
        # rational C_ell + Galois: one zero at k!=0 => all zero => C_ell
        # constant in s; with C_ell(0)=0 => cell colour-dead:
        require(nzk in (0, 12), "per-cell all-or-nothing violated")
    log("    per cell: either all 12 nonzero colours fire or the cell is"
        " colour-dead")
    log(f"    colour-dead cells: {dead if dead else 'NONE'}")
    log("")

    log("[10] scope notes (MISTAKE-281/283 discipline)")
    log("    * TYPED row (THM-2309 (25)); NOT an asserted scalar cover; no"
        " physical")
    log("      current, no row exclusion; LRC(14) and JC(2) remain OPEN.")
    log("    * W is ONE pairing on ONE common base (the physical circle"
        " T): the")
    log("      collision-colour density and the owner-clock probe are"
        " evaluated at the")
    log("      SAME y inside a single integral; no product of separately")
    log("      marginalized controls occurs anywhere (MISTAKE-281 rule).")
    log("    * The array W is NOT the THM-2512 lawful-interaction defect"
        " d_A and is not")
    log("      claimed to be; (T1)-(T3) are tested as exact facts about"
        " THIS array on")
    log("      this row's r = 3 packet at the packet clock, nothing more.")
    log("    * The owner-clock label ell is carried, never summed, exactly"
        " as HYP-9032")
    log("      obstacle (3) requires; the raw k-sum law shows why summing"
        " would kill it.")
    log("    * All decisions are integer/Fraction/Z[zeta_13]-exact;"
        " floats never enter.")
    log("")
    log(f"{elapsed()} all checks passed")


if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    main()
