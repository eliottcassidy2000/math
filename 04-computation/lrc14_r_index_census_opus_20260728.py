#!/usr/bin/env python3
r"""r-index census across the lawful (j, sigma, K) configuration space of the
canonical typed row THM-2309 (25) -- the THM-2471 (23)-(25) collision index.

TASK (baton 3, 2026-07-28).  THM-2471 (23)-(24) define, from the THM-2306
data on the canonical typed row

    w = (H,q1..q5,c1,c2,c3) = (1,14,27,40,53,66,13,2197,742586),

    e = 1_{E_j},   f = 1_{Q_{j,sigma}} P^K 1_{E_j},   c = c_j,
    A = P_c f,     B = P_c e,
    I_s = int_T (P^s A)(P^s B),      r = min{s >= 1 : I_s > 0}.

Stage 1 (04-computation/lrc14_ancestry_realization_stage_opus_20260727.py,
machinery reused verbatim here) decided r = 3 at the single configuration
(j=1, sigma={a}, K=2).  This script censuses r across every configuration
the files license and reports whether the r = 5 window (the unique affine
invariant theta = t - 2u of the THM-2471 (44)-(47) strict-profile table)
opens at ANY lawful configuration.

LICENSING LEDGER (each item cited; checked exactly in the run):
  owners:
    THM-2309 (25) types the row with owner j = the depth-one blocker,
    (w_j, w_a, w_b) = (13, 13^3, 2*13^5).  THM-2349 (2)/(23a) licenses a
    SECOND shallow owner only on repeated-first rows (nu13(c_2) = 1); this
    row has strict profile nu13(c1,c2,c3) = (1,3,5), so j = 2 is NOT
    lawful.  Diagnostic: E_2 = Q_{1,{a}} (THM-2305 (3) literally) has
    positive measure, but measure positivity is not a license -- the
    THM-2349 second-owner license is profile-typed.  Census owner: j = 1.
  word strata (THM-2305 fixed order):
    THM-2305 (3)-(4): sigma in {{a},{b},{a,b}} with a = c_2, b = c_3.
    Under the scalar-cover reading these partition R_1 = A_0 \ D_c1; the
    row is TYPED (MISTAKE-281), so the residual EMPTY-WORD set
    A_0 \ (D_c1 u D_c2 u D_c3) may carry positive measure.  It is NOT a
    THM-2305 canonical terminal word, and THM-2306 Sec. 1 / THM-2471 (23)
    restrict Q to the three canonical words; the empty word is therefore
    NOT lawful for the census.  It is computed as an explicitly flagged
    DIAGNOSTIC only (it does satisfy the bare (41a) premises).
  clocks:
    THM-2471 (23) prescribes K = lambda+1 = 2 (packet clock R = 13^K).
    THM-2306 Sec. 6 item 4: the disjoint-support proof uses only
      (41a)  E subset D_c,   Q subset D_c^c,   int_Q P^K 1_E > 0,
    at ANY fixed clock; THM-2471 Sec. 3: "The same proof works at any
    fixed positive-return clock allowed by THM-2306 Section 6."  Hence
    K = 3, 4 are lawful exactly when rho^(K) = int_Q P^K 1_E > 0,
    computed exactly per configuration below.  (THM-2349 Sec. 3 licenses
    delayed clocks k >= 2 by BV mixing as a sufficient condition; here
    the exact computed positivity is the license actually used.)

CENSUS SPACE: (j, sigma, K) in {1} x {{a},{b},{a,b}} x {2,3,4}, nine
lawful candidates, each contingent on its exactly computed rho^(K) > 0;
plus DIAGNOSTIC rows (empty word at K = 2,3,4, and the aggregate R_1).

REDUCTION (proved in the Stage-1 docstring; self-tested below).  With
m = 13^(s+1), R = 13^K, the transfer-adjoint chain gives

    I_s = int_T 1_E(u) 1_Q({R u}) Phi_s({R u}) du,
    Phi_s(x) = (1/m) sum_{j<m} 1_E(x + j/m),

and, by one adjoint step in the other order (the Stage-1 route (f),
audited there against the packet-walk engine),

    I_s = int_T 1_Q(x) (N_R(x)/R) Phi_s(x) dx,
    N_R(x) = sum_{b<R} 1_E((x+b)/R)  (the multiplicity fold of R*E).

The y-side form is the MAIN engine here: at K = 4 the packet clock
R = 28561 breaks the Stage-1 single-wrap F-walk precondition
(R * maxlen(E) < T_DEN), while N_R is a 2|E|-event fold at any R.

CONVENTIONS (each an explicit logged choice; files override summaries):
  * grid: unions of half-open [A,B) with integer endpoints on the
    T_DEN = 182*lcm(w) grid; Phi_s profiles on g*T_DEN, g = 13^max(0,s-5);
    no floats enter any decision;
  * P_a h(y) = (1/a) sum_{b<a} h((y+b)/a)  (THM-2471 (2));
  * E_1, Q_1,{a} built exactly as in the audited Stage-1/owner-loop
    scripts (anchor-checked); Q_1,{b}, Q_1,{a,b}, the empty-word residual
    and R_1 built by the same comb machinery;
  * the full word Q_{j,sigma} of THM-2471 (23) is used (source-safe factor
    RETAINED; no THM-2449 source deletion);
  * danger d(y) = 1_{||y|| < 1/14}; windows via the 182/91-denominator
    combs of the audited scripts;
  * K is the PACKET clock of THM-2471 (23)/(28) (ancestry clock R = 13^K),
    not a response clock.

VERIFICATION MATRIX (independent routes, all exact):
  (a) toy brute-force self-test on a hostile little model (grid 405,
      base p = 3): direct rational computation of int (P_m f)(P_m e)
      == the y-side N/Phi reduction, for ALL packet clocks
      R' in {3,9,27} x m in {3,9,27}, including wrap cases (R'*len > 1);
      the toy N-profile also cross-checks the production count_fold;
  (b) anchors: mu(E_1), mu(Q_1,{a}), interval counts, rho^(2)(Q_1,{a}),
      and the full Stage-1 I_s sequence and r = 3 at (j=1,{a},K=2)
      against the audited Stage-1 artifact;
  (c) structural owner-normalized disjointness: E_1 subset D_c1 and
      Q n D_c1 = empty interval-wise for every stratum, forcing I_0 = 0
      at every K (THM-2306 normalization premise);
  (d) support-fold route (NO shared machinery with the prefix engine):
      supp P_m h = union fold; I_s = 0 iff the supports of P_m f and
      P_m e share no positive measure -- decides zero/positive for EVERY
      (config, s) independently of the exact-value engine;
  (e) correlation route for s = 0,1 (replica-dichotomy prefix identity):
      I_s = (1/m) sum_j I_R(E, Q n (E - j/m)), exact I_R prefix sweeps,
      run for all four strata at all three clocks;
  (f) packet-walk route (the Stage-1 MAIN engine, generalized to R):
      F = E n (R.)^-1 Q on the R*T_DEN grid, I_s = int_F Phi_s({Ru}) du,
      run wherever the single-wrap precondition holds AND the expected
      F interval count ~ |Q|*R*mu(E) stays below a 5e6 materialization
      cap (the route is a redundant 4th check), at s in {r-1, r};
      also mu(F) == rho there;
  (g) full-circle Phi controls at every (K,s), and the exact partition
      additivity control I_s(R_1) = sum_parts I_s at every computed s.

SCOPE (MISTAKE-281/283 discipline, hard):
  * the row is TYPED (THM-2309 (25)); NOT an asserted scalar cover;
  * no physical current is claimed, no row exclusion, no LRC(14)/JC(2)
    progress beyond what THM-2471/THM-2306/THM-2305/THM-2309/THM-2349
    literally license on this one row;
  * every I_s is one integral of a triple product at one common base
    point (paired on one base BEFORE any marginalization);
  * the deep-root sidecar branches are THM-2471 (44)-(47) arithmetic
    applied to the exact computed r; no intertwiner is constructed and
    none is claimed.

Script: 04-computation/lrc14_r_index_census_opus_20260728.py
Output: 05-knowledge/results/lrc14_r_index_census_opus_20260728.out
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
# Row data and grid (identical to the audited Stage-1 script)
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
K_CANON = LAMBDA + 1
require(K_CANON == 2, "canonical packet clock K = lambda+1 = 2")

LCM_W = 1
for v in W:
    LCM_W = LCM_W * v // gcd(LCM_W, v)
T_DEN = 182 * LCM_W
require(T_DEN == 297836897838480, "T_DEN anchor")
require(nu13(T_DEN) == 6, "13-adic valuation of T_DEN must be 6")


# --------------------------------------------------------------------------
# Interval machinery (verbatim from the audited Stage-1 script; build_set
# generalized ONLY by the no-'in' branch needed for the all-safe residual)
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
    """pattern: dict index -> 'in' | 'out' | 'gout'; ell: 9-tuple of ints.
    Generalization vs Stage-1: with no 'in' label, start from the full
    circle (needed for R_1 and the all-safe empty-word residual)."""
    ins = [i for i, mmode in pattern.items() if mmode == "in"]
    if ins:
        start = min(ins, key=lambda i: W[i])
        iv = in_comb(start, ell[start])
    else:
        start = None
        iv = [(0, T_DEN)]
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
PAT_QB = {0: "gout", 1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
          6: "out", 7: "out", 8: "in"}
PAT_QAB = {0: "gout", 1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
           6: "out", 7: "in", 8: "in"}
PAT_QE = {0: "gout", 1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
          6: "out", 7: "out", 8: "out"}
PAT_R1 = {0: "gout", 1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
          6: "out"}
ZELL = (0,) * 9


# --------------------------------------------------------------------------
# Generic exact helpers (verbatim from Stage-1)
# --------------------------------------------------------------------------
def rotate(iv, sh, grid):
    """The set - sh/grid, i.e. intervals [A-sh, B-sh) mod grid, sorted."""
    out = []
    for a, b in iv:
        a2 = (a - sh) % grid
        b2 = a2 + (b - a)
        if b2 <= grid:
            out.append((a2, b2))
        else:
            out.append((a2, grid))
            out.append((0, b2 - grid))
    out.sort()
    return out


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


def measure(iv):
    return sum(b - a for a, b in iv)


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


def union_fold(iv, mult, grid):
    """Union set {mult*x mod 1 : x in iv} on [0,grid), sorted merged.
    This is exactly supp(P_mult h) for h with supp h = iv (nonneg h)."""
    pieces = []
    for a, b in iv:
        L = mult * (b - a)
        if L >= grid:
            return [(0, grid)]
        a2 = (mult * a) % grid
        b2 = a2 + L
        if b2 <= grid:
            pieces.append((a2, b2))
        else:
            pieces.append((a2, grid))
            pieces.append((0, b2 - grid))
    pieces.sort()
    out = []
    for a, b in pieces:
        if out and a <= out[-1][1]:
            if b > out[-1][1]:
                out[-1] = (out[-1][0], b)
        else:
            out.append((a, b))
    return out


def count_fold(iv, mult, grid):
    """Multiplicity fold N(y) = sum_{k in Z} 1_{mult*[A,B)}(y + k*grid) on
    [0,grid), as a sorted piecewise-constant profile.
    Returns (starts, vals, cums, total)."""
    const = 0
    events = []
    ap = events.append
    for a, b in iv:
        L = mult * (b - a)
        const += L // grid
        rem = L % grid
        if rem:
            st = (mult * a) % grid
            en = st + rem
            if en <= grid:
                ap((st, 1))
                ap((en, -1))
            else:
                ap((st, 1))
                ap((grid, -1))
                ap((0, 1))
                ap((en - grid, -1))
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
    require(all(v >= 0 for v in vals), "negative fold multiplicity")
    cums = [0] * len(starts)
    for i in range(1, len(starts)):
        cums[i] = cums[i - 1] + vals[i - 1] * (starts[i] - starts[i - 1])
    total = cums[-1] + vals[-1] * (grid - starts[-1])
    require(total == mult * measure(iv), "fold mass not conserved")
    return starts, vals, cums, total


# --------------------------------------------------------------------------
# Exact prefix identity of the replica-dichotomy script (routes (e), rho)
# --------------------------------------------------------------------------
def make_prefix(iv):
    starts = [a for a, _ in iv]
    lens = [b - a for a, b in iv]
    pref = [0] * (len(iv) + 1)
    for i, ln in enumerate(lens):
        pref[i + 1] = pref[i] + ln
    return starts, lens, pref


def sweep_acc(iv, Rred, starts, lens, pref):
    """Return (acc_r, acc_phi) = sums of (rB-rA) and (Phi(rB)-Phi(rA))."""
    br = bisect_right
    Tl = T_DEN
    acc_r = 0
    acc_p = 0
    for a, b in iv:
        ra = a * Rred % Tl
        rb = b * Rred % Tl
        acc_r += rb - ra
        i = br(starts, ra) - 1
        if i >= 0:
            d = ra - starts[i]
            L = lens[i]
            pa = pref[i] + (d if d < L else L)
        else:
            pa = 0
        i = br(starts, rb) - 1
        if i >= 0:
            d = rb - starts[i]
            L = lens[i]
            pb = pref[i] + (d if d < L else L)
        else:
            pb = 0
        acc_p += pb - pa
    return acc_r, acc_p


def IR_exact(E_iv, lenE_, Q_iv, R):
    """I_R(E,Q) = int 1_E(x) 1_Q({Rx}) dx, exact."""
    if not Q_iv:
        return Fraction(0)
    starts, lens, pref = make_prefix(Q_iv)
    LQ = pref[-1]
    r, p = sweep_acc(E_iv, R % T_DEN, starts, lens, pref)
    num = R * lenE_ - r
    require(num % T_DEN == 0, "floor-count not integral")
    return Fraction(LQ * (num // T_DEN) + p, R * T_DEN)


# --------------------------------------------------------------------------
# Toy brute-force self-tests (route (a))
# --------------------------------------------------------------------------
def toy_test():
    """Hostile little model: circle grid 405 = 3^4*5, base p = 3, sets with
    wraps and unequal pieces.  For EVERY toy packet clock R' in {3,9,27}
    and every m in {3,9,27}: direct rational computation of
    int (P_m f')(P_m e'), f' = 1_Q' P_{R'} 1_E', is compared against the
    y-side N/Phi reduction used as the main engine below.  The toy fold
    N_{R'} is also compared cellwise against the production count_fold
    (which must handle the R'=27 wrap case L > grid)."""
    Gr = 405
    Ep = [(0, 7), (50, 61), (100, 130), (200, 203), (390, 405)]
    Qp = [(10, 40), (120, 190), (300, 370)]

    def ind(iv, pos):
        return any(a <= pos < b for a, b in iv)

    for Rp in (3, 9, 27):
        # f' = 1_Q' * P_Rp 1_E' per Gr-cell (integer counts, units 1/Rp)
        fv = []
        for c in range(Gr):
            if not ind(Qp, c):
                fv.append(0)
                continue
            cnt = 0
            for b in range(Rp):
                q = c + b * Gr          # preimage cell on the Rp*Gr grid
                if ind(Ep, q // Rp):
                    cnt += 1
            fv.append(cnt)
        ev = [1 if ind(Ep, c) else 0 for c in range(Gr)]
        # toy multiplicity fold N_{Rp}(c) == fv-count without the Q factor
        nv = []
        for c in range(Gr):
            nv.append(sum(1 for b in range(Rp)
                          if ind(Ep, (c + b * Gr) // Rp)))
        # production count_fold cross-check (integer breakpoints)
        st, vl, cm, tot = count_fold(Ep, Rp, Gr)
        for c in range(Gr):
            i = bisect_right(st, c) - 1
            require(vl[i] == nv[c],
                    f"count_fold != brute fold at R'={Rp}, cell {c}")
        require(tot == Rp * measure(Ep), "toy fold mass")

        for m in (3, 9, 27):
            def Pm(h):
                out = []
                for c in range(Gr):
                    sacc = 0
                    for b in range(m):
                        q = c + b * Gr  # preimage cell on the m*Gr grid
                        sacc += h[q // m]
                    out.append(sacc)
                return out

            Pf = Pm(fv)                 # units 1/(m*Rp)
            Pe = Pm(ev)                 # units 1/m
            I_brute = Fraction(sum(x * y for x, y in zip(Pf, Pe)),
                               Gr * m * m * Rp)
            # y-side: I = sum_c 1_Q(c) N(c) PhiCnt(c) / (Gr*Rp*m)
            require(Gr % m == 0, "toy grid must resolve m")
            stp = Gr // m
            acc = 0
            for c in range(Gr):
                if not ind(Qp, c):
                    continue
                if not nv[c]:
                    continue
                phic = sum(1 for j in range(m)
                           if ind(Ep, (c + j * stp) % Gr))
                acc += nv[c] * phic
            I_y = Fraction(acc, Gr * Rp * m)
            require(I_brute == I_y,
                    f"toy y-side mismatch at R'={Rp}, m={m}")
    log("    toy brute-force == y-side N/Phi reduction at all")
    log("    R' in {3,9,27} x m in {3,9,27} (wrap cases included): PASS")
    log("    toy count_fold cellwise == brute multiplicity fold: PASS")


# --------------------------------------------------------------------------
# Phi_s prefix engine (verbatim from Stage-1)
# --------------------------------------------------------------------------
def phi_engine(s, E_iv, lenE_):
    """Exact machinery for Phi_s(x) = (1/m) sum_{j<m} 1_E(x + j/m),
    m = 13^(s+1), profile on the G = g*T_DEN grid, g = 13^max(0, s+1-6).
    query(P) = m*G * int_0^{P/G} Phi_s(y) dy, valid for any P >= 0."""
    m = 13 ** (s + 1)
    g = 13 ** max(0, (s + 1) - 6)
    G = g * T_DEN
    Tm = G // m
    require(Tm * m == G, "grid does not resolve the translate step")
    Eg = [(a * g, b * g) for a, b in E_iv]
    starts, vals, cums, total = count_fold(Eg, 1, Tm)
    lenEg = lenE_ * g
    require(total == lenEg, "Phi fold mass mismatch")

    def query(P):
        q, r0 = divmod(P, Tm)
        i = bisect_right(starts, r0) - 1
        return q * lenEg + cums[i] + vals[i] * (r0 - starts[i])

    return m, g, G, query


# --------------------------------------------------------------------------
# Main y-side engine:  I_s = int 1_Q (N_R/R) Phi_s dx
# --------------------------------------------------------------------------
def yside_Is(Q_iv, nfold, R, phi):
    m, g, G, query = phi
    n_st, n_v, _, _ = nfold
    npieces = len(n_st)
    acc = 0
    for qa, qb in Q_iv:
        i = bisect_right(n_st, qa) - 1
        while True:
            pa = n_st[i]
            pb = n_st[i + 1] if i + 1 < npieces else T_DEN
            a = qa if qa > pa else pa
            b = qb if qb < pb else pb
            if a < b and n_v[i]:
                acc += n_v[i] * (query(b * g) - query(a * g))
            if pb >= qb:
                break
            i += 1
    return Fraction(acc, R * m * G)


# --------------------------------------------------------------------------
# Packet-walk route (the Stage-1 MAIN engine, generalized to clock R)
# --------------------------------------------------------------------------
def packet_F(E_iv, Q_iv, R):
    """F = E n (R.)^-1 Q on the R*T_DEN grid, single-wrap walk.
    Returns None when the precondition R*(B-A) < T_DEN fails."""
    for A, B in E_iv:
        if R * (B - A) >= T_DEN:
            return None
    Qs = [a for a, _ in Q_iv]
    nQ = len(Q_iv)
    if nQ == 0:
        return []
    F = []
    for A, B in E_iv:
        LA = R * A
        sA = LA % T_DEN
        span = R * (B - A)
        sE = sA + span
        LAoff = LA - sA
        idx = bisect_right(Qs, sA) - 1
        off = 0
        if idx < 0:
            idx = nQ - 1
            off = -T_DEN
        while True:
            qa0, qb0 = Q_iv[idx]
            qa = qa0 + off
            if qa >= sE:
                break
            qb = qb0 + off
            if qb > sA:
                lo = sA if qa < sA else qa
                hi = sE if qb > sE else qb
                if hi > lo:
                    F.append((LAoff + lo, LAoff + hi))
            idx += 1
            if idx == nQ:
                idx = 0
                off += T_DEN
    return F


def froute_Is(F, R, phi):
    m, g, G, query = phi
    acc = 0
    for a, b in F:
        acc += query(b * g) - query(a * g)
    return Fraction(acc, R * m * G)


# --------------------------------------------------------------------------
# Deep-root sidecar branch (THM-2471 (44)-(47); c_j = 13, C = 2*13^5)
# --------------------------------------------------------------------------
def sidecar_branch(r):
    d = 13 ** r
    C = 2 * 13 ** 5
    delta = gcd(C, d)
    d0 = d // delta
    if d0 > 1:
        require(r >= 6, "d0>1 must mean r>=6 in the strict profile")
        return f"r>=6 SHEET: a mod 13^{r - 5} indispensable"
    h = C // d
    if h % 13 == 0:
        require(r <= 4, "13|h must mean r<=4 in the strict profile")
        return f"BASE-ONLY (h = C/d = {h} = 2*13^{nu13(h)}, 13 | h)"
    require(r == 5 and h == 2, "13 coprime h must mean r=5, h=2")
    return "r=5 WINDOW: unique affine invariant theta = t - 2u"


# ==========================================================================
def main():
    log("r-index census across the lawful (j, sigma, K) space of the"
        " canonical typed row")
    log("script=04-computation/lrc14_r_index_census_opus_20260728.py")
    log("machine=opus  script-date=2026-07-28  run-date=2026-07-27")
    log("")
    log("[1] row, definitions, and logged conventions")
    log(f"    w=(H,q1..q5,c1,c2,c3)={W}   profile nu13(c1,c2,c3)=(1,3,5)")
    log(f"    owner j=1: c = c_1 = {C1}, lambda = {LAMBDA};"
        f" canonical packet clock K = lambda+1 = {K_CANON}")
    log("    e = 1_(E_1);  f = 1_(Q_1,sigma) P^K 1_(E_1)   (THM-2471 (23),")
    log("    full word, source-safe factor RETAINED; K generalized per")
    log("    THM-2306 Sec.6 item 4 / THM-2471 Sec.3 positive-return clocks)")
    log("    A = P_13 f, B = P_13 e;  I_s = int (P^s A)(P^s B)"
        " = int (P_(13^(s+1)) f)(P_(13^(s+1)) e)")
    log("    r = min{s>=1 : I_s > 0}   (THM-2471 (24))")
    log("    P_a h(y) = (1/a) sum_(b<a) h((y+b)/a); half-open [A,B) intervals")
    log(f"    on the T_DEN = 182*lcm(w) = {T_DEN} grid")
    log("    MAIN ENGINE (y-side reduction, = Stage-1 route (f), self-tested):")
    log("      I_s = int_T 1_Q(x) (N_R(x)/R) Phi_s(x) dx,   R = 13^K,")
    log("      N_R = multiplicity fold of R*E,")
    log("      Phi_s(x) = (1/m) sum_(j<m) 1_E(x + j/m),  m = 13^(s+1)")
    log("")

    log("[2] licensing census from the files (exact checks)")
    log("    OWNERS: THM-2309 (25) types (w_j,w_a,w_b) = (13, 13^3, 2*13^5):")
    log("      owner j=1 is the depth-one blocker; a = c_2, b = c_3.")
    log("      THM-2349 (2)/(23a): a SECOND shallow owner is licensed only on")
    log(f"      repeated-first rows (nu13(c_2)=1); here nu13(c_2) ="
        f" {nu13(C2)} (strict")
    log("      profile 1<3<5)  =>  j = 2 is NOT LAWFUL.  Census owner: j=1"
        " ONLY.")
    log("    WORDS (THM-2305 (3), fixed order): sigma in {{a},{b},{a,b}}:")
    log("      Q_1,{a}  = A_0 n D_c2 \\ (D_c1 u D_c3)   (= E_2 literally)")
    log("      Q_1,{b}  = A_0 n D_c3 \\ (D_c1 u D_c2)   (= E_3 literally)")
    log("      Q_1,{a,b}= A_0 n D_c2 n D_c3 \\ D_c1")
    log("      EMPTY WORD: not a THM-2305 canonical stratum (partition (4)")
    log("      has exactly three cells under the cover reading; THM-2306")
    log("      Sec.1 and THM-2471 (23) restrict Q to the three canonical")
    log("      words)  =>  NOT LAWFUL; computed below as DIAGNOSTIC only.")
    log("    CLOCKS: K = 2 canonical (THM-2471 (23)); K = 3, 4 lawful iff")
    log("      rho^(K) = int_Q P^K 1_E > 0 (THM-2306 (41a) + Sec.6 item 4;")
    log("      THM-2471 Sec.3), decided exactly in [5].")
    log("")

    log("[3] machinery self-test (hostile toy model, brute force)")
    toy_test()
    log("")

    log("[4] exact interval sets, anchors, and the typed-row word partition")
    E = build_set(PAT_E, ZELL)
    lenE = check_intervals(E, T_DEN)
    require(Fraction(lenE, T_DEN) == Fraction(1882176, 28589561),
            "mu(E_1) anchor")
    require(len(E) == 57072, "E_1 interval count anchor")
    log(f"    E_1: intervals={len(E)}  measure={Fraction(lenE, T_DEN)}: PASS")
    maxlenE = max(b - a for a, b in E)
    log(f"    maxlen(E_1) = {maxlenE}/T_DEN"
        f"  (single-wrap F-walk lawful iff 13^K*maxlen < T_DEN)")

    QA = build_set(PAT_QA, ZELL)
    lenQA = check_intervals(QA, T_DEN)
    require(Fraction(lenQA, T_DEN) == Fraction(143103830843, 5727632650740),
            "mu(Q_1,{a}) anchor")
    require(len(QA) == 22478, "Q_1,{a} interval count anchor")
    log(f"    Q_1,{{a}}:   intervals={len(QA)}"
        f"  measure={Fraction(lenQA, T_DEN)}: PASS (Stage-1 anchor)")

    QB = build_set(PAT_QB, ZELL)
    lenQB = check_intervals(QB, T_DEN)
    log(f"    Q_1,{{b}}:   intervals={len(QB)}"
        f"  measure={Fraction(lenQB, T_DEN)}")
    QAB = build_set(PAT_QAB, ZELL)
    lenQAB = check_intervals(QAB, T_DEN)
    log(f"    Q_1,{{a,b}}: intervals={len(QAB)}"
        f"  measure={Fraction(lenQAB, T_DEN)}")
    QE = build_set(PAT_QE, ZELL)
    lenQE = check_intervals(QE, T_DEN)
    log(f"    empty-word residual A_0\\(Dc1uDc2uDc3): intervals={len(QE)}"
        f"  measure={Fraction(lenQE, T_DEN)}")
    R1 = build_set(PAT_R1, ZELL)
    lenR1 = check_intervals(R1, T_DEN)
    log(f"    R_1 = A_0\\Dc1: intervals={len(R1)}"
        f"  measure={Fraction(lenR1, T_DEN)}")

    # typed-row partition: R_1 = Qa u Qb u Qab u empty, exactly
    require(lenQA + lenQB + lenQAB + lenQE == lenR1,
            "word partition measure additivity")
    require(merge_touch(QA + QB + QAB + QE) == merge_touch(R1),
            "word partition interval-exact union")
    for X, Y, nm in ((QA, QB, "a/b"), (QA, QAB, "a/ab"), (QA, QE, "a/e"),
                     (QB, QAB, "b/ab"), (QB, QE, "b/e"), (QAB, QE, "ab/e")):
        require(intersect_lists(X, Y) == [], f"strata {nm} not disjoint")
    log("    partition R_1 = Q_{a} u Q_{b} u Q_{a,b} u empty-residual:")
    log("    interval-exact, pairwise disjoint, measures additive: PASS")
    log("    NOTE (typed row, MISTAKE-281): the empty-word residual has")
    log(f"    POSITIVE measure {Fraction(lenQE, T_DEN)} ~"
        f" {float(Fraction(lenQE, T_DEN)):.6e};")
    log("    the scalar-cover identity THM-2305 (4) fails on this typed row")
    log("    by exactly that mass.  This does not make the empty word a")
    log("    lawful stratum; it quantifies the typed/cover gap.")
    log("    DIAGNOSTIC (owner j=2): E_2 = Q_1,{a} literally (THM-2305 (3)),")
    log(f"    mu(E_2) = {Fraction(lenQA, T_DEN)} > 0, yet j=2 stays"
        " UNLICENSED:")
    log("    the THM-2349 second-owner license is profile-typed"
        " (repeated-first")
    log("    only), and measure positivity is not a license.")
    log("")

    log("[5] structural I_0 = 0 and the exact returns rho^(K) (licensing)")
    D13 = in_comb(OWNER, 0)             # D_c1 = {||13x|| < 1/14}
    require(measure(intersect_lists(E, D13)) == lenE, "E_1 not inside D_c1")
    for Q_iv, nm in ((QA, "Q_{a}"), (QB, "Q_{b}"), (QAB, "Q_{a,b}"),
                     (QE, "empty"), (R1, "R_1")):
        require(intersect_lists(Q_iv, D13) == [], f"{nm} meets D_c1")
    log("    E_1 subset D_c1 and every stratum n D_c1 = empty"
        " (interval-exact):")
    log("    A = P_13 f and B = P_13 e are disjoint at EVERY (sigma,K)")
    log("    => I_0 = 0 structurally (THM-2306 normalization premise).")
    log("")
    SETS = [("{a}", QA, lenQA), ("{b}", QB, lenQB), ("{a,b}", QAB, lenQAB),
            ("empty*", QE, lenQE), ("R_1*", R1, lenR1)]
    KLIST = (2, 3, 4)
    rho = {}
    for K in KLIST:
        R = 13 ** K
        parts = 0
        for nm, Q_iv, _ln in SETS:
            rv = IR_exact(E, lenE, Q_iv, R)
            rho[(nm, K)] = rv
            if nm != "R_1*":
                parts += rv
            tag = "" if rv == 0 else f"  ~{float(rv):.6e}"
            log(f"    K={K}: rho^({K})({nm}) = int_Q P^{K} 1_E = {rv}{tag}")
        require(parts == rho[("R_1*", K)], f"rho additivity fails at K={K}")
        log(f"         additivity sum(strata+empty) == rho(R_1): PASS")
    require(rho[("{a}", 2)] == Fraction(21376087, 17907461390),
            "rho^(2)(Q_{a}) != stored Stage-1 stratum measure")
    log("    anchor rho^(2)({a}) == Stage-1 mu(F) = 21376087/17907461390:"
        " PASS")
    log("    LICENSE TABLE (lawful census configs = canonical word AND")
    log("    rho^(K) > 0; empty word and R_1 are diagnostics, marked *):")
    lawful = []
    for nm, Q_iv, _ln in SETS:
        for K in KLIST:
            pos = rho[(nm, K)] > 0
            if nm in ("{a}", "{b}", "{a,b}") and pos:
                lawful.append((nm, K))
                verdict = "LAWFUL"
            elif nm in ("{a}", "{b}", "{a,b}"):
                verdict = "NOT LAWFUL (zero return)"
            else:
                verdict = "DIAGNOSTIC ONLY"
            log(f"      (j=1, sigma={nm}, K={K}): rho>0 ="
                f" {pos}  =>  {verdict}")
    log("")

    log("[6] census: exact I_s per (sigma, K) until first positive"
        " (cap s<=7)")
    S_CAP = 7
    phi_cache = {}

    def phi(s):
        if s not in phi_cache:
            phi_cache[s] = phi_engine(s, E, lenE)
        return phi_cache[s]

    Ivals = {}          # (nm, K) -> {s: Fraction}
    r_idx = {}          # (nm, K) -> int or None
    nfolds = {}
    for K in KLIST:
        R = 13 ** K
        nfolds[K] = count_fold(E, R, T_DEN)
        log(f"    {elapsed()} K={K}: N_R fold ready"
            f" ({len(nfolds[K][0])} pieces, mass check PASS)")
    for K in KLIST:
        R = 13 ** K
        nf = nfolds[K]
        for nm, _q, _l in SETS:
            Ivals[(nm, K)] = {}
            r_idx[(nm, K)] = None
        for s in range(0, S_CAP + 1):
            m, g, G, query = phi(s)
            ctrl = query(R * T_DEN * g) - query(0)
            require(Fraction(ctrl, m * G) == R * Fraction(lenE, T_DEN),
                    f"full-circle Phi control failed at K={K}, s={s}")
            tot_parts = 0
            for nm, Q_iv, _ln in SETS:
                Is = yside_Is(Q_iv, nf, R, phi(s))
                require(Is >= 0, "I_s negative (impossible)")
                Ivals[(nm, K)][s] = Is
                if nm != "R_1*":
                    tot_parts += Is
                if r_idx[(nm, K)] is None and s >= 1 and Is > 0:
                    r_idx[(nm, K)] = s
            require(tot_parts == Ivals[("R_1*", K)][s],
                    f"I_s additivity fails at K={K}, s={s}")
            require(Ivals[("R_1*", K)][s - 1] == 0 if s == 1 else True,
                    "I_0 != 0")
            shown = "  ".join(
                f"{nm}:{'0' if Ivals[(nm, K)][s] == 0 else 'POS'}"
                for nm, _q, _l in SETS)
            log(f"    {elapsed()} K={K} s={s} m=13^{s + 1}:  {shown}"
                f"   [Phi control + additivity PASS]")
            if all(r_idx[(nm, K)] is not None
                   for nm, _q, _l in SETS):
                break
        for nm, _q, _l in SETS:
            require(Ivals[(nm, K)][0] == 0,
                    f"I_0 != 0 at ({nm}, K={K}) contradicts [5]")
    log("")

    log("[7] independent verification matrix")
    log("    (d) support-fold route (union folds only, no prefix engine):")
    for K in KLIST:
        R = 13 ** K
        UR = union_fold(E, R, T_DEN)
        for nm, Q_iv, _ln in SETS:
            Fsup = intersect_lists(Q_iv, UR)
            smax = max(Ivals[(nm, K)].keys())
            pat_eng = []
            pat_sup = []
            for s in range(0, smax + 1):
                m = 13 ** (s + 1)
                Se = union_fold(E, m, T_DEN)
                Sf = union_fold(Fsup, m, T_DEN)
                ol = measure(intersect_lists(Se, Sf))
                pat_sup.append(ol > 0)
                pat_eng.append(Ivals[(nm, K)][s] > 0)
                require((ol > 0) == (Ivals[(nm, K)][s] > 0),
                        f"support route contradicts engine at"
                        f" ({nm},K={K},s={s})")
            log(f"        K={K} {nm}: zero/positive pattern"
                f" s=0..{smax} agrees with engine: PASS")
    log("        (nonnegative step functions: null support overlap <=> zero")
    log("         integral; positive overlap => strictly positive integral)")
    log("    (e) correlation route s=0,1 (replica-dichotomy prefix"
        " identity),")
    log("        four strata x three clocks (R_1 covered by additivity):")
    for s in (0, 1):
        m = 13 ** (s + 1)
        require(T_DEN % m == 0, "correlation route needs m | T_DEN")
        sh_unit = T_DEN // m
        tots = {}
        for j in range(m):
            Esh = rotate(E, j * sh_unit, T_DEN)
            for nm, Q_iv, _ln in SETS:
                if nm == "R_1*":
                    continue
                Qp = intersect_lists(Q_iv, Esh)
                for K in KLIST:
                    tots[(nm, K)] = tots.get((nm, K), Fraction(0)) + \
                        IR_exact(E, lenE, Qp, 13 ** K)
        for nm, _q, _l in SETS:
            if nm == "R_1*":
                continue
            for K in KLIST:
                tv = tots[(nm, K)] / m
                require(tv == Ivals[(nm, K)][s],
                        f"correlation route mismatch at ({nm},K={K},s={s})")
        log(f"        {elapsed()} s={s}: all 12 (sigma,K) values"
            f" == engine: PASS")
    log("    (f) packet-walk route (Stage-1 main engine, generalized R),")
    log("        at s in {r-1, r} where the single-wrap precondition"
        " holds:")
    F_SIZE_CAP = 5_000_000
    for K in KLIST:
        R = 13 ** K
        if R * maxlenE >= T_DEN:
            log(f"        K={K}: precondition 13^K*maxlen(E) < T_DEN FAILS"
                f" (R*maxlen = {R * maxlenE} >= T_DEN); route skipped --")
            log("        the y-side engine is exactly the audited Stage-1")
            log("        route (f) and needs no wrap bound.")
            continue
        for nm, Q_iv, _ln in SETS:
            # expected F interval count ~ |Q| * R * mu(E) + |E|; the walk
            # materializes F, so cap it (the route is a redundant 4th check;
            # engine/support/correlation/additivity already decide all).
            est = len(Q_iv) * R * lenE // T_DEN + len(E)
            if est > F_SIZE_CAP:
                log(f"        K={K} {nm}: skipped (expected F size"
                    f" ~{est} intervals exceeds cap {F_SIZE_CAP};"
                    " three other exact routes cover this config)")
                continue
            F = packet_F(E, Q_iv, R)
            require(F is not None, "precondition checked above")
            lenF = check_intervals(F, R * T_DEN)
            require(Fraction(lenF, R * T_DEN) == rho[(nm, K)],
                    f"mu(F) != rho at ({nm},K={K})")
            rr = r_idx[(nm, K)]
            if rr is None:
                log(f"        K={K} {nm}: mu(F) == rho^({K}): PASS"
                    " (r undecided at cap, no I_r check)")
                continue
            for s in (rr - 1, rr):
                Iw = froute_Is(F, R, phi(s))
                require(Iw == Ivals[(nm, K)][s],
                        f"packet-walk mismatch at ({nm},K={K},s={s})")
            log(f"        K={K} {nm}: mu(F) == rho^({K}), and packet-walk"
                f" I_{rr - 1}, I_{rr} == engine: PASS")
    log("    verification matrix complete: engine values confirmed by the")
    log("    support route (every s), correlation route (s=0,1), and the")
    log("    packet-walk route (where lawful), plus partition additivity")
    log("    and full-circle controls at every step.")
    log("")

    log("[8] THE CENSUS TABLE  (j=1, sigma, K) -> r, exact first-positive"
        " I_r")
    log("    (columns: config | lawful | r | I_r exact | 169*I_r | deep-root"
        " branch)")
    stage1_anchor = Fraction(9926558757352, 109707098520974955)
    require(Ivals[("{a}", 2)][3] == stage1_anchor,
            "Stage-1 I_3 anchor mismatch at ({a}, K=2)")
    require(r_idx[("{a}", 2)] == 3, "Stage-1 r=3 anchor mismatch")
    log("    anchor: (sigma={a}, K=2) reproduces the audited Stage-1"
        " r = 3 and")
    log(f"    I_3 = {stage1_anchor}: PASS")
    log("")
    window_open = []
    for nm, _q, _l in SETS:
        for K in KLIST:
            rr = r_idx[(nm, K)]
            isl = (nm, K) in lawful
            tagl = "LAWFUL " if isl else ("DIAG   " if nm in
                                          ("empty*", "R_1*") else "UNLAWFL")
            if rr is None:
                log(f"    (j=1, sigma={nm:6s}, K={K}) {tagl}:"
                    f"  r >= {S_CAP + 1} (cap; undecided)")
                continue
            Ir = Ivals[(nm, K)][rr]
            br = sidecar_branch(rr)
            if isl and rr == 5:
                window_open.append((nm, K))
            log(f"    (j=1, sigma={nm:6s}, K={K}) {tagl}:  r = {rr}")
            log(f"        I_r = {Ir}  ~{float(Ir):.6e}")
            log(f"        169*I_r = {169 * Ir}")
            log(f"        zeros: I_1..I_{rr - 1} = 0 exactly;"
                f"  branch: {br}")
    log("")

    log("[9] DECISION VALUE and pattern report")
    if window_open:
        log(f"    r = 5 WINDOW OPENS at lawful configs: {window_open}")
    else:
        log("    The r = 5 window (unique-intertwiner case theta = t - 2u)")
        log("    does NOT open at any lawful configuration: every lawful")
        log("    (j=1, sigma, K) census value satisfies r <= 4 or r >= 6 as")
        log("    printed above -- see the exact table.")
    rvals = {}
    for nm, _q, _l in SETS:
        if nm in ("empty*", "R_1*"):
            continue
        for K in KLIST:
            if (nm, K) in lawful and r_idx[(nm, K)] is not None:
                rvals[(nm, K)] = r_idx[(nm, K)]
    log("    lawful r-table: " + "  ".join(
        f"({nm},K={K})->r={v}" for (nm, K), v in sorted(rvals.items())))
    all3 = all(v == 3 for v in rvals.values())
    if all3:
        log("    PATTERN: r = 3 = nu13(c_2) at EVERY lawful configuration;")
        log("    r tracks the middle-blocker valuation and moves with")
        log("    NEITHER the word sigma NOR the packet clock K in the")
        log("    lawful range censused here.")
    else:
        byK = {}
        bys = {}
        for (nm, K), v in rvals.items():
            byK.setdefault(K, set()).add(v)
            bys.setdefault(nm, set()).add(v)
        log(f"    PATTERN: r varies.  by clock K: "
            + "  ".join(f"K={K}:{sorted(vs)}" for K, vs in sorted(byK.items())))
        log(f"             by word sigma: "
            + "  ".join(f"{nm}:{sorted(vs)}" for nm, vs in sorted(bys.items())))
        log("    (compare nu13(c_2) = 3, nu13(c_3) = 5; exact table in [8]")
        log("    is the record -- this paragraph is description, not claim)")
    for K in KLIST:
        if ("empty*", K) in r_idx and r_idx[("empty*", K)] is not None:
            log(f"    diagnostic (NOT lawful): empty word at K={K} has"
                f" r = {r_idx[('empty*', K)]}")
    log("")

    log("[10] scope notes (MISTAKE-281/283 discipline)")
    log("    * TYPED row (THM-2309 (25)); NOT an asserted scalar cover; the")
    log("      empty-word residual mass in [4] is the exact typed/cover gap")
    log("      on this row and licenses nothing.")
    log("    * Census scope: owner j=1 only (THM-2349 profile-typed"
        " license);")
    log("      words = the three THM-2305 canonical strata; clocks K=2,3,4")
    log("      via THM-2306 (41a)/Sec.6-item-4 + THM-2471 Sec.3 with the")
    log("      return positivity computed exactly.  Nothing beyond these")
    log("      configurations is claimed lawful, and nothing about other")
    log("      rows, physical currents, row exclusion, LRC(14) or JC(2)")
    log("      follows from this census.")
    log("    * Every I_s is one integral of a triple product at one common")
    log("      base point (paired on one base BEFORE marginalization).")
    log("    * No intertwiner was constructed and none is claimed; the")
    log("      deep-root branches are THM-2471 (44)-(47) arithmetic applied")
    log("      to the exact computed r values, nothing more.")
    log("    * All decisions are integer/Fraction-exact on the T_DEN and")
    log("      13^max(0,s-5)*T_DEN grids; floats appear only in ~displays.")
    log("")
    log(f"{elapsed()} all checks passed")


if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    main()
