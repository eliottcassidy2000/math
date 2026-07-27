#!/usr/bin/env python3
"""Two-stage ancestry/Boolean realization test on the canonical typed row:
STAGE 1 -- the THM-2471 (23)-(24) first-collision index r, decided exactly.

TASK (THM-2512 Sec. 5 bridge test + THM-2471 Secs. 6-7).  On the canonical
typed row THM-2309 (25),

    w = (H,q1..q5,c1,c2,c3) = (1,14,27,40,53,66,13,2197,742586),
    owner j = 1 (c_1 = 13, lambda = nu13(c_1) = 1, K = lambda+1 = 2),
    word sigma = {a}  (Q_{1,{a}} = A_0 n D_{c2} \\ (D_{c1} u D_{c3})),

take the exact THM-2471 (23) data

    e = 1_{E_1},          f = 1_{Q_{1,{a}}} P^2 1_{E_1},
    c = c_1 = 13,         A = P_c f,   B = P_c e,

and compute the collision index of THM-2471 (24):

    I_s = int_T (P^s A)(P^s B),        r = min{s >= 1 : I_s > 0},

with P = P_13 the mult-by-13 transfer operator, so P^s A = P_{13^(s+1)} f
and P^s B = P_{13^(s+1)} e.  Branch per the strict-profile deep-root sidecar
THM-2471 (44)-(47) (c_j = 13, C = 2*13^5, d = 13^r -- exactly this row):

    r <= 4 :  deep root is BASE-ONLY (degenerate sub-case; STOP),
    r  = 5 :  unique affine invariant theta = t - 2u (proceed to Stage 2),
    r >= 6 :  sheet residue a mod 13^(r-5) indispensable (DECIDED NEGATIVE
              structurally; STOP).

REDUCTION (all steps exact; logged here because they carry the decision).
With m = 13^(s+1):

  (i)  adjoint of the transfer operator:  int (P_m f)(P_m e)
       = int f(x) (P_m e)({m x}) dx;
  (ii) preimage identity: the m preimages of {m x} under y -> {m y} are
       exactly {x + j/m : j = 0..m-1}, hence
       (P_m e)({m x}) = Phi_s(x),   Phi_s(x) = (1/m) sum_j 1_E(x + j/m);
  (iii) f = 1_Q * P_169 1_E and one more adjoint step give

       I_s = int_T 1_E(u) * 1_Q({169 u}) * Phi_s({169 u}) du.

So each I_s is ONE interval-overlap mass: the packet set
F = E n (169.)^{-1}Q (an explicit interval list on the NN = 169*T_DEN grid,
the same F_{0,0} as the audited owner-loop-drift script) integrated against
the 1/m-periodic fold of E, evaluated through an exact prefix/cumulative
identity (the same prefix-count mechanism as the audited replica-dichotomy
script, here applied to the multiplicity fold of E instead of a Boolean Q).

CONVENTIONS (each an explicit logged choice; files override summaries):
  * grid: all base sets are unions of half-open intervals [A,B) with integer
    endpoints on the T_DEN = 182*lcm(w) grid; the packet F lives on
    NN = 169*T_DEN; no floats enter any decision;
  * P_a h(y) = (1/a) sum_{b=0}^{a-1} h((y+b)/a)  (THM-2471 (2));
  * E_1 = build_set(PAT_E): guard-safe (||Hx|| > 1/7), five unit safes,
    c_1 DANGER d(13x) (owner 'in'), c_2-safe, c_3-safe -- identical to the
    audited owner-loop-drift script (anchor-checked);
  * Q_{1,{a}} = build_set(PAT_QA): guard-safe, unit safes, c_1-safe,
    D_{c2} 'in', c_3-safe -- the FULL word of THM-2471 (23) (its source-safe
    factor RETAINED; the source-deleted T_a of THM-2449 (15) is NOT used
    here, because (23) prescribes f = 1_Q P^K e with Q = Q_{j,sigma});
  * danger d(y) = 1_{||y|| < 1/14}, windows via the 182/91-denominator combs
    of the audited scripts;
  * K = lambda + 1 = 2, so the ancestry clock inside f is R = 13^K = 169
    (THM-2471 (28)); this is the packet clock, NOT a response clock.

VERIFICATION MATRIX (independent routes, all exact):
  (a) toy brute-force self-test of the reduction (grid 405, m = 3,9,27):
      direct rational computation of int (P_m f)(P_m e) on the fine grid
      == the Phi/prefix reduction, on a hostile little model;
  (b) anchors: mu(E_1), mu(Q_{1,{a}}), mu(F) and interval counts against the
      audited twist-variance/owner-loop artifacts;
  (c) structural owner-normalized disjointness: E_1 subset D_{c1} and
      Q n D_{c1} = empty interval-wise, forcing supp(P_13 f) n supp(P_13 e)
      = empty, i.e. I_0 = 0 (THM-2306 normalization premise);
  (d) support-fold route (NO shared machinery with the prefix engine):
      supp P_m h = {m x mod 1 : x in supp h}, computed as a union fold;
      I_s = 0 iff the two supports meet in measure zero, and (nonneg step
      functions) a positive-measure overlap forces I_s > 0.  This decides
      zero/positive for EVERY s independently of the exact-value engine;
  (e) correlation route for s = 0,1 (machinery of the audited
      replica-dichotomy script): I_s = (1/m) sum_j I_169(E, Q n (E - j/m))
      via the exact prefix-count identity, 13 + 169 lag sweeps;
  (f) y-side route for s = r-1, r: I_s = int 1_Q (N_169/169) Phi_s dx with
      N_169 the multiplicity fold of 169*E -- a different decomposition of
      the same integral (piece walk on the T_DEN grid instead of the
      169-substituted F-walk);
  (g) full-circle control at every s: replacing F by the full circle must
      return mu(E) exactly.

SCOPE (MISTAKE-281 discipline, hard):
  * the row is TYPED (THM-2309 (25)); it is NOT an asserted scalar cover;
  * no physical current is claimed, no row exclusion, no LRC(14) progress
    beyond what THM-2471/THM-2512 literally license on this one row;
  * full support is never noncancellation; every claim below is an exact
    computed identity about THIS row's packet;
  * all pairings below live on one common base before any marginalization
    (each I_s is a single integral of a triple product at one point u).

Script: 04-computation/lrc14_ancestry_realization_stage_opus_20260727.py
Output: 05-knowledge/results/lrc14_ancestry_realization_stage_opus_20260727.out
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
NN = RPKT * T_DEN
require(T_DEN == 297836897838480, "T_DEN anchor")
require(NN == 50334435734703120, "NN anchor")
require(nu13(T_DEN) == 6, "13-adic valuation of T_DEN must be 6")


# --------------------------------------------------------------------------
# Interval machinery (verbatim from the audited owner-loop-drift script)
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
# Generic exact helpers
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
    Returns (starts, vals, cums, total): value vals[i] on
    [starts[i], starts[i+1]); cums[i] = integral up to starts[i];
    total = integral over [0,grid) = mult * measure(iv)."""
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
# Exact prefix identity of the replica-dichotomy script (for route (e))
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
    starts, lens, pref = make_prefix(Q_iv)
    LQ = pref[-1]
    r, p = sweep_acc(E_iv, R % T_DEN, starts, lens, pref)
    num = R * lenE_ - r
    require(num % T_DEN == 0, "floor-count not integral")
    return Fraction(LQ * (num // T_DEN) + p, R * T_DEN)


# --------------------------------------------------------------------------
# Toy brute-force self-test of the reduction (route (a))
# --------------------------------------------------------------------------
def toy_test():
    """Hostile little model: circle grid 405 = 3^4*5, base p = 3, packet
    clock R' = 9, sets E',Q' with wraps and unequal pieces.  Checks
    int (P_m f')(P_m e') computed by direct rational averaging on the fine
    m-grid == the Phi/preimage reduction, for m = 3, 9, 27."""
    Gr = 405
    Rp = 9
    Ep = [(0, 7), (50, 61), (100, 130), (200, 203), (390, 405)]
    Qp = [(10, 40), (120, 190), (300, 370)]

    def ind(iv, pos):
        return any(a <= pos < b for a, b in iv)

    for m in (3, 9, 27):
        # f' = 1_Q' * P_Rp 1_E' per Gr-cell (values cnt/Rp)
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

        def Pm(h):
            out = []
            for c in range(Gr):
                sacc = 0
                for b in range(m):
                    q = c + b * Gr      # preimage cell on the m*Gr grid
                    sacc += h[q // m]
                out.append(sacc)
            return out

        Pf = Pm(fv)                     # units 1/(m*Rp)
        Pe = Pm(ev)                     # units 1/m
        I_brute = Fraction(sum(x * y for x, y in zip(Pf, Pe)),
                           Gr * m * m * Rp)
        # reduction: I = int 1_E(u) 1_Q({Rp u}) Phi({Rp u}) du on NN2 = Rp*Gr
        NN2 = Rp * Gr
        require(NN2 % m == 0, "toy grid must resolve m")
        stp = NN2 // m
        acc = 0
        for u in range(NN2):
            if not ind(Ep, u // Rp):
                continue
            x = (Rp * u) % NN2
            for q in range(x, x + Rp):
                qq = q % NN2
                if not ind(Qp, qq // Rp):
                    continue
                acc += sum(1 for j in range(m)
                           if ind(Ep, ((qq + j * stp) % NN2) // Rp))
        I_red = Fraction(acc, NN2 * Rp * m)
        require(I_brute == I_red, f"toy self-test mismatch at m={m}")
    log("    toy brute-force == Phi/preimage reduction at m=3,9,27: PASS")


# --------------------------------------------------------------------------
# Phi_s prefix engine
# --------------------------------------------------------------------------
def phi_engine(s, E_iv, lenE_):
    """Exact machinery for Phi_s(x) = (1/m) sum_{j<m} 1_E(x + j/m),
    m = 13^(s+1).  Phi_s has period 1/m; its one-period profile is the
    multiplicity fold of E into [0, Tm) on the G = g*T_DEN grid,
    g = 13^max(0, s+1-6) (the translate step j*G/m must be integral).
    Returns (m, g, G, query) with

        query(P) = m*G * int_0^{P/G} Phi_s(y) dy
                 = (P // Tm) * lenE*g + PrefN(P mod Tm),

    valid for any P >= 0 (P/G may exceed 1: Phi_s is 1/m-periodic)."""
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


# ==========================================================================
def main():
    log("Two-stage ancestry/Boolean realization test on the canonical typed"
        " row -- STAGE 1")
    log("script=04-computation/lrc14_ancestry_realization_stage_opus_20260727.py")
    log("machine=opus  date=2026-07-27")
    log("")
    log("[1] row, packet, and logged conventions")
    log(f"    w=(H,q1..q5,c1,c2,c3)={W}   profile nu13(c1,c2,c3)=(1,3,5)")
    log(f"    owner j=1: c = c_1 = {C1}, lambda = {LAMBDA}, K = lambda+1 ="
        f" {KCLK}, packet clock R = 13^K = {RPKT}")
    log("    word sigma={a}: Q_1,{a} = A_0 n D_c2 \\ (D_c1 u D_c3), full word")
    log("    (source-safe factor RETAINED: THM-2471 (23) uses Q = Q_(j,sigma))")
    log("    e = 1_(E_1);  f = 1_(Q_1,{a}) P^2 1_(E_1)   (THM-2471 (23))")
    log("    A = P_13 f, B = P_13 e;  I_s = int (P^s A)(P^s B)"
        " = int (P_(13^(s+1)) f)(P_(13^(s+1)) e)")
    log("    r = min{s>=1 : I_s > 0}   (THM-2471 (24))")
    log("    P_a h(y) = (1/a) sum_(b<a) h((y+b)/a); half-open [A,B) intervals")
    log(f"    on the T_DEN = 182*lcm(w) = {T_DEN} grid; NN = 169*T_DEN")
    log("    REDUCTION (proved in the docstring, self-tested below):")
    log("      I_s = int_T 1_E(u) 1_Q({169u}) Phi_s({169u}) du,")
    log("      Phi_s(x) = (1/m) sum_(j<m) 1_E(x + j/m),  m = 13^(s+1)")
    log("")

    log("[2] machinery self-test (hostile toy model, brute force)")
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
    log(f"    Q_1,{{a}}: intervals={len(Q)}  measure={Fraction(lenQ, T_DEN)}:"
        " PASS")
    log("    anchors vs the audited twist-variance/owner-loop artifacts: PASS")
    log("")

    log("[4] packet set F = E n (169.)^(-1) Q on the NN grid (the x-side of f)")
    Qs = [a for a, _ in Q]
    nQ = len(Q)
    F = []
    for A, B in E:
        LA = 169 * A
        sA = LA % T_DEN
        span = 169 * (B - A)
        require(span < T_DEN, "E interval too long for single-wrap walk")
        sE = sA + span
        LAoff = LA - sA
        idx = bisect_right(Qs, sA) - 1
        off = 0
        if idx < 0:
            idx = nQ - 1
            off = -T_DEN
        while True:
            qa0, qb0 = Q[idx]
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
    lenF = check_intervals(F, NN)
    require(Fraction(lenF, NN) == Fraction(21376087, 17907461390),
            "mu(F) != stored stratum measure(E_1 n T^-2 Q_1,{a})")
    log(f"    F: intervals={len(F)}  measure={Fraction(lenF, NN)}")
    log("    == stored measure(E_1 n T^-2 Q_1,{a}) of the audited artifacts:"
        " PASS")
    log(f"    {elapsed()}")
    log("")

    log("[5] structural owner-normalized disjointness (I_0 = 0 premise)")
    D13 = in_comb(OWNER, 0)             # D_c1 = {||13x|| < 1/14}
    require(measure(intersect_lists(E, D13)) == lenE, "E_1 not inside D_c1")
    require(intersect_lists(Q, D13) == [], "Q meets D_c1")
    log("    E_1 subset D_c1 and Q_1,{a} n D_c1 = empty (interval-exact).")
    log("    Multiplication by 13 maps each D_c1 window onto (-1/14,1/14), so")
    log("    supp(P_13 e) subset {||y||<1/14} while supp(P_13 f) avoids it:")
    log("    A = P_13 f and B = P_13 e are disjoint => I_0 = 0 structurally")
    log("    (THM-2306 owner normalization; checked numerically below too).")
    log("")

    log("[6] STAGE 1: exact I_s until the first positive (THM-2471 (24))")
    S_CAP = 9
    Ivals = {}
    r_idx = None
    for s in range(0, S_CAP):
        m, g, G, query = phi_engine(s, E, lenE)
        # (g) full-circle control: F -> [0,NN) must return mu(E)
        ctrl = query(NN * g) - query(0)
        require(Fraction(ctrl, m * G) == 169 * Fraction(lenE, T_DEN),
                f"full-circle control failed at s={s}")
        acc = 0
        for a, b in F:
            acc += query(b * g) - query(a * g)
        Is = Fraction(acc, 169 * m * G)
        require(Is >= 0, "I_s negative (impossible)")
        Ivals[s] = Is
        tag = "" if Is == 0 else f"  ~{float(Is):.6e}"
        log(f"    {elapsed()} s={s}  m=13^{s + 1}  I_s = {Is}{tag}"
            f"   [full-circle control PASS]")
        if s >= 1 and Is > 0:
            r_idx = s
            break
    require(r_idx is not None, "no collision found below the cap")
    require(Ivals[0] == 0, "I_0 != 0 contradicts [5]")
    R_COLL = r_idx
    log(f"    FIRST POSITIVE: r = {R_COLL}   "
        f"(I_1..I_{R_COLL - 1} = 0, I_{R_COLL} > 0)")
    log(f"    I_r = {Ivals[R_COLL]}")
    log(f"    169*I_r = {169 * Ivals[R_COLL]}   (the (14) service mass M)")
    log("")

    log("[7] independent verification matrix")
    # (d) support-fold route: engine-independent zero/positivity for every s
    log("    (d) support-fold route (union folds only, no prefix engine):")
    U169 = union_fold(E, 169, T_DEN)
    Fsup = intersect_lists(Q, U169)     # supp f = Q n supp(P_169 1_E)
    log(f"        supp(P_169 1_E): {len(U169)} iv, measure "
        f"{Fraction(measure(U169), T_DEN)};  supp f: {len(Fsup)} iv")
    for s in range(0, R_COLL + 1):
        m = 13 ** (s + 1)
        require(T_DEN % m == 0, "support fold needs m | T_DEN")
        Se = union_fold(E, m, T_DEN)          # supp P_m e
        Sf = union_fold(Fsup, m, T_DEN)       # supp P_m f
        ol = measure(intersect_lists(Se, Sf))
        verdict = "ZERO" if ol == 0 else "POSITIVE"
        log(f"        s={s}: mu(supp P_m f n supp P_m e) = "
            f"{Fraction(ol, T_DEN)}  => I_s {verdict}")
        require((ol == 0) == (Ivals[s] == 0),
                f"support route contradicts engine at s={s}")
    log("        (nonnegative step functions: null support overlap <=> zero")
    log("         integral; positive overlap => strictly positive integral)")
    log("        support route == engine zero/positive pattern: PASS")
    # (e) correlation route for s = 0, 1
    log("    (e) correlation route (replica-dichotomy prefix identity):")
    for s in (0, 1):
        m = 13 ** (s + 1)
        sh_unit = T_DEN // m
        tot = Fraction(0)
        for j in range(m):
            Esh = rotate(E, j * sh_unit, T_DEN)
            Qp = intersect_lists(Q, Esh)
            tot += IR_exact(E, lenE, Qp, 169)
        tot /= m
        require(tot == Ivals[s], f"correlation route mismatch at s={s}")
        log(f"        {elapsed()} s={s}: (1/m) sum_j I_169(E, Q n (E - j/m))"
            f" = {tot} == engine I_{s}: PASS")
    # (f) y-side route for s = r-1, r
    log("    (f) y-side route  I_s = int 1_Q (N_169/169) Phi_s dx :")
    n_st, n_v, n_c, n_tot = count_fold(E, 169, T_DEN)
    require(n_tot == 169 * lenE, "N_169 fold mass")
    npieces = len(n_st)
    for s in (R_COLL - 1, R_COLL):
        m, g, G, query = phi_engine(s, E, lenE)
        acc = 0
        for qa, qb in Q:
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
        Iy = Fraction(acc, 169 * m * G)
        require(Iy == Ivals[s], f"y-side route mismatch at s={s}")
        log(f"        {elapsed()} s={s}: y-side I_{s} = {Iy} == engine: PASS")
    log("    verification matrix complete: 3 independent routes agree with")
    log("    the engine on every decided value.")
    log("")

    log("[8] minimality premise and deep-root sidecar (THM-2471 (25)-(26),"
        " (44)-(47))")
    d_coll = C1 * 13 ** (R_COLL - 1)
    log(f"    d = c*13^(r-1) = 13^{R_COLL} = {d_coll};  U = P_d f, V = P_d e")
    log(f"    I_(r-1) = 0 (exact, three routes) => UV = 0 a.e.  and")
    log(f"    int (P_13 U)(P_13 V) = I_r = {Ivals[R_COLL]} > 0:")
    log("    Sections 1-2 of THM-2471 apply with (m,I) = (d, I_r); note")
    log(f"    nu13(d*q) = lambda + r - 1 = {LAMBDA + R_COLL - 1} for the (27)"
        " landing valuation.")
    log("    Deep-root sidecar at the DEEPEST speed C = c_3:")
    CDEEP = C3
    delta = gcd(CDEEP, d_coll)
    C0 = CDEEP // delta
    d0 = d_coll // delta
    log(f"      C = {CDEEP} = 2*13^5,  delta = gcd(C,d) = {delta},"
        f"  C_0 = {C0},  d_0 = {d0}")
    if d0 > 1:
        branch = "r>=6-type: sheet residue a mod d_0 indispensable"
        require(R_COLL >= 6, "d_0>1 must mean r>=6 in the strict profile")
    else:
        h = CDEEP // d_coll
        log(f"      d_0 = 1, h = C/d = {h} = 2*13^{nu13(h)}")
        if h % 13 == 0:
            branch = "BASE-ONLY: 13 | h, deep phase independent of u"
            require(R_COLL <= 4, "13|h must mean r<=4 in the strict profile")
        else:
            branch = "unique affine invariant theta = t - h*u"
            require(R_COLL == 5, "13 coprime h must mean r=5")
    log(f"      branch: {branch}")
    log(f"    exact agreement with the THM-2471 (47) strict-profile table at"
        f" r = {R_COLL}.")
    log("")

    log("[9] DECISION (task branch on r)")
    if R_COLL == 5:
        log("    r = 5: Stage 2 would proceed -- NOT the case here.")
    elif R_COLL >= 6:
        log("    r >= 6: DECIDED NEGATIVE structurally -- NOT the case here.")
    else:
        log(f"    r = {R_COLL} <= 4  =>  DEGENERATE BASE-ONLY SUB-CASE; STOP.")
        log("")
        log("    Content of the outcome (exact scope):")
        log(f"    * The canonical typed row's first owner-normalized collision"
            f" is at depth")
        log(f"      r = {R_COLL}: d = 13^3 = 2197 (= the middle blocker scale"
            f" c_2; the collision")
        log("      depth coincides with nu13(c_2) = 3, two levels above the"
            " deepest leg).")
        log("    * By THM-2471 (44)-(47), at this collision the deep root of"
            " the c_3 probe")
        log("      is BASE-ONLY: h = C/d = 338 = 2*13^2 is divisible by 13,"
            " so the deep")
        log("      phase C*X_(u,a) = (y+u)*h/13 + C_0*a ... reduces mod 1 to a"
            " function")
        log("      of the base y alone -- it carries NO collision-root (u)"
            " information")
        log("      and NO sheet (a) information.")
        log("    * Stage 2 (the THM-2512 Sec. 5 bridge test in the form"
            " prescribed here,")
        log("      a toothpick contraction with v slaved through the unique"
            " affine")
        log("      invariant theta = t - 2u) requires r = 5; at r = 3 the"
            " slaving map is")
        log("      undefined (there is no u-dependence to slave: theta would"
            " degenerate")
        log("      to t itself).  The two-stage test therefore terminates at"
            " Stage 1.")
        log("    * This is NOT the r >= 6 structural negative (no claim that"
            " an")
        log("      ancestry-point-independent intertwiner fails to exist), and"
            " NOT the")
        log("      r = 5 live case.  It is the degenerate sub-case: on THIS"
            " row the deep")
        log("      root is trivially ancestry-point-independent because it is"
            " base-only;")
        log("      any lawful owner/deep coupling must run through the base"
            " coordinate y")
        log("      alone, and the u-channel pairing against the collision"
            " colours J(k)")
        log("      carries zero deep-side coupling at the actual first"
            " collision.")
    log("")

    log("[10] scope notes (MISTAKE-281 discipline)")
    log("    * TYPED row (THM-2309 (25)); NOT an asserted scalar cover; no")
    log("      physical current, no row exclusion; LRC(14) remains OPEN.")
    log("    * Every I_s is one integral of a triple product at one common")
    log("      base point u (no separate marginalizations were multiplied).")
    log("    * No intertwiner was constructed and none is claimed; the")
    log("      base-only statement is THM-2471 (45)-(47) arithmetic plus the")
    log("      exact computed r, nothing more.")
    log("    * All decisions are integer/Fraction-exact on the T_DEN, NN and")
    log("      13^max(0,s-5)*T_DEN grids; floats appear only in ~displays.")
    log("")
    log(f"{elapsed()} all checks passed")


if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    main()
