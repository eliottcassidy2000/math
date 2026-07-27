#!/usr/bin/env python3
"""THM-2512 replica dichotomy decided on one concrete live row (exact).

TASK (cheapest decisive step of the live-row transplant analysis): on the
canonical typed row THM-2309 (25), for ONE fixed word/parity clock class rho,
compute the lawful THM-2449 response table A^R EXACTLY at distinct lawful
clocks in the class, extract M_rho and C_rho from the exact clock law
A^R = M_rho + C_rho/R (THM-2449 (31), THM-2512 (26)), and run the two finite
additive tests of THM-2449 (33) / THM-2512 (27)-(28):

    d_(M_rho) = 0 ?     and     d_(C_rho) = 0 ?

where d_A(s,ell) = I_(ell,s) is the transposed doubly centred ANOVA
interaction (THM-2512 (1)-(5)).  By THM-2512 (17)/(27)-(28): if either is
nonzero, ALL 5,184 primitive cut coefficients fire at every sufficiently
large lawful clock in the class (at most one exception); if both vanish,
every clock in the class is on the replica branch (typed no-go for this
host packet; cf. THM-2456's uniform-offset hostile which keeps the replica
branch physically realizable).

ROW (THM-2309 (25), conventions of 04-computation/lrc14_169_twist_variance_
opus_20260726.py):  w = (H,q1..q5,c1,c2,c3) = (1,14,27,40,53,66,13,2197,742586),
guard H=1 (odd 13-unit), valuation profile nu13(c1,c2,c3) = (1,3,5),
owner j=1 (c_1=13), targets a=c_2=2197, b=c_3=742586 (deepest comb leg).
Danger d(y) = 1_{||y||<1/14}; safe g = 1-d (THM-2407 (4)); guard window
||H x|| <= 1/7 removed (C_H = {||Hx||>1/7}).

LAWFUL TABLE (THM-2449 Sec. 4, eqs (14)-(16)):
    Delta_r(x)     = d(c_3 x - r/13)                     (deep probe)
    d_(j,ell)(x)   = d(c_1 x - ell/7)                    (7 septimal source phases)
    U_(s,0)(x)     = source-deleted present packet: product of the eight
                     retained safe factors (guard-safe, five unit safes,
                     c_2-safe, c_3-safe) with the c_1 factor DELETED
                     (THM-2407 Sec. 1), shifted by the lawful target action
                     s*eta at t=0, with the THM-2403/THM-2309 dipole
                     eta = e_{c_2} - e_{q_1}   (star+graft: u0=q_5, graft
                     q_1 -> a; graft q_2 -> b carries t and is unshifted at
                     t=0).  Convention: a factor v with covector coefficient
                     theta_v is the safe complement of d(v x - theta_v s/13);
                     so c_2-safe = 1 - d(c_2 x - s/13) and
                     q_1-safe = 1 - d(q_1 x + s/13).
    word           = sigma = {a} (same stratum as THM-2541; first positive-
                     measure delayed word in THM-2305's fixed order):
                     Q_{1,{a}} = A_0 n D_{c_2} \\ (D_{c_1} u D_{c_3}).
                     Delete its source-safe factor (the \\ D_{c_1} part):
                     T_a = A_0 n D_{c_2} \\ D_{c_3};
                     Q^eps_ell(y) = T_a(y) * (1 - d(c_1 y - eps*ell/7))
                     (THM-2449 (15), the mandatory THM-2409 skew-diagonal).
    H^R_ell(r,s,0) = int_T U_(s,0)(x) d_(1,ell)(x) Delta_r(x) Q^eps_ell(Rx) dx
    A^R_(ell,s)    = sum_r H^R_ell(r,s,0)
                   = int_T U_(s,0) d_(1,ell) (2 - d(13 c_3 x)) Q^eps_ell(Rx) dx,
    using the exact root count sum_r Delta_r = 2 - d(13 c_3 x)
    (THM-2403 (26)); the r=0 term vanishes because the c_3-safe factor is
    retained (THM-2449 (17)).

CLOCK CLASS (THM-2449 (25)-(31)):  R = 13^k, eps = R mod 7 = (-1)^k.  All
present-side sets live on the grid Z/T_DEN with T_DEN = 182*lcm(w)
= 13^K * D_0, K = 6, gcd(D_0,13)=1.  The exact covariance (28) applies to
clocks R = 13^K N, N ≡ N' (mod D_0); hence the guaranteed class is
    rho:  sigma = {a},  eps = +1,  k ≡ 6 (mod ord_{D_0}(13)),  k >= 6.
We compute A^R at the two distinct lawful clocks k1 = 6 and
k2 = 6 + ord_{D_0}(13) IN THE SAME CLASS, extract (M_rho, C_rho) from the
two-point linear system, then VERIFY the law at a third same-class clock
k3 = 6 + 2*ord, verify M_rho against the exact product formula
M_rho = sum_r mu(F n Delta_r) * mu(Q^eps_ell) implied by (28), and report
whether the law already extends down to the sub-threshold clock k = 2
(the THM-2541 clock; not guaranteed by (28) since k < K).

EXACTNESS:  every set is a finite union of half-open intervals with integer
endpoints on the T_DEN grid; every mass is a Fraction; the dilated overlap
uses the exact prefix identity
    R*T*I_R(E,Q) = L_Q * (R*lenE - sum_i (rB_i - rA_i))/T
                   + sum_i (Phi(rB_i) - Phi(rA_i)),
    rX = (R mod T)*X mod T,  Phi = cumulative Q-length,
which makes astronomically large clocks (13^16386) exactly computable.
No floats anywhere in the decision.

SCALAR-COVER CAVEAT (THM-2541 scope guardrail): the row is typed, NOT an
asserted scalar cover, so THM-2449 (18)'s anchor A^R_(ell,0)=a_R 1_(ell=0)
is not guaranteed here; its empirical status is logged.  The ANOVA
interaction tests (THM-2512 (1)-(5)) need no anchor; only the identification
of the zero branch with the literal replica FORM (10) does.

Script: 04-computation/lrc14_replica_dichotomy_typed_row_opus_20260727.py
Output: 05-knowledge/results/lrc14_replica_dichotomy_typed_row_opus_20260727.out
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
# Row data and grid
# --------------------------------------------------------------------------
W = (1, 14, 27, 40, 53, 66, 13, 13**3, 2 * 13**5)
GUARD, OWNER, TA, TB = 0, 6, 7, 8
UNIT_IDX = (1, 2, 3, 4, 5)
C1, C2, C3 = W[OWNER], W[TA], W[TB]

LCM_W = 1
for v in W:
    LCM_W = LCM_W * v // gcd(LCM_W, v)
T = 182 * LCM_W
require(T == 297836897838480, "T_DEN mismatch vs row-conventions script")


def nu13(n):
    v = 0
    while n % 13 == 0:
        n //= 13
        v += 1
    return v


K13 = nu13(T)
require(K13 == 6, "13-adic valuation of T must be 6")
D0 = T // 13**6
require(D0 % 13 != 0, "D0 must be a 13-unit")

# multiplicative order of 13 mod D0  (THM-2449 (30))
x, ORD = 13 % D0, 1
while x != 1:
    x = x * 13 % D0
    ORD += 1
    require(ORD < 10**7, "order search overflow")

# --------------------------------------------------------------------------
# Interval machinery (integer endpoints on the T grid; half-open [a,b))
# --------------------------------------------------------------------------
def make_comb(v, PD, lo, hi):
    """Sorted disjoint intervals of  U_n ((lo+PD n)/(PD v),(hi+PD n)/(PD v))
    intersected with [0,1), on the T grid.  Requires 0 < hi-lo < PD."""
    require(T % (PD * v) == 0, f"grid does not resolve comb (v={v},PD={PD})")
    U = T // (PD * v)
    wlen = (hi - lo) * U
    base = (lo % PD) * U
    out = []
    for n in range(v):
        s = base + n * PD * U
        e = s + wlen
        if e <= T:
            out.append((s, e))
        else:
            out.append((s, T))
            out.append((0, e - T))
    out.sort()
    return out


def subtract_comb(iv, wsp, PD, lo, hi):
    """Subtract the periodic windows ((lo+PD*n)/(PD*wsp),(hi+PD*n)/(PD*wsp))."""
    require(T % (PD * wsp) == 0, f"grid does not resolve comb (v={wsp},PD={PD})")
    U = T // (PD * wsp)
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


def intersect_comb(iv, wsp, PD, lo, hi):
    """Intersect with the windows = subtract their complement windows."""
    return subtract_comb(iv, wsp, PD, hi, lo + PD)


def check_intervals(iv):
    last = -1
    for a, b in iv:
        require(0 <= a < b <= T and a >= last, "interval list corrupt")
        last = b
    return sum(b - a for a, b in iv)


# --------------------------------------------------------------------------
# Present-side builders
# --------------------------------------------------------------------------
def build_F(ell, s):
    """F_(ell,s) = U_(s,0) n {d(c_1 x - ell/7)}  as an interval list.
    13-shift windows: theta/13 center -> PD=182 numerators (14*theta -+ 13).
    7-shift windows:  ell/7 center   -> PD=182 numerators (26*ell  -+ 13)."""
    iv = make_comb(C1, 182, 26 * ell - 13, 26 * ell + 13)   # d_(1,ell)
    iv = subtract_comb(iv, W[GUARD], 91, -13, 13)           # guard-safe ||Hx||>1/7
    iv = subtract_comb(iv, W[1], 182, -14 * s - 13, -14 * s + 13)  # q1-safe, theta=-s
    iv = subtract_comb(iv, W[2], 182, -13, 13)              # q2-safe (t=0)
    iv = subtract_comb(iv, W[3], 182, -13, 13)              # q3-safe
    iv = subtract_comb(iv, W[4], 182, -13, 13)              # q4-safe
    iv = subtract_comb(iv, W[5], 182, -13, 13)              # q5-safe (= restored q_*)
    iv = subtract_comb(iv, C2, 182, 14 * s - 13, 14 * s + 13)      # c2-safe, theta=+s
    iv = subtract_comb(iv, C3, 182, -13, 13)                # c3-safe (t=0)
    return iv


def build_FD(F_iv):
    """F n D_{13 c_3},  D = {||13 c_3 x|| < 1/14}  (double-cover deficiency set:
    sum_r Delta_r = 2 - 1_D, THM-2403 (26))."""
    return intersect_comb(F_iv, 13 * C3, 14, -1, 1)


def build_word_Ta():
    """T_a = word Q_{1,{a}} with its source-safe factor deleted:
    C_H n (all unit safes) n D_{c_2} n c_3-safe."""
    iv = make_comb(C2, 182, -13, 13)                        # D_{c_2}
    iv = subtract_comb(iv, W[GUARD], 91, -13, 13)
    for i in UNIT_IDX:
        iv = subtract_comb(iv, W[i], 182, -13, 13)
    iv = subtract_comb(iv, C3, 182, -13, 13)
    return iv


# --------------------------------------------------------------------------
# Exact dilated overlap via the prefix identity
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
    Tl = T
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


def IR_from_acc(R, lenE, LQ, acc_r, acc_p):
    """Exact I_R(E,Q) = (LQ*floor-sum + Phi-sum)/(R*T)."""
    num = R * lenE - acc_r
    require(num % T == 0, "floor-count not integral")
    return Fraction(LQ * (num // T) + acc_p, R * T)


# --------------------------------------------------------------------------
# ANOVA interaction (THM-2512 (1)-(5))
# --------------------------------------------------------------------------
def interaction(X):
    """X: 7x13 Fractions.  Return I with I[ell][s]=X-R_ell-C_s+mu."""
    R7 = [sum(row, Fraction(0)) / 13 for row in X]
    C13 = [sum(X[l][s] for l in range(7)) / 7 for s in range(13)]
    mu = sum(R7, Fraction(0)) / 7
    return [[X[l][s] - R7[l] - C13[s] + mu for s in range(13)] for l in range(7)]


def lcm2(a, b):
    return a * b // gcd(a, b)


# ==========================================================================
def main():
    log("THM-2512 replica dichotomy on the canonical typed row -- exact decision")
    log("script=04-computation/lrc14_replica_dichotomy_typed_row_opus_20260727.py")
    log("machine=opus  date=2026-07-27")
    log("")
    log("[1] row and lawful-table conventions (logged choices)")
    log(f"    w=(H,q1..q5,c1,c2,c3)={W}   profile nu13(c1,c2,c3)=(1,3,5)")
    log(f"    owner j=1 (c_1={C1}); targets a=c_2={C2}, b=c_3={C3}")
    log("    dipole (THM-2403 (16), THM-2309 star+graft u0=q5, q1->a, q2->b):")
    log("      eta = e_{c2} - e_{q1};  at t=0 only s*eta acts;")
    log("      factor shift convention: safe(v x - theta_v s/13),")
    log("      i.e. c2-safe = 1-d(c2 x - s/13),  q1-safe = 1-d(q1 x + s/13)")
    log("    U_(s,0) = 8 retained safe factors, c_1 factor DELETED (THM-2407 Sec.1)")
    log("    d_(1,ell)(x) = d(13x - ell/7);  Delta_r(x) = d(c3 x - r/13)")
    log("    sum_r Delta_r = 2 - d(13 c3 x) exactly (THM-2403 (26)); r=0 term = 0")
    log("    word sigma={a} (same stratum as THM-2541, first positive in THM-2305")
    log("    order); T_a = Q_{1,{a}} minus its source-safe factor;")
    log("    Q^eps_ell(y) = T_a(y)*(1 - d(13y - eps*ell/7))  (THM-2449 (15))")
    log("")

    log("[2] clock-class arithmetic (THM-2449 (25)-(31))")
    log(f"    T = 182*lcm(w) = {T} = 13^{K13} * D0,  D0 = {D0}")
    log(f"    ord_(D0)(13) = {ORD}")
    k1, k2, k3 = 6, 6 + ORD, 6 + 2 * ORD
    log(f"    class rho: sigma={{a}}, eps=+1 (even k), k == 6 (mod {ORD}), k >= 6")
    log(f"    extraction clocks: k1={k1}, k2={k2}  (same class, guaranteed by (28))")
    log(f"    verification clocks: k3={k3} (same class), k0=2 (sub-threshold, even)")
    for kk in (2, k1, k2, k3):
        require(kk % 2 == 0, "clock parity must be even (eps=+1)")
    require(pow(13, k2 - 6, D0) == 1 and pow(13, k3 - 6, D0) == 1,
            "k2,k3 not in the class of k1")
    log("    eps = R mod 7 = (-1)^k = +1 for all four clocks: PASS")
    log("")

    # -- word side ----------------------------------------------------------
    log("[3] word sets (eps=+1)")
    Ta = build_word_Ta()
    lenTa = check_intervals(Ta)
    log(f"    T_a (source-safe factor deleted): intervals={len(Ta)}  "
        f"measure={Fraction(lenTa, T)}  ~{lenTa / T:.6f}")
    Gl = []
    Gtab = []
    LQ = []
    for ell in range(7):
        g = subtract_comb(Ta, C1, 182, 26 * ell - 13, 26 * ell + 13)
        ln = check_intervals(g)
        Gl.append(g)
        Gtab.append(make_prefix(g))
        LQ.append(ln)
        log(f"    Q^+_{ell} = T_a * (1-d(13y-{ell}/7)): intervals={len(g)}  "
            f"measure={Fraction(ln, T)}  ~{ln / T:.6f}")
    # cross-check vs THM-2541 conventions script: Q^+_0 = Q_{1,{a}}
    require(Fraction(LQ[0], T) == Fraction(143103830843, 5727632650740),
            "Q^+_0 measure != stored Q_1,{a} measure (convention drift)")
    require(len(Gl[0]) == 22478, "Q^+_0 interval count != stored 22478")
    log("    cross-check: Q^+_0 == Q_1,{a} of lrc14_169_twist_variance (measure")
    log("    143103830843/5727632650740, 22478 intervals): PASS")
    log(f"    {elapsed()}")
    log("")

    # -- clocks -------------------------------------------------------------
    R_of = {}
    Rred_of = {}
    for kk in (2, k1, k2, k3):
        R = pow(13, kk)
        R_of[kk] = R
        Rred_of[kk] = R % T
    CLKS = (2, k1, k2, k3)

    # -- main table computation --------------------------------------------
    log("[4] exact tables A^R (91 cells; F and F n D swept at 4 clocks)")
    A = {kk: [[None] * 13 for _ in range(7)] for kk in CLKS}
    lenF_tab = [[0] * 13 for _ in range(7)]
    lenFD_tab = [[0] * 13 for _ in range(7)]
    nF_tot = nFD_tot = 0
    for ell in range(7):
        starts, lens, pref = Gtab[ell]
        for s in range(13):
            F = build_F(ell, s)
            lenF = check_intervals(F)
            FD = build_FD(F)
            lenFD = sum(b - a for a, b in FD)
            lenF_tab[ell][s] = lenF
            lenFD_tab[ell][s] = lenFD
            nF_tot += len(F)
            nFD_tot += len(FD)
            for kk in CLKS:
                r1, p1 = sweep_acc(F, Rred_of[kk], starts, lens, pref)
                r2, p2 = sweep_acc(FD, Rred_of[kk], starts, lens, pref)
                IF = IR_from_acc(R_of[kk], lenF, LQ[ell], r1, p1)
                IFD = IR_from_acc(R_of[kk], lenFD, LQ[ell], r2, p2)
                A[kk][ell][s] = 2 * IF - IFD
            del F, FD
        log(f"    {elapsed()} ell={ell} done "
            f"(row |F| avg {sum(lenF_tab[ell]) // 13}/T mass, intervals so far "
            f"F:{nF_tot} FD:{nFD_tot})")
    log("")

    # -- sanity checks ------------------------------------------------------
    log("[5] sanity checks")
    # (a) full-circle control: I_R([0,1), Q) == mu(Q) at every clock
    full = [(0, T)]
    for kk in CLKS:
        r0, p0 = sweep_acc(full, Rred_of[kk], *Gtab[0])
        I0 = IR_from_acc(R_of[kk], T, LQ[0], r0, p0)
        require(I0 == Fraction(LQ[0], T), f"full-circle control failed at k={kk}")
    log("    (a) I_R(full circle, Q^+_0) == mu(Q^+_0) at all 4 clocks: PASS")
    # (b) E_1 cross-check: F_(0,0) == E_1 of the conventions script
    F00 = build_F(0, 0)
    require(Fraction(check_intervals(F00), T) == Fraction(1882176, 28589561),
            "mu(F_(0,0)) != stored E_1 measure")
    require(len(F00) == 57072, "F_(0,0) interval count != stored 57072")
    log("    (b) F_(0,0) == E_1 (measure 1882176/28589561, 57072 intervals): PASS")
    # (c) stored k=2 stratum overlap: I_169(E_1, Q_1,{a}) (word factor at ell=0)
    r0, p0 = sweep_acc(F00, Rred_of[2], *Gtab[0])
    I169 = IR_from_acc(R_of[2], check_intervals(F00), LQ[0], r0, p0)
    require(I169 == Fraction(21376087, 17907461390),
            "I_169(E_1,Q_1,{a}) != stored stratum measure")
    log("    (c) I_169(F_(0,0),Q^+_0) == stored measure(E_1 n T^-2 Q_1,{a})")
    log("        = 21376087/17907461390: PASS  (independent machinery check)")
    # (d) Delta_0 vanishing and the root-count identity on two cells
    for (el, sc) in ((0, 0), (3, 7)):
        Fc = build_F(el, sc)
        lenFc = sum(b - a for a, b in Fc)
        FDc = build_FD(Fc)
        lenFDc = sum(b - a for a, b in FDc)
        d0 = intersect_comb(Fc, C3, 182, -13, 13)
        require(sum(b - a for a, b in d0) == 0, "F n Delta_0 != empty")
        tot = 0
        sumIR = Fraction(0)
        st, ln_, pf = Gtab[el]
        for r in range(13):
            Fr = intersect_comb(Fc, C3, 182, 14 * r - 13, 14 * r + 13)
            lr = sum(b - a for a, b in Fr)
            tot += lr
            rr, pp = sweep_acc(Fr, Rred_of[k1], st, ln_, pf)
            sumIR += IR_from_acc(R_of[k1], lr, LQ[el], rr, pp)
        require(tot == 2 * lenFc - lenFDc, "sum_r mu(F n Delta_r) != 2muF-muFD")
        require(sumIR == A[k1][el][sc], "sum_r I_R != weighted table entry")
        log(f"    (d) cell (ell={el},s={sc}): F n Delta_0 = empty; "
            f"sum_r mu(F n Delta_r) == 2*mu(F)-mu(FD); sum_r I_R == A entry: PASS")
    # (e) partition of the seven source phases at s=0
    U0 = [(0, T)]
    U0 = subtract_comb(U0, W[GUARD], 91, -13, 13)
    for i in UNIT_IDX:
        U0 = subtract_comb(U0, W[i], 182, -13, 13)
    U0 = subtract_comb(U0, C2, 182, -13, 13)
    U0 = subtract_comb(U0, C3, 182, -13, 13)
    require(sum(lenF_tab[ell][0] for ell in range(7)) == check_intervals(U0),
            "sum_ell mu(F_(ell,0)) != mu(U_(0,0))")
    log("    (e) sum_ell mu(F_(ell,0)) == mu(U_(0,0)) (septimal partition): PASS")
    log(f"    {elapsed()}")
    log("")

    # -- extraction ---------------------------------------------------------
    log("[6] extraction of M_rho, C_rho from the two same-class clocks (k1,k2)")
    R1, R2, R3 = R_of[k1], R_of[k2], R_of[k3]
    M = [[(R1 * A[k1][l][s] - R2 * A[k2][l][s]) / (R1 - R2)
          for s in range(13)] for l in range(7)]
    Cm = [[R1 * (A[k1][l][s] - M[l][s]) for s in range(13)] for l in range(7)]
    okM = all(M[l][s] == (2 * Fraction(lenF_tab[l][s], T)
                          - Fraction(lenFD_tab[l][s], T)) * Fraction(LQ[l], T)
              for l in range(7) for s in range(13))
    require(okM, "M_rho != product-measure formula (THM-2449 (28) limit)")
    log("    M_rho == [sum_r mu(F n Delta_r)] * mu(Q^+_ell) entrywise: PASS")
    ok3 = all(R3 * (A[k3][l][s] - M[l][s]) == Cm[l][s]
              for l in range(7) for s in range(13))
    require(ok3, "third same-class clock violates A = M + C/R")
    log(f"    third same-class clock k3={k3}: A^R == M_rho + C_rho/R exact: PASS")
    ok2 = all(A[2][l][s] == M[l][s] + Cm[l][s] / R_of[2]
              for l in range(7) for s in range(13))
    log(f"    sub-threshold clock k=2 (R=169, NOT guaranteed by (28)): "
        f"law holds exactly = {ok2}")
    log("")

    denM = 1
    denC = 1
    l1M = Fraction(0)
    l1C = Fraction(0)
    for l in range(7):
        for s in range(13):
            denM = lcm2(denM, M[l][s].denominator)
            denC = lcm2(denC, Cm[l][s].denominator)
            l1M += abs(M[l][s])
            l1C += abs(Cm[l][s])
    log("[7] exact class matrices")
    log(f"    M_rho: common denominator {denM}")
    log(f"           l1 norm = {l1M}  ~{float(l1M):.6f}")
    log(f"    C_rho: common denominator {denC}")
    log(f"           l1 norm = {l1C}  ~{float(l1C):.6f}")
    log(f"    lcm(den M, den C) = {lcm2(denM, denC)}")
    log("    M_rho row ell=0 (s=0..3): "
        + ", ".join(str(M[0][s]) for s in range(4)))
    log("    C_rho row ell=0 (s=0..3): "
        + ", ".join(str(Cm[0][s]) for s in range(4)))
    log("")

    # -- anchor status ------------------------------------------------------
    log("[8] owner-anchor status at the untwisted target column s=0 (THM-2449 (18);")
    log("    NOT guaranteed here -- typed row is not an asserted scalar cover)")
    for kk in (k1,):
        col = [A[kk][l][0] for l in range(7)]
        log(f"    k={kk}: A_(ell,0) = "
            + ", ".join(f"{c}" if c.denominator < 10**12 else f"~{float(c):.3e}"
                        for c in col))
    anchor_M = all(M[l][0] == 0 for l in range(1, 7)) and M[0][0] > 0
    anchor_C = all(Cm[l][0] == 0 for l in range(1, 7))
    log(f"    M_rho anchor (M_(0,0)>0, M_(ell,0)=0 for ell!=0): {anchor_M}")
    log(f"    C_rho anchor column zero off owner: {anchor_C}")
    log(f"    M_(0,0) = {M[0][0]}  (positive owner overlap at large clocks: "
        f"{M[0][0] > 0})")
    log("")

    # -- DECISION -----------------------------------------------------------
    log("[9] DECISION: the two finite additive tests (THM-2449 (33), THM-2512 (27))")
    IM = interaction(M)
    IC = interaction(Cm)
    wM = [(l, s) for l in range(7) for s in range(13) if IM[l][s] != 0]
    wC = [(l, s) for l in range(7) for s in range(13) if IC[l][s] != 0]
    log(f"    d_(M_rho) == 0 ?  {'YES' if not wM else 'NO'}   "
        f"(nonzero interaction entries: {len(wM)}/91)")
    log(f"    d_(C_rho) == 0 ?  {'YES' if not wC else 'NO'}   "
        f"(nonzero interaction entries: {len(wC)}/91)")
    if wM:
        l, s = wM[0]
        log(f"    witness d_(M_rho): (ell,s)=({l},{s})  "
            f"I_(ell,s) = {IM[l][s]}  ~{float(IM[l][s]):.3e}")
        mn = min(abs(IM[l][s]) for (l, s) in wM)
        log(f"    min nonzero |I(M)| = {mn}  ~{float(mn):.3e}")
    if wC:
        l, s = wC[0]
        vv = IC[l][s]
        log(f"    witness d_(C_rho): (ell,s)=({l},{s})  "
            f"I_(ell,s) = {vv}  ~{float(vv):.3e}")
    log("")
    if wM or wC:
        log("    VERDICT: NON-REPLICA BRANCH on the canonical typed row.")
        log("    By THM-2512 (27)-(28) with the exact class law A^R=M_rho+C_rho/R:")
        log("    at every sufficiently large lawful clock in the class rho")
        log(f"    (k == 6 mod {ORD}, k >= 6, eps=+1, sigma={{a}}), with at most ONE")
        log("    exceptional clock, d_(A^R) != 0, hence ALL 5,184 primitive cut")
        log("    coefficients Psi^A_(tau,a0)(alpha,beta) fire (THM-2512 (17)-(19)),")
        log("    and at least 294 toothpick components / 3,528 Theta entries are")
        log("    nonzero (THM-2512 (20)).  The first three transplant laws hold on")
        log("    a live-branch array.")
        log("    The missing l1/denominator floor of THM-2512 Sec. 2 for THIS row:")
        log(f"      common denominator of M_rho:      {denM}")
        log(f"      common denominator of C_rho:      {denC}")
        log(f"      common denominator (both, lcm):   {lcm2(denM, denC)}")
        log(f"      l1(M_rho) = {l1M}")
        log(f"      l1(C_rho) = {l1C}")
        log("    (Any nonzero interaction entry of Den*M_rho resp. Den*C_rho is an")
        log("    integer multiple of 1/91, so |nonzero d| >= 1/(91*Den).)")
        if wM and not anchor_M:
            log("    NOTE: the delta anchor (THM-2449 (18)) did NOT hold exactly on")
            log("    this typed row (see [8]); the non-replica conclusion (17)/(19)")
            log("    needs only rationality + d_A != 0, so it is unaffected; only")
            log("    the replica-FORM normalization (10) would have needed the")
            log("    anchor.")
    else:
        log("    VERDICT: REPLICA BRANCH -- typed no-go for host packet 1.")
        log("    Both finite matrices have zero interaction, so by THM-2449 (33)")
        log("    every clock in the class is mixed-zero: the canonical typed row's")
        log("    lawful owner table sits on the delta-plus-six-replicas boundary")
        log("    (physically realizable per THM-2456's uniform-offset hostile).")
        log("    The whole THM-2512 cut bundle vanishes identically on this class.")
        log("    REDIRECT: the THM-2535/THM-2542 clock-boundary system and its")
        log("    bidegree-(1,1) 2-cell.")
    log("")
    log("[10] scope notes")
    log("    * Typed row, not an asserted scalar cover; no scalar row is removed;")
    log("      LRC(14) remains open.  The decision is about THIS row's lawful")
    log("      THM-2449 response table on the class rho only.")
    log("    * The tables are exact rationals (Fraction interval intersection on")
    log("      the T grid); no floats enter any decision.")
    log("    * Clock-law scope: guaranteed for k >= 6 in the class by THM-2449")
    log("      (25)-(28); the k=2 status above is reported, not assumed.")
    log("")
    log(f"{elapsed()} all checks passed")


if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    main()
