#!/usr/bin/env python3
"""THM-2368 Sec. 8 cheapest decisive step, executed as THM-2365 eq (19b):
the first-hostile successor-mass drift test on the canonical typed row.

OBJECT.  Canonical typed row THM-2309 (25) (same row/word THM-2541 certified):

    w = (H, q1..q5, c1, c2, c3) = (1, 14, 27, 40, 53, 66, 13, 2197, 742586),
    owner j = 1 (c1 = 13, depth 1), targets a = c2 = 13^3, c = c3 = 2*13^5
    (deepest), word sigma = {a} (first positive in THM-2305 order), clock
    k = 2, R = 13^k = 169.

Source-refined owner packet (THM-2461 eq (9), from THM-2306; THM-2461 eq (7)
proves the packet equals its own source-atom refinement -- E_1 sits in the
single source atom tau_src, so no further refinement factor exists):

    f_sigma(y) = 1_{Q_{1,{a}}}(y) * (P^2 1_{E_1})(y),
    x-side representation  F(x) = 1_{E_1}(x) * 1_{Q_{1,{a}}}(169 x),

which is exactly THM-2365 Sec. 7's delayed-word packet F_R = E * Q(R x) at
R = 169.  Lawful target action (THM-2365 (5), minus-shift convention;
concrete lawful covector basis of THM-2365 Sec. 1 / THM-2367 (22), the same
grafts q1 -> a=c2, q2 -> c=c3 as the certified twist bank [THM-2541]):

    eta = e_{c2} - e_{q1},    ell = e_{c3} - e_{q2},
    E_{s,t}(x) = prod_i I_i(w_i x - (s eta_i + t ell_i)/13),
    F_{s,t}(x) = E_{s,t}(x) * 1_{Q_{1,{a}}}(169 x)      (word twist-neutral:
                 the word translations R*shift/13 = 13*shift are integers,
                 THM-2365 (5) / THM-2334 (44)).

Deep probe and tensors (THM-2365 (6)-(7), c = c3 the deepest blocker):

    Delta_r(x) = d(c3 x - r/13),        d(y) = 1_{||y|| < 1/14},
    H_F(r,s,t) = int F_{s,t}(x) Delta_r(x) dx     (packet / owner-loop tensor,
                                                   = THM-2365 H_R at R=169),
    H_E(r,s,t) = int E_{s,t}(x) Delta_r(x) dx     (bare-owner tensor, the
                                                   THM-2365 (32) object).

DECISION (THM-2365 (19b)): the 169 successor masses

    S(s,t) = sum_r H(r,s,t) = int F_{s,t}(x) (2 - d(13 c3 x)) dx

are exact rationals.  If they are NOT all equal, the tensor is not circulant
(if H(r,s,t) = G(r-t) then S is constant), hence by THM-2365 (17)-(18) the
drift energy D_H > 0 -- an exact certificate.  Constancy of S would only be
a first hostile test (THM-2365 after (19b)); we therefore ALSO compute the
full 13^3 tensors and test the rigid circulant law H(r,s,t) = G(r-t)
directly (THM-2365 (32) escalation), plus the exact drift energies
D_H = (1/13^3) sum |H - PH|^2 with PH the (16) target-action average,
cross-checked against the (19a) group-average identity.

Zero shape (THM-2367 Sec. 4): a zero outcome would be H(r,s,t) = A*J(r-t),
i.e. the sharp circulant object -- a structured law, not a failure.

SCOPE (THM-2541 guardrail): the row is TYPED, not asserted to be a scalar
cover.  A positive drift here is the first drift-positive instance feeding
the proved chain THM-2365 (31) + THM-2466 (32) + THM-2471 (48) on THIS
row's packet; it is NOT a universal bypass, removes no scalar row, and
leaves LRC(14) open.

Everything decision-bearing is exact integer/Fraction arithmetic on the
common denominators T_DEN (bare) and NN = 169*T_DEN (packet); no floats
enter any decision.  Interval machinery is reused verbatim from the audited
04-computation/lrc14_169_twist_variance_opus_20260726.py (in_comb,
subtract_comb, build_set, and the x_sweep walking pattern).

Script: 04-computation/lrc14_owner_loop_drift_typed_row_opus_20260727.py
Output: 05-knowledge/results/lrc14_owner_loop_drift_typed_row_opus_20260727.out
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
        raise AssertionError(msg)


# --------------------------------------------------------------------------
# Row data and scales (identical to lrc14_169_twist_variance_opus_20260726.py)
# --------------------------------------------------------------------------
W = (1, 14, 27, 40, 53, 66, 13, 13**3, 2 * 13**5)
GUARD, OWNER, TA, TB = 0, 6, 7, 8
UNIT_IDX = (1, 2, 3, 4, 5)
K_CLOCK = 2
RDIL = 13**K_CLOCK  # 169
C3 = W[TB]


def nu13(n):
    v = 0
    while n % 13 == 0:
        n //= 13
        v += 1
    return v


require(W[GUARD] % 2 == 1 and W[GUARD] % 13 != 0, "guard must be odd 13-unit")
require(all(W[i] % 13 != 0 for i in UNIT_IDX), "units must be 13-units")
require(len({W[i] for i in UNIT_IDX}) == 5, "units distinct")
PROFILE = (nu13(W[OWNER]), nu13(W[TA]), nu13(W[TB]))
require(PROFILE == (1, 3, 5), f"profile {PROFILE}")
require(RDIL % 13 == 0, "13 | R: word translations are target-neutral")

LCM_W = 1
for v in W:
    LCM_W = LCM_W * v // gcd(LCM_W, v)
T_DEN = 182 * LCM_W          # common denominator of all base breakpoints
NN = RDIL * T_DEN            # common denominator after the k=2 word pullback
require(all(T_DEN % (182 * v) == 0 for v in W), "T_DEN divisibility")
require(T_DEN == 297836897838480, "T_DEN anchor (twist-variance artifact)")
require(NN == 50334435734703120, "NN anchor (twist-variance artifact)")

# Deep-probe / successor comb geometry.
#   Probe Delta_r: windows centered at (13n+r)/(13*c3), half-width 1/(14*c3).
#   Successor d(13 c3 x): windows centered at n/(13*c3), half-width 1/(14*13*c3).
# On the T_DEN scale the center spacing is PP_T and both half-widths are
# integers; on the NN scale everything is multiplied by 169.
L0 = LCM_W // C3
require(LCM_W % C3 == 0, "c3 | LCM_W")
PP_T = T_DEN // (13 * C3)            # = 14*L0: spacing of window centers
require(PP_T == 14 * L0 and T_DEN == 13 * C3 * PP_T, "comb grid exact")
HPROBE_T = 13 * L0                   # probe half-width  = T_DEN/(14*c3)
HSUCC_T = L0                         # successor half-width = T_DEN/(14*13*c3)
PP_N, HPROBE_N, HSUCC_N = 169 * PP_T, 169 * HPROBE_T, 169 * HSUCC_T
require((13 * C3) % 13 == 0, "center count per circle divisible by 13")


# --------------------------------------------------------------------------
# Exact interval machinery (verbatim from the audited twist-variance script)
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
# Comb accumulation: probe windows [c*PP - hp, c*PP + hp] attributed to
# r = c mod 13, successor windows [c*PP - hs, c*PP + hs] (same centers).
# Since hs < hp and window width 2*hp = (13/7)*PP < 13*PP, scanning centers
# c with c*PP + hp > a and c*PP - hp < b covers every contributing window.
# --------------------------------------------------------------------------
def accumulate_interval(a, b, PP, hp, hs, probe_acc):
    """Add overlaps of [a,b] with the 13 probe combs and return successor
    overlap contribution.  Pure integers."""
    succ = 0
    c = (a - hp) // PP + 1
    cend = (b + hp) // PP
    while c <= cend:
        ctr = c * PP
        lo = ctr - hp
        hi = ctr + hp
        ov = (b if b < hi else hi) - (a if a > lo else lo)
        if ov > 0:
            probe_acc[c % 13] += ov
            lo2 = ctr - hs
            hi2 = ctr + hs
            ov2 = (b if b < hi2 else hi2) - (a if a > lo2 else lo2)
            if ov2 > 0:
                succ += ov2
        c += 1
    return succ


def selftest_combs():
    """Validate accumulate_interval against an independent unit-cell count.
    Windows have integer endpoints, so the overlap with an integer interval
    equals the number of unit cells [x, x+1) whose (half-integer) midpoint
    lies inside a window -- and a midpoint never hits a window boundary."""
    PP, hp, hs = 14, 13, 1          # same 14/13/1 shape as the real combs
    top = 13 * 7 * PP               # one full circle: 91 centers, 91 % 13 == 0
    require((top // PP) % 13 == 0, "selftest circle must close mod 13")
    cases = ((0, 37), (5, 200), (100, 1274), (0, top), (1270, top),
             (13, 15), (27, 29), (170, 171), (0, 1), (top - 1, top))
    for (a, b) in cases:
        acc = [0] * 13
        s = accumulate_interval(a, b, PP, hp, hs, acc)
        ref = [0] * 13
        rs = 0
        for x in range(a, b):
            m2 = 2 * x + 1          # midpoint*2 of cell [x, x+1)
            # nearest centers: c*PP within hp of m2/2  <=>  |m2 - 2*c*PP| < 2*hp
            for c in range((x - hp) // PP - 1, (x + 1 + hp) // PP + 2):
                d2 = m2 - 2 * c * PP
                if -2 * hp < d2 < 2 * hp:
                    ref[c % 13] += 1
                    if -2 * hs < d2 < 2 * hs:
                        rs += 1
        require(acc == ref and s == rs, f"selftest comb mismatch on {(a, b)}")
        # multiplicity identity sum_r probe = 2*len - succ on the full circle
        if (a, b) == (0, top):
            require(sum(acc) == 2 * (b - a) - s, "selftest identity")
    log("    comb accumulator self-test vs independent unit-cell count: PASS")


# --------------------------------------------------------------------------
# Per-twist computation
# --------------------------------------------------------------------------
def twist_shift_vector(s, t):
    """THM-2365 (5): factor i is shifted by -(s*eta_i + t*ell_i)/13, and the
    interval machinery uses the +e_i/13 convention, so e_i = -(s eta + t ell).
    eta = e_{TA} - e_{q1} (index 7 minus index 1),
    ell = e_{TB} - e_{q2} (index 8 minus index 2)."""
    e = [0] * 9
    e[1] = s % 13            # graft q1 moves by +s/13
    e[2] = t % 13            # graft q2 moves by +t/13
    e[TA] = (-s) % 13        # target blocker c2 moves by -s/13
    e[TB] = (-t) % 13        # target blocker c3 moves by -t/13
    return tuple(e)


def process_twist(s, t, Q_iv, Q_starts, nQ):
    """Return (lenE, HE[13], succE, lenF, HF[13], succF) integers.
    HE/succE on the T_DEN scale, HF/succF on the NN scale."""
    E_iv = build_set(PAT_E, twist_shift_vector(s, t))
    lenE = 0
    HE = [0] * 13
    succE = 0
    lenF = 0
    HF = [0] * 13
    succF = 0
    ppn, hpn, hsn = PP_N, HPROBE_N, HSUCC_N
    ppt, hpt, hst = PP_T, HPROBE_T, HSUCC_T
    tden = T_DEN
    for A, B in E_iv:
        lenE += B - A
        succE += accumulate_interval(A, B, ppt, hpt, hst, HE)
        # word refinement: walk Q intervals across the 169x image (x_sweep
        # pattern); piece [LAoff+lo, LAoff+hi] on the NN scale.
        LA = 169 * A
        sA = LA % tden
        span = 169 * (B - A)
        require(span < tden, "E interval too long for single-wrap walk")
        sE = sA + span
        LAoff = LA - sA
        idx = bisect_right(Q_starts, sA) - 1
        off = 0
        if idx < 0:
            idx = nQ - 1
            off = -tden
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
                    a = LAoff + lo
                    b = LAoff + hi
                    lenF += b - a
                    succF += accumulate_interval(a, b, ppn, hpn, hsn, HF)
            idx += 1
            if idx == nQ:
                idx = 0
                off += tden
    return lenE, HE, succE, lenF, HF, succF


# --------------------------------------------------------------------------
# Tensor post-processing (THM-2365 (16)-(19a), (32))
# --------------------------------------------------------------------------
def tensor_report(name, tab, scale_den, unit_name):
    """tab[(r,s,t)] integer masses on scale scale_den.  Returns dict of
    results; prints the full analysis."""
    log(f"    --- tensor {name} (integer masses in units of 1/{unit_name}) ---")
    # nonnegativity + diagonal-plane zero (THM-2365 (8)-(9))
    require(all(v >= 0 for v in tab.values()), f"{name}: negative mass")
    diag_ok = all(tab[(t, s, t)] == 0 for s in range(13) for t in range(13))
    log(f"    diagonal-plane zero H(t,s,t)=0 for all 169 (s,t): "
        f"{'PASS' if diag_ok else 'FAIL'}")
    require(diag_ok, f"{name}: THM-2365 (9) diagonal zero failed -- "
                     "sign-convention error")
    # (19b) successor masses
    S = {(s, t): sum(tab[(r, s, t)] for r in range(13))
         for s in range(13) for t in range(13)}
    s00 = S[(0, 0)]
    witnesses = [(s, t) for s in range(13) for t in range(13)
                 if S[(s, t)] != s00]
    # circulant test (THM-2365 (18)/(32))
    G = [tab[(d % 13, 0, 0)] for d in range(13)]
    viol = [(r, s, t) for r in range(13) for s in range(13) for t in range(13)
            if tab[(r, s, t)] != G[(r - t) % 13]]
    # exact drift energy (17): PH(r,s,t) = psi(r-t),
    # psi(d) = (1/169) sum_{u,v} H(d+v,u,v)
    num = [sum(tab[((d + v) % 13, u, v)] for u in range(13) for v in range(13))
           for d in range(13)]
    dnum = 0
    for r in range(13):
        for s in range(13):
            for t in range(13):
                diff = 169 * tab[(r, s, t)] - num[(r - t) % 13]
                dnum += diff * diff
    D = Fraction(dnum, 13**3 * 169**2 * scale_den**2)
    # (19a) group-average cross-check
    dnum2 = 0
    for u in range(13):
        for v in range(13):
            acc = 0
            for r in range(13):
                for s in range(13):
                    for t in range(13):
                        diff = tab[(r, s, t)] - tab[((r + v) % 13,
                                                     (s + u) % 13,
                                                     (t + v) % 13)]
                        acc += diff * diff
            dnum2 += acc
    D2 = Fraction(dnum2, 2 * 169 * 13**3 * scale_den**2)
    require(D == D2, f"{name}: (17) vs (19a) drift-energy mismatch")
    log(f"    (17)-vs-(19a) exact drift-energy identity: PASS")
    log(f"    successor masses S(s,t): {len(set(S.values()))} distinct "
        f"values among 169")
    log(f"    S nonconstancy witnesses: {len(witnesses)} / 168 nonzero twists")
    log(f"    circulant-law violations H(r,s,t) != G(r-t): {len(viol)} / 2197")
    log(f"    exact drift energy D_H = {D}")
    log(f"                          ~ {float(D):.6e}")
    return {"S": S, "witnesses": witnesses, "viol": viol, "D": D, "G": G}


def main():
    log("THM-2365 (19b) first-hostile owner-loop drift test on the canonical"
        " typed row")
    log("script=04-computation/lrc14_owner_loop_drift_typed_row_opus_20260727.py")
    log("machine=opus  date=2026-07-27")
    log("")
    log("[1] row, packet, and conventions (choices logged; files override"
        " summaries)")
    log(f"    w=(H,q1..q5,c1,c2,c3)={W}   valuation profile {PROFILE}")
    log(f"    owner j=1, word sigma={{a}} (Q_1,{{a}} = A_0 n D_c2 \\"
        f" (D_c1 u D_c3)), clock k=2, R=169")
    log("    packet: f = 1_Q P^2 1_(E_1)  (THM-2461 (9)/THM-2306);"
        " x-side F(x)=1_(E_1)(x)*1_Q(169x)")
    log("      = THM-2365 Sec.7 F_R at R=169.  THM-2461 (7): E_1 lies in the"
        " single source atom,")
    log("      so F IS the source-refined packet (no extra refinement factor"
        " exists).")
    log("    covectors (THM-2365 Sec.1 concrete basis / THM-2367 (22); grafts"
        " q1->c2, q2->c3):")
    log("      eta = e_c2 - e_q1,  ell = e_c3 - e_q2")
    log("    target action (THM-2365 (5), minus convention): factor i shifted"
        " by -(s*eta_i+t*ell_i)/13")
    log("      => c2: -s/13, q1: +s/13, c3: -t/13, q2: +t/13; word factors"
        " shifted by -13*(s,t) in Z => unchanged")
    log("    probe: Delta_r(x)=d(c3 x - r/13), c3=742586 the deepest blocker"
        " (THM-2365 (6))")
    log(f"    scales: T_DEN={T_DEN} (bare), NN=169*T_DEN={NN} (packet)")
    log(f"    comb grid: center spacing PP_T={PP_T}, probe half-width"
        f" {HPROBE_T}, successor half-width {HSUCC_T}")
    log("    tensor naming: H_F = packet (word-refined, THM-2365 H_R at"
        " R=169); H_E = bare owner (THM-2365 (32)).")
    log("      The task brief labels the packet drift 'D_(H_E)'; in THM-2365's"
        " own notation the packet tensor")
    log("      is H_R.  Both tensors are computed below, so both readings are"
        " decided.")
    log("")

    log("[0] machinery self-test")
    selftest_combs()
    log("")

    log("[2] exact interval sets and anchors")
    E0 = build_set(PAT_E, ZELL)
    lenE0 = check_intervals(E0, T_DEN)
    mu_E = Fraction(lenE0, T_DEN)
    log(f"    E_1: intervals={len(E0)}  measure={mu_E}")
    require(mu_E == Fraction(1882176, 28589561),
            "E_1 measure anchor vs twist-variance artifact")
    Q_iv = build_set(PAT_QA, ZELL)
    lenQ = check_intervals(Q_iv, T_DEN)
    log(f"    Q_1,{{a}}: intervals={len(Q_iv)}  measure={Fraction(lenQ, T_DEN)}")
    require(Fraction(lenQ, T_DEN) == Fraction(143103830843, 5727632650740),
            "Q_1,{a} measure anchor vs twist-variance artifact")
    Q_starts = [a for a, _ in Q_iv]
    nQ = len(Q_iv)
    log("    anchors vs 05-knowledge/results/lrc14_169_twist_variance_opus_"
        "20260726.out: PASS")
    log("")

    log("[3] main sweep: 169 lawful twists x (13 probes + successor), bare"
        " and packet")
    HE_tab = {}
    HF_tab = {}
    SE = {}
    SF = {}
    muE = {}
    muF = {}
    for s in range(13):
        for t in range(13):
            lenE, HE, succE, lenF, HF, succF = process_twist(
                s, t, Q_iv, Q_starts, nQ)
            # THM-2365 (19b) successor identity, both scales, exact:
            require(sum(HE) == 2 * lenE - succE,
                    f"(19b) identity fails for H_E at {(s, t)}")
            require(sum(HF) == 2 * lenF - succF,
                    f"(19b) identity fails for H_F at {(s, t)}")
            for r in range(13):
                HE_tab[(r, s, t)] = HE[r]
                HF_tab[(r, s, t)] = HF[r]
            SE[(s, t)] = sum(HE)
            SF[(s, t)] = sum(HF)
            muE[(s, t)] = lenE
            muF[(s, t)] = lenF
        log(f"    {elapsed()} s={s} done")
    require(muF[(0, 0)] == Fraction(21376087, 17907461390) * NN,
            "stratum measure anchor vs twist-variance artifact")
    log("    untwisted stratum measure mu(F_0,0) = "
        f"{Fraction(muF[(0, 0)], NN)}  (anchor PASS)")
    log("    rowwise (19b) identity sum_r H = 2*mu - successor: PASS at all"
        " 169 twists, both tensors")
    log("")

    log("[4] PACKET tensor H_F: the task's decision object")
    resF = tensor_report("H_F (source-refined owner packet, R=169)",
                         HF_tab, NN, "NN")
    log("")
    log("[5] BARE-OWNER tensor H_E: the THM-2365 (32) escalation object")
    resE = tensor_report("H_E (bare owner, word deleted)", HE_tab, T_DEN,
                         "T_DEN")
    log("")

    log("[6] DECISION (THM-2365 (19b) first-hostile test on the packet)")
    wF = resF["witnesses"]
    if wF:
        s0 = SF[(0, 0)]
        ws, wt = wF[0]
        log("    the 169 successor masses S(s,t) of the source-refined owner"
            " packet are NOT all equal.")
        log(f"    witness pair: (s,t)=(0,0) vs (s,t)=({ws},{wt})")
        log(f"      S(0,0)      = {Fraction(s0, NN)}")
        log(f"                  = {s0}/NN")
        log(f"      S({ws},{wt})"
            f"      = {Fraction(SF[(ws, wt)], NN)}")
        log(f"                  = {SF[(ws, wt)]}/NN")
        log(f"      difference  = {Fraction(SF[(ws, wt)] - s0, NN)}  (exact,"
            " nonzero)")
        log("    CERTIFICATE: nonconstant S => H_F not circulant => by"
            " THM-2365 (17)-(18)")
        log(f"    D_(H_F) = {resF['D']} > 0 (owner-loop drift POSITIVE on this"
            " row's packet).")
    else:
        log("    the 169 packet successor masses are constant -- first"
            " hostile test inconclusive")
        log("    (constancy of S is NOT a proof of zero drift); escalating"
            " per THM-2365 (32):")
        if resF["viol"]:
            r, s, t = resF["viol"][0]
            log(f"    circulant law FAILS: H_F({r},{s},{t}) = "
                f"{Fraction(HF_tab[(r, s, t)], NN)} != G_F({(r - t) % 13}) = "
                f"{Fraction(resF['G'][(r - t) % 13], NN)}")
            log(f"    => D_(H_F) = {resF['D']} > 0.")
        else:
            log("    H_F(r,s,t) = G_F(r-t) EXACTLY on all 2197 cells: the"
                " packet tensor is the rigid")
            log("    circulant object (THM-2365 (18), the (26)"
                " inverse-character line; THM-2367 Sec.4 shape).")
    log("")

    log("[7] bare-owner branch (THM-2365 (31)-(32))")
    wE = resE["witnesses"]
    if wE:
        s0 = SE[(0, 0)]
        ws, wt = wE[0]
        log("    the 169 BARE successor masses S_E(s,t) are also nonconstant.")
        log(f"    witness pair: (s,t)=(0,0) vs ({ws},{wt}):")
        log(f"      S_E(0,0)   = {Fraction(s0, T_DEN)}")
        log(f"      S_E({ws},{wt})   = {Fraction(SE[(ws, wt)], T_DEN)}")
        log(f"    => D_(H_E) = {resE['D']} > 0: BARE-OWNER drift positive."
            "  By THM-2365 (31),")
        log("    at every sufficiently large delayed clock some exact 91-unit"
            " m and X have a")
        log("    nonzero target fibre on this row.")
    else:
        if resE["viol"]:
            r, s, t = resE["viol"][0]
            log("    S_E constant but circulant law fails: "
                f"H_E({r},{s},{t}) != G_E({(r - t) % 13});"
                f" D_(H_E) = {resE['D']} > 0.")
        else:
            log("    H_E(r,s,t) = G_E(r-t) EXACTLY: the bare owner is the"
                " sharp circulant object")
            log("    (THM-2365 (30)-(31) second branch; THM-2367 Sec.4"
                " shape).  Zero here is a law,")
            log("    not a failure; the packet result above stands on its"
                " own.")
    log("")

    log("[8] exact successor-mass tables (numerators on the common"
        " denominators; rows s, cols t)")
    log("    S_F(s,t)*NN:")
    for s in range(13):
        log("      " + " ".join(str(SF[(s, t)]) for t in range(13)))
    log("    S_E(s,t)*T_DEN:")
    for s in range(13):
        log("      " + " ".join(str(SE[(s, t)]) for t in range(13)))
    log("")

    log("[9] scope (THM-2541 guardrail, unchanged)")
    log("    * The row is the canonical TYPED row THM-2309 (25): typed, NOT"
        " asserted to be a")
    log("      scalar cover.  Nothing here excludes a scalar row; the ledger"
        " stays 165 and")
    log("      LRC(14) remains OPEN.")
    log("    * A positive drift is the first drift-positive instance feeding"
        " the proved chain")
    log("      THM-2365 (31) [target landing] + THM-2466 (32) [delayed-word"
        " drift service")
    log("      retention] + THM-2471 (48) [weighted root service] on THIS"
        " row's packet -- not")
    log("      a universal bypass.")
    log("    * A zero would have been the sharp circulant object H = G(r-t)"
        " (THM-2365 (18),")
    log("      THM-2367 Sec.4), i.e. the precise functional form of the"
        " inverse-character")
    log("      hostile -- reported as a law, not a failure.")
    log("")
    log(f"{elapsed()} all checks passed")


if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    main()
