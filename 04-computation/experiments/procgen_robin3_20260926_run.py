#!/usr/bin/env python3
"""procgen_robin3_20260926_run -- runner for the robin3 lane (session collatz-procgen-20260922).

Theorem R_2: N_m(L) <= 2.1287 A_(m+2)(L) and Theorem R_3: N_m(L) <= 1.1632 A_(m+3)(L), for ALL m >= 2 and L >= 1
(HYP-9142 with shifts 2 and 3; THM-4513 has shift 15; the conjectured shift 1 stays open).
Proof architecture (note procgen_robin3_20260926_shift_one.md): robin2 Proposition 1 + Lemma 2 + Lemma 3, with
  * EM (eventual monotonicity) from a NEW criterion (FKG on the top-avoiding sublattice): (n-e)/(e+1) <= 1 - P_e(TH),
    exact for 2 <= m < m_E (all landing Sturmian factors), analytic for m >= m_E (Lemma LR + Ville + Hoeffding);
  * TB (top-hit bound on the last k1 steps) from a NEW decoupled bound B' (exact survival sequences) for m <= m_X, and
    robin2 regime A + a sharpened bridge regime D (Lemma TOP instead of Hoeffding union bounds) for m > m_X.
Shift 2 adds the refined EM criterion EM** (hitting lower bounds, Lemmas LR-, TOP-LOW) and the room-weighted regime D'.
Diagnostics for shifts 0, 1, 2.
Every claim goes through check(); a failure aborts. Output: stdout only.
"""
import io
import math
import os
import random
import resource
import sys
import time
from contextlib import redirect_stdout
from fractions import Fraction as F

import numpy as np
from mpmath import mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import procgen_robin2_20260926_lib as R2        # noqa: E402  read-only
import procgen_robin3_20260926_lib as L         # noqa: E402

T_START = time.time()
NCHECK = 0


def check(cond, msg):
    global NCHECK
    if not cond:
        print("[FAIL] " + msg, flush=True)
        raise SystemExit("CHECK FAILED: " + msg)
    NCHECK += 1
    print("[OK] " + msg, flush=True)


def section(title):
    print("\n=== %s   (t=%.1fs)" % (title, time.time() - T_START), flush=True)


def up(x, d=5):
    """x rounded UP to d decimals (for printing upper bounds)"""
    return math.ceil(x * 10 ** d - 1e-12) / 10 ** d


# ------------------------------------------------------------------------------------------------ parameters
C0 = 3
V0 = F(11, 100)          # EM bridge length n1(m) = ceil((m+1)/V0); TB window k1(m) = n1(m) + 1
M_E = 80                 # EM: exact for 2 <= m < M_E, analytic for m >= M_E
M_X = 500                # TB: B' exact for 2 <= m <= M_X, regimes A + D for m > M_X
V_MIN = F(105, 1000)     # regime D class boundaries
V_SMALL = F(8, 100)
J_ROOM = 19              # room-weighted regime D': number of endpoint levels j = 0..J used
# shift 2 (Theorem R_2)
V0_2 = F(8, 100)          # window parameter for m >= M_E2 (EM** needs it)
V0_2S = F(1, 10)          # window parameter for 25 <= m < M_E2 (exact EM; shorter windows, better B')
M_E2 = 150
M_X2 = 1000
V_MIN2 = F(75, 1000)
V_SMALL2 = F(5, 100)
J_ROOM2 = 40
J_EM2 = 12               # EM**: number of leading-zero levels j = 1..J with a hitting lower bound
T0_EM2 = 600             # EM**: time horizon of the hitting lower bound
THETAS = np.linspace(0.20, 1.095, 180)
SMALL_R = {1: F(1156186, 10 ** 6), 2: F(1047346, 10 ** 6), 3: F(1014780, 10 ** 6)}   # robin2 Theorem R_small, m <= 24
TMAX = 400000
DL = R2.deltas(TMAX)
DB = bytes(DL)
random.seed(20260926)
mp.dps = 50

section("0. parameters and constants")
print("c0=%d v0=%s m_E=%d m_X=%d v_min=%s v_small=%s ; c=%.12f d=%.12f mu=%.12f" % (C0, V0, M_E, M_X, V_MIN, V_SMALL,
                                                                                  L.C, L.D, L.MU))
check(abs(L.C - L.lo(L.I_C)) < 1e-15 and (L.I_C.b - L.I_C.a) < 1e-30,
      "interval constant c = log_3 2 has width < 1e-30 and agrees with the float value to 1e-15")
rob = mp.sqrt(2 / mp.pi) * mp.exp(-mp.mpf(1) / 6)
check(mp.mpf('0.6753') < rob, "Robbins constant sqrt(2/pi) e^(-1/6) = %s > 0.6753" % mp.nstr(rob, 8))

# ------------------------------------------------------------------------------------------------ 1
section("1. model: robin2 counters vs the orchestrator's independent counters (read-only)")
src = open(os.path.join(HERE, "procgen_robin2_20260926_orchestrator_check.py")).read().split("# 1. N_1 = 0")[0]
ORCH = {}
with redirect_stdout(io.StringIO()):
    exec(compile(src, "orchestrator_counters", "exec"), ORCH)
okm = True
for m in (2, 4, 7):
    Nr = R2.robin_counts(300, m, DL)
    No = ORCH["N_barrier"](m, 300)
    Ar = R2.dirichlet_counts(300, m + 3, DL)
    Ao = ORCH["A_hard"](m + 3, 300)
    okm &= (Nr == No) and (Ar == Ao)
check(okm, "N_m(L), A_(m+3)(L) of the robin2 library equal the orchestrator's counters (m = 2, 4, 7; L <= 300)")

# ------------------------------------------------------------------------------------------------ 2
section("2. lemma sanity checks (exact or exact-float; not used by the proof)")
# 2a: Lemma A (sublattice FKG) and the EM criterion chain, exhaustive small bridges
import itertools   # noqa: E402
cases = bad = 0
for trial in range(400):
    n = random.randint(8, 15)
    s = random.randint(2, 20000)
    while not (DL[s - 2] == 0 and DL[s - 1] == 1):
        s += 1
    word = DL[s:s + n]
    m = random.randint(2, 6)
    M = m + random.randint(0, 3)
    e = sum(word) + 1 - m
    if e < 1 or e + 1 > n:
        continue

    def paths(x0, ee):
        out = []
        for A in itertools.combinations(range(n), ee):
            bits = [0] * n
            for a in A:
                bits[a] = 1
            out.append((A, L.vpath(x0, bits, word)))
        return out
    Pm = paths(m, e)
    Tm = [(A, V) for A, V in Pm if all(v <= M for v in V[1:n])]
    Cm = [(A, V) for A, V in Tm if all(v >= 1 for v in V[1:])]
    if not Cm:
        continue
    cases += 1
    Ef = F(sum(A[0] for A, V in Pm), len(Pm))
    ET = F(sum(A[0] for A, V in Tm), len(Tm))
    EC = F(sum(A[0] for A, V in Cm), len(Cm))
    PT = F(len(Tm), len(Pm))
    bad += not (Ef == F(n - e, e + 1) and EC <= ET <= Ef / PT)
    for _ in range(10):
        (A1, V1), (A2, V2) = random.choice(Tm), random.choice(Tm)
        bj, bm = [0] * n, [0] * n
        for a in (min(x, y) for x, y in zip(A1, A2)):
            bj[a] = 1
        for a in (max(x, y) for x, y in zip(A1, A2)):
            bm[a] = 1
        bad += (L.vpath(m, bj, word) != [max(x, y) for x, y in zip(V1, V2)]
                or L.vpath(m, bm, word) != [min(x, y) for x, y in zip(V1, V2)])
    Cm1 = sum(1 for A, V in paths(m - 1, e + 1) if all(1 <= v <= M for v in V[1:n]) and V[n] >= 1)
    bad += Cm1 > sum(A[0] for A, V in Cm)
check(cases > 300 and bad == 0,
      "2a Lemma A: on %d exhaustive small bridges: join/meet = pointwise max/min of paths (sublattice), "
      "E[f] = (n-e)/(e+1), E_C[f] <= E_T[f] <= E[f]/P(T), and |C_(m-1)| <= sum_(C_m) f (injection)" % cases)
# 2b: Lemma LR on random hitting configurations (exact rationals, 50-digit comparison)
nlr = 0
worst = 0.0
for trial in range(4000):
    n = random.randint(64, 700)
    e = random.randint(int(0.4 * n) + 1, int(0.6 * n) - 1)
    p = F(e, n)
    y = random.uniform(2, 5)
    t = random.randint(1, n // 2)
    i = math.floor(y + L.C * t) + 1
    delta = i - p * t
    if not (delta >= 2 and delta / (n - t) <= F(1, 4) and e - i >= 1 and (n - e) - (t - i) >= 1):
        continue
    LRv = L.exact_LR(n, e, t, i)
    ratio = (mp.mpf(LRv.numerator) / LRv.denominator) / mp.exp(mp.mpf(t) / n + mp.mpf(1) / (6 * n))
    worst = max(worst, float(ratio))
    nlr += 1
check(nlr > 2000 and worst <= 1.0,
      "2b Lemma LR (delta >= 2): exact LR_t(i) <= exp(t/n + 1/(6n)) on %d hitting configurations (max ratio %.4f)"
      % (nlr, worst))
# 2b': Lemma LR- (lower bound) on random configurations, exact rationals
nlr = 0
okl = True
for trial in range(3000):
    n = random.randint(200, 1500)
    e = random.randint(int(0.45 * n) + 1, int(0.6 * n) - 1)
    p = F(e, n)
    t = random.randint(1, n // 3)
    y = random.uniform(1, 8)
    i = math.floor(y + L.C * t) + 1
    if not (1 <= i <= e - 1 and 0 <= t - i <= (n - e) - 1):
        continue
    delta = float(i - p * t)
    if delta <= 0:
        continue
    pf = float(p)
    lb = L.lr_lower_log(n, t, delta, pf, pf)
    if lb == -math.inf:
        continue
    LRv = L.exact_LR(n, e, t, i)
    okl &= mp.log(mp.mpf(LRv.numerator) / LRv.denominator) >= lb
    nlr += 1
check(okl and nlr > 1000, "2b' Lemma LR-: exact log LR_t(i) >= the lower bound on %d random configurations" % nlr)


def bridge_top_hit(word, x0, x1, top):
    """float: P(V > top at some time 1..n-1) for the uniform bridge x0 -> x1 (top only; normalised DP).
    The array covers every site a path can reach (V in [x0 - n, max(top, x0) + n]), so the total count is exact."""
    n = len(word)
    off = n + 2
    size = max(top, x0) + n + 3 + off
    fa = np.zeros(size)
    fb = np.zeros(size)
    fa[x0 + off] = fb[x0 + off] = 1.0
    for i, d in enumerate(word):
        ga = np.zeros(size)
        gb = np.zeros(size)
        if d == 1:
            ga[:-1] += fa[1:]
            ga += fa
            gb[:-1] += fb[1:]
            gb += fb
        else:
            ga += fa
            ga[1:] += fa[:-1]
            gb += fb
            gb[1:] += fb[:-1]
        if i < n - 1:
            gb[top + 1 + off:] = 0.0
        z = ga.sum()
        fa, fb = ga / z, gb / z
    return 1.0 - fb[x1 + off] / fa[x1 + off]


def top_bound(n, e, y, v):
    eta = L.eta_plus_float(e / n, 1.0 / n)
    return math.exp(1 / (6 * n)) * math.exp(-eta * y) + math.exp(-v * v * n - 4 * v * y) / (1 - math.exp(-4 * v * v))


# 2c: Lemma TOP vs exact bridge top-hit probabilities (EM bridges and TB bridges)
worst = 0.0
cnt = 0
for m in (80, 120):
    M = m + C0
    n = L.n1_of(m, V0)
    words, comp = L.landing_words(DB, n, 250000)
    for w, s in random.sample(list(words.items()), 20):
        e = sum(w) + 1 - m
        a = m - (L.C * s) % 1.0
        b = 1 - (L.C * (s + n)) % 1.0
        worst = max(worst, bridge_top_hit(list(w), m, 1, M) / top_bound(n, e, M - a, (a - b) / n))
        cnt += 1
for m in (80, 150):
    M = m + C0
    k1 = L.n1_of(m, V0) + 1
    kb = L.k_b(m)
    for _ in range(8):
        u = random.randint(1, 100000)
        k = random.randint(kb + 1, k1)
        word = DL[u:u + k]
        a = m - (L.C * u) % 1.0
        for e in range(0, k + 1, max(1, k // 50)):
            z = a + e - L.C * k
            v = (a - z) / k
            if 0.05 <= v <= 0.2 and z > 0:
                worst = max(worst, bridge_top_hit(word, m, m + e - sum(word), M) / top_bound(k, e, M - a, v))
                cnt += 1
check(worst < 1.0, "2c Lemma TOP: exact bridge top-hit probability <= bound on %d EM/TB bridges (max ratio %.3f)"
      % (cnt, worst))
# 2c': Lemma TOP-LOW (hitting lower bound used in EM**) against exact probabilities for bridges after j leading zeros
okt = True
cnt = 0
worst = 9.0
for m in (150,):
    M = m + 2
    n = L.n1_of(m, V0_2)
    words, comp = L.landing_words(DB, n, 250000)
    for w, s0 in random.sample(list(words.items()), 6):
        e = sum(w) + 1 - m
        for j in (1, 3, 6):
            rest = list(w[j:])
            x0 = m - sum(w[:j])            # site after j leading zeros
            ex = bridge_top_hit(rest, x0, 1, M)
            cdf = L.iid_hit_cdf_lower(L.C - float(V0_2), 2, j, 400)
            nprime = n - j
            p_act = e / nprime
            Ls = []
            for t in range(1, 401):
                ll = L.lr_lower_log(nprime, t, 2 + L.C * (j + 1) + float(V0_2) * t + 1, p_act, p_act)
                Ls.append(math.exp(ll) if ll > -math.inf else 0.0)
            for i2 in range(1, len(Ls)):
                Ls[i2] = min(Ls[i2], Ls[i2 - 1])
            lowb = sum((Ls[t - 1] - (Ls[t] if t < 400 else 0.0)) * cdf[t] for t in range(1, 401))
            okt &= lowb <= ex
            worst = min(worst, ex / lowb if lowb > 0 else 9.0)
            cnt += 1
check(okt, "2c' Lemma TOP-LOW: exact top-hit probability of the bridge after j leading zeros >= lower bound on %d cases "
      "(min exact/bound %.3f)" % (cnt, worst))
# 2d: Lemma B' vs the exact top-hit ratio q = 1 - V_k(m)/H_k(m) at sampled start times
okb = True
rows = []
for m in (25, 60, 150):
    M = m + C0
    k1 = L.n1_of(m, V0) + 1
    qB = L.tb_Bprime(m, C0, k1, DL, THETAS)[0]
    qmax = 0.0
    for _ in range(12):
        u = random.randint(1, 100000)
        cap = m + 200
        fH = np.zeros(cap + 3)
        fV = np.zeros(M + 3)
        fH[m] = fV[m] = 1.0
        for k in range(1, k1 + 1):
            d = DL[u + k - 1]
            gH = np.zeros(cap + 3)
            gV = np.zeros(M + 3)
            if d == 1:
                gH[:-1] += fH[1:]
                gH += fH
                gV[:-1] += fV[1:]
                gV += fV
            else:
                gH += fH
                gH[1:] += fH[:-1]
                gV += fV
                gV[1:] += fV[:-1]
            gH *= 0.5
            gV *= 0.5
            gH[0] = gV[0] = 0.0
            qmax = max(qmax, 1.0 - gV[1:].sum() / gH.sum())
            gV[M + 1:] = 0.0
            fH, fV = gH, gV
    okb &= qmax <= qB
    rows.append("m=%d: sampled q %.4f <= q_B' %.4f" % (m, qmax, qB))
check(okb, "2d Lemma B' dominates the exact top-hit ratio at sampled phases: " + "; ".join(rows))
# 2e: Lemma DIP (Hoeffding union) vs exact dip probabilities of small bridges
okd = True
nd = 0
for trial in range(200):
    k = random.randint(40, 200)
    u = random.randint(1, 50000)
    word = DL[u:u + k]
    h = random.randint(8, 30)
    a = h - (L.C * u) % 1.0
    e = random.randint(max(0, math.ceil(L.C * k - a + 0.5)), max(0, math.floor(L.C * k)))
    z = a + e - L.C * k
    if z <= 0.5 or (a - z) <= 0:
        continue
    v = (a - z) / k
    x1 = h + e - sum(word)
    # exact dip probability: bridge from h to x1 with a bottom kill (V >= 1 at times 1..k-1)
    f = {h: 1}
    g_all = {h: 1}
    for i2, d in enumerate(word):
        nf, ng = {}, {}
        for x, c in f.items():
            for x2 in (x - d, x + 1 - d):
                if x2 >= 1 or i2 == k - 1:
                    nf[x2] = nf.get(x2, 0) + c
        for x, c in g_all.items():
            for x2 in (x - d, x + 1 - d):
                ng[x2] = ng.get(x2, 0) + c
        f, g_all = nf, ng
    pdip = 1 - F(f.get(x1, 0), g_all[x1])
    bound = sum(math.exp(-2 * (z + v * s2) ** 2 / s2) for s2 in range(1, k))
    if float(pdip) > bound + 1e-12:
        print("  DIP counterexample? k=%d h=%d e=%d u=%d a=%.4f z=%.4f v=%.5f pdip=%.6g bound=%.6g"
              % (k, h, e, u, a, z, v, float(pdip), bound))
    okd &= float(pdip) <= bound + 1e-12
    nd += 1
check(okd and nd > 50, "2e Lemma DIP: exact bridge dip probability <= sum_s exp(-2(z+vs)^2/s) on %d bridges" % nd)

# ------------------------------------------------------------------------------------------------ 3
section("3. EM, exact part: 2 <= m < m_E, all landing Sturmian factors of length n1(m) = ceil((m+1)/v0)")
worst = (9.0, None)
okx = True
tab3 = []
for m in range(2, M_E):
    M = m + C0
    n = L.n1_of(m, V0)
    words, complete = L.landing_words(DB, n, 150000)
    Pm, Pm1 = L.em_kernel_check(list(words), m, M)
    err = 2.1 * n * L.U_ROUND
    ratio = float((Pm / Pm1).min())
    okx &= complete and bool((Pm >= Pm1 * (1 + err)).all()) and bool((Pm1 > 0.5 * 2.0 ** -n).all())
    if ratio < worst[0]:
        worst = (ratio, m)
    if m in (2, 3, 5, 10, 20, 40, 60, 79):
        tab3.append("m=%d:n1=%d,%d factors,min %.4f" % (m, n, len(words), ratio))
print("  " + "; ".join(tab3))
check(okx, "3: for 2 <= m < %d and every landing factor, P_(s,s+n1)(m,1) >= P_(s,s+n1)(m-1,1) > 0 (float-rigorous, "
      "relative error <= 1.01 n u; factor sets complete); smallest ratio %.4f at m = %d" % (M_E, worst[0], worst[1]))

# ------------------------------------------------------------------------------------------------ 4
section("4. EM, analytic part: m >= m_E (Lemma A + Lemma TOP), interval arithmetic at m_E + monotone terms")
r = L.em_analytic(M_E, C0, V0)
print("  n1(m_E)=%d  Ef<=%.6f  Top<=%.6f (eta_lo=%.6f, a lower bound)  Late<=%.4e  total<=%.6f  (upper bounds rounded up)" % (
    r['n'], up(L.hi(r['Ef']), 6), up(L.hi(r['top']), 6), r['eta'], up(L.hi(r['late']), 7), up(L.hi(r['total']), 6)))
check(r['ok'], "4: for every m >= %d and landing time s: (n-e)/(e+1) + P_e(TH) <= %.6f < 1 (Lemma LR side conditions "
      "hold), hence P(m,1) >= P(m-1,1)" % (M_E, up(L.hi(r['total']), 6)))
cons = [L.hi(L.em_analytic(mm, C0, V0)['total']) for mm in (M_E, 2 * M_E, 10 * M_E, 100 * M_E)]
check(all(x >= y - 1e-15 for x, y in zip(cons, cons[1:])),
      "4b (consistency) the certified EM total is nonincreasing along m = m_E, 2m_E, 10m_E, 100m_E: "
      + ", ".join("%.6f" % x for x in cons))

# ------------------------------------------------------------------------------------------------ 5
section("5. top-hit bound (TB) on the last k1(m) = n1(m)+1 steps")
qB = {}
t5 = time.time()
for m in range(2, M_X + 1):
    k1 = L.n1_of(m, V0) + 1
    qB[m] = L.tb_Bprime(m, C0, k1, DL, THETAS)[0]
print("  B' (%.1fs, rounded up): " % (time.time() - t5) + ", ".join("m=%d:%.4f" % (m, up(qB[m], 4)) for m in
                                                       (2, 3, 4, 5, 8, 12, 24, 25, 50, 100, 200, 300, 400, 500)))
q25 = max(qB[m] for m in range(25, M_X + 1))
q2 = max(qB.values())
check(all(v < 1 for v in qB.values()),
      "5a regime B' (exact survival, float-rigorous): q <= %.5f for 25 <= m <= %d, q <= %.5f for 2 <= m <= %d"
      % (up(q25), M_X, up(q2), M_X))
dw = L.tb_regime_Dw(M_X, C0, V0, V_MIN, V_SMALL, J_ROOM)
d = dw['base']
qA = L.tb_regime_A(C0)
print("  regime D' at m_X: k_b=%d k1=%d  Top_L*c/v_min<=%.5f  Phi<=%.4f (r_lo=%.4f, J=%d)  Q_L'<=%.5f  Q_H<=%.5f "
      "(dip<=%.3e, z*=%.2f)  R_B<=%.3e  Delta=%.5f  q_D'<=%.5f ; regime A: q <= 4/(3*3^c0) = %.5f" % (
          d['kb'], d['k1'], L.hi(d['QL']), L.hi(dw['Phi']), L.lo(dw['rlo']), J_ROOM, L.hi(dw['QLw']), L.hi(d['QH']),
          L.hi(d['dip']), L.lo(d['zstar']), L.hi(d['RB']), L.lo(d['Delta']), L.hi(dw['qDw']), float(qA)))
check(dw['ok'], "5b regimes A + D' (room-weighted): for every m > %d, start time u >= 1 and 1 <= k <= k1(m): "
      "q <= max(%.5f, %.5f)" % (M_X, up(float(qA)), up(L.hi(dw['qDw']))))
cons = [L.hi(L.tb_regime_Dw(mm, C0, V0, V_MIN, V_SMALL, J_ROOM)['qDw']) for mm in (M_X, 2 * M_X, 10 * M_X)]
check(all(x >= y - 1e-12 for x, y in zip(cons, cons[1:])),
      "5c (consistency) q_D' nonincreasing along m = m_X, 2m_X, 10m_X: " + ", ".join("%.5f" % x for x in cons))

# ------------------------------------------------------------------------------------------------ 6
section("6. assembly: Theorem R_3")
qD = max(float(qA), L.hi(dw['qDw']))
G_small = SMALL_R[3]
G_mid = 1.0 / (1.0 - q25)
G_big = 1.0 / (1.0 - qD)
G3 = max(float(G_small), G_mid, G_big)
G3_self = max(1.0 / (1.0 - q2), G_big)
print("  Gamma_3(m) (rounded up): m <= 24: %.6f (robin2 Theorem R_small); 25 <= m <= %d: %.6f (B'); m > %d: %.6f (A + D')" % (
    float(G_small), M_X, up(G_mid, 6), M_X, up(G_big, 6)))
print("  Gamma_3 <= %.6f ; self-contained (B' also for m <= 24): %.6f" % (up(G3, 6), up(G3_self, 6)))
GAMMA3 = F(11632, 10000)
check(G3 <= float(GAMMA3) and G3_self <= 2.8251,
      "6: Theorem R_3: N_m(L) <= %s A_(m+3)(L) for all m >= 2, L >= 1 (constant <= %.6f; <= %.4f without R_small)"
      % (float(GAMMA3), up(G3, 6), up(G3_self, 4)))

# ------------------------------------------------------------------------------------------------ 7
section("7. exact verification of Theorem R_3 on a finite range and the observed ratio")
okv = True
worst = F(0)
arg = None
for m in range(2, 41):
    Lm = 800 if m <= 20 else 500
    N = R2.robin_counts(Lm, m, DL)
    A = R2.dirichlet_counts(Lm, m + 3, DL)
    for Lx in range(1, Lm + 1):
        okv &= N[Lx] * GAMMA3.denominator <= GAMMA3.numerator * A[Lx]
        if A[Lx] and N[Lx] * worst.denominator > worst.numerator * A[Lx]:
            worst = F(N[Lx], A[Lx])
            arg = (m, Lx)
check(okv, "7: N_m(L) <= %s A_(m+3)(L) exactly for 2 <= m <= 40, L <= 800/500; observed max ratio 1 + %.3e at "
      "(m,L) = %s" % (float(GAMMA3), float(worst - 1), arg))

# ------------------------------------------------------------------------------------------------ 8
section("8. Theorem R_2: shift 2 for all m (EM exact + EM**; B' + regimes A, D')")
# 8a EM exact at shift 2, 25 <= m < M_E2, bridge length n = ceil((m+1)/V0_2)
okx = True
worst = (9.0, None)
t8 = time.time()
for m in range(25, M_E2):
    n = L.n1_of(m, V0_2S)
    words, complete = L.landing_words(DB, n, 300000)
    Pm, Pm1 = L.em_kernel_check(list(words), m, m + 2)
    err = 2.1 * n * L.U_ROUND
    okx &= complete and bool((Pm >= Pm1 * (1 + err)).all()) and bool((Pm1 > 0.5 * 2.0 ** -n).all())
    ratio = float((Pm / Pm1).min())
    if ratio < worst[0]:
        worst = (ratio, m)
check(okx, "8a EM at shift 2, exact part: for 25 <= m < %d and every landing factor of length n1(m) = 10(m+1), "
      "P(m,1) >= P(m-1,1) > 0 (float-rigorous); smallest ratio %.4f at m = %d (%.1fs)"
      % (M_E2, worst[0], worst[1], time.time() - t8))
# 8b EM** (refined criterion) at m >= M_E2
r2 = L.em_refined(M_E2, 2, V0_2, J_EM2, T0_EM2)
print("  EM** at m_E=%d: n1=%d h0<=%.5f  h_1..h_4 >= %s  value <= %.5f" % (
    M_E2, r2['n'], up(r2['h0']), ", ".join("%.4f" % (math.floor(h * 1e4) / 1e4) for h in r2['hs'][:4]), up(r2['value'])))
check(r2['ok'], "8b EM at shift 2, analytic part (Proposition EM**): E_T[f] <= %.5f < 1 for every m >= %d and landing "
      "time s" % (up(r2['value']), M_E2))
cons = [L.em_refined(mm, 2, V0_2, J_EM2, T0_EM2)['value'] for mm in (M_E2, 2 * M_E2, 10 * M_E2)]
check(all(x >= y - 1e-12 for x, y in zip(cons, cons[1:])),
      "8c (consistency) the EM** value is nonincreasing along m = m_E, 2m_E, 10m_E: " + ", ".join("%.5f" % x for x in cons))
# 8d B' at shift 2, 25 <= m <= M_X2
t8 = time.time()
qB2 = {m: L.tb_Bprime(m, 2, L.n1_of(m, V0_2S if m < M_E2 else V0_2) + 1, DL, THETAS)[0] for m in range(25, M_X2 + 1)}
q2a = max(qB2[m] for m in range(25, M_E2))
q2b = max(qB2[m] for m in range(M_E2, M_X2 + 1))
q2max = max(q2a, q2b)
print("  B' at shift 2 (%.1fs, rounded up): " % (time.time() - t8) + ", ".join("m=%d:%.4f" % (m, up(qB2[m], 4)) for m in
                                                                (25, 50, 100, 149, 150, 300, 500, 1000)))
check(q2max < 1, "8d TB at shift 2, Lemma B': q <= %.5f for 25 <= m < %d (k1 = 10(m+1)+1) and q <= %.5f for "
      "%d <= m <= %d (k1 = ceil((m+1)/%.2f)+1)" % (up(q2a), M_E2, up(q2b), M_E2, M_X2, float(V0_2)))
# 8e regimes A + D' at shift 2, m > M_X2
dw2 = L.tb_regime_Dw(M_X2, 2, V0_2, V_MIN2, V_SMALL2, J_ROOM2)
qA2 = L.tb_regime_A(2)
print("  regime D' at shift 2, m_X=%d: Phi<=%.4f (r_lo=%.4f)  Q_L'<=%.5f  Q_H<=%.5f  R_B<=%.3e  q_D'<=%.5f ; "
      "regime A: %.5f" % (M_X2, L.hi(dw2['Phi']), L.lo(dw2['rlo']), L.hi(dw2['QLw']), L.hi(dw2['base']['QH']),
                          L.hi(dw2['base']['RB']), L.hi(dw2['qDw']), float(qA2)))
check(dw2['ok'], "8e TB at shift 2, regimes A + D': for every m > %d, u >= 1, k <= k1(m): q <= max(%.5f, %.5f)"
      % (M_X2, up(float(qA2)), up(L.hi(dw2['qDw']))))
cons = [L.hi(L.tb_regime_Dw(mm, 2, V0_2, V_MIN2, V_SMALL2, J_ROOM2)['qDw']) for mm in (M_X2, 2 * M_X2, 10 * M_X2)]
check(all(x >= y - 1e-12 for x, y in zip(cons, cons[1:])),
      "8f (consistency) q_D' at shift 2 nonincreasing along m = m_X, 2m_X, 10m_X: " + ", ".join("%.5f" % x for x in cons))
# 8g assembly
qbig2 = max(float(qA2), L.hi(dw2['qDw']))
G2_small = float(SMALL_R[2])
G2_mid = 1.0 / (1.0 - q2max)       # max over 25 <= m <= M_X2 (both window choices)
G2_big = 1.0 / (1.0 - qbig2)
G2 = max(G2_small, G2_mid, G2_big)
GAMMA2 = F(21285, 10000)
print("  Gamma_2(m) (rounded up): m <= 24: %.6f (R_small); 25 <= m < %d: %.6f, %d <= m <= %d: %.6f (B'); m > %d: %.6f (A + D')" % (
    G2_small, M_E2, up(1 / (1 - q2a), 6), M_E2, M_X2, up(1 / (1 - q2b), 6), M_X2, up(G2_big, 6)))
check(G2 <= float(GAMMA2), "8g Theorem R_2: N_m(L) <= %s A_(m+2)(L) for all m >= 2, L >= 1 (constant <= %.6f)"
      % (float(GAMMA2), up(G2, 6)))
# 8h exact verification
okv = True
worst = F(0)
arg = None
for m in range(2, 41):
    Lm = 800 if m <= 20 else 500
    N = R2.robin_counts(Lm, m, DL)
    A = R2.dirichlet_counts(Lm, m + 2, DL)
    for Lx in range(1, Lm + 1):
        okv &= N[Lx] * GAMMA2.denominator <= GAMMA2.numerator * A[Lx]
        if A[Lx] and N[Lx] * worst.denominator > worst.numerator * A[Lx]:
            worst = F(N[Lx], A[Lx])
            arg = (m, Lx)
check(okv, "8h: N_m(L) <= %s A_(m+2)(L) exactly for 2 <= m <= 40, L <= 800/500; observed max ratio 1 + %.3e at (m,L) = %s"
      % (float(GAMMA2), float(worst - 1), arg))

# ------------------------------------------------------------------------------------------------ 9
section("9. diagnostics: why the method stops at shift 2 (EMPIRICAL / method statements, float)")
# 9a shift 1: EM bridge-crossing onset at sampled landing phases, incl. the near-worst phase {cs} -> 2c-1
land = [s for s in range(2, 150000) if DL[s - 2] == 0 and DL[s - 1] == 1]
near = sorted(land, key=lambda s: (L.C * s) % 1.0)[:3]
rows = []
onsets = []
for m in (10, 20, 40, 80):
    M = m + 1
    desc = (m - 1) / L.MU
    best_on, ev = 0, 9.0
    for s in near + random.sample(land, 5):
        fa = np.zeros(M + 3)
        fb = np.zeros(M + 3)
        fa[m - 1] = fb[m] = 1.0
        cross = None
        nmax = int(12 * desc) + 300
        for i in range(nmax):
            dd = DL[s + i]
            for X in (fa, fb):
                g = np.zeros(M + 3)
                if dd == 1:
                    g[:-1] += X[1:]
                    g += X
                else:
                    g += X
                    g[1:] += X[:-1]
                g[0] = 0.0
                g[M + 1:] = 0.0
                X[:] = g
            z = fa.sum() + fb.sum()
            fa /= z
            fb /= z
            if cross is None and fa[1] > 0 and fb[1] >= fa[1]:
                cross = i + 1
        dn = DL[s + nmax]
        ratio = (2 * fb[1:M + 1].sum() - (fb[1] if dn else 0)) / (2 * fa[1:M + 1].sum() - (fa[1] if dn else 0))
        best_on = max(best_on, cross / desc)
        ev = min(ev, ratio)
    onsets.append(best_on)
    rows.append("m=%d: onset %.2f desc, eventual ratio %.4f" % (m, best_on, ev))
print("  shift 1: " + "; ".join(rows))
check(all(x < y for x, y in zip(onsets, onsets[1:])),
      "9a (EMPIRICAL) shift 1: the EM bridge-crossing onset / descent time grows with m (%s); the worst-phase "
      "eventual ratio V_k(m)/V_k(m-1) is only ~1.01-1.06" % ", ".join("%.2f" % x for x in onsets))
# 9b shift 2: the analytic EM criterion (d+v0)/(c-v0) + exp(-eta_+(c-v0) (2c+1)) exceeds 1 for every v0 in (0, mu)
vals = []
for v in np.linspace(0.005, L.MU - 0.002, 250):
    eta = L.eta_plus_float(L.C - v, 0.0)
    vals.append((L.D + v) / (L.C - v) + math.exp(-eta * (2 + 2 * L.C - 1)))
check(min(vals) > 1.0, "9b (method) shift 2: the simple criterion E[f] + e^(-eta y) of Proposition EM* is >= %.4f > 1 "
      "for every v0 in (0, mu) (grid of 250); shift 2 needs the refined criterion EM** (section 8b)" % min(vals))
# 9c shift 0: N_m/A_m grows (exact)
N3 = R2.robin_counts(2000, 3, DL)
A3 = R2.dirichlet_counts(2000, 3, DL)
r1000 = F(N3[1000], A3[1000])
r2000 = F(N3[2000], A3[2000])
check(r1000 > 300 and r2000 > 10 ** 5, "9c (FINITE-EXACT) shift 0: N_3(L)/A_3(L) = %.4g at L=1000 and %.4g at L=2000 "
      "(geometric growth; shift 0 is false, EMPIRICAL)" % (float(r1000), float(r2000)))

# 9d shift 1: the decoupled bound B' exceeds 1 on EM-sized windows (window ~3.6 descent times)
rows = []
q1 = []
for m in (10, 20, 40, 80):
    q = L.tb_Bprime(m, 1, L.n1_of(m, F(36, 1000)) + 1, DL, THETAS)[0]
    q1.append(q)
    rows.append("m=%d:%.3f" % (m, q))
check(all(x > 1 for x in q1), "9d (method) shift 1: with the window k1 = ceil((m+1)/0.036)+1 (~3.6 descent times, needed "
      "for EM at shift 1) Lemma B' gives q_B' > 1: " + ", ".join(rows))
# 9e shift 2 (heuristic, iid values; not a proof): the refined criterion sum_j p0^j (1 - e^{-eta(y+cj+d)}) / (1 - e^{-eta y})
rows = []
for v in (0.04, 0.06, 0.08, 0.09, 0.10):
    eta = L.eta_plus_float(L.C - v, 0.0)
    y = 2 + 2 * L.C - 1
    p0 = L.D + v
    num = sum(p0 ** j * (1 - math.exp(-eta * (y + L.C * j + L.D))) for j in range(1, 400))
    rows.append("v=%.2f:%.4f" % (v, num / (1 - math.exp(-eta * y))))
print("  9e (heuristic) shift 2 refined EM criterion, iid values: " + ", ".join(rows))

# ------------------------------------------------------------------------------------------------ summary
section("summary")
ru = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
print("Theorem R_3: N_m(L) <= %s A_(m+3)(L) for all m >= 2, L >= 1." % float(GAMMA3))
print("Theorem R_2: N_m(L) <= %s A_(m+2)(L) for all m >= 2, L >= 1." % float(GAMMA2))
print("peak RSS (ru_maxrss, bytes on macOS) = %d  (%.0f MB)" % (ru, ru / 2 ** 20))
print("checks passed: %d   wall time %.1fs" % (NCHECK, time.time() - T_START))
print("ALL CHECKS PASSED")
