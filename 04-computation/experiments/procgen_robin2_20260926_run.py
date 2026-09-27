#!/usr/bin/env python3
"""procgen_robin2_20260926_run -- re-checks every finite claim of
05-knowledge/results/procgen_robin2_20260926_constant_robin.md  (lane robin2, session collatz-procgen-20260922).

Main theorem checked here:  N_m(L) <= GAMMA * A_{m+K}(L) for all m >= 2, L >= 1, with K = 15 and the explicit GAMMA
assembled in section 7.  Every claim goes through check(); a failure aborts.  Output: stdout only.
Usage:  python3 -u procgen_robin2_20260926_run.py > 05-knowledge/results/procgen_robin2_20260926.out
"""
import sys, os, math, time, itertools, random, resource
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import procgen_robin2_20260926_lib as R

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
    print("\n== %s ==  (t=%.1fs)" % (title, time.time() - T_START), flush=True)


C, D, MU = R.C, R.D, R.MU
# ------------------------------------------------------------------------------------------------ parameters
K = 15            # shift c0
W = 1e-3          # EM margin: n1(m) = ceil((m+1)/(mu - W))
M_E = 26          # eventual monotonicity: exact for 2 <= m < M_E, analytic for m >= M_E
M_X = 60          # top-hit window: exact worst-phase survival for m <= M_X
M_2 = 2550        # top-hit window: crude closed form for M_X < m < M_2, bridge bound for m >= M_2
V_MIN = 0.11      # bridge bound speed threshold
M_BIG = 10 ** 6   # bridge bound checked exactly for M_2 <= m <= M_BIG, monotone majorant beyond


def n1(m):
    """ceil((m+1)/(mu-W)), certified with a 1e-9 guard against rounding"""
    x = (m + 1) / (MU - W)
    n = math.ceil(x)
    if n - x < 1e-9 or x - math.floor(x) < 1e-9:
        raise SystemExit("n1 too close to an integer boundary at m=%d" % m)
    return n


section("0. constants")
theta_star = math.log(3.0)
check(abs((math.exp(theta_star * D) + math.exp(-theta_star * C)) / 2 - 1) < 1e-14,
      "theta* = ln 3 solves (e^{theta d} + e^{-theta c})/2 = 1  (3^d = 3/2, 3^-c = 1/2)")
print("c = %.12f  d = %.12f  mu = c - 1/2 = %.12f" % (C, D, MU))
print("parameters: K=%d W=%g M_E=%d M_X=%d M_2=%d V_MIN=%g M_BIG=%d" % (K, W, M_E, M_X, M_2, V_MIN, M_BIG))

# ------------------------------------------------------------------------------------------------ 1
section("1. V-coordinate model = the pairpeak/robin definitions of N_m(L) and A_M(L)")
import procgen_pairpeak_20260926_lib as PP   # read-only reuse (unchanged)
DL = R.deltas(80000)
FL = R.floors(80000)
LM = 160
flPP = PP.floor_table(400)
ok = True
for m in range(2, 15):
    ok &= (R.robin_counts(LM, m, DL) == PP.refl_counts_all_L(LM, m, flPP, 400))
for M in range(1, 30):
    ok &= (R.dirichlet_counts(LM, M, DL) == PP.hard_counts_all_L(LM, M, flPP, 400))
check(ok, "robin_counts/dirichlet_counts equal PP.refl_counts_all_L/hard_counts_all_L for m<=14, M<=29, L<=%d" % LM)
check(all(not (DL[t] == 0 and DL[t + 1] == 0) for t in range(80000)),
      "delta has no two consecutive zeros (c > 1/2), t < 80000")

# ------------------------------------------------------------------------------------------------ 2
section("2. structural lemmas: exact sanity checks on small cases")


def robin_distributions(Lmax, m, dl):
    """list over t of dicts site -> Robin weighted count R_t(site) from (0,0)"""
    f = {0: 1}
    out = [dict(f)]
    for t in range(Lmax):
        d = dl[t]
        g = {}
        for v, y in f.items():
            if v == m:
                g[m - 1] = g.get(m - 1, 0) + 2 * y
            else:
                for x in (v - d, v + 1 - d):
                    if x >= 1:
                        g[x] = g.get(x, 0) + y
        f = g
        out.append(dict(f))
    return out


# 2a Lemma 1 (hybrid telescoping identity)
ok = True
for m, Kk in [(4, 1), (5, 2), (7, 3)]:
    M = m + Kk
    Lmax = 70
    Rd = robin_distributions(Lmax, m, DL)
    N = R.robin_counts(Lmax, m, DL)
    A = R.dirichlet_counts(Lmax, M, DL)
    for L in range(1, Lmax + 1):
        tot = 0
        for t in range(L):
            rt = Rd[t].get(m, 0)
            if rt:
                j = L - t - 1
                Vlo = R.free_counts(t + 1, m - 1, j, M, DL)[j]
                Vhi = R.free_counts(t + 1, m, j, M, DL)[j]
                tot += rt * (Vlo - Vhi)
        ok &= (N[L] - A[L] == tot)
check(ok, "Lemma 1: N_m(L) - A_M(L) = sum_t R_t(m) [V^(t+1)_(L-t-1)(m-1) - V^(t+1)_(L-t-1)(m)] exactly (3 cases, L<=70)")

# 2b Lemma 2: W <= H and V/H nonincreasing in x
ok1 = ok2 = True
for m, Kk in [(6, 2), (9, 3), (12, 5)]:
    M = m + Kk
    for u in [1, 4, 9, 13, 30, 57]:
        kmax = 120
        H = {x: R.free_counts(u, x, kmax, None, DL) for x in range(1, M + 1)}
        V = {x: R.free_counts(u, x, kmax, M, DL) for x in range(1, M + 1)}
        for x in range(1, m + 1):
            if x == m and not (DL[u - 1] == 0 and DL[u] == 1):
                continue
            Wc = R.robin_from(u, x, kmax, m, DL)
            ok1 &= all(Wc[k] <= H[x][k] for k in range(kmax + 1))
        for k in range(1, kmax + 1):
            ok2 &= all(V[x + 1][k] * H[x][k] <= V[x][k] * H[x + 1][k] for x in range(1, M))
check(ok1, "Lemma 2a: W^(u)_k(x) <= H^(u)_k(x) (Robin <= half-line), sampled (m,K,u), k<=120")
check(ok2, "Lemma 2b: V^(u)_k(x)/H^(u)_k(x) nonincreasing in x on [1,M], sampled, k<=120")

# 2c Lemma 3: TP2 of the confined kernel and the propagation P(m,y) >= P(m-1,y)
def kernel_full(s, n, M, dl):
    """dict x -> list over y of P_{s,s+n}(x,y)"""
    res = {}
    for x in range(1, M + 1):
        f = [0] * (M + 3)
        f[x] = 1
        for i in range(n):
            d = dl[s + i]
            g = [0] * (M + 3)
            for v in range(1, M + 1):
                if f[v]:
                    g[v - d] += f[v]
                    g[v + 1 - d] += f[v]
            g[0] = 0; g[M + 1] = 0; g[M + 2] = 0
            f = g
        res[x] = f
    return res


ok = True
for M in [5, 8]:
    for s in [1, 3, 10]:
        for n in [1, 2, 5, 13, 30]:
            P = kernel_full(s, n, M, DL)
            for x in range(1, M):
                for y in range(1, M):
                    ok &= (P[x][y] * P[x + 1][y + 1] >= P[x][y + 1] * P[x + 1][y])
check(ok, "Lemma 3: TP2 of P_{s,s+n}(x,y) on adjacent pairs (M in {5,8}, several s, n)")
ok = True
cnt = 0
for m, Kk in [(4, 2), (6, 3), (9, 3)]:
    M = m + Kk
    for s in range(1, 200):
        if not (DL[s - 2] == 0 and DL[s - 1] == 1):
            continue
        kmax = 200
        P = R.kernel_to_site1(s, [m - 1, m], 1, M, DL)
        # find first n with P(m,1) >= P(m-1,1) > 0, then check V_k(m) >= V_k(m-1) for all k >= n+1 up to kmax
        first = None
        for n in range(1, 120):
            P = R.kernel_to_site1(s, [m - 1, m], n, M, DL)
            if P[m - 1] > 0 and P[m] >= P[m - 1]:
                first = n
                break
        if first is None:
            continue
        Vlo = R.free_counts(s, m - 1, kmax, M, DL)
        Vhi = R.free_counts(s, m, kmax, M, DL)
        ok &= all(Vhi[k] >= Vlo[k] for k in range(first + 1, kmax + 1))
        cnt += 1
check(ok and cnt > 30, "Lemma 3: P_{s,s+n}(m,1) >= P_{s,s+n}(m-1,1) > 0 implies V^(s)_k(m) >= V^(s)_k(m-1) for all k >= n+1 "
      "(%d landing times, k <= 200)" % cnt)

# 2d Lemma 4 ingredients, exhaustive on small bridges
def vpath(s, x0, bits, dl):
    v = x0
    out = []
    for i, b in enumerate(bits):
        v = v + b - dl[s + i]
        out.append(v)
    return out


bad = 0
cases = 0
for m, Kk in [(3, 2), (4, 2), (5, 3)]:
    M = m + Kk
    for s in range(1, 26):
        for n in range(4, 15):
            e = 1 - m + (FL[s + n] - FL[s])
            if e < 0 or e + 1 > n:
                continue
            tot = sf = sb = sfb = sumCf = 0
            for Aset in itertools.combinations(range(n), e):
                bits = [0] * n
                for p in Aset:
                    bits[p] = 1
                pth = vpath(s, m, bits, DL)
                f = Aset[0] if e > 0 else n
                bok = all(v >= 1 for v in pth)
                tok = all(v <= M for v in pth[:-1])
                tot += 1; sf += f; sb += bok; sfb += f * bok
                if bok and tok:
                    sumCf += f
            prev = 0
            for Aset in itertools.combinations(range(n), e + 1):
                bits = [0] * n
                for p in Aset:
                    bits[p] = 1
                pth = vpath(s, m - 1, bits, DL)
                if pth[-1] == 1 and all(v >= 1 for v in pth) and all(v <= M for v in pth[:-1]):
                    prev += 1
            cases += 1
            if Fraction(sf, tot) != Fraction(n - e, e + 1):
                bad += 1
            if sfb * tot > sf * sb:
                bad += 1
            if prev > sumCf:
                bad += 1
check(bad == 0 and cases > 500,
      "Lemma 4: E[f]=(n-e)/(e+1), FKG E[f 1_BOK] <= E[f]P(BOK), injection P^n(m-1,1) <= sum_C f  (%d bridge cases, exhaustive)" % cases)

# 2e Lemma 5 (real cycle lemma), exhaustive over rotations for random multisets
random.seed(20260926)
bad = 0
for trial in range(4000):
    n = random.randint(2, 32)
    kk = random.randint(0, n)
    xs = [C] * kk + [C - 1] * (n - kk)
    random.shuffle(xs)
    T = sum(xs)
    if T <= 1e-12:
        continue
    good = 0
    for j in range(n):
        ps = 0.0
        okr = True
        for i in range(n):
            ps += xs[(j + i) % n]
            if ps <= 1e-12:
                okr = False
                break
        good += okr
    if good < math.ceil(T / C - 1e-9):
        bad += 1
check(bad == 0, "Lemma 5 (real cycle lemma): >= ceil(T/c) good rotations for 4000 random {c, -d} sequences")

# 2f Lemmas 5-6 on exact bridge probabilities
bad = 0
mn = 1e9
for m, Kk in [(6, 3), (10, 5), (20, 8)]:
    M = m + Kk
    for s in range(1, 70, 5):
        for n in range(int(0.9 * (m - 1) / MU), int(1.3 * (m - 1) / MU), 3):
            e = 1 - m + (FL[s + n] - FL[s])
            if e < 0 or e > n:
                continue
            a = m - (C * s - FL[s])
            b = 1 - (C * (s + n) - FL[s + n])
            tot = math.comb(n, e)
            pbok = R.bridge_counts(s, m, 1, n, M, DL, True, False) / tot
            pth = 1 - R.bridge_counts(s, m, 1, n, M, DL, False, True) / tot
            lb = math.ceil((a - b) / C) / n
            v = (a - b) / n
            ub = sum(math.exp(-2 * (M - a + v * t) ** 2 / min(t, n - t)) for t in range(1, n))
            bad += (pbok < lb - 1e-12) + (pth > ub + 1e-12)
            mn = min(mn, pbok / lb)
check(bad == 0, "Lemmas 5-6 on exact bridges: P(BOK) >= ceil(T/c)/n and P(TH) <= Hoeffding sum (min ratio P(BOK)/LB = %.3f)" % mn)

# ------------------------------------------------------------------------------------------------ 3
section("3. eventual monotonicity (EM), exact part: 2 <= m < M_E, all landing Sturmian factors of length n1(m)")
ok = True
nf_tot = 0
for m in range(2, M_E):
    M = m + K
    n = n1(m)
    lf, complete = R.landing_factors(n, DL, 79000 - n - 5)
    ok &= complete
    worst = None
    for word in lf:
        P = R.kernel_to_site1_word(word, [m - 1, m], M)
        ok &= (P[m - 1] > 0 and P[m] >= P[m - 1])
        r = Fraction(P[m], P[m - 1])
        worst = r if worst is None or r < worst else worst
        nf_tot += 1
    print("  m=%2d n1=%3d landing factors=%3d  min P(m,1)/P(m-1,1) = %.4f" % (m, n, len(lf), float(worst)))
check(ok, "EM exact: for 2<=m<%d and every landing factor, P_{s,s+n1}(m,1) >= P_{s,s+n1}(m-1,1) > 0 "
      "(%d factors; factor sets certified complete: n+3 factors of length n+2 found)" % (M_E, nf_tot))

# ------------------------------------------------------------------------------------------------ 4
section("4. EM, analytic part: m >= M_E")
v_star = (MU - W) * (M_E - 2) / (M_E + 2)
G_star = R.G_upper(K, v_star)
lhs = (0.5 - W) / (0.5 + W) + 2 * C * G_star / v_star
print("  v_*(M_E) = %.6f   G(K,v_*) <= %.4e   (1/2-W)/(1/2+W) = %.6f   2cG/v_* = %.4e   sum = %.6f"
      % (v_star, G_star, (0.5 - W) / (0.5 + W), 2 * C * G_star / v_star, lhs))
check(lhs < 1, "EM criterion (1/2-W)/(1/2+W) + 2cG(K,v_*)/v_* < 1 at m = M_E, hence for all m >= M_E (monotone in m)")
# sanity: the ceiling property m+1 <= (mu-W) n1(m) used in the proof
check(all((MU - W) * n1(m) >= m + 1 for m in range(2, 200000)), "n1(m)(mu - W) >= m + 1 for 2 <= m < 2e5 (ceiling, with guard)")

# ------------------------------------------------------------------------------------------------ 5
section("5. top-hit bound (TB): q_k = P(TH'|BOK) for the free walk from site m, k <= k1(m) = n1(m)+1")
Q3K = 3.0 ** (-K)
# 5a: k <= k_b(m): Kolmogorov gives P(BOK) >= 3/4, Doob gives P(TH') <= 3^-K
q_a = Q3K / 0.75
print("  5a (k <= k_b(m), all m): q <= 3^-K/(3/4) = %.3e" % q_a)
check(q_a < 1e-6, "5a: q <= (4/3) 3^-K on k <= k_b(m)")
# 5b: m <= M_X, window k in (k_b, k1]: exact worst-phase survival beta_m(k) = H^(0)_k(m-1)/2^k
q_b = 0.0
arg_b = None
for m in range(2, M_X + 1):
    k1 = n1(m) + 1
    kb = R.k_b(m)
    Hc = R.free_counts(0, m - 1, k1, None, DL)
    for k in range(kb + 1, k1 + 1):
        beta = Fraction(Hc[k], 2 ** k)
        qq = Q3K / float(beta)
        if qq > q_b:
            q_b, arg_b = qq, (m, k, float(beta))
print("  5b (m <= %d): max q <= 3^-K/beta_m(k) = %.3e at (m,k,beta) = %s" % (M_X, q_b, arg_b))
check(q_b < 1e-5, "5b: exact worst-phase survival bound, q <= %.2e for 2 <= m <= %d on the window" % (q_b, M_X))
# 5c: M_X < m < M_2: closed form at k = k1(m) (survival is nonincreasing in k, so this covers the whole window)
q_c = 0.0
arg_c = None
for m in range(M_X + 1, M_2):
    k1 = n1(m) + 1
    best = 0.0
    for Dd in range(3, 30):
        j0 = FL[k1] + Dd - m + 2          # U_k > c k + D - m + 1  <=>  U_k >= j0
        lb = R.binom_tail_lower(k1, j0) - 3.0 ** (-Dd)
        best = max(best, lb)
    qq = Q3K / best if best > 0 else float('inf')
    if qq > q_c:
        q_c, arg_c = qq, (m, k1, best)
print("  5c (%d < m < %d): max q <= 3^-K / [P_lb(U_k1 >= j0) - 3^-D] = %.3e at (m,k1,P_lb) = %s" % (M_X, M_2, q_c, arg_c))
check(q_c < 1e-4, "5c: closed-form survival bound, q <= %.2e for %d < m < %d on the window" % (q_c, M_X, M_2))

# 5d: m >= M_2: bridge decomposition. q <= h + (c/v_min) R(m)
h = 2 * C * R.G_upper(K, V_MIN) / V_MIN
KLv = R.KL_half(C - V_MIN)
print("  5d: v_min=%.3f  h = 2cG(K,v_min)/v_min = %.4e   KL(c-v_min||1/2) = %.5f" % (V_MIN, h, KLv))
check(C - V_MIN > 0.5, "5d: c - v_min > 1/2 (Chernoff applies)")
qd_max = 0.0
arg_d = None
ok = True
for m in range(M_2, M_BIG + 1):
    k1 = n1(m) + 1
    kb = R.k_b(m)
    ymax = max(abs(MU - (m - 2) / k1), abs(m / (kb + 1) - MU))
    klu = R.KL_half_upper(ymax)
    if not (klu < KLv and m - 2 >= V_MIN * k1):
        ok = False
        break
    Rm = math.sqrt(k1) / 0.6753 * math.exp(-(kb + 1) * (KLv - klu))
    qd = h + C / V_MIN * Rm
    if qd > qd_max:
        qd_max, arg_d = qd, (m, k1, kb, ymax, Rm)
check(ok, "5d: for M_2 <= m <= M_BIG: KL(1/2+y_max) < KL(c-v_min) and m-2 >= v_min k1 (e* is a good endpoint)")
print("  5d: max over M_2<=m<=M_BIG of h + (c/v_min)R(m) = %.4e at (m,k1,kb,ymax,R) = %s" % (qd_max, arg_d))
# majorant for m >= M_BIG: y_max <= Y (monotone pieces), k1 <= (m+1)/(mu-W)+2, kb+1 >= (m-1-sqrt((m-1)/mu))/mu
mb = M_BIG
Y = max(MU - (mb - 2) / ((mb + 1) / (MU - W) + 2), mb * MU / (mb - 1 - math.sqrt((mb - 1) / MU)) - MU)
Y = max(Y, W)          # the first piece decreases to W from above; W itself is also an upper bound for large m
klY = R.KL_half_upper(Y)


def Rmaj(m):
    return math.sqrt((m + 1) / (MU - W) + 2) / 0.6753 * math.exp(-((m - 1 - math.sqrt((m - 1) / MU)) / MU) * (KLv - klY))


Rb = Rmaj(M_BIG)
# d/dm log Rmaj = 1/(2(m+1+2(mu-W))) - (KLv-klY)(1 - 1/(2 sqrt(mu(m-1))))/mu < 0 for m >= M_BIG:
dlog = 1.0 / (2 * (M_BIG + 1)) - (KLv - klY) * (1 - 1 / (2 * math.sqrt(MU * (M_BIG - 1)))) / MU
check(klY < KLv and dlog < 0 and h + C / V_MIN * Rb < qd_max + 1e-12 and V_MIN / (MU - W) < 1
      and M_BIG - 2 >= V_MIN * ((M_BIG + 1) / (MU - W) + 2),
      "5d: for m > M_BIG the monotone majorant gives R(m) <= Rmaj(M_BIG) = %.3e (Y = %.4e, dlog/dm < 0; "
      "m-2 >= v_min k1 persists since v_min/(mu-W) < 1)" % (Rb, Y))
q_d = max(qd_max, h + C / V_MIN * Rb)
check(q_d < 0.01, "5d: bridge bound q <= %.4e for all m >= M_2" % q_d)

# ------------------------------------------------------------------------------------------------ 6
section("6. assembly: Gamma = 1/(1 - max q)")
q_max = max(q_a, q_b, q_c, q_d)
GAMMA_F = 1.0 / (1.0 - q_max)
GAMMA = Fraction(10039, 10000)          # stated constant 1.0039 >= 1/(1-q_max)
print("  max q = %.6e  ->  1/(1-max q) = %.7f ;  stated GAMMA = %s = %.4f" % (q_max, GAMMA_F, GAMMA, float(GAMMA)))
check(float(GAMMA) >= GAMMA_F, "Theorem R_15: N_m(L) <= %s * A_(m+15)(L) for all m >= 2, L >= 1 (all ingredients checked)" % GAMMA)

# ------------------------------------------------------------------------------------------------ 7
section("7. exact verification of the proved inequality and the observed ratio profile")
ok = True
worst15 = Fraction(0)
prof = {1: [], 2: [], 3: [], 15: []}
for m in range(2, 61):
    Lm = 2500 if m <= 20 else (1500 if m <= 40 else 400)
    N = R.robin_counts(Lm, m, DL)
    for c0 in (1, 2, 3, 15):
        A = R.dirichlet_counts(Lm, m + c0, DL)
        rmax = Fraction(0)
        argL = None
        for L in range(1, Lm + 1):
            r = Fraction(N[L], A[L])
            if r > rmax:
                rmax, argL = r, L
            if c0 == 15:
                ok &= (N[L] * GAMMA.denominator <= GAMMA.numerator * A[L])
        prof[c0].append((m, rmax, argL))
        if c0 == 15:
            worst15 = max(worst15, rmax)
check(ok, "N_m(L) <= %s A_(m+15)(L) exactly for 2<=m<=20, L<=2500; 21<=m<=40, L<=1500; 41<=m<=60, L<=400 "
      "(max observed ratio 1 + %.3e)" % (GAMMA, float(worst15 - 1)))
for c0 in (1, 2, 3, 15):
    best = max(prof[c0], key=lambda z: z[1])
    print("  c0=%2d: max_{m<=60,L} N_m(L)/A_(m+c0)(L) = 1 + %.4e at (m,L) = (%d,%d)" % (c0, float(best[1] - 1), best[0], best[2]))
print("  profile (c0=1): " + ", ".join("m=%d:%.6f" % (m, float(r)) for (m, r, L) in prof[1][:19]))
print("  profile (c0=2): " + ", ".join("m=%d:%.7f" % (m, float(r)) for (m, r, L) in prof[2][:19]))
print("  profile (c0=15), excess over 1 where positive: " + (", ".join("m=%d:%.2e" % (m, float(r - 1)) for (m, r, L) in prof[15] if r > 1) or "none"))

# ------------------------------------------------------------------------------------------------ 8
section("8. consistency checks (not used by the proof)")
# 8a: the worst-phase survival beta_m(k) = H^(0)_k(m-1)/2^k is below H^(u)_k(m)/2^k for sampled phases u
ok = True
for m in (3, 5, 9):
    kk = n1(m) + 1
    beta = R.free_counts(0, m - 1, kk, None, DL)
    for u in range(1, 400, 7):
        Hu = R.free_counts(u, m, kk, None, DL)
        ok &= all(Hu[k] >= beta[k] for k in range(kk + 1))
check(ok, "8a: H^(u)_k(m) >= H^(0)_k(m-1) (worst phase = S-height m-1) for sampled u, m in {3,5,9}")
# 8b: EM directly at sampled landing times for some m >= M_E, K = 15
ok = True
cnt = 0
for m in (26, 30, 40):
    M = m + K
    kk = n1(m) + 1
    for s in range(2, 3000):
        if not (DL[s - 2] == 0 and DL[s - 1] == 1):
            continue
        if cnt % 5 == 0:
            Vlo = R.free_counts(s, m - 1, kk + 60, M, DL)
            Vhi = R.free_counts(s, m, kk + 60, M, DL)
            ok &= all(Vhi[k] >= Vlo[k] for k in range(kk, kk + 61))
        cnt += 1
        if cnt > 400:
            break
check(ok, "8b: V^(s)_k(m) >= V^(s)_k(m-1) for k in [k1, k1+60] at sampled landing times, m in {26,30,40}, K=15")
# 8c: exact top-hit ratio sup_{u,k<=k1} H_k(m)/V_k(m) at sampled phases vs the proved Gamma
worst = 1.0
for m in (5, 10, 20, 40):
    M = m + K
    kk = n1(m) + 1
    for u in range(1, 300, 3):
        Hu = R.free_counts(u, m, kk, None, DL)
        Vu = R.free_counts(u, m, kk, M, DL)
        worst = max(worst, max(Hu[k] / Vu[k] for k in range(1, kk + 1)))
print("  8c: sampled max H_k(m)/V_k(m) over u, k<=k1(m), m in {5,10,20,40}: %.3e above 1 (proved bound: %.3e)"
      % (worst - 1, GAMMA_F - 1))
check(worst <= GAMMA_F, "8c: sampled top-hit ratios lie below the proved Gamma")

# 8d: Proposition 1 at small shifts (sampled Gamma*, observed EM onset): N <= Gamma* A on L <= 700
ok = True
rows = []
for m, Kk in [(6, 3), (8, 3), (10, 4), (14, 5)]:
    M = m + Kk
    nn = 0
    zts = [t for t in range(1, 1500) if DL[t - 1] == 0 and DL[t] == 1]
    for t in zts[:120]:
        s_ = t + 1
        first = None
        for n in range(1, 401):
            P = R.kernel_to_site1(s_, [m - 1, m], n, M, DL)
            if P[m - 1] > 0 and P[m] >= P[m - 1]:
                first = n
                break
        nn = max(nn, first)
    kk = nn + 1
    Gs = 1.0
    for u in range(1, 400):
        for x in range(1, m + 1):
            if x == m and not (DL[u - 1] == 0 and DL[u] == 1):
                continue
            Wv = R.robin_from(u, x, kk, m, DL)
            Vv = R.free_counts(u, x, kk, M, DL)
            Gs = max(Gs, max(Wv[k] / Vv[k] for k in range(1, kk + 1)))
    N = R.robin_counts(700, m, DL)
    A = R.dirichlet_counts(700, M, DL)
    worst = max(N[L] / A[L] for L in range(1, 701))
    ok &= worst <= Gs + 1e-12
    rows.append("(m=%d,K=%d: k1=%d, Gamma*~%.6f, max N/A=%.6f)" % (m, Kk, kk, Gs, worst))
print("  8d: " + " ".join(rows))
check(ok, "8d: Proposition 1 bound N <= Gamma* A holds with sampled Gamma* at small shifts (L <= 700)")

# ------------------------------------------------------------------------------------------------ 9
section("9. small shifts c0 in {1,2,3}, bounded m: exact constants valid for ALL L (Proposition 1 + Lemma 3)")
DB = bytes(DL)


def factor_set(n, scan):
    seen = {}
    for s_ in range(scan):
        w_ = DB[s_:s_ + n]
        if w_ not in seen:
            seen[w_] = s_
    return seen


def em_onset(word, m, M):
    fa = [0] * (M + 3); fb = [0] * (M + 3); fa[m - 1] = 1; fb[m] = 1
    for i, d in enumerate(word):
        ga = [0] * (M + 3); gb = [0] * (M + 3)
        for v in range(1, M + 1):
            if fa[v]:
                ga[v - d] += fa[v]; ga[v + 1 - d] += fa[v]
            if fb[v]:
                gb[v - d] += fb[v]; gb[v + 1 - d] += fb[v]
        for g in (ga, gb):
            g[0] = 0; g[M + 1] = 0; g[M + 2] = 0
        fa, fb = ga, gb
        if fa[1] > 0 and fb[1] >= fa[1]:
            return i + 1
    return None


def counts_word(word, x, top, robin_m=None):
    """k -> V (top wall, relaxed last step) or, if robin_m is given, the Robin count W"""
    cap = robin_m if robin_m is not None else top
    f = [0] * (cap + 3); f[x] = 1; out = [1]
    for d in word:
        g = [0] * (cap + 3)
        for v in range(1, cap + 1):
            y = f[v]
            if not y:
                continue
            if robin_m is not None and v == robin_m:
                g[v - 1] += 2 * y
            else:
                g[v - d] += y; g[v + 1 - d] += y
        g[0] = 0
        if robin_m is None:
            out.append(sum(g[1:])); g[cap + 1] = 0; g[cap + 2] = 0
        else:
            out.append(sum(g))
        f = g
    return out


SMALL_M = 24
small_table = {}
ok = True
for c0 in (1, 2, 3):
    for m in range(2, SMALL_M + 1):
        M = m + c0
        Nmax = int(6 * (m - 1) / MU) + 60
        fs = factor_set(Nmax + 2, 79000 - Nmax)
        ok &= (len(fs) == Nmax + 3)
        n1s = 0
        for w_ in fs:
            if w_[0] == 0 and w_[1] == 1:
                cr = em_onset(w_[2:], m, M)
                ok &= (cr is not None)
                n1s = max(n1s, cr if cr is not None else 10 ** 9)
        k1s = n1s + 1
        fs2 = factor_set(k1s + 1, 79000 - k1s)
        ok &= (len(fs2) == k1s + 2)
        G = Fraction(1)
        for w_ in fs2:
            prev, word = w_[0], w_[1:]
            for x in range(1, m + 1):
                if x == m and not (prev == 0 and word[0] == 1):
                    continue
                Wv = counts_word(word, x, M, robin_m=m)
                Vv = counts_word(word, x, M)
                for k in range(1, k1s + 1):
                    if Wv[k] * G.denominator > G.numerator * Vv[k]:
                        G = Fraction(Wv[k], Vv[k])
        small_table[(c0, m)] = (n1s, G)
        # sanity: the implied inequality on L <= 600
        N = R.robin_counts(600, m, DL)
        A = R.dirichlet_counts(600, M, DL)
        ok &= all(N[L] * G.denominator <= G.numerator * A[L] for L in range(1, 601))
for c0 in (1, 2, 3):
    Gmax = max(small_table[(c0, m)][1] for m in range(2, SMALL_M + 1))
    print("  c0=%d: max_{2<=m<=%d} Gamma*(m) = %.6f ;  " % (c0, SMALL_M, float(Gmax)) +
          ", ".join("m=%d:%.4f(n1=%d)" % (m, float(small_table[(c0, m)][1]), small_table[(c0, m)][0])
                    for m in range(2, SMALL_M + 1, 2)))
check(ok, "9: for c0 in {1,2,3} and 2<=m<=%d: EM onset over all landing factors and Gamma* over all phases computed "
      "exactly (factor sets complete); hence N_m(L) <= Gamma*(m) A_(m+c0)(L) for ALL L" % SMALL_M)

section("summary")
ru = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
print("peak RSS (ru_maxrss, bytes on macOS) = %d  (%.0f MB)" % (ru, ru / 2 ** 20))
print("checks passed: %d   wall time %.1fs" % (NCHECK, time.time() - T_START))
print("ALL CHECKS PASSED")
