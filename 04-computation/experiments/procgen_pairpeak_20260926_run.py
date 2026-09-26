#!/usr/bin/env python3
"""procgen_pairpeak_20260926_run -- runner for the pairing-peak lane (note:
05-knowledge/results/procgen_pairpeak_20260926_pairing_peak_price.md).  Every printed claim is a check(...) that aborts
on failure.  Deterministic apart from the [time] lines.

Sections
  0  constants
  1  coupling lemma of a single flip (exact chain, merge law, merge rates of bad sources)
  2  bijection lemma for prefix-determined flip rules; slope survival = actual non-descent; private realization
  3  Theorem A (failure bound of the reflected barrier) against exact DP counts
  4  the private price: 2 rho^peak <= pi_L <= DP bounds; elementary exponent gamma; sampled private barrier cost
  5  the Robin inequality N_m(L) <= K0 A_(m+1)(L) (exact to L = 180, float to L = 1024); conditional bound
  6  the consistent construction (C program): valid sections of members of P_L, densities, interference
  7  one-flip domination statistics
"""
import math
import os
import random
import re
import resource
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import procgen_pairpeak_20260926_lib as P
import procgen_peak_20260926_lib as PK

check = P.check
T0 = time.time()


def rss_mb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1024 * 1024) if sys.platform == "darwin" else r / 1024


def stamp(label):
    print("[time] %s: %.1fs elapsed, peak RSS %.0f MB" % (label, time.time() - T0, rss_mb()), flush=True)


def header(s):
    print("\n" + "=" * 100 + "\n" + s + "\n" + "=" * 100, flush=True)


# ------------------------------------------------------------------------------------------------------------------
header("0. Constants")
c, s2, lam, eta = P.C, P.SIG2, P.LAM, P.ETA
a_ln3 = math.log(3)
b_blk = math.log(1 / P.RHO_B) * s2 / 2
GAMMA = 1.5 * (2 * b_blk) ** (1 / 3) * a_ln3 ** (2 / 3)
kap3 = PK.kappa(3)[0]
print("   c = log_3 2 = %.6f  sigma^2 = %.6f  lambda = %.6f  1-H(c) = %.6f  rho_b = 1-(1-lambda^2)/14 = %.6f" % (c, s2, lam, eta, P.RHO_B))
print("   elementary exponent gamma = (3/2) (ln(1/rho_b) sigma^2)^(1/3) (ln 3)^(2/3) = %.5f ; Mogulskii kappa_3 = %.5f" % (GAMMA, kap3))
check(abs(eta - 0.050044) < 1e-6, "1 - H(log_3 2) = 0.050044")
check(abs(P.RHO_B - 0.953013) < 1e-6 and lam ** 2 < 1, "rho_b = 0.953013 < 1")
check(abs(GAMMA - 0.35742) < 1e-4, "gamma = 0.3574")

# ------------------------------------------------------------------------------------------------------------------
header("1. Coupling lemma of a single flip: y_t = 3^(a_t) x_t + b_t; merge iff the Collatz word of y_1 begins 1^r 0 0")
rng = random.Random(20260926)
nmerge = 0
for _ in range(3000):
    v = rng.getrandbits(400) | 1
    mt, a0, steps = P.chain_verify(v, 250)
    fo = P.first_opportunity_merge(v)
    # merge time r+3 exactly when the word begins 1^r 0 0
    y = (3 * v + 1) // 2
    r = 0
    while y & 1:
        y = P.T(y)
        r += 1
    if fo:
        nmerge += 1
        assert mt == r + 3, (v, mt, r)
check(True, "for 3000 random odd v < 2^400 the identity y_t = 3^(a_t) x_t + b_t (four-case chain from (a,b) = (1,2)) and "
      "parity(x_t) = parity(y_t) + b_t hold at every step until a_t = 0 (up to 250 steps); whenever the word of y_1 "
      "begins 1^r 0 0 the orbits merge exactly at time r+3 (%d merges)" % nmerge)
K = 16
cnt = 0
for v in range(1, 1 << K, 2):
    # the first K-1 letters of the word of y_1 are fixed by v mod 2^K; use a large lift
    vv = v + (1 << (K + 60)) * 12345
    yy = (3 * vv + 1) // 2
    w = []
    for _ in range(K - 1):
        w.append(yy & 1)
        yy = P.T(yy)
    r = 0
    while r < len(w) and w[r] == 1:
        r += 1
    if r + 2 <= K - 1 and w[r] == 0 and w[r + 1] == 0:
        cnt += 1
check(cnt == 2 ** (K - 2) - 1, "exact merge law: among the 2^%d odd classes v mod 2^%d, exactly 2^%d - 1 have a word of y_1 beginning "
      "1^r00 with r+2 <= %d (Haar probability 1/2 in the limit; under the tilted measure Q it is 1-c = %.4f)" % (K - 1, K, K - 2, K - 1, 1 - c))

# merge rates of bad sources: flip at the odd point before the peak (the peak-lane candidate)
print("   single flip at the highest odd point of a bad source (L steps; n random ~ 2^90, rejection sampling):")
for L in (16, 32, 48):
    rng2 = random.Random(1000 + L)
    tot = ok = mer = 0
    while tot < 1500:
        n = rng2.getrandbits(90) | (1 << 89)
        p = [n]
        for _ in range(L):
            p.append(P.T(p[-1]))
        if not all(p[j] > n for j in range(1, L + 1)):
            continue
        tot += 1
        k = max((j for j in range(L) if p[j] & 1), key=lambda j: p[j])
        x = (p[k] - 1) >> 1
        res = None
        for t in range(k + 1, L + 1):
            if x < n:
                res = "ok"
                break
            if x == p[t]:
                res = "merge"
                break
            x = P.T(x)
        ok += (res == "ok")
        mer += (res == "merge")
    print("     L=%d: rescued %.3f, re-merged with the bad Collatz orbit %.3f, other failure %.3f" % (L, ok / tot, mer / tot, 1 - (ok + mer) / tot))
    check(mer / tot > 0.35, "L=%d: the single high flip re-merges with the undecided orbit for more than 35%% of bad sources (EMPIRICAL)" % L)
stamp("section 1")

# ------------------------------------------------------------------------------------------------------------------
header("2. Bijection lemma for the reflected barrier Pi^(m); slope survival = non-descent; private realization")
for L in (8, 10, 12, 14):
    M = L + 8
    fl = P.floor_table(M)
    for m in (2, 3, 4, 5):
        Kbig = (random.Random(L * 10 + m).getrandbits(50) | 1) << (L + 40)
        words = set()
        surv = 0
        agree = True
        for r in range(1 << L):
            n = r + Kbig
            w, s, flips, desc, cost, nfd = P.barrier_run(n, L, m, fl, M)
            words.add(tuple(w))
            surv += s
            if s == desc:
                agree = False
        exact = P.refl_count(L, m, fl, M)
        check(len(words) == 1 << L and surv == exact and agree,
              "L=%d m=%d: the 2^%d classes give 2^%d distinct modified words; survivors %d = DP count N_m(L) = %d; slope survival <=> no actual descent" % (L, m, L, L, surv, exact))
# private realization: the member with exactly the flipped pairs of Pi^(m) reproduces the path (large n)
L, m = 14, 3
M = L + 8
fl = P.floor_table(M)
Kbig = 987654321 << (L + 40)
bad_real = 0
for r in range(1 << L):
    n = r + Kbig
    w, s, flips, desc, cost, nfd = P.barrier_run(n, L, m, fl, M)
    F = set((x + 1) // 2 for (t, x) in flips)
    x = n
    xs = [n]
    for t in range(L):
        i = (x + 1) // 2
        if i in F:
            x = (x - 1) // 2 if (x & 1) else 3 * (x // 2)
        else:
            x = P.T(x)
        xs.append(x)
    y = n
    ys = [n]
    U = 0
    for t in range(L):
        ft = 0 if t == 0 else int(fl[M + t])
        thr = (m - 1) if t == 0 else (m - 1 + ft + 1)
        if (y & 1) and U >= thr:
            y = (y - 1) // 2
        elif y & 1:
            y = (3 * y + 1) // 2
            U += 1
        else:
            y //= 2
        ys.append(y)
    if xs != ys:
        bad_real += 1
check(bad_real == 0, "L=14, m=3, all 2^14 classes (n ~ 2^80): the pairing member with exactly the pairs flipped by Pi^(m) reproduces the Pi^(m) path")
stamp("section 2")

# ------------------------------------------------------------------------------------------------------------------
header("3. Theorem A: N_m(L)/2^L <= 2^(-(1-H)L) rho_b^floor(L/n_m), n_m = ceil(2(m+1)^2/sigma^2)")
LMAX = 180
M = LMAX + 8
fl = P.floor_table(M)
NM = {}
worst = 0.0
for m in range(1, LMAX + 2):
    NM[m] = P.refl_counts_all_L(LMAX, m, fl, M)
    for L in range(1, LMAX + 1):
        lhs = NM[m][L] / 2 ** L
        rhs = P.block_bound(L, m)
        worst = max(worst, lhs / rhs)
check(worst <= 1.0, "exact DP, 1 <= L <= %d, 1 <= m <= L+1: N_m(L)/2^L <= block bound (max ratio %.4f)" % (LMAX, worst))
nb = [PK.rho_exact(3, L) for L in (20, 60, 180)]
check(NM[25][20] == nb[0] and NM[65][60] == nb[1] and NM[LMAX + 1][LMAX] <= nb[2],
      "barrier above the reachable range: N_m(L) = |Bad_L| (L = 20, m = 25; L = 60, m = 65); N_m(L) <= |Bad_L| always")
check(NM[1][LMAX] == 0, "m = 1 (barrier at level 0): no class survives")
print("   exact N_m(L)/2^L versus |Bad_L|/2^L:")
for L in (60, 120, 180):
    rho = PK.rho_exact(3, L) / 2 ** L
    print("     L=%d  rho_L=%.4e  " % (L, rho) + "  ".join("m=%d:%.3e" % (m, NM[m][L] / 2 ** L) for m in (2, 3, 4, 5, 6, 8)))
for L in (256, 512, 1024):
    Mf = L + 8
    flf = P.floor_table(Mf)
    wr = 0.0
    for m in range(2, 41):
        pf, _ = P.refl_float(L, m, flf, Mf)
        wr = max(wr, pf / P.block_bound(L, m))
    check(wr <= 1.0, "float DP, L=%d, 2 <= m <= 40: N_m(L)/2^L <= block bound (max ratio %.3e)" % (L, wr))
stamp("section 3")

# ------------------------------------------------------------------------------------------------------------------
header("4. The private price: 2 rho^peak_L <= pi_L <= single/multi-scale DP bounds; elementary exponent gamma")
rows = []
for L in (32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024):
    Mf = 2 * L + 8
    flf = P.floor_table(Mf)
    rho = PK.rho_float(3, L, flf, Mf)
    rp = PK.rho_peak_float(3, L, smax=min(70.0, 3 * L ** (1 / 3) + 12))
    rpk = rp["hi"]
    fails = {}
    EJ = {}
    for m in range(2, min(L, 60)):
        fails[m], EJ[m] = P.refl_float(L, m, flf, Mf)
    single = min(2 * 3 ** (1 - m) * min(L * rho, EJ[m]) + 2 * fails[m] for m in fails)
    mtop = max(fails)
    multi = 2 * L * (3 ** (1 - mtop) * rho + sum(3 ** (1 - m) * fails[m + 1] for m in range(2, mtop))) + 2 * fails[2]
    rows.append((L, rho, rpk, single, multi))
    check(2 * rp["lo"] <= min(single, multi) + 1e-300, "L=%d: lower bound 2 rho^peak = %.3e <= upper bound %.3e" % (L, 2 * rpk, min(single, multi)))
print("      L   rho_L        2rho^peak    single-scale bound  /rho   multi-scale bound  /rho    (2rho^peak)/rho")
for (L, rho, rpk, single, multi) in rows:
    print("   %5d  %.4e  %.4e  %.4e  %8.4f  %.4e  %8.4f   %.3e" % (L, rho, 2 * rpk, single, single / rho, multi, multi / rho, 2 * rpk / rho))
check(rows[-1][4] / rows[-1][1] < 0.01 and rows[-3][4] / rows[-3][1] < 1.0,
      "the rigorous multi-scale bound is below rho_L at L = 512 and below 0.01 rho_L at L = 1024 (float DP)")
# elementary asymptotic bound: 2^((1-H)L) * B_L^elem  versus exp(-gamma L^(1/3))
print("   elementary bound (Theorem A inserted in the multi-scale bound), normalised: ln(2^((1-H)L) B_L) / L^(1/3):")
vals = []
for L in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7, 10 ** 8, 10 ** 9, 10 ** 10):
    terms = []
    for m in range(2, int(3 * L ** (1 / 3)) + 40):
        terms.append(math.log(2 * L) + (1 - m) * a_ln3 + P.block_bound_log(L, m + 1) + eta * L * math.log(2))
    lse = max(terms) + math.log(sum(math.exp(t - max(terms)) for t in terms))
    vals.append((L, lse / L ** (1 / 3)))
    print("     L=%.0e : %.4f   (target -gamma = %.4f)" % (L, lse / L ** (1 / 3), -GAMMA))
check(vals[-1][1] < -0.33 and vals[-1][1] > -GAMMA and all(vals[i + 1][1] < vals[i][1] for i in range(len(vals) - 1)),
      "the normalised elementary bound decreases monotonically toward -gamma (it is -%.4f at L = 1e10)" % -vals[-1][1])
# sampled private cost of the multi-scale barrier strategy
print("   sampled private cost of the multi-scale barrier strategy (uniform bad classes, n ~ 2^(L+80)):")
for L, K in ((32, 1500), (64, 1500), (128, 700)):
    Mf = 2 * L + 8
    flf = P.floor_table(Mf)
    sample, nbad = P.bad_word_sampler(L, flf, Mf)
    rho = nbad / 2 ** L
    rp = PK.rho_peak_float(3, L, smax=min(70.0, 3 * L ** (1 / 3) + 12))["hi"]
    rngs = random.Random(7 * L)
    tot = 0.0
    pk = 0.0
    fails = 0
    nfl = []
    for _ in range(K):
        w = sample(rngs)
        n = P.n_from_word(w, rngs)
        x = n
        mx = n
        for t in range(L - 1):
            x = P.T(x)
            mx = max(mx, x)
        pk += 2.0 * n / mx
        best = None
        for m in range(45, 1, -1):
            _, s, flips, desc, cost, nfd = P.barrier_run(n, L, m, flf, Mf)
            if desc:
                best = (cost, nfd)
                break
        if best is None:
            tot += 2.0
            fails += 1
        else:
            tot += best[0]
            nfl.append(best[1])
    est = rho * tot / K
    print("     L=%d: rho_L=%.3e 2rho^peak=%.3e (sample %.3e); private barrier cost %.3e = %.2f x 2rho^peak = %.4f x rho_L; fallbacks %d/%d; mean flips %.2f"
          % (L, rho, 2 * rp, rho * pk / K, est, est / (2 * rp), est / rho, fails, K, sum(nfl) / max(1, len(nfl))))
    check(abs(rho * pk / K / (2 * rp) - 1) < 0.1, "L=%d: the sample reproduces 2 rho^peak within 10%%" % L)
    check(est / (2 * rp) < 12, "L=%d: sampled private barrier cost < 12 x 2 rho^peak (EMPIRICAL)" % L)
stamp("section 4")

# ------------------------------------------------------------------------------------------------------------------
header("5. The Robin inequality N_m(L) <= K0 A_(m+K)(L) and the conditional private bound")
worst1 = (0.0, None)
worst2 = (0.0, None)
for m in range(1, LMAX + 2):
    A1 = P.hard_counts_all_L(LMAX, m + 1, fl, M)
    A2 = P.hard_counts_all_L(LMAX, m + 2, fl, M)
    for L in range(1, LMAX + 1):
        if A1[L]:
            r1 = NM[m][L] / A1[L]
            if r1 > worst1[0]:
                worst1 = (r1, (L, m))
        elif NM[m][L]:
            worst1 = (float("inf"), (L, m))
        if A2[L]:
            r2 = NM[m][L] / A2[L]
            if r2 > worst2[0]:
                worst2 = (r2, (L, m))
check(worst1[0] <= 1.0033, "exact, 1 <= L <= %d, all m: N_m(L) <= 1.0033 A_(m+1)(L) (max %.6f at (L,m) = %s)" % (LMAX, worst1[0], worst1[1]))
check(worst2[0] <= 1.0001, "exact, 1 <= L <= %d, all m: N_m(L) <= 1.0001 A_(m+2)(L) (max %.7f at %s); the ratio does exceed 1" % (LMAX, worst2[0], worst2[1]))
check(worst1[0] > 1.0, "the ratio N_m/A_(m+1) exceeds 1 somewhere (so K0 = 1 is false for K = 1)")
for L in (256, 512, 1024):
    Mf = L + 8
    flf = P.floor_table(Mf)
    wr = 0.0
    for m in range(2, 61):
        pf, _ = P.refl_float(L, m, flf, Mf)
        ha = P.hard_float(L, m + 1, flf, Mf)
        if ha > 0:
            wr = max(wr, pf / ha)
    check(wr <= 1.01, "float DP, L=%d, 2 <= m <= 60: N_m(L) <= 1.01 A_(m+1)(L) (max ratio %.5f)" % (L, wr))
print("   conditional private bound (Theorem C with K = 1, K0 = 1.0033): pi_L <= K0 3 (27 L + 18) rho^peak_L + 2L 2^(1-L) rho^peak_L")
stamp("section 5")

# ------------------------------------------------------------------------------------------------------------------
header("6. The consistent construction (single member of the pairing family), valid on n <= N")
src = os.path.join(HERE, "procgen_pairpeak_20260926_greedy.c")
tmpd = tempfile.mkdtemp(prefix="pairpeak_")
exe = os.path.join(tmpd, "greedy")
subprocess.run(["cc", "-O2", "-o", exe, src, "-lm"], check=True)


def run_c(L, N, W, mode):
    out = subprocess.run([exe, str(L), str(N), str(W), str(mode)], check=True, capture_output=True, text=True).stdout
    d = {}
    for line in out.splitlines():
        for k, v in re.findall(r"(\w+)=([-\d.e+]+)", line):
            d[k] = float(v) if ("." in v or "e" in v) else int(v)
        if line.startswith("SCALES"):
            d["scales"] = line[7:]
        if line.startswith("IDEAL"):
            d["ideal"] = line[6:]
    return d


peak_lane = {8: 0.031274, 12: 0.017886, 16: 0.006902, 20: 0.004658, 24: 0.002284, 28: 0.001542, 32: 0.000908}
thm4475 = {8: 0.068052, 12: 0.049190, 16: 0.028426, 20: 0.022694}
for L, ref in sorted(thm4475.items()):
    d = run_c(L, 10 ** 6, 12, 2)
    check(abs(d["density"] - ref) < 2e-6 and d["stuck"] == 0 and d["verifyfail"] == 0,
          "MODE 2 (THM-4475's G_L) L=%d: density %.6f = the peak lane's re-implementation %.6f (simulator check)" % (L, d["density"], ref))
print("\n   MODE 0 (all certificates frozen), N = 10^6, scales W = 12..2 then THM-4475 fallback:")
print("     L   density    2rho_L     2rho^peak   /2rho^peak  /2rho_L  peak-lane  rescues barrier  fallback  degraded  private-cost  riders(late)")
res0 = {}
for L in list(range(8, 41, 2)) + [48, 56, 64]:
    d = run_c(L, 10 ** 6, 12, 0)
    Mf = 2 * L + 8
    flf = P.floor_table(Mf)
    rho = PK.rho_float(3, L, flf, Mf)
    rpk = PK.rho_peak_float(3, L, smax=min(70.0, 3 * L ** (1 / 3) + 12))["hi"]
    res0[L] = (d, rho, rpk)
    fbk = d["A"] + d["B"] + d["F"] + d["P"] + d["dfs"]
    check(d["stuck"] == 0 and d["verifyfail"] == 0,
          "L=%d: no stuck n; every 3 <= n <= 10^6 falls below itself within L steps (valid section of a member of P_L)" % L)
    print("    %3d  %.6f  %.6f  %.6f  %6.3f     %.4f  %s  %6d  %6d  %6d  %6d    %.6f  %d(%d)" % (
        L, d["density"], 2 * rho, 2 * rpk, d["density"] / (2 * rpk), d["density"] / (2 * rho),
        ("%.6f" % peak_lane[L]) if L in peak_lane else "   -    ", d["rescues"], d["barrier"], fbk, d["degraded"],
        d["ideal_cost_density"], d["riders"], d["late"]))
check(all(res0[L][0]["harmed"] == 0 for L in res0), "every rescued n is Collatz-bad (no harm to Collatz-good n) in all MODE 0 runs")
check(res0[40][0]["density"] / (2 * res0[40][2]) < 5 and res0[40][0]["density"] / (2 * res0[40][1]) < 0.1,
      "L=40: density < 5 x 2rho^peak and < 0.1 x 2rho_L (FINITE-EXACT on n <= 10^6)")
print("\n   interference in MODE 0: degraded = rescues whose actual best scale is below the private best scale; "
      "first deviation of the actual barrier path from the private one:")
for L in (16, 24, 32, 40, 48, 64):
    d = res0[L][0]
    print("     L=%d: degraded %d of %d rescues (%.1f%%): blocked (a desired flip point already frozen) %d, high flip %d, low flip %d; "
          "private needs a low rescue: %d" % (L, d["degraded"], d["rescues"], 100 * d["degraded"] / d["rescues"], d["blocked"], d["highflip"],
                                              d["lowflip"], d["ideal_needs_low"]))
check(all(res0[L][0]["blocked"] >= 0.95 * res0[L][0]["degraded"] for L in (16, 24, 32, 40)),
      "L = 16..40: at least 95% of the degraded rescues are blocked flips (the desired point lies on an earlier frozen certificate)")
print("\n   expensive (fallback A/F/B/P/dfs) rescues per undecided source, MODE 0 versus MODE 1 (default paths of Collatz-good n"
      " not frozen; diagnostic, NOT a valid member):")
print("     L   rho_L*N   MODE0: fallback  /(rho_L N)   degraded | MODE1: fallback  /(rho_L N)  degraded  density    verify-failures")
for L in (16, 24, 32, 40, 48, 64):
    d = run_c(L, 10 ** 6, 12, 1)
    d0, rho, rpk = res0[L]
    f0 = d0["A"] + d0["B"] + d0["F"] + d0["P"] + d0["dfs"]
    f1 = d["A"] + d["B"] + d["F"] + d["P"] + d["dfs"]
    print("    %3d  %7.0f   %6d        %.4f      %5d   |  %6d        %.4f     %5d    %.6f   %d" % (
        L, rho * 1e6, f0, f0 / (rho * 1e6), d0["degraded"], f1, f1 / (rho * 1e6), d["degraded"], d["density"], d["verifyfail"]))
    check(d["degraded"] < res0[L][0]["degraded"] and (f1 <= f0 or L < 24),
          "L=%d: not freezing Collatz-good paths reduces degraded rescues (%d < %d)%s" % (L, d["degraded"], res0[L][0]["degraded"],
                                                                                    "" if L < 24 else " and fallbacks (%d <= %d)" % (f1, f0)))
    check(d["verifyfail"] > 0 or L >= 64, "L=%d: MODE 1 breaks some earlier n (verify failures %d): freezing is needed without a domination argument" % (L, d["verifyfail"]))
print("\n   MODE 0 at N = 10^7:")
for L in (16, 32):
    d = run_c(L, 10 ** 7, 12, 0)
    check(d["stuck"] == 0 and d["verifyfail"] == 0 and d["harmed"] == 0,
          "L=%d N=10^7: valid (no stuck n, all n <= 10^7 descend within L, no harm); density %.6f (N=10^6: %.6f)" % (L, d["density"], res0[L][0]["density"]))
stamp("section 6")

# ------------------------------------------------------------------------------------------------------------------
header("7. One-flip domination (1FD): T^(j-1)((v-1)/2) <= T^j(v) for 1 <= j <= L")
for L in (16, 32, 64):
    rngd = random.Random(99 + L)
    nv = 0
    tot = 20000
    for _ in range(tot):
        v = rngd.getrandbits(400)
        v = v - (v % 12) + 11
        if P.one_fd_violation(v, L):
            nv += 1
    print("     L=%d: 1FD fails for %.3f of random v = 11 mod 12" % (L, nv / tot))
    check(0.15 < nv / tot < 0.5, "L=%d: one-flip domination fails for a constant fraction of barrier-type points (EMPIRICAL)" % L)
stamp("section 7")
print("\nALL CHECKS PASSED")
