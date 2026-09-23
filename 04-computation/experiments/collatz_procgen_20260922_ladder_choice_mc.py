#!/usr/bin/env python3
"""collatz_procgen_20260922_ladder_choice_mc.py -- does choice collapse the q = 5 exceptional set?

Uses collatz_procgen_20260922_ladder_eq_dfs.c (exact branch-and-bound per residue class; an independent
code path for lane one's DP program exceptional_general.c):
  1. exhaustive counts at m = 16, 20 for six games, compared with the DP program (must be identical);
  2. Monte Carlo over uniform random classes mod 2^m, m <= 64: every sampled class is decided EXACTLY
     (certificate of precision m exists or not); only the sampling is random (fixed seeds, splitmix64).
     Controls: the MC agrees with the exact DP fraction at m = 26 and, for the no-choice game, with the
     exact ballot fraction f_64(5) (collatz_procgen_20260922_ladder_ballot.py);
  3. first-moment indices: rho_q = min 2^(t-1)(1+q^-t) (no choice) and r_q = min 2^(t-1)/(1-q^-t)
     (full choice E_q), and the analytic positivity bound (1-K_q(t))/2 of the ladder note, Theorem 3.
"""
import importlib.util
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

import mpmath

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("ballot", os.path.join(HERE, "collatz_procgen_20260922_ladder_ballot.py"))
ballot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ballot)

MC = re.compile(r"mc m=(\d+) samples=(\d+) exceptional=(\d+) frac=(\S+) \+- (\S+)")
EX = re.compile(r"exhaustive m=(\d+) exceptional=(\d+)")
DP = re.compile(r"mode=(\d) m=\s*(\d+) exceptional=(\d+)")


def main():
    quick = "--quick" in sys.argv
    tmp = tempfile.mkdtemp(prefix="ladder_mc_")
    dfs = os.path.join(tmp, "eq_dfs")
    dp = os.path.join(tmp, "excgen")
    subprocess.run(["clang", "-O3", "-o", dfs, os.path.join(HERE, "collatz_procgen_20260922_ladder_eq_dfs.c"), "-lm"], check=True)
    subprocess.run(["clang", "-O3", "-o", dp, os.path.join(HERE, "collatz_procgen_20260922_exceptional_general.c"), "-lm"], check=True)
    fails = 0
    print("=" * 78)
    print("LADDER-MC. exact per-class search (eq_dfs.c) vs DP, then Monte Carlo to m=64")
    print("=" * 78)
    games = [("q=5 E_5", 5, 1, []), ("q=5 no choice", 5, 0, []), ("q=3 E", 3, 1, []),
             ("q=3 S={6 mod 8}", 3, 2, [3, 6]), ("q=3 S={2 mod 4}", 3, 2, [2, 2]), ("q=7 E_7", 7, 1, [])]
    for name, q, mode, extra in games:
        out = subprocess.run([dp, str(q), "1", "20", str(mode)] + [str(e) for e in extra],
                             capture_output=True, text=True, check=True).stdout
        dpv = {int(m): int(c) for _, m, c in DP.findall(out)}
        exv = {}
        for m in (16, 20):
            o = subprocess.run([dfs, str(q), str(mode), "exhaustive", str(m)] + [str(e) for e in extra],
                               capture_output=True, text=True, check=True).stdout
            exv[m] = int(EX.search(o).group(2))
        ok = all(exv[m] == dpv[m] for m in (16, 20))
        fails += not ok
        print(f"  {name:<18} exhaustive DFS m=16,20: {exv[16]:>7} {exv[20]:>8}   DP: {dpv[16]:>7} {dpv[20]:>8}   "
              f"{'IDENTICAL' if ok else 'MISMATCH'}")

    def mc(q, mode, m, n, seed):
        o = subprocess.run([dfs, str(q), str(mode), "mc", str(m), str(n), str(seed)],
                           capture_output=True, text=True, check=True).stdout
        g = MC.search(o)
        return int(g.group(3)), int(g.group(2)), float(g.group(4)), float(g.group(5))

    print()
    print("q=5: level-m exceptional Haar mass f_m (no choice: exact ballot; E_5: exact DP for m<=26, MC beyond)")
    N5 = ballot.exact_counts(5, 64)
    exact_E5 = {26: 4457200 / 2 ** 26}
    sched = [(26, 100000), (32, 100000), (40, 100000), (48, 100000), (56, 100000), (64, 100000)]
    if quick:
        sched = [(26, 20000), (32, 20000), (40, 20000), (48, 20000)]
    prev = None
    for m, n in sched:
        k, n_, f, se = mc(5, 1, m, n, 1000 + m)
        f0 = N5[m] / 2 ** m
        loc = ""
        if prev is not None and k > 0 and prev[1] > 0:
            loc = f"  local exponent log2(count)/m-slope: {1 + math.log2(f / prev[1]) / (m - prev[0]):.3f}"
        extra = f"  (exact DP {exact_E5[m]:.5f})" if m in exact_E5 else ""
        # Clopper-Pearson-style 95% upper bound for small counts (Poisson approx.)
        ub = (k + 1.96 * math.sqrt(k) + 2) / n_ if k < 50 else f + 1.96 * se
        print(f"  m={m:>2}: no-choice f_m={f0:.5f}   E_5 MC {k}/{n_} = {f:.5f} +- {se:.5f} (95% ub {ub:.5f}){extra}{loc}")
        if m in exact_E5 and abs(f - exact_E5[m]) > 4 * se:
            fails += 1
            print("    MC/DP disagreement at m=26")
        prev = (m, f)
    k, n_, f, se = mc(5, 0, 64, 20000 if quick else 100000, 777)
    f0 = N5[64] / 2 ** 64
    okc = abs(f - f0) < 4 * se
    fails += not okc
    print(f"  control: no-choice MC at m=64: {f:.5f} +- {se:.5f} vs exact ballot f_64(5) = {f0:.5f}: {'CONSISTENT' if okc else 'INCONSISTENT'}")
    print(f"  mu_5 = lim f_m (no choice) = 0.176025784562 (LADDER-C);  E_5 upper bound (rigorous, exact DP): "
          f"Haar(Bad(E_5)) <= f_26(E_5) = {exact_E5[26]:.5f}")

    print()
    print("q = 7, 9, 11, 13: E_q exceptional fraction (MC) vs no-choice f_m (exact) and mu_q")
    qsched = [32, 48] if quick else [32, 48, 64]
    for q in (7, 9, 11, 13):
        Nq = ballot.exact_counts(q, 64)
        row = []
        for m in qsched:
            k, n_, f, se = mc(q, 1, m, 5000 if quick else 20000, 2000 + 100 * q + m)
            row.append(f"m={m}: {f:.4f}+-{se:.4f} (no-choice {Nq[m]/2**m:.4f})")
        print(f"  q={q:>2}: " + "   ".join(row))

    print()
    print("first-moment indices (min over t on a grid; any fixed t gives a valid bound):")
    print("  rho_q = min_t 2^(t-1)(1+q^-t)   [no choice: < 1 iff drift > 0 iff q >= 5]")
    print("  r_q   = min_t 2^(t-1)/(1-q^-t)  [full choice: expected number of descending certificates per bit]")
    ts = [mpmath.mpf(i) / 2000 for i in range(1, 3000)]
    for q in (3, 5, 7, 9, 11, 13, 21, 31, 41, 101):
        rho = min(mpmath.power(2, t - 1) * (1 + mpmath.power(q, -t)) for t in ts)
        rq = min(mpmath.power(2, t - 1) / (1 - mpmath.power(q, -t)) for t in ts)
        t0 = mpmath.mpf("0.713")
        r0 = mpmath.power(2, t0 - 1) / (1 - mpmath.power(q, -t0))
        K = (mpmath.power(2, t0) * mpmath.power(q, -t0) / ((1 - mpmath.power(q, -2 * t0)) * (1 - r0))) if r0 < 1 else mpmath.inf
        lb = max(mpmath.mpf(0), (1 - K) / 2) if K != mpmath.inf else mpmath.mpf(0)
        print(f"  q={q:>3}: rho_q={mpmath.nstr(rho, 6):>8}  r_q={mpmath.nstr(rq, 6):>8}  K_q(0.713)={mpmath.nstr(K, 6):>8}  "
              f"=> Theorem 3 lower bound Haar(Bad(E_q)) >= {mpmath.nstr(lb, 5)}")
    shutil.rmtree(tmp)
    print(f"LADDER-MC TOTAL: {'ALL CHECKS PASS' if not fails else 'SOME CHECK FAILED'}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
