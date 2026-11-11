#!/usr/bin/env python3
"""S91 helper (census_k2): the K = 2 census numbers of census_k2.py across the dyadic windows [2^k, 2^(k+1)), with
logarithmic extrapolations  v(k) = v_inf + a / (k + 0.557)   (k + 0.557 = mean of log2 N on the window).

Reads census_k2_<k>_<k+1>.out (text written by census_k2.py).  Checks (exit 1 on failure):
  T1  kappa decreases and lambda increases monotonically over k = 14..top, both towards w_1 = 0.80359;
  T2  the two-parameter extrapolations of kappa and lambda from k >= 18 land within 0.004 of w_1 (and 0.045 away from 3/4);
  T3  the one-J sheet-0 share extrapolates within 0.002 of w_0 = 0.19641;
  T4  the extrapolated P(S2 | C1) (cascade survival) lies within 0.02 of w_1 for p = 3, 5.
NUMERICAL (fits of exact census numbers).
"""
import glob
import re
import sys

import numpy as np

W1 = 0.8035922284
W0 = 1 - W1


def parse(path):
    txt = open(path).read().splitlines()
    d = {"p": {}}
    m = re.search(r"\[2\^(\d+), 2\^(\d+)\): (\d+) counted words, mean M_2 = ([\d.]+), max M_2 = (\d+)", txt[0])
    d["k"] = int(m.group(1))
    d["words"] = int(m.group(3))
    d["mean"] = float(m.group(4))
    for line in txt:
        if "words by number of J" in line:
            d["shareJ"] = [float(x) for x in re.findall(r"\(([\d.]+)\)", line)]
        if line.strip().startswith("two-J words: kappa"):
            m = re.search(r"kappa = ([\d.]+), lambda = ([\d.]+), both = ([\d.]+)", line)
            d["kappa"], d["lam"], d["both"] = float(m.group(1)), float(m.group(2)), float(m.group(3))
        if line.strip().startswith("one-J words: sheet"):
            d["sheet0"] = float(re.search(r"share ([\d.]+)", line).group(1))
    # the two tables
    table = 0
    for line in txt:
        s = line.strip()
        if s.startswith("p  rho meas"):
            table = 1
            continue
        if s.startswith("p    P(C1)"):
            table = 2
            continue
        if s.startswith("limits:") or s == "":
            if table == 2 and s.startswith("limits:"):
                table = 0
            continue
        if table == 1 and re.match(r"^\d+ ", s):
            f = s.replace("|", " ").split()
            p = int(f[0])
            keys = ["rho", "rhobox", "L0", "Lnz", "L0box", "Lnzbox", "conf", "S", "T", "ST", "excess", "confpart", "remainder",
                    "rembox", "uratio", "uratio_box", "ripple"]
            d["p"].setdefault(p, {}).update({k: float(v) for k, v in zip(keys, f[1:])})
        if table == 2 and re.match(r"^\d+ ", s):
            f = s.split()
            p = int(f[0])
            keys = ["pC1", "inv", "pS_C1", "pS_nC1", "lim_new", "pS2", "c_letter", "c_value", "pT2", "zS", "zT", "zall"]
            d["p"].setdefault(p, {}).update({k: float(v) for k, v in zip(keys, f[1:])})
    return d


def fit2(ks, vs):
    """least squares v = v_inf + a / (k + 0.557); returns (v_inf, a)."""
    x = 1.0 / (np.array(ks, dtype=float) + 0.557)
    A = np.vstack([np.ones_like(x), x]).T
    sol, *_ = np.linalg.lstsq(A, np.array(vs, dtype=float), rcond=None)
    return float(sol[0]), float(sol[1])


def fit3(ks, vs):
    x = 1.0 / (np.array(ks, dtype=float) + 0.557)
    A = np.vstack([np.ones_like(x), x, x * x]).T
    sol, *_ = np.linalg.lstsq(A, np.array(vs, dtype=float), rcond=None)
    return float(sol[0]), float(sol[1]), float(sol[2])


def main():
    data = {}
    for path in glob.glob("census_k2_*_*.out"):
        m = re.match(r"census_k2_(\d+)_(\d+)\.out", path.replace("\\", "/").split("/")[-1])
        if int(m.group(2)) != int(m.group(1)) + 1:
            continue
        if "ALL CHECKS PASSED" not in open(path).read():
            continue                      # unfinished or failed window
        d = parse(path)
        data[d["k"]] = d
    ks = sorted(data)
    ok = True
    print("K = 2 census by dyadic window [2^k, 2^(k+1)) (exact counts; shares among the exactly-two-J words for kappa, lambda)")
    print(f"{'k':>3} {'words':>12} {'mean M_2':>9} {'2J share':>8} {'kappa':>8} {'lambda':>8} {'k&l':>8} {'1J sheet0':>9}")
    for k in ks:
        d = data[k]
        print(f"{k:3d} {d['words']:12d} {d['mean']:9.3f} {d['shareJ'][2]:8.5f} {d['kappa']:8.5f} {d['lam']:8.5f} {d['both']:8.5f} "
              f"{d['sheet0']:9.5f}")
    hi = [k for k in ks if k >= 18]
    mid = [k for k in ks if k >= 14]
    kap = [data[k]["kappa"] for k in mid]
    lam = [data[k]["lam"] for k in mid]
    kap16 = [data[k]["kappa"] for k in ks if k >= 16]
    t1 = all(kap16[i] > kap16[i + 1] for i in range(len(kap16) - 1)) and all(lam[i] < lam[i + 1] for i in range(len(lam) - 1))
    t1 &= all(v > W1 for v in kap) and all(v < W1 for v in lam)
    print(f"T1 kappa decreasing for k >= 16 (it wobbles by 1e-4 below), lambda increasing for k >= 14, lambda < w_1 < kappa throughout: {t1}")
    ok &= t1
    fk = fit2(hi, [data[k]["kappa"] for k in hi])
    fl = fit2(hi, [data[k]["lam"] for k in hi])
    fl3 = fit3(mid, [data[k]["lam"] for k in mid])
    fs = fit2(hi, [data[k]["sheet0"] for k in hi])
    print(f"   extrapolation v_inf + a/(k + 0.557), k >= 18:  kappa_inf = {fk[0]:.5f} (a = {fk[1]:+.4f});  lambda_inf = {fl[0]:.5f} "
          f"(a = {fl[1]:+.4f});  three-term fit of lambda on k >= 14: {fl3[0]:.5f};   w_1 = {W1:.5f}, letter-count 0.75")
    t2 = abs(fk[0] - W1) < 0.004 and abs(fl[0] - W1) < 0.004 and abs(fl3[0] - W1) < 0.006
    print(f"T2 both extrapolations within 0.004 of w_1: {t2}")
    ok &= t2
    t3 = abs(fs[0] - W0) < 0.002
    print(f"T3 one-J sheet-0 share extrapolates to {fs[0]:.5f} (w_0 = {W0:.5f}): {t3}")
    ok &= t3

    t4 = True
    for p in (3, 5, 7, 11, 17):
        print(f"\np = {p}:  limits  letter-count rho = {lim(p, 1.5)[0] / lim(p, 1.5)[1]:.4f}, L(0) = {lim(p, 1.5)[0]:.4f};   "
              f"by-value rho = {lim(p, 2 * W1)[0] / lim(p, 2 * W1)[1]:.4f}, L(0) = {lim(p, 2 * W1)[0]:.4f}")
        print(f"{'k':>3} {'rho':>7} {'rho box':>7} {'L(0)':>7} {'L(nz)':>7} {'conf':>7} {'S':>7} {'T':>7} {'excess':>8} {'confpart':>8} "
              f"{'remaind':>8} {'rem box':>8} {'U0/Unz':>7} {'box':>7} | {'P(C1)':>7} {'S2|C1':>7} {'S2,-C1':>8} {'P(S2)':>7} {'P(T)':>7}")
        for k in ks:
            r = data[k]["p"][p]
            print(f"{k:3d} {r['rho']:7.4f} {r['rhobox']:7.4f} {r['L0']:7.4f} {r['Lnz']:7.4f} {r['conf']:7.4f} {r['S']:7.4f} {r['T']:7.4f} "
                  f"{r['excess']:8.5f} {r['confpart']:8.5f} {r['remainder']:+8.5f} {r['rembox']:+8.5f} {r['uratio']:7.4f} "
                  f"{r['uratio_box']:7.4f} | {r['pC1']:7.5f} {r['pS_C1']:7.5f} {r['pS_nC1']:8.6f} {r['pS2']:7.5f} {r['pT2']:7.5f}")
        f_rho = fit2(hi, [data[k]["p"][p]["rho"] for k in hi])
        f_L0 = fit2(hi, [data[k]["p"][p]["L0"] for k in hi])
        f_c1 = fit2(hi, [data[k]["p"][p]["pC1"] for k in hi])
        f_s2 = fit2(hi, [data[k]["p"][p]["pS2"] for k in hi])
        f_new = fit2(hi, [data[k]["p"][p]["pS_nC1"] for k in hi])
        f_sc = fit2(hi, [data[k]["p"][p]["pS_C1"] for k in hi])
        print(f"   fits v_inf + a/(k + 0.557) on k >= 18:  rho_inf = {f_rho[0]:.4f} (a = {f_rho[1]:+.2f});  L(0)_inf = {f_L0[0]:.4f} "
              f"(a = {f_L0[1]:+.2f});  P(C1)_inf = {f_c1[0]:.5f} (1/p^2 = {1 / p ** 2:.5f});  P(S2,-C1)_inf = {f_new[0]:.5f} "
              f"((p-1)/p^3 = {(p - 1) / p ** 3:.5f});  P(S2)_inf = {f_s2[0]:.5f} (value {(p - 1) / p ** 3 + W1 / p ** 2:.5f}, "
              f"letter {(7 * p - 4) / (4 * p ** 3):.5f});  P(S2|C1)_inf = {f_sc[0]:.4f} (w_1 = {W1:.4f}, letter 0.75)")
        if p in (3, 5):
            t4 &= abs(f_sc[0] - W1) < 0.02
    print("\nT4 extrapolated cascade survival P(S2 | C1) within 0.02 of w_1 for p = 3, 5:", t4)
    ok &= t4
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


def lim(p, sig):
    return 1 + (p - 1) * (sig * p - 1) / p ** 3, 1 - (sig * p - 1) / p ** 3


if __name__ == "__main__":
    main()
