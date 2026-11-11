#!/usr/bin/env python3
"""S91 helper (census_k2), step 3: the two-root-change census M_2 of HYP-9506 by value, with every counted word
classified as confined or not (rc_census.py), and the first-order mechanism tested at K = 2.

Usage:  python3 census_k2.py LOG_LO LOG_HI [nthreads]       window = odd N in [2^LOG_LO, 2^LOG_HI)
Output: text on stdout and the counters in census_k2_<LO>_<HI>.npz (read by k2_trend.py).

For each prime p (3, 5, 7, 11, 13, 17, 23, 31):
  rho_p = mean(M_2 | p | N) / mean of the class means of the nonzero classes     (the statistic of eval_s15.py)
  confined = S or T (see rc_census.py);  a confined word lies in its distinguished class c_g exactly.
  groups g = (number of J, shift type st);  c_g = 0 for st = 0, 1 and c_g = -E != 0 for st = 2.
  FIRST-ORDER LAW (box identity of S91, CONTEXT section 2 item 3, applied at the last root change with the MEASURED
  confined shares of the same window; it predicts the p class values from the group shares):
        P(N = r | g) = (1 - D_g)/p + D_g [r = c_g],     D_g = pi_S|g + pi_T|g - 1/p^2
  (for the J-free words the first block is the single state (1, 1): D_g = (pi_T p^2 - 1)/(p^2 - 1)).
  Hence   P(0) - mean_{r != 0} P(r) = D^0 - D^x/(p - 1),   D^0 = sum over groups with c_g = 0 of P(g) D_g, D^x the rest;
  the unconfined words alone:  U(0) - mean U(r != 0) = - sum_g eps_g P(g) (1/p^2 - pi_ST|g),  eps = 1 or -1/(p - 1).
Cascade quantities of the two-J words A J B J C (the by-value analogues of THM-4916's letter-count numbers):
  kappa = share with k1 = 1 (A ends in O or B starts with Z),  lambda = share with l2 = 1 (B ends in O or C starts with Z)
  letter-count: 3/4, 3/4;  by-value limit claimed here: w_1 = 0.80359 for both.
"""
import math
import sys
import time

import numpy as np

import rc_census as rc

PRIMES = [3, 5, 7, 11, 13, 17, 23, 31]
W1 = 0.8035922284
W0 = 1 - W1
# rc_lf.c output lf_2p24.txt (kind-pasteur S15), window [2^23, 2^24), K = 2
REF_2P24 = {"mean": 264.6185, "maxM": 3792,
            5: [380.3774, 233.8882, 234.2145, 233.6335, 240.9789],
            7: [317.1457, 253.0765, 259.7396, 253.4783, 256.5729, 250.5001, 261.8165]}
REF_RHO = {5: 1.6140, 17: 1.5338, 11: 1.3739, 7: 1.2395, 13: 1.1710, 23: 1.1240, 31: 1.1020, 3: 1.8565}


def limits(p):
    """(L(0), L(nonzero)) in the letter-count limit (THM-4916) and in the by-value limit claimed by this helper."""
    letter = (1 + (3 * p - 2) * (p - 1) / (2 * p ** 3), 1 - (3 * p - 2) / (2 * p ** 3))
    sig = 2 * W1
    value = (1 + (p - 1) * (sig * p - 1) / p ** 3, 1 - (sig * p - 1) / p ** 3)
    return letter, value


def analyse(res, lo, hi, out):
    A, B = 1 << lo, 1 << hi
    ok = True
    tot = res["tot"]
    nodd = (B - A) // 2
    words = res["nwords"]
    allcnt = res["allcnt"]           # [nJ, parity, z]
    allp = res["allp"]               # [prime, nJ, st]
    cnt = res["cnt"]                 # [prime, nJ, S, T, C1, st]
    cnt2 = res["cnt2"]               # [k1, l2]
    byJ = allcnt.sum(axis=(1, 2))
    mean = words / nodd
    out(f"K = 2 census, odd N in [2^{lo}, 2^{hi}): {words} counted words, mean M_2 = {mean:.4f}, max M_2 = {res['maxM']}, "
        f"{res['ntasks']} tasks")
    out(f"  words by number of J: 0: {int(byJ[0])} ({byJ[0] / words:.5f}), 1: {int(byJ[1])} ({byJ[1] / words:.5f}), "
        f"2: {int(byJ[2])} ({byJ[2] / words:.5f})")
    out(f"  confinement lemma (every word with S or T has value = -[c = 0] E mod p): violations {res['lemma_violations']};  "
        f"16-bit counters exact: {res['overflow_ok']}")
    ok &= res["lemma_violations"] == 0 and res["overflow_ok"] and int(allcnt.sum()) == words
    two = int(cnt2.sum())
    ok &= two == int(byJ[2])
    kappa = cnt2[1, :].sum() / two
    lam = cnt2[:, 1].sum() / two
    out(f"  two-J words: kappa = {kappa:.5f}, lambda = {lam:.5f}, both = {cnt2[1, 1] / two:.5f} (product {kappa * lam:.5f});  "
        f"letter-count 0.75, by-value limit {W1:.5f}")
    one = int(byJ[1])
    sheet0_1 = int(allcnt[1, 0, 0])
    out(f"  one-J words: sheet s = 0 share {sheet0_1 / one:.5f}  (by-value limit {W0:.5f}, letter-count 0.25)")
    rows = {}
    out("")
    out(f"  {'p':>3} {'rho meas':>9} {'rho box':>8} {'L(0)':>7} {'L(nz)':>7} {'L0 box':>7} {'Lnz box':>7} | {'conf':>7} {'S':>7} {'T':>7} "
        f"{'S&T':>8} | {'excess':>8} {'conf part':>9} {'remainder':>9} {'rem box':>8} | {'U0/Unz':>7} {'box':>7} {'ripple':>7}")
    for j, p in enumerate(PRIMES):
        means, n = rc.class_means(tot, A, B, p)
        L = means / mean
        Lnz = float(np.mean(L[1:]))
        rho = L[0] / Lnz
        conf = res["confhist"][j, :p].astype(np.float64)
        ncf = conf.sum()
        cmeans = conf / n
        umeans = means - cmeans
        u0 = umeans[0]
        unz = float(np.mean(umeans[1:]))
        ripple = float(np.max(np.abs(L[1:] - Lnz)) / Lnz)
        # shares
        nS = cnt[j, :, 1].sum()
        nT = cnt[j, :, :, 1].sum()
        nST = cnt[j, :, 1, 1].sum()
        ok &= int(nS + nT - nST) == int(ncf)
        # box prediction from the group shares
        P0 = 0.0
        Pc0 = 0.0          # confined mass with c_g = 0
        Pcx = 0.0
        rem = 0.0          # predicted U(0) - mean U(nz)
        unc_tot = 0.0
        unc0 = 0.0
        for nJ in range(3):
            for st in range(3):
                g = allp[j, nJ, st]
                if g == 0:
                    continue
                Pg = g / words
                piS = cnt[j, nJ, 1, :, :, st].sum() / g
                piT = cnt[j, nJ, :, 1, :, st].sum() / g
                piST = cnt[j, nJ, 1, 1, :, st].sum() / g
                cf = piS + piT - piST
                if nJ == 0:
                    D = (piT * p * p - 1) / (p * p - 1)
                else:
                    D = piS + piT - 1 / p ** 2
                zero = st != 2
                P0 += Pg * ((1 - D) / p + (D if zero else 0.0))
                if zero:
                    Pc0 += Pg * cf
                    u_dist = (1 - D) / p + D - cf        # unconfined mass in the distinguished class 0
                    unc0 += Pg * u_dist
                else:
                    Pcx += Pg * cf
                    unc0 += Pg * (1 - D) / p
                unc_tot += Pg * (1 - cf)
        L0box = p * P0
        Lnzbox = p * (1 - P0) / (p - 1)
        rhobox = L0box / Lnzbox
        uratio_box = unc0 / ((unc_tot - unc0) / (p - 1))
        excess = (L[0] - Lnz) / p
        # measured confined part: class-0 confined share minus the mean confined share of a nonzero class
        cshare = conf / ncf * (ncf / words)
        confpart = cshare[0] - float(np.mean(cshare[1:]))
        remainder = excess - confpart
        rembox = (unc0 - (unc_tot - unc0) / (p - 1))
        out(f"  {p:3d} {rho:9.4f} {rhobox:8.4f} {L[0]:7.4f} {Lnz:7.4f} {L0box:7.4f} {Lnzbox:7.4f} | {ncf / words:7.4f} {nS / words:7.4f} "
            f"{nT / words:7.4f} {nST / words:8.5f} | {excess:8.5f} {confpart:9.5f} {remainder:+9.5f} {rembox:+8.5f} | "
            f"{u0 / unz:7.4f} {uratio_box:7.4f} {ripple:7.4f}")
        rows[p] = dict(rho=rho, rhobox=rhobox, L0=L[0], Lnz=Lnz, conf=ncf / words, S=nS / words, T=nT / words, excess=excess,
                       confpart=confpart, remainder=remainder, rembox=rembox, uratio=u0 / unz, uratio_box=uratio_box,
                       L0box=L0box, Lnzbox=Lnzbox, classes=L.copy())
    out("")
    out("  cascade of the two-J words (shares among the two-J words of the window):")
    out(f"  {'p':>3} {'P(C1)':>8} {'1/p^2':>8} {'P(S2|C1)':>9} {'P(S2,-C1)':>10} {'(p-1)/p^3':>10} {'P(S2)':>8} {'letter c':>9} {'value c':>9} "
        f"{'P(T)':>8} {'P(c_g=0|S2)':>12} {'P(c_g=0|T)':>11} {'P(c_g=0)':>9}")
    for j, p in enumerate(PRIMES):
        c2 = cnt[j, 2]               # [S, T, C1, st]
        g2 = allp[j, 2].sum()
        pC1 = c2[:, :, 1].sum() / g2
        pS = c2[1].sum() / g2
        pS_C1 = c2[1, :, 1].sum() / max(c2[:, :, 1].sum(), 1)
        pS_nC1 = c2[1, :, 0].sum() / g2
        pT = c2[:, 1].sum() / g2
        z_S = c2[1, :, :, :2].sum() / max(c2[1].sum(), 1)
        z_T = c2[:, 1, :, :2].sum() / max(c2[:, 1].sum(), 1)
        z_all = allp[j, 2, :2].sum() / g2
        out(f"  {p:3d} {pC1:8.5f} {1 / p ** 2:8.5f} {pS_C1:9.5f} {pS_nC1:10.6f} {(p - 1) / p ** 3:10.6f} {pS:8.5f} "
            f"{(7 * p - 4) / (4 * p ** 3):9.5f} {(p - 1) / p ** 3 + W1 / p ** 2:9.5f} {pT:8.5f} {z_S:12.5f} {z_T:11.5f} {z_all:9.5f}")
        rows[p].update(pC1=pC1, pS2=pS, pS_C1=pS_C1, pS_nC1=pS_nC1, pT2=pT, zS=z_S, zT=z_T, zall=z_all)
    out("")
    out("  limits:   p   letter-count L(0), L(nz), rho        by-value L(0), L(nz), rho        measured rho")
    for p in PRIMES:
        (a0, a1), (b0, b1) = limits(p)
        out(f"          {p:3d}   {a0:.4f} {a1:.4f} {a0 / a1:.4f}               {b0:.4f} {b1:.4f} {b0 / b1:.4f}               {rows[p]['rho']:.4f}")
    return ok, rows, dict(kappa=kappa, lam=lam, mean=mean, words=words, byJ=byJ, sheet0_1=sheet0_1 / one)


def main():
    lo = int(sys.argv[1])
    hi = int(sys.argv[2])
    nthreads = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    lines = []

    def out(s):
        print(s, flush=True)
        lines.append(s)

    t0 = time.time()
    res = rc.census(2, 1 << lo, 1 << hi, PRIMES, nthreads=nthreads)
    t1 = time.time()
    ok, rows, glob = analyse(res, lo, hi, out)
    out(f"  census time {t1 - t0:.1f} s on {nthreads} threads")
    if (lo, hi) == (23, 24):
        tot = res["tot"]
        rep = abs(glob["mean"] - REF_2P24["mean"]) < 5e-5 and res["maxM"] == REF_2P24["maxM"]
        for p in (5, 7):
            means, n = rc.class_means(tot, 1 << lo, 1 << hi, p)
            rep &= bool(np.all(np.abs(means - np.array(REF_2P24[p])) < 5e-5))
        for p, r in REF_RHO.items():
            rep &= abs(rows[p]["rho"] - r) < 5e-5
        out(f"  reproduction of the packet (lf_2p24.txt: overall mean 264.6185, max 3792, class means mod 5 and mod 7 to 4 decimals; "
            f"rho_p of HYP-9506 (B) and of census_mechanism_2p24.txt for p = 3): {rep}")
        ok &= rep
    np.savez_compressed(f"census_k2_{lo}_{hi}.npz", confhist=res["confhist"], cnt=res["cnt"], allcnt=res["allcnt"],
                        allp=res["allp"], cnt2=res["cnt2"], nwords=res["nwords"], maxM=res["maxM"],
                        primes=np.array(PRIMES),
                        classes=np.array([np.pad(rows[p]["classes"], (0, 31 - p)) for p in PRIMES]))
    out("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
