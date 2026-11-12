# v_sumT.py - measurement of sum_eps T'(eps) of Prop. 8 (the quantity the lane left untested) in balanced cells
# t' = ceil(w/2), s = floor(w/2) (the regime w ~ n/2 of hypothesis (I_g)), w = 6..12.  Exact integer counts.
# T'(eps) = #{(m, m') in ({0,1}^Z)^2 : sum_j (m_j - m'_j) lambda_j in Lam+},
# Lam+ = {dc X - dP Xi + th : |dc| < 2^s, |dP| < 2^(t'-1), |th| <= 2}.
import sys, math, time
import numpy as np
from numba import njit

@njit(cache=True)
def kap(d):
    if d == 1:
        return 2
    return 1 << (d - 1)

@njit(cache=True)
def sumT(tp, s, w):
    eps = np.zeros(w, dtype=np.int64)
    lam = np.zeros(w, dtype=np.int64)
    zpos = np.zeros(w, dtype=np.int64)
    dcmax = (1 << s) - 1
    dPmax = (1 << (tp - 1)) - 1
    Bmax = (dcmax + dPmax) * (1 << w) + 8
    tab = np.zeros(2 * Bmax + 1, dtype=np.uint8)
    total = 1
    for i in range(w):
        total *= 3
    sT = 0; sdiag = 0; soff_small = 0; npat = 0
    for code in range(total):
        c = code
        for i in range(w):
            eps[i] = (c % 3) - 1
            c //= 3
        if eps[0] == 0 or eps[w - 1] == 0:
            continue
        npat += 1
        z = 0; X = 0; Xi = 0
        for i in range(w):
            X += eps[i] * (1 << i)
            Xi += eps[i] * (1 << (w - 1 - i))
            if eps[i] == 0:
                zpos[z] = i; z += 1
        for a in range(z):
            j = zpos[a]
            t = 0
            for i in range(w):
                if eps[i] != 0:
                    t += kap(abs(i - j)) * eps[i]
            lam[a] = t
        for dc in range(-dcmax, dcmax + 1):
            for dP in range(-dPmax, dPmax + 1):
                tab[dc * X - dP * Xi + Bmax] = 1
        T = 0
        tot3 = 1
        for a in range(z):
            tot3 *= 3
        for dcode in range(tot3):
            c2 = dcode
            v = 0; wt = 1
            for a in range(z):
                d = (c2 % 3) - 1
                c2 //= 3
                if d == 0:
                    wt *= 2
                else:
                    v += d * lam[a]
            hit = False
            for th in range(-2, 3):
                u = v - th
                if -Bmax <= u <= Bmax and tab[u + Bmax]:
                    hit = True
                    break
            if hit:
                T += wt
        for dc in range(-dcmax, dcmax + 1):
            for dP in range(-dPmax, dPmax + 1):
                tab[dc * X - dP * Xi + Bmax] = 0
        sT += T
        sdiag += 1 << z
    return sT, sdiag, npat

def main():
    w0 = int(sys.argv[1]); w1 = int(sys.argv[2])
    prev = None
    for w in range(w0, w1 + 1):
        tp = (w + 1) // 2; s = w // 2
        t0 = time.time()
        sT, sdiag, npat = sumT(tp, s, w)
        x = sT ** (1.0 / w)
        print("w=%2d cell (t'=%d, s=%d): sum_eps T' = %d  diagonal = %d (4^(w-1))  off-diagonal = %d  (sum T')^(1/w) = %.4f  ratio to previous w = %s  sum T'/(16/3)^w = %.4e  sum T'/(6/sqrt2)^w = %.4f  CS bound/4^w = (3^w sumT')^(1/2)/4^w = %.4f   (%.1fs)" % (
            w, tp, s, sT, sdiag, sT - sdiag, x, ("%.3f" % (sT / prev)) if prev else "-", sT / (16 / 3) ** w, sT / (6 / math.sqrt(2)) ** w, math.sqrt(3.0 ** w * sT) / 4.0 ** w, time.time() - t0))
        prev = sT
        sys.stdout.flush()

if __name__ == "__main__":
    main()
