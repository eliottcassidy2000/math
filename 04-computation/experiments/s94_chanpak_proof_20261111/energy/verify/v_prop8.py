# v_prop8.py - brute-force check of Prop. 8 (|M(eps)|^2 <= T'(eps)) and measurement of sum_eps T'(eps);
# also an independent recount of "odd window pairs admitting a solution in the cell" (lane column distinct_window_pairs).
# Usage: python v_prop8.py tp s w [tp s w ...]
import sys, math, time
import numpy as np
from numba import njit

@njit(cache=True)
def eta_table(n):
    N = 1 << n
    q = np.zeros(N, dtype=np.int64)
    e = np.zeros(N, dtype=np.int64)
    for h in range(1, N):
        p = h >> 1
        if h & 1:
            q[h] = q[p] + (p + 1) // 2
        else:
            q[h] = q[p]
    for h in range(N):
        e[h] = 2 * h + 1 - 2 * q[h]
    return e

@njit(cache=True)
def admissible(e, tp, s, w):
    # adm[x, x'] = 1 iff some P in [2^(tp-1), 2^tp), Z in [0, 2^s) has eta(P|x|Z) = eta(P|x'|Z); odd pairs only
    W = 1 << w
    adm = np.zeros((W, W), dtype=np.uint8)
    for x in range(W):
        for xp in range(W):
            if ((x ^ xp) & 1) == 0 or (((x ^ xp) >> (w - 1)) & 1) == 0:
                continue
            found = False
            for P in range(1 << (tp - 1), 1 << tp):
                base = P << (w + s)
                for Z in range(1 << s):
                    if e[base | (x << s) | Z] == e[base | (xp << s) | Z]:
                        found = True
                        break
                if found:
                    break
            if found:
                adm[x, xp] = 1
    return adm

@njit(cache=True)
def kap(d):
    if d == 1:
        return 2
    return 1 << (d - 1)

@njit(cache=True)
def prop8(adm, tp, s, w):
    # loop over odd patterns eps in {-1,0,1}^w; returns sumM, sumT, sum sqrt T, violations, npatterns
    eps = np.zeros(w, dtype=np.int64)
    lam = np.zeros(w, dtype=np.int64)
    zpos = np.zeros(w, dtype=np.int64)
    sumM = 0; sumT = 0; sumsq = 0.0; viol = 0; npat = 0; diag = 0
    total = 1
    for i in range(w):
        total *= 3
    dcmax = (1 << s) - 1
    dPmax = (1 << (tp - 1)) - 1
    for code in range(total):
        c = code
        for i in range(w):
            eps[i] = (c % 3) - 1
            c //= 3
        if eps[0] == 0 or eps[w - 1] == 0:
            continue
        npat += 1
        z = 0
        X = 0; Xi = 0; xbase = 0; xpbase = 0
        for i in range(w):
            X += eps[i] * (1 << i)
            Xi += eps[i] * (1 << (w - 1 - i))
            if eps[i] == 1:
                xbase |= 1 << i
            elif eps[i] == -1:
                xpbase |= 1 << i
            else:
                zpos[z] = i; z += 1
        for a in range(z):
            j = zpos[a]
            t = 0
            for i in range(w):
                if eps[i] != 0:
                    t += kap(abs(i - j)) * eps[i]
            lam[a] = t
        # |M(eps)|
        M = 0
        for mm in range(1 << z):
            com = 0
            for a in range(z):
                if (mm >> a) & 1:
                    com |= 1 << zpos[a]
            if adm[xbase | com, xpbase | com]:
                M += 1
        # Lam+ table
        B = dcmax * abs(X) + dPmax * abs(Xi) + 2
        tab = np.zeros(2 * B + 1, dtype=np.uint8)
        for dc in range(-dcmax, dcmax + 1):
            for dP in range(-dPmax, dPmax + 1):
                v = dc * X - dP * Xi
                for th in range(-2, 3):
                    tab[v + th + B] = 1
        # T'(eps) by enumeration of delta in {-1,0,1}^z with weight 2^(#zeros)
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
            if -B <= v <= B and tab[v + B]:
                T += wt
        sumM += M; sumT += T; sumsq += math.sqrt(T)
        diag += 1 << z
        if M * M > T:
            viol += 1
    return sumM, sumT, sumsq, viol, npat, diag

def main():
    args = [int(a) for a in sys.argv[1:]]
    for k in range(0, len(args), 3):
        tp, s, w = args[k:k + 3]
        n = tp + s + w
        t0 = time.time()
        e = eta_table(n)
        adm = admissible(e, tp, s, w)
        nadm = int(adm.sum())
        sumM, sumT, sumsq, viol, npat, diag = prop8(adm, tp, s, w)
        print("cell (t'=%d, s=%d), w=%d, n=%d: ordered odd window pairs admitting a solution = %d (unordered %d) of %d odd pairs ; sum_eps |M(eps)| = %d (equal: %s)" % (
            tp, s, w, n, nadm, nadm // 2, 4 ** (w - 1), sumM, sumM == nadm))
        print("   Prop. 8: patterns=%d violations of |M|^2 <= T' = %d ; sum_eps T'(eps) = %d ; diagonal sum 2^z = %d (4^(w-1) = %d)" % (npat, viol, sumT, diag, 4 ** (w - 1)))
        print("   sum T' as x^w: x = %.4f  (4, 6/sqrt2 = 4.243, 16/3 = 5.333, 6) ; bound (3^w sum T')^(1/2) = %.3e vs 4^w = %.3e vs actual %d ; sum sqrt(T') = %.3e" % (
            sumT ** (1.0 / w), math.sqrt(3.0 ** w * sumT), 4.0 ** w, nadm, sumsq))
        print("   (%.1fs)" % (time.time() - t0))
        sys.stdout.flush()

if __name__ == "__main__":
    main()
