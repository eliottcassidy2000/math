# Auditor numerics-prereg: the census (4) of every_base.py with its two defects repaired, printed in the format of
# every_base.out, next to the stored rows.  The function below is the session's census_b with exactly two edits:
#   (a) the table is the full one (the caller passes q for h < 2^22 instead of q[: 2^21]);
#   (b) the skip "if X1 * X1 + coef * y1 * y1 > X: ... continue" is removed; the loop over h1 stops on the monotone
#       bound (1 + b (h1 // 2))^2 + coef y1^2 > X only (psi(h) >= h // 2, so X(h) >= 1 + b (h // 2) for every h).
# It is compared with my independent census (my_every_base.py).
import sys
import numpy as np
from numba import njit
from aud_common import psi_table
from my_every_base import census_b as my_census_b, PRIMES


@njit
def census_b_fixed(q, b, X, K, primes, cnt, div, digit2):
    H = q.shape[0]
    bad = 0
    for s in range(1, K + 1):
        coef = b << (s - 1)
        for h1 in range(H):
            y1 = 2 * h1 + 1
            X1 = 1 + b * (h1 - q[h1])
            if (1 + b * (h1 // 2)) * (1 + b * (h1 // 2)) + coef * y1 * y1 > X:
                break
            for h2 in range(h1, H):
                y2 = 2 * h2 + 1
                if (1 + b * (h2 // 2)) * X1 + coef * y1 * y2 > X:
                    break
                X2 = 1 + b * (h2 - q[h2])
                N = X1 * X2 + coef * y1 * y2
                if N > X:
                    continue
                w = 1 if h2 == h1 else 2
                if N % b != 1:
                    bad += 1
                digit2[((N - 1) // b) % b] += w
                for j in range(primes.shape[0]):
                    p = primes[j]
                    cnt[j, N % p] += w
                    if h2 != h1 and X1 % p == 0 and y1 % p == 0:
                        div[j, 0] += w
                        if N % p == 0:
                            div[j, 1] += w
    return bad


def main():
    K = 24
    X = 1 << K
    NB = 22
    psi = psi_table(1 << NB)
    h = np.arange((1 << NB) + 1, dtype=np.int64)
    q = (h - psi).astype(np.int64)[: 1 << NB]
    primes = np.array([3, 5, 7, 11, 13], dtype=np.int64)
    stored = {}
    for ln in open("C:/tmp/wt-s92/04-computation/experiments/s92_prime_powers_20261111/every_base.out"):
        if ln.startswith("    b = ") and "triples" in ln:
            stored[int(ln.split(":")[0].split("=")[1])] = ln.rstrip()
    print("(4) sheet s >= 1, N_b <= 2^24, repaired census (format of every_base.out); 'stored' = the row of every_base.out")
    worst = 0.0
    for b in (2, 3, 4, 5, 6, 10):
        cnt = np.zeros((len(primes), 16), np.int64)
        div = np.zeros((len(primes), 2), np.int64)
        digit2 = np.zeros(b, np.int64)
        bad = census_b_fixed(q, b, X, K, primes, cnt, div, digit2)
        tot = int(cnt[0].sum())
        mine = int(my_census_b(psi, b, X, 1, PRIMES)[4])
        cells = []
        for j, p in enumerate(primes.tolist()):
            f = p * cnt[j, :p] / tot
            if b % p == 0:
                cells.append("p=%d: N = 1 (mod p)" % p)
            else:
                fs = div[j, 0] / tot
                pred0, predr = 1 + (p - 1) * fs, 1 - fs
                rem = max(abs(f[0] - pred0), abs(f[1:] - predr).max())
                worst = max(worst, rem)
                cells.append("p=%d: %.3f [%.3f, %.3f] F p^2=%.2f rem=%.3f" % (p, f[0], f[1:].min(), f[1:].max(), fs * p * p, rem))
        d2 = b * digit2 / digit2.sum()
        print("    b = %2d: %9d triples (my independent census: %d, %s); N = 1 (mod b): %s; second digit flat to %.3f; %s" % (
            b, tot, mine, "equal" if tot == mine else "DIFFERENT", bad == 0, abs(d2 - 1).max(), "; ".join(cells)))
        print("      stored: " + stored[b].strip())
    print("largest remainder over the six bases and the primes not dividing b (repaired census, session's definition): %.4f" % worst)


if __name__ == "__main__":
    main()
