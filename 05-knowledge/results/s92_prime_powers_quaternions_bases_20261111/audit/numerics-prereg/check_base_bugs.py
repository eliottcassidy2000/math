# Auditor numerics-prereg: localisation of two enumeration defects in the base-b censuses of the session.
#   every_base.py, census_b (item (4) of every_base.out):  (i) the table q is cut at h < 2^21;
#       (ii) a state h1 is skipped when N(h1, h1) > X, although N(h1, h2) < N(h1, h1) is possible for h2 > h1
#       (X(h) = 1 + b psi(h) is not monotone in h).
#   bridges_exact.py, S10 (base_head_census): the table psi has only 2^12 entries, so the census is the box
#       h1 <= h2 < 4096, not "N_beta <= 2^22".
# Everything below calls the session's OWN functions (imported from the copies in session_copies/) with their own
# arguments and with corrected arguments, and compares with (a) a pure-python brute force, (b) my census.
import sys
import time
import numpy as np

sys.path.insert(0, "session_copies")
import every_base as EB          # noqa: E402
import bridges_exact as BE       # noqa: E402
from aud_common import psi_table  # noqa: E402
from my_every_base import census_b as my_census_b, PRIMES  # noqa: E402


def brute_python(psi, b, X, sheet):
    """ordered triples with N_b <= X on one sheet; trivial bounds only."""
    tot = 0
    s = 0 if sheet == 0 else 1
    while True:
        e = 1 if s == 0 else 0
        if 1 + (b * ((1 << s) - e)) // 2 > X:
            break
        h1 = 0
        while True:
            y1 = 2 * h1 + 1
            if (b * ((y1 << s) - e)) // 2 > X:
                break
            X1 = 1 + b * int(psi[h1])
            h2 = 0
            while True:
                y2 = 2 * h2 + 1
                c3 = ((y1 << s) * y2 - e) // 2
                if b * c3 > X:
                    break
                if X1 * (1 + b * int(psi[h2])) + b * c3 <= X:
                    tot += 1
                h2 += 1
            h1 += 1
        if sheet == 0:
            break
        s += 1
    return tot


def session_every_base(q, b, K):
    primes = np.array([3, 5, 7, 11, 13], dtype=np.int64)
    cnt = np.zeros((len(primes), 16), np.int64)
    div = np.zeros((len(primes), 2), np.int64)
    digit2 = np.zeros(b, np.int64)
    EB.census_b(q, b, 1 << K, K, primes, cnt, div, digit2)
    return int(cnt[0].sum())


def main():
    NB = 24
    psi = psi_table(1 << NB)
    h = np.arange((1 << NB) + 1, dtype=np.int64)
    q = (h - psi).astype(np.int64)          # q(h) = h - psi(h)
    print("A. every_base.py census_b, sheet s >= 1, N_b <= 2^K: number of ordered triples")
    print("   K   b   python brute force   my census   session census_b (table 2^21, as in the script)   session census_b (full table 2^24)")
    for K in (14, 16):
        for b in (2, 3, 5, 10):
            bf = brute_python(psi, b, 1 << K, 1)
            mine = my_census_b(psi, b, 1 << K, 1, PRIMES)[4]
            s21 = session_every_base(q[: 1 << 21], b, K)
            sfull = session_every_base(q[: 1 << 24], b, K)
            print("  %2d  %2d   %12d   %12d   %12d   %12d   %s" % (K, b, bf, mine, s21, sfull,
                  "session census misses %d triples (defect (ii))" % (bf - sfull) if bf != sfull else "equal"))
    for K in (20, 24):
        for b in (2, 3, 4, 5, 6, 10):
            mine = my_census_b(psi, b, 1 << K, 1, PRIMES)[4]
            s21 = session_every_base(q[: 1 << 21], b, K)
            sfull = session_every_base(q[: 1 << 24], b, K)
            print("  %2d  %2d   %12s   %12d   %12d   %12d   missing with the script's table: %d (%.3f %%); with the full table: %d" % (
                K, b, "-", mine, s21, sfull, mine - s21, 100.0 * (mine - s21) / mine, mine - sfull))
    print()
    print("B. bridges_exact.py S10 base_head_census, sheet s = 0, N_beta <= 2^22: number of ordered triples")
    print("   beta   my census (all triples)   session function, table 2^12 (as in the script)   session function, table 2^22")
    X = 1 << 22
    for beta in (2, 3, 5, 10):
        primes = np.array([p for p in (3, 5, 7, 11, 13) if beta % p], np.int64)
        mine = my_census_b(psi, beta, X, 0, PRIMES)[4]
        out = []
        for size in (1 << 12, 1 << 22):
            cnt = np.zeros((len(primes), 2, 16), np.int64)
            exact = np.zeros(len(primes), np.int64)
            tot = BE.base_head_census(psi[:size].astype(np.int64), beta, X, primes, cnt, exact)
            cells = []
            for j, p in enumerate(primes):
                p = int(p)
                c0 = (-beta * pow(2, -1, p)) % p
                D1, Dp = int(cnt[j, 0].sum()), int(cnt[j, 1].sum())
                rem = max(abs(int(cnt[j, 0, c] + cnt[j, 1, c]) - (D1 / p + (Dp if c == c0 else 0))) for c in range(p)) * p / tot
                fac0 = (cnt[j, 0, c0] + cnt[j, 1, c0]) * p / tot
                cells.append("p=%d: class %d factor %.3f, F p^2 = %.2f, rem %.3f" % (p, c0, fac0, Dp / tot * p * p, rem))
            out.append((int(tot), "; ".join(cells), int(exact.sum())))
        print("   %2d   %10d   %10d (%.1f %% of the triples)   %10d" % (beta, mine, out[0][0], 100.0 * out[0][0] / mine, out[1][0]))
        print("        table 2^12 (stored output): " + out[0][1])
        print("        table 2^22 (all triples):   " + out[1][1] + "   [exactness violations %d]" % out[1][2])
    print()
    print("C. brute-force check of my sheet-0 census at small size: beta, python brute force, my census")
    for beta in (2, 3, 5, 10):
        bf = brute_python(psi, beta, 1 << 14, 0)
        mine = my_census_b(psi, beta, 1 << 14, 0, PRIMES)[4]
        print("   %2d   %8d   %8d   %s" % (beta, bf, mine, "equal" if bf == mine else "DIFFERENT"))


if __name__ == "__main__":
    main()
