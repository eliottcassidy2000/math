"""Merge-time tails of the base-p pair chain (session opus-2026-10-08-S22).
Note: 05-knowledge/results/zp_rank_one_coalescence_20261008.md (section 5).  Uses istep from zp_chain_checks.py.

Start (0, 1): Haar y and y + 1.  Digits are fresh uniform (Fact D), so the chain is simulated exactly from random digits.
For p = 2, 3, 7, 11, 13, 17, 23: 4800 paths to T = 10^5; q(T) = P(no merge by T) and sqrt(T) q(T) (THM-4610 5(c):
q >= c T^(-1/2); the T^(-1/2) upper rate is not claimed).  Also q at T_N = ln(10^7)/|Lambda_p|, the typical number of steps
an integer near 10^7 takes to descend, to compare with q_eff of zp_basins_survey.py (S2).  Also the mean flip rate (2/p per step away from p | A) and the number of
flips by T, so that q can be read against the flip clock.  NUMERICAL.  Prints ALL CHECKS PASSED.  Runtime about 10 minutes."""
import sys, os, math, random, time
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zp_chain_checks as Z          # noqa: E402

FAIL = []
TS = (10, 100, 1000, 10 ** 4, 10 ** 5)


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


def job(args):
    p, seed, n_paths, Tmax = args
    rng = random.Random(seed)
    out = []
    for _ in range(n_paths):
        k, A = 0, 1; ab = None; flips = 0
        for n in range(Tmax):
            if k == 0 and A == 0:
                ab = n; break
            k2, A = Z.istep(p, k, A, rng.randrange(p))
            flips += k2 != k
            k = k2
        out.append((ab, flips))
    return out


if __name__ == '__main__':
    t0 = time.time()
    res = {}
    with Pool(12) as pool:
        for p in (2, 3, 7, 11, 13, 17, 23):
            chunks = pool.map(job, [(p, 7919 * p + w, 400, max(TS)) for w in range(12)])
            r = [x for c in chunks for x in c]
            N = len(r)
            row = []
            for T in TS:
                q = sum(1 for ab, _ in r if ab is None or ab > T) / N
                row.append((T, round(q, 5), round(math.sqrt(T) * q, 3)))
            alive = [f for ab, f in r if ab is None]
            lam = (1 / p) * math.log(1 / p) + ((p - 1) / p) * math.log((p + 1) / p)
            TN = round(math.log(10 ** 7) / abs(lam))
            qTN = sum(1 for ab, _ in r if ab is None or ab > TN) / N
            row.append(('T_N', TN, round(qTN, 4)))
            res[p] = row
            print(f'     p = {p}: (T, q, sqrt(T) q) = {row}; unmerged at 10^5: {len(alive)} of {N}, '
                  f'their mean flips {sum(alive) / max(1, len(alive)):.0f}  [{time.time() - t0:.0f}s]', flush=True)
    check(all(res[p][-2][1] > 0 and res[p][-2][1] < res[p][0][1] for p in res),
          'q(T) decreases and stays positive to T = 10^5 for p = 2, 3, 7, 11, 13, 17, 23 (NUMERICAL)')
    print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')
