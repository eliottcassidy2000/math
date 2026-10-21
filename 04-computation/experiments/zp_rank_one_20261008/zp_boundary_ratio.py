"""Basin boundaries against the Haar merge tail, at every scale N (session opus-2026-10-08-S22; added after the audit).
Note: 05-knowledge/results/zp_rank_one_coalescence_20261008.md (section 3); HYP-9245 (a).

For p = 3, 11, 17, 23 and N = 10^3 .. 10^7:
  * second-basin share s(N) and boundary density b(N) = #{n < N : n, n+1 in different basins}/N (basin_array of
    zp_basins_survey.py, n <= 10^7);
  * the Haar no-merge probability q_p(T_N) from the start (0, 1) at T_N = round(ln N/|Lambda_p|), from 10^5 chain paths
    (12 workers), with its binomial standard error;
  * the ratio R(N) = b(N)/(2 s(N)(1 - s(N)) q_p(T_N)) and its standard error from q alone.
NUMERICAL.  Prints ALL CHECKS PASSED.  Runtime about 30 seconds."""
import sys, os, math, random, time
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zp_chain_checks as Z          # noqa: E402
import zp_basins_survey as S         # noqa: E402

FAIL = []
NS = [10 ** j for j in range(3, 8)]
PATHS = 100000


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


def lam(p):
    return (1 / p) * math.log(1 / p) + ((p - 1) / p) * math.log((p + 1) / p)


def job(args):
    p, seed, n, Tmax = args
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        k, A = 0, 1; ab = None
        for t in range(Tmax + 1):
            if k == 0 and A == 0:
                ab = t; break
            k, A = Z.istep(p, k, A, rng.randrange(p))
        out.append(ab)
    return out


if __name__ == '__main__':
    t0 = time.time()
    rows = {}
    with Pool(12) as pool:
        for p in (3, 11, 17, 23):
            TN = {N: round(math.log(N) / abs(lam(p))) for N in NS}
            Tmax = max(TN.values())
            chunks = pool.map(job, [(p, 104729 * p + w, PATHS // 12, Tmax) for w in range(12)])
            ab = [x for c in chunks for x in c]
            M = len(ab)
            b = S.basin_array(p, 10 ** 7)
            line = []
            for N in NS:
                q = sum(1 for a in ab if a is None or a > TN[N]) / M
                se = math.sqrt(q * (1 - q) / M)
                s = b[1:N + 1].count(2) / N
                bd = sum(1 for n in range(1, N) if b[n] != b[n + 1]) / N
                R = bd / (2 * s * (1 - s) * q)
                line.append((N, TN[N], round(s, 5), round(bd, 5), round(q, 4), round(se, 4), round(R, 3), round(R * se / q, 3)))
            rows[p] = line
            print(f'     p = {p}: (N, T_N, s, boundary, q_p(T_N), se(q), ratio R, se(R)) = {line}  [{time.time() - t0:.0f}s]', flush=True)
    r3 = [x[6] for x in rows[3]]
    check(all(r3[i] > r3[i + 1] for i in range(len(r3) - 1)) and all(0.9 < rows[p][-1][6] < 1.3 for p in rows),
          'for p = 3 the ratio R(N) falls strictly from N = 10^3 to 10^7; at N = 10^7 it lies in (0.9, 1.3) for p = 3, 11, 17, 23 '
          '(NUMERICAL; whether R(N) tends to 1 is not resolved)')
    print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')
