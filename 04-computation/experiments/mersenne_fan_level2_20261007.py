"""Level 2 below the t = 23 fan, pushed to depth 2^LOGW (session opus-2026-10-07-S21).
Note: 05-knowledge/results/mersenne_line_barriers_20261007.md.  Uses par/stepf from mersenne_fan_escape_20261007.py.
The 2084 chains D <= 4000 of the fan bottom M_99708993677 that are not absorbed at 19,000,765 form one cluster from step
4,474,989 on (mersenne_fan_block_20261007.py (K)).  Here the D = 1911 chain is run from step 0 to 2^LOGW - 8 with coins
streamed from the fast parity vector, and the D = 4000 chain alongside it until the two collapse; absorption at time n would
certify M_99708993677 ~> M_(99708993677 - D) for every D in the cluster.  Absorption is checked at every Terras time
<= 2^LOGW - 9.  Prints ALL CHECKS PASSED (the check is only that the run is internally consistent: the two chains collapse,
as in the block script) and reports the outcome.  LOGW = 29 takes about 21 minutes."""
import sys, os, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mersenne_fan_escape_20261007 as F          # noqa: E402
from flint import fmpz                              # noqa: E402

LOGW = int(sys.argv[1]) if len(sys.argv) > 1 else 29
E0 = 99708993677
FAIL = []


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


if __name__ == '__main__':
    t0 = time.time()
    N = (1 << LOGW) - 8
    V, a, c = F.par(F.xres(E0, N + 8), N)
    del c
    nb = (N + 7) // 8
    Vb = int(V).to_bytes(nb, 'little')
    del V
    print(f'parity vector of the fan bottom to 2^{LOGW}: {time.time() - t0:.0f}s, weight {a}', flush=True)
    k1, a1 = -1911, 1 - F.P3[1911]
    k2, a2 = -4000, 1 - F.P3[4000]
    collapsed = None; absorbed = None; kmin = k1; kmax = k1
    n = 0
    BLK = 1 << 20                                     # bytes per block = 8M coins
    for off in range(0, nb, BLK):
        coins = np.unpackbits(np.frombuffer(Vb[off:off + BLK], dtype=np.uint8), bitorder='little').tolist()
        for b in coins:
            if n >= N:
                break
            if k1 == 0 and a1 == 0:
                absorbed = n
                break
            k1, a1 = F.stepf(k1, a1, b)
            if collapsed is None:
                k2, a2 = F.stepf(k2, a2, b)
                if (k2, a2) == (k1, a1):
                    collapsed = n + 1
            if k1 < kmin:
                kmin = k1
            if k1 > kmax:
                kmax = k1
            n += 1
        if absorbed is not None or n >= N:
            break
        if (off // BLK) % 8 == 7:
            print(f'  step {n}: level {k1}, range [{kmin}, {kmax}]  [{time.time() - t0:.0f}s]', flush=True)
    print(f'D = 4000 chain collapsed into the D = 1911 chain at step {collapsed}')
    print(f'outcome: ' + (f'ABSORBED at Terras time {absorbed}' if absorbed is not None else f'not absorbed within {N} steps') + f'; level range [{kmin}, {kmax}]')
    check(collapsed is not None and collapsed <= (1 << 27), 'the D = 4000 chain collapses into the D = 1911 chain before 2^27 (consistent with the companion script)')
    print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')
