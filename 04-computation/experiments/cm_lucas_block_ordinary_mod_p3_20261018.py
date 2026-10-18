#!/usr/bin/env python3
"""Third-order Lucas factorization at ordinary primes of CM points (S43, refining HYP-9310/9311).

Probe (cm_lucas_block_mod_p3_probe_20261018.out): (B1') sum_{k=p}^{2p-1} t_k == t_1 sum_{k<p} t_k held mod p^3 at every
split prime and failed at every inert prime of the CM field for eight parameters checked (one exception: A1, m = 1,
p = 17).
PRE-REGISTERED (written before this run): for all 33 CM parameters of the four families and all primes
7 <= p <= PMAX not dividing 6m:
  (O)  if p is ordinary (F_p = sum_{k<p} t_k is a p-adic unit), then (B1') holds modulo p^3;
  (S)  if p is supersingular (F_p == 0 mod p), (B1') holds mod p^2 but fails mod p^3 except at sporadic primes
       (expected frequency about 1/p).
Usage: python cm_lucas_block_ordinary_mod_p3_20261018.py [PMAX]
Session opus-2026-10-17-S43.
"""
import sys, json, os
from math import comb

FAM = {
    'A1': lambda k: comb(2 * k, k) ** 3,
    'A2': lambda k: comb(2 * k, k) ** 2 * comb(3 * k, k),
    'A3': lambda k: comb(2 * k, k) ** 2 * comb(4 * k, 2 * k),
    'A4': lambda k: comb(2 * k, k) * comb(3 * k, k) * comb(6 * k, 3 * k),
}
SURV = {
    'A1': [-512, -64, -8, 1, 16, 64, 256, 4096],
    'A2': [-8640, -1728, -192, -27, 8, 64, 108, 216, 1458, 3375],
    'A3': [-82944, -12288, -3969, -1024, -144, 81, 256, 648, 2304, 20736],
    'A4': [-32768, -3375, 1728, 8000, 54000],
}


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    primes = [p for p in range(7, PMAX + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    out = {}
    totO = totOfail = totS = totShold = 0
    for fam, Af in FAM.items():
        A = [Af(k) for k in range(2 * PMAX)]
        res = {}
        for m in SURV[fam]:
            O_fail, S_hold, nO, nS = [], [], 0, 0
            for p in primes:
                if (6 * m) % p == 0:
                    continue
                P3 = p ** 3
                im = pow(m % P3, -1, P3)
                pw = [1] * (2 * p)
                for k in range(1, 2 * p):
                    pw[k] = pw[k - 1] * im % P3
                S0 = sum(A[k] % P3 * pw[k] for k in range(p)) % P3
                S1 = sum(A[k] % P3 * pw[k] for k in range(p, 2 * p)) % P3
                d = (S1 - A[1] * im * S0) % P3
                assert d % (p * p) == 0, (fam, m, p)
                if S0 % p:
                    nO += 1
                    if d:
                        O_fail.append(p)
                else:
                    nS += 1
                    if d == 0:
                        S_hold.append(p)
            res[str(m)] = dict(ordinary=nO, ordinary_fail_mod_p3=O_fail, supersingular=nS, supersingular_hold_mod_p3=S_hold)
            totO += nO
            totOfail += len(O_fail)
            totS += nS
            totShold += len(S_hold)
            print('%s m = %7d: ordinary %3d (mod p^3 failures %s) | supersingular %3d (mod p^3 holds at %s)'
                  % (fam, m, nO, O_fail if O_fail else 'none', nS, S_hold), flush=True)
        out[fam] = res
    print('TOTAL: ordinary prime-instances %d, mod p^3 failures %d; supersingular %d, sporadic mod p^3 holds %d'
          % (totO, totOfail, totS, totShold))
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'cm_lucas_block_ordinary_mod_p3_20261018.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
