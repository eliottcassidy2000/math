"""Fast torsion / spurious-factor census without full Smith forms (for n = 7, 8 with k = 3).

Per class: rational rank r (FLINT), ranks of M and of the unglued N mod 2 and 3 (spur_2, spur_3, and the p = 2 ghost
connecting rank), torsion decided as follows:
  * torsion if rank_Fp(M) < r for some prime p <= 97 (certificate: that p);
  * otherwise, if M is nonsingular, torsion-free iff |det M| = 1 (FLINT det);
  * otherwise an r x r minor on independent rows/columns (chosen mod a large prime) with |det| = 1 certifies
    torsion-free (the gcd of maximal minors divides it); if none is found the full FLINT Smith form decides.
usage: python3 audit_fast.py n k [--hc LIST]
"""
import sys
import random
import flint
from runlib import (tournaments, adj_from_bits, hp_table, ham_cycles, is_strong, run_matrix, rank_mod_p, rank_Q,
                    unglued, sector, snf_diag)
from audit_census import connecting_rank

PRIMES = [p for p in range(2, 98) if all(p % q for q in range(2, int(p ** 0.5) + 1))]
BIGP = 2305843009213693951


def independent_rows(M, p):
    rows = []
    basis = []
    for i, row in enumerate(M):
        cand = basis + [[x % p for x in row]]
        if rank_mod_p(cand, p) > len(basis):
            basis = cand
            rows.append(i)
    return rows


def torsion_status(M, r):
    for p in PRIMES:
        if rank_mod_p(M, p) < r:
            return True, 'rank drop mod %d' % p
    n = len(M)
    if r == n:
        d = int(flint.fmpz_mat(M).det())
        return abs(d) != 1, 'det=%d' % d if abs(d) < 10 ** 6 else 'det has %d digits' % len(str(abs(d)))
    rnd = random.Random(1)
    for attempt in range(6):
        perm = list(range(n))
        if attempt:
            rnd.shuffle(perm)
        Mp = [M[i] for i in perm]
        R = [perm[i] for i in independent_rows(Mp, BIGP)][:r]
        sub = [M[i] for i in R]
        T = [list(col) for col in zip(*sub)]
        cperm = list(range(n))
        if attempt:
            rnd.shuffle(cperm)
        C = [cperm[i] for i in independent_rows([T[j] for j in cperm], BIGP)][:r]
        minor = [[M[i][j] for j in C] for i in R]
        d = int(flint.fmpz_mat(minor).det())
        if abs(d) == 1:
            return False, 'unit maximal minor'
    d = snf_diag(M)
    return any(x > 1 for x in d), 'full SNF'


def main():
    n, k = int(sys.argv[1]), int(sys.argv[2])
    hc_filter = None
    if '--hc' in sys.argv:
        hc_filter = set(int(x) for x in sys.argv[sys.argv.index('--hc') + 1].split(','))
    stats = {}
    for bits in tournaments(n):
        A = adj_from_bits(n, bits)
        hp = hp_table(A)
        hc = ham_cycles(A, hp)
        if hc_filter is not None and hc not in hc_filter:
            continue
        strong = is_strong(A)
        tup, M = run_matrix(A, k, hp)
        r = rank_Q(M)
        N = unglued(M, k)
        sp2 = rank_mod_p(M, 2) - rank_mod_p(N, 2)
        sp3 = rank_mod_p(M, 3) - rank_mod_p(N, 3)
        conn = connecting_rank(sector(N, k, [0]), sector(N, k, [1]))
        assert conn == sp2
        tor, how = torsion_status(M, r)
        print('%s %s HC=%d rank=%d/%d spur2=%d spur3=%d torsion=%s (%s)'
              % (bits, 'strong' if strong else 'nonstr', hc, r, len(M), sp2, sp3, tor, how), flush=True)
        key = ('strong' if strong else 'nonstr', 'HC=0' if hc == 0 else 'HC=1' if hc == 1 else 'HC>=2',
               'torsion' if tor else 'tf')
        stats[key] = stats.get(key, 0) + 1
        if sp2 or sp3:
            stats['spur>0'] = stats.get('spur>0', 0) + 1
        hk = 'HC=%d' % hc if hc <= 3 else None
        if hk:
            stats[(hk, 'torsion' if tor else 'tf')] = stats.get((hk, 'torsion' if tor else 'tf'), 0) + 1
    print('SUMMARY n=%d k=%d %s' % (n, k, sorted(stats.items(), key=str)))


if __name__ == '__main__':
    main()
