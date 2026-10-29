"""Two-run plus and sign sectors versus Hamiltonian cycles (spot check of the producer's out/s2_vs_hc_n8.out).

P_2[{a1<a2},{b1<b2}] = M[(a1,a2),(b1,b2)] + M[(a1,a2),(b2,b1)],  S_2 = the same with a minus sign, M from the subset DP.
Smith forms by FLINT. usage: python3 audit_s2.py n [n ...]
"""
import sys
import itertools
from collections import Counter
from runlib import tournaments, adj_from_bits, hp_table, ham_cycles, run_matrix, snf_diag


def main():
    for n in [int(x) for x in sys.argv[1:]]:
        c = Counter()
        exc = []
        for bits in tournaments(n):
            A = adj_from_bits(n, bits)
            hp = hp_table(A)
            hc = ham_cycles(A, hp)
            if hc == 0:
                continue
            tup, M = run_matrix(A, 2, hp)
            idx = {t: i for i, t in enumerate(tup)}
            pairs = list(itertools.combinations(range(n), 2))
            P = [[M[idx[a]][idx[b]] + M[idx[a]][idx[(b[1], b[0])]] for b in pairs] for a in pairs]
            S = [[M[idx[a]][idx[b]] - M[idx[a]][idx[(b[1], b[0])]] for b in pairs] for a in pairs]
            dP, dS = snf_diag(P), snf_diag(S)
            tP, tS = any(d > 1 for d in dP), any(d > 1 for d in dS)
            sS = 0 in dS
            key = ('HC=1' if hc == 1 else 'HC>=2', 'S torsion' if tS else 'S tf', 'S singular' if sS else 'S nonsing',
                   'P torsion' if tP else 'P tf')
            c[key] += 1
            if (hc == 1) == tS:
                exc.append((bits, hc, sorted(Counter(dS).items())))
        print('n=%d strong classes: %s' % (n, sorted(c.items())))
        print('   violations of [S_2 torsion-free <=> unique HC]: %d %s' % (len(exc), exc[:5]), flush=True)


if __name__ == '__main__':
    main()
