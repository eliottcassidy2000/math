"""Fast two-run observer by subset dynamic programming (numpy), for whole censuses at n = 7, 8.

hp[S][s,t] = number of Hamiltonian paths of T[S] from s to t;  M_2[(s1,s2),(t1,t2)] = sum_{S1} hp[S1][s1,t1] hp[V-S1][s2,t2].
Cross-checked against the permutation-cutting construction of runs.py for every class with n <= 6.
Records, per class: strong?, number of Hamiltonian cycles, Smith form of M_2 (PARI), and the field "unglued=", which
is the concatenation of the Smith forms of the plus and minus sectors P, S (not re-sorted into a Smith form of N;
spur_2 is unaffected) -- audit data-file note.
usage: python3 fast_runs.py n [--check]
"""
import sys
import itertools
import numpy as np
from runs import tournaments, tour_from_bits, run_matrix, gp_matrix, gp_batch, parse_vec, snf_summary, fmt_snf


def hp_table(A):
    n = len(A)
    An = np.array(A, dtype=np.int64)
    hp = {}
    for v in range(n):
        M = np.zeros((n, n), dtype=np.int64)
        M[v, v] = 1
        hp[1 << v] = M
    for size in range(2, n + 1):
        for comb in itertools.combinations(range(n), size):
            S = 0
            for v in comb:
                S |= 1 << v
            M = np.zeros((n, n), dtype=np.int64)
            for t in comb:
                prev = hp[S & ~(1 << t)]
                M[:, t] = prev @ An[:, t]
            # rows outside S are zero automatically (prev rows outside S are zero)
            hp[S] = M
    return hp


def m2_fast(A):
    n = len(A)
    hp = hp_table(A)
    full = (1 << n) - 1
    T = np.zeros((n, n, n, n), dtype=np.int64)    # T[s1,s2,t1,t2]
    for S1 in range(1, full):
        H1 = hp[S1]
        H2 = hp[full & ~S1]
        T += np.einsum('ac,bd->abcd', H1, H2)
    return T, hp


def ham_cycles(A, hp):
    n = len(A)
    full = (1 << n) - 1
    # cycles through vertex 0: paths 0 -> t covering V with t -> 0
    return int(sum(hp[full][0, t] for t in range(1, n) if A[t][0]))


def tuples2(n):
    out = []
    for b in itertools.combinations(range(n), 2):
        out.append((b[0], b[1]))
        out.append((b[1], b[0]))
    return out


def main():
    n = int(sys.argv[1])
    check = '--check' in sys.argv
    tup = tuples2(n)
    recs = []
    cmds = []
    for bits in tournaments(n):
        A = tour_from_bits(n, bits)
        T, hp = m2_fast(A)
        M = [[int(T[s[0], s[1], t[0], t[1]]) for t in tup] for s in tup]
        if check:
            _, M0 = run_matrix(A, 2, tup)
            assert M == M0, bits
        m = len(tup) // 2
        P = [[M[2 * a][2 * b] + M[2 * a][2 * b + 1] for b in range(m)] for a in range(m)]
        S = [[M[2 * a][2 * b] - M[2 * a][2 * b + 1] for b in range(m)] for a in range(m)]
        full = (1 << n) - 1
        strong = all(hp_reach for hp_reach in [True])
        # strong iff a Hamiltonian cycle exists (Camion), n >= 3
        hc = ham_cycles(A, hp)
        recs.append((bits, hc))
        cmds.append('matsnf(%s)' % gp_matrix(M))
        cmds.append('concat(matsnf(%s),matsnf(%s))' % (gp_matrix(P), gp_matrix(S)))
    out = []
    step = 60
    for i in range(0, len(cmds), step):
        out += gp_batch(cmds[i:i + step])
    stats = {}
    for c, (bits, hc) in enumerate(recs):
        snf = snf_summary(parse_vec(out[2 * c]))
        un = snf_summary(parse_vec(out[2 * c + 1]))
        tf = all(d in (0, 1) for d in snf)
        key = ('HC=0' if hc == 0 else 'HC=1' if hc == 1 else 'HC>=2', 'torsion-free' if tf else 'torsion')
        stats[key] = stats.get(key, 0) + 1
        a2 = sum(mm for d, mm in snf.items() if d and d % 2 == 0)
        b2 = sum(mm for d, mm in un.items() if d and d % 2 == 0)
        print('%s HC=%d M2=%s unglued=%s spur2=%d' % (bits, hc, fmt_snf(snf), fmt_snf(un), b2 - a2))
    print('SUMMARY n=%d (two-run observer):' % n, ', '.join('%s/%s=%d' % (k[0], k[1], v) for k, v in sorted(stats.items())))


if __name__ == '__main__':
    main()
