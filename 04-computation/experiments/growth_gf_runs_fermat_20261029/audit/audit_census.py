"""Independent census for THM-4846 (E) / HYP-9445 (1),(2): per class, Smith forms of M_k (FLINT), glue losses,
ghost-lemma connecting rank, boundary killing, Hamiltonian cycles.

usage: python3 audit_census.py n k [--hc LIST]   (LIST: comma-separated HC counts to restrict to, e.g. 1,2,3)
Writes one line per class to stdout and a SUMMARY line.
"""
import sys
import itertools
from runlib import (tournaments, adj_from_bits, hp_table, ham_cycles, ham_paths, is_strong, run_matrix,
                    plus_sector, snf_diag, snf_summary, fmt, rank_mod_p, unglued, sector, ndiv)


def boundary_vectors(n, k, sets):
    idx = {s: i for i, s in enumerate(sets)}
    out = []
    for Q in itertools.combinations(range(n), k + 1):
        v = [0] * len(sets)
        for pos, y in enumerate(Q):
            face = tuple(x for x in Q if x != y)
            v[idx[face]] += -1 if pos % 2 else 1
        out.append(v)
    return out


def connecting_rank(P, S):
    r = len(P)
    B = []
    for i in range(r):
        row = []
        for j in range(r):
            d = S[i][j] - P[i][j]
            assert d % 2 == 0
            row.append(d // 2)
        B.append(row)
    block = [P[i] + [0] * r for i in range(r)] + [B[i] + P[i] for i in range(r)]
    return rank_mod_p(block, 2) - 2 * rank_mod_p(P, 2)


def main():
    n, k = int(sys.argv[1]), int(sys.argv[2])
    hc_filter = None
    if '--hc' in sys.argv:
        hc_filter = set(int(x) for x in sys.argv[sys.argv.index('--hc') + 1].split(','))
    classes = tournaments(n)
    if '--slice' in sys.argv:
        a, b = map(int, sys.argv[sys.argv.index('--slice') + 1].split(':'))
        classes = classes[a:b]
    stats = {}
    for bits in classes:
        A = adj_from_bits(n, bits)
        hp = hp_table(A)
        hc = ham_cycles(A, hp)
        strong = is_strong(A)
        assert strong == (hc > 0 or n == 1)
        if hc_filter is not None and hc not in hc_filter:
            continue
        H = ham_paths(hp, n)
        sc = ''.join(map(str, sorted(sum(r) for r in A)))
        tup, M = run_matrix(A, k, hp)
        sm = snf_summary(snf_diag(M))
        tf = all(d in (0, 1) for d, _ in sm)
        line = '%s %s %s H=%d HC=%d M=%s' % (bits, sc, 'strong' if strong else 'nonstr', H, hc, fmt(sm))
        rec = {'strong': strong, 'hc': hc, 'tf': tf}
        if k in (2, 3):
            N = unglued(M, k)
            sn = snf_summary(snf_diag(N))
            sp2 = ndiv(sn, 2) - ndiv(sm, 2)
            sp3 = ndiv(sn, 3) - ndiv(sm, 3)
            # rank form of the spurious counts
            r2 = rank_mod_p(M, 2) - rank_mod_p(N, 2)
            r3 = rank_mod_p(M, 3) - rank_mod_p(N, 3)
            assert (sp2, sp3) == (r2, r3)
            P = sector(N, k, [0])
            S = sector(N, k, [1])
            # plus sector from the chart = direct unordered start-set x end-set counts
            sets, Pd = plus_sector(A, k, hp)
            assert P == Pd, 'plus sector mismatch'
            conn = connecting_rank(P, S)
            assert conn == sp2, ('ghost lemma fails', bits, conn, sp2)
            loss = sm != sn
            line += ' N=%s spur2=%d spur3=%d conn2=%d gluediff=%s' % (fmt(sn), sp2, sp3, conn, loss)
            rec.update(spur2=sp2, spur3=sp3, loss=loss)
            if k == 3 or k == 2:
                kills = all(all(sum(P[i][j] * v[j] for j in range(len(v))) == 0 for i in range(len(P)))
                            for v in boundary_vectors(n, k, sets))
                line += ' kills=%s' % kills
                rec['kills'] = kills
            # torsion odd-only?
            rec['oddonly'] = (not tf) and all(d in (0, 1) or d % 2 == 1 for d, _ in sm)
        print(line, flush=True)
        key = ('strong' if strong else 'nonstr', 'HC=0' if hc == 0 else 'HC=1' if hc == 1 else 'HC>=2',
               'tf' if tf else 'torsion')
        stats[key] = stats.get(key, 0) + 1
        for nm in ('loss', 'kills', 'oddonly'):
            if rec.get(nm):
                stats[nm] = stats.get(nm, 0) + 1
        if rec.get('spur2', 0) or rec.get('spur3', 0):
            stats['spur>0'] = stats.get('spur>0', 0) + 1
        if strong and hc >= 2 and rec.get('oddonly'):
            stats['HC>=2 odd-only torsion'] = stats.get('HC>=2 odd-only torsion', 0) + 1
    print('SUMMARY n=%d k=%d classes=%d %s' % (n, k, len(classes), sorted(stats.items(), key=str)))


if __name__ == '__main__':
    main()
