"""Census of k-run response matrices over all tournaments of order n.

For each isomorphism class (gentourng) and each run count k this records
  * the ordered observer M_k (rank over Q, Smith form for k <= 3 when requested),
  * the plus sector P_k (unordered start set x end set path-cover counts; P_k[a,b] counts unordered
    covers by k paths with start set a and end set b), its rank, and whether P_k kills every boundary of a
    (k+1)-subset (the (k-1)-chains on k-sets, run-symmetric lifts),
  * for k = 2, 3 the unglued (sector-split) Smith form and the spurious p-factor counts
    spur_p = rank_Fp(M_k) - rank_Fp(N_k)  (equal to the excess of invariant factors divisible by p).

usage: python3 census_runs.py n k [--snf]
"""
import sys
import itertools
from runs import (tournaments, tour_from_bits, run_matrix, ordered_tuples, unglued_k3, gp_matrix,
                  gp_batch, parse_vec, snf_summary, fmt_snf)


def scores(A):
    return ''.join(str(s) for s in sorted(sum(r) for r in A))


def is_strong(A):
    n = len(A)
    reach = [{i} for i in range(n)]
    changed = True
    while changed:
        changed = False
        for i in range(n):
            for j in range(n):
                if A[i][j] and not reach[j] <= reach[i]:
                    reach[i] |= reach[j]
                    changed = True
    return all(len(r) == n for r in reach)


def ham_paths(A):
    n = len(A)
    return sum(1 for w in itertools.permutations(range(n)) if all(A[w[i]][w[i + 1]] for i in range(n - 1)))


def plus_sector(A, k):
    """Unordered covers: start set x end set counts (k-subsets in lexicographic order)."""
    n = len(A)
    sets = list(itertools.combinations(range(n), k))
    idx = {s: i for i, s in enumerate(sets)}
    P = [[0] * len(sets) for _ in sets]
    for w in itertools.permutations(range(n)):
        breaks = [i for i in range(n - 1) if not A[w[i]][w[i + 1]]]
        if len(breaks) > k - 1:
            continue
        free = [i for i in range(n - 1) if A[w[i]][w[i + 1]]]
        for extra in itertools.combinations(free, k - 1 - len(breaks)):
            cuts = sorted(breaks + list(extra))
            starts = [w[0]] + [w[c + 1] for c in cuts]
            if starts != sorted(starts):      # one ordering of the k paths per unordered cover
                continue
            ends = tuple(sorted([w[c] for c in cuts] + [w[n - 1]]))
            P[idx[tuple(starts)]][idx[ends]] += 1
    return sets, P


def boundaries(n, k, sets):
    """Boundaries of the (k+1)-subsets as (k-1)-chains on k-sets (columns of the returned list)."""
    idx = {s: i for i, s in enumerate(sets)}
    out = []
    for Q in itertools.combinations(range(n), k + 1):
        v = [0] * len(sets)
        for pos, y in enumerate(Q):
            face = tuple(x for x in Q if x != y)
            v[idx[face]] += -1 if pos % 2 else 1
        out.append(v)
    return out


def mat_vec_zero(P, v):
    return all(sum(P[i][j] * v[j] for j in range(len(v))) == 0 for i in range(len(P)))


def main():
    n, k = int(sys.argv[1]), int(sys.argv[2])
    want_snf = '--snf' in sys.argv
    classes = tournaments(n)
    recs = []
    cmds = []
    for bits in classes:
        A = tour_from_bits(n, bits)
        sets, P = plus_sector(A, k)
        B = boundaries(n, k, sets)
        kills = all(mat_vec_zero(P, v) for v in B)
        rec = {'bits': bits, 'scores': scores(A), 'strong': is_strong(A), 'H': ham_paths(A), 'kills': kills}
        cmds.append('matrank(%s)' % gp_matrix(P))
        rec['cmd0'] = len(cmds) - 1
        if not want_snf:
            tuples, M = run_matrix(A, k)
            cmds.append('matrank(%s)' % gp_matrix(M))
            rec['cmd_rankM'] = len(cmds) - 1
        if want_snf:
            tuples, M = run_matrix(A, k)
            cmds.append('matsnf(%s)' % gp_matrix(M))
            rec['cmd_snf'] = len(cmds) - 1
            if k == 3:
                N, Pp, S, L0, L1 = unglued_k3(M, len(tuples) // 6)
                cmds.append('matsnf(%s)' % gp_matrix(N))
                rec['cmd_un'] = len(cmds) - 1
                for nm, X in (('P', Pp), ('S', S), ('L', L0)):
                    cmds.append('matsnf(%s)' % gp_matrix(X))
                    rec['cmd_' + nm] = len(cmds) - 1
            if k == 2:
                # S_2 chart on each unordered pair: (u, v) = (x_ab + x_ba, x_ab - x_ba)
                m = len(tuples) // 2
                Fb = [[0] * (2 * m) for _ in range(2 * m)]
                for o in range(m):
                    Fb[2 * o][2 * o] = Fb[2 * o][2 * o + 1] = 1
                    Fb[2 * o + 1][2 * o] = 1
                    Fb[2 * o + 1][2 * o + 1] = -1
                Hb = [[x for x in row] for row in Fb]       # F^-1 = F / 2
                from runs import matmul
                X = matmul(matmul(Fb, M), Hb)
                N = [[v // 2 for v in row] for row in X]
                assert all(v % 2 == 0 for row in X for v in row)
                cmds.append('matsnf(%s)' % gp_matrix(N))
                rec['cmd_un'] = len(cmds) - 1
        recs.append(rec)
    out = []
    step = 40
    for i in range(0, len(cmds), step):
        out += gp_batch(cmds[i:i + step])
    print('n=%d k=%d classes=%d  (plus sector = unordered start-set x end-set covers; boundary space dim C(%d,%d)=%d)'
          % (n, k, len(classes), n - 1, k, len(boundaries(n, k, list(itertools.combinations(range(n), k))))
             and __import__('math').comb(n - 1, k)))
    for rec in recs:
        line = '%-22s %-8s %-6s H=%-4d rankP=%-4s killsBdry=%-5s' % (
            rec['bits'], rec['scores'], 'strong' if rec['strong'] else 'nonstr', rec['H'], out[rec['cmd0']],
            rec['kills'])
        if 'cmd_rankM' in rec:
            line += ' rankM=%s/%d' % (out[rec['cmd_rankM']], len(run_matrix(tour_from_bits(n, rec['bits']), k)[0]))
        if want_snf:
            snf = snf_summary(parse_vec(out[rec['cmd_snf']]))
            line += ' M=%s' % fmt_snf(snf)
            if 'cmd_un' in rec:
                un = snf_summary(parse_vec(out[rec['cmd_un']]))
                sp = {}
                for p in (2, 3):
                    a = sum(m for d, m in snf.items() if d and d % p == 0)
                    b = sum(m for d, m in un.items() if d and d % p == 0)
                    sp[p] = b - a
                line += ' unglued=%s spur2=%d spur3=%d' % (fmt_snf(un), sp[2], sp[3])
            for nm in ('P', 'S', 'L'):
                if 'cmd_' + nm in rec:
                    line += ' %s=%s' % (nm, fmt_snf(snf_summary(parse_vec(out[rec['cmd_' + nm]]))))
        print(line)


if __name__ == '__main__':
    main()
