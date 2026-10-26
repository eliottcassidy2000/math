# Audit check A: TT5(P_q), TT6(P_q) by two algorithms that differ from the lane's (the lane uses the unique SOURCE
# recursively: TT_k = q * TT_(k-1)(T[QR]), then out-degree sums inside nested out-neighbourhoods).
#   (1) Source-SINK decomposition: a transitive k-set has a unique source s and sink t, s -> t, and its middle k-2
#       vertices are a transitive set inside M(s,t) = N+(s) cap N-(t).  Aut(P_q) contains x -> a x + b (a a square),
#       which is transitive on the q(q-1)/2 arcs, so TT_k = q(q-1)/2 * TT_(k-2)(T[M(0,1)]),
#       M(0,1) = {x : chi(x) = 1, chi(1-x) = 1}.  Inside M we use TT3 = sum_(s->t) (A^2)_st and
#       TT4 = sum_(s->t) C((A^2)_st, 2)  (the source-sink formula once more), A the adjacency of T[M].
#   (2) Plain depth-first enumeration of all transitive chains v1 -> v2 -> ... -> vk (all pairs i<j arcs) with
#       Python-integer bitsets over ALL vertices, no symmetry used (small q only).
# Prime powers use the audit's own GF(p^r) (audit_gf.py; a different irreducible polynomial from the lane's).
# Compares with the lane's paley_tt5_data.out, paley_tt6_data.out, paley_prime_powers.out.  Audit 2026-10-25.
import re, sys, time
import numpy as np
from sympy import primerange
sys.path.insert(0, '.')
from audit_gf import GF

LANE = '..'


def paley_M_adjacency(F):
    q = F.q
    els = np.arange(q)
    one_minus = F.sub(np.full(q, 1), els)
    M = els[(F.chi[els] == 1) & (F.chi[one_minus] == 1)]
    assert len(M) == (q - 3) // 4
    D = F.sub(M[None, :], M[:, None])           # D[s, t] = t - s
    A = (F.chi[D] == 1).astype(np.int64)
    assert (A + A.T + np.eye(len(M), dtype=np.int64) == 1).all()
    return A


def tt56_source_sink(F):
    q = F.q
    A = paley_M_adjacency(F)
    A2 = A @ A
    tt3 = int((A * A2).sum())
    tt4 = int((A * (A2 * (A2 - 1) // 2)).sum())
    arcs = q * (q - 1) // 2
    return arcs * tt3, arcs * tt4


def paley_out_bitsets(F):
    q = F.q
    els = np.arange(q)
    out = []
    for v in range(q):
        diff = F.sub(els, np.full(q, v))
        nb = np.nonzero(F.chi[diff] == 1)[0]
        m = 0
        for u in nb:
            m |= 1 << int(u)
        out.append(m)
    return out


def dfs_count(out, k):
    q = len(out)
    def rec(d, mask):
        if d == k - 1:
            return bin(mask).count('1')
        tot = 0
        m = mask
        while m:
            low = m & -m
            u = low.bit_length() - 1
            tot += rec(d + 1, mask & out[u])
            m ^= low
        return tot
    return sum(rec(1, out[v]) for v in range(q))


def load(path, pat):
    d = {}
    for l in open(path):
        m = re.match(pat, l)
        if m:
            d[int(m.group(1))] = int(m.group(2))
    return d


if __name__ == '__main__':
    lane5 = load(f'{LANE}/paley_tt5_data.out', r'p=(\d+) TT5=(\d+)')
    lane6 = load(f'{LANE}/paley_tt6_data.out', r'p=(\d+) TT6=(\d+)')
    lanepp = {}
    for l in open(f'{LANE}/paley_prime_powers.out'):
        m = re.match(r'q = \d+\^\d+ = (\d+): TT5 = (\d+) .*TT6 = (\d+) ', l)
        if m:
            lanepp[int(m.group(1))] = (int(m.group(2)), int(m.group(3)))

    t0 = time.time()
    print('(1) source-sink formula vs the lane (primes)')
    mism5 = mism6 = 0
    mine5, mine6 = {}, {}
    for p in primerange(3, 1500):
        if p % 4 != 3:
            continue
        F = GF(p, 1)
        a, b = tt56_source_sink(F)
        mine5[p] = a
        if p <= 751:
            mine6[p] = b
        if p in lane5 and lane5[p] != a:
            mism5 += 1; print('  TT5 MISMATCH', p, a, lane5[p])
        if p in lane6 and lane6[p] != b:
            mism6 += 1; print('  TT6 MISMATCH', p, b, lane6[p])
    print(f'  TT5: compared at {sum(1 for p in lane5 if p in mine5)} primes (lane file has {len(lane5)}), mismatches {mism5}')
    print(f'  TT6: compared at {sum(1 for p in lane6 if p in mine6)} primes (lane file has {len(lane6)}), mismatches {mism6}')
    print('  small q (not in the lane data files, inserted there as zeros):',
          {p: (mine5[p], mine6[p]) for p in (3, 7, 11)})
    print('  zeros of TT6 among primes <= 751:', [p for p in mine6 if mine6[p] == 0])
    print('  zeros of TT5 among primes < 1500:', [p for p in mine5 if mine5[p] == 0])

    print('(1b) prime powers (audit GF(p^r)) vs paley_prime_powers.out')
    pp = {}
    for (p, r) in [(3, 3), (3, 5), (7, 3), (11, 3), (3, 7)]:
        F = GF(p, r)
        a, b = tt56_source_sink(F)
        pp[F.q] = (a, b)
        lp = lanepp.get(F.q)
        print(f'  q = {p}^{r} = {F.q}: TT5 = {a}, TT6 = {b}; lane {lp}; {"AGREE" if lp == (a, b) else "DISAGREE"}  '
              f'(f = {F.f})')

    print('(2) depth-first enumeration of transitive chains (no symmetry) vs (1)')
    for q, (p, r) in [(7, (7, 1)), (11, (11, 1)), (19, (19, 1)), (23, (23, 1)), (27, (3, 3)), (31, (31, 1)),
                      (43, (43, 1)), (47, (47, 1)), (59, (59, 1)), (67, (67, 1)), (71, (71, 1))]:
        F = GF(p, r)
        out = paley_out_bitsets(F)
        c5, c6 = dfs_count(out, 5), dfs_count(out, 6)
        ref = (mine5[q], mine6[q]) if r == 1 else pp[q]
        print(f'  q = {q}: DFS TT5 = {c5}, TT6 = {c6}; source-sink {ref}; {"AGREE" if (c5, c6) == ref else "DISAGREE"}')
    print(f'time {time.time() - t0:.1f}s')

    # write my values for later checks
    with open('audit_tt_values.out', 'w') as fh:
        for p in sorted(mine5):
            fh.write(f'p={p} TT5={mine5[p]}' + (f' TT6={mine6[p]}' if p in mine6 else '') + '\n')
        for q, (a, b) in pp.items():
            fh.write(f'q={q} TT5={a} TT6={b}\n')
