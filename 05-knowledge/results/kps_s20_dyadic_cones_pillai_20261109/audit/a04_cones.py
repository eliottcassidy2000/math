# a04: THM-5109 (C) checks by brute force.
#  (a) Q(X => Y) = Q(X) Q(Y) on random pairs; (b) the source over X has triple (H(X), 0, Q(X));
#  (c) explicit cone realisations H(X) + 2^K for X = C4, C4=>C4, C4=>C19, C19=>C4, C217 (codes from the certified lists),
#      built as one-arc joins of (source over X, source) with TT_m at entry p = K - j + 1, and counted by brute force;
#  (d) e_v = O - E (alternating path sum) on random rooted tournaments (the identity the census inputs used).
import random
from a_lib import *
rng = random.Random(2026)
bad = 0; cnt = 0
for _ in range(60):
    X = rand_tour(rng.randrange(0, 6), rng); Y = rand_tour(rng.randrange(0, 6), rng)
    r = run([X, Y, ordinal_sum(X, Y)])
    cnt += 1
    if r[0]['Q'] * r[1]['Q'] != r[2]['Q']: bad += 1
print('(a) Q multiplicative over ordinal sums: %d random pairs, failures %d' % (cnt, bad))
bad = 0; cnt = 0
for _ in range(60):
    X = rand_tour(rng.randrange(0, 8), rng)
    T = ordinal_sum(TT(1), X)
    rX, rT = run([X, T])
    cnt += 1
    if (rT['H'], rT['e'][0], rT['R'][0]) != (rX['H'], 0 if len(X) else 1, rX['Q']): bad += 1
print('(b) source over X has (h, e, R) = (H(X), 0, Q(X)) (X nonempty): %d random X, failures %d' % (cnt, bad))

C4 = from_code(34, 4); C19 = from_code(17057, 6); C19b = from_code(6658, 6); C217 = from_code(144179468, 8)
C223 = from_code(147063080, 8); C599 = from_code(72166211, 8)
tests = [('C4', C4), ('C4=>C4', ordinal_sum(C4, C4)), ('C4=>C19', ordinal_sum(C4, C19)),
         ('C19b=>C4', ordinal_sum(C19b, C4)), ('C217', C217), ('C223', C223), ('C599', C599),
         ('TT2=>C4=>TT1', ordinal_sum(ordinal_sum(TT(2), C4), TT(1)))]
for name, X in tests:
    T = ordinal_sum(TT(1), X); rT = run([T])[0]
    h, e, R = rT['H'], rT['e'][0], rT['R'][0]
    assert ispow2(R), (name, R)
    j = R.bit_length() - 1
    outs = []
    for K in range(j + 1, j + 4):
        p = K - j + 1
        if len(T) + p > 16: break
        N = one_arc_join(T, 0, TT(p), p - 1)
        HN = run([N])[0]['H']
        outs.append((K, HN, HN == h + 2 ** K))
    print('(c) %-13s H=%d Q=2^%d : ' % (name, h, j) + ', '.join('K=%d: H(N)=%d %s' % (K, HN, 'ok' if ok else 'FAIL') for K, HN, ok in outs))
# (d) e = O - E
bad = 0; cnt = 0
for _ in range(25):
    n = rng.randrange(1, 8); A = rand_tour(n, rng)
    r = run([A])[0]
    for v in range(n):
        alt = 0
        for P in paths_from(A, v, list(range(n))):
            rest = [w for w in range(n) if w not in P]
            alt += (1 if len(P) % 2 else -1) * naive_H(A, rest)
        cnt += 1
        if alt != r['e'][v]: bad += 1
print('(d) e_v = O - E: %d ports, failures %d' % (cnt, bad))
