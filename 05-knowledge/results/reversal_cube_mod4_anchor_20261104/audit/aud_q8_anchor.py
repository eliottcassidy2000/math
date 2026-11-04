"""Theorem B on the S188 prime-eight endpoint-obstruction quotient Q8 (p(Q8)=0, coarse-blind with equal
interiors) on ACTUAL graphs: TT_k child at j (u = b_1), other children TT2 (source) and independent random
2-vertex interiors (target).  Predictions: reversed L_u = 1; unchanged k=2: L_u = 1; unchanged k=3: Lambda' = 1."""
import random, sys, time
from aud_lib import (battery_queries, coef, compose, converse, expand_query, half_contrast, Lfun, rand_tour,
                     run_jobs, subq, transitive_from_order)
from aud_actual_checks import kernel_of
Q8 = [238, 4, 240, 118, 195, 146, 34, 74]
rng = random.Random(8)
m = 8
fails = 0
t0 = time.time()
for k, js in ((2, range(8)), (3, (0, 3, 5, 7))):
    for j in js:
        sizes = [2] * m
        sizes[j] = k
        Aj = transitive_from_order(list(range(k)), k)   # b_1 = local 0
        src = [rand_tour(2, rng) for _ in range(m)]
        src[j] = Aj
        W, blocks = compose(Q8, src)
        Z, kb = kernel_of(Q8, j, sizes, blocks, 0, list(range(1, k)), True)
        u, beta, n = j, m, m + 1
        for typ, Ajt in (('unchanged', Aj), ('reversed', converse(Aj))):
            tgt = [rand_tour(2, rng) for _ in range(m)]
            tgt[j] = Ajt
            Wp, _ = compose(converse(Q8), tgt)
            vs = [u] if (typ == 'reversed' or k == 2) else [u, beta]
            qs = sorted(set(battery_queries(n, vs)) | set(subq((1 << u) | (1 << beta), 1 << beta)))
            act = [expand_query(q, kb) for q in qs]
            r = run_jobs([('H', len(W), W, None, act), ('H', len(Wp), Wp, None, act)])
            half = {q: half_contrast(r[0][i], r[1][i]) for i, q in enumerate(qs)}
            if typ == 'reversed' or k == 2:
                val, name = Lfun(half, Z, u), 'L_u'
            else:
                val = Lfun(half, Z, u) ^ Lfun(half, Z, beta) ^ coef(half, (1 << u) | (1 << beta), 1 << beta)
                name = "Lambda'"
            ok = val == 1
            fails += not ok
            print('k=%d j=%d %-9s %s=%d (pred 1) witness=%s  [%d queries, %.0fs]' % (
                k, j, typ, name, val, any(half.values()), len(qs), time.time() - t0), flush=True)
print('FAILURES:', fails)
