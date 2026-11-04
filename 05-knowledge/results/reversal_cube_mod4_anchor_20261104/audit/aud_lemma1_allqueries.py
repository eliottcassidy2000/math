"""Lemma 1 on ALL child-uniform queries (not just the low-degree battery): actual (F_W-F_W')/2 mod 2 minus the
weighted model delta_Z = g_Z - g_{Z^op} is one constant over all 2^(2n) coarse queries.  Random kernels n=3,4,5
(weights = [child size 1]), source/target interiors independent random, sizes 1-3; also anchored kernels
(Theorem B shape: u named, beta = A_j - u contracted) with all 2^(2(m+1)) anchored queries."""
import random
from aud_lib import compose, converse, expand_query, half_contrast, rand_tour, run_jobs, transitive_from_order
from aud_actual_checks import kernel_of
rng = random.Random(41)
cases = fails = 0
for trial in range(24):
    n = rng.choice([3, 4, 5])
    Z = rand_tour(n, rng)
    while True:
        sizes = [rng.choice([1, 2, 3]) for _ in range(n)]
        if sum(sizes) <= 11:
            break
    w = sum(1 << i for i in range(n) if sizes[i] == 1)
    W, blocks = compose(Z, [rand_tour(s, rng) for s in sizes])
    Wp, _ = compose(converse(Z), [rand_tour(s, rng) for s in sizes])
    qs = [(x & ((1 << n) - 1), x >> n) for x in range(1 << (2 * n))]
    act = [expand_query(q, blocks) for q in qs]
    r = run_jobs([('H', len(W), W, None, act), ('H', len(Wp), Wp, None, act),
                  ('C', n, Z, w, qs), ('C', n, converse(Z), w, qs)])
    vals = {half_contrast(r[0][i], r[1][i]) ^ r[2][i] ^ r[3][i] for i in range(len(qs))}
    cases += 1
    fails += len(vals) != 1
print('coarse: %d instances (all 2^(2n) queries each), failures %d' % (cases, fails))
cases2 = fails2 = 0
for trial in range(16):
    m = rng.choice([2, 3])
    Q = rand_tour(m, rng)
    j = rng.randrange(m)
    k = rng.choice([2, 3, 4])
    while True:
        sizes = [rng.choice([1, 2, 3]) for _ in range(m)]
        sizes[j] = k
        if sum(sizes) <= 11:
            break
    src = [rand_tour(s, rng) for s in sizes]
    src[j] = transitive_from_order(list(range(k)), k)
    W, blocks = compose(Q, src)
    Zk, kb = kernel_of(Q, j, sizes, blocks, 0, list(range(1, k)), True)
    n = m + 1
    w = (1 << j) | ((1 << m) if k == 2 else 0) | sum(1 << i for i in range(m) if i != j and sizes[i] == 1)
    for typ in ('unchanged', 'reversed'):
        tgt = [rand_tour(s, rng) for s in sizes]
        tgt[j] = src[j] if typ == 'unchanged' else converse(src[j])
        Wp, _ = compose(converse(Q), tgt)
        # target kernel: Q^op with u-beta arc as in A'_j
        Zt = converse(Zk)
        if typ == 'unchanged':  # restore u->beta
            Zt[m] &= ~(1 << j); Zt[j] |= 1 << m
        qs = [(x & ((1 << n) - 1), x >> n) for x in range(1 << (2 * n))]
        act = [expand_query(q, kb) for q in qs]
        r = run_jobs([('H', len(W), W, None, act), ('H', len(Wp), Wp, None, act),
                      ('C', n, Zk, w, qs), ('C', n, Zt, w, qs)])
        vals = {half_contrast(r[0][i], r[1][i]) ^ r[2][i] ^ r[3][i] for i in range(len(qs))}
        cases2 += 1
        fails2 += len(vals) != 1
print('anchored kernel (source AND target contraction): %d instances (all 2^(2(m+1)) queries), failures %d' % (cases2, fails2))
