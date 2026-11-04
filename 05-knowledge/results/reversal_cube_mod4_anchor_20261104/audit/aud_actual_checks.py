"""Referee checks of THM-4864 on ACTUAL named tournaments (exact Hamilton counts; no contracted model is
used to produce the tested values -- the cycle engine is only used for the Lemma 1 cross-check).

A1  Proposition 4 through actual blow-ups: weight-0 vertices realized as 2- or 3-vertex children with
    random interiors, target = Z^op[independent random same-size interiors]; L_v(actual half) == w(v).
    Also Lemma 1: actual half-contrast == const + delta_Z (cycle engine) on every battery query.
A2  Theorem A: random Q (m=3..6) with a singleton child z and other children of sizes 2-3 (variant: also
    further singletons), target interiors independent random; L_z(actual) == 1, witness query residue 2.
A3  Theorem B: TT_k child (k=2,3,4), u = b_1 or b_k, targets unchanged / reversed / random transitive with u
    at the same end / at the opposite end; other children random (hard case: sizes 2-3; non-hard: also
    singletons), other target interiors independent random.  Predicted functional values + witness.
A4  Ambient lift: ports copy the core's exterior incidences; half parity of every battery query unchanged.
A5  Restoration control: a second changed node N (target N^op, equal profile) elsewhere in the ambient.
A6  Generalisation probe (NOT claimed by the packet): A_j non-transitive with u a source/sink of A_j and A'_j.
A7  m = 1 control (missing hypothesis m >= 2 in Theorem B).
"""
import random
import sys
import time

from aud_lib import (arc, battery_queries, coef, compose, converse, expand_query, ham_start_end, half_contrast,
                     Lfun, rand_tour, run_jobs, subq, transitive_from_order)

T0 = time.time()
FAIL = []
COUNT = {}
MAXN = 14  # max old names (augmented order <= 16)


def tick(key, ok, info=None):
    c = COUNT.setdefault(key, [0, 0])
    c[0] += 1
    if not ok:
        c[1] += 1
        if len(FAIL) < 40:
            FAIL.append((key, info))


def rand_interior(s, rng):
    return rand_tour(s, rng)


def contrast_table(W, Wp, kblocks, qs):
    act = [expand_query(q, kblocks) for q in qs]
    r = run_jobs([('H', len(W), W, None, act), ('H', len(Wp), Wp, None, act)])
    half = {}
    for i, q in enumerate(qs):
        half[q] = half_contrast(r[0][i], r[1][i])
    return half, r


# ----------------------------------------------------------------------------- A1
def a1(rng):
    for n in (3, 4, 5, 6, 7, 8):
        for _ in range({3: 30, 4: 30, 5: 25, 6: 20, 7: 15, 8: 12}[n]):
            while True:
                w = rng.randrange(1 << n)
                sizes = [1 if (w >> v) & 1 else rng.choice([2, 2, 3]) for v in range(n)]
                if sum(sizes) <= 13:
                    break
                sizes = [1 if (w >> v) & 1 else 2 for v in range(n)]
                if sum(sizes) <= 13:
                    break
            Z = rand_tour(n, rng)
            src = [rand_interior(s, rng) for s in sizes]
            tgt = [rand_interior(s, rng) for s in sizes]
            W, blocks = compose(Z, src)
            Wp, blocks2 = compose(converse(Z), tgt)
            assert blocks == blocks2
            qs = sorted(battery_queries(n, range(n)))
            half, _ = contrast_table(W, Wp, blocks, qs)
            for v in range(n):
                tick('A1 actual blow-up: L_v(half) == w(v)  [Prop 4 + Lemma 1], n=%d' % n,
                     Lfun(half, Z, v) == ((w >> v) & 1), (Z, w, sizes, v))
            r = run_jobs([('C', n, Z, w, qs), ('C', n, converse(Z), w, qs)])
            diffs = {half[q] ^ r[0][i] ^ r[1][i] for i, q in enumerate(qs)}
            tick('A1 Lemma 1: actual half - delta_Z is query-independent', len(diffs) == 1, (Z, w, sizes))


# ----------------------------------------------------------------------------- A2
def a2(rng):
    for trial in range(80):
        m = rng.choice([3, 4, 5, 6])
        Q = rand_tour(m, rng)
        z = rng.randrange(m)
        extra_singletons = trial % 4 == 3
        while True:
            sizes = [rng.choice([1, 2, 3] if extra_singletons else [2, 3]) for _ in range(m)]
            sizes[z] = 1
            if sum(sizes) <= 13:
                break
        src = [rand_interior(s, rng) for s in sizes]
        tgt = [rand_interior(s, rng) for s in sizes]
        W, blocks = compose(Q, src)
        Wp, _ = compose(converse(Q), tgt)
        sing = [i for i in range(m) if sizes[i] == 1]
        qs = sorted(battery_queries(m, sing))
        half, raw = contrast_table(W, Wp, blocks, qs)
        key = 'A2 Theorem A: L_z(actual coarse half) == 1' + (' [several singletons]' if extra_singletons else '')
        for zz in sing:
            tick(key, Lfun(half, Q, zz) == 1, (Q, sizes, zz))
        wit = [q for q in qs if half[q] == 1]
        tick('A2 Theorem A: some coarse battery query has F_W-F_W\' = 2 mod 4', len(wit) > 0, (Q, sizes))
        tick('A2 battery size <= 4(m+1) for L_z', len(battery_queries(m, [z])) <= 4 * (m + 1), m)
        # Lemma 1 cross-check with the weighted model
        wq = sum(1 << i for i in range(m) if sizes[i] == 1)
        r = run_jobs([('C', m, Q, wq, qs), ('C', m, converse(Q), wq, qs)])
        tick('A2 Lemma 1: actual coarse half - delta_Q is query-independent',
             len({half[q] ^ r[0][i] ^ r[1][i] for i, q in enumerate(qs)}) == 1, (Q, sizes))


# ----------------------------------------------------------------------------- Theorem B instances
def pj(Q, j):
    S, E = ham_start_end(Q)
    return (S[j] + E[j]) % 2


def build_B(rng, m, k, hard, force_p0, u_end, nontransitive=False):
    while True:
        Q = rand_tour(m, rng)
        j = rng.randrange(m)
        if not force_p0 or m < 3 or pj(Q, j) == 0:
            break
    while True:
        sizes = [rng.choice([2, 3]) if hard else rng.choice([1, 2, 3]) for _ in range(m)]
        sizes[j] = k
        if not hard and not any(sizes[i] == 1 for i in range(m) if i != j):
            sizes[rng.choice([i for i in range(m) if i != j])] = 1
        if sum(sizes) <= MAXN:
            break
    # A_j on local indices 0..k-1 with a random physical order b_1..b_k
    order = list(range(k))
    rng.shuffle(order)
    u_loc = order[0] if u_end == 'b1' else order[-1]
    beta_loc = [x for x in range(k) if x != u_loc]
    if nontransitive:
        # u is the source (b1) / sink (bk) of A_j; beta interior random (any tournament)
        Aj = [0] * k
        sub = rand_tour(k - 1, rng)
        for a_, x in enumerate(beta_loc):
            for b_, y in enumerate(beta_loc):
                if sub[a_] >> b_ & 1:
                    Aj[x] |= 1 << y
        for x in beta_loc:
            if u_end == 'b1':
                Aj[u_loc] |= 1 << x
            else:
                Aj[x] |= 1 << u_loc
    else:
        Aj = transitive_from_order(order, k)
    src = [rand_interior(s, rng) for s in sizes]
    src[j] = Aj
    return Q, j, sizes, src, order, u_loc, beta_loc


def targets_for(rng, k, Aj, order, u_loc, beta_loc, u_end, nontransitive):
    T = {}
    if nontransitive:
        for typ in ('same', 'flip'):
            A = [0] * k
            sub = rand_tour(k - 1, rng)
            for a_, x in enumerate(beta_loc):
                for b_, y in enumerate(beta_loc):
                    if sub[a_] >> b_ & 1:
                        A[x] |= 1 << y
            u_src = (u_end == 'b1') == (typ == 'same')
            for x in beta_loc:
                if u_src:
                    A[u_loc] |= 1 << x
                else:
                    A[x] |= 1 << u_loc
            T['nontransitive-' + typ] = (A, typ == 'same')
        return T
    T['unchanged'] = (Aj, True)
    T['reversed'] = (converse(Aj), False)
    bo = list(beta_loc)
    rng.shuffle(bo)
    same_order = [u_loc] + bo if u_end == 'b1' else bo + [u_loc]
    T['random transitive, u same end'] = (transitive_from_order(same_order, k), True)
    bo2 = list(beta_loc)
    rng.shuffle(bo2)
    flip_order = bo2 + [u_loc] if u_end == 'b1' else [u_loc] + bo2
    T['random transitive, u other end'] = (transitive_from_order(flip_order, k), False)
    return T


def kernel_of(Q, j, sizes, blocks, u_loc, beta_loc, u_src):
    """kernel index: j = u, m = beta, others children.  Returns (Z, kblocks)."""
    m = len(Q)
    n = m + 1
    Z = [0] * n
    for i in range(m):
        for i2 in range(m):
            if i != i2 and arc(Q, i, i2):
                Z[i] |= 1 << i2
    for i in range(m):
        if i != j:
            if arc(Q, j, i):
                Z[m] |= 1 << i
            else:
                Z[i] |= 1 << m
    if u_src:
        Z[j] |= 1 << m
    else:
        Z[m] |= 1 << j
    kb = [list(b) for b in blocks] + [[blocks[j][x] for x in beta_loc]]
    kb[j] = [blocks[j][u_loc]]
    return Z, kb


def run_B_instance(rng, m, k, hard, force_p0, u_end, nontransitive=False, lift=None):
    Q, j, sizes, src, order, u_loc, beta_loc = build_B(rng, m, k, hard, force_p0, u_end, nontransitive)
    W, blocks = compose(Q, src)
    Zk, kb = kernel_of(Q, j, sizes, blocks, u_loc, beta_loc, u_end == 'b1')
    n = m + 1
    u, beta = j, m
    sing = [i for i in range(m) if i != j and sizes[i] == 1]
    vs = [u, beta] + sing
    qs = set(battery_queries(n, vs))
    qs |= set(subq((1 << u) | (1 << beta), (1 << u) | (1 << beta)))  # extra cubics + coarse a_j/c_j
    qs = sorted(qs)
    p = pj(Q, j)
    tg = targets_for(rng, k, src[j], order, u_loc, beta_loc, u_end, nontransitive)
    results = {}
    for name, (Ajt, same) in tg.items():
        tgt = [rand_interior(s, rng) for s in sizes]
        tgt[j] = Ajt
        Wp, b2 = compose(converse(Q), tgt)
        assert b2 == blocks
        half, raw = contrast_table(W, Wp, kb, qs)
        results[name] = (same, half, Wp)
        tag = '%s k=%s' % ('u=b1' if u_end == 'b1' else 'u=bk', '2' if k == 2 else '>=3')
        pre = 'A6 [generalisation] ' if nontransitive else 'A3 '
        if not same:
            tick(pre + 'reversed-type: L_u == 1 (%s)' % tag, Lfun(half, Zk, u) == 1, (Q, j, sizes, name))
        else:
            if not hard:
                for z in sing:
                    tick(pre + 'unchanged-type non-hard: L_z == 1 (%s)' % tag, Lfun(half, Zk, z) == 1,
                         (Q, j, sizes, name))
            elif k >= 3:
                lam = Lfun(half, Zk, u) ^ Lfun(half, Zk, beta) ^ coef(half, (1 << u) | (1 << beta), 1 << beta)
                tick(pre + "unchanged-type hard k>=3: Lambda' == 1+p_j(Q) (%s)" % tag, lam == 1 ^ p,
                     (Q, j, sizes, name, p))
                if p == 1:
                    # coarse a_j coefficient: query (u+beta) vs (0): Moebius over the 1-dim coarse cube
                    cj = half[((1 << u) | (1 << beta), 0)] ^ half[(0, 0)]
                    tick(pre + 'unchanged-type hard k>=3, p_j=1: coarse [a_j] == 1 (%s)' % tag, cj == 1,
                         (Q, j, sizes, name))
            else:
                tick(pre + 'unchanged-type hard k=2: L_u == 1 (%s)' % tag, Lfun(half, Zk, u) == 1,
                     (Q, j, sizes, name))
        if not same and not hard:
            for z in sing:
                tick(pre + 'reversed-type non-hard: L_z == 1 (%s)' % tag, Lfun(half, Zk, z) == 1, (Q, j, sizes, name))
        tick(pre + 'some query of the stated form has F_W-F_W\' = 2 mod 4', any(half[q] for q in qs),
             (Q, j, sizes, name))
        COUNT.setdefault('A3 p_j(Q) distribution (hard, unchanged-type, k>=3) p=0 / p=1', [0, 0])
        if hard and same and k >= 3 and not nontransitive:
            COUNT['A3 p_j(Q) distribution (hard, unchanged-type, k>=3) p=0 / p=1'][p] += 1
        if lift is not None and not nontransitive:
            lift(W, Wp, blocks, kb, qs, half, name, m)
    return Q, j, sizes


def make_ambient_lift(rng, restoration):
    def lift(W, Wp, blocks, kb, qs, half_core, name, m):
        N0 = len(W)
        e = rng.choice([1, 2]) if restoration else rng.choice([1, 2, 3])
        nN = 3 if restoration else 0
        if N0 + e + nN > MAXN:
            return
        # ambient quotient on {M, N?, e1..ee}: random tournament R; index 0 = M, 1 = N (if any)
        r = 1 + (1 if restoration else 0) + e
        R = rand_tour(r, rng)
        Nint = rand_tour(3, rng) if restoration else None
        parts_src = [W] + ([Nint] if restoration else []) + [[0]] * e
        parts_tgt = [Wp] + ([converse(Nint)] if restoration else []) + [[0]] * e
        T, bl = compose(R, parts_src)
        U, bl2 = compose(R, parts_tgt)
        assert bl == bl2
        core = bl[0]
        assert core == list(range(N0))
        ext = [p for p in range(len(T)) if p >= N0]
        # ports imitate M on the exterior: a_t = c_t = [t -> M]
        extmask = 0
        for t in ext:
            if arc(T, t, 0):  # t -> first core vertex == t -> M (module)
                extmask |= 1 << t
        act = []
        for q in qs:
            a, c = expand_query(q, kb)
            act.append((a | extmask, c | extmask))
        rr = run_jobs([('H', len(T), T, None, act), ('H', len(U), U, None, act)])
        ok = all(half_contrast(rr[0][i], rr[1][i]) == half_core[q] for i, q in enumerate(qs))
        key = ('A5 restoration (2nd changed node N -> N^op) + lift' if restoration else 'A4 ambient lift') + \
              ': half parity equals bare core at every battery query'
        tick(key, ok, (name, len(T)))
    return lift


def a3(rng):
    configs = []
    for m in (2, 3, 4, 5, 6):
        for k in (2, 3, 4):
            for hard in (True, False):
                if not hard and m < 2:
                    continue
                for u_end in ('b1', 'bk'):
                    configs.append((m, k, hard, u_end))
    for (m, k, hard, u_end) in configs:
        reps = 3 if m <= 4 else 2
        for r in range(reps):
            run_B_instance(rng, m, k, hard, force_p0=(r % 2 == 0), u_end=u_end)


def a4a5(rng):
    lift = make_ambient_lift(rng, restoration=False)
    rest = make_ambient_lift(rng, restoration=True)
    for t in range(14):
        m = rng.choice([2, 3, 4])
        k = rng.choice([2, 3])
        run_B_instance(rng, m, k, hard=rng.random() < 0.7, force_p0=True, u_end=rng.choice(['b1', 'bk']),
                       lift=lift if t % 2 == 0 else rest)


def a6(rng):
    for t in range(24):
        m = rng.choice([2, 3, 4, 5])
        k = rng.choice([3, 4])
        run_B_instance(rng, m, k, hard=rng.random() < 0.6, force_p0=(t % 2 == 0), u_end=rng.choice(['b1', 'bk']),
                       nontransitive=True)


def a7(rng):
    # m = 1: W = A_1 (TT_k), unchanged target W' = A_1 -> identical tournaments, no separation possible.
    for k in (2, 3, 4):
        A = transitive_from_order(list(range(k)), k)
        qs = [(a, c) for a in range(1 << k) for c in range(1 << k)]
        r = run_jobs([('H', k, A, None, qs), ('H', k, A, None, qs)])
        tick('A7 m=1, unchanged: F_W == F_W\' at every query (Theorem B needs m>=2)', r[0] == r[1], k)
        rr = run_jobs([('H', k, A, None, qs), ('H', k, converse(A), None, qs)])
        sep = any(half_contrast(rr[0][i], rr[1][i]) for i in range(len(qs)))
        print('A7 m=1 reversed k=%d: some named query separates mod 4: %s' % (k, sep))


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 4864
    rng = random.Random(seed)
    for name, fn in (('A7', a7), ('A1', a1), ('A2', a2), ('A3', a3), ('A4/A5', a4a5), ('A6', a6)):
        t = time.time()
        fn(rng)
        print('%s done (%.1fs)' % (name, time.time() - t), flush=True)
    print('\nseed', seed)
    for key in sorted(COUNT):
        if key.startswith('A3 p_j'):
            print('%-100s p=0: %d  p=1: %d' % (key, COUNT[key][0], COUNT[key][1]))
        else:
            print('%-100s cases=%5d  failures=%d' % (key, COUNT[key][0], COUNT[key][1]))
    print('FAILURES:', FAIL if FAIL else 'none')
    print('total time %.1fs' % (time.time() - T0))


if __name__ == '__main__':
    main()
