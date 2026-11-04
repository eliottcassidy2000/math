"""Referee model-level checks for THM-4864 (independent code; nothing imported from the packet).

M1  formula (1) == brute-force odd-weight port-cycle parity (all queries)
M2  Lemma 2 coefficient formulas vs brute-force Moebius coefficients
M3  converse duality g_{Z^op}(a,c) = g_Z(1+c,1+a)
M4  Lemma 3 endpoint parity + its consequence
M5  Proposition 4 (weighted E1) by brute-force cycles: exhaustive n=3,4; random n=5..9; n=2 boundary
M6  Theorem B kernel algebra in the weighted model: tau-orbit equality on ALL monomials,
    k=2 pure-monomial correction, functional predictions (u = b_1 and u = b_k, both extra-monomial forms)
"""
import random
import sys
import time

from aud_lib import (all_paths, all_tours, arc, battery_queries, coef, converse, ham_start_end,
                     Lfun, rand_tour, run_jobs, wt)

T0 = time.time()
FAIL = []
COUNT = {}


def tick(key, ok, info=None):
    c = COUNT.setdefault(key, [0, 0])
    c[0] += 1
    if not ok:
        c[1] += 1
        if len(FAIL) < 30:
            FAIL.append((key, info))


def all_queries(n):
    return [(a, c) for c in range(1 << n) for a in range(1 << n)]


# ----------------------------------------------------------------------------- formula (1)
def formula1_table(Z, w, qs):
    """Evaluate the packet's formula (1) at every query in qs (as functions, literals 1+a = 1-a)."""
    P = all_paths(Z)
    mask = {p: sum(1 << v for v in p) for p in P}
    single = {}   # (kind, s, e) -> parity
    for p in P:
        s, e = p[0], p[-1]
        if wt(p, w) == 0:
            single[('xa', s, e)] = single.get(('xa', s, e), 0) ^ 1
            single[('yc', s, e)] = single.get(('yc', s, e), 0) ^ 1
        else:
            single[('xy', s, e)] = single.get(('xy', s, e), 0) ^ 1
    pair = {}
    for p in P:
        for q in P:
            if mask[p] & mask[q]:
                continue
            if (wt(p, w) + wt(q, w)) % 2 == 1:
                key = (p[0], p[-1], q[0], q[-1])
                pair[key] = pair.get(key, 0) ^ 1
    out = []
    for a, c in qs:
        A = lambda v: (a >> v) & 1
        C = lambda v: (c >> v) & 1
        t = 0
        for (kind, s, e), par in single.items():
            if not par:
                continue
            if kind == 'xa':
                t ^= (1 - A(s)) & A(e)
            elif kind == 'yc':
                t ^= (1 - C(s)) & C(e)
            else:
                t ^= (1 - C(s)) & A(e)
        for (s, e, s2, e2), par in pair.items():
            if par:
                t ^= (1 - A(s)) & C(e) & (1 - C(s2)) & A(e2)
        out.append(t)
    return out


def m1(rng):
    jobs, meta = [], []
    for n in (1, 2, 3, 4):
        for Z in all_tours(n):
            for w in range(1 << n):
                jobs.append(('C', n, Z, w, all_queries(n)))
                meta.append(('M1 formula(1)==bruteforce, exhaustive n=%d (all queries)' % n, Z, w))
    for _ in range(40):
        n = 5
        Z = rand_tour(n, rng)
        w = rng.randrange(1 << n)
        jobs.append(('C', n, Z, w, all_queries(n)))
        meta.append(('M1 formula(1)==bruteforce, random n=5 (all 1024 queries)', Z, w))
    res = run_jobs(jobs)
    for (mode, n, Z, w, qs), r, (key, _, _) in zip(jobs, res, meta):
        tick(key, formula1_table(Z, w, qs) == r, (Z, w))


# ----------------------------------------------------------------------------- Lemma 2
def m2(rng):
    jobs, data = [], []
    for n in (3, 4, 5, 6):
        reps = {3: 60, 4: 60, 5: 40, 6: 6}[n]
        for _ in range(reps):
            Z = rand_tour(n, rng)
            w = rng.randrange(1 << n)
            qs = sorted(battery_queries(n, range(n)))
            jobs.append(('C', n, Z, w, qs))
            data.append((n, Z, w, qs))
    res = run_jobs(jobs)
    for (n, Z, w, qs), r in zip(data, res):
        g = dict(zip(qs, r))
        P = all_paths(Z)
        mask = {p: sum(1 << v for v in p) for p in P}
        N = {}
        for p in P:
            if len(p) >= 2 and wt(p, w) == 1:
                N[(p[0], p[-1])] = N.get((p[0], p[-1]), 0) ^ 1
        Pi = {}
        for p in P:
            for q in P:
                if mask[p] & mask[q]:
                    continue
                if (wt(p, w) + wt(q, w)) % 2 == 1:
                    Pi[(p[-1], q[-1])] = Pi.get((p[-1], q[-1]), 0) ^ 1
        for v in range(n):
            for u in range(n):
                if u == v:
                    continue
                pred = N.get((v, u), 0) ^ Pi.get((v, u), 0)
                tick('M2 [a_u c_v]g_Z = N_odd(v~>u)+Pi(->v,->u)', coef(g, 1 << u, 1 << v) == pred, (Z, w, u, v))
                cub = 0
                for p in P:
                    if p[-1] == u and not (mask[p] >> v & 1) and (wt(p, w) + ((w >> v) & 1)) % 2 == 1:
                        cub ^= 1
                tick('M2 [a_u a_v c_v]g_Z = #{P in Z-v ending at u, w(P)+w(v) odd}',
                     coef(g, (1 << u) | (1 << v), 1 << v) == cub, (Z, w, u, v))


# ----------------------------------------------------------------------------- duality
def m3(rng):
    jobs, data = [], []
    for n in (2, 3, 4, 5):
        for _ in range(30):
            Z = rand_tour(n, rng)
            w = rng.randrange(1 << n)
            qs = all_queries(n)
            jobs += [('C', n, Z, w, qs), ('C', n, converse(Z), w, qs)]
            data.append((n, qs))
    res = run_jobs(jobs)
    for i, (n, qs) in enumerate(data):
        g = dict(zip(qs, res[2 * i]))
        gop = dict(zip(qs, res[2 * i + 1]))
        full = (1 << n) - 1
        tick('M3 g_{Z^op}(a,c)=g_Z(1+c,1+a)', all(gop[(a, c)] == g[(full ^ c, full ^ a)] for (a, c) in qs))


# ----------------------------------------------------------------------------- Lemma 3
def induced(R, idx):
    sub = [0] * len(idx)
    for a_, i in enumerate(idx):
        for b_, j in enumerate(idx):
            if R[i] >> j & 1:
                sub[a_] |= 1 << b_
    return sub


def m4(rng):
    for n in range(1, 9):
        for _ in range({8: 8, 7: 15}.get(n, 30)):
            R = rand_tour(n, rng)
            S, E = ham_start_end(R)
            P = all_paths(R)
            for v in range(n):
                ends = sum(1 for p in P if p[-1] == v)
                starts = sum(1 for p in P if p[0] == v)
                tick('M4 #paths ending at v == S_v (mod 2)', ends % 2 == S[v] % 2, (R, v))
                tick('M4 #paths starting at v == E_v (mod 2)', starts % 2 == E[v] % 2, (R, v))
            if n <= 6:
                cache = {}
                for Smask in range(1, 1 << n):
                    idx = [i for i in range(n) if Smask >> i & 1]
                    cache[Smask] = (idx, ham_start_end(induced(R, idx)))
                for v in range(n):
                    tot = 0
                    for Smask in range(1, (1 << n) - 1):
                        if Smask >> v & 1:
                            idx, (Ss, Es) = cache[Smask]
                            k = idx.index(v)
                            tot += Ss[k] + Es[k]
                    tick('M4 sum_{v in S proper} eps_v(S) == 0', tot % 2 == 0, (R, v))


# ----------------------------------------------------------------------------- Proposition 4
def m5(rng):
    inst = []
    for n in (3, 4):
        for Z in all_tours(n):
            for w in range(1 << n):
                inst.append((Z, w, 'M5 weighted E1 L_v=w(v), exhaustive n=%d (all Z, all w, all v)' % n))
    for n in (5, 6, 7, 8, 9):
        reps = {5: 600, 6: 300, 7: 200, 8: 120, 9: 60}[n]
        for _ in range(reps):
            Z = rand_tour(n, rng)
            inst.append((Z, rng.randrange(1 << n), 'M5 weighted E1 L_v=w(v), random n=%d (all v)' % n))
    jobs, qss = [], []
    for Z, w, key in inst:
        n = len(Z)
        qs = sorted(battery_queries(n, range(n)))
        qss.append(qs)
        jobs += [('C', n, Z, w, qs), ('C', n, converse(Z), w, qs)]
    res = run_jobs(jobs)
    for i, (Z, w, key) in enumerate(inst):
        qs = qss[i]
        d = {q: res[2 * i][t] ^ res[2 * i + 1][t] for t, q in enumerate(qs)}
        for v in range(len(Z)):
            tick(key, Lfun(d, Z, v) == ((w >> v) & 1), (Z, w, v))
    # n = 2 boundary
    fails = tot = 0
    for Z in all_tours(2):
        for w in range(4):
            qs = sorted(battery_queries(2, range(2)))
            r = run_jobs([('C', 2, Z, w, qs), ('C', 2, converse(Z), w, qs)])
            d = {q: r[0][t] ^ r[1][t] for t, q in enumerate(qs)}
            for v in range(2):
                tot += 1
                u = 1 - v
                pred = ((w >> u) & 1) ^ ((w >> v) & 1)  # g(v) at n=2
                tick('M5 n=2 boundary: L_v == g(v) = w(u)+w(v)', Lfun(d, Z, v) == pred, (Z, w, v))
                fails += Lfun(d, Z, v) != ((w >> v) & 1)
    print('M5 n=2: L_v != w(v) in %d of %d cases (boundary is real)' % (fails, tot))


# ----------------------------------------------------------------------------- Theorem B kernel algebra (model)
def moebius(vals, nv):
    f = list(vals)
    for i in range(nv):
        bit = 1 << i
        for x in range(1 << nv):
            if x & bit:
                f[x] ^= f[x ^ bit]
    return f


def kernel(Q, j, u_is_source, wts_other, k):
    """Kernel on m+1 vertices: index j = u, index m = beta, others = children.  Returns (Z, w)."""
    m = len(Q)
    n = m + 1
    beta = m
    Z = [0] * n
    for i in range(m):
        for i2 in range(m):
            if i != i2 and arc(Q, i, i2):
                Z[i] |= 1 << i2
    for i in range(m):
        if i != j:
            if arc(Q, j, i):
                Z[beta] |= 1 << i
            else:
                Z[i] |= 1 << beta
    if u_is_source:
        Z[j] |= 1 << beta
    else:
        Z[beta] |= 1 << j
    w = 1 << j
    if k == 2:
        w |= 1 << beta
    for i in range(m):
        if i != j and wts_other[i]:
            w |= 1 << i
    return Z, w


def flip_ub(Z, j, beta):
    Z2 = list(Z)
    if Z[j] >> beta & 1:
        Z2[j] &= ~(1 << beta)
        Z2[beta] |= 1 << j
    else:
        Z2[beta] &= ~(1 << j)
        Z2[j] |= 1 << beta
    return Z2


def m6(rng):
    trials = []
    for trial in range(500):
        m = rng.choice([2, 3, 4, 5])
        Q = rand_tour(m, rng)
        j = rng.randrange(m)
        k = rng.choice([2, 3])  # only w_beta = [k=2] enters the kernel
        u_src = rng.random() < 0.5
        wo = [0] * m if rng.random() < 0.6 else [rng.randrange(2) for _ in range(m)]
        wo[j] = 0
        trials.append((m, Q, j, k, u_src, wo))
    jobs, data = [], []
    for (m, Q, j, k, u_src, wo) in trials:
        Z, w = kernel(Q, j, u_src, wo, k)
        n = m + 1
        Zop = converse(Z)
        Zu = flip_ub(Zop, j, m)  # Q^op with the same u-beta arc as Z (unchanged type)
        qs = [(x & ((1 << n) - 1), x >> n) for x in range(1 << (2 * n))]
        jobs += [('C', n, Z, w, qs), ('C', n, Zop, w, qs), ('C', n, Zu, w, qs)]
        data.append((Z, w, n))
    res = run_jobs(jobs)
    for t, (m, Q, j, k, u_src, wo) in enumerate(trials):
        Z, w, n = data[t]
        gZ, gZop, gZu = res[3 * t], res[3 * t + 1], res[3 * t + 2]
        nv = 2 * n
        u, beta = j, m
        mD = moebius([gZ[i] ^ gZop[i] for i in range(1 << nv)], nv)
        mH = moebius([gZ[i] ^ gZu[i] for i in range(1 << nv)], nv)

        def tau(x):
            def sw(y):
                bu, bb = (y >> u) & 1, (y >> beta) & 1
                y &= ~((1 << u) | (1 << beta))
                return y | (bb << u) | (bu << beta)
            return sw(x & ((1 << n) - 1)) | (sw(x >> n) << n)
        ok = True
        for x in range(1 << nv):
            tx = tau(x)
            if tx < x:
                continue
            if (mD[x] ^ (mD[tx] if tx != x else 0)) != (mH[x] ^ (mH[tx] if tx != x else 0)):
                ok = False
                break
        tick('M6 tau-orbit sums: unchanged half == delta_Z on EVERY orbit', ok, (Q, j, k, u_src, wo))
        hard = all(wo[i] == 0 for i in range(m) if i != j)
        if k == 2 and hard:
            corr = moebius([gZop[i] ^ gZu[i] for i in range(1 << nv)], nv)
            mixed = any(corr[x] for x in range(1 << nv) if (x & ((1 << n) - 1)) and (x >> n))
            tick('M6 k=2 hard: g_{Z^op}-g_{Z_u} has no mixed monomial', not mixed, (Q, j, u_src))

        def cm(mob, A, C):
            return mob[A | (C << n)]

        def L(mob, v):
            s = 0
            for t2 in range(n):
                if t2 == v:
                    continue
                s ^= cm(mob, 1 << t2, 1 << v)
                if Z[v] >> t2 & 1:
                    s ^= cm(mob, (1 << t2) | (1 << v), 1 << v)
            return s
        tick('M6 reversed: L_u = 1', L(mD, u) == 1, (Q, j, k, u_src, wo))
        S, E = ham_start_end(Q)
        p = (S[j] + E[j]) % 2
        if hard:
            if k >= 3:
                lit = L(mH, u) ^ L(mH, beta) ^ cm(mH, (1 << u) | (1 << beta), 1 << beta)
                mp = L(mH, u) ^ L(mH, beta) ^ (cm(mH, (1 << u) | (1 << beta), 1 << beta) if u_src
                                               else cm(mH, (1 << u) | (1 << beta), 1 << u))
                tag = 'u=b_1' if u_src else 'u=b_k'
                tick("M6 unchanged k>=3 hard (%s): Lambda' literal [a_u a_b c_b] = 1+p_j" % tag, lit == 1 ^ p,
                     (Q, j, u_src, p))
                tick("M6 unchanged k>=3 hard (%s): Lambda' with missing tau-partner = 1+p_j" % tag, mp == 1 ^ p,
                     (Q, j, u_src, p))
                # p_j = 1 branch: coarse [a_j] of the fully contracted model must be 1
                if p == 1:
                    Qop = converse(Q)
                    wq = 0
                    qsq = [(0, 0), (1 << j, 0)]
                    r = run_jobs([('C', m, Q, wq, qsq), ('C', m, Qop, wq, qsq)])
                    cj = (r[0][0] ^ r[1][0]) ^ (r[0][1] ^ r[1][1])
                    tick('M6 hard, p_j=1: coarse [a_j] delta_Q = 1', cj == 1, (Q, j))
            else:
                tick('M6 unchanged k=2 hard: L_u = 1', L(mH, u) == 1, (Q, j, u_src))
        else:
            for z in range(m):
                if z != j and wo[z]:
                    tick('M6 non-hard: L_z = 1 (unchanged target)', L(mH, z) == 1, (Q, j, z))
                    tick('M6 non-hard: L_z = 1 (reversed target)', L(mD, z) == 1, (Q, j, z))


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20261104
    rng = random.Random(seed)
    for name, fn in (('M1', m1), ('M2', m2), ('M3', m3), ('M4', m4), ('M5', m5), ('M6', m6)):
        t = time.time()
        fn(rng)
        print('%s done (%.1fs)' % (name, time.time() - t), flush=True)
    print('\nseed', seed)
    for k in sorted(COUNT):
        print('%-92s cases=%6d  failures=%d' % (k, COUNT[k][0], COUNT[k][1]))
    print('FAILURES:', FAIL if FAIL else 'none')
    print('total time %.1fs' % (time.time() - T0))


if __name__ == '__main__':
    main()
