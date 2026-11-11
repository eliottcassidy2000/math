"""verify_s22.py - master verification of THM-5130 (kind-pasteur-2026-11-10-S22). Pure Python, no input files.
Every block recounts paths literally (subset dynamic programme of towers.py, or the transfer-matrix count of
tower_build.py, which shares no code with it) and compares with the stated formulas.  Prints PASS / FAIL per block."""
import itertools, random
from towers import tables, H, state, gadget, apply, level_up, TT, rand_t
from tower_build import build, adjacency, count_apex, start_state, step, tower_state
rng = random.Random(20261110)
allok = True


def report(name, ok, extra=""):
    global allok
    allok &= bool(ok)
    print(("PASS  " if ok else "FAIL  ") + name + ("  " + extra if extra else ""))


def adj(m, code):
    A = [[0] * m for _ in range(m)]
    b = 0
    for i in range(m):
        for j in range(i + 1, m):
            if code >> b & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            b += 1
    return A


def labelled(m):
    k = m * (m - 1) // 2
    for code in range(1 << k):
        yield adj(m, code)


def Qmat(A):
    n = len(A)
    full = (1 << n) - 1
    en, st, hs = tables(A)
    Q = [[0] * n for _ in range(n)]
    for X in range(1, full):
        Y = full ^ X
        for e in range(n):
            if X >> e & 1 and en[X][e]:
                for f in range(n):
                    if Y >> f & 1:
                        Q[e][f] += en[X][e] * st[Y][f]
    U = [sum(st[X][o] * hs[full ^ X] for X in range(1, full + 1) if X >> o & 1) for o in range(n)]
    return hs[full], Q, U


def components(A):
    n = len(A)
    reach = [[i == j or bool(A[i][j]) for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
    comp = []
    seen = set()
    for i in range(n):
        if i in seen:
            continue
        c = [j for j in range(n) if reach[i][j] and reach[j][i]]
        seen |= set(c)
        comp.append(c)
    return comp, reach


# V1 closing law
bad = 0
cnt = 0
for _ in range(80):
    n = rng.randint(1, 5)
    A = rand_t(n, rng)
    O = [i for i in range(n) if rng.random() < 0.5]
    r, M, h, E = state(A, O)
    for t in (1, 2, 3):
        cnt += 1
        bad += (H(level_up(A, O, TT(t))) != r + 2 ** t * M)
report("V1 closing law H(B => TT_t => v, v -> O) = r + 2^t M", bad == 0, "%d cases" % cnt)

# V2 gadget law
bad = 0
cnt = 0
for _ in range(300):
    n = rng.randint(1, 4)
    A = rand_t(n, rng)
    O = [i for i in range(n) if rng.random() < 0.6]
    k = rng.randint(1, 4)
    AX = rand_t(k, rng)
    P = [i for i in range(k) if rng.random() < 0.6]
    cnt += 1
    bad += (state(level_up(A, O, AX), [n + p for p in P]) != apply(gadget(AX, P)["G"], state(A, O)))
report("V2 gadget law s(B+, P) = G(X, P) s(B, O), random bases and gadgets", bad == 0, "%d cases" % cnt)

# V3 transitive pairs: constants for every distance
ok = True
for delta in (1, 2, 3, 4, 5):
    for i in range(1, 4):
        for m in range(0, 3):
            a = i + delta + m
            g = gadget(TT(a), [i - 1, i - 1 + delta])
            rX = -1 if delta == 1 else 1 - 2 ** (delta - 1)
            kap = 2 ** (a - delta) * (2 ** delta - 1 - 3 ** (delta - 1))
            G = [[rX, kap, 0, 0],
                 [(2 ** delta + 1) * 2 ** m, 2 * (3 ** delta + 2 ** delta) * 2 ** (i - 1) * 3 ** m,
                  (2 ** delta + 1) * 2 ** m, -(1 if i == 1 else 0)],
                 [1, 2 ** a, 0, 0],
                 [(1 if m == 0 else 0), (2 ** delta + 1) * 2 ** (i - 1), -(1 if m == 0 else 0), 0]]
            ok &= (g["G"] == G)
report("V3 transitive pair at distance delta <= 5: r_X, kappa_X = 2^(a-delta)(2^delta - 1 - 3^(delta-1)), whole matrix", ok)

# V4 sink and source
bad = 0
for _ in range(60):
    n = rng.randint(1, 5)
    A = rand_t(n, rng)
    O = [i for i in range(n) if rng.random() < 0.5]
    r, M, h, E = state(A, O)
    S = [row + [1] for row in A] + [[0] * (n + 1)]
    bad += (state(S, O) != (r, 2 * M, h, 0))
    R = [[0] + [1] * n] + [[0] + row for row in A]
    bad += (state(R, [o + 1 for o in O]) != (r, M, h, E))
report("V4 a sink doubles M and kills E; a source changes nothing", bad == 0)

# V5 structure lemmas on all labelled tournaments of order <= 5
g1 = g2 = g3 = g4 = True
unit = 0
for m in range(1, 6):
    for A in labelled(m):
        h, Q, U = Qmat(A)
        comp, reach = components(A)
        cid = {v: k for k, c in enumerate(comp) for v in c}
        for o in range(m):
            g1 &= (U[o] == h + sum(Q[x][o] for x in range(m) if A[o][x]))
        for e in range(m):
            for f in range(e + 1, m):
                if cid[e] != cid[f]:
                    g2 &= (Q[e][f] + Q[f][e] >= h)
                B2 = [row[:] for row in A]
                B2[e][f], B2[f][e] = A[f][e], A[e][f]
                h2, Q2, U2 = Qmat(B2)
                g4 &= (Q2[e][f] == Q[e][f] and Q2[f][e] == Q[f][e])
        for Om in range(3, 1 << m):
            O = [i for i in range(m) if Om >> i & 1]
            if len(O) < 2:
                continue
            w = sum(Q[e][f] for e in O for f in O if e != f)
            if h - w == 1:
                unit += 1
                cs = {cid[o] for o in O}
                g3 &= (len(cs) == 1 and all(len(c) == 1 for k, c in enumerate(comp) if k not in cs))
                if len(cs) == 1:
                    Qc = comp[cs.pop()]
                    idx = {v: k for k, v in enumerate(Qc)}
                    AQ = [[A[x][y] for y in Qc] for x in Qc]
                    rq, Mq, hq, Eq = state(AQ, [idx[o] for o in O])
                    after = sum(1 for c in comp if c[0] not in Qc and reach[Qc[0]][c[0]])
                    g3 &= (rq == 1 and sum(U[o] for o in O) == 2 ** after * Mq)
report("V5a U_B(o) = H(B) + sum over out-neighbours x of Q_B(x, o)  (all labelled tournaments of order <= 5)", g1)
report("V5b W_ef >= H(B) for e, f in different strong components", g2)
report("V5c unit sets with |O| >= 2: one strong component, all others single vertices, M = 2^(sinks after) M_core", g3,
       "%d unit sets" % unit)
report("V5d W_ef does not depend on the orientation of the arc ef (so r(B) - r(B') = H(B) - H(B'))", g4)
bad = 0
for _ in range(200):
    n = rng.randint(1, 5)
    K = rand_t(n, rng)
    y = rng.randrange(n)
    D = [[0] * (n + 1) for _ in range(n + 1)]
    for i in range(n):
        for j in range(n):
            D[i][j] = K[i][j]
        D[i][n] = K[i][y]
        D[n][i] = K[y][i]
    D[y][n] = 1
    D[n][y] = 0
    r, M, h, E = state(D, [y, n])
    bad += (r != H(D) - 2 * H(K))
report("V5e twin lemma r(K[TT_2 at y], twins) = H(K[TT_2 at y]) - 2 H(K)", bad == 0)

# V6 explicit families by the transfer-matrix count
ok = True
cnt = 0
for d in (1, 2):
    for i in range(1, 4):
        for m in range(0, 3):
            for t in (1, 2, 3):
                N, q, beats = build((i + d + m, [i, i + d]), [], t)
                cnt += 1
                ok &= (count_apex(N, q, adjacency(N, q, beats)) == (2 ** d + 1) * 2 ** (m + t) - 1)
report("V6a one closer on a transitive pair: H = (2^delta + 1) 2^(m+t) - 1", ok, "%d tournaments" % cnt)
ok = True
cnt = 0
for n in (1, 2, 3):
    for p in range(1, n + 1):
        for d in (1, 2):
            for i in (1, 2, 3):
                for m in (0, 1, 2):
                    for t in (1, 2):
                        M1 = (2 ** d + 1) * 2 ** (m + 1) + 2 * (3 ** d + 2 ** d) * 2 ** (i - 1 + n - p) * 3 ** m \
                            - (1 if (i == 1 and p == n) else 0)
                        N, q, beats = build((n, [p]), [(d, i, m)], t)
                        cnt += 1
                        ok &= (count_apex(N, q, adjacency(N, q, beats)) == 2 ** t * M1 - 1)
report("V6b two closers from a fan: H = 2^t [(2^d+1) 2^(m+1) + 2 (3^d+2^d) 2^(i-1+n-p) 3^m - [i=1][p=n]] - 1", ok,
       "%d tournaments" % cnt)
ok = True
cnt = 0
for d1 in (1, 2):
    for i1 in (1, 2):
        for m1 in (0, 1, 2):
            for d2 in (1, 2):
                for i2 in (1, 2, 3):
                    for m2 in (0, 1, 2):
                        for t in (1, 2):
                            M2 = 2 * (3 ** d2 + 2 ** d2) * (2 ** d1 + 1) * 2 ** (m1 + i2 - 1) * 3 ** m2 \
                                - (1 if (i2 == 1 and m1 == 0) else 0)
                            N, q, beats = build((i1 + d1 + m1, [i1, i1 + d1]), [(d2, i2, m2)], t)
                            cnt += 1
                            ok &= (count_apex(N, q, adjacency(N, q, beats)) == 1 + 2 ** t * M2)
report("V6c two closers from a transitive pair: H = 1 + 2^t [2 (2^d1+1)(3^d2+2^d2) 2^(m1+i2-1) 3^m2 - [i2=1][m1=0]]", ok,
       "%d tournaments" % cnt)

# V7 recursion against two independent counts, random small towers; three-term recurrence
bad = 0
cnt = 0
for _ in range(80):
    n0 = rng.randint(1, 3)
    P0 = sorted(rng.sample(range(1, n0 + 1), rng.randint(0, 1)))
    depth = rng.randint(0, 2)
    levels = [(rng.randint(1, 2), rng.randint(1, 2), rng.randint(0, 1)) for _ in range(depth)]
    t = rng.randint(1, 2)
    N, q, beats = build((n0, P0), levels, t)
    A = adjacency(N, q, beats)
    if len(A) > 15:
        continue
    r, M, h, E = tower_state((n0, P0), levels)
    cnt += 1
    bad += not (H(A) == count_apex(N, q, A) == r + 2 ** t * M)
report("V7a towers: subset count = transfer-matrix count = r + 2^t M", bad == 0, "%d towers" % cnt)
rec = True
for _ in range(300):
    st = (rng.randint(1, 4), [1])
    levels = [(rng.randint(1, 2), rng.randint(1, 5), rng.randint(0, 3)) for _ in range(4)]
    S = [start_state(st)]
    for l in levels:
        S.append(step(S[-1], l))
    for k in range(1, 4):
        d, i, m = levels[k]
        ak = sum(levels[k - 1])
        rec &= (S[k + 1][1] == (2 ** d + 1) * 2 ** (m + ak) * S[k - 1][1]
                + 2 * (3 ** d + 2 ** d) * 2 ** (i - 1) * 3 ** m * S[k][1] - (S[k][3] if i == 1 else 0))
report("V7b three-term recurrence M_(k+1) = (2^d+1) 2^(m+a_k) M_(k-1) + 2 (3^d+2^d) 2^(i-1) 3^m M_k - [i=1] E_k", rec)

# V8 Mersenne line
A8 = adj(8, 14705281)
O8 = [1, 3, 4]
ok = (state(A8, O8) == (-1, 1024, 217, 64))
for t in (1, 2, 3, 4):
    ok &= (H(level_up(A8, O8, TT(t))) == 2 ** (10 + t) - 1)
report("V8a order-8 pair (code 14705281, out-set {1,3,4}): state (-1, 1024, 217, 64); H = 2^(10+t) - 1 on 9+t vertices, t <= 4", ok)
wit = {1: (1, 0), 2: (3, 2), 4: (5, 40), 5: (6, 6241), 6: (8, 135528708), 7: (7, 657682), 8: (8, 34259078),
       9: (8, 105437317), 10: (9, 3359819025)}
ok = all(H(adj(m, code)) == 2 ** s - 1 for s, (m, code) in wit.items())
report("V8b witnesses for 2^s - 1, s = 1, 2, 4, ..., 10 (orders <= 9)", ok)

# V9 neutral gadgets among all labelled tournaments of order <= 5
ok = True
neutral = 0
trans = 0
for m in range(1, 6):
    for A in labelled(m):
        n = m
        full = (1 << n) - 1
        en, st, hs = tables(A)
        pi2 = sum(hs[X] * hs[full ^ X] for X in range(full + 1))
        h = hs[full]
        Q2 = [[0] * n for _ in range(n)]
        Q3 = [[0] * n for _ in range(n)]
        for X in range(1, full + 1):
            rest = full ^ X
            Y = rest
            while Y:
                Z = rest ^ Y
                for e in range(n):
                    if X >> e & 1 and en[X][e]:
                        for f in range(n):
                            if Y >> f & 1 and st[Y][f]:
                                if Z == 0:
                                    Q2[e][f] += en[X][e] * st[Y][f]
                                else:
                                    Q3[e][f] += en[X][e] * st[Y][f] * hs[Z]
                Y = (Y - 1) & rest
        for Pm in range(1, 1 << n):
            P = [i for i in range(n) if Pm >> i & 1]
            w = sum(Q2[e][f] for e in P for f in P if e != f)
            w3 = sum(Q3[e][f] for e in P for f in P if e != f)
            if pi2 - 2 * w - 2 * w3 == 0:
                neutral += 1
                ok &= (len(P) == 2 and h - w == -h)
                if h == 1:
                    trans += 1
report("V9 gadgets with kappa = 0 (all labelled tournaments of order <= 5): always |P| = 2 and r_X = -H(X)", ok,
       "%d neutral labelled gadgets, %d with H(X) = 1" % (neutral, trans))
print("ALL PASS" if allok else "SOME FAILED")

# V10 explicit small pairs behind the rays through 13, 19, 27, 41 and 1039 (states recounted literally)
pairs = [((3, 5, 3), (1, 6), 1, 13), ((3, 0, 5), (-1, 5), 2, 19), ((3, 5, 1), (3, 3), 3, 27), ((5, 642, 3), (1, 20), 1, 41),
         ((7, 790536, 6), (-1, 130), 3, 1039), ((6, 4608, 38), (-17, 66), 4, 1039), ((6, 20612, 34), (15, 128), 3, 1039)]
ok = True
for (m, code, Om), (r0, M0), t, N in pairs:
    A = adj(m, code)
    O = [i for i in range(m) if Om >> i & 1]
    r, M, h, E = state(A, O)
    ok &= ((r, M) == (r0, M0) and r + 2 ** t * M == N and H(level_up(A, O, TT(t))) == N)
report("V10 explicit pairs: 13 = 1 + 2*6, 19 = -1 + 4*5, 27 = 3 + 8*3, 41 = 1 + 2*20, 1039 = -1 + 8*130 = -17 + 16*66 = 15 + 8*128", ok)
print("ALL PASS (with V10)" if allok else "SOME FAILED")
