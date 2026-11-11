"""aud_f_structure.py - independent checks of THM-5130 (F), of the remark "15 lies on no ray of offset +1 or -1",
of H3 (i) at order 6 and of a few explicit states quoted in the write-up.

F   all labelled tournaments of order <= 5, Q-matrix by plain permutation enumeration (aud_core.brute_Q logic):
    F1  U_B(o) = H(B) + sum over out-neighbours x of o of Q_B(x, o)
    F2a Q_B(f, e) = prod H * U_{C_a}(e) * prod Pi_2 * V_{C_b}(f) * prod H  >= H(B)   (e in an earlier component)
        and Q_B(e, f) = [b = a+1] prod H * End_{C_a}(e) * Start_{C_b}(f) * prod H   (auditor's companion formula)
    F2b r(B, O) = r(C_a, O) prod_{i != a} H(C_i);  U_B(O) = U_{C_a}(O) prod_{i<a} H(C_i) prod_{i>a} Pi_2(C_i)
    F2c unit sets with |O| >= 2: one component, others single vertices, M = 2^(after) M_core
    F3  twins: W = 2 H(K);  F4: W_ef independent of the arc ef
S   all labelled pairs of order <= 6 (subset recursion): slopes of seeds / anti-seeds; 7 is not a seed slope;
    1, 2, 4, 8 are not anti-seed slopes; small-slope lists; M >= delta_dom with equality iff B transitive (order 6).
X   explicit states quoted in the write-up.
"""
import itertools, sys, time
from aud_core import *

allok = True


def report(name, ok, extra=""):
    global allok
    allok &= bool(ok)
    print(("PASS  " if ok else "FAIL  ") + name + ("  " + extra if extra else ""))
    sys.stdout.flush()


def brute_all(A):
    """H, Q[e][f], U[o], V[o], Start[o], End[o], Pi2 by permutation enumeration."""
    n = len(A)
    H = 0
    Q = [[0] * n for _ in range(n)]
    U = [0] * n
    V = [0] * n
    St = [0] * n
    En = [0] * n
    Pi = 0
    for p in itertools.permutations(range(n)):
        br = [k for k in range(1, n) if not A[p[k - 1]][p[k]]]
        if len(br) >= 2:
            continue
        if not br:
            H += 1
            St[p[0]] += 1
            En[p[-1]] += 1
            Pi += n + 1
            for k in range(0, n):
                U[p[k]] += 1               # P1 = p[:k] possibly empty, P2 = p[k:] starts at p[k]
            for k in range(1, n + 1):
                V[p[k - 1]] += 1           # P1 = p[:k] ends at p[k-1], P2 possibly empty
            for k in range(1, n):
                Q[p[k - 1]][p[k]] += 1
        else:
            k = br[0]
            Pi += 1
            U[p[k]] += 1
            V[p[k - 1]] += 1
            Q[p[k - 1]][p[k]] += 1
    return H, Q, U, V, St, En, Pi


def components(A):
    """strong components in their linear order (earlier => later)."""
    n = len(A)
    reach = [[i == j or bool(A[i][j]) for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if reach[i][k] and reach[k][j]:
                    reach[i][j] = True
    comps = []
    seen = set()
    for i in range(n):
        if i in seen:
            continue
        c = [j for j in range(n) if reach[i][j] and reach[j][i]]
        seen |= set(c)
        comps.append(c)
    comps.sort(key=lambda c: -sum(1 for j in range(n) if reach[c[0]][j]))   # earlier components reach more
    return comps


t0 = time.time()
f1 = f2a = f2a2 = f2b = f2c = f4 = True
units = 0
cache = {}
for n in range(1, 6):
    for code in range(1 << (n * (n - 1) // 2)):
        A = adj_from_code(n, code)
        H, Q, U, V, St, En, Pi = brute_all(A)
        cache[(n, code)] = (H, Q, U)
        for o in range(n):
            f1 &= (U[o] == H + sum(Q[x][o] for x in range(n) if A[o][x]))
        comps = components(A)
        data = []
        for c in comps:
            sub = [[A[x][y] for y in c] for x in c]
            data.append(brute_all(sub))
        cid = {v: k for k, c in enumerate(comps) for v in c}
        loc = {v: c.index(v) for c in comps for v in c}
        # check order: all arcs go from earlier to later component
        assert all(A[x][y] for a in range(len(comps)) for b in range(a + 1, len(comps)) for x in comps[a] for y in comps[b])
        assert H == eval("*".join(str(d[0]) for d in data))
        for e in range(n):
            for f in range(n):
                a, b = cid[e], cid[f]
                if a < b:
                    val = 1
                    for i, d in enumerate(data):
                        if i < a or i > b:
                            val *= d[0]
                        elif a < i < b:
                            val *= d[6]
                    val *= data[a][2][loc[e]] * data[b][3][loc[f]]
                    f2a &= (Q[f][e] == val and val >= H)
                    val2 = 0
                    if b == a + 1:
                        val2 = data[a][5][loc[e]] * data[b][4][loc[f]]
                        for i, d in enumerate(data):
                            if i < a or i > b:
                                val2 *= d[0]
                    f2a2 &= (Q[e][f] == val2)
        for Om in range(1, 1 << n):
            O = [i for i in range(n) if Om >> i & 1]
            W = sum(Q[e][f] for e in O for f in O if e != f)
            M = sum(U[o] for o in O)
            cs = {cid[o] for o in O}
            if len(cs) >= 2:
                f2a &= (H - W <= 0)
            if len(cs) == 1:
                a = cs.pop()
                d = data[a]
                Wc = sum(d[1][loc[e]][loc[f]] for e in O for f in O if e != f)
                Mc = sum(d[2][loc[o]] for o in O)
                rc = d[0] - Wc
                prodH = 1
                before = 1
                after = 1
                for i, dd in enumerate(data):
                    if i != a:
                        prodH *= dd[0]
                    if i < a:
                        before *= dd[0]
                    if i > a:
                        after *= dd[6]
                f2b &= (H - W == rc * prodH and M == Mc * before * after)
                if H - W == 1 and len(O) >= 2:
                    units += 1
                    nafter = sum(len(c) for i, c in enumerate(comps) if i > a)
                    f2c &= (all(len(c) == 1 for i, c in enumerate(comps) if i != a) and rc == 1 and M == 2 ** nafter * Mc)
            elif H - W == 1 and len(O) >= 2:
                f2c = False
        # F4
        for e in range(n):
            for f in range(e + 1, n):
                B2 = [row[:] for row in A]
                B2[e][f], B2[f][e] = A[f][e], A[e][f]
                b = 0
                code2 = 0
                for i in range(n):
                    for j in range(i + 1, n):
                        if B2[i][j]:
                            code2 |= 1 << b
                        b += 1
                if (n, code2) in cache:
                    H2, Q2, U2 = cache[(n, code2)]
                    f4 &= (Q2[e][f] == Q[e][f] and Q2[f][e] == Q[f][e])
report("F1 U_B(o) = H(B) + sum_{o -> x} Q_B(x, o)", f1)
report("F2a Q_B(f, e) product formula and >= H(B); r <= 0 when O meets two components", f2a)
report("F2a' (companion) Q_B(e, f) = [adjacent components] prod H * End * Start", f2a2)
report("F2b r and U factor over the components when O lies in one component", f2b)
report("F2c unit sets with |O| >= 2", f2c, "%d unit sets" % units)
report("F4 W_ef does not depend on the orientation of the arc ef", f4, "%.0fs" % (time.time() - t0))
# F3 twins
import random
rng = random.Random(3)
ok = True
for _ in range(150):
    n = rng.randint(1, 5)
    K = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                K[i][j] = 1
            else:
                K[j][i] = 1
    y = rng.randrange(n)
    D = [[0] * (n + 1) for _ in range(n + 1)]
    for i in range(n):
        for j in range(n):
            D[i][j] = K[i][j]
        D[i][n] = K[i][y]
        D[n][i] = K[y][i]
    D[y][n] = 1
    D[n][y] = 0
    H2, Q2, U2, V2, S2, E2, P2 = brute_all(D)
    HK = brute_H(K)
    ok &= (Q2[y][n] == HK and Q2[n][y] == HK and brute_state(D, [y, n])[0] == H2 - 2 * HK)
report("F3 twin lemma: Q(e,f) = Q(f,e) = H(K), r = H(K_2) - 2 H(K)", ok)

# ---------------------------------------------------------------- S: slopes of unit pairs, order <= 6
t0 = time.time()
seed_slopes = set()
anti_slopes = set()
seed_min = {}
bad_dom = eq_nontrans = trans_noneq = 0
pairs6 = 0
seedO_sizes = set()
for n in range(1, 7):
    for code in range(1 << (n * (n - 1) // 2)):
        A = adj_from_code(n, code)
        endat, startat, tot = dp_tables(A)
        full = (1 << n) - 1
        h = tot[full]
        U = [0] * n
        Q = [[0] * n for _ in range(n)]
        for S in range(0, full):
            R = full ^ S
            tS = tot[S]
            for o in range(n):
                if (R >> o) & 1:
                    s = startat[R][o]
                    if s:
                        U[o] += tS * s
                        if S:
                            eS = endat[S]
                            for e in range(n):
                                if (S >> e) & 1 and eS[e]:
                                    Q[e][o] += eS[e] * s
        outm = [sum(1 << j for j in range(n) if A[i][j]) for i in range(n)]
        for Om in range(1, 1 << n):
            O = [i for i in range(n) if Om >> i & 1]
            W = sum(Q[e][f] for e in O for f in O if e != f)
            M = sum(U[o] for o in O)
            r = h - W
            if r == 1:
                seed_slopes.add(M)
                seedO_sizes.add(len(O))
            elif r == -1:
                anti_slopes.add(M)
            if n == 6:
                delta = 0
                sub = Om
                while sub:
                    un = 0
                    for o in O:
                        if sub >> o & 1:
                            un |= outm[o]
                    delta += 2 ** bin(un & ~Om).count("1")
                    sub = (sub - 1) & Om
                pairs6 += 1
                bad_dom += (M < delta)
                if M == delta and h > 1:
                    eq_nontrans += 1
                if h == 1 and M != delta:
                    trans_noneq += 1
print("   seed slopes <= 40 (order <= 6):", sorted(x for x in seed_slopes if x <= 40))
print("   anti-seed slopes <= 40 (order <= 6):", sorted(x for x in anti_slopes if x <= 40))
print("   sizes of unit sets (order <= 6):", sorted(seedO_sizes))
report("S1 no seed of slope 7 (nor 3, 5, 10), no anti-seed of slope 1, 2, 4, 8 (nor 11) among all labelled pairs of order <= 6",
       not ({3, 5, 7, 10} & seed_slopes) and not ({1, 2, 4, 8, 11} & anti_slopes))
small_seed = sorted(x for x in seed_slopes if x < 10)
report("S2 seed slopes below 10 are 1, 2, 4, 6, 8 (powers of two and 6); anti-seed slopes below 9 are 3, 5, 6",
       small_seed == [1, 2, 4, 6, 8] and sorted(x for x in anti_slopes if x < 9) == [3, 5, 6])
report("S3 H3 (i) at order 6: M >= delta_dom for every labelled pair, equality iff B transitive",
       bad_dom == 0 and eq_nontrans == 0 and trans_noneq == 0, "%d pairs, %.0fs" % (pairs6, time.time() - t0))

# ---------------------------------------------------------------- X: explicit states
checks = [((6, 20612, 34), (15, 128, 31)), ((3, 5, 3), (1, 6, 3)), ((3, 0, 5), (-1, 5, 1)), ((5, 642, 3), (1, 20, 9)),
          ((7, 790536, 6), (-1, 130, 39)), ((6, 4608, 38), (-17, 66, 3)), ((6, 17410, 18), (1, 64, 15)),
          ((4, 34, 12), (1, 16, 5)), ((5, 642, 10), (1, 32, 9))]
ok = True
for (m, code, Om), exp in checks:
    A = adj_from_code(m, code)
    O = [i for i in range(m) if Om >> i & 1]
    st = brute_state(A, O)
    good = (st[:3] == exp)
    ok &= good
    print("   (order %d, code %d, mask %d): state (r, M, h, E) = %s %s" % (m, code, Om, st, "" if good else "MISMATCH, expected %s" % (exp,)))
report("X explicit states quoted in the write-up (first three coordinates)", ok)
print("   note: (6, 20612, mask 34) has H(B) = 31 and slope 128 with offset 15: 'no pair of order <= 9 with H(B) = 31 and slope 64 or 128' is false without the condition r = 1")
print("ALL PASS" if allok else "SOME FAILED")
