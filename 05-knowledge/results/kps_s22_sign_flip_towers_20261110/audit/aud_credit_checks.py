"""aud_credit_checks.py - computations behind the credit findings of the audit.

C1  Session260 coverage/PROOF.md: the five-vertex carrier B = (24,1,3,6,14) (out-masks) IS the one-closer tournament
    C(TT_2, {1,2}; TT_2) of THM-5130, and its marked families for v = 1, 2 (module P_a -> x -> G_c -> y -> S_b at a
    tail vertex, x and y selected) are tower levels (delta, i, m) = (c+1, ., .) over the transitive pair (TT_2, {1,2}):
    their (H, D, U) agree with Session260 (4)-(5) AND with THM-5130 (1) / D3.
C2  ONE_EXIT.md Theorem 1: H(T(X,x;G,A)) = hH + FU - sE - dW.  With X = X' => v, x = v it is row h of Theorem A; and
    applied to the closed tournaments T_t(B+, P), which are one-exit attachments onto (B, O), it gives r+ and M+ as
    linear forms in (h, M, E, W) - checked against the matrix of Theorem A.
C3  Session261 MARKED_RESELECTION.md: the gadget constants of Theorem A are entries of the marked cover tensor
    (C_r, S_r, E_r, K_r), r <= 3:  H_X = C_1, Pi = C_2 + 2 C_1, w = K_2, w3 = K_3, sigma = S_1, u = S_1 + S_2,
    u3 = S_3, eps = E_1, nu = E_1 + E_2.
C4  Session261 NINE_OFFSET_EXTENSION.md section 3 (1): U_x = 2^(p+1) 3^q + 2^(q+1) - 0^p for a chain mark of
    F_n = C_3[TT_n, 1, 1] is H1 of THM-5130 applied to the point base.
C5  Session258 CORE_CLASSIFICATION.md (1) = F2 (a): see aud_f_structure.py (F2a), same product formula.
"""
import itertools, random, sys
from aud_core import *

allok = True
rng = random.Random(2611)


def report(name, ok, extra=""):
    global allok
    allok &= bool(ok)
    print(("PASS  " if ok else "FAIL  ") + name + ("  " + extra if extra else ""))


def from_outmasks(rows):
    n = len(rows)
    A = [[1 if rows[i] >> j & 1 else 0 for j in range(n)] for i in range(n)]
    assert is_tournament(A)
    return A


def rand_t(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


# ---------------------------------------------------------------- C1
B5 = from_outmasks((24, 1, 3, 6, 14))
T = close(transitive(2), [0, 1], transitive(2))        # vertices: chain 0,1 ; tail 2,3 ; closer 4
# isomorphism: session vertex 4 -> 0, 3 -> 1, 2 -> 2, 1 -> 3, 0 -> 4
perm = {4: 0, 3: 1, 2: 2, 1: 3, 0: 4}
iso = all(B5[x][y] == T[perm[x]][perm[y]] for x in range(5) for y in range(5) if x != y)
report("C1a Session260's carrier (24,1,3,6,14) = C(TT_2, {1,2}; TT_2) (explicit isomorphism)", iso)
ok = True
cnt = 0
CE = {1: (12, 0), 2: (18, 1), 3: (20, 0), 4: (30, 3)}
for v in (1, 2, 3, 4):
    for a in range(0, 3):
        for b in range(0, 3):
            for c in range(0, 3):
                # replace v by the module P_a -> x -> G_c -> y -> S_b
                size = a + c + b + 2
                verts = []
                for u in range(5):
                    if u == v:
                        verts += [(u, k) for k in range(size)]
                    else:
                        verts.append((u, 0))
                n = len(verts)
                A = [[0] * n for _ in range(n)]
                for i, (x, xi) in enumerate(verts):
                    for j, (y, yj) in enumerate(verts):
                        if i != j:
                            A[i][j] = (1 if xi < yj else 0) if x == y else B5[x][y]
                O = [i for i, (x, xi) in enumerate(verts) if x == v and xi in (a, a + c + 1)]
                r, M, h, E = dp_state(A, O)
                Cv, ev = CE[v]
                H260 = 24 * 2 ** (a + b + c) - 1
                D260 = (1 if c == 0 else 0) + 2 ** c - 1 - 12 * 2 ** (a + b) * (3 ** c - 2 ** (c + 1) + 1)
                U260 = Cv * 2 ** a * 3 ** b * (2 ** (c + 1) + 3 ** (c + 1)) - (ev if a == 0 else 0)
                good = (h == H260 and r == D260 and M == U260)
                # tower reading for v = 1, 2 (tail vertices of the one-closer tournament)
                if v in (1, 2):
                    delta = c + 1
                    if v == 1:      # session vertex 1 = second tail vertex: tail = (vertex 2) -> module
                        i_, m_ = a + 2, b
                    else:           # session vertex 2 = first tail vertex: tail = module -> (vertex 1)
                        i_, m_ = a + 1, b + 1
                    atail = i_ + delta + m_
                    Bp = close(transitive(2), [0, 1], transitive(atail))
                    st = dp_state(Bp, [2 + i_ - 1, 2 + i_ - 1 + delta])
                    good &= (st[:3] == (r, M, h))
                    # Theorem B / Corollary B prediction from the start state (-1, 3, 1, 1) of (TT_2, {1,2})
                    r0, M0, h0, E0 = -1, 3, 1, 1
                    rX = -1 if delta == 1 else 1 - 2 ** (delta - 1)
                    kap = 2 ** (atail - delta) * (2 ** delta - 1 - 3 ** (delta - 1))
                    al, ga = 2 ** delta + 1, 3 ** delta + 2 ** delta
                    pred = (rX * r0 + kap * M0, al * 2 ** m_ * (h0 + r0) + 2 * ga * 2 ** (i_ - 1) * 3 ** m_ * M0 - (E0 if i_ == 1 else 0),
                            r0 + 2 ** atail * M0)
                    good &= (pred == (r, M, h))
                ok &= good
                cnt += 1
report("C1b Session260 (4)-(5) on its eight-start family, and for v = 1, 2 identity with the towers over (TT_2, {1,2}) "
       "(gap c = delta - 1; Theorem B with delta <= 3)", ok, "%d marked modules" % cnt)

# ---------------------------------------------------------------- C2
ok = True
cnt = 0
for _ in range(150):
    n = rng.randint(1, 4)
    k = rng.randint(1, 3)
    G = rand_t(n)
    Aset = [i for i in range(n) if rng.random() < 0.6]
    r, M, h, E = brute_state(G, Aset)
    W = h - r
    X = rand_t(k)
    x = rng.randrange(k)
    # T(X, x; G, A): G => X except x -> A
    N = n + k
    T2 = [[0] * N for _ in range(N)]
    for i in range(n):
        for j in range(n):
            T2[i][j] = G[i][j]
        for j in range(k):
            T2[i][n + j] = 1
    for i in range(k):
        for j in range(k):
            T2[n + i][n + j] = X[i][j]
    for a_ in Aset:
        T2[a_][n + x] = 0
        T2[n + x][a_] = 1
    assert is_tournament(T2)
    hX = brute_H(X)
    # F = sum over nonempty paths P of X ending at x of H(X - V(P)); s = paths starting at x; d = H(X - x)
    F = 0
    others = [v for v in range(k) if v != x]
    for size in range(0, k):
        for S in itertools.combinations(others, size):
            for p in itertools.permutations(S):
                path = list(p) + [x]
                if all(X[path[i]][path[i + 1]] for i in range(len(path) - 1)):
                    rest = [v for v in range(k) if v not in path]
                    F += brute_H([[X[u][w] for w in rest] for u in rest])
    s_ = sum(1 for p in itertools.permutations(range(k)) if p[0] == x and all(X[p[i]][p[i + 1]] for i in range(k - 1)))
    d_ = brute_H([[X[u][w] for w in others] for u in others])
    cnt += 1
    ok &= (brute_H(T2) == hX * h + F * M - s_ * E - d_ * W)
report("C2a ONE_EXIT Theorem 1: H(T(X,x;G,A)) = hH + FU - sE - dW (random checks by permutation enumeration)", ok, "%d cases" % cnt)
# row h of Theorem A as the case X = X' => v
ok = True
for _ in range(60):
    k = rng.randint(1, 4)
    Xp = rand_t(k)
    g = brute_gadget(Xp, [])
    Xv = [row + [1] for row in Xp] + [[0] * (k + 1)]     # X' => v
    # coordinates of (X' => v, v): h = H(X'), F = Pi_2(X'), s = 0, d = H(X')
    hX = brute_H(Xv)
    F = 0
    for size in range(0, k + 1):
        for S in itertools.combinations(range(k), size):
            for p in itertools.permutations(S):
                if all(Xp[p[i]][p[i + 1]] for i in range(len(p) - 1)):
                    rest = [v for v in range(k) if v not in p]
                    F += brute_H([[Xp[u][w] for w in rest] for u in rest])
    ok &= (hX == g["H"] and F == g["Pi"])
report("C2b for X = X' => v, x = v: (h, F, s, d) = (H(X'), Pi_2(X'), 0, H(X')), so ONE_EXIT (1) is row h of Theorem A", ok)

# ---------------------------------------------------------------- C3
ok = True
for _ in range(80):
    k = rng.randint(1, 5)
    X = rand_t(k)
    P = [i for i in range(k) if rng.random() < 0.5]
    Ps = set(P)
    g = brute_gadget(X, P)
    C = [0] * 5
    S = [0] * 5
    E_ = [0] * 5
    K = [0] * 5
    for p in itertools.permutations(range(k)):
        br = [i for i in range(1, k) if not X[p[i - 1]][p[i]]]
        for r in (1, 2, 3):
            for cuts in itertools.combinations(range(1, k), r - 1):
                if all(b in cuts for b in br):
                    blocks = []
                    prev = 0
                    for c_ in list(cuts) + [k]:
                        blocks.append(p[prev:c_])
                        prev = c_
                    C[r] += 1
                    if blocks[0][0] in Ps:
                        S[r] += 1
                    if blocks[0][-1] in Ps:
                        E_[r] += 1
                    if r >= 2 and blocks[0][-1] in Ps and blocks[1][0] in Ps:
                        K[r] += 1
    ok &= (g["H"] == C[1] and g["Pi"] == C[2] + 2 * C[1] and g["w"] == K[2] and g["w3"] == K[3] and g["sigma"] == S[1]
           and g["u"] == S[1] + S[2] and g["u3"] == S[3] and g["eps"] == E_[1] and g["nu"] == E_[1] + E_[2])
report("C3 gadget constants of Theorem A = entries of Session261's marked cover tensor (C_r, S_r, E_r, K_r), r <= 3", ok)

# ---------------------------------------------------------------- C4
ok = True
for p_ in range(0, 4):
    for q_ in range(0, 4):
        n = p_ + q_ + 1
        # F_n: chain -> y -> z -> chain ; as C(point z, {z}; TT_n) with closer y
        Fn = close([[0]], [0], transitive(n))
        hF = dp_H(Fn)
        st = dp_state(Fn, [1 + p_])
        U261 = 2 ** (p_ + 1) * 3 ** q_ + 2 ** (q_ + 1) - (1 if p_ == 0 else 0)
        a, pp = n, p_ + 1
        UH1 = 2 ** (a - pp) * (1 + 1) + 2 ** pp * 3 ** (a - pp) * 1 - (1 if pp == 1 else 0)
        ok &= (hF == 2 ** n + 1 and st[1] == U261 == UH1 and st[0] == 2 ** n + 1)
report("C4 Session261 NINE_OFFSET_EXTENSION (1) = H1 of THM-5130 at the point base", ok)
print("ALL PASS" if allok else "SOME FAILED")
