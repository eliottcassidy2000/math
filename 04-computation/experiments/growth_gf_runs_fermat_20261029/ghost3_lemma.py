"""Mod-3 ghost lemma for the three-run observer (HYP-9445 (4)).

3-locally the S_3 Fourier image is {(u, v, W) : W' = Q^-1 W Q = [[w11, w12], [w21, w22]], w21 = 0, w11 = u, w22 = v
(mod 3)} with Q = [[1, 1], [1, 0]] (columns: the invariant line e = (1,1) and f = (1,0)).  M acts by u -> P u,
v -> S v and W' -> L' W' (left multiplication, columns separately), with L' = [[L11, L12], [L21, L22]],
L21 = 0, L11 = P, L22 = S (mod 3).  In the Z_(3)-basis (u, v, y11, w12, y21, y22) of the image
(w11 = u + 3 y11, w21 = 3 y21, w22 = v + 3 y22), M reduces mod 3 to two block-lower-triangular pieces
    I  on (u, y21, y11):  [[P, 0, 0], [D2, S, 0], [D1, L12, P]],
    II on (v, w12, y22):  [[S, 0, 0], [L12, P, 0], [D3, D2, S]],
D1 = (L11 - P)/3, D2 = L21/3, D3 = (L22 - S)/3.  Hence
    rank_F3 M = rank I + rank II,   spur_3 = rank I + rank II - rank P - rank S - 2 rank [[P, L12], [0, S]],
and when S is invertible mod 3 the right side is rank(E : ker P -> coker P) with E = D1 - L12 S^-1 D2.
Checked against direct mod-3 ranks of M and N for every class, n = 4, 5, 6.
"""
from runs import tournaments, tour_from_bits, run_matrix, unglued_k3


def rank_mod(M, p):
    M = [[x % p for x in row] for row in M]
    r = 0
    cols = len(M[0]) if M else 0
    for c in range(cols):
        piv = next((i for i in range(r, len(M)) if M[i][c]), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [x * inv % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        r += 1
    return r


def blocks(N, r):
    """Extract P, S and the standard-sector operator L on the (a, c) column of W, in the adapted basis."""
    P = [[N[6 * i][6 * j] for j in range(r)] for i in range(r)]
    S = [[N[6 * i + 1][6 * j + 1] for j in range(r)] for i in range(r)]
    # W = [[a, b], [c, d]] at coordinates 2, 3, 4, 5; left multiplication acts on the column (a, c) = (2, 4)
    La = lambda i, j, s, t: N[6 * i + s][6 * j + t]
    L = {}
    for i in range(r):
        for j in range(r):
            L[(i, j)] = [[La(i, j, 2, 2), La(i, j, 2, 4)], [La(i, j, 4, 2), La(i, j, 4, 4)]]
    # adapted basis: L' = Q^-1 L Q, Q = [[1,1],[1,0]], Q^-1 = [[0,1],[1,-1]]
    def conj(m):
        q, qi = [[1, 1], [1, 0]], [[0, 1], [1, -1]]
        t = [[sum(qi[x][y] * m[y][z] for y in range(2)) for z in range(2)] for x in range(2)]
        return [[sum(t[x][y] * q[y][z] for y in range(2)) for z in range(2)] for x in range(2)]
    Lp = {key: conj(m) for key, m in L.items()}
    blk = lambda s, t: [[Lp[(i, j)][s][t] for j in range(r)] for i in range(r)]
    return P, S, blk(0, 0), blk(0, 1), blk(1, 0), blk(1, 1)


def assemble(rows):
    out = []
    for row in rows:
        for k in range(len(row[0])):
            out.append(sum((blk[k] for blk in row), []))
    return out


zero = lambda r: [[0] * r for _ in range(r)]
for n in (4, 5, 6):
    bad = nonzero = 0
    for bits in tournaments(n):
        A = tour_from_bits(n, bits)
        tuples, M = run_matrix(A, 3)
        r = len(tuples) // 6
        N, _, _, _, _ = unglued_k3(M, r)
        P, S, L11, L12, L21, L22 = blocks(N, r)
        for X, Y in ((L11, P), (L22, S)):
            assert all((X[i][j] - Y[i][j]) % 3 == 0 for i in range(r) for j in range(r))
        assert all(x % 3 == 0 for row in L21 for x in row)
        D1 = [[(L11[i][j] - P[i][j]) // 3 for j in range(r)] for i in range(r)]
        D2 = [[L21[i][j] // 3 for j in range(r)] for i in range(r)]
        D3 = [[(L22[i][j] - S[i][j]) // 3 for j in range(r)] for i in range(r)]
        Z = zero(r)
        I = assemble([[P, Z, Z], [D2, S, Z], [D1, L12, P]])
        II = assemble([[S, Z, Z], [L12, P, Z], [D3, D2, S]])
        Lp = assemble([[P, L12], [Z, S]])
        rP, rS = rank_mod(P, 3), rank_mod(S, 3)
        formula_M = rank_mod(I, 3) + rank_mod(II, 3)
        rM, rN = rank_mod(M, 3), rank_mod(N, 3)
        assert rN == rP + rS + 2 * rank_mod(Lp, 3)
        spur3 = rM - rN
        if formula_M != rM:
            bad += 1
        if spur3:
            nonzero += 1
    print('n=%d: classes where rank_F3 M != rank I + rank II: %d; classes with spur_3 > 0: %d' % (n, bad, nonzero))
