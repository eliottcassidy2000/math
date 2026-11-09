#!/usr/bin/env python3
"""AUDIT: THM-5070 (2) apex bridges and (3) the tower, by brute force (packing polynomials from a1_blocks.py).

(2) T = TT_k + u (u beats chain vertex i iff bit i of W, most significant first), B = TT_l + w (chain vertex j beats
    w iff bit j of W', least significant first). Glue T (.) B at (u, w). Check, for x = 2..8:
    - I(Omega(T (.) B), x) equals the displayed formula;
    - E, O of (T, u) equal ceil(W/2), 1 + floor(W/2); q_e, q_o of (B, w) equal ceil(W'/2), 1 + floor(W'/2);
    - a(W) (number of odd cycles of T) satisfies the recurrence, and a(W') is the odd-cycle count of B.
(3) M_T = [[1,x,0],[0,0,1],[1,1,0]]: det(tI - M_T) = t^3 - t^2 - t - (x-1) on a grid; (t-2)(t^2+t+1) at x = 3;
    M_T (3,1,2) = 2 (3,1,2); iterating the one-vertex block gives the tower (each vertex beats only its predecessor),
    and brute-force H_3 of the tower 1..14 equals (3*2^n + c_n)/7 with c_n = 1, -5, 4 (period 3).
"""
import random
import sys

sys.path.insert(0, __import__("os").path.dirname(__file__))
from a1_blocks import glue, state_polys, ev, odd_cycle_masks  # noqa: E402


def apex_T(k, W):
    n = k + 1
    A = [[0] * n for _ in range(n)]
    for i in range(1, n):
        for j in range(i + 1, n):
            A[i][j] = 1
    for i in range(1, n):
        b = W >> (k - i) & 1
        if b:
            A[0][i] = 1
        else:
            A[i][0] = 1
    return A


def apex_B(l, Wp):
    n = l + 1  # vertex 0 = w, chain 1..l
    A = [[0] * n for _ in range(n)]
    for i in range(1, n):
        for j in range(i + 1, n):
            A[i][j] = 1
    for j in range(1, n):
        if Wp >> (j - 1) & 1:
            A[j][0] = 1
        else:
            A[0][j] = 1
    return A


def a_rec(W):
    if W == 0:
        return 0
    if W % 2:
        return a_rec(W // 2)
    return a_rec(W // 2) + (W // 2 + 1) // 2


def q_counts(B, w, x):
    """(q_e, q_o): paths of B ending at w, weighted by I(Omega(B - V(Q)), x), by reversing the arcs."""
    n = len(B)
    R = [[B[j][i] for j in range(n)] for i in range(n)]
    h, E, O = state_polys(n, R, w)   # paths from w in the converse = paths to w in B; converse keeps packings
    return ev(E, x), ev(O, x)


def main():
    rng = random.Random(5070)
    fails = checks = 0
    for trial in range(400):
        k, l = rng.randint(1, 6), rng.randint(1, 6)
        W, Wp = rng.randrange(1 << k), rng.randrange(1 << l)
        T, B = apex_T(k, W), apex_B(l, Wp)
        A, _ = glue(T, 0, B, 0, 0)
        hT, ET, OT = state_polys(len(T), T, 0)
        hB = state_polys(len(B), B, 0)[0]
        hA = state_polys(len(A), A, 0)[0]
        aT, aB = sum(odd_cycle_masks(len(T), T).values()), sum(odd_cycle_masks(len(B), B).values())
        for x in range(2, 9):
            checks += 1
            ce, fl = (W + 1) // 2, W // 2
            cep, flp = (Wp + 1) // 2, Wp // 2
            formula = (1 + x * a_rec(W)) * (1 + x * a_rec(Wp)) + x * (ce * (1 + flp) + (1 + fl) * cep)
            qe, qo = q_counts(B, 0, x)
            ok = (ev(hA, x) == formula and ev(ET, x) == ce and ev(OT, x) == 1 + fl and qe == cep and qo == 1 + flp
                  and aT == a_rec(W) and aB == a_rec(Wp) and ev(hT, x) == 1 + x * aT and ev(hB, x) == 1 + x * aB)
            if not ok:
                fails += 1
                if fails < 5:
                    print("APEX FAIL", k, W, l, Wp, x, ev(hA, x), formula)
    print(f"(2) apex bridges: {checks} checks (400 random pairs, k, l <= 6, x = 2..8), {fails} failures")

    # (3) characteristic polynomial on a grid
    bad = 0
    for x in range(-3, 6):
        for t in range(-4, 6):
            M = [[t - 1, -x, 0], [0, t, -1], [-1, -1, t]]
            det = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                   + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
            if det != t ** 3 - t ** 2 - t - (x - 1):
                bad += 1
            if x == 3 and det != (t - 2) * (t * t + t + 1):
                bad += 1
    MT = [[1, 3, 0], [0, 0, 1], [1, 1, 0]]
    v = [3, 1, 2]
    Mv = [sum(MT[r][c] * v[c] for c in range(3)) for r in range(3)]
    print(f"(3) det(tI - M_T) = t^3 - t^2 - t - (x - 1) on a 9 x 10 grid, (t-2)(t^2+t+1) at x = 3: {bad} failures;"
          f" M_T (3,1,2) = {Mv}")
    # tower: vertex i beats only i-1
    bad = 0
    state = (1, 0, 1)
    for nv in range(1, 15):
        A = [[0] * nv for _ in range(nv)]
        for i in range(nv):
            for j in range(nv):
                if i != j:
                    A[i][j] = 1 if (i == j + 1 or j >= i + 2) else 0
        h, E, O = state_polys(nv, A, nv - 1)
        c = [1, -5, 4][(nv - 1) % 3]
        pred = (3 * 2 ** nv + c) // 7
        if (ev(h, 3), ev(E, 3), ev(O, 3)) != state or ev(h, 3) != pred or (3 * 2 ** nv + c) % 7:
            bad += 1
            print("TOWER FAIL", nv, ev(h, 3), state, pred)
        state = (state[0] + 3 * state[1], state[2], state[0] + state[1])
    print(f"(3) tower 1..14: brute-force (h, E, O) at x = 3 equals M_T^(n-1) (1,0,1) and h = (3*2^n + c_n)/7: {bad} failures")


if __name__ == "__main__":
    main()
