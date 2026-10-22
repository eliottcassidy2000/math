# Fast version of the C_3-substitution path-cover recursion (see ternary_recursion.py).
# The generating function of single words (start block, length) in the three block variables is
#     W = (e1 + e2 + 3 e3) / (1 - e3),     e_i elementary symmetric in (xA, xB, xC),
# so the ordered j-cover count is c_j(C_3[A,B,C]) = E[W^j], E(xA^a xB^b xC^c) = cA_a cB_b cC_c (a, b, c >= 1).
# One DP step multiplies the occurrence array by the 7-term numerator and divides by (1 - e3) (a diagonal
# prefix sum).  Big integers via numpy object arrays.  mac-mini-2026-10-22-sierp.
import sys, math
import numpy as np

def step(P, M):
    """P -> P * (e1 + e2 + 3 e3) / (1 - e3), truncated to exponents <= M"""
    R = np.zeros_like(P)
    R[1:, :, :] += P[:-1, :, :]
    R[:, 1:, :] += P[:, :-1, :]
    R[:, :, 1:] += P[:, :, :-1]
    R[1:, 1:, :] += P[:-1, :-1, :]
    R[:, 1:, 1:] += P[:, :-1, :-1]
    R[1:, :, 1:] += P[:-1, :, :-1]
    R[1:, 1:, 1:] += 3 * P[:-1, :-1, :-1]
    # divide by (1 - e3): Q[m] = sum_(t >= 0) R[m - t(1,1,1)], by a doubling (Hillis-Steele) scan along the diagonal
    Q = R
    s = 1
    while s <= M:
        Q2 = Q.copy()
        Q2[s:, s:, s:] += Q[:M + 1 - s, :M + 1 - s, :M + 1 - s]
        Q = Q2
        s *= 2
    return Q

def covers(b, J):
    """b[m] = ordered m-covers of the (common) block, m = 0..M (b[0] unused).  Returns c[0..J] of C_3[B, B, B]."""
    M = len(b) - 1
    bv = np.array([0] + [int(x) for x in b[1:]], dtype=object)
    P = np.zeros((M + 1, M + 1, M + 1), dtype=object); P[0, 0, 0] = 1
    c = [0] * (J + 1)
    for j in range(1, J + 1):
        P = step(P, M)
        # contraction sum_m P[m] b[m1] b[m2] b[m3]
        t = np.tensordot(P, bv, axes=([2], [0]))
        t = np.tensordot(t, bv, axes=([1], [0]))
        c[j] = int(np.dot(t, bv))
    return c

def H_only(b):
    """H(C_3[B,B,B]) = 3 sum_(q>=0) (b_(q+1) b_q^2 + b_(q+1)^2 b_q + b_(q+1)^3), b_0 = 0"""
    M = len(b) - 1
    bb = list(b) + [0]
    bb[0] = 0
    return 3 * sum(bb[q + 1] * bb[q] ** 2 + bb[q + 1] ** 2 * bb[q] + bb[q + 1] ** 3 for q in range(M))

def v3(x):
    v = 0
    while x % 3 == 0: x //= 3; v += 1
    return v

if __name__ == '__main__':
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    b = [1, 1]                       # a single vertex
    for k in range(1, K + 1):
        n = 3 ** k
        Hk = H_only(b)
        aut = (3 ** k - 1) // 2
        print(f"k={k} n={n}: H(S_k) = {Hk}")
        print(f"   v_3(H) = {v3(Hk)}  v_3|Aut| = {aut}  (H/|Aut|) mod 3 = {(Hk // 3 ** aut) % 3}  H*2^(n-1)/n! = {Hk * 2 ** (n - 1) / math.factorial(n):.6f}", flush=True)
        if k < K:
            b = covers(b, n)
            assert b[1] == Hk, (k, b[1], Hk)
            assert b[n] == math.factorial(n)
            print(f"   ordered covers c_j(S_k) computed for j <= {n}; c_n = n! checked", flush=True)
