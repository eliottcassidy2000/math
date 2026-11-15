#!/usr/bin/env python3
"""Wave 38: the triangle trilinear form and its Z[sqrt -2] specialisation.

Cyclic join cyc(A, B, C): blocks A => B => C => A (every arc from A to B, B to C, C to A).
Port vector of a block X: P_k(X) = number of ordered k-tuples of disjoint nonempty paths covering X (P_1 = H).

CLAIM (proved in the README, wave 38):
    H(cyc(A, B, C)) = sum over admissible visit-count triples (a, b, c) of mult(a, b, c) P_a(A) P_b(B) P_c(C),
where a triple is admissible if, for some k >= 1, it is (k, k, k) [mult 3], a rotation of (k+1, k, k) [mult 1] or a
rotation of (k+1, k+1, k) [mult 1].  So H(cyc) is a symmetric trilinear form T(v_A, v_B, v_C).
With C = one vertex it is the bilinear form v_A^T [[3,1],[1,1]] v_B = Re(beta_A conj(beta_B)), beta = (H+P2) + H sqrt(-2)
(THM-5222).  The script
  (1) checks the claim against brute-force subset DP on random blocks,
  (2) prints the Gram matrices G_C = (T(e_a, e_b, v_C))_{a,b} for small C and their determinants,
  (3) checks the full symmetry H(cyc(A,B,C)) = H(cyc(B,A,C)),
  (4) prints the cubic T(v, v, v) for blocks with P_3 = 0.
Usage: python3 w38_triangle_form.py HQ_BINARY
"""
import random
import sys
from fractions import Fraction
from itertools import product

src = open(__file__.replace("w38_triangle_form.py", "w38_two_apex_additivity.py")).read().rsplit("main()", 1)[0]
sys.argv = sys.argv[:2]
ns = {"__name__": "helpers"}
exec(compile(src, "helpers", "exec"), ns)
hq_batch, gentourng, ports = ns["hq_batch"], ns["gentourng"], ns["ports"]
rng = random.Random(5242)


def mult(a, b, c):
    t = (a, b, c)
    lo, hi = min(t), max(t)
    if lo < 1:
        return 0
    if lo == hi:
        return 3
    if hi == lo + 1:
        return 1
    return 0


def T(vA, vB, vC):
    s = 0
    for a, x in enumerate(vA, 1):
        for b, y in enumerate(vB, 1):
            for c, z in enumerate(vC, 1):
                m = mult(a, b, c)
                if m:
                    s += m * x * y * z
    return s


def cyc(A, B, C):
    blocks = [A, B, C]
    offs = [0, len(A), len(A) + len(B)]
    n = sum(map(len, blocks))
    out = [0] * n
    for i, X in enumerate(blocks):
        for u in range(len(X)):
            mask = 0
            for w in range(len(X)):
                if (X[u] >> w) & 1:
                    mask |= 1 << (offs[i] + w)
            nxt = (i + 1) % 3
            for w in range(len(blocks[nxt])):
                mask |= 1 << (offs[nxt] + w)
            out[offs[i] + u] = mask
    return out


def rand_tour(k):
    out = [0] * k
    for i in range(k):
        for j in range(i + 1, k):
            if rng.random() < 0.5:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
    return out


# (1) the trilinear law against brute force
tests = []
for _ in range(400):
    A, B, C = (rand_tour(rng.randint(1, 4)) for _ in range(3))
    tests.append((A, B, C))
res = hq_batch([cyc(A, B, C) for A, B, C in tests])
bad = 0
for (A, B, C), (h, q) in zip(tests, res):
    if h != T(ports(A), ports(B), ports(C)):
        bad += 1
print(f"(1) trilinear law H(cyc(A,B,C)) = T(v_A, v_B, v_C): {len(tests)} random triples (blocks of order 1-4), "
      f"{bad} failures")

# (3) symmetry under swapping two blocks (converse cyclic order)
res2 = hq_batch([cyc(B, A, C) for A, B, C in tests])
asym = sum(1 for (h1, _), (h2, _) in zip(res, res2) if h1 != h2)
print(f"(3) H(cyc(A,B,C)) = H(cyc(B,A,C)) on all {len(tests)} triples: {asym} failures")


# (2) Gram matrices
def det(M):
    M = [[Fraction(x) for x in row] for row in M]
    n = len(M)
    d = Fraction(1)
    for i in range(n):
        p = next((r for r in range(i, n) if M[r][i] != 0), None)
        if p is None:
            return 0
        if p != i:
            M[i], M[p] = M[p], M[i]
            d = -d
        d *= M[i][i]
        for r in range(i + 1, n):
            f = M[r][i] / M[i][i]
            for c in range(i, n):
                M[r][c] -= f * M[i][c]
    return d


print("(2) Gram matrices G_C[a][b] = T(e_a, e_b, v_C), a, b = 1..|C|+1:")
smalls = [("point", [0])]
for k in (2, 3, 4):
    for i, X in enumerate(gentourng(k)):
        smalls.append((f"order{k}#{i}", X))
for name, X in smalls:
    v = ports(X)
    K = len(X) + 1
    G = [[T([1 if t == a else 0 for t in range(1, K + 1)], [1 if t == b else 0 for t in range(1, K + 1)], v)
          for b in range(1, K + 1)] for a in range(1, K + 1)]
    print(f"  C = {name}, ports {v}: det G_C = {det(G)}; G_C = {G}")

# (4) the cubic on blocks with P_3 = 0 (orders 1, 2)
print("(4) cubic T(v,v,v) for v = (P1, P2): 3 (P1 + P2)(P1^2 + P2^2)?",
      all(T([x, y], [x, y], [x, y]) == 3 * (x + y) * (x * x + y * y) for x in range(1, 9) for y in range(0, 9)))
