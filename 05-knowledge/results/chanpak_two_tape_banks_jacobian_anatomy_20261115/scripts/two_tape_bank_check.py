# Literal check of the two-tape three-count bank.
# Seeds (block masks, head pattern):  S0 = ([22,44,24,49,34,5], 23),  S1 = ([38,28,24,49,33,6], 39),  S2 = ([38,28,24,49,33,6], 23).
# Claim: for every upper tape word and every lower tape word the three tournaments have counts K+90, K+92, K+94 with
#   K = 432 A B + 174 A + 303 B + 45 - 30 u1 - 15 z1 - 90 (J(A) + J(B)),
#   A = sum over upper contacts of 2^(height-1), B = sum over lower NON-contacts of 2^(depth-1), u1 = A mod 2, z1 = B mod 2,
#   J = pair energy of the binary digits: sum over pairs of set bits at distance g of c(g), c(1) = 1, c(g) = 2^(g-2).
import random
from master_formula_check import block_stats, brute, formula
def mat(masks):
    m = len(masks); return [[1 if masks[a] >> b & 1 else 0 for b in range(m)] for a in range(m)]
seeds = [([22, 44, 24, 49, 34, 5], 23), ([38, 28, 24, 49, 33, 6], 39), ([38, 28, 24, 49, 33, 6], 23)]
for masks, X in seeds:
    C = mat(masks); m = 6; h, W, G, lam, kap = block_stats(C, m)
    ao = sum(lam[s][1] for s in range(m) if not X >> s & 1); ae = sum(lam[s][0] for s in range(m) if not X >> s & 1)
    bo = sum(kap[s][1] for s in range(m) if X >> s & 1); be = sum(kap[s][0] for s in range(m) if X >> s & 1)
    cut = sum(W[a][b][0] for a in range(m) for b in range(m) if a != b and (X >> a & 1) and not (X >> b & 1))
    print("seed", masks, X, ": h =", h[63], "G =", G, "a =", (ao, ae), "b =", (bo, be), "cut =", cut, "; H(C + x) =", h[63] + 2 * cut)
def J(n):
    b = [i for i in range(n.bit_length()) if n >> i & 1]
    return sum(1 if y - x == 1 else 2 ** (y - x - 2) for i, x in enumerate(b) for y in b[i + 1:])
rng = random.Random(5); bad = 0; tests = 0
for trial in range(40):
    i = rng.randint(0, 4); j = rng.randint(0, 4)
    U = {hh for hh in range(1, i + 1) if rng.random() < 0.5}; L = {d for d in range(1, j + 1) if rng.random() < 0.5}
    A = sum(2 ** (hh - 1) for hh in U); B = sum(2 ** (d - 1) for d in range(1, j + 1) if d not in L)
    K = 432 * A * B + 174 * A + 303 * B + 45 - 30 * (A & 1) - 15 * (B & 1) - 90 * (J(A) + J(B))
    vals = [brute(mat(masks), 6, i, j, U, X, L) for masks, X in seeds]
    tests += 1
    if vals != [K + 90, K + 92, K + 94]: bad += 1; print("MISMATCH", i, j, U, L, vals, K)
print("two-tape bank: %d random tape pairs checked by brute force (orders up to 15), mismatches: %d" % (tests, bad))
print("example: A = 5, B = 6 ->", 432 * 30 + 174 * 5 + 303 * 6 + 45 - 30 - 90 * (J(5) + J(6)) + 90, "and the next two odd numbers")
