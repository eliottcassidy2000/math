# Literal check of a compensated (line) bank: five seeds with h = 41, G_e + G_o = 1424 whose signatures (a, b) lie on a
# line of direction (1, -2); for every exterior with alpha = 2 beta (and u1 = z1 = 0) the five counts are consecutive odd numbers.
from master_formula_check import block_stats, brute
def mat(masks):
    m = len(masks); return [[1 if masks[u] >> v & 1 else 0 for v in range(m)] for u in range(m)]
B1 = [0, 5, 17, 7, 107, 79, 143, 63]; B2 = [2, 4, 25, 19, 35, 79, 31, 127]
seeds = [(B1, 187), (B1, 219), (B1, 69), (B1, 37), (B2, 159)]
for masks, X in seeds:
    C = mat(masks); m = 8; h, W, G, lam, kap = block_stats(C, m)
    a = sum(lam[s][0] + lam[s][1] for s in range(m) if not X >> s & 1); b = sum(kap[s][0] + kap[s][1] for s in range(m) if X >> s & 1)
    cut = sum(W[u][v][0] for u in range(m) for v in range(m) if u != v and (X >> u & 1) and not (X >> v & 1))
    print("seed", masks, X, ": h =", h[-1], "G =", G[0] + G[1], "(a, b) =", (a, b), "cut =", cut, " q a + p b =", 2 * a + b)
def tapes(alpha, beta):
    """upper tape with u1 = 0 and given alpha; lower tape with z1 = 0 (depth 1 a contact) and given beta."""
    U = set(); hh = 2; x = alpha
    while x: 
        if x & 1: U.add(hh)
        x >>= 1; hh += 1
    i = max(U) if U else 0
    Z = set(); d = 2; x = beta
    while x:
        if x & 1: Z.add(d)
        x >>= 1; d += 1
    j = max(Z) if Z else 0
    L = {d for d in range(1, j + 1) if d not in Z}
    return i, j, U, L
ok = True
for beta in range(0, 5):
    alpha = 2 * beta; i, j, U, L = tapes(alpha, beta)
    if i + j + 9 > 17: continue
    vals = [brute(mat(masks), 8, i, j, U, X, L) for masks, X in seeds]
    good = all(vals[k + 1] - vals[k] == 2 for k in range(4)); ok &= good
    print("alpha = %d, beta = %d (tapes %d + %d, order %d): counts %s %s" % (alpha, beta, i, j, i + j + 9, vals, "consecutive" if good else "NOT consecutive"))
# off the line the bank breaks
i, j, U, L = tapes(1, 1); print("off the line (alpha, beta) = (1, 1):", [brute(mat(masks), 8, i, j, U, X, L) for masks, X in seeds])
print("ALL OK" if ok else "FAILURE")
