"""Verifier of THM-5135 (S92 audit follow-up): Remark 2 (two kernels with 3 x 3 visit matrices, determinants 32
and -6, and the two 2 x 2 visit matrices quoted in (4)) and Remark 3 (collapse modulo 3 of every such kernel).
Own code: walks are enumerated directly; H(K[A, B]) is counted on the substituted tournament.
Run: python v5_visit_matrices_rank3.py > v5_visit_matrices_rank3.out
"""
import itertools
import random

FAILS = []


def check(cond, msg):
    print(("  ok   " if cond else "  FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def from_bits(m, bits):
    """pairs (i, j), i < j, in lexicographic order; bit k set: i -> j, otherwise j -> i (the audit's encoding)."""
    adj = [[0] * m for _ in range(m)]
    for k, (i, j) in enumerate((i, j) for i in range(m) for j in range(i + 1, m)):
        if (bits >> k) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj


def ham_table(adj):
    n = len(adj)
    full = 1 << n
    dp = [[0] * n for _ in range(full)]
    hp = [0] * full
    hp[0] = 1
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, full):
        tot = 0
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if c:
                tot += c
                for w in range(n):
                    if adj[v][w] and not (mask >> w) & 1:
                        dp[mask | (1 << w)][w] += c
        hp[mask] = tot
    return hp


def cover_numbers(adj, imax):
    """P_i = number of ordered i-tuples of disjoint non-empty paths covering the tournament, i = 1..imax."""
    n = len(adj)
    hp = ham_table(adj)
    full = (1 << n) - 1
    # f[i][mask] = number of ordered i-tuples of disjoint non-empty paths covering exactly mask
    f = [dict() for _ in range(imax + 1)]
    f[0] = {0: 1}
    for i in range(1, imax + 1):
        cur = {}
        for mask, c in f[i - 1].items():
            rest = full & ~mask
            sub = rest
            while sub:
                cur[mask | sub] = cur.get(mask | sub, 0) + c * hp[sub]
                sub = (sub - 1) & rest
        f[i] = cur
    return [0] + [f[i].get(full, 0) for i in range(1, imax + 1)]


def visit_matrix(adj, x, y, vmax):
    """N[i][j] = number of walks visiting every vertex other than x, y exactly once, x exactly i and y exactly j times."""
    m = len(adj)
    others = [v for v in range(m) if v not in (x, y)]
    N = [[0] * (vmax + 2) for _ in range(vmax + 2)]

    def rec(v, used, cx, cy):
        if len(used) == len(others):
            N[cx][cy] += 1
        for w in range(m):
            if adj[v][w]:
                if w == x:
                    if cx <= vmax:
                        rec(w, used, cx + 1, cy)
                elif w == y:
                    if cy <= vmax:
                        rec(w, used, cx, cy + 1)
                elif w not in used:
                    rec(w, used | {w}, cx, cy)

    for v in range(m):
        if v == x:
            rec(v, frozenset(), 1, 0)
        elif v == y:
            rec(v, frozenset(), 0, 1)
        else:
            rec(v, frozenset([v]), 0, 0)
    return N


def substitute(adjK, x, y, A, B):
    m = len(adjK)
    blocks = []
    for v in range(m):
        if v == x:
            blocks.append([("A", i) for i in range(len(A))])
        elif v == y:
            blocks.append([("B", i) for i in range(len(B))])
        else:
            blocks.append([("K", v)])
    verts = [u for blk in blocks for u in blk]
    owner = {}
    for v, blk in enumerate(blocks):
        for u in blk:
            owner[u] = v
    n = len(verts)
    adj = [[0] * n for _ in range(n)]
    for a, u in enumerate(verts):
        for b, w in enumerate(verts):
            if a == b:
                continue
            if owner[u] == owner[w]:
                T = A if u[0] == "A" else B
                adj[a][b] = T[u[1]][w[1]]
            else:
                adj[a][b] = adjK[owner[u]][owner[w]]
    return adj


def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))


# modules: all tournaments of orders 1, 2, 3 up to isomorphism, and three of order 4
MODULES = [from_bits(1, 0), from_bits(2, 1), from_bits(3, 7), from_bits(3, 5), from_bits(4, 63), from_bits(4, 0b101011), from_bits(4, 0b010110)]
MODINFO = []
for T in MODULES:
    P = cover_numbers(T, 4)
    MODINFO.append((T, P))
print("modules (order, P_1 = H, P_2, P_3, P_4):", [(len(T), P[1], P[2], P[3], P[4]) for T, P in MODINFO])
check(all(P[2] % 2 == 0 and P[3] % 6 == 0 and P[4] % 24 == 0 for T, P in MODINFO), "P_m is divisible by m! (m = 2, 3, 4) on these modules")
check(all((P[2] - ((-1) ** len(T) + 1) * P[1]) % 3 == 0 for T, P in MODINFO), "P_2 = ((-1)^n + 1) H (mod 3) on these modules")


def test_kernel(name, adjK, x, y, want=None, want_det=None, vmax=3):
    m = len(adjK)
    N = visit_matrix(adjK, x, y, vmax)
    more = any(N[i][j] for i in range(vmax + 2) for j in range(vmax + 2) if i > vmax or j > vmax)
    G = [[N[i][j] for j in range(1, vmax + 1)] for i in range(1, vmax + 1)]
    ok = True
    consts = {}
    okc = True
    for (A, PA) in MODINFO:
        for (B, PB) in MODINFO:
            if m - 2 + len(A) + len(B) > 10:
                continue
            Hsub = ham_table(substitute(adjK, x, y, A, B))[-1]
            pred = sum(N[i][j] * PA[i] * PB[j] for i in range(1, vmax + 2) for j in range(1, vmax + 2) if i < len(PA) and j < len(PB))
            ok &= (Hsub == pred)
            key = (len(A) % 2, len(B) % 2)
            hh = (PA[1] * PB[1]) % 3
            if hh:
                c = (Hsub * pow(hh, -1, 3)) % 3
                if key in consts and consts[key] != c:
                    okc = False
                consts.setdefault(key, c)
            else:
                okc &= (Hsub % 3 == 0)
    msg = f"{name}: visit matrix {G}"
    if vmax == 3:
        msg += f", determinant {det3(G)}"
    else:
        msg += f", determinant {G[0][0] * G[1][1] - G[0][1] * G[1][0]}"
    good = ok and not more
    if want is not None:
        good &= (G == want)
    if want_det is not None:
        good &= ((det3(G) if vmax == 3 else G[0][0] * G[1][1] - G[0][1] * G[1][0]) == want_det)
    check(good, msg + f"; no walk with more visits: {not more}; H(K[A, B]) = sum N_ij P_i(A) P_j(B) on the modules: {ok}")
    check(okc, f"   modulo 3: H(K[A, B]) = c H(A) H(B) with c depending only on the parities of |A|, |B|: constants {dict(sorted(consts.items()))}")
    return G


print("[1] Remark 2: the two kernels with three visits (encoding of the audit: order 4 bits 4 ports (0, 3); order 5 bits 76 ports (0, 1))")
test_kernel("order 4", from_bits(4, 4), 0, 3, want=[[5, 3, 0], [3, 7, 2], [0, 2, 2]], want_det=32)
test_kernel("order 5", from_bits(5, 76), 0, 1, want=[[15, 13, 2], [13, 14, 4], [2, 4, 2]], want_det=-6)

print("[2] (4): kernels of orders 4 and 5 with the 2 x 2 visit matrices [[5, 2], [2, 1]] and [[11, 6], [6, 3]] (searched here)")
found = {}
for m, target in ((4, [[5, 2], [2, 1]]), (5, [[11, 6], [6, 3]])):
    for bits in range(1 << (m * (m - 1) // 2)):
        adjK = from_bits(m, bits)
        hit = None
        for x in range(m):
            for y in range(m):
                if x != y:
                    N = visit_matrix(adjK, x, y, 2)
                    if [[N[1][1], N[1][2]], [N[2][1], N[2][2]]] == target and not any(
                            N[i][j] for i in range(4) for j in range(4) if i > 2 or j > 2):
                        hit = (bits, x, y)
                        break
            if hit:
                break
        if hit:
            found[m] = hit
            break
check(4 in found and 5 in found, f"found (bits, x, y): {found}")
for m, target in ((4, [[5, 2], [2, 1]]), (5, [[11, 6], [6, 3]])):
    bits, x, y = found[m]
    test_kernel(f"order {m}, bits {bits}, ports ({x}, {y})", from_bits(m, bits), x, y, want=target, want_det=(1 if m == 4 else -3), vmax=2)

print("[3] Remark 3: random kernels of orders 3, 4, 5 with every number of visits up to 4")
random.seed(3)
allok = True
cnt = 0
for m in (3, 4, 5):
    for _ in range(12 if m < 5 else 8):
        bits = random.randrange(1 << (m * (m - 1) // 2))
        adjK = from_bits(m, bits)
        x, y = random.sample(range(m), 2)
        N = visit_matrix(adjK, x, y, 4)
        consts = {}
        for (A, PA) in MODINFO:
            for (B, PB) in MODINFO:
                if m - 2 + len(A) + len(B) > 9:
                    continue
                Hsub = ham_table(substitute(adjK, x, y, A, B))[-1]
                pred = sum(N[i][j] * PA[i] * PB[j] for i in range(1, 5) for j in range(1, 5))
                allok &= (Hsub == pred)
                key = (len(A) % 2, len(B) % 2)
                hh = (PA[1] * PB[1]) % 3
                if hh:
                    c = (Hsub * pow(hh, -1, 3)) % 3
                    if key in consts and consts[key] != c:
                        allok = False
                    consts.setdefault(key, c)
                else:
                    allok &= (Hsub % 3 == 0)
        cnt += 1
check(allok, f"{cnt} random (kernel, ports): substitution formula exact and the collapse modulo 3 holds")

print("ALL CHECKS PASSED" if not FAILS else f"{len(FAILS)} FAILURES: {FAILS}")
