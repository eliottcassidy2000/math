#!/usr/bin/env python3
"""
Independent adversarial re-verification (fresh code, written cold from definitions).

TASK A — HYP-2312 (max-H => |Pf|=1) at n=6:
  Enumerate ALL 2^15 labeled tournaments on 6 vertices.
  H(T) = number of directed Hamiltonian paths (counted via subset DP).
  Find max H, collect all maximizers, group into isomorphism classes
  (canonical form = min bit-signature over all 720 relabelings).
  For each class compute det(S) with S = A - A^T (skew) and det(I + 2A),
  exact integer arithmetic; |Pf| = sqrt(det S) (n even).

TASK B — HYP-2398 mod-4 score law, fresh check:
  S = A - A^T  (S_ij = +1 if i->j, -1 if j->i, 0 diagonal)
  s_i = out-degree (score) of vertex i.
  adj(S) = exact integer adjugate via cofactors (Bareiss integer determinant).
  Odd n (5,7,9; 200 random samples each):
      adj(S)_ij == (-1)^(s_i + s_j)  (mod 4)  for ALL i,j ?
  Even n (6,8; 200 random samples each):
      u = adj(S) @ 1 ;  u_i == +(-1)^(s_i) (mod 4) universally?
      Also test global-flip alternative u_i == -(-1)^(s_i) (mod 4),
      and the opposite orientation convention S' = A^T - A.

No numpy float determinants anywhere; all arithmetic exact over Z.
"""

import random
from itertools import permutations, combinations

# ---------------------------------------------------------------- utilities

def bareiss_det(M):
    """Exact integer determinant (fraction-free Bareiss). M: list of lists of ints."""
    n = len(M)
    if n == 0:
        return 1
    M = [row[:] for row in M]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            piv = next((r for r in range(k + 1, n) if M[r][k] != 0), None)
            if piv is None:
                return 0
            M[k], M[piv] = M[piv], M[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def adjugate(M):
    """Exact integer adjugate via cofactors: adj(M)[i][j] = (-1)^(i+j) * minor(j,i)."""
    n = len(M)
    adj = [[0] * n for _ in range(n)]
    for j in range(n):          # row deleted from M
        rows = [r for r in range(n) if r != j]
        for i in range(n):      # col deleted from M
            cols = [c for c in range(n) if c != i]
            minor = [[M[r][c] for c in cols] for r in rows]
            adj[i][j] = (-1) ** (i + j) * bareiss_det(minor)
    return adj


# tournament encodings: pairs in fixed order, bit b set <=> i -> j for pair (i,j), i<j
def pairs_of(n):
    return list(combinations(range(n), 2))


def adj_matrix(bits, n, prs):
    A = [[0] * n for _ in range(n)]
    for b, (i, j) in enumerate(prs):
        if (bits >> b) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1
    return A


def ham_path_count(A, n):
    """Count directed Hamiltonian paths via DP over (subset, last vertex)."""
    out = [0] * n
    for i in range(n):
        m = 0
        for j in range(n):
            if A[i][j]:
                m |= 1 << j
        out[i] = m
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for last in range(n):
            c = row[last]
            if c == 0 or not (mask >> last) & 1:
                continue
            ext = out[last] & ~mask
            while ext:
                nb = ext & -ext
                v = nb.bit_length() - 1
                dp[mask | nb][v] += c
                ext ^= nb
    full = (1 << n) - 1
    return sum(dp[full])


def canonical(bits, n, prs, idx):
    """Canonical form: min bit-signature over all relabelings."""
    best = None
    arc = {}
    for b, (i, j) in enumerate(prs):
        if (bits >> b) & 1:
            arc[(i, j)] = 1
            arc[(j, i)] = 0
        else:
            arc[(i, j)] = 0
            arc[(j, i)] = 1
    for perm in permutations(range(n)):
        sig = 0
        for b, (i, j) in enumerate(prs):
            pi, pj = perm[i], perm[j]
            if arc[(pi, pj)]:
                sig |= 1 << b
        if best is None or sig < best:
            best = sig
    return best


# ---------------------------------------------------------------- TASK A

def task_A():
    n = 6
    prs = pairs_of(n)
    idx = {p: b for b, p in enumerate(prs)}
    N = 1 << len(prs)
    print("=" * 72)
    print("TASK A: exhaustive n=6 census of H-maximizers (fresh code)")
    print("=" * 72)
    Hmax = -1
    maximizers = []
    hist = {}
    for bits in range(N):
        A = adj_matrix(bits, n, prs)
        H = ham_path_count(A, n)
        hist[H] = hist.get(H, 0) + 1
        if H > Hmax:
            Hmax = H
            maximizers = [bits]
        elif H == Hmax:
            maximizers.append(bits)
    print(f"max H over all {N} labeled tournaments on 6 vertices: {Hmax}")
    print(f"number of labeled maximizers: {len(maximizers)}")

    classes = {}
    for bits in maximizers:
        c = canonical(bits, n, prs, idx)
        classes.setdefault(c, []).append(bits)
    print(f"number of ISO CLASSES attaining H = {Hmax}: {len(classes)}")
    for c, members in sorted(classes.items()):
        A = adj_matrix(c, n, prs)
        S = [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]
        detS = bareiss_det(S)
        I2A = [[(1 if i == j else 0) + 2 * A[i][j] for j in range(n)] for i in range(n)]
        detI2A = bareiss_det(I2A)
        scores = tuple(sorted(sum(A[i]) for i in range(n)))
        pf = round(detS ** 0.5)
        assert pf * pf == detS, f"det(S)={detS} not a perfect square!"
        assert detS == detI2A, f"det(S)={detS} != det(I+2A)={detI2A}"
        Hc = ham_path_count(A, n)
        print(f"  class canon={c:5d}: labeled count={len(members):4d}, "
              f"scores={scores}, H={Hc}, det(S)={detS}, det(I+2A)={detI2A}, |Pf|={pf}")
    return Hmax, classes


# ---------------------------------------------------------------- TASK B

def random_tournament(n, rng):
    prs = pairs_of(n)
    bits = rng.getrandbits(len(prs))
    return adj_matrix(bits, n, prs)


def task_B():
    print()
    print("=" * 72)
    print("TASK B: mod-4 score law for adj(S), S = A - A^T (fresh code)")
    print("=" * 72)
    rng = random.Random(20260611)

    # --- odd n
    for n in (5, 7, 9):
        viol = 0
        checked = 0
        for t in range(200):
            A = random_tournament(n, rng)
            s = [sum(A[i]) for i in range(n)]
            S = [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]
            adjS = adjugate(S)
            for i in range(n):
                for j in range(n):
                    want = 1 if (s[i] + s[j]) % 2 == 0 else -1
                    if (adjS[i][j] - want) % 4 != 0:
                        viol += 1
                    checked += 1
        print(f"odd n={n}: 200 random tournaments, {checked} (i,j) entries checked, "
              f"violations of adj(S)_ij == (-1)^(s_i+s_j) (mod 4): {viol}")

    # --- even n
    for n in (6, 8):
        plus_ok = 0
        minus_only = 0
        mixed = 0
        viol_entries = 0
        flip_conv_plus = 0   # opposite convention S' = A^T - A, + sign
        for t in range(200):
            A = random_tournament(n, rng)
            s = [sum(A[i]) for i in range(n)]
            S = [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]
            adjS = adjugate(S)
            u = [sum(adjS[i][j] for j in range(n)) for i in range(n)]
            okp = all((u[i] - (1 if s[i] % 2 == 0 else -1)) % 4 == 0 for i in range(n))
            okm = all((u[i] + (1 if s[i] % 2 == 0 else -1)) % 4 == 0 for i in range(n))
            if okp:
                plus_ok += 1
            elif okm:
                minus_only += 1
            else:
                mixed += 1
                viol_entries += sum(
                    1 for i in range(n)
                    if (u[i] - (1 if s[i] % 2 == 0 else -1)) % 4 != 0)
            # opposite orientation convention
            S2 = [[-S[i][j] for j in range(n)] for i in range(n)]
            adjS2 = adjugate(S2)
            u2 = [sum(adjS2[i][j] for j in range(n)) for i in range(n)]
            if all((u2[i] - (1 if s[i] % 2 == 0 else -1)) % 4 == 0 for i in range(n)):
                flip_conv_plus += 1
        print(f"even n={n}: 200 random tournaments; "
              f"u=adj(S)@1, S=A-A^T:")
        print(f"    u_i == +(-1)^(s_i) (mod 4) for ALL i:  {plus_ok}/200")
        print(f"    only global-flip  -(-1)^(s_i) works:   {minus_only}/200")
        print(f"    neither (mixed signs / true violation): {mixed}/200  "
              f"(violating entries: {viol_entries})")
        print(f"    [convention check] S'=A^T-A, + sign holds: {flip_conv_plus}/200")

    # --- bonus: exhaustive even n=6 sign check (cheap, 32768)
    n = 6
    prs = pairs_of(n)
    plus_all = True
    bad = 0
    for bits in range(1 << len(prs)):
        A = adj_matrix(bits, n, prs)
        s = [sum(A[i]) for i in range(n)]
        S = [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]
        adjS = adjugate(S)
        u = [sum(adjS[i][j] for j in range(n)) for i in range(n)]
        if not all((u[i] - (1 if s[i] % 2 == 0 else -1)) % 4 == 0 for i in range(n)):
            plus_all = False
            bad += 1
    print(f"exhaustive n=6 (all 32768): + sign universal: {plus_all} "
          f"(failures: {bad})")


if __name__ == "__main__":
    task_A()
    task_B()
