"""AUDIT: relation between Theorem Q (3) of the note and the cycle-cut bridge (2.3) of Session130
(05-knowledge/results/session130_blue_black_theorems_20261031/hamiltonian/PROOFS.md, Theorem 2):

    q = sum over u in S with u -> rho(u) of b_(u, rho(u))   (mod 2),

S = out-neighbours of the fixed vertex z, B = T - z, W_uv = number of ordered pairs of non-empty directed paths
partitioning V(B), the first ending at u and the second starting at v, b_uv = (W_uv + W_vu)/2.
In the notation of the note: q = F, the sum runs over the pairs of D_1, u = a_i^+, rho(u) = a_i^-.
Checked on every labelled tournament invariant under the fixed anti-involution, orders 3, 5, 7:
 (i)   (2.3) itself (F by direct enumeration of invariant odd cycles);
 (ii)  Theorem Q (3): F = sum over D_1 of P_Z(i), P_Z(i) = number of paths of the fold Z ending at i;
 (iii) whether the two weights agree termwise modulo 2: b_(a_i^+, a_i^-) = P_Z(i) for every pair i in D_1.
"""
import itertools, sys


def all_antiinvariant(m):
    n = 2 * m + 1
    psi = [0] + [0] * (2 * m)
    for i in range(m):
        psi[1 + 2 * i] = 2 + 2 * i
        psi[2 + 2 * i] = 1 + 2 * i
    free = []
    for i in range(m):
        free.append(("z", i))
        free.append(("in", i))
    for i in range(m):
        for j in range(i + 1, m):
            free.append(("pp", i, j))
            free.append(("pm", i, j))
    for bits in range(1 << len(free)):
        A = [[0] * n for _ in range(n)]
        for k, f in enumerate(free):
            b = bits >> k & 1
            if f[0] == "z":
                a, c = 1 + 2 * f[1], 2 + 2 * f[1]
                if b:
                    A[0][c] = 1; A[a][0] = 1
                else:
                    A[0][a] = 1; A[c][0] = 1
            elif f[0] == "in":
                a, c = 1 + 2 * f[1], 2 + 2 * f[1]
                if b:
                    A[c][a] = 1
                else:
                    A[a][c] = 1
            elif f[0] == "pp":
                a, c = 1 + 2 * f[1], 1 + 2 * f[2]
                if b:
                    A[a][c] = 1; A[psi[c]][psi[a]] = 1
                else:
                    A[c][a] = 1; A[psi[a]][psi[c]] = 1
            else:
                a, c = 1 + 2 * f[1], 2 + 2 * f[2]
                if b:
                    A[a][c] = 1; A[psi[c]][psi[a]] = 1
                else:
                    A[c][a] = 1; A[psi[a]][psi[c]] = 1
        yield A, psi


def invariant_odd_cycles(A, psi):
    n = len(A)
    cnt = 0
    for L in range(3, n + 1, 2):
        for S in itertools.combinations(range(n), L):
            v0 = S[0]
            for p in itertools.permutations(S[1:]):
                seq = (v0,) + p
                if all(A[seq[i]][seq[(i + 1) % L]] for i in range(L)):
                    succ = {seq[i]: seq[(i + 1) % L] for i in range(L)}
                    if all(succ.get(psi[w]) == psi[u] for u, w in succ.items()):
                        cnt += 1
    return cnt


def path_counts(A, verts):
    """end[S][v], start[S][v]: Hamiltonian paths of the sub-tournament on the vertex set S (bitmask over verts)."""
    k = len(verts)
    end = [[0] * k for _ in range(1 << k)]
    for i in range(k):
        end[1 << i][i] = 1
    for S in range(1, 1 << k):
        for i in range(k):
            x = end[S][i]
            if x:
                for j in range(k):
                    if not S >> j & 1 and A[verts[i]][verts[j]]:
                        end[S | 1 << j][j] += x
    start = [[0] * k for _ in range(1 << k)]
    for i in range(k):
        start[1 << i][i] = 1
    for S in range(1, 1 << k):
        for i in range(k):
            x = start[S][i]
            if x:
                for j in range(k):
                    if not S >> j & 1 and A[verts[j]][verts[i]]:
                        start[S | 1 << j][j] += x
    return end, start


for m in (1, 2, 3):
    tot = bad23 = badQ = termwise_bad = terms = 0
    for A, psi in all_antiinvariant(m):
        n = 2 * m + 1
        tot += 1
        F = invariant_odd_cycles(A, psi)
        verts = list(range(1, n))
        k = len(verts)
        end, start = path_counts(A, verts)
        full = (1 << k) - 1
        def W(u, v):
            iu, iv = verts.index(u), verts.index(v)
            s = 0
            for S in range(1, full):
                if S >> iu & 1 and not S >> iv & 1:
                    s += end[S][iu] * start[full ^ S][iv]
            return s
        plus = [a for a in verts if A[0][a]]
        minus = [psi[a] for a in plus]
        D1 = [A[plus[i]][minus[i]] for i in range(m)]
        Z = [[1 if (i != j and (A[plus[i]][plus[j]] + A[plus[i]][minus[j]]) == 1) else 0 for j in range(m)] for i in range(m)]
        # paths of Z ending at i
        P = [0] * m
        def rec(last, used):
            P[last] += 1
            for j in range(m):
                if not used >> j & 1 and Z[last][j]:
                    rec(j, used | 1 << j)
        for i in range(m):
            rec(i, 1 << i)
        s130 = 0
        sQ = 0
        for i in range(m):
            if D1[i]:
                w1, w2 = W(plus[i], minus[i]), W(minus[i], plus[i])
                assert (w1 + w2) % 2 == 0
                b = (w1 + w2) // 2
                s130 += b
                sQ += P[i]
                terms += 1
                if (b - P[i]) % 2:
                    termwise_bad += 1
        if (F - s130) % 2:
            bad23 += 1
        if (F - sQ) % 2:
            badQ += 1
    print("order %d: %d tournaments; violations of Session130 (2.3): %d; of Theorem Q (3): %d; pairs in D_1 examined: %d, with b != P_Z (mod 2): %d" % (
        2 * m + 1, tot, bad23, badQ, terms, termwise_bad))
