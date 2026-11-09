#!/usr/bin/env python3
"""AUDIT (independent): port-chain block matrices by a BLACK-BOX linear solve, no case analysis.

For a block (B, w, u') we glue it after several small ported tournaments (T, u) (every T-vertex beats every
B-vertex except w -> u; new port u'), compute the old and new states (h, E, O) by brute force (DFS cycle
enumeration + memoised packing polynomial + DFS path enumeration), solve N from three test tournaments with
independent states, and then CHECK N on further test tournaments (this tests that a transfer matrix exists at all).
Everything is done with packing polynomials alpha(T) = (alpha_0, alpha_1, ...), so all bases x = 2..8 at once.

Writes OUTDIR/mats_x{X}.txt: lines "kmin | 9 entries | k rows w u'" (one line per distinct matrix, kmin = least
block order realising it), and prints counts of distinct matrices for k <= 4 and k <= 5.
Usage: a1_blocks.py OUTDIR [KMAX]
"""
import itertools
import sys
from fractions import Fraction

XS = range(2, 9)


def odd_cycle_masks(n, adj):
    """dict mask -> number of directed odd cycles (length >= 3) on exactly that vertex set (DFS from min)."""
    cnt = {}
    for s in range(n):
        stack = [(s, 1 << s, 1)]
        while stack:
            v, mask, ln = stack.pop()
            if ln >= 3 and ln % 2 == 1 and adj[v][s]:
                cnt[mask] = cnt.get(mask, 0) + 1
            for t in range(s + 1, n):
                if not mask >> t & 1 and adj[v][t]:
                    stack.append((t, mask | 1 << t, ln + 1))
    return cnt


def packing_poly_factory(n, adj):
    cnt = odd_cycle_masks(n, adj)
    bymin = {}
    for m, c in cnt.items():
        low = (m & -m).bit_length() - 1
        bymin.setdefault(low, []).append((m, c))
    memo = {0: (1,)}

    def add(p, q, c=1, shift=0):
        r = list(p) + [0] * max(0, len(q) + shift - len(p))
        for i, a in enumerate(q):
            r[i + shift] += c * a
        return tuple(r)

    def P(S):
        if S in memo:
            return memo[S]
        low = (S & -S).bit_length() - 1
        res = P(S & ~(1 << low))
        for (m, c) in bymin.get(low, []):
            if m & ~S == 0:
                res = add(res, P(S & ~m), c, 1)
        memo[S] = res
        return res

    return P


def state_polys(n, adj, port):
    """(h, E, O) as packing polynomials (tuples of coefficients)."""
    P = packing_poly_factory(n, adj)
    full = (1 << n) - 1

    def add(p, q):
        r = list(p) + [0] * max(0, len(q) - len(p))
        for i, a in enumerate(q):
            r[i] += a
        return tuple(r)

    E, O = (0,), (0,)
    stack = [(port, 1 << port, 1)]
    while stack:
        v, mask, ln = stack.pop()
        p = P(full & ~mask)
        if ln % 2:
            O = add(O, p)
        else:
            E = add(E, p)
        for t in range(n):
            if not mask >> t & 1 and adj[v][t]:
                stack.append((t, mask | 1 << t, ln + 1))
    return P(full), E, O


def ev(p, x):
    return sum(a * x ** i for i, a in enumerate(p))


def glue(T, u, B, w, up):
    n, k = len(T), len(B)
    A = [[0] * (n + k) for _ in range(n + k)]
    for i in range(n):
        for j in range(n):
            A[i][j] = T[i][j]
        for j in range(k):
            A[i][n + j] = 1
    for i in range(k):
        for j in range(k):
            A[n + i][n + j] = B[i][j]
    A[u][n + w] = 0
    A[n + w][u] = 1
    return A, n + up


def tour_from_bits(k, bits):
    A = [[0] * k for _ in range(k)]
    t = 0
    for i in range(k):
        for j in range(i + 1, k):
            if bits >> t & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            t += 1
    return A


def canon(A):
    k = len(A)
    best = None
    for p in itertools.permutations(range(k)):
        code = tuple(A[p[i]][p[j]] for i in range(k) for j in range(k))
        if best is None or code < best:
            best = code
    return best


def classes(k):
    seen = {}
    for bits in range(1 << (k * (k - 1) // 2)):
        A = tour_from_bits(k, bits)
        c = canon(A)
        if c not in seen:
            seen[c] = A
    return list(seen.values())


# test ported tournaments: the first three have independent states for every x >= 2 (det = 1 - x)
TESTS = [
    ([[0]], 0),
    ([[0, 1], [0, 0]], 0),
    ([[0, 1, 0], [0, 0, 1], [1, 0, 0]], 0),
    ([[0, 0], [1, 0]], 0),                      # same state as the first, different tournament
    ([[0, 1, 1], [0, 0, 1], [0, 0, 0]], 0),     # transitive, port = source
    ([[0, 1, 1], [0, 0, 1], [0, 0, 0]], 2),     # transitive, port = sink
    (tour_from_bits(4, 0b101101), 1),
    (tour_from_bits(4, 0b011010), 3),
]


def solve3(S, Sp):
    """N with N S = Sp, S and Sp 3x3 (columns = states), exact rationals."""
    M = [[Fraction(S[r][c]) for c in range(3)] for r in range(3)]
    # inverse of M by Gauss-Jordan
    inv = [[Fraction(int(r == c)) for c in range(3)] for r in range(3)]
    for col in range(3):
        piv = next(r for r in range(col, 3) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        inv[col], inv[piv] = inv[piv], inv[col]
        d = M[col][col]
        M[col] = [a / d for a in M[col]]
        inv[col] = [a / d for a in inv[col]]
        for r in range(3):
            if r != col and M[r][col] != 0:
                f = M[r][col]
                M[r] = [a - f * b for a, b in zip(M[r], M[col])]
                inv[r] = [a - f * b for a, b in zip(inv[r], inv[col])]
    return [[sum(Fraction(Sp[r][i]) * inv[i][c] for i in range(3)) for c in range(3)] for r in range(3)]


def main():
    outdir = sys.argv[1]
    kmax = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    test_states = [state_polys(len(T), T, u) for (T, u) in TESTS]
    mats = {x: {} for x in XS}
    nblocks = 0
    bad = 0
    allf = open(f"{outdir}/blocks_all.txt", "w")
    for k in range(1, kmax + 1):
        for B in classes(k):
            for w in range(k):
                for up in range(k):
                    nblocks += 1
                    rows = "".join(str(B[i][j]) for i in range(k) for j in range(k))
                    line = [f"{k} {rows} {w} {up}"]
                    new = []
                    for (T, u) in TESTS:
                        A, port = glue(T, u, B, w, up)
                        new.append(state_polys(len(A), A, port))
                    for x in XS:
                        old = [[ev(p, x) for p in s] for s in test_states]
                        nw = [[ev(p, x) for p in s] for s in new]
                        S = [[old[c][r] for c in range(3)] for r in range(3)]
                        Sp = [[nw[c][r] for c in range(3)] for r in range(3)]
                        N = solve3(S, Sp)
                        ok = all(a.denominator == 1 and a >= 0 for row in N for a in row)
                        for t in range(3, len(TESTS)):
                            pred = [sum(N[r][c] * old[t][c] for c in range(3)) for r in range(3)]
                            if pred != nw[t]:
                                ok = False
                        if not ok:
                            bad += 1
                            print("NOT A NONNEGATIVE INTEGER TRANSFER", k, B, w, up, x, N)
                            continue
                        key = tuple(int(a) for row in N for a in row)
                        line.append(" ".join(map(str, key)))
                        if key not in mats[x]:
                            mats[x][key] = (k, f"{k} {rows} {w} {up}")
                    allf.write(" | ".join(line) + "\n")
        print(f"k <= {k}: blocks {nblocks}; distinct matrices per base x=2..8:",
              [len(mats[x]) for x in XS], flush=True)
    allf.close()
    for x in XS:
        with open(f"{outdir}/mats_x{x}.txt", "w") as f:
            for key in sorted(mats[x]):
                kmin, rep = mats[x][key]
                f.write(f"{kmin} | {' '.join(map(str, key))} | {rep}\n")
        k4 = sum(1 for v in mats[x].values() if v[0] <= 4)
        print(f"x={x}: distinct matrices k<=4: {k4}, k<={kmax}: {len(mats[x])}")
    print("linearity/integrality/nonnegativity failures:", bad)


if __name__ == "__main__":
    main()
