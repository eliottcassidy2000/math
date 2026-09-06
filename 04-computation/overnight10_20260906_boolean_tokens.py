#!/usr/bin/env python3
"""Exact native Boolean triangle sectors and their path-token spectrum.

No repository imports, floating point, or external packages. Literal triangle
XOR/BFS, integer characteristic polynomials, and cyclotomic subset-sum products
are separate finite controls for the companion all-parameter proof.
"""
from collections import Counter
from itertools import combinations, combinations_with_replacement
from math import comb, factorial
import sys

sys.stdout.reconfigure(newline="\n")
GATES = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def trim(a):
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def pdiv(a, b):
    a = a[:]
    q = [0] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b) and a != [0]:
        shift = len(a) - len(b)
        x, rem = divmod(a[-1], b[-1])
        need(rem == 0, "integral monic polynomial division")
        q[shift] = x
        for j, y in enumerate(b):
            a[j + shift] -= x * y
        trim(a)
    return trim(q), a


PHI = {}


def cyclotomic(m):
    if m not in PHI:
        a = [-1] + [0] * (m - 1) + [1]
        for d in range(1, m):
            if m % d == 0:
                a, rem = pdiv(a, cyclotomic(d))
                need(rem == [0], "cyclotomic exact factorization")
        PHI[m] = a
    return PHI[m]


def ring_reduce(a, phi):
    a = a[:] + [0] * max(0, len(phi) - 1 - len(a))
    for j in range(len(a) - 1, len(phi) - 2, -1):
        x = a[j]
        for k in range(len(phi)):
            a[j - len(phi) + 1 + k] -= x * phi[k]
    return tuple(a[:len(phi) - 1])


def ring_add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def ring_mul(a, b, phi):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return ring_reduce(out, phi)


def spectral_polynomial(N, c):
    """Product over sums of c distinct 2*cos(pi*j/(N+1)), exactly."""
    h = 2 * (N + 1)
    phi = cyclotomic(h)
    zero = (0,) * (len(phi) - 1)
    one = (1,) + zero[1:]
    modes = {}
    for j in range(1, N + 1):
        x = [0] * h
        x[j] += 1
        x[h - j] += 1
        modes[j] = ring_reduce(x, phi)
    char = [one]  # increasing powers of the eigenvalue indeterminate
    zero_sets = []
    for S in combinations(range(1, N + 1), c):
        value = zero
        for j in S:
            value = ring_add(value, modes[j])
        if value == zero:
            zero_sets.append(S)
        out = [zero] * (len(char) + 1)
        for i, x in enumerate(char):
            out[i + 1] = ring_add(out[i + 1], x)
            minus = tuple(-y for y in ring_mul(x, value, phi))
            out[i] = ring_add(out[i], minus)
        char = out
    need(all(all(x == 0 for x in a[1:]) for a in char),
         "cyclotomic characteristic polynomial descends to Z")
    return [a[0] for a in reversed(char)], zero_sets


def matrix_charpoly(A):
    """Integer Faddeev-LeVerrier recurrence; exact division gates retained."""
    n = len(A)
    B = [[int(i == j) for j in range(n)] for i in range(n)]
    rows = [[(j, x) for j, x in enumerate(row) if x] for row in A]
    out = [1]
    for k in range(1, n + 1):
        AB = [[sum(x * B[j][z] for j, x in rows[i])
               for z in range(n)] for i in range(n)]
        ck, rem = divmod(-sum(AB[i][i] for i in range(n)), k)
        need(rem == 0, "integer characteristic recurrence")
        out.append(ck)
        for i in range(n):
            AB[i][i] += ck
        B = AB
    need(all(x == 0 for row in B for x in row), "Cayley-Hamilton residual")
    return out


def token_matrix(N, c):
    states = list(combinations(range(1, N + 1), c))
    ids = {x: i for i, x in enumerate(states)}
    A = [[0] * len(states) for _ in states]
    for i, S in enumerate(states):
        for x in S:
            for y in (x - 1, x + 1):
                if 1 <= y <= N and y not in S:
                    T = tuple(sorted(set(S) - {x} | {y}))
                    A[i][ids[T]] += 1
    need(all(x in (0, 1) for row in A for x in row), "token edges are Boolean")
    return states, A


def graph_from_lengths(n, lengths):
    edges = set()
    start = 0
    for m in lengths:
        vertices = list(range(start, start + m))
        for i in range(m):
            edges.add(tuple(sorted((vertices[i], vertices[(i + 1) % m]))))
        start += m
    need(start <= n, "literal representative fits vertex budget")
    return edges


def cycle_type(n, edges):
    adjacency = [set() for _ in range(n)]
    for u, v in edges:
        adjacency[u].add(v)
        adjacency[v].add(u)
    if any(len(row) not in (0, 2) for row in adjacency):
        return None
    seen = set()
    sizes = []
    for u in range(n):
        if u in seen or not adjacency[u]:
            continue
        seen.add(u)
        queue = [u]
        for x in queue:
            for y in adjacency[x]:
                if y not in seen:
                    seen.add(y)
                    queue.append(y)
        sizes.append(len(queue))
    return tuple(sorted(sizes))


def orbit_size(n, S):
    den = factorial(n - sum(S))
    for s, multiplicity in Counter(S).items():
        den *= (2 * s) ** multiplicity * factorial(multiplicity)
    value, rem = divmod(factorial(n), den)
    need(rem == 0, "orbit volume integer")
    return value


def literal_sector(n, c, M):
    states = [S for S in combinations_with_replacement(range(3, M + 1), c)
              if sum(S) <= n]
    ids = {x: i for i, x in enumerate(states)}
    matrix = [[0] * len(states) for _ in states]
    trials = 0
    for i, S in enumerate(states):
        F = graph_from_lengths(n, S)
        need(cycle_type(n, F) == S, "literal initial cycle classification")
        for triangle in combinations(range(n), 3):
            H = F ^ set(combinations(triangle, 2))
            T = cycle_type(n, H)
            trials += 1
            if T in ids:
                matrix[i][ids[T]] += 1
        expected = Counter()
        for s, multiplicity in Counter(S).items():
            for change in (-1, 1):
                t = s + change
                Tlist = list(S)
                Tlist.remove(s)
                Tlist.append(t)
                T = tuple(sorted(Tlist))
                if T not in ids:
                    continue
                count = multiplicity * s
                if change == 1:
                    count *= n - sum(S)
                expected[T] += count
        need(matrix[i] == [expected[T] for T in states],
             "literal multiplicity versus insertion/deletion classification")
    boolean = [[int(x > 0) for x in row] for row in matrix]
    volumes = [orbit_size(n, S) for S in states]
    N = c + M - 3
    token_states, A = token_matrix(N, c)
    token_ids = {x: i for i, x in enumerate(token_states)}
    mapped = [tuple(s - 3 + i for i, s in enumerate(S, 1)) for S in states]
    charge_cap = n - 3 * c + c * (c + 1) // 2
    need(set(mapped) == {S for S in token_states if sum(S) <= charge_cap},
         "exact vertex-budget cap in token coordinates")
    for i, S in enumerate(mapped):
        for j, T in enumerate(mapped):
            need(boolean[i][j] == A[token_ids[S]][token_ids[T]],
                 "literal Boolean sector equals induced token cap")
            need(volumes[i] * matrix[i][j] == volumes[j] * matrix[j][i],
                 "literal multiplicity retains detailed balance")
    need(trials == len(states) * comb(n, 3), "literal triangle universe")
    return states, boolean, matrix, trials


def main():
    print("NATIVE BOOLEAN TRIANGLE SECTORS: EXACT CONTROLS")
    total_trials = 0
    for n, c, M in [(6, 1, 6), (8, 2, 4), (10, 2, 5),
                    (12, 3, 4), (24, 3, 8), (7, 2, 4), (9, 2, 5)]:
        states, A, weighted, trials = literal_sector(n, c, M)
        total_trials += trials
        cap = "inactive" if n >= c * M else "active"
        print(f"literal n={n} c={c} M={M}: states={len(states)} "
              f"triangles={trials} cap={cap} PASS")
        if (n, c, M) == (8, 2, 4):
            need(states == [(3, 3), (3, 4), (4, 4)], "minimal repeated-cycle order")
            need(weighted == [[0, 12, 0], [4, 0, 3], [0, 8, 0]],
                 "minimal Boolean/multiplicity mismatch")
            print("  multiplicity matrix=" + repr(weighted))
        if (n, c, M) == (7, 2, 4):
            need(len(states) == 2 and A == [[0, 1], [1, 0]],
                 "resource-cap hostile: P2 instead of P3")
    print(f"literal complete universe: {total_trials} triangle toggles")

    cases = 0
    for N in range(1, 9):
        for c in range(1, N + 1):
            states, A = token_matrix(N, c)
            direct = matrix_charpoly(A)
            spectral, zero_sets = spectral_polynomial(N, c)
            need(direct == spectral, "matrix versus cyclotomic spectral product")
            nullity = 0
            for x in reversed(direct):
                if x:
                    break
                nullity += 1
            need(nullity == len(zero_sets), "exact spectral zero multiplicity")
            need(len(states) == comb(N, c), "token state universe")
            cases += 1
            if (N, c) == (4, 2):
                need(direct == [1, 0, -6, 0, 5, 0, 0], "six-state exact spectrum")
                print("c=2,M=5: chi=lambda^2*(lambda^2-1)*(lambda^2-5)")
            if (N, c) == (8, 3):
                parity = Counter(sum(S) % 2 for S in states)
                need(parity == Counter({0: 28, 1: 28}), "balanced sector parity")
                need(zero_sets == [(1, 5, 7), (2, 4, 8)], "two exact resonant zero modes")
                print("c=3,M=8: 56 states; parity=28+28; nullity=2; "
                      "mode sets=" + repr(zero_sets))
    print(f"matrix/cyclotomic spectra: all {cases} cases 1<=c<=N<=8 PASS")

    for h in (3, 5, 7, 9):
        N = 3 * h - 1
        modulus = 2 * (N + 1)
        phi = cyclotomic(modulus)
        witnesses = []
        for j in range(1, h):
            S = (j, 2 * h - j, 2 * h + j)
            value = [0] * modulus
            for k in S:
                value[k] += 1
                value[modulus - k] += 1
            need(all(x == 0 for x in ring_reduce(value, phi)),
                 "three-cosine zero-mode family")
            witnesses.append(S)
        need(len(set(witnesses)) == h - 1, "distinct resonance modes")
        need(sum((-1) ** sum(S) for S in combinations(range(1, N + 1), 3)) == 0,
             "zero parity index in infinite-family controls")
    print("balanced resonance family h=3,5,7,9: exact cyclotomic controls PASS")

    states, A = token_matrix(4, 2)
    token_L = [[-x for x in row] for row in A]
    compound_L = [[-x for x in row] for row in A]
    for i, S in enumerate(states):
        token_L[i][i] = sum(A[i])
        compound_L[i][i] = sum(1 if x in (1, 4) else 2 for x in S)
        occupied_edges = sum(x + 1 in S for x in S)
        need(compound_L[i][i] - token_L[i][i] == 2 * occupied_edges,
             "exact exclusion diagonal potential")
    lap_char = matrix_charpoly(token_L)
    compound_char = matrix_charpoly(compound_L)
    need(lap_char != compound_char, "Laplacian free-sum hostile")
    print("P4 two-token Laplacian chi coefficients=" + repr(lap_char))
    print("additive one-particle Laplacian chi coefficients=" + repr(compound_char))
    need(sum(-1 for u, v in combinations(range(4), 2) if v == u + 1) == -3,
         "multiplicative exterior adjacency has trace -3, not native trace zero")

    F = graph_from_lengths(5, (3,))
    four_cycle = {(0, 1), (1, 3), (3, 4), (0, 4)}
    need(cycle_type(5, F ^ four_cycle) == (5,),
         "four-cycle layer adds C3-C5 skip edge")
    F = graph_from_lengths(4, (3,))
    need(cycle_type(4, F ^ set(combinations(range(3), 2))) == (),
         "sector is not ambient-invariant: C3 exits to empty")
    print("hostiles: active vertex cap, multiplicity, multiplicative compound, "
          "Laplacian potential, longer cycle, ambient boundary PASS")
    print(f"PASS {GATES} always-active gates")


if __name__ == "__main__":
    main()
