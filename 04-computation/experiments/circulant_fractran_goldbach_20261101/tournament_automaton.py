#!/usr/bin/env python3
"""mac-mini-2026-11-01-lenses: the GF(2) automaton of a tournament - a translation from tournaments to (linear)
cellular automata.

For a tournament T with adjacency matrix A (A[u][v] = 1 iff u -> v), let Phi_T : GF(2)^V -> GF(2)^V, x -> A x: each
vertex takes the parity of the bits of the vertices it beats. For a circulant Circ(Z_q, S) this is the additive
cellular automaton on the ring Z_q with neighbourhood S: (Phi x)_i = sum_{s in S} x_{i+s} mod 2.

Facts used / tested:
  * an isomorphism T -> T' conjugates Phi_T and Phi_T' (permutation matrices), so the dynamics of Phi_T are
    isomorphism invariants;
  * T^op has matrix A^T, and A, A^T are similar over any field, so Phi_{T^op} is conjugate to Phi_T: the automaton
    is converse-blind and factors through the MERGED metagraph (classes up to converse);
  * circulant: rank_GF(2) A = q - deg gcd(a(x), x^q - 1), a(x) = sum_{s in S} x^s (q odd);
  * Paley P_q (q prime = 3 mod 4, S = quadratic residues): over GF(2^m), a(alpha)^2 = a(alpha^2); if q = 7 mod 8
    then 2 is a residue, a(alpha) in GF(2) and a(alpha) + a(alpha^n) = 1 (n a non-residue), so a vanishes on exactly
    one residue class of roots and rank = (q+1)/2, the image being a binary quadratic-residue code (CITED:
    MacWilliams-Sloane ch. 16); if q = 3 mod 8, a(alpha) is a root of y^2 + y + 1, never 0, and Phi is invertible.

PRE-REGISTERED (before the first run):
  A1  For every tournament class on n <= 6 vertices (2, 4, 12, 56 classes) the dynamical invariants of Phi_T
      (rank, nilpotent-part size, cycle multiset of the periodic part on all 2^n states) coincide for T and T^op,
      and are unchanged under 20 random relabellings of each class representative.
  A2  For every circulant tournament of prime order q <= 31 (all connection sets with 1 in S), Gaussian rank over
      GF(2) equals q - deg gcd(a(x), x^q - 1).
  A3  Paley tournaments, q in {3, 7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83}: rank = q when q = 3 mod 8 and
      (q+1)/2 when q = 7 mod 8; the image for q = 7 has minimum distance 3 (the Hamming [7,4] code) and for q = 23
      minimum distance 7 (the binary Golay [23,12] code).
  A4  Every circulant tournament of order q = 1 mod 4 has the all-ones vector in its kernel (rank <= q - 1).
  A5  (exploratory, no prediction) the invariants of the ten merged n = 5 nodes, and whether they single out the
      three pure-blue nodes.
usage: python3 tournament_automaton.py
"""
import itertools
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = []


def check(name, ok, detail=''):
    print(f'  [{"PASS" if ok else "FAIL"}] {name}' + (f'  {detail}' if detail else ''))
    if not ok:
        FAIL.append(name)


def rank_gf2(rows, n):
    """rows: list of int bitmasks"""
    rows = list(rows)
    r = 0
    for c in range(n):
        piv = next((i for i in range(r, len(rows)) if rows[i] >> c & 1), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i] >> c & 1:
                rows[i] ^= rows[r]
        r += 1
    return r


def apply(rows, x):
    """y = A x over GF(2); rows[u] = bitmask of out-neighbours of u"""
    y = 0
    for u, row in enumerate(rows):
        if bin(row & x).count('1') & 1:
            y |= 1 << u
    return y


def dynamics(rows, n):
    """functional graph of x -> Ax on 2^n states: rank, nilpotent depth, multiset of cycle lengths"""
    f = [apply(rows, x) for x in range(1 << n)]
    # periodic points: iterate image until stable
    img = set(range(1 << n))
    depth = 0
    while True:
        nxt = {f[x] for x in img}
        if nxt == img:
            break
        img = nxt
        depth += 1
    seen, cycles = set(), []
    for x in sorted(img):
        if x in seen:
            continue
        L, y = 0, x
        while True:
            seen.add(y)
            y = f[y]
            L += 1
            if y == x:
                break
        cycles.append(L)
    return rank_gf2(rows, n), depth, tuple(sorted(cycles)), len(img)


def adj_rows(A):
    return [sum(1 << v for v in range(len(A)) if A[u][v]) for u in range(len(A))]


def classes(n):
    """all tournament classes on n vertices, as adjacency matrices. Every tournament has a Hamiltonian path
    (Redei), so it suffices to orient the non-path pairs of the fixed path n-1 -> ... -> 0 (the tiling model);
    classes are identified by the least adjacency string over all relabellings."""
    perms = list(itertools.permutations(range(n)))
    tiles = [(i, j) for i in range(n) for j in range(i + 2, n)]
    seen = {}
    for bits in range(1 << len(tiles)):
        A = [[0] * n for _ in range(n)]
        for i in range(n - 1):
            A[i + 1][i] = 1
        for k, (i, j) in enumerate(tiles):
            if bits >> k & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        key = min(tuple(A[p[i]][p[j]] for i in range(n) for j in range(n)) for p in perms)
        if key not in seen:
            seen[key] = A
    return list(seen.values())


def polymod_gcd(a, b):
    """gcd of GF(2) polynomials as int bitmasks"""
    while b:
        while a and a.bit_length() >= b.bit_length():
            a ^= b << (a.bit_length() - b.bit_length())
        a, b = b, a
    return a


def circ_rows(q, S):
    return [sum(1 << ((u + s) % q) for s in S) for u in range(q)]


def qr(q):
    return {(x * x) % q for x in range(1, q)}


def min_distance_of_image(rows, n):
    # image of A = span of the columns; build generator from a basis of the column space
    cols = [sum((rows[u] >> v & 1) << u for u in range(n)) for v in range(n)]
    basis = []
    for c in cols:
        x = c
        for b in basis:
            x = min(x, x ^ b)
        if x:
            basis.append(x)
    best = n
    for k in range(1, 1 << len(basis)):
        w = 0
        for i, b in enumerate(basis):
            if k >> i & 1:
                w ^= b
        best = min(best, bin(w).count('1'))
    return len(basis), best


def main():
    random.seed(20261101)
    print('A1  invariants of Phi_T on all classes n <= 6: converse-blind and relabelling-invariant')
    ok1 = True
    store = {}
    for n in (3, 4, 5, 6):
        cl = classes(n)
        for A in cl:
            inv = dynamics(adj_rows(A), n)
            At = [[A[j][i] for j in range(n)] for i in range(n)]
            if dynamics(adj_rows(At), n) != inv:
                ok1 = False
            for _ in range(20):
                p = list(range(n))
                random.shuffle(p)
                B = [[A[p[i]][p[j]] for j in range(n)] for i in range(n)]
                if dynamics(adj_rows(B), n) != inv:
                    ok1 = False
        invs = {}
        for A in cl:
            invs.setdefault(dynamics(adj_rows(A), n), 0)
            invs[dynamics(adj_rows(A), n)] += 1
        store[n] = invs
        print(f'    n={n}: {len(cl)} classes, {len(invs)} distinct automaton types; types (rank, depth, cycles, '
              f'#periodic): {sorted(invs.items(), key=lambda kv: -kv[1])[:6]}')
    check('A1: converse-blind and relabelling-invariant on all classes n <= 6', ok1)
    print('A2  circulants of prime order q <= 31: Gaussian rank = q - deg gcd(a, x^q - 1)')
    ok2, cnt = True, 0
    for q in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        pairs = [(x, q - x) for x in range(2, (q + 1) // 2)]
        for choice in itertools.product((0, 1), repeat=len(pairs)):
            S = [1] + [p[c] for p, c in zip(pairs, choice)]
            a = sum(1 << s for s in S)
            g = polymod_gcd((1 << q) | 1, a)
            if rank_gf2(circ_rows(q, S), q) != q - (g.bit_length() - 1):
                ok2 = False
            cnt += 1
    check('A2: rank formula on all circulant tournaments of prime order <= 31', ok2, f'{cnt} connection sets')
    print('A3  Paley tournaments: rank by q mod 8; Hamming and Golay images')
    ok3, rows3 = True, []
    for q in (3, 7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83):
        r = rank_gf2(circ_rows(q, sorted(qr(q))), q)
        want = q if q % 8 == 3 else (q + 1) // 2
        rows3.append((q, q % 8, r, want))
        ok3 &= r == want
    dim7, d7 = min_distance_of_image(circ_rows(7, sorted(qr(7))), 7)
    dim23, d23 = min_distance_of_image(circ_rows(23, sorted(qr(23))), 23)
    print(f'    (q, q mod 8, rank, predicted): {rows3}')
    print(f'    image of P_7: dimension {dim7}, minimum distance {d7}; image of P_23: dimension {dim23}, '
          f'minimum distance {d23}')
    check('A3: Paley rank = q (q = 3 mod 8), (q+1)/2 (q = 7 mod 8); P_7 -> [7,4,3], P_23 -> [23,12,7]',
          ok3 and (dim7, d7) == (4, 3) and (dim23, d23) == (12, 7))
    print('A4  q = 1 mod 4: all-ones vector in the kernel')
    ok4 = True
    for q in (5, 9, 13, 17, 21, 25, 29):
        pairs = [(x, q - x) for x in range(2, (q + 1) // 2)]
        for choice in itertools.product((0, 1), repeat=len(pairs)):
            S = [1] + [p[c] for p, c in zip(pairs, choice)]
            if apply(circ_rows(q, S), (1 << q) - 1) != 0:
                ok4 = False
    check('A4: all-ones in the kernel for every circulant of order 5, 9, ..., 29', ok4)
    print('A5  exploratory: the ten merged n = 5 nodes (owner explorer labels)')
    M = json.load(open(os.path.join(HERE, '..', 'groups_of_five_primegame_20261031', 'merged_metagraph_n5.json')))
    for node in M['merged_nodes']:
        c = M['classes'][node['members'][0]]
        n = 5
        A = [[int(c['sig'][i * n + j]) for j in range(n)] for i in range(n)]
        inv = dynamics(adj_rows(A), n)
        print(f'    node {node["members"]}: type {c["type"]}, scores {c["scores"]}, H {c["H"]}: rank {inv[0]}, '
              f'nilpotent depth {inv[1]}, cycles {inv[2]}, periodic states {inv[3]}')
    print()
    print('ALL PRE-REGISTERED CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
    return 0 if not FAIL else 1


if __name__ == '__main__':
    sys.exit(main())
