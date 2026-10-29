"""(1) Trivial block of the one-arc observer: unordered covers by n-1 paths, start set V-y x end set V-x, equal A[x][y]
(so in complement coordinates the block is A^T), for every class n = 3..7; its sign block is (-1)^(|x-y|-1) A[x][y].
(2) For T5 the mod-2 kernel of the three-run plus sector is the cut space of K5 (the span of the vertex stars mod 2),
and the mod-2 connecting map B = (S-P)/2 maps it isomorphically onto the cokernel."""
import itertools
from census_runs import plus_sector
from runs import tournaments, tour_from_bits, circulant, run_matrix, unglued_k3
from ghost_lemma import rank2

for n in range(3, 8):
    for bits in tournaments(n):
        A = tour_from_bits(n, bits)
        sets, P = plus_sector(A, n - 1)
        idx = {s: i for i, s in enumerate(sets)}
        for x in range(n):
            for y in range(n):
                if x == y:
                    continue
                a = idx[tuple(v for v in range(n) if v != y)]
                b = idx[tuple(v for v in range(n) if v != x)]
                assert P[a][b] == A[x][y]
            assert P[idx[tuple(v for v in range(n) if v != x)]][idx[tuple(v for v in range(n) if v != x)]] == 0
    print('n=%d: trivial block of M_{n-1} = A^T in complement coordinates, all %d classes' % (n, len(tournaments(n))))
# sign block: signed sums over the end orders
for n in range(3, 7):
    for bits in tournaments(n):
        A = tour_from_bits(n, bits)
        tup = list(itertools.permutations(range(n), n - 1))
        idx = {t: i for i, t in enumerate(tup)}
        _, M = run_matrix(A, n - 1, tup)
        def perm_sign(p):
            s, p = 1, list(p)
            for i in range(len(p)):
                while p[i] != i:
                    j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
            return s
        for y in range(n):
            for x in range(n):
                if x == y:
                    continue
                start = tuple(v for v in range(n) if v != y)
                endset = sorted(v for v in range(n) if v != x)
                tot = 0
                for order in itertools.permutations(endset):
                    perm = [endset.index(v) for v in order]
                    tot += perm_sign(perm) * M[idx[start]][idx[order]]
                assert tot == (-1) ** (abs(x - y) - 1) * A[x][y], (n, bits, x, y, tot)
    print('n=%d: sign block = (-1)^(|x-y|-1) A[x][y]' % n)
# (2) T5
A = circulant(5, {1, 2})
tuples, M = run_matrix(A, 3)
N, P, S, L0, L1 = unglued_k3(M, len(tuples) // 6)
sets = list(itertools.combinations(range(5), 3))
# plus-sector index = end SETS (columns) -- star at w: triangles Q-y for Q = V - w, as indicator mod 2
stars = []
for w in range(5):
    v = [0] * len(sets)
    for y in range(5):
        if y == w: continue
        tri = tuple(sorted(set(range(5)) - {w, y}))
        v[sets.index(tri)] = 1
    stars.append(v)
r = len(P)
print('T5: rank_F2 P =', rank2(P), '; rank_F2 of the 5 star indicators =', rank2(stars),
      '; stars in ker P mod 2:', all(sum(P[i][j] * v[j] for j in range(r)) % 2 == 0 for v in stars for i in range(r)))
B = [[(S[i][j] - P[i][j]) // 2 for j in range(r)] for i in range(r)]
block = [P[i] + [0] * r for i in range(r)] + [B[i] + P[i] for i in range(r)]
print('T5: connecting rank (spur_2) =', rank2(block) - 2 * rank2(P))
