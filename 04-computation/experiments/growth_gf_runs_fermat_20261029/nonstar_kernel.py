"""HYP-9445 (3): the non-star kernel of the second permanental compound at order 6.
For each class whose compound kernel is not spanned by vertex stars, compute an exact rational kernel basis, reduce
modulo the star span, and report the minimal supports (as pair sets) of the quotient representatives."""
import itertools, sympy
from runs import tournaments, tour_from_bits
from compound_kernel import compound, star_vectors, rank_Q
n = 6
pairs = list(itertools.combinations(range(n), 2))
for bits in tournaments(n):
    A = tour_from_bits(n, bits)
    _, Pc = compound(A)
    K = sympy.Matrix(Pc).nullspace()
    if not K:
        continue
    stars = []
    for w in range(n):
        stars += star_vectors(A, w, pairs)
    span = rank_Q(stars) if stars else 0
    if span == len(K):
        continue
    sc = ''.join(map(str, sorted(sum(r) for r in A)))
    strong = all(len(x) == n for x in [set()] ) if False else None
    # supports of a reduced basis (row-reduce kernel basis and list minimal-support vectors)
    Kb = sympy.Matrix.hstack(*K).T.rref()[0]
    supports = []
    for i in range(Kb.rows):
        v = Kb.row(i)
        supp = [pairs[j] for j in range(len(pairs)) if v[j] != 0]
        supports.append((len(supp), supp))
    supports.sort()
    print(bits, sc, 'ker', len(K), 'stars', span, 'smallest supports:', [s for s in supports[:2]])
