"""Star law for the second permanental compound (THM-4846 (B), HYP-9445 (3)), recomputed independently.

Pc[alpha, beta] = per A[beta, alpha] on 2-subsets (complement coordinates of the (n-2)-run plus sector).
Star kernel at w := kernel of Pc restricted to the columns {w,y}, y != w (computed directly, exact rational rank via
FLINT; no case analysis).  The star law holds iff the sum of the star kernels is the whole kernel.
Also reported: the producer's criterion (compound_kernel.py star_kernel_dim), which imposes C(u) = 0 for every u
even at a sink, to measure its effect.
For n <= 7, Pc is cross-checked against the direct (n-2)-run plus sector of runlib (subset DP).
"""
import sys
import itertools
import flint
from runlib import tournaments, adj_from_bits, is_strong, plus_sector, hp_table


def compound(A):
    n = len(A)
    pairs = list(itertools.combinations(range(n), 2))
    Pc = [[A[b[0]][a[0]] * A[b[1]][a[1]] + A[b[0]][a[1]] * A[b[1]][a[0]] for b in pairs] for a in pairs]
    return pairs, Pc


def nullspace_basis(cols_matrix):
    """Rational nullspace of an integer matrix (list of rows) as integer column vectors."""
    X = flint.fmpz_mat(cols_matrix)
    Nm, nullity = X.nullspace()
    vecs = []
    for j in range(nullity):
        vecs.append([int(Nm[i, j]) for i in range(Nm.nrows())])
    for v in vecs:
        assert all(sum(r[c] * v[c] for c in range(len(v))) == 0 for r in cols_matrix)
    return vecs


def producer_star_dim(A, w):
    n = len(A)
    others = [y for y in range(n) if y != w]
    N = [u for u in range(n) if A[w][u]]
    eqs = []
    for u in range(n):
        if u not in N:
            eqs.append([A[y][u] for y in others])
    if len(N) == 2:
        eqs.append([A[y][N[0]] + A[y][N[1]] for y in others])
    elif len(N) >= 3:
        for u in N:
            eqs.append([A[y][u] for y in others])
    if not eqs:
        return len(others)
    return len(others) - flint.fmpz_mat(eqs).rank()


def main():
    for n in [int(x) for x in sys.argv[1:]]:
        pairs = list(itertools.combinations(range(n), 2))
        stats = {}
        prod_fail = 0
        examples = []
        for bits in tournaments(n):
            A = adj_from_bits(n, bits)
            _, Pc = compound(A)
            if n <= 7:
                sets, P = plus_sector(A, n - 2, hp_table(A))
                comp = {s: tuple(v for v in range(n) if v not in s) for s in sets}
                pidx = {p: i for i, p in enumerate(pairs)}
                for i, a in enumerate(sets):
                    for j, b in enumerate(sets):
                        assert P[i][j] == Pc[pidx[comp[a]]][pidx[comp[b]]]
            r = flint.fmpz_mat(Pc).rank()
            kdim = len(pairs) - r
            starvecs = []
            stardims = []
            for w in range(n):
                cols = [i for i, p in enumerate(pairs) if w in p]
                sub = [[Pc[i][j] for j in cols] for i in range(len(pairs))]
                basis = nullspace_basis(sub)
                stardims.append(len(basis))
                for b in basis:
                    v = [0] * len(pairs)
                    for c, x in zip(cols, b):
                        v[c] = x
                    starvecs.append(v)
            span = flint.fmpz_mat(starvecs).rank() if starvecs else 0
            law = span == kdim
            strong = is_strong(A)
            sinks = [v for v in range(n) if sum(A[v]) == 0]
            # producer's span (dimension only via its per-vertex criterion is not enough; replicate its vectors)
            pd = [producer_star_dim(A, w) for w in range(n)]
            key = ('strong' if strong else 'nonstr', 'singular' if kdim else 'nonsing', 'law' if law else 'FAIL')
            stats[key] = stats.get(key, 0) + 1
            if pd != stardims:
                prod_fail += 1
                if len(examples) < 5:
                    examples.append((bits, 'sinks', sinks, 'true star dims', stardims, 'producer dims', pd,
                                     'kernel', kdim, 'span', span))
        fails = sum(v for k, v in stats.items() if k[2] == 'FAIL')
        print('n=%d: star-law failures %d; breakdown %s' % (n, fails, sorted(stats.items())))
        print('     classes where the producer star-kernel criterion differs from the direct one: %d; e.g. %s'
              % (prod_fail, examples[:3]))


if __name__ == '__main__':
    main()
