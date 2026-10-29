"""The (n-2)-run unordered observer is the second permanental compound of the tournament matrix.

With k = n-2 paths a spanning cover has exactly two arcs, which match the two non-end vertices bijectively onto
the two non-start vertices (no cycle is possible with two arcs in a tournament).  Hence the plus sector
P_{n-2}[a, b] (unordered covers with start set a, end set b) equals per A[V-b, V-a]; in complement coordinates
alpha = V-a, beta = V-b it is Pc[alpha, beta] = per A[beta, alpha].

For every class of order n this script
  * cross-checks Pc against direct path-cover enumeration (n <= 7),
  * computes rank_Q(Pc) (exact elimination mod two 61-bit primes; agreement is required),
  * finds the 'star' kernel: vectors supported on the pairs through one vertex w, v = sum_y c(y)[{w,y}];
    Pc v = 0 iff C = c^T A[V-w, .] vanishes off N+(w) and, on N+(w), C(u1) = -C(u2) when |N+(w)| = 2,
    C = 0 when |N+(w)| >= 3 (no condition on N+(w) when |N+(w)| = 1, none at all at a sink),
  * reports whether the whole kernel is spanned by star vectors (the T5 phenomenon: tetrahedron boundaries).
usage: python3 compound_kernel.py n
"""
import sys
import itertools
from runs import tournaments, tour_from_bits

PR = (2305843009213693951, 4611686018427387847)   # 2^61-1 and a 62-bit prime


def rank_mod(M, p):
    M = [[x % p for x in row] for row in M]
    r = 0
    rows, cols = len(M), len(M[0]) if M else 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i][c]:
                piv = i
                break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [x * inv % p for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        r += 1
    return r


def rank_Q(M):
    r1, r2 = rank_mod(M, PR[0]), rank_mod(M, PR[1])
    if r1 != r2:
        raise RuntimeError('rank disagreement')
    return r1


def compound(A):
    n = len(A)
    pairs = list(itertools.combinations(range(n), 2))
    Pc = [[A[b[0]][a[0]] * A[b[1]][a[1]] + A[b[0]][a[1]] * A[b[1]][a[0]] for b in pairs] for a in pairs]
    return pairs, Pc


def direct_plus(A):
    """Unordered (n-2)-path covers: start set x end set, indexed by complements (non-starts x non-ends)."""
    n = len(A)
    k = n - 2
    pairs = list(itertools.combinations(range(n), 2))
    idx = {p: i for i, p in enumerate(pairs)}
    P = [[0] * len(pairs) for _ in pairs]
    for w in itertools.permutations(range(n)):
        breaks = [i for i in range(n - 1) if not A[w[i]][w[i + 1]]]
        if len(breaks) > k - 1:
            continue
        free = [i for i in range(n - 1) if A[w[i]][w[i + 1]]]
        for extra in itertools.combinations(free, k - 1 - len(breaks)):
            cuts = sorted(breaks + list(extra))
            starts = [w[0]] + [w[c + 1] for c in cuts]
            if starts != sorted(starts):
                continue
            ends = [w[c] for c in cuts] + [w[n - 1]]
            nonstart = tuple(sorted(set(range(n)) - set(starts)))
            nonend = tuple(sorted(set(range(n)) - set(ends)))
            P[idx[nonstart]][idx[nonend]] += 1
    return P


def star_kernel_dim(A, w):
    """Dimension of the space of c on V-w giving a kernel vector of the compound."""
    n = len(A)
    others = [y for y in range(n) if y != w]
    N = [u for u in range(n) if A[w][u]]
    if not N:
        return len(others)          # at a sink A[w, .] = 0, so every star vector is in the kernel (audit fix)
    eqs = []
    for u in range(n):
        if u not in N:
            eqs.append([A[y][u] for y in others])
    if len(N) == 2:
        eqs.append([A[y][N[0]] + A[y][N[1]] for y in others])
    elif len(N) >= 3:
        for u in N:
            eqs.append([A[y][u] for y in others])
    return len(others) - rank_Q(eqs)


def star_vectors(A, w, pairs):
    """Rational basis of the star kernel at w, as vectors on the pairs (exact, sympy nullspace)."""
    import sympy
    n = len(A)
    others = [y for y in range(n) if y != w]
    N = [u for u in range(n) if A[w][u]]
    eqs = []
    if N:                           # at a sink there is no constraint (audit fix)
        for u in range(n):
            if u not in N:
                eqs.append([A[y][u] for y in others])
        if len(N) == 2:
            eqs.append([A[y][N[0]] + A[y][N[1]] for y in others])
        elif len(N) >= 3:
            for u in N:
                eqs.append([A[y][u] for y in others])
    if eqs:
        basis = sympy.Matrix(eqs).nullspace()
    else:
        basis = [sympy.eye(len(others))[:, i] for i in range(len(others))]
    idx = {p: i for i, p in enumerate(pairs)}
    out = []
    for b in basis:
        den = sympy.ilcm(*[x.q for x in b]) if len(b) else 1
        v = [0] * len(pairs)
        for coef, y in zip(b, others):
            v[idx[tuple(sorted((w, y)))]] = int(coef * den)
        out.append(v)
    return out


def main():
    n = int(sys.argv[1])
    classes = tournaments(n)
    pairs = list(itertools.combinations(range(n), 2))
    stats = {}
    special = []
    for bits in classes:
        A = tour_from_bits(n, bits)
        _, Pc = compound(A)
        if n <= 7:
            D = direct_plus(A)
            if D != Pc:
                raise RuntimeError('compound != direct enumeration for ' + bits)
        r = rank_Q(Pc)
        kerdim = len(pairs) - r
        sinks = [v for v in range(n) if sum(A[v]) == 0]
        sources = [v for v in range(n) if sum(A[v]) == n - 1]
        sd = [star_kernel_dim(A, w) for w in range(n)]
        # span of all star vectors: build them as columns via nullspace mod p per vertex is avoided;
        # instead compare kerdim with the rank of the union of star kernels using the dual description:
        # the star kernels at different vertices intersect only in 0 except through vectors on a common
        # pair {w,y}; we record the per-vertex dimensions and the kernel dimension.
        stars = []
        for w in range(n):
            stars += star_vectors(A, w, pairs)
        for v in stars:
            if any(sum(Pc[i][j] * v[j] for j in range(len(v))) for i in range(len(Pc))):
                raise RuntimeError('star vector not in kernel')
        span = rank_Q(stars) if stars else 0
        strong = True
        reach = [{i} for i in range(n)]
        for _ in range(n):
            for i in range(n):
                for j in range(n):
                    if A[i][j]:
                        reach[i] |= reach[j]
        strong = all(len(x) == n for x in reach)
        key = ('strong' if strong else 'nonstrong', kerdim > 0, sum(sd) > 0, span == kerdim)
        stats[key] = stats.get(key, 0) + 1
        if strong and kerdim > 0:
            special.append((bits, ''.join(map(str, sorted(sum(r_) for r_ in A))), kerdim, span, sd))
    print('n=%d classes=%d  (compound order C(n,2)=%d; checked against direct enumeration: %s)'
          % (n, len(classes), len(pairs), n <= 7))
    for key in sorted(stats):
        print('  %-9s singular=%-5s starkernel=%-5s stars_span_kernel=%-5s : %d' % (key[0], key[1], key[2], key[3], stats[key]))
    print('strong classes with singular compound (bits, scores, kernel dim, span of stars, star-kernel dims by vertex):')
    for s in special[:60]:
        print('  ', s)
    if len(special) > 60:
        print('   ... %d more' % (len(special) - 60))


if __name__ == '__main__':
    main()
