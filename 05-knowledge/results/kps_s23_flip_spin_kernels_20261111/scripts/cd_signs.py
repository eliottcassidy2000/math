"""cd_signs.py - E3 of the S23 pre-registration: sign tournaments of the Cayley-Dickson algebras.

Standard doubling (a, b)(c, d) = (ac - conj(d) b, d a + b conj(c)), conj(a, b) = (conj(a), -b).
Basis e_0 .. e_(2^n - 1), e_(i + 2^(n-1)) = (0, e_i); e_i e_j = g(i, j) e_(i xor j).
Sign tournament CD_n on the nonzero indices: i -> j iff g(i, j) = +1 (the units anticommute, so g(i,j) = -g(j,i)).

Checks: lines {a, b, a^b} cyclic; regular; doubly regular (sign table Hadamard); H; zero-sum quadrangles
{a,b}|{c,d} (a^b = c^d) with an even number of arcs each way (two-term zero divisors (e_a + e_b)(e_c +- e_d) = 0,
verified by multiplying in the algebra); the affine F_2 system 'no even quadrangle' on PG(n-1, 2)."""
import itertools, sys


def cd_mul(x, y):
    """x, y lists of length 2^n (integer coordinates)."""
    n = len(x)
    if n == 1:
        return [x[0] * y[0]]
    h = n // 2
    a, b, c, d = x[:h], x[h:], y[:h], y[h:]
    conj = lambda v: [v[0]] + [-t for t in v[1:]]
    ac = cd_mul(a, c)
    db = cd_mul(conj(d), b)
    da = cd_mul(d, a)
    bc = cd_mul(b, conj(c))
    return [p - q for p, q in zip(ac, db)] + [p + q for p, q in zip(da, bc)]


def sign_table(n):
    N = 1 << n
    g = [[0] * N for _ in range(N)]
    for i in range(N):
        ei = [0] * N
        ei[i] = 1
        for j in range(N):
            ej = [0] * N
            ej[j] = 1
            p = cd_mul(ei, ej)
            k = i ^ j
            assert all(p[t] == 0 for t in range(N) if t != k) and abs(p[k]) == 1
            g[i][j] = p[k]
    return g


def ham_paths(adj):
    n = len(adj)
    out = [sum(1 << j for j in range(n) if adj[i][j]) for i in range(n)]
    dp = [dict() for _ in range(1 << n)]
    cnt = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        cnt[1 << v][v] = 1
    full = (1 << n) - 1
    for s in range(1, full + 1):
        row = cnt[s]
        for v in range(n):
            c = row[v]
            if c:
                nx = out[v] & ~s
                while nx:
                    u = (nx & -nx).bit_length() - 1
                    nx &= nx - 1
                    cnt[s | 1 << u][u] += c
    return sum(cnt[full])


def analyse(n, g, label, do_h=True):
    N = 1 << n
    V = list(range(1, N))
    assert all(g[i][i] == -1 for i in V) and all(g[i][j] == -g[j][i] for i in V for j in V if i != j)
    arc = lambda i, j: g[i][j] == 1
    lines = sorted({tuple(sorted((a, b, a ^ b))) for a in V for b in V if a != b})
    cyc = sum(1 for (a, b, c) in lines if (arc(a, b) and arc(b, c) and arc(c, a)) or (arc(b, a) and arc(c, b) and arc(a, c)))
    scores = sorted(sum(1 for j in V if j != i and arc(i, j)) for i in V)
    common = {}
    for i in V:
        for j in V:
            if i < j:
                c = sum(1 for k in V if k != i and k != j and arc(i, k) and arc(j, k))
                common[c] = common.get(c, 0) + 1
    rows_orth = all(sum(g[i][k] * g[j][k] for k in range(N)) == 0 for i in range(N) for j in range(i))
    # zero-sum quadrangle pairings
    even = 0
    total = 0
    verified = 0
    quads_even = set()
    pairs_in = set()
    for fsum in V:
        pairs = sorted({tuple(sorted((a, a ^ fsum))) for a in V if a != fsum})
        for (p, q) in itertools.combinations(pairs, 2):
            total += 1
            a, b = p
            c, d = q
            fw = sum(1 for u in (a, b) for v in (c, d) if arc(u, v))
            if fw % 2 == 0:
                even += 1
                quads_even.add(tuple(sorted((a, b, c, d))))
                pairs_in.add(p)
                pairs_in.add(q)
                # verify in the algebra: (e_a + e_b)(e_c + s e_d) = 0 for one sign s
                ok = False
                for s in (1, -1):
                    x = [0] * N
                    y = [0] * N
                    x[a] = 1
                    x[b] = 1
                    y[c] = 1
                    y[d] = s
                    if all(t == 0 for t in cd_mul(x, y)) if label.startswith("CD") else False:
                        ok = True
                verified += ok
    print("%s: %d vertices, %d lines, cyclic lines %d, scores %s" % (label, N - 1, len(lines), cyc, sorted(set(scores))))
    print("   common out-neighbour counts over pairs: %s ; full sign table Hadamard: %s" % (sorted(common.items()), rows_orth))
    print("   zero-sum quadrangle pairings: %d, even ones: %d (algebra check passed on %d), distinct 4-sets %d, pairs {a,b} involved %d"
          % (total, even, verified, len(quads_even), len(pairs_in)))
    if do_h:
        adj = [[1 if (i != j and arc(V[i], V[j])) else 0 for j in range(N - 1)] for i in range(N - 1)]
        print("   H =", ham_paths(adj))
    return even


def affine_system(n, extra_cyclic_lines=False):
    """Unknown x_{uv} (u < v nonzero): arc u -> v iff x = 0.  Equation per quadrangle pairing {a,b}|{c,d}:
    number of arcs from {a,b} to {c,d} is odd.  Returns (consistent, dimension of the solution space, one solution)."""
    N = 1 << n
    V = list(range(1, N))
    var = {}
    for u in V:
        for v in V:
            if u < v:
                var[(u, v)] = len(var)
    nv = len(var)
    rows = []

    def term(u, v):
        # indicator of arc u -> v as affine form: (mask, const):  u<v: 1 + x_uv ; u>v: x_vu
        if u < v:
            return (1 << var[(u, v)], 1)
        return (1 << var[(v, u)], 0)
    for fsum in V:
        pairs = sorted({tuple(sorted((a, a ^ fsum))) for a in V if a != fsum})
        for (p, q) in itertools.combinations(pairs, 2):
            mask, const = 0, 0
            for u in p:
                for v in q:
                    mk, c = term(u, v)
                    mask ^= mk
                    const ^= c
            rows.append((mask, const ^ 1))        # sum of indicators = 1
    if extra_cyclic_lines:
        lines = sorted({tuple(sorted((a, b, a ^ b))) for a in V for b in V if a != b})
        for (a, b, c) in lines:                   # cyclic iff [a->b] = [b->c] = [c->a]
            for (p, q) in (((a, b), (b, c)), ((b, c), (c, a))):
                m1, c1 = term(*p)
                m2, c2 = term(*q)
                rows.append((m1 ^ m2, c1 ^ c2))
    # Gaussian elimination over F_2
    piv = {}
    consistent = True
    for mask, const in rows:
        while mask:
            hb = mask.bit_length() - 1
            if hb in piv:
                pm, pc = piv[hb]
                mask ^= pm
                const ^= pc
            else:
                piv[hb] = (mask, const)
                break
        else:
            if const:
                consistent = False
    return consistent, nv - len(piv), nv, len(rows), len(piv)


if __name__ == "__main__":
    for n in (2, 3, 4):
        g = sign_table(n)
        analyse(n, g, "CD_%d" % n)
    g5 = sign_table(5)
    analyse(5, g5, "CD_5", do_h=False)
    for n in (3, 4):
        for ex in (False, True):
            print("affine system on PG(%d,2), %s: consistent %s, solution dimension %d (unknowns %d, equations %d, rank %d)"
                  % ((n - 1, "with all lines cyclic" if ex else "quadrangles only") + affine_system(n, ex)))
