"""a5_cd.py (AUDIT, independent) - items 5 and 8: Cayley-Dickson algebras by my own doubling.

Doubling (README 3.5): (a, b)(c, d) = (a c + g conj(d) b, d a + b conj(c)), conj(a, b) = (conj a, -b), g = -1 for the
standard algebras (g = +1 at a level gives a split form).  Vectors are integer lists; nothing is taken from the
author's code.

A. basis products e_i e_j = sign * e_(i xor j); anticommutation; squares.
B. Proposition 3.2 inside the algebra: for ALL ordered quadruples of distinct units (a, b, c, d) and both signs,
   (e_a + e_b)(e_c + s e_d) = 0 exactly when a^b = c^d and the number of arcs {a,b} -> {c,d} is even; counts.
C. Block recursion, all lines cyclic, Hadamard, doubly regular.
D. H(CD_4), automorphisms of CD_3, CD_4; planes of the sedenion tournament; even 4-sets for n = 4, 5, 6.
E. Generalised (split) forms of dimension 4, 8, 16, 32: line rule, units of square -1, cyclic lines, even 4-sets."""
import itertools, sys, os, random
from aud_common import ham_c, ham_py, HERE

out = []


def P(s):
    print(s)
    out.append(s)
    sys.stdout.flush()


def conj(x):
    if len(x) == 1:
        return x[:]
    h = len(x) // 2
    return conj(x[:h]) + [-t for t in x[h:]]


def add(x, y):
    return [a + b for a, b in zip(x, y)]


def neg(x):
    return [-a for a in x]


def mul(x, y, gam):
    """gam = list of doubling parameters, gam[-1] used at the top level."""
    if len(x) == 1:
        return [x[0] * y[0]]
    h = len(x) // 2
    a, b, c, d = x[:h], x[h:], y[:h], y[h:]
    g = gam[-1]
    sub = gam[:-1]
    first = add(mul(a, c, sub), [g * t for t in mul(conj(d), b, sub)])
    second = add(mul(d, a, sub), mul(b, conj(c), sub))
    return first + second


def basis(i, N):
    v = [0] * N
    v[i] = 1
    return v


def sign_table(gam):
    n = len(gam)
    N = 1 << n
    g = [[0] * N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            p = mul(basis(i, N), basis(j, N), gam)
            nz = [(k, v) for k, v in enumerate(p) if v]
            assert len(nz) == 1 and nz[0][0] == i ^ j and abs(nz[0][1]) == 1, (i, j, nz)
            g[i][j] = nz[0][1]
    return g


def tournament(g):
    N = len(g)
    return [set(j - 1 for j in range(1, N) if j != i and g[i][j] == 1) for i in range(1, N)]


def c3_count(arc, S):
    """number of cyclic triangles inside the vertex set S (arc(u, v) True iff u -> v)."""
    c = 0
    for u, v, w in itertools.combinations(S, 3):
        if (arc(u, v) and arc(v, w) and arc(w, u)) or (arc(v, u) and arc(w, v) and arc(u, w)):
            c += 1
    return c


def zero_sum_quads(N):
    qs = set()
    for a in range(1, N):
        for b in range(a + 1, N):
            for c in range(b + 1, N):
                d = a ^ b ^ c
                if d > c:
                    qs.add((a, b, c, d))
    return sorted(qs)


def lines(N):
    return sorted({tuple(sorted((a, b, a ^ b))) for a in range(1, N) for b in range(1, N) if a != b})


tables = {}
for n in range(1, 7):
    tables[n] = sign_table([-1] * n)
P("A. standard Cayley-Dickson algebras, dimensions 2..64: e_i e_j = +- e_(i xor j) for all basis pairs (checked by the doubling itself)")
for n in range(2, 7):
    g = tables[n]
    N = 1 << n
    anti = all(g[i][j] == -g[j][i] for i in range(1, N) for j in range(1, N) if i != j)
    sq = all(g[i][i] == -1 for i in range(1, N))
    unit = all(g[0][j] == 1 and g[j][0] == 1 for j in range(N))
    P("   n = %d: anticommutation of distinct imaginary units %s, squares -1: %s, e_0 is the unit: %s" % (n, anti, sq, unit))

# ---- B. Proposition 3.2 in the algebra
P("B. zero products (e_a + e_b)(e_c + s e_d), a, b, c, d distinct imaginary units:")
for n in (3, 4, 5):
    g = tables[n]
    N = 1 << n
    arc = lambda u, v: g[u][v] == 1
    gam = [-1] * n
    zero_pairings = set()      # (frozenset{a,b}, frozenset{c,d}) ordered: left factor, right factor
    zero_products = 0          # (unordered {a,b}, unordered {c,d}, s) with product zero
    mism = 0
    outside = 0
    tested = 0
    both_signs = 0
    full_mult_checked = 0
    for a in range(1, N):
        for b in range(a + 1, N):
            for c in range(1, N):
                if c in (a, b):
                    continue
                for d in range(c + 1, N):
                    if d in (a, b):
                        continue
                    tested += 1
                    zs = []
                    for s in (1, -1):
                        # product by bilinearity from the basis products computed by the doubling
                        vec = {}
                        for (u, v, coef) in ((a, c, 1), (a, d, s), (b, c, 1), (b, d, s)):
                            vec[u ^ v] = vec.get(u ^ v, 0) + coef * g[u][v]
                        z = all(val == 0 for val in vec.values())
                        if n <= 4 or (a * 31 + b * 17 + c * 7 + d) % 97 == 0:
                            # full multiplication of the two vectors by the doubling formula
                            x = add(basis(a, N), basis(b, N))
                            y = add(basis(c, N), [s * t for t in basis(d, N)])
                            z2 = all(t == 0 for t in mul(x, y, gam))
                            assert z2 == z
                            full_mult_checked += 1
                        if z:
                            zs.append(s)
                    arcs_ab_cd = sum(1 for u in (a, b) for v in (c, d) if arc(u, v))
                    pred = (a ^ b == c ^ d) and (arcs_ab_cd % 2 == 0)
                    if bool(zs) != pred:
                        mism += 1
                    if zs:
                        zero_products += len(zs)
                        if len(zs) == 2:
                            both_signs += 1
                        if a ^ b != c ^ d:
                            outside += 1
                        zero_pairings.add((frozenset((a, b)), frozenset((c, d))))
    unordered = {frozenset(p) for p in zero_pairings}
    sym = all((q, p) in zero_pairings for (p, q) in zero_pairings)
    quads = {tuple(sorted(set(p) | set(q))) for (p, q) in zero_pairings}
    pairs_involved = {p for (p, q) in zero_pairings} | {q for (p, q) in zero_pairings}
    # parity classification of all zero-sum 4-sets
    Q = zero_sum_quads(N)
    even = [q for q in Q if c3_count(arc, q) % 2 == 0]
    hist = {}
    for q in Q:
        k = c3_count(arc, q)
        hist[k] = hist.get(k, 0) + 1
    per_quad = {}
    for pr in unordered:
        q = tuple(sorted(set().union(*pr)))
        per_quad[q] = per_quad.get(q, 0) + 1
    P("   n = %d: (unordered {a,b}) x (unordered {c,d}) tested %d (x 2 signs; full vector multiplication on %d of them);" % (n, tested, full_mult_checked))
    P("      mismatches with the criterion [a^b = c^d and even number of arcs {a,b} -> {c,d}]: %d; zero products outside zero-sum 4-sets: %d; both signs zero: %d" % (mism, outside, both_signs))
    P("      zero products with left factor e_a + e_b: %d (ordered pairings), symmetric in the two factors: %s; unordered pairings {a,b}|{c,d}: %d" % (len(zero_pairings), sym, len(unordered)))
    P("      zero-sum 4-sets %d, by number of cyclic triangles %s; even 4-sets %d; 4-sets carrying zero pairings %d (all even: %s; exactly 2 pairings each: %s); pairs {a,b} involved %d"
      % (len(Q), sorted(hist.items()), len(even), len(quads), set(quads) == set(even), set(per_quad.values()) == {2} if per_quad else True, len(pairs_involved)))

# ---- C. block recursion, lines, Hadamard, doubly regular
P("C. sign tables G_n (order 2^n, with the unit):")
for n in range(1, 6):
    G = tables[n]
    G1 = tables[n + 1]
    N = 1 << n
    D = [1] + [-1] * (N - 1)
    ok = True
    for i in range(N):
        for j in range(N):
            ok &= G1[i][j] == G[i][j]
            ok &= G1[i][j + N] == G[j][i]
            ok &= G1[i + N][j] == G[i][j] * D[j]
            ok &= G1[i + N][j + N] == -G[j][i] * D[j]
    P("   block recursion G_(n+1) = [[G, G^T], [G D, -G^T D]] for n = %d -> %d: %s" % (n, n + 1, ok))
for n in range(2, 7):
    G = tables[n]
    N = 1 << n
    arc = lambda u, v: G[u][v] == 1
    L = lines(N)
    cyc = sum(1 for l in L if c3_count(arc, l) == 1)
    had = all(sum(G[i][k] * G[j][k] for k in range(N)) == (N if i == j else 0) for i in range(N) for j in range(N))
    hadT = all(sum(G[k][i] * G[k][j] for k in range(N)) == (N if i == j else 0) for i in range(N) for j in range(N))
    scores = sorted({sum(1 for j in range(1, N) if j != i and arc(i, j)) for i in range(1, N)})
    common = sorted({sum(1 for k in range(1, N) if k != i and k != j and arc(i, k) and arc(j, k)) for i in range(1, N) for j in range(i + 1, N)})
    P("   n = %d: lines %d, cyclic %d; G G^T = N I: %s, G^T G = N I: %s; scores %s; common out-neighbours of pairs %s (N/4 - 1 = %d)"
      % (n, len(L), cyc, had, hadT, scores, common, N // 4 - 1))

# ---- D. H, automorphisms, planes, even 4-sets
def automorphisms(T):
    n = len(T)
    adj = [[1 if j in T[i] else 0 for j in range(n)] for i in range(n)]
    count = 0
    img = [-1] * n
    used = [False] * n

    def rec(k):
        nonlocal count
        if k == n:
            count += 1
            return
        for v in range(n):
            if used[v]:
                continue
            if all(adj[i][k] == adj[img[i]][v] for i in range(k)):
                img[k] = v
                used[v] = True
                rec(k + 1)
                used[v] = False
        img[k] = -1
    rec(0)
    return count


P("D. the sign tournaments CD_n:")
for n in (2, 3, 4):
    T = tournament(tables[n])
    P("   CD_%d: H = %d, automorphisms %d" % (n, ham_c([T])[0], automorphisms(T)))
# is CD_3 the Paley tournament on 7 points?  (the unique doubly regular tournament of order 7; check H = 189 and
# double regularity above; explicit isomorphism test with QR_7)
T3 = tournament(tables[3])
paley = [set((i + q) % 7 for q in (1, 2, 4)) for i in range(7)]
iso = False
for p in itertools.permutations(range(7)):
    if all((p[j] in paley[p[i]]) == (j in T3[i]) for i in range(7) for j in range(7) if i != j):
        iso = True
        break
P("   CD_3 isomorphic to the Paley tournament QR_7: %s" % iso)


def planes(N):
    """hyperplanes through... all 3-dimensional subspaces minus 0 (as point sets of size 7)."""
    res = set()
    for a in range(1, N):
        for b in range(a + 1, N):
            for c in range(b + 1, N):
                if c == a ^ b:
                    continue
                res.add(tuple(sorted({a, b, c, a ^ b, a ^ c, b ^ c, a ^ b ^ c})))
    return sorted(res)


for n in (4, 5):
    G = tables[n]
    N = 1 << n
    arc = lambda u, v: G[u][v] == 1
    PL = planes(N)
    info = {}
    for pl in PL:
        qs = [q for q in zero_sum_quads(N) if set(q) <= set(pl)] if n == 4 else None
        sub = [set(pl.index(v) for v in pl if v != u and arc(u, v)) for u in pl]
        h = ham_py(sub)
        # even 4-sets of the plane = complements of its lines
        ev = 0
        for a, b in itertools.combinations(pl, 2):
            c = a ^ b
            if a < b < c:
                comp = [v for v in pl if v not in (a, b, c)]
                if c3_count(arc, comp) % 2 == 0:
                    ev += 1
        key = (ev, h)
        info[key] = info.get(key, 0) + 1
    P("   CD_%d: planes %d by (even 4-sets in the plane, H of the plane): %s" % (n, len(PL), sorted(info.items())))
for n in (4, 5, 6):
    G = tables[n]
    N = 1 << n
    arc = lambda u, v: G[u][v] == 1
    Q = zero_sum_quads(N)
    ev = sum(1 for q in Q if c3_count(arc, q) % 2 == 0)

    def gauss(m, k, q=2):
        num = 1
        den = 1
        for i in range(k):
            num *= q ** (m - i) - 1
            den *= q ** (i + 1) - 1
        return num // den
    P("   CD_%d: zero-sum 4-sets %d, even %d; 42 [n-1 choose 3]_2 = %d" % (n, len(Q), ev, 42 * gauss(n - 1, 3)))

# ---- E. generalised forms
P("E. generalised Cayley-Dickson algebras (all 2^n parameter choices):")
for n in (2, 3, 4, 5):
    N = 1 << n
    res = {}
    for gam in itertools.product((-1, 1), repeat=n):
        g = sign_table(list(gam))
        anti = all(g[i][j] == -g[j][i] for i in range(1, N) for j in range(1, N) if i != j)
        arc = lambda u, v: g[u][v] == 1
        L = lines(N)
        rule = True
        cyc = 0
        for (a, b, c) in L:
            negs = [u for u in (a, b, c) if g[u][u] == -1]
            k = c3_count(arc, (a, b, c))
            cyc += k
            if k == 1:
                rule &= (len(negs) == 3)
            else:
                # transitive: exactly one unit of square -1 and it is the middle vertex (score 1 inside the line)
                mid = [u for u in (a, b, c) if sum(1 for v in (a, b, c) if v != u and arc(u, v)) == 1]
                rule &= (len(negs) == 1 and mid == negs)
        negset = {0} | {u for u in range(1, N) if g[u][u] == -1}
        subgroup = all((u ^ v) in negset for u in negset for v in negset)
        ev = None
        if n in (3, 4):
            Q = zero_sum_quads(N)
            ev = sum(1 for q in Q if c3_count(arc, q) % 2 == 0)
        h = ham_py(tournament(g)) if n == 3 else None
        key = (anti, rule, subgroup, len(negset) - 1, cyc, ev, h)
        res.setdefault(key, []).append(gam)
    for key, lst in sorted(res.items(), key=lambda kv: str(kv[0])):
        P("   dim %d: anticommuting %s, line rule %s, {square -1} + 0 a subgroup %s, units of square -1: %d, cyclic lines %d, even 4-sets %s, H %s : %d parameter choices"
          % ((N,) + key + (len(lst),)))
open(os.path.join(HERE, "a5_cd.out"), "w").write("\n".join(out) + "\n")
