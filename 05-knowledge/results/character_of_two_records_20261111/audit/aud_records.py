"""AUDIT: section 6 (records; THM-5137 (R)), recomputed independently by brute force.

 R1 (b): mean of the Redei-Berge function U_T over all labelled tournaments on n <= 5 vertices, in the power-sum
         basis (U_T is computed from the Grinberg-Stanley odd-cycle formula: sets of disjoint directed odd cycles);
         q_n and the hook Schur functions: numerical identity at random integer points.
 R2:     mean weighted number of odd-cycle packings by max-ordered type; the identity sum 2^l/pi_u = 2, n <= 16.
 R3:     tournaments by the sizes of their pieces (cut the sequence of strong components before every record of the
         component maxima), graphs by the sizes of their components, both listed by increasing largest vertex;
         l_n = #{vertex n reaches every vertex}, n <= 6 here (n = 7, 8: aud_ln.c); formula (b).
 H check: U_T at p_k -> 1 equals the number of Hamiltonian paths (sanity of the cycle formula), every tournament.
"""
import itertools, random
from fractions import Fraction
from math import factorial, comb


def tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        A = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A


def graphs(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        adj = [0] * n
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
        yield adj


def odd_cycles(A):
    """All directed odd cycles as (vertex mask, length)."""
    n = len(A)
    res = []
    for L in range(3, n + 1, 2):
        for S in itertools.combinations(range(n), L):
            v0 = S[0]
            for p in itertools.permutations(S[1:]):
                seq = (v0,) + p
                if all(A[seq[i]][seq[(i + 1) % L]] for i in range(L)):
                    res.append((sum(1 << v for v in S), L))
    return res


def packings(cycles):
    """All sets of pairwise disjoint cycles (as lists of (mask, len))."""
    res = []

    def rec(i, used, cur):
        if i == len(cycles):
            res.append(list(cur))
            return
        rec(i + 1, used, cur)
        m, L = cycles[i]
        if not (m & used):
            cur.append(cycles[i])
            rec(i + 1, used | m, cur)
            cur.pop()

    rec(0, 0, [])
    return res


def ham(A):
    n = len(A)
    c = 0
    for p in itertools.permutations(range(n)):
        if all(A[p[i]][p[i + 1]] for i in range(n - 1)):
            c += 1
    return c


def z(lam):
    r = 1
    for k in set(lam):
        m = lam.count(k)
        r *= k ** m * factorial(m)
    return r


def compositions_odd(n):
    if n == 0:
        yield ()
        return
    for k in range(1, n + 1, 2):
        for rest in compositions_odd(n - k):
            yield (k,) + rest


def compositions(n):
    if n == 0:
        yield ()
        return
    for k in range(1, n + 1):
        for rest in compositions(n - k):
            yield (k,) + rest


def partial_sums(alpha):
    s, out = 0, []
    for a in alpha:
        s += a
        out.append(s)
    return out


def pi_u(alpha):
    r = 1
    for s in partial_sums(alpha):
        r *= s
    return r


print("== R1 (b), R2: means over all labelled tournaments")
for n in range(1, 6):
    by_partition = {}
    by_composition = {}
    count = 0
    hfail = 0
    for A in tournaments(n):
        count += 1
        cyc = odd_cycles(A)
        tot = 0
        for S in packings(cyc):
            w = 2 ** len(S)
            tot += w
            used = 0
            blocks = []
            for m, L in S:
                used |= m
                blocks.append((m.bit_length() - 1, L))  # (largest vertex, size)
            for v in range(n):
                if not used >> v & 1:
                    blocks.append((v, 1))
            blocks.sort()
            alpha = tuple(L for _, L in blocks)
            lam = tuple(sorted(alpha, reverse=True))
            by_partition[lam] = by_partition.get(lam, 0) + w
            by_composition[alpha] = by_composition.get(alpha, 0) + w
        if tot != ham(A):
            hfail += 1
    ok1 = True
    for lam, tot in by_partition.items():
        mean = Fraction(tot, count)
        pred = Fraction(factorial(n) * 2 ** len(lam), 2 ** n * z(lam))
        ok1 &= (mean == pred)
    lams_pred = set(tuple(sorted(a, reverse=True)) for a in compositions_odd(n))
    ok1 &= set(by_partition) == lams_pred
    ok2 = True
    rows = []
    for alpha in sorted(compositions_odd(n)):
        mean = Fraction(by_composition.get(alpha, 0), count)
        pred = Fraction(factorial(n) * 2 ** len(alpha), 2 ** n * pi_u(alpha))
        ok2 &= (mean == pred)
        rows.append("%s: %s" % (alpha, mean))
    print("n = %d: %d tournaments; U_T(1,0,...) = H(T) fails: %d; mean [p_lambda] = n! 2^(l-n)/z_lambda for all lambda: %s; "
          "mean by max-ordered type = n! 2^(l-n)/pi_u: %s" % (n, count, hfail, ok1, ok2))
    print("     " + "; ".join(rows))

print("== identity sum over compositions of n into odd parts of 2^l/pi_u = 2")
for n in range(1, 17):
    s = sum(Fraction(2 ** len(a), pi_u(a)) for a in compositions_odd(n))
    print("n = %2d: %s" % (n, s), end="; ")
print()

print("== q_n = sum over odd-part partitions of 2^l p_lambda/z_lambda = 2 * (sum of hook Schur functions), random integer points")
random.seed(5137)


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def schur_at(lam, xs):
    """Bialternant formula with exact rationals."""
    N = len(xs)
    lam = list(lam) + [0] * (N - len(lam))
    def det(M):
        M = [row[:] for row in M]
        d = Fraction(1)
        for c in range(N):
            p = None
            for r in range(c, N):
                if M[r][c] != 0:
                    p = r
                    break
            if p is None:
                return Fraction(0)
            if p != c:
                M[c], M[p] = M[p], M[c]
                d = -d
            d *= M[c][c]
            for r in range(c + 1, N):
                f = M[r][c] / M[c][c]
                for k in range(c, N):
                    M[r][k] -= f * M[c][k]
        return d
    num = det([[Fraction(x) ** (lam[j] + N - 1 - j) for j in range(N)] for x in xs])
    den = det([[Fraction(x) ** (N - 1 - j) for j in range(N)] for x in xs])
    return num / den


allok = True
for n in range(1, 8):
    for trial in range(3):
        N = n + 1
        xs = random.sample(range(1, 40), N)
        # q_n from the generating function prod (1 + x t)/(1 - x t)
        series = [Fraction(1)] + [Fraction(0)] * n
        for x in xs:
            # multiply by (1 + x t)
            new = series[:]
            for k in range(1, n + 1):
                new[k] += x * series[k - 1]
            # divide by (1 - x t): cumulative
            for k in range(1, n + 1):
                new[k] += x * new[k - 1]
            series = new
        qn = series[n]
        p = lambda k: sum(Fraction(x) ** k for x in xs)
        ps = Fraction(0)
        for lam in partitions(n):
            if all(k % 2 for k in lam):
                term = Fraction(2 ** len(lam), z(lam))
                for k in lam:
                    term *= p(k)
                ps += term
        hooks = sum(schur_at((n - k,) + (1,) * k, xs) for k in range(n))
        if not (qn == ps == 2 * hooks):
            allok = False
            print("FAIL", n, xs)
print("n = 1..7, three random points each (n + 1 variables): q_n = power-sum form = 2 * hooks: %s" % allok)

print("== R3: pieces of tournaments and components of graphs by max-ordered type")
A001187 = [1, 1, 1, 4, 38, 728, 26704, 1866256, 251548592]
strong = {}
lvals = {}
for n in range(1, 7):
    by_t, by_g = {}, {}
    ln = 0
    sn = 0
    for A in tournaments(n):
        out = [sum(1 << j for j in range(n) if A[i][j]) for i in range(n)]
        # reach sets
        def reach(v):
            r = 1 << v
            while True:
                nr = r
                for u in range(n):
                    if r >> u & 1:
                        nr |= out[u]
                if nr == r:
                    return r
                r = nr
        R = [reach(v) for v in range(n)]
        full = (1 << n) - 1
        if R[n - 1] == full:
            ln += 1
        if all(r == full for r in R):
            sn += 1
        # strong components in the order of the condensation: sort vertices by decreasing size of reach set
        comps = {}
        for v in range(n):
            comps.setdefault(R[v], []).append(v)
        order = sorted(comps.items(), key=lambda kv: -bin(kv[0]).count("1"))
        pieces = []
        best = -1
        for _, vs in order:
            mx = max(vs)
            if mx > best:
                pieces.append(len(vs))
                best = mx
            else:
                pieces[-1] += len(vs)
        alpha = tuple(pieces)
        by_t[alpha] = by_t.get(alpha, 0) + 1
    for adj in graphs(n):
        seen = 0
        comps = []
        for v in range(n):
            if not seen >> v & 1:
                r = 1 << v
                while True:
                    nr = r
                    for u in range(n):
                        if r >> u & 1:
                            nr |= adj[u]
                    if nr == r:
                        break
                    r = nr
                seen |= r
                comps.append((r.bit_length() - 1, bin(r).count("1")))
        comps.sort()
        alpha = tuple(s for _, s in comps)
        by_g[alpha] = by_g.get(alpha, 0) + 1
    lvals[n] = ln
    strong[n] = sn
    okc = True
    for alpha in compositions(n):
        pred = 1
        for s, a in zip(partial_sums(alpha), alpha):
            pred *= comb(s - 1, a - 1) * A001187[a]
        if by_t.get(alpha, 0) != pred or by_g.get(alpha, 0) != pred:
            okc = False
    print("n = %d: l_n = %d (A001187: %d); strong labelled tournaments S_n = %d; piece types of tournaments = component types of graphs = "
          "prod C(s_j - 1, a_j - 1) l_(a_j) for all %d compositions: %s" % (n, ln, A001187[n], sn, 2 ** (n - 1), okc))
    if n == 5:
        print("     type (1,3,1) at n = 5: tournaments %d, graphs %d" % (by_t.get((1, 3, 1), 0), by_g.get((1, 3, 1), 0)))
    if n == 4:
        print("     n = 4 tournaments by type: %s" % sorted(by_t.items()))

print("== formula (b): l_n = sum_j C(n-1, j-1) S_j 2^C(n-j, 2), with S_j from T = 1/(1 - S)")
# strong labelled tournaments from the EGF identity
NMAX = 12
T = [Fraction(2 ** comb(n, 2), factorial(n)) for n in range(NMAX + 1)]
# S = 1 - 1/T
inv = [Fraction(0)] * (NMAX + 1)
inv[0] = Fraction(1)
for n in range(1, NMAX + 1):
    inv[n] = -sum(T[k] * inv[n - k] for k in range(1, n + 1))
S = [0] + [int(-inv[n] * factorial(n)) for n in range(1, NMAX + 1)]
# connected graphs from log
logT = [Fraction(0)] * (NMAX + 1)
for n in range(1, NMAX + 1):
    logT[n] = T[n] - sum(Fraction(k, n) * logT[k] * T[n - k] for k in range(1, n))
conn = [0] + [int(logT[n] * factorial(n)) for n in range(1, NMAX + 1)]
okb = True
for n in range(1, NMAX + 1):
    ln = sum(comb(n - 1, j - 1) * S[j] * 2 ** comb(n - j, 2) for j in range(1, n + 1))
    if ln != conn[n]:
        okb = False
print("S_n (n = 1..8): %s" % S[1:9])
print("connected labelled graphs from log of the series (n = 1..8): %s" % conn[1:9])
print("formula (b) equals the connected-graph numbers for n <= %d: %s; brute-force S_n agree for n <= 6: %s; brute-force l_n agree: %s" % (
    NMAX, okb, all(strong[n] == S[n] for n in strong), all(lvals[n] == conn[n] for n in lvals)))
