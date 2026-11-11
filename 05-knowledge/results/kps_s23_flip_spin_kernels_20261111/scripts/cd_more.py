"""cd_more.py - second round on sign tournaments of twisted group algebras on Z_2^n (S23, after E3).

A. generalized Cayley-Dickson doubling with parameters g_1..g_n (g = -1 standard, g = +1 split):
   (a, b)(c, d) = (ac + g conj(d) b, d a + b conj(c)).  Lines: cyclic iff definite; split lines transitive with the
   definite unit in the middle.
B. zero-sum 4-sets: number of cyclic triangles inside (parity even <=> two-term zero divisors, exactly 2 pairings).
C. planes of PG(3,2) in the sedenions: bad 4-sets per plane, H per plane.
D. F_2 linear algebra on PG(3,2): largest sets of planes all of whose 4-sets can be odd at once; least number of
   even 4-sets over all tournaments (local search upper bound, plane bound lower bound).
E. counts for n = 6."""
import itertools, random, sys
from cd_signs import ham_paths


def gcd_mul(x, y, gs):
    n = len(x)
    if n == 1:
        return [x[0] * y[0]]
    h = n // 2
    g = gs[-1]
    rest = gs[:-1]
    a, b, c, d = x[:h], x[h:], y[:h], y[h:]
    conj = lambda v: [v[0]] + [-t for t in v[1:]]
    ac = gcd_mul(a, c, rest)
    db = gcd_mul(conj(d), b, rest)
    da = gcd_mul(d, a, rest)
    bc = gcd_mul(b, conj(c), rest)
    return [p + g * q for p, q in zip(ac, db)] + [p + q for p, q in zip(da, bc)]


def table(gs):
    n = len(gs)
    N = 1 << n
    g = [[0] * N for _ in range(N)]
    for i in range(N):
        ei = [0] * N
        ei[i] = 1
        for j in range(N):
            ej = [0] * N
            ej[j] = 1
            p = gcd_mul(ei, ej, gs)
            k = i ^ j
            assert all(p[t] == 0 for t in range(N) if t != k) and abs(p[k]) == 1
            g[i][j] = p[k]
    return g


def c3_of(arc, S):
    c = 0
    for (x, y, z) in itertools.combinations(S, 3):
        if (arc(x, y) and arc(y, z) and arc(z, x)) or (arc(y, x) and arc(z, y) and arc(x, z)):
            c += 1
    return c


def quads(n):
    N = 1 << n
    out = set()
    for a in range(1, N):
        for b in range(a + 1, N):
            for c in range(b + 1, N):
                d = a ^ b ^ c
                if d > c:
                    out.add((a, b, c, d))
    return sorted(out)


def report(gs, label, do_h=True):
    n = len(gs)
    N = 1 << n
    g = table(gs)
    V = range(1, N)
    sq = [g[i][i] for i in range(N)]
    anti = all(g[i][j] == -g[j][i] for i in V for j in V if i != j)
    arc = lambda i, j: g[i][j] == 1
    lines = sorted({tuple(sorted((a, b, a ^ b))) for a in V for b in V if a != b})
    cyc = 0
    rule_ok = True
    for L in lines:
        neg = [x for x in L if sq[x] == -1]
        is_c = c3_of(arc, L) == 1
        cyc += is_c
        if len(neg) == 3:
            rule_ok &= is_c
        else:
            # transitive with the definite unit in the middle
            ok = (not is_c) and len(neg) == 1
            if ok:
                m = neg[0]
                o = [x for x in L if x != m]
                ok = (arc(o[0], m) and arc(m, o[1])) or (arc(o[1], m) and arc(m, o[0]))
            rule_ok &= ok
    Q = quads(n)
    hist = {}
    for q in Q:
        k = c3_of(arc, q)
        hist[k] = hist.get(k, 0) + 1
    scores = sorted(set(sum(1 for j in V if j != i and arc(i, j)) for i in V))
    print("%s: squares -1 on %d units, +1 on %d; anticommuting %s; lines %d, cyclic %d; line rule holds: %s; scores %s"
          % (label, sum(1 for i in V if sq[i] == -1), sum(1 for i in V if sq[i] == 1), anti, len(lines), cyc, rule_ok, scores))
    print("   zero-sum 4-sets by number of cyclic triangles inside:", sorted(hist.items()), " even:", sum(v for k, v in hist.items() if k % 2 == 0))
    if do_h:
        adj = [[1 if (i != j and arc(i + 1, j + 1)) else 0 for j in range(N - 1)] for i in range(N - 1)]
        print("   H =", ham_paths(adj))
    return g


def planes(n):
    N = 1 << n
    out = set()
    for a in range(1, N):
        for b in range(a + 1, N):
            for c in range(b + 1, N):
                if c == a ^ b:
                    continue
                span = {0}
                for v in (a, b, c):
                    span |= {s ^ v for s in span}
                out.add(tuple(sorted(span - {0})))
    return sorted(out)


class F2System:
    def __init__(self, n):
        self.N = 1 << n
        self.var = {}
        for u in range(1, self.N):
            for v in range(u + 1, self.N):
                self.var[(u, v)] = len(self.var)

    def term(self, u, v):
        if u < v:
            return (1 << self.var[(u, v)], 1)       # [u -> v] = 1 + x_uv   (x = 1 iff v -> u)
        return (1 << self.var[(v, u)], 0)

    def quad_rows(self, q):
        a, b, c, d = q
        rows = []
        for (p, r) in (((a, b), (c, d)), ((a, c), (b, d))):   # two of the three pairings are independent
            mask, const = 0, 0
            for u in p:
                for v in r:
                    mk, cc = self.term(u, v)
                    mask ^= mk
                    const ^= cc
            rows.append((mask, const ^ 1))
        return rows

    @staticmethod
    def solve(rows):
        piv = {}
        ok = True
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
                    ok = False
        return ok, len(piv)


if __name__ == "__main__":
    report([-1, -1], "quaternions (-1,-1)")
    report([-1, 1], "split quaternions (-1,+1)")
    report([1, 1], "split quaternions (+1,+1)")
    report([-1, -1, -1], "octonions")
    report([-1, -1, 1], "split octonions (-1,-1,+1)")
    report([-1, 1, 1], "split octonions (-1,+1,+1)")
    g4 = report([-1, -1, -1, -1], "sedenions")
    report([-1, -1, -1, 1], "split sedenions (-1,-1,-1,+1)")
    # C. planes in the sedenions
    arc4 = lambda i, j: g4[i][j] == 1
    P = planes(4)
    Q4 = quads(4)
    res = []
    for pl in P:
        s = set(pl)
        qs = [q for q in Q4 if set(q) <= s]
        bad = sum(1 for q in qs if c3_of(arc4, q) % 2 == 0)
        adj = [[1 if (i != j and arc4(pl[i], pl[j])) else 0 for j in range(7)] for i in range(7)]
        res.append((bad, ham_paths(adj), c3_of(arc4, pl)))
    print("sedenion planes (even 4-sets, H, cyclic triangles):", sorted(res))
    # D. linear algebra
    S = F2System(4)
    rows_by_plane = []
    for pl in P:
        s = set(pl)
        rr = []
        for q in Q4:
            if set(q) <= s:
                rr += S.quad_rows(q)
        rows_by_plane.append(rr)
    best = 0
    count_by_size = {}
    for mask in range(1, 1 << 15):
        rr = []
        for i in range(15):
            if mask >> i & 1:
                rr += rows_by_plane[i]
        ok, rk = S.solve(rr)
        if ok:
            k = bin(mask).count("1")
            count_by_size[k] = count_by_size.get(k, 0) + 1
            best = max(best, k)
    print("PG(3,2): sets of planes whose 4-sets can all be odd at once, by size:", sorted(count_by_size.items()), " largest:", best)
    # all-odd planes of the sedenions
    odd_planes = [i for i, pl in enumerate(P) if all(c3_of(arc4, q) % 2 == 1 for q in Q4 if set(q) <= set(pl))]
    print("sedenions: planes with all 4-sets odd:", len(odd_planes))
    # E. counts for n = 5, 6
    for n in (5, 6):
        g = table([-1] * n)
        arc = lambda i, j: g[i][j] == 1
        Q = quads(n)
        ev = sum(1 for q in Q if c3_of(arc, q) % 2 == 0)
        pl = planes(n) if n == 5 else None
        print("n = %d: zero-sum 4-sets %d, even %d (zero-divisor pairings %d)" % (n, len(Q), ev, 2 * ev))
        if pl:
            cnt = {}
            for p_ in pl:
                s = set(p_)
                bad = 0
                for a, b, c in itertools.combinations(p_, 3):
                    d = a ^ b ^ c
                    if d in s and d > c and len({a, b, c, d}) == 4 and (a ^ b) != c:
                        bad += c3_of(arc, (a, b, c, d)) % 2 == 0
                cnt[bad] = cnt.get(bad, 0) + 1
            print("   planes of PG(4,2) by number of even 4-sets:", sorted(cnt.items()))
