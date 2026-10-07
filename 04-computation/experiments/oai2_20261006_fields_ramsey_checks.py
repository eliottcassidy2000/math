#!/usr/bin/env python3
"""Independent checks of the session readers' claims on openai/math #158/#172 (plane colouring, Euclidean Ramsey) and
#164/#189 (Hindman, cycle-clique) against our Hadwiger-Nelson and Erdos-592 fronts (mac-mini, 2026-10-06).

  1. kappa(q) := chi(Cay(F_(q^2), mu_(q+1))), the residue graphs of Lemma R (exact, CaDiCaL via pysat):
     kappa = 4, 3, 4, 4, 4, 4, 3, 5, 4 for q = 2, 3, 4, 5, 7, 8, 9, 11, 16.
     Cay(F_49, mu_8) is the 7x7 knight torus G_7 of THM-4552 (alpha*mu_8 = the circle x^2 + y^2 = 5 = the knight set).
  2. Lemma R in the field Q(sqrt-3, sqrt-19) at the prime over 3 (and Q(sqrt-3, sqrt-7)): every unit vector
     u = z / conj(z) (Hilbert 90; all unit vectors arise so) has 3-integral coordinates in the basis
     {1, sqrt-3, sqrt-D, sqrt-3 sqrt-D}, and its residue a + c*s (s^2 = -D = -1 in F_9) lies in mu_4; hence the residue
     colouring of Cay(F_9, mu_4) = K_3 x K_3 3-colours the whole plane. Random tests (FINITE-EXACT sanity check of a
     PROVED lemma).
  3. Erdos 592 on [t]^2: the FULLY gap-determined game (edge iff y - x in E) dies at t = 4 (Q_inv(2,3) SAT,
     Q_inv(2,4) UNSAT), while the ROW-invariant game of THM-453 F/G (R_a = R, B_(a,a') = B_(a'-a)) dies at t = 5
     (SAT at 4, UNSAT at 5). The two notions differ; THM-470's "Finv = the translation-invariant game of THM-453 F" is
     a mislabel (Finv is the fully gap-determined algebra).
  4. Schur triples with 3-smooth summands: 2-colourings of [N] avoiding monochromatic {a, b, a + b} (a, b 3-smooth):
     least forcing N = 5 (a = b allowed) and 13 (a != b); and the parity colouring Omega mod 2 of the 3-smooth numbers
     has no monochromatic Schur triple inside S below 10^12 (Gersonides: the primitive solutions are 1+1, 1+2, 1+3, 1+8).
Run: python3 oai2_20261006_fields_ramsey_checks.py   (about 1-3 min)
"""
import itertools, math, random
from fractions import Fraction as Fr
from pysat.solvers import Solver
from pysat.card import CardEnc

OK = True


def check(cond, msg):
    global OK
    print(("  ok   " if cond else "  FAIL ") + msg, flush=True)
    OK &= bool(cond)


# ---------- finite fields GF(p^m) as integers 0..p^m-1 (base-p digit vectors) ----------
def gf(p, m):
    """returns (elements count, add, mul) tables for GF(p^m) using a monic irreducible of degree m."""
    def polymulmod(a, b, mod):
        res = [0] * (2 * m - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    res[i + j] = (res[i + j] + x * y) % p
        for d in range(2 * m - 2, m - 1, -1):     # reduce by x^m = -sum mod[k] x^k
            c = res[d]
            if c:
                res[d] = 0
                for k in range(m):
                    res[d - m + k] = (res[d - m + k] - c * mod[k]) % p
        return res[:m]

    def to_vec(x):
        return [(x // p ** k) % p for k in range(m)]

    def to_int(v):
        return sum(c * p ** k for k, c in enumerate(v))

    q = p ** m
    for tail in itertools.product(range(p), repeat=m):
        mod = list(tail)                           # x^m + mod[m-1] x^(m-1) + ... + mod[0]
        if m == 1:
            mod = [0]
            break
        # irreducible iff x^(p^m) = x and x^(p^(m/r)) != x for prime r | m ... simple test: no roots/factors by brute force
        ok = True
        # brute-force: the multiplicative group of F_p[x]/(f) has order q-1 iff some element has order q-1 -> field
        # cheaper: f irreducible iff gcd tests; here sizes are tiny, so test that every nonzero element is invertible
        elems = [to_vec(x) for x in range(1, q)]
        one = to_vec(1)
        prods = set()
        a = to_vec(p if m > 1 else 1)              # the class of x
        # x has an inverse and the ring has no zero divisors: check x*y != 0 for all nonzero x, y (q <= 256: fine)
        for x in elems:
            for y in elems:
                if all(c == 0 for c in polymulmod(x, y, mod)):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            break
    add = [[to_int([(u + v) % p for u, v in zip(to_vec(x), to_vec(y))]) for y in range(q)] for x in range(q)]
    mul = [[to_int(polymulmod(to_vec(x), to_vec(y), mod)) if m > 1 else (x * y) % p for y in range(q)] for x in range(q)]
    return q, add, mul


def prime_power(q):
    for p in range(2, q + 1):
        if q % p == 0:
            m = 0
            r = q
            while r % p == 0:
                r //= p
                m += 1
            return (p, m) if r == 1 else None


def kappa_graph(q):
    p, m = prime_power(q)
    Q, add, mul = gf(p, 2 * m)
    def power(x, e):
        r = 1
        for _ in range(e):
            r = mul[r][x]
        return r
    mu = [x for x in range(1, Q) if power(x, q + 1) == 1]
    assert len(mu) == q + 1
    neg = [next(y for y in range(Q) if add[x][y] == 0) for x in range(Q)]
    edges = set()
    for x in range(Q):
        for u in mu:
            y = add[x][u]
            edges.add((min(x, y), max(x, y)))
    return Q, sorted(edges), mu, add, mul


def colorable(nv, edges, k, fix=None):
    var = lambda v, c: v * k + c + 1
    s = Solver(name="cadical153")
    top = nv * k
    for v in range(nv):
        s.add_clause([var(v, c) for c in range(k)])
        for c1 in range(k):
            for c2 in range(c1 + 1, k):
                s.add_clause([-var(v, c1), -var(v, c2)])
    for a, b in edges:
        for c in range(k):
            s.add_clause([-var(a, c), -var(b, c)])
    if fix:
        for v, c in fix:
            s.add_clause([var(v, c)])
    r = s.solve()
    model = s.get_model() if r else None
    s.delete()
    if r:
        col = [next(c for c in range(k) if model[var(v, c) - 1] > 0) for v in range(nv)]
        assert all(col[a] != col[b] for a, b in edges)
    return r


def clique_fix(nv, edges, size):
    """a clique of the given size (greedy search) to pin colours (symmetry breaking)."""
    adj = [set() for _ in range(nv)]
    for a, b in edges:
        adj[a].add(b); adj[b].add(a)
    best = []
    def grow(cl, cand):
        nonlocal best
        if len(cl) > len(best):
            best = cl[:]
        if len(best) >= size:
            return
        for v in sorted(cand):
            grow(cl + [v], cand & adj[v])
            if len(best) >= size:
                return
    grow([0], adj[0])
    return best


print("1. kappa(q) = chi(Cay(F_(q^2), mu_(q+1)))")
expect = {2: 4, 3: 3, 4: 4, 5: 4, 7: 4, 8: 4, 9: 3, 11: 5, 16: 4}
got = {}
for q in sorted(expect):
    Q, E, mu, add, mul = kappa_graph(q)
    k = expect[q]
    cl = clique_fix(Q, E, k - 1)
    fix = [(v, i) for i, v in enumerate(cl[: k - 1])]
    up = colorable(Q, E, k, fix)
    down = colorable(Q, E, k - 1, fix[: k - 1] if len(cl) >= k - 1 else None)
    got[q] = (up, down)
    print(f"   q={q:2d}: |V| = {Q}, degree {q+1}, clique found {len(cl)}; {k}-colourable: {up}; {k-1}-colourable: {down}", flush=True)
check(all(u and not d for u, d in got.values()), "kappa(q) = " + ", ".join(f"{expect[q]} (q={q})" for q in sorted(expect)))

# knight torus
Q, E, mu, add, mul = kappa_graph(7)                      # F_49 = F_7[x]/(f); find the knight set as alpha*mu_8
knight = {((a % 7), (b % 7)) for a, b in [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]}
circle = {(a, b) for a in range(7) for b in range(7) if (a * a + b * b) % 7 == 5}
check(knight == circle, "the 7x7 knight set is the circle x^2 + y^2 = 5 in F_7^2 (THM-4552), i.e. alpha*mu_8 in F_7(i) "
      "with N(alpha) = 5; so G_7 = Cay(F_49, mu_8) up to the additive automorphism z -> alpha z, and chi(G_7) = kappa(7) = 4")

print("2. Lemma R in Q(sqrt-3, sqrt-D) at the prime over 3 (D = 7, 19, 43, 67, 163: N = (D+1)/4 = 2 mod 3)")


class BQ:
    """a + b r3 + c rD + d r3 rD with r3^2 = -3, rD^2 = -D (rationals)."""
    def __init__(s, D, a, b=0, c=0, d=0):
        s.D, s.v = D, (Fr(a), Fr(b), Fr(c), Fr(d))

    def __add__(s, o):
        return BQ(s.D, *[x + y for x, y in zip(s.v, o.v)])

    def __mul__(s, o):
        a, b, c, d = s.v
        e, f, g, h = o.v
        D = s.D
        # basis products: r3^2 = -3, rD^2 = -D, (r3 rD)^2 = 3D, r3*rD = w, r3*w = -3 rD, rD*w = -D r3
        return BQ(D,
                  a * e - 3 * b * f - D * c * g + 3 * D * d * h,
                  a * f + b * e - D * c * h - D * d * g,
                  a * g + c * e - 3 * b * h - 3 * d * f,
                  a * h + d * e + b * g + c * f)

    def conj(s):            # complex conjugation: r3 -> -r3, rD -> -rD, r3 rD -> r3 rD
        a, b, c, d = s.v
        return BQ(s.D, a, -b, -c, d)

    def inv(s):
        # multiply by the three Galois conjugates: sigma3 (r3 -> -r3), sigmaD (rD -> -rD), and both
        a, b, c, d = s.v
        s3 = BQ(s.D, a, -b, c, -d)
        sD = BQ(s.D, a, b, -c, -d)
        sb = BQ(s.D, a, -b, -c, d)
        num = s3 * sD * sb
        n = (s * num).v
        assert n[1] == n[2] == n[3] == 0 and n[0] != 0
        return BQ(s.D, *[x / n[0] for x in num.v])


rng = random.Random(5)
good = True
for D in (7, 19, 43, 67, 163):
    for trial in range(400):
        z = BQ(D, *[rng.randint(-6, 6) for _ in range(4)])
        if all(x == 0 for x in z.v):
            continue
        u = z * z.conj().inv()
        nrm = u * u.conj()
        good &= nrm.v == (1, 0, 0, 0)
        a, b, c, d = u.v
        good &= all(x.denominator % 3 != 0 for x in u.v)
        # residue a + c s in F_9 = F_3[s]/(s^2 + 1): s^2 = -D = -1 mod 3 since D = 1 mod 3
        ra, rc = int(a.numerator * pow(a.denominator, -1, 3)) % 3, int(c.numerator * pow(c.denominator, -1, 3)) % 3
        # (ra + rc s)^4 = 1 in F_9 ?
        def m9(x, y):
            return ((x[0] * y[0] - x[1] * y[1]) % 3, (x[0] * y[1] + x[1] * y[0]) % 3)
        x2 = m9((ra, rc), (ra, rc))
        good &= m9(x2, x2) == (1, 0)
check(good, "2000 unit vectors u = z/conj(z): N(u) = 1, coordinates 3-integral, residue in mu_4 subset F_9 (so chi <= kappa(3) = 3)")

print("3. Erdos 592 on [t]^2: fully gap-determined vs row-invariant")


def lexpos(d):
    for x in d:
        if x:
            return x > 0
    return False


def binary_subgrids(t, n):
    """leaf tuples (lex-sorted) of every binary subgrid of [t]^n."""
    def rec(prefix, level):
        if level == n:
            yield [tuple(prefix)]
            return
        for a, b in itertools.combinations(range(t), 2):
            for L in rec(prefix + [a], level + 1):
                for R in rec(prefix + [b], level + 1):
                    yield L + R
    # the above re-chooses children independently per node: implement as product over nodes
    def node(prefix, level):
        if level == n:
            return [[tuple(prefix)]]
        out = []
        for a, b in itertools.combinations(range(t), 2):
            for L in node(prefix + [a], level + 1):
                for R in node(prefix + [b], level + 1):
                    out.append(L + R)
        return out
    return node([], 0)


def solve_game(t, n, mode):
    pts = list(itertools.product(range(t), repeat=n))
    var = {}
    def v(key):
        if key not in var:
            var[key] = len(var) + 1
        return var[key]
    def edge(x, y):
        if y < x:
            x, y = y, x
        if mode == "full":
            return v(("g",) + tuple(b - a for a, b in zip(x, y)))
        # row-invariant (n = 2): R(c, c') within a row, B_g(c, c') across rows
        if x[0] == y[0]:
            return v(("R", x[1], y[1]))
        return v(("B", y[0] - x[0], x[1], y[1]))
    s = Solver(name="cadical153")
    for x, y, z in itertools.combinations(sorted(pts), 3):
        s.add_clause([-edge(x, y), -edge(y, z), -edge(x, z)])
    for leaves in binary_subgrids(t, n):
        s.add_clause([edge(a, b) for a, b in itertools.combinations(leaves, 2)])
    r = s.solve()
    s.delete()
    return r


full = {t: solve_game(t, 2, "full") for t in (3, 4)}
row = {t: solve_game(t, 2, "row") for t in (4, 5)}
print(f"   fully gap-determined: SAT at t=3: {full[3]}, t=4: {full[4]};  row-invariant: SAT at t=4: {row[4]}, t=5: {row[5]}")
check(full[3] and not full[4] and row[4] and not row[5],
      "fully gap-determined cutoff 4 (Q_inv(2,3) SAT, Q_inv(2,4) UNSAT); row-invariant cutoff 5 (THM-453 G): different games")

print("4. Schur triples with 3-smooth summands")


def smooth(x):
    while x % 2 == 0:
        x //= 2
    while x % 3 == 0:
        x //= 3
    return x == 1


def least_forcing(allow_equal):
    for N in range(2, 30):
        S = [a for a in range(1, N + 1) if smooth(a)]
        trip = [(a, b, a + b) for a in S for b in S if (a < b or (allow_equal and a == b)) and a + b <= N]
        s = Solver(name="cadical153")
        for a, b, c in trip:
            s.add_clause([a, b, c])
            s.add_clause([-a, -b, -c])
        r = s.solve()
        s.delete()
        if not r:
            return N
    return None


f_eq, f_neq = least_forcing(True), least_forcing(False)
check(f_eq == 5 and f_neq == 13, f"least N forcing a monochromatic {{a, b, a+b}} (a, b 3-smooth) in every 2-colouring of [N]: "
      f"{f_eq} with a = b allowed, {f_neq} with a != b")
S = sorted(2 ** i * 3 ** j for i in range(41) for j in range(26) if 2 ** i * 3 ** j <= 10 ** 12)
Sset = set(S)
def Om(x):
    c = 0
    while x % 2 == 0:
        x //= 2; c += 1
    while x % 3 == 0:
        x //= 3; c += 1
    return c
prim = sorted({(a // math.gcd(a, b), b // math.gcd(a, b)) for a in S for b in S if a <= b and a + b in Sset})
mono = [(a, b) for a in S for b in S if a < b and a + b in Sset and Om(a) % 2 == Om(b) % 2 == Om(a + b) % 2]
check(prim == [(1, 1), (1, 2), (1, 3), (1, 8)] and not mono,
      f"3-smooth solutions of x + y = z below 10^12 reduce to {prim}; the colouring Omega mod 2 has no monochromatic Schur triple in S")
print("ALL CHECKS PASSED" if OK else "SOME CHECK FAILED")
