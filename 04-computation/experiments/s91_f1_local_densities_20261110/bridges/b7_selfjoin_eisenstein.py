#!/usr/bin/env python3
"""S91 bridges (7): self-joins of F1 and a search for the Eisenstein form x^2 + xy + y^2 in the join calculus.

Rooted state of the one-apex tournament R(h): (x, y) = (2 psi(h) + 1, 2h + 1) (b2_1039.py checks 1a, 1b).

A. Self-joins.  For h1 = h2 = h the F1 value is x^2 + 2^s y^2 - [s = 0]: the norm form of Z[sqrt(-2^s)] at the state.
   Which of 3, 5, 7, 13, 19, 27, 1039 are self-join values; the 2-adic obstruction (never 7 mod 8).

B. Core gluings (the search universe U1).  Left piece R(h1) = chain C1 + apex a, right piece R(h2)^op = chain C2 +
   apex b, and a set Z of extra vertices.  Arcs: C1 -> everything outside the left piece, everything outside the right
   piece -> C2, and an ARBITRARY tournament K on the core {a, b} u Z.  (F1 is the case K = a -> z_1 -> ... -> z_s -> b
   transitive except b -> a.)  For such a gluing the Hamiltonian path count is a polynomial
        N = c0 x1 x2 + c1 x1 y2 + c2 y1 x2 + c3 y1 y2 + c4 x1 + c5 y1 + c6 x2 + c7 y2 + c8
   with integer coefficients depending only on K (sum over the orders of the core vertices; proof sketch below),
   so on the diagonal h1 = h2 it is the binary quadratic form  c0 x^2 + (c1 + c2) x y + c3 y^2  plus lower terms.
   The script (i) checks the polynomial against the explicit tournaments (bitmask DP) for every K with |Z| <= 2 and
   36 pairs (h1, h2), and for 300 random K with |Z| = 3 and 60 with |Z| = 4; (ii) enumerates every K with |Z| <= 4
   (2 + 8 + 64 + 1024 + 32768 labelled cores) and lists the discriminants D = (c1 + c2)^2 - 4 c0 c3; the Eisenstein
   form needs D = -3 m^2.

   Why a polynomial: a chain vertex of C1 can only be preceded by a lower chain vertex or by a, so C1 is met in at
   most two runs, one at the very start of the path and one right after a; the pairs of runs number x1 (a is the
   first core vertex) or y1 (otherwise), and exactly one of them has an empty run after a.  Symmetrically for C2, b.

C. Middle blocks (the search universe U2).  K = (a => M => b, b -> a) for a tournament M: then
        N = H(M) x1 x2 + D(M) y1 y2,     D(M) = sum over U of H(M[U]) H(M - U)   (ordered two-bin covers),
   the rooted form of THM-5124 (4) (two-hub suspension; D is Session197's count).  Diagonal form H(M) x^2 + D(M) y^2,
   a norm form of Q(sqrt(-H(M) D(M))).  The script lists the squarefree kernels of H(M) D(M) over ALL tournaments M
   on n <= 8 vertices (completions of a Hamiltonian path) and says whether 3 occurs.

D. Three-fold cyclic gluing of one-apex pieces (apexes a1 -> a2 -> a3 -> a1, otherwise R(h1) => R(h2) => R(h3)^op):
   the case M = R(h2) of C.  On the full diagonal the count is x^3 + D(R(h)) y^2; the script tests whether D(R(h))
   is a polynomial of degree <= 2 in (x, y) (it is not: no binary form, so no factor x^2 -+ xy + y^2 of x^3 +- y^3).
"""
import itertools
import math
import sys

import numpy as np
from numba import njit, prange

RESULTS = []


def check(name, cond, detail=""):
    RESULTS.append(bool(cond))
    print(("PASS  " if cond else "FAIL  ") + name + (("   [" + detail + "]") if detail else ""), flush=True)


def psi(h):
    if h == 0:
        return 0
    if h % 2 == 0:
        return psi(h // 2) + h // 2
    k = h // 2
    return psi(k) + k // 2 + 1


def xy(h):
    return 2 * psi(h) + 1, 2 * h + 1


@njit(cache=False)
def ham_paths(adj, n):
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), np.int64)
    for i in range(n):
        dp[1 << i, i] = 1
    for mask in range(1, full + 1):
        for last in range(n):
            c = dp[mask, last]
            if c == 0:
                continue
            for nxt in range(n):
                if (mask >> nxt) & 1 == 0 and adj[last, nxt]:
                    dp[mask | (1 << nxt), nxt] += c
    tot = 0
    for last in range(n):
        tot += dp[full, last]
    return tot


def one_apex(h):
    u = 2 * h
    m = u.bit_length()
    adj = np.zeros((m + 1, m + 1), np.uint8)
    for i in range(m):
        for j in range(i + 1, m):
            adj[i, j] = 1
        if (u >> (m - 1 - i)) & 1:
            adj[m, i] = 1
        else:
            adj[i, m] = 1
    return adj, m


OWNER_ODD = [3, 5, 7, 13, 19, 27, 1039]

# ---------------------------------------------------------------- A. self-joins
def selfjoin_solutions(N):
    out = []
    s = 0
    while (1 << s) <= N + 1:
        h = 0
        while True:
            x, y = xy(h)
            if (y * y << s) > N + 1:
                break
            if x * x + (y * y << s) - (1 if s == 0 else 0) == N:
                out.append((h, s))
            h += 1
        s += 1
    return out


sj = {N: selfjoin_solutions(N) for N in OWNER_ODD + [21]}
check("A1 self-join values among 3, 5, 7, 13, 19, 21, 27, 1039: exactly 3 = (h, s) = (0, 1), 5 = (0, 2), 27 = (1, 1)",
      sj[3] == [(0, 1)] and sj[5] == [(0, 2)] and sj[27] == [(1, 1)] and all(not sj[N] for N in (7, 13, 19, 21, 1039)), str(sj))
res = {s: set() for s in range(0, 13)}
for h in range(1 << 12):
    x, y = xy(h)
    for s in range(0, 13):
        res[s].add((x * x + (y * y << s) - (1 if s == 0 else 0)) % 8)
check("A2 self-join values are 1 (s = 0), 3 (s = 1), 5 (s = 2), 1 (s >= 3) modulo 8: never 7 mod 8, and 7 = 1039 = 7 mod 8",
      res[0] == {1} and res[1] == {3} and res[2] == {5} and all(res[s] == {1} for s in range(3, 13)) and 7 % 8 == 7 and 1039 % 8 == 7)
vals = sorted({x * x + (y * y << s) - (1 if s == 0 else 0) for h in range(64) for (x, y) in [xy(h)] for s in range(0, 14)
               if x * x + (y * y << s) <= 4097})
print("      self-join values <= 4096:", vals[:40], "... (%d values)" % len(vals))
check("A3 at a diagonal state x = y the s = 1 self-join is 3 x^2: 3 (h = 0), 27 (h = 1), 75, 243, 867, 3267 (h = 2, 4, 8, 16)",
      [xy(h)[0] ** 2 + 2 * xy(h)[1] ** 2 for h in (0, 1, 2, 4, 8, 16)] == [3, 27, 75, 243, 867, 3267])


# ---------------------------------------------------------------- B. core gluings
def pair_index(r):
    idx = np.zeros((r, r), np.int64)
    k = 0
    for i in range(r):
        for j in range(i + 1, r):
            idx[i, j] = k
            idx[j, i] = k
            k += 1
    return idx, k


@njit(cache=False)
def core_poly(code, r, perms, pidx):
    """coefficients (x1x2, x1y2, y1x2, y1y2, x1, y1, x2, y2, 1) of the core K coded by 'code' (core vertex 0 = a, 1 = b)."""
    K = np.zeros((r, r), np.int64)
    for i in range(r):
        for j in range(i + 1, r):
            if (code >> pidx[i, j]) & 1:
                K[i, j] = 1
            else:
                K[j, i] = 1
    c = np.zeros(9, np.int64)
    for t in range(perms.shape[0]):
        ok = True
        posa = -1
        posb = -1
        for i in range(r):
            if perms[t, i] == 0:
                posa = i
            if perms[t, i] == 1:
                posb = i
        for i in range(r - 1):
            u = perms[t, i]
            v = perms[t, i + 1]
            if u != 0 and v != 1 and K[u, v] == 0:
                ok = False
                break
        if not ok:
            continue
        ax = 1 if posa == 0 else 0
        bx = 1 if posb == r - 1 else 0
        prod = (0 if ax == 1 else 2) + (0 if bx == 1 else 1)
        ia = 4 if ax == 1 else 5
        ib = 6 if bx == 1 else 7
        if posb == posa + 1:
            c[prod] += 1
            c[8] -= K[1, 0]
        else:
            da = 0
            if posa < r - 1:
                da = K[perms[t, posa + 1], 0]
            db = 0
            if posb > 0:
                db = K[1, perms[t, posb - 1]]
            c[prod] += 1
            c[ia] -= db
            c[ib] -= da
            c[8] += da * db
    return c


def eval_poly(c, h1, h2):
    x1, y1 = xy(h1)
    x2, y2 = xy(h2)
    return int(c[0]) * x1 * x2 + int(c[1]) * x1 * y2 + int(c[2]) * y1 * x2 + int(c[3]) * y1 * y2 + int(c[4]) * x1 \
        + int(c[5]) * y1 + int(c[6]) * x2 + int(c[7]) * y2 + int(c[8])


def glued(code, r, pidx, h1, h2):
    """explicit tournament: C1, a, Z, b, C2."""
    a1, m1 = one_apex(h1)
    a2, m2 = one_apex(h2)
    a2 = a2.T.copy()
    z = r - 2
    n = m1 + 1 + z + 1 + m2
    apos, bpos = m1, m1 + 1 + z
    gl = [apos, bpos] + [m1 + 1 + i for i in range(z)]          # core index -> global index
    adj = np.zeros((n, n), np.uint8)
    for i in range(m1):                                          # C1 internal and C1 - a
        for j in range(m1):
            adj[i, j] = a1[i, j]
        adj[i, apos] = a1[i, m1]
        adj[apos, i] = a1[m1, i]
    off2 = bpos + 1
    for i in range(m2):                                          # C2 internal and C2 - b
        for j in range(m2):
            adj[off2 + i, off2 + j] = a2[i, j]
        adj[off2 + i, bpos] = a2[i, m2]
        adj[bpos, off2 + i] = a2[m2, i]
    for i in range(m1):                                          # C1 -> everything outside the left piece
        for j in range(m1 + 1, n):
            adj[i, j] = 1
    for j in range(m2):                                          # everything outside the right piece -> C2
        for i in range(m1, bpos):
            adj[i, off2 + j] = 1
    for i in range(r):
        for j in range(i + 1, r):
            if (code >> pidx[i, j]) & 1:
                adj[gl[i], gl[j]] = 1
            else:
                adj[gl[j], gl[i]] = 1
    return adj, n


PERMS = {r: np.array(list(itertools.permutations(range(r))), dtype=np.int64) for r in range(2, 7)}
PIDX = {r: pair_index(r) for r in range(2, 7)}
HS = [0, 1, 2, 3, 5, 6]
ok = True
ncheck = 0
for r in (2, 3, 4):
    pidx, npairs = PIDX[r]
    for code in range(1 << npairs):
        c = core_poly(code, r, PERMS[r], pidx)
        for h1 in HS:
            for h2 in HS:
                adj, n = glued(code, r, pidx, h1, h2)
                ncheck += 1
                if ham_paths(adj, n) != eval_poly(c, h1, h2):
                    ok = False
rng = np.random.default_rng(20261110)
for r, nsamp, hs in ((5, 300, [0, 1, 2, 5]), (6, 60, [0, 1, 3])):
    pidx, npairs = PIDX[r]
    for code in rng.integers(0, 1 << npairs, size=nsamp):
        code = int(code)
        c = core_poly(code, r, PERMS[r], pidx)
        for h1 in hs:
            for h2 in hs:
                adj, n = glued(code, r, pidx, h1, h2)
                ncheck += 1
                if ham_paths(adj, n) != eval_poly(c, h1, h2):
                    ok = False
check("B1 the core polynomial equals the Hamiltonian path count of the explicit gluing (%d tournaments: every core with |Z| <= 2, "
      "random cores with |Z| = 3, 4)" % ncheck, ok)
# F1 cores
okf = True
for s in range(0, 5):
    r = s + 2
    pidx, npairs = PIDX[r]
    code = 0
    order = [0] + list(range(2, r)) + [1]                        # a, z_1, ..., z_s, b
    for i in range(r):
        for j in range(i + 1, r):
            u, v = order[i], order[j]
            fwd = not (u == 0 and v == 1)                        # the arc a -> b is reversed
            lo, hi = min(u, v), max(u, v)
            bit = 1 if ((u == lo) == fwd) else 0
            code |= bit << pidx[lo, hi]
    c = core_poly(code, r, PERMS[r], pidx)
    okf &= list(c) == [1, 0, 0, 2 ** s, 0, 0, 0, 0, -1 if s == 0 else 0]
check("B2 the F1 cores give x1 x2 + 2^s y1 y2 - [s = 0] for s = 0..4", okf)

print("B. census of cores: |Z|, number of cores, distinct diagonal forms (c0, c1 + c2, c3), discriminants")
census = {}
eis = {}
all_disc = {}
for r in (2, 3, 4, 5, 6):
    pidx, npairs = PIDX[r]
    forms = {}
    for code in range(1 << npairs):
        c = core_poly(code, r, PERMS[r], pidx)
        key = (int(c[0]), int(c[1] + c[2]), int(c[3]))
        forms.setdefault(key, []).append((code, tuple(int(v) for v in c)))
    census[r] = forms
    discs = {}
    for (al, be, ga), lst in forms.items():
        D = be * be - 4 * al * ga
        discs.setdefault(D, []).append((al, be, ga))
    all_disc[r] = discs
    neg = sorted(D for D in discs if D < 0)
    eis[r] = sorted(D for D in neg if (-D) % 3 == 0 and math.isqrt((-D) // 3) ** 2 == (-D) // 3)
    print("      |Z| = %d: %6d cores, %4d distinct forms, %4d distinct discriminants, %d negative; Eisenstein type (D = -3 m^2): %s"
          % (r - 2, 1 << npairs, len(forms), len(discs), len(neg), eis[r]))


def kernel(n):
    """squarefree kernel of a positive integer."""
    k, d = 1, 2
    while d * d <= n:
        e = 0
        while n % d == 0:
            n //= d
            e += 1
        if e % 2:
            k *= d
        d += 1
    return k * n


fields = {r: sorted({kernel(-D) for D in all_disc[r] if D < 0}) for r in all_disc}
for r in fields:
    print("      |Z| = %d: imaginary quadratic fields Q(sqrt(-d)) of the definite forms, d = %s" % (r - 2, fields[r][:40])
          + (" ... (%d values)" % len(fields[r]) if len(fields[r]) > 40 else ""))
check("B3 |Z| = 0: forms x^2 (transitive) and x^2 + y^2 (b -> a): discriminants 0, -4", sorted(all_disc[2]) == [-4, 0])
f1 = sorted(census[3])
check("B4 |Z| = 1 (all 8 triangles on a, b, z): diagonal forms (1,0,0), (1,0,2), (1,1,0), (1,1,2), (1,2,1), (1,2,3); "
      "discriminants 0, 1, -7, -8: x^2 + xy + 2y^2 (the norm form of Q(sqrt(-7))) occurs, x^2 + xy + y^2 does not",
      f1 == [(1, 0, 0), (1, 0, 2), (1, 1, 0), (1, 1, 2), (1, 2, 1), (1, 2, 3)] and sorted(all_disc[3]) == [-8, -7, 0, 1], str(f1))
def reduce_form(a, b, c):
    """Gauss reduction of a positive definite form a x^2 + b x y + c y^2: returns (a, b, c) with -a < b <= a <= c, b >= 0 if a = c."""
    while True:
        t = (a - b) // (2 * a)
        b, c = b + 2 * a * t, a * t * t + b * t + c
        if a > c:
            a, b, c = c, -b, a
            continue
        if a == c and b < 0:
            b = -b
        return a, b, c


def arcs_of(code, r):
    pidx, _ = PIDX[r]
    names = "ab" + "pqrs"[:r - 2]
    out = []
    for i in range(r):
        for j in range(i + 1, r):
            out.append("%s->%s" % ((names[i], names[j]) if (code >> pidx[i, j]) & 1 else (names[j], names[i])))
    return ", ".join(out)


def confirm_diagonal(code, r, c, hs=(0, 1, 2, 3, 5)):
    pidx, _ = PIDX[r]
    good = True
    for h in hs:
        adj, n = glued(code, r, pidx, h, h)
        x, y = xy(h)
        good &= ham_paths(adj, n) == c[0] * x * x + (c[1] + c[2]) * x * y + c[3] * y * y + (c[4] + c[6]) * x + (c[5] + c[7]) * y + c[8]
    return good


check("B5 reduction routine: (3, 9, 7) -> (1, 1, 1), (1, 2, 4) -> (1, 0, 3), (3, 9, 19) -> (3, 3, 13), (1, 1, 2) -> (1, 1, 2)",
      reduce_form(3, 9, 7) == (1, 1, 1) and reduce_form(1, 2, 4) == (1, 0, 3) and reduce_form(3, 9, 19) == (3, 3, 13)
      and reduce_form(1, 1, 2) == (1, 1, 2))
print("      disc -7 cores at |Z| = 1 (N as coefficients of x1x2, x1y2, y1x2, y1y2, x1, y1, x2, y2, 1):")
ok7 = True
for code, c in census[3][(1, 1, 2)]:
    print("        %s   N = %s" % (arcs_of(code, 3), c))
    ok7 &= confirm_diagonal(code, 3, c, hs=(0, 1, 2, 3, 5, 6, 9, 10))
check("B6 the two disc -7 cores (z a source or a sink of the core triangle, b -> a): explicit tournaments have "
      "x^2 + xy + 2y^2 - x - y - 1 Hamiltonian paths on the diagonal (h = 0, 1, 2, 3, 5, 6, 9, 10)",
      ok7 and all(c[0] == 1 and c[1] + c[2] == 1 and c[3] == 2 and c[4] + c[6] == -1 and c[5] + c[7] == -1 and c[8] == -1
                  for _, c in census[3][(1, 1, 2)]))
print("      Eisenstein-type forms (D = -3 m^2): |Z|, form (c0, c1 + c2, c3), D, reduced form, number of cores, example")
eis_rows = {}
okE = True
for r in (4, 5, 6):
    rows_ = []
    for key in sorted(census[r]):
        al, be, ga = key
        D = be * be - 4 * al * ga
        if D < 0 and (-D) % 3 == 0 and math.isqrt((-D) // 3) ** 2 == (-D) // 3:
            red = reduce_form(al, be, ga)
            g = math.gcd(math.gcd(red[0], red[1]), red[2])
            multiple = (red[0] // g, red[1] // g, red[2] // g) == (1, 1, 1)
            rows_.append((key, D, red, multiple, len(census[r][key])))
            code, c = census[r][key][0]
            if r <= 5:
                okE &= confirm_diagonal(code, r, c)
                print("        |Z| = %d: %s D = %d reduced %s%s; %d cores; e.g. %s   N = %s"
                      % (r - 2, key, D, red, " = %d * Eisenstein" % g if multiple else "", len(census[r][key]), arcs_of(code, r), c))
    eis_rows[r - 2] = rows_
    if r == 6:
        print("        |Z| = 4: %d Eisenstein-type forms; reduced forms and numbers of forms: %s"
              % (len(rows_), sorted({(row[2], sum(1 for q_ in rows_ if q_[2] == row[2])) for row in rows_})))
check("B7 explicit confirmation (h = 0, 1, 2, 3, 5 on the diagonal) of one core for every Eisenstein-type form with |Z| <= 3", okE)
exact = {zz: [row for row in eis_rows[zz] if row[0] == (1, 1, 1)] for zz in eis_rows}
proper = {zz: [row[0] for row in eis_rows[zz] if row[3]] for zz in eis_rows}
print("      forms equivalent to a multiple of x^2 + xy + y^2:", proper)
check("B8 U1 answer: the literal form x^2 + xy + y^2 occurs for NO core with |Z| <= 4; no Eisenstein-type form at |Z| <= 1; "
      "x^2 + 3y^2 (order Z[sqrt(-3)]) first at |Z| = 2; a form properly equivalent to x^2 + xy + y^2 first at |Z| = 3: "
      "3x^2 + 9xy + 7y^2 = N(x + y, x + 2y) with Z a 3-cycle",
      all(not exact[zz] for zz in exact) and not eis[2] and not eis[3] and proper[2] == [] and (3, 9, 7) in proper[3]
      and all((1, 1, 1) != k_ for r in census for k_ in census[r]))
check("B9 identity 3x^2 + 9xy + 7y^2 = X^2 + XY + Y^2 at (X, Y) = (x + y, x + 2y), tested on a grid",
      all(3 * x * x + 9 * x * y + 7 * y * y == (x + y) ** 2 + (x + y) * (x + 2 * y) + (x + 2 * y) ** 2 for x in range(-6, 7) for y in range(-6, 7)))
cyc = True
for code, c in census[5][(3, 9, 7)]:
    pidx, _ = PIDX[5]
    zarcs = [(code >> pidx[i, j]) & 1 for i, j in ((2, 3), (2, 4), (3, 4))]
    cyc &= zarcs in ([1, 0, 1], [0, 1, 0])                       # p->q, r->p, q->r  or  q->p, p->r, r->q
    cyc &= (code >> pidx[0, 1]) & 1 == 1                          # a -> b (no back arc)
check("B10 all %d cores with the form 3x^2 + 9xy + 7y^2 have Z a cyclic triangle and a -> b" % len(census[5][(3, 9, 7)]), cyc)
lower_free = [(r - 2, key) for r in (4, 5, 6) for key, D, red, mult, cnt in eis_rows[r - 2]
              if any(c[4:] == (0, 0, 0, 0, 0) for _, c in census[r][key])]
print("      Eisenstein-type forms realised with no lower-order terms:", lower_free)


# ---------------------------------------------------------------- C. middle blocks
@njit(parallel=True, cache=False)
def block_census(n, nthreads):
    pi = np.zeros(64, np.int64)
    pj = np.zeros(64, np.int64)
    m = 0
    for i in range(n):
        for j in range(i + 2, n):
            pi[m] = i
            pj[m] = j
            m += 1
    total = 1 << m
    full = (1 << n) - 1
    Hs = np.zeros(total, np.int64)
    Cs = np.zeros(total, np.int64)
    chunk = total // nthreads + 1
    for t in prange(nthreads):
        out = np.zeros(n, np.int64)
        dp = np.zeros((1 << n, n), np.int64)
        paths = np.zeros(1 << n, np.int64)
        for code in range(t * chunk, min(total, (t + 1) * chunk)):
            for i in range(n):
                out[i] = 0
            for i in range(n - 1):
                out[i] |= 1 << (i + 1)
            for b in range(m):
                if (code >> b) & 1:
                    out[pi[b]] |= 1 << pj[b]
                else:
                    out[pj[b]] |= 1 << pi[b]
            for mask in range(1, full + 1):
                for last in range(n):
                    dp[mask, last] = 0
            for i in range(n):
                dp[1 << i, i] = 1
            paths[0] = 1
            for mask in range(1, full + 1):
                tot = 0
                for last in range(n):
                    c = dp[mask, last]
                    if c == 0:
                        continue
                    tot += c
                    avail = out[last] & (full ^ mask)
                    while avail:
                        low = avail & (-avail)
                        nxt = 0
                        while (low >> nxt) != 1:
                            nxt += 1
                        dp[mask | low, nxt] += c
                        avail ^= low
                paths[mask] = tot
            cc = 0
            for mask in range(0, full + 1):
                cc += paths[mask] * paths[full ^ mask]
            Hs[code] = paths[full]
            Cs[code] = cc
    return Hs, Cs


@njit(cache=False)
def h_and_d(adj, n):
    """(H(M), D(M)) for the digraph adj on n vertices: D = sum over U of paths(U) paths(complement), paths(empty) = 1."""
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), np.int64)
    paths = np.zeros(1 << n, np.int64)
    for i in range(n):
        dp[1 << i, i] = 1
    paths[0] = 1
    for mask in range(1, full + 1):
        tot = 0
        for last in range(n):
            c = dp[mask, last]
            if c == 0:
                continue
            tot += c
            for nxt in range(n):
                if (mask >> nxt) & 1 == 0 and adj[last, nxt]:
                    dp[mask | (1 << nxt), nxt] += c
        paths[mask] = tot
    d = 0
    for mask in range(full + 1):
        d += paths[mask] * paths[full ^ mask]
    return paths[full], d


def block_join(h1, M, s, h2):
    """R(h1) => M => (transitive chain of s vertices) => R(h2)^op, with the single back arc b -> a between the apexes."""
    a1, m1 = one_apex(h1)
    a2, m2 = one_apex(h2)
    a2 = a2.T.copy()
    k = M.shape[0]
    n = (m1 + 1) + k + s + (m2 + 1)
    adj = np.zeros((n, n), np.uint8)
    for i in range(n):
        for j in range(i + 1, n):
            adj[i, j] = 1                                        # everything forward
    adj[:m1 + 1, :m1 + 1] = a1
    adj[m1 + 1:m1 + 1 + k, m1 + 1:m1 + 1 + k] = M
    off = m1 + 1 + k + s
    adj[off:, off:] = a2
    ap1, ap2 = m1, off + m2
    adj[ap1, ap2] = 0
    adj[ap2, ap1] = 1
    return adj, n


def decode_path_code(n, code):
    prs = [(i, j) for i in range(n) for j in range(i + 2, n)]
    adj = np.zeros((n, n), np.uint8)
    for i in range(n - 1):
        adj[i, i + 1] = 1
    for b_, (i, j) in enumerate(prs):
        if (code >> b_) & 1:
            adj[i, j] = 1
        else:
            adj[j, i] = 1
    return adj


print("C. middle blocks M: pairs (H(M), D(M)) over all tournaments on n vertices; squarefree kernel d of H(M) D(M)")
pairs_by_n = {1: {(1, 2)}}
arrays = {}
for n in range(2, 9):
    Hs, Cs = block_census(n, 3)
    arrays[n] = (Hs, Cs)
    pairs_by_n[n] = {(int(a_), int(b_)) for a_, b_ in zip(Hs.tolist(), Cs.tolist())} if n <= 7 else \
        {(int(a_), int(b_)) for a_, b_ in np.unique(np.stack([Hs, Cs], axis=1), axis=0)}
kern = {}
three, six = [], []
for n in sorted(pairs_by_n):
    ks = sorted({kernel(hh * cc) for hh, cc in pairs_by_n[n]})
    kern[n] = ks
    hits3 = sorted((hh, cc) for hh, cc in pairs_by_n[n] if kernel(hh * cc) == 3)
    hits6 = sorted((hh, cc) for hh, cc in pairs_by_n[n] if kernel(hh * cc) == 6)
    three.extend((n, hh, cc) for hh, cc in hits3)
    six.extend((n, hh, cc) for hh, cc in hits6)
    print("      n = %d: %5d distinct pairs; kernels: %s%s; kernel 3: %s; kernel 6: %s"
          % (n, len(pairs_by_n[n]), ks[:14], " ..." if len(ks) > 14 else "", hits3[:5], hits6[:5]))
check("C1 transitive blocks give (1, 2^n); the 3-cycle (3, 12); the strong 4-tournament (5, 32) (ratio 32/5 of the S82 note)",
      all((1, 2 ** n) in pairs_by_n[n] for n in range(1, 9)) and (3, 12) in pairs_by_n[3] and (5, 32) in pairs_by_n[4])
check("C2 order 5: the ratios D/H include 36/5, 100/11, 28/3, 64/5 (S82 note)",
      all(any(cc * q_ == p_ * hh for hh, cc in pairs_by_n[5]) for p_, q_ in ((36, 5), (100, 11), (28, 3), (64, 5))))
# the middle-block formula as a special case of B
okm = True
nm_checked = 0
for r in (3, 4, 5, 6):
    z = r - 2
    pidx, npairs = PIDX[r]
    for sub in range(1 << (z * (z - 1) // 2)):
        code = 0
        k = 0
        M = np.zeros((z, z), np.uint8)
        for i in range(2, r):
            for j in range(i + 1, r):
                if (sub >> k) & 1:
                    code |= 1 << pidx[i, j]
                    M[i - 2, j - 2] = 1
                else:
                    M[j - 2, i - 2] = 1
                k += 1
        for i in range(2, r):
            code |= 1 << pidx[0, i]                              # a -> M;  the bits of (b, z) stay 0: z -> b;  bit (a, b) = 0: b -> a
        c = core_poly(code, r, PERMS[r], pidx)
        hM, dM = h_and_d(M, z)
        nm_checked += 1
        okm &= list(c) == [int(hM), 0, 0, int(dM), 0, 0, 0, 0, 0]
check("C3 for K = (a => M => b, b -> a) the polynomial is H(M) x1 x2 + D(M) y1 y2 (all %d labelled M on 1..4 vertices)" % nm_checked, okm)
check("C4 a sink or source added to M doubles D and keeps H: so padding by s vertices gives H(M) x1 x2 + 2^s D(M) y1 y2, "
      "and the field is Q(sqrt(-3)) iff the kernel of H(M) D(M) 2^s is 3",
      all((hh, 2 * cc) in pairs_by_n[n + 1] for n in range(1, 8) for hh, cc in pairs_by_n[n]))
first3 = min(n for n, _, _ in three) if three else None
first6 = min(n for n, _, _ in six) if six else None
check("C5 U2 answer: kernel 3 first occurs at n = %s and kernel 6 at n = %s; no block on <= 5 vertices has kernel 3 or 6" % (first3, first6),
      first3 == 7 and first6 == 6 and not any(kernel(hh * cc) in (3, 6) for n in range(1, 6) for hh, cc in pairs_by_n[n]),
      "kernel 6 at n = 6: %s; kernel 3 at n = 7: %s" % ([(hh, cc) for n, hh, cc in six if n == 6], [(hh, cc) for n, hh, cc in three if n == 7]))
# explicit Eisenstein-type middle block: the first 6-vertex block with kernel 6, padded by one vertex
Hs6, Cs6 = arrays[6]
code6 = next(code for code in range(len(Hs6)) if kernel(int(Hs6[code]) * int(Cs6[code])) == 6)
M6 = decode_path_code(6, code6)
hM6, dM6 = (int(v) for v in h_and_d(M6, 6))
masks6 = [int(sum(1 << j for j in range(6) if M6[i, j])) for i in range(6)]
scores6 = sorted(bin(m_).count("1") for m_ in masks6)
print("      explicit block M6: out-neighbour masks %s (bit j of masks[i] set iff i -> j), scores %s, (H, D) = (%d, %d), H D = %d = 6 * %d^2"
      % (masks6, scores6, hM6, dM6, hM6 * dM6, math.isqrt(hM6 * dM6 // 6)))
okx = True
for h in (0, 1, 2):
    adj, n = block_join(h, M6, 1, h)
    x, y = xy(h)
    okx &= ham_paths(adj, n) == hM6 * x * x + 2 * dM6 * y * y
check("C6 explicit tournaments R(h) => M6 => z => R(h)^op (back arc between the apexes) have H(M6) x^2 + 2 D(M6) y^2 Hamiltonian paths "
      "(h = 0, 1, 2), a form of discriminant -3 m^2", okx and kernel(hM6 * 2 * dM6) == 3,
      "%d x^2 + %d y^2" % (hM6, 2 * dM6))
P7 = np.array([[1 if ((j - i) % 7) in (1, 2, 4) else 0 for j in range(7)] for i in range(7)], dtype=np.uint8)
hP7, dP7 = (int(v) for v in h_and_d(P7, 7))
adj, n = block_join(1, P7, 1, 1)
check("C7 Paley P_7: (H, D) = (189, 2016), H D = 6 * 252^2; with one padding vertex the join form is 189 x^2 + 4032 y^2 = 63 (3 x^2 + (8y)^2); "
      "explicit count at h = 1: 37989", (hP7, dP7) == (189, 2016) and 189 * 2016 == 6 * 252 ** 2 and 189 * 4032 == 3 * 504 ** 2
      and 189 == 63 * 3 and 4032 == 63 * 64 and int(ham_paths(adj, n)) == 189 * 9 + 4032 * 9 == 37989)
print("      kernels d of H(M) D(M) for n <= 8 that are <= 30:", sorted({k_ for n in kern for k_ in kern[n] if k_ <= 30}))
seven = sorted({(n, hh, cc) for n in pairs_by_n for hh, cc in pairs_by_n[n] if kernel(hh * cc) in (7, 14)})
print("      blocks with kernel 7 or 14 (field Q(sqrt(-7)) after padding):", seven[:6])
nineteen = sorted({(n, hh, cc) for n in pairs_by_n for hh, cc in pairs_by_n[n] if kernel(hh * cc) in (19, 38)})
print("      blocks with kernel 19 or 38:", nineteen[:6])

# ---------------------------------------------------------------- D. three-fold cyclic gluing of one-apex pieces
rows = []
for h in range(0, 512):
    adj, m = one_apex(h)
    hR, dR = h_and_d(adj, m + 1)
    x, y = xy(h)
    assert int(hR) == x
    rows.append((h, x, y, int(dR)))
okd = True
for h1, h2 in ((1, 1), (2, 1), (1, 3), (0, 2), (3, 2), (2, 10)):
    adjM, mM = one_apex(h2)
    adj, n = block_join(h1, adjM, 0, h1)
    x1, y1 = xy(h1)
    x2, _ = xy(h2)
    okd &= ham_paths(adj, n) == x2 * x1 * x1 + rows[h2][3] * y1 * y1
check("D1 three-fold cyclic gluing (apexes a1 -> a2 -> a3 -> a1): H = x(h2) x1 x3 + D(R(h2)) y1 y3 on explicit tournaments (6 cases, h3 = h1)", okd)
print("      h, x, y, D(R(h)) for h < 12:", rows[:12])
import sympy
sub = rows[:64]
Mat = sympy.Matrix([[1, x, y, x * x, x * y, y * y] for _, x, y, _ in sub])
aug = Mat.row_join(sympy.Matrix([d for _, _, _, d in sub]))
check("D2 D(R(h)) is NOT a polynomial of degree <= 2 in (x, y) on h < 64 (rank test): the fully diagonal three-fold count "
      "x^3 + D(R(h)) y^2 is not a binary form, so x^3 + y^3 = (x + y)(x^2 - xy + y^2) has no analogue here", aug.rank() > Mat.rank(),
      "rank %d vs %d" % (Mat.rank(), aug.rank()))
same_x = {}
for h, x, y, d in rows:
    same_x.setdefault(x, set()).add(d)
multi = sorted((x, sorted(v)) for x, v in same_x.items() if len(v) > 1)
check("D3 D(R(h)) is not a function of x alone (same x, different D)", len(multi) > 0, str(multi[:3]))
eis_h = [(h, x, d, kernel(x * d)) for h, x, y, d in rows if kernel(x * d) in (3, 6)]
print("      one-apex middle blocks R(h), h < 512, with kernel of x D in {3, 6}: %d; first: %s" % (len(eis_h), eis_h[:6]))
k19 = [(h, x, d) for h, x, y, d in rows if kernel(x * d) in (19, 38)][:4]
print("      one-apex middle blocks with kernel 19 or 38 (first four):", k19)

if all(RESULTS):
    print("ALL CHECKS PASSED (%d checks)" % len(RESULTS))
    sys.exit(0)
print("CHECK FAILED (%d of %d)" % (RESULTS.count(False), len(RESULTS)))
sys.exit(1)
