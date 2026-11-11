#!/usr/bin/env python3
"""S91 bridges (4): the strip automaton of J. Allikvere, arXiv:2607.28813 (multiset dimension of P_3 [x] P_n),
rebuilt from the paper's local conditions, and the helper's unverified claims about its spectrum.

Paper's objects.  Column word t_x = (p_x, m_x, q_x) in {0,1}^3 (landmarks in rows 0, 1, 2), t_0 = t_(n+1) = 000.
A sequence t_1 .. t_n is ADMISSIBLE if at every column x
   (A) (p_x, p_(x-1) + p_(x+1)) != (q_x, q_(x-1) + q_(x+1)),
   (B) p_x != m_x  or  q_(x-1) + q_x + q_(x+1) >= 1,
   (C) q_x != m_x  or  p_(x-1) + p_x + p_(x+1) >= 1.
Min-plus transfer:  V_1(000, v) = |v|,  (T V)(v, w) = min over u with (A, B, C) at the middle of (u, v, w) of V(u, v) + |w|.

Checked here, exactly (integers and rationals only):
  1. the paper: V_5 = V_4 + 1; V_4 finite at 62 of 64 states; W(3) = 3, W(4) = 4; 184 tight transitions; 320
     minimum-weight 4-prefixes; 12 accepting states; a recurrent core L of 19 states, strongly connected; the second
     components of L are 000, 001, 010, 100, 101; 106 minimum 4-prefixes end in L; A w = 3 w with the entries of the
     paper's appendix (9 ten times, 6 three times, 5 twice, 2 four times);
  2. the conditions are necessary, and the paper's resolving counts: among the n-subsets of P_3 [x] P_n the multiset
     resolving ones number 0, 32, 320 for n = 5, 6, 7 (exhaustive over all C(3n, n) subsets) and each is admissible;
     among the weight-8 admissible sequences 1328 resolve;
  3. helper claim (d1): the five column types are the independent sets of the path P_3 and L is 19 of the 25 pairs;
  4. helper claim (d2): the characteristic polynomial of the core is x^14 (x-3)(x-1)(x+1)(x^2+x+2)
     (two independent exact computations: sympy, and Bareiss determinants at 20 integers + Lagrange interpolation);
  5. helper claim (d3): the number a(n) of weight-n admissible sequences (direct dynamic programme on the
     conditions, no automaton; brute force for n <= 7) satisfies a(n) = 2a(n-1) + a(n-2) + 6a(n-3) and
     42 a(n) = 32*3^n - 9 V_n - 21 U_n; the automaton path counts agree for 5 <= n <= 200;
  6. helper claim (d4): for the Lucas sequences of x^2 + x + 2, V_n^2 + 7 U_n^2 = 2^(n+2), and |U_n| = 1 exactly
     for n = 1, 2, 3, 5, 13 among n <= 20000 (for all n this is the Ramanujan-Nagell theorem, Nagell 1948: CITED).
"""
import itertools
import sys
from fractions import Fraction

import sympy

RESULTS = []


def check(name, cond, detail=""):
    RESULTS.append(bool(cond))
    print(("PASS  " if cond else "FAIL  ") + name + (("   [" + detail + "]") if detail else ""), flush=True)


INF = 10 ** 9
COLS = [(p, m, q) for p in (0, 1) for m in (0, 1) for q in (0, 1)]
ZERO = (0, 0, 0)
nm = {t: "%d%d%d" % t for t in COLS}


def legal(u, v, w):
    """conditions (A), (B), (C) at the middle column v of the consecutive columns u, v, w."""
    (pu, mu, qu), (pv, mv, qv), (pw, mw, qw) = u, v, w
    if (pv, pu + pw) == (qv, qu + qw):
        return False
    if pv == mv and qu + qv + qw == 0:
        return False
    if qv == mv and pu + pv + pw == 0:
        return False
    return True


LEGAL = {(u, v, w): legal(u, v, w) for u in COLS for v in COLS for w in COLS}
wt = {t: sum(t) for t in COLS}


def transfer(V):
    out = {}
    for v in COLS:
        for w in COLS:
            best = INF
            for u in COLS:
                if V[(u, v)] < INF and LEGAL[(u, v, w)]:
                    best = min(best, V[(u, v)] + wt[w])
            out[(v, w)] = best
    return out


V = {1: {(u, v): (wt[v] if u == ZERO else INF) for u in COLS for v in COLS}}
for L in range(2, 8):
    V[L] = transfer(V[L - 1])


def Wmin(n):
    return min(V[n][(u, v)] for u in COLS for v in COLS if V[n][(u, v)] < INF and LEGAL[(u, v, ZERO)])


V4 = V[4]
finite = [s for s in V4 if V4[s] < INF]
check("1a V_5 = V_4 + 1 on all 64 states (and V_6 = V_5 + 1, V_7 = V_6 + 1)",
      all((V[5][s] == V4[s] + 1) if V4[s] < INF else V[5][s] >= INF for s in V4)
      and all((V[6][s] == V[5][s] + 1) if V[5][s] < INF else V[6][s] >= INF for s in V4)
      and all((V[7][s] == V[6][s] + 1) if V[6][s] < INF else V[7][s] >= INF for s in V4))
check("1b V_4 is finite at 62 of 64 states with values 3..8", len(finite) == 62 and sorted({V4[s] for s in finite}) == [3, 4, 5, 6, 7, 8],
      "infinite at %s" % [nm[s[0]] + "," + nm[s[1]] for s in V4 if V4[s] >= INF])
check("1c W(3) = 3, W(4) = 4, W(5) = 5, W(6) = 6, W(7) = 7", [Wmin(n) for n in (3, 4, 5, 6, 7)] == [3, 4, 5, 6, 7])
tight = [((u, v), (v, w)) for u in COLS for v in COLS for w in COLS
         if V4[(u, v)] < INF and V4[(v, w)] < INF and LEGAL[(u, v, w)] and wt[w] + V4[(u, v)] - V4[(v, w)] == 1]
check("1d number of tight transitions = 184", len(tight) == 184, str(len(tight)))
prefixes = [t for t in itertools.product(COLS, repeat=4)
            if LEGAL[(ZERO, t[0], t[1])] and LEGAL[(t[0], t[1], t[2])] and LEGAL[(t[1], t[2], t[3])]
            and sum(wt[c] for c in t) == V4[(t[2], t[3])]]
check("1e number of minimum-weight admissible 4-prefixes = 320", len(prefixes) == 320, str(len(prefixes)))
accept = [s for s in finite if V4[s] == 4 and LEGAL[(s[0], s[1], ZERO)]]
check("1f number of accepting states (V_4 = 4, legal last column) = 12", len(accept) == 12, str(len(accept)))
succ, pred = {}, {}
for a, b in tight:
    succ.setdefault(a, set()).add(b)
    pred.setdefault(b, set()).add(a)


def closure(start, nb):
    seen = set(start)
    stack = list(start)
    while stack:
        s = stack.pop()
        for t in nb.get(s, ()):
            if t not in seen:
                seen.add(t)
                stack.append(t)
    return seen


core = sorted(closure({(t[2], t[3]) for t in prefixes}, succ) & closure(set(accept), pred))
check("1g the recurrent core L has 19 states", len(core) == 19, str(len(core)))
idx = {s: i for i, s in enumerate(core)}
n19 = len(core)
Aadj = [[0] * n19 for _ in range(n19)]
for a, b in tight:
    if a in idx and b in idx:
        Aadj[idx[a]][idx[b]] = 1
strong = all(closure({s}, {a: {b for b in succ.get(a, ()) if b in idx} for a in core}) >= set(core) for s in core)
check("1h L is strongly connected (within L) and contains loops", strong and any(Aadj[i][i] for i in range(n19)),
      "loops at %s; %d edges" % ([nm[core[i][0]] + "," + nm[core[i][1]] for i in range(n19) if Aadj[i][i]], sum(map(sum, Aadj))))
types = sorted({s[1] for s in core})
FIVE = [(0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0), (1, 0, 1)]
check("1i second components of L are exactly 000, 001, 010, 100, 101; all 12 accepting states lie in L",
      types == FIVE and all(s in idx for s in accept))
pre_in = [t for t in prefixes if (t[2], t[3]) in idx]
check("1j 106 minimum 4-prefixes end in L, and they use only the five types",
      len(pre_in) == 106 and all(c in FIVE for t in pre_in for c in t), str(len(pre_in)))
M = sympy.Matrix(Aadj)
ns = (M - 3 * sympy.eye(n19)).nullspace()
vec = ns[0]
den = sympy.ilcm(*[sympy.Rational(x).q for x in vec])
vec = [int(x * den) for x in vec]
g = 0
for x in vec:
    g = sympy.igcd(g, x)
vec = [x // g for x in vec]
if vec[0] < 0:
    vec = [-x for x in vec]
check("1k A w = 3 w has a one-dimensional solution space; w is positive with entries 9 (x10), 6 (x3), 5 (x2), 2 (x4)",
      len(ns) == 1 and sorted(vec) == sorted([9] * 10 + [6] * 3 + [5] * 2 + [2] * 4), str(sorted(vec)))


# ---------------------------------------------------------------- 2. necessity and resolving counts
def multiset_resolving(n, S):
    """S: set of (x, y) in the king strip P_3 [x] P_n; Chebyshev distance."""
    seen = set()
    for x in range(n):
        for y in range(3):
            key = tuple(sorted(max(abs(x - a), abs(y - b)) for a, b in S))
            if key in seen:
                return False
            seen.add(key)
    return True


def columns_of(n, S):
    return [tuple(1 if (x, y) in S else 0 for y in range(3)) for x in range(n)]


def admissible(seq):
    ext = [ZERO] + list(seq) + [ZERO]
    return all(LEGAL[(ext[i - 1], ext[i], ext[i + 1])] for i in range(1, len(ext) - 1))


res_counts = {}
necessity = True
for n in (5, 6, 7):
    verts = [(x, y) for x in range(n) for y in range(3)]
    cnt = 0
    for S in itertools.combinations(verts, n):
        if multiset_resolving(n, S):
            cnt += 1
            necessity &= admissible(columns_of(n, set(S)))
    res_counts[n] = cnt
check("2a multiset resolving n-subsets of P_3 [x] P_n (all C(3n, n) subsets): 0, 32, 320 for n = 5, 6, 7", res_counts == {5: 0, 6: 32, 7: 320},
      str(res_counts))
check("2b every one of them is an admissible sequence (the conditions are necessary)", necessity)
small = True
for n in (3, 4):
    verts = [(x, y) for x in range(n) for y in range(3)]
    for k in range(1, 3 * n + 1):
        for S in itertools.combinations(verts, k):
            if multiset_resolving(n, S) and not admissible(columns_of(n, set(S))):
                small = False
check("2c necessity of (A), (B), (C) for ALL resolving sets (any size) of P_3 [x] P_3 and P_3 [x] P_4", small)


# ---------------------------------------------------------------- 5 (first part). a(n) directly from the conditions
def count_weight_n(nmax):
    """a(n) = number of admissible sequences of length n and weight exactly n, 1 <= n <= nmax (no automaton)."""
    out = {}
    dp = {}
    for v in COLS:
        dp[(ZERO, v, wt[v])] = dp.get((ZERO, v, wt[v]), 0) + 1
    for n in range(1, nmax + 1):
        out[n] = sum(c for (u, v, w0), c in dp.items() if w0 == n and LEGAL[(u, v, ZERO)])
        nd = {}
        for (u, v, w0), c in dp.items():
            for w in COLS:
                if LEGAL[(u, v, w)]:
                    w1 = w0 + wt[w]
                    if w1 <= nmax:                 # weights never decrease; only weight <= nmax can matter
                        key = (v, w, w1)
                        nd[key] = nd.get(key, 0) + c
        dp = nd
    return out


NMAX = 200
a = count_weight_n(NMAX)


def brute_count(n):
    cnt = 0
    res = 0

    def rec(seq, weight):
        nonlocal cnt
        if weight > n:
            return
        if len(seq) == n:
            if weight == n and LEGAL[(seq[-2] if n >= 2 else ZERO, seq[-1], ZERO)]:
                cnt += 1
            return
        for w in COLS:
            if len(seq) >= 1:
                u = seq[-2] if len(seq) >= 2 else ZERO
                if not LEGAL[(u, seq[-1], w)]:
                    continue
            rec(seq + [w], weight + wt[w])

    rec([], 0)
    return cnt


check("5a a(5..8) = 188, 556, 1660, 5004 and a(10) = 44972 (paper), by the direct dynamic programme",
      [a[n] for n in (5, 6, 7, 8)] == [188, 556, 1660, 5004] and a[10] == 44972, str([a[n] for n in range(1, 11)]))
check("5b brute-force enumeration of all admissible sequences agrees for n = 1..8", all(brute_count(n) == a[n] for n in range(1, 9)))
# resolving among weight-8 admissible sequences


def all_weight_n(n):
    out = []

    def rec(seq, weight):
        if weight > n:
            return
        if len(seq) == n:
            if weight == n and LEGAL[(seq[-2], seq[-1], ZERO)]:
                out.append(tuple(seq))
            return
        for w in COLS:
            if len(seq) >= 1:
                u = seq[-2] if len(seq) >= 2 else ZERO
                if not LEGAL[(u, seq[-1], w)]:
                    continue
            rec(seq + [w], weight + wt[w])

    rec([], 0)
    return out


seq8 = all_weight_n(8)
res8 = sum(1 for seq in seq8 if multiset_resolving(8, {(x, y) for x in range(8) for y in range(3) if seq[x][y]}))
check("2d among the 5004 weight-8 admissible sequences exactly 1328 are multiset resolving", len(seq8) == 5004 and res8 == 1328,
      "%d of %d" % (res8, len(seq8)))
# automaton path counts (Theorem 16.1)
entry = [0] * n19
for t in pre_in:
    entry[idx[(t[2], t[3])]] += 1
final = [1 if s in accept else 0 for s in core]
vecl = entry[:]
auto = {}
for n in range(4, NMAX + 1):
    auto[n] = sum(x * y for x, y in zip(vecl, final))
    vecl = [sum(vecl[i] * Aadj[i][j] for i in range(n19)) for j in range(n19)]
check("5c automaton path counts (106 entries, tight transitions in L, 12 accepting states) equal a(n) for 5 <= n <= 200",
      all(auto[n] == a[n] for n in range(5, NMAX + 1)), "n = 4: automaton %d, direct %d" % (auto[4], a[4]))
rec_from = min(n0 for n0 in range(4, 60) if all(a[n] == 2 * a[n - 1] + a[n - 2] + 6 * a[n - 3] for n in range(n0, NMAX + 1)))
check("5d a(n) = 2 a(n-1) + a(n-2) + 6 a(n-3) for all n0 <= n <= 200 with n0 = %d (helper: n >= 6)" % rec_from, rec_from <= 6,
      "fails at n = %s" % [n for n in range(4, rec_from) if a[n] != 2 * a[n - 1] + a[n - 2] + 6 * a[n - 3]])
U = [0, 1]
Vl = [2, -1]
for n in range(2, 20001):
    U.append(-U[-1] - 2 * U[-2])
    Vl.append(-Vl[-1] - 2 * Vl[-2])
cf_from = min(n0 for n0 in range(1, 60) if all(42 * a[n] == 32 * 3 ** n - 9 * Vl[n] - 21 * U[n] for n in range(n0, NMAX + 1)))
check("5e 42 a(n) = 32*3^n - 9 V_n - 21 U_n for all n0 <= n <= 200 with n0 = %d (helper: n >= 3)" % cf_from, cf_from <= 3)
check("5f x^3 - 2x^2 - x - 6 = (x - 3)(x^2 + x + 2); a(n)/3^n -> 16/21",
      sympy.expand((sympy.Symbol("x") - 3) * (sympy.Symbol("x") ** 2 + sympy.Symbol("x") + 2)
                   - (sympy.Symbol("x") ** 3 - 2 * sympy.Symbol("x") ** 2 - sympy.Symbol("x") - 6)) == 0
      and abs(Fraction(a[200], 3 ** 200) - Fraction(16, 21)) < Fraction(1, 10 ** 30))

# ---------------------------------------------------------------- 3. column types
indep = sorted(t for t in COLS if not (t[0] and t[1]) and not (t[1] and t[2]))
missing = sorted((u, v) for u in FIVE for v in FIVE if (u, v) not in idx)
check("3a the five column types are the independent sets of the column path P_3", indep == FIVE)
check("3b L consists of 19 of the 25 ordered pairs of types", all(s[0] in FIVE and s[1] in FIVE for s in core) and len(missing) == 6,
      "the six absent pairs: %s" % [nm[u] + "," + nm[v] for u, v in missing])

# ---------------------------------------------------------------- 4. characteristic polynomial
x = sympy.Symbol("x")
cp = sympy.factor(M.charpoly(x).as_expr())
target = x ** 14 * (x - 3) * (x - 1) * (x + 1) * (x ** 2 + x + 2)
check("4a sympy: charpoly of the core = x^14 (x-3)(x-1)(x+1)(x^2+x+2)", sympy.expand(cp - target) == 0, str(cp))


def bareiss_det(mat):
    n = len(mat)
    Mx = [row[:] for row in mat]
    sign, prev = 1, 1
    for k in range(n - 1):
        if Mx[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if Mx[i][k] != 0), None)
            if sw is None:
                return 0
            Mx[k], Mx[sw] = Mx[sw], Mx[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                Mx[i][j] = (Mx[i][j] * Mx[k][k] - Mx[i][k] * Mx[k][j]) // prev
        prev = Mx[k][k]
    return sign * Mx[n - 1][n - 1]


pts = list(range(-9, 11))
vals = [bareiss_det([[(t if i == j else 0) - Aadj[i][j] for j in range(n19)] for i in range(n19)]) for t in pts]
tvals = [int(target.subs(x, t)) for t in pts]
check("4b Bareiss determinants det(tI - A) at t = -9..10 (20 points determine a monic polynomial of degree 19) equal the target",
      vals == tvals)


def matmul(P, Q):
    return [[sum(P[i][k] * Q[k][j] for k in range(n19)) for j in range(n19)] for i in range(n19)]


P = [[int(i == j) for j in range(n19)] for i in range(n19)]
tr_ok = True
for j in range(1, 31):
    P = matmul(P, Aadj)
    tr_ok &= (sum(P[i][i] for i in range(n19)) == 3 ** j + 1 + (-1) ** j + Vl[j])
check("4c tr(A^j) = 3^j + 1 + (-1)^j + V_j for j = 1..30 (third independent confirmation of the nonzero spectrum)", tr_ok)
rank = M.rank()
check("4d rank of A is %d: the eigenvalue 0 has algebraic multiplicity 14 and geometric multiplicity %d" % (rank, n19 - rank), rank < n19)

# ---------------------------------------------------------------- 6. Lucas sequence and Ramanujan-Nagell
check("6a V_n^2 + 7 U_n^2 = 2^(n+2) for n <= 2000 (norm of ((-1 + sqrt(-7))/2)^n is 2^n)",
      all(Vl[n] ** 2 + 7 * U[n] ** 2 == 2 ** (n + 2) for n in range(0, 2001)))
ones = [n for n in range(1, 20001) if abs(U[n]) == 1]
check("6b |U_n| = 1 for n <= 20000 exactly at n = 1, 2, 3, 5, 13", ones == [1, 2, 3, 5, 13], str(ones))
check("6c the corresponding |V_n| = 1, 3, 5, 11, 181 solve x^2 + 7 = 2^(n+2) with exponents 3, 4, 5, 7, 15 (Ramanujan-Nagell list)",
      [abs(Vl[n]) for n in ones] == [1, 3, 5, 11, 181] and all(Vl[n] ** 2 + 7 == 2 ** (n + 2) for n in ones))
check("6d the indices 1, 2, 3, 5, 13 contain 2, 3, 5, 13 of the owner's list but not 7, 19, 27; discriminant of x^2 + x + 2 is -7",
      1 - 8 == -7 and set(ones) & {2, 3, 5, 7, 13, 19, 27, 1039} == {2, 3, 5, 13})

if all(RESULTS):
    print("ALL CHECKS PASSED (%d checks)" % len(RESULTS))
    sys.exit(0)
print("CHECK FAILED (%d of %d)" % (RESULTS.count(False), len(RESULTS)))
sys.exit(1)
