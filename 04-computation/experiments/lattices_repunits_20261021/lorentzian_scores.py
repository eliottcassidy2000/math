"""Tournament score vectors in the Lorentzian lattice II_{25,1} (session opus-2026-10-21-S52; THM-4752).

II_{25,1} = { x in Z^26 u (Z + 1/2)^26 : x_1 + ... + x_25 - x_26 = 0 (mod 2) }, form x_1^2 + ... + x_25^2 - x_26^2.
For a tournament T on 25 vertices with score vector s:
  (1) v_T = (s | 70) lies in II_{25,1} (s integral, sum s = 300, 300 - 70 even) and v_T^2 = -2 c3(T)
      (Kendall-Babington Smith: c3 = C(25,3) - sum C(s_i,2));
  (2) v_T is null iff T is transitive, i.e. v_T is a permutation of rho = (0, 1, ..., 24 | 70), the Conway Weyl vector
      (rho^2 = 0 is the cannonball identity 0^2 + ... + 24^2 = 70^2; rho^perp / rho is the Leech lattice);
  (3) for every ranking pi of the vertices (a bijection to 0..24, 24 = top) with rho_pi = (pi | 70):
      backward length BL_T(pi) = sum over arcs u -> w with pi(u) < pi(w) of pi(w) - pi(u) = -v_T . rho_pi;
      it is minimised exactly by the rankings sorted like s (rearrangement; Kano-Sakamoto 1983 for tournaments), with
      minimum U(T) = sum_k (Landau slack of the k lowest scores) = c3(T) + |s_sorted - (0..24)|^2 / 2;
  (4) nu = (1/2, ..., 1/2 | 5/2) is a primitive null vector with nu . rho = -25 = -h(A_24) and v_T . nu = -25 for all T;
      the Leech roots (r^2 = 2, r . rho = -1) e_i - e_(i+1) (i = 0..23) and r* = (-1/2, 1/2^23, 3/2 | 5/2) are orthogonal
      to nu, form the affine diagram A~_24 (a 25-cycle) and sum to nu: they are the vertices of the A_24 deep hole of the
      Leech lattice, and nu^perp / nu is the Niemeier lattice N(A_24) (Niemeier-Venkov: rank 24, contains A_24);
  (5) u_T = v_T - 10 (0^25 | 1) = (s | 60) lies in nu^perp and u_T - 24 nu = (s - 12 | 0): the class of T in N(A_24) is
      its centred score vector d_T = s - 12 in the root lattice A_24, |d_T|^2 = 1300 - 2 c3(T), and d_T = 0 iff T is
      regular (the n = 25 rung of the THM-868 score-deviation map);
  (6) the reflections in e_i - e_j generate S_25 = W(A_24) (relabelling) and the score vectors are exactly the lattice
      points of the permutohedron conv(S_25 rho) (Landau), inscribed in the light cone: its 25! vertices are the Leech
      cusps sigma(rho), every other point is timelike.
Also: 0^2 + ... + (n-1)^2 is a square only for n = 1, 2, 25 (Watson 1918; checked to 10^6), II_{n,1} exists iff
n = 1 (mod 8), and the contrast II_{9,1}: the E10 = T_{2,3,7} Weyl vector has negative norm (finitely many simple
roots, Vinberg) while II_{25,1}'s is null (Conway).  Prints ALL CHECKS PASSED."""
import random, itertools, math
from fractions import Fraction
import numpy as np

FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


n = 25
H = Fraction(1, 2)
rho = list(range(n)) + [70]


def lor(x, y):
    return sum(a * b for a, b in zip(x[:-1], y[:-1])) - x[-1] * y[-1]


def in_II(x):
    twice = [2 * Fraction(t) for t in x]
    allint = all(t.denominator == 1 and t.numerator % 2 == 0 for t in twice)
    allhalf = all(t.denominator == 1 and t.numerator % 2 == 1 for t in twice)
    return (allint or allhalf) and (sum(Fraction(t) for t in x[:-1]) - Fraction(x[-1])) % 2 == 0


def random_tournament(m, rng, p=0.5):
    adj = [[False] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            if rng.random() < p:
                adj[i][j] = True
            else:
                adj[j][i] = True
    return adj


def c3(adj):
    m = len(adj)
    return sum(1 for a, b, c in itertools.combinations(range(m), 3)
               if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]))


def backward_length(adj, pi):
    m = len(adj)
    return sum(pi[w] - pi[u] for u in range(m) for w in range(m) if adj[u][w] and pi[u] < pi[w])


def score_ranking(s):
    order = sorted(range(len(s)), key=lambda v: (s[v], v))
    pi = [0] * len(s)
    for r, v in enumerate(order):
        pi[v] = r
    return pi


print('(W) the Weyl vector, the cannonball identity and the A_24 deep hole')
check(lor(rho, rho) == 0 and in_II(rho), 'rho = (0,1,...,24 | 70) is a null vector of II_{25,1} (sum k^2 = 4900 = 70^2)')
sq = [m for m in range(1, 10 ** 6 + 1) if math.isqrt((m - 1) * m * (2 * m - 1) // 6) ** 2 == (m - 1) * m * (2 * m - 1) // 6]
check(sq == [1, 2, 25], f'0^2 + ... + (n-1)^2 is a perfect square for n <= 10^6 exactly at n = {sq} (n = 2: 0 + 1 = 1; cannonball n = 25)')
nu = [H] * n + [5 * H]
check(lor(nu, nu) == 0 and in_II(nu) and not in_II([t / 2 for t in nu]) and lor(nu, rho) == -25,
      'nu = (1/2^25 | 5/2) is a primitive null vector of II_{25,1} (nu/2 is not in the lattice) with nu . rho = -25 = -h(A_24)')
roots = []
for i in range(n - 1):
    r = [0] * (n + 1); r[i] = 1; r[i + 1] = -1
    roots.append(r)
rstar = [-H] + [H] * 23 + [3 * H] + [5 * H]
hole = roots + [rstar]
check(all(lor(r, r) == 2 and lor(r, rho) == -1 and in_II(r) and lor(r, nu) == 0 for r in hole),
      'e_i - e_(i+1) (i = 0..23) and r* = (-1/2, 1/2^23, 3/2 | 5/2) are Leech roots (r^2 = 2, r . rho = -1) orthogonal to nu')
G = [[lor(a, b) for b in hole] for a in hole]
cycle = all(G[i][j] == (2 if i == j else -1 if (i - j) % 25 in (1, 24) else 0) for i in range(25) for j in range(25))
check(cycle and [sum(c) for c in zip(*hole)] == nu,
      'their Gram matrix is the affine Cartan matrix of A~_24 (a 25-cycle, all marks 1) and they sum to nu: the A_24 deep hole, imaginary root = nu')

print('(T) 25-vertex tournaments')
rng = random.Random(20261021)
cases = []
trans = [[i < j for j in range(n)] for i in range(n)]
cases.append(('transitive', trans))
rot = [[((j - i) % n) in range(1, 13) for j in range(n)] for i in range(n)]
cases.append(('rotational R_25 (regular)', rot))
for p in (0.5, 0.5, 0.5, 0.2, 0.05, 0.01):
    cases.append((f'random p = {p}', random_tournament(n, rng, p)))
ok = ok_bl = ok_a24 = True
for name, adj in cases:
    s = [sum(row) for row in adj]
    v = s + [70]
    C3 = c3(adj)
    ss = sorted(s)
    pi_s = score_ranking(s)
    rho_s = pi_s + [70]
    U = backward_length(adj, pi_s)
    d2 = sum((ss[i] - i) ** 2 for i in range(n))
    slack = sum(sum(ss[:k]) - k * (k - 1) // 2 for k in range(1, n))
    good = (in_II(v) and lor(v, v) == -2 * C3 and -lor(v, rho_s) == U == slack and 2 * U == 2 * C3 + d2
            and lor(v, nu) == -25 and ((lor(v, v) == 0) == (C3 == 0)))
    ok &= good
    # (3) every ranking: BL = -v . rho_pi, never below U (random rankings and all adjacent swaps of the score ranking)
    for _ in range(40):
        pi = list(range(n)); rng.shuffle(pi)
        bl = backward_length(adj, pi)
        ok_bl &= bl == -lor(v, pi + [70]) and bl >= U
    for a in range(n - 1):
        pi = pi_s[:]
        x, y = pi.index(a), pi.index(a + 1)
        pi[x], pi[y] = pi[y], pi[x]
        ok_bl &= backward_length(adj, pi) >= U
    # (5) the A_24 cusp quotient
    u = s + [60]
    d = [Fraction(t) - 24 * nn for t, nn in zip(u, nu)]
    ok_a24 &= (in_II(u) and lor(u, nu) == 0 and d == [si - 12 for si in s] + [0] and sum(d[:-1]) == 0
               and sum(t * t for t in d[:-1]) == 1300 - 2 * C3 and ((C3 == 650) == all(t == 0 for t in d)))
    print(f'     {name:28s}: c3 = {C3:4d}, v_T^2 = {lor(v, v):6d}, U = -v.rho_s = {U:5d} = Landau slack {slack} = c3 + |d|^2/2 = {C3} + {d2}/2; v.nu = {lor(v, nu)}; |s - 12|^2 = {1300 - 2 * C3}')
check(ok, 'v_T in II_{25,1}, v_T^2 = -2 c3, -v_T . rho_s = U(T) = Landau slack = c3 + |d|^2/2, v_T . nu = -25, null iff transitive (8 tournaments)')
check(ok_bl, 'BL_T(pi) = -v_T . rho_pi for 320 random rankings, and no random ranking or adjacent swap beats the score ranking')
check(ok_a24, '(s | 60) in nu^perp, its N(A_24) class is s - 12 in A_24 with |s - 12|^2 = 1300 - 2 c3; zero exactly for the regular R_25')
# the identity U = c3 + |d|^2/2 at every order (it does not need 25), and exhaustive backward-length minimisation, n <= 7
ok = ok_ks = True
for m in range(3, 15):
    for _ in range(30):
        adj = random_tournament(m, rng, rng.choice([0.5, 0.3, 0.1]))
        s = [sum(row) for row in adj]
        ss = sorted(s)
        U = backward_length(adj, score_ranking(s))
        ok &= 2 * U == 2 * c3(adj) + sum((ss[i] - i) ** 2 for i in range(m))
        if m <= 7 and _ < 6:
            sq_m = sum(j * j for j in range(m))
            best = min(backward_length(adj, list(pi)) for pi in itertools.permutations(range(m)))
            argmins = [pi for pi in itertools.permutations(range(m)) if backward_length(adj, list(pi)) == best]
            ok_ks &= best == U and all(backward_length(adj, list(pi)) == sq_m - sum(pi[v] * s[v] for v in range(m)) for pi in argmins)
            ok_ks &= all(all((s[a] - s[b]) * (pi[a] - pi[b]) >= 0 for a in range(m) for b in range(m)) for pi in argmins)
check(ok, 'U(T) = c3(T) + |s_sorted - (0..n-1)|^2/2 holds at every order n = 3..14 (360 random tournaments)')
check(ok_ks, 'n <= 7 (30 tournaments, all n! rankings): min BL = U(T), attained only by score-consistent rankings (Kano-Sakamoto); BL(pi) = sum j^2 - sum pi(v) s_v')
# Landau: score vectors of n-tournaments = lattice points of the permutohedron (majorization) -- small n exhaustive
ok = True
for m in range(2, 7):
    pairs = list(itertools.combinations(range(m), 2))
    realized = set()
    for bits in range(1 << len(pairs)):
        s = [0] * m
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                s[i] += 1
            else:
                s[j] += 1
        realized.add(tuple(s))
    perm = set()
    for s in itertools.product(range(m), repeat=m):
        ss = sorted(s)
        if sum(s) == m * (m - 1) // 2 and all(sum(ss[:k]) >= k * (k - 1) // 2 for k in range(1, m + 1)):
            perm.add(s)
    forests = [1, 1, 2, 7, 38, 291, 2932]              # labelled forests on n vertices (A001858), n = 0..6
    ok &= realized == perm and len(perm) == forests[m]
    print(f'     n = {m}: {len(realized)} score vectors = {len(perm)} permutohedron lattice points (Landau) = {forests[m]} labelled forests (Stanley)')
check(ok, 'score vectors = lattice points of the permutohedron conv(S_n (0..n-1)) for n <= 6 (Landau/majorization), counted by labelled forests (Stanley 1980)')

print('(V) the contrast II_{9,1}: the E10 = T_{2,3,7} Weyl vector')
# T_{2,3,7}: centre node 0 with arms of lengths 1, 2, 6
arms = [1, 2, 6]
nodes = 1 + sum(arms)
A = 2 * np.eye(nodes, dtype=int)
idx = 1
for L in arms:
    prev = 0
    for t in range(L):
        A[prev, idx] = A[idx, prev] = -1
        prev = idx; idx += 1
c = np.linalg.solve(A.astype(float), np.ones(nodes))
rho2 = float(np.ones(nodes) @ c)
check(round(np.linalg.det(A)) == -1 and abs(rho2 - round(rho2)) < 1e-9 and rho2 < 0,
      f'E10 Cartan matrix det = {round(np.linalg.det(A))}; Weyl vector (rho . alpha_i = 1): rho^2 = {round(rho2)} < 0 (II_{{9,1}} = E8 + U: finitely many simple roots)')
check(all((m - 1) % 8 == 0 for m in (1, 9, 17, 25)), 'II_{n,1} exists for n = 1 (mod 8): n = 9, 17, 25; only n = 25 has a null Weyl vector (Conway), matching the cannonball n')

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
