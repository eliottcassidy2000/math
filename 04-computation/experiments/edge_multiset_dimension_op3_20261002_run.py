#!/usr/bin/env python3
"""edge_multiset_dimension_op3_20261002_run.py -- Allikvere's Open Problem 3: one estimate from d = 10.

Companion runner of 05-knowledge/results/edge_multiset_dimension_op3_20261002.md.  Pure Python 3, exact integer and
rational arithmetic throughout (square roots are bounded from above by integer square roots).

Setting: a uniformly random subset S of V(Q_d) (density 1/2).  For a pair type t of edge pairs, its level graph has
levels 0..T and edge weights N_ab (the number of vertices s with {d(e,s), d(f,s)} = {a, b}).  The forest lemma bounds
P(H_e = H_f) by the product of beta(N) over the edges of any forest, beta(N) = C(N, floor(N/2)) / 2^N.  The potential
forest of t lets every level pick its heaviest edge to a strictly more central level.  V_d is the resulting union
bound with antipodal halving.  The note proves V_(d+1) <= rho_d V_d + A_(d+1) (transfer lemma); this runner checks

  S1  the closed-form cells and type counts against brute force over representative edge pairs (d = 3..8);
  S2  the union bound itself: potential forests against maximum-weight (Kruskal) forests, d = 6..12, and the exact
      value V_10 < 0.6603 (the base of the induction);
  S3  the transfer lemma at the level of sub-cells, for every type and 3 <= d <= 40: the explicit sub-cell maps are
      injective, land on the image edge and do not decrease sizes; targets stay strictly more central; the explicit
      extra chooser is not an image and its edge weight equals the closed form N*(t');
  S4  rho_d <= 1/2 for 10 <= d <= 40 (rational upper bounds) and the crude bound for d >= 41;
  S5  the antipodal terms A_d, d >= 11, and the conclusion V_d <= 0.677 for every d >= 10.
usage: python3 edge_multiset_dimension_op3_20261002_run.py [DMAX_TRANSFER=40]
"""
import sys
import time
from math import comb, isqrt, log, log2, pi, sqrt
from fractions import Fraction as Fr

T0 = time.time()
CHECKS = 0


def ok(cond, msg):
    global CHECKS
    CHECKS += 1
    if not cond:
        print('FAILED:', msg)
        sys.exit(1)


def section(title):
    print()
    print('=' * 78)
    print(title)
    print('=' * 78)
    sys.stderr.write('[%8.1fs] %s\n' % (time.time() - T0, title.split(':')[0]))


# ----------------------------------------------------------------------------------------------------------------
# types, sub-cells, weights
# Parallel type ('P', h, nu): h + nu = d - 1, 1 <= h.  Sub-cell (c, j): 2 C(nu,c) C(h,j) vertices, levels
#   (c + j) -> (c + h - j).   Crossing type ('X', h, nu): h + nu = d - 2, 0 <= h.  Sub-cell (x, y, c, j):
#   C(nu,c) C(h,j) vertices, levels (x + c + j) -> (y + c + h - j).
# ----------------------------------------------------------------------------------------------------------------
def top_of(t):
    k, h, nu = t
    return h + nu if k == 'P' else h + nu + 1


def subcells(t):
    """yields (frm, to, size, label)"""
    k, h, nu = t
    if k == 'P':
        for c in range(nu + 1):
            for j in range(h + 1):
                yield c + j, c + h - j, 2 * comb(nu, c) * comb(h, j), ('P', c, j)
    else:
        for x in (0, 1):
            for y in (0, 1):
                for c in range(nu + 1):
                    for j in range(h + 1):
                        yield x + c + j, y + c + h - j, comb(nu, c) * comb(h, j), ('X', x, y, c, j)


def size_of(t, label):
    k, h, nu = t
    if label[0] == 'P':
        _, c, j = label
        return 2 * comb(nu, c) * comb(h, j) if (0 <= c <= nu and 0 <= j <= h) else 0
    _, x, y, c, j = label
    return comb(nu, c) * comb(h, j) if (0 <= c <= nu and 0 <= j <= h) else 0


def ends_of(t, label):
    k, h, nu = t
    if label[0] == 'P':
        _, c, j = label
        return c + j, c + h - j
    _, x, y, c, j = label
    return x + c + j, y + c + h - j


_W = {}


def weights(t):
    if t not in _W:
        W = {}
        for a, b, s, _ in subcells(t):
            if a != b:
                key = (min(a, b), max(a, b))
                W[key] = W.get(key, 0) + s
        _W[t] = W
    return _W[t]


def types(d):
    return [('P', h, d - 1 - h) for h in range(1, d)] + [('X', h, d - 2 - h) for h in range(0, d - 1)]


def count(d, t):
    k, h, nu = t
    return d * 2 ** (d - 2) * comb(d - 1, h) if k == 'P' else d * (d - 1) * 2 ** (d - 1) * comb(d - 2, h)


def wt(t):
    """antipodal halving: the antipodal type ('P', d-1, 0) has weight 1, every other type 1/2"""
    k, h, nu = t
    return Fr(1) if (k == 'P' and nu == 0) else Fr(1, 2)


# ----------------------------------------------------------------------------------------------------------------
# centrality and potential forests
# b is strictly more central than a (top T) iff |2b - T| < |2a - T|, or the two are mirror images and b is the upper.
# ----------------------------------------------------------------------------------------------------------------
def more_central(b, a, T):
    db, da = abs(2 * b - T), abs(2 * a - T)
    return db < da or (db == da and b > a)


_POT = {}


def potential(t):
    """returns (T, adj, choice) with choice[a] = (b, N_ab) for every level a having a strictly more central neighbour"""
    if t not in _POT:
        T = top_of(t)
        adj = {a: {} for a in range(T + 1)}
        for (a, b), w in weights(t).items():
            adj[a][b] = w
            adj[b][a] = w
        ch = {}
        for a in range(T + 1):
            best = None
            for b, w in adj[a].items():
                if more_central(b, a, T) and (best is None or w > best[1]):
                    best = (b, w)
            if best is not None:
                ch[a] = best
        _POT[t] = (T, adj, ch)
    return _POT[t]


_BETA = {}


def beta(N):
    if N not in _BETA:
        _BETA[N] = Fr(comb(N, N // 2), 2 ** N)
    return _BETA[N]


def beta_upper(N, exact_up_to=4000):
    """a rational upper bound for beta(N): exact for small N, else sqrt(2/(pi N)) < sqrt(2/(3.14159 N))"""
    if N <= exact_up_to:
        return beta(N)
    K = 10 ** 15
    x = Fr(2 * 100000, 314159 * N)
    return Fr(isqrt(x.numerator * K * K // x.denominator) + 1, K)


def is_forest(edges, n):
    par = list(range(n))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for a, b in edges:
        ra, rb = f(a), f(b)
        if ra == rb:
            return False
        par[ra] = rb
    return True


def B_pot(t):
    T, adj, ch = potential(t)
    B = Fr(1)
    for b, w in ch.values():
        B *= beta(w)
    return B


def V_exact(d):
    return sum((wt(t) * count(d, t) * B_pot(t) for t in types(d)), Fr(0))


def kruskal_bound(d, halving):
    tot = Fr(0)
    for t in types(d):
        W = weights(t)
        T = top_of(t)
        par = list(range(T + 1))

        def f(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x
        B = Fr(1)
        for w, (a, b) in sorted(((w, e) for e, w in W.items()), reverse=True):
            ra, rb = f(a), f(b)
            if ra != rb:
                par[ra] = rb
                B *= beta(w)
        tot += (wt(t) if halving else 1) * count(d, t) * B
    return tot


# ----------------------------------------------------------------------------------------------------------------
# the transfer: predecessor rule, level map, sub-cell maps, extra chooser
# ----------------------------------------------------------------------------------------------------------------
def pred(t2):
    """predecessor of a non-antipodal type t2 at d+1, and the step used"""
    k, h2, nu2 = t2
    if k == 'P':
        return 'nu', ('P', h2, nu2 - 1)
    if h2 >= nu2 and h2 >= 1:
        return 'h', ('X', h2 - 1, nu2)
    return 'nu', ('X', h2, nu2 - 1)


def succ_type(t, step):
    k, h, nu = t
    return (k, h, nu + 1) if step == 'nu' else (k, h + 1, nu)


def phi(a, T):
    return a if 2 * a < T else a + 1


def psi(b, a, T, step):
    lower = 2 * a < T
    if step == 'nu':
        return b if lower else b + 1
    return b + 1 if lower else b


def map_label(label, lower, from_a, step):
    """the sub-cell map of the transfer lemma; from_a: the sub-cell is oriented from the chooser a"""
    if label[0] == 'P':
        _, c, j = label
        if step == 'nu':
            return ('P', c, j) if lower else ('P', c + 1, j)
        bump = (not from_a) if lower else from_a
        return ('P', c, j + 1) if bump else ('P', c, j)
    _, x, y, c, j = label
    if step == 'nu':
        return ('X', x, y, c, j) if lower else ('X', x, y, c + 1, j)
    bump = (not from_a) if lower else from_a
    return ('X', x, y, c, j + 1) if bump else ('X', x, y, c, j)


def extra_chooser(t, step):
    """(a*, b*) at d+1: a* is not the image of a chooser of t, b* is strictly more central (cases A1-B2, X)"""
    k, h, nu = t
    T = top_of(t)
    if k == 'X' or h % 2 == 1:
        return (T // 2, T // 2 + 1) if T % 2 == 0 else ((T + 3) // 2, (T + 1) // 2)
    return (T // 2 + 2, T // 2) if T % 2 == 0 else ((T - 1) // 2, (T + 3) // 2)


def N_star(t2):
    k, h2, nu2 = t2
    if k == 'P':
        return 4 * comb(nu2, nu2 // 2) * comb(h2, (h2 + 1) // 2 - 1)
    if h2 % 2 == 1:
        return 2 * comb(nu2 + 1, (nu2 + 1) // 2) * comb(h2, (h2 - 1) // 2)
    return 2 * comb(nu2, nu2 // 2) * comb(h2 + 1, h2 // 2)


def transfer_check(d):
    """S3 for one d; returns the number of chooser edges checked"""
    nedges = 0
    for t2 in types(d + 1):
        if t2[0] == 'P' and t2[2] == 0:
            continue
        step, t = pred(t2)
        ok(t in types(d) and succ_type(t, step) == t2, ('predecessor', d, t2))
        T, adj, ch = potential(t)
        T2, adj2, ch2 = potential(t2)
        ok(T2 == T + 1, 'top')
        images = set()
        for a, (b, w) in ch.items():
            ok(2 * a != T, ('centre chose', t))
            lower = 2 * a < T
            a2, b2 = phi(a, T), psi(b, a, T, step)
            ok(more_central(b2, a2, T2), ('target not more central', t2, a, b))
            seen = set()
            tot_old = 0
            for frm, to, s, lab in subcells(t):
                if {frm, to} != {a, b} or frm == to or s == 0:
                    continue
                tot_old += s
                lab2 = map_label(lab, lower, frm == a, step)
                f2, g2 = ends_of(t2, lab2)
                want = (a2, b2) if frm == a else (b2, a2)
                ok((f2, g2) == want, ('image sub-cell on the wrong edge', t2, lab, lab2))
                ok(size_of(t2, lab2) >= s, ('image sub-cell smaller', t2, lab, lab2))
                ok(lab2 not in seen, ('sub-cell map not injective', t2, lab2))
                seen.add(lab2)
            ok(tot_old == w, ('edge weight', t, a, b))
            ok(adj2[a2].get(b2, 0) >= w, ('image edge lighter', t2))
            ok(a2 not in images, 'level map not injective')
            images.add(a2)
            nedges += 1
        a_s, b_s = extra_chooser(t, step)
        ok(a_s not in images, ('extra chooser is an image', t2))
        ok(more_central(b_s, a_s, T2), ('extra target not more central', t2))
        ok(adj2[a_s].get(b_s, 0) == N_star(t2), ('extra edge weight', t2, adj2[a_s].get(b_s, 0), N_star(t2)))
        # consequence: B(t2) <= B(t) beta(N*)
        if d <= 13:
            ok(B_pot(t2) <= B_pot(t) * beta(N_star(t2)), ('B(t2) <= B(t) beta(N*)', t2))
    return nedges


def rho_upper(d):
    """rational upper bound for rho_d and the predecessor attaining it"""
    load = {}
    for t2 in types(d + 1):
        if t2[0] == 'P' and t2[2] == 0:
            continue
        step, t = pred(t2)
        R = Fr(wt(t2) * count(d + 1, t2), wt(t) * count(d, t))
        load[t] = load.get(t, Fr(0)) + R * beta_upper(N_star(t2))
    t = max(load, key=load.get)
    return load[t], t


def A_upper(d):
    """rational upper bound for the antipodal term A_d = d 2^(d-2) prod_{j < (d-1)/2} beta(4 C(d-1, j))"""
    L = d - 1
    B = Fr(d * 2 ** (d - 2))
    for j in range((L + 1) // 2):
        B *= beta_upper(4 * comb(L, j))
    return B


def main():
    dmax_tr = int(sys.argv[1]) if len(sys.argv) > 1 else 40

    section('S1: cells and type counts against brute force')
    for d in range(3, 9):
        E = d * 2 ** (d - 1)
        tot = 0
        for t in types(d):
            k, h, nu = t
            v = (1 << h) - 1
            e = (0, 1 << (d - 1))
            f = (v, v | (1 << (d - 1))) if k == 'P' else (v, v | (1 << (d - 2)))
            W = {}
            for s in range(1 << d):
                da = min(bin(e[0] ^ s).count('1'), bin(e[1] ^ s).count('1'))
                db = min(bin(f[0] ^ s).count('1'), bin(f[1] ^ s).count('1'))
                if da != db:
                    key = (min(da, db), max(da, db))
                    W[key] = W.get(key, 0) + 1
            ok(W == weights(t), ('weights', d, t))
            tot += count(d, t)
        ok(tot == comb(E, 2), ('type counts', d))
        print('d=%d: %d types; closed-form edge weights match brute force; counts add up to C(E,2) = %d'
              % (d, len(types(d)), comb(E, 2)))
    for d in range(9, 41):
        ok(sum(count(d, t) for t in types(d)) == comb(d * 2 ** (d - 1), 2), ('type counts', d))
    print('type counts add up to C(E,2) for 9 <= d <= 40')

    section('S2: the union bound V_d (potential forests, antipodal halving) and the base value V_10')
    for d in range(6, 13):
        for t in types(d):
            T, adj, ch = potential(t)
            ok(is_forest([(a, b) for a, (b, w) in ch.items()], T + 1), ('potential forest', d, t))
    print('potential forests are forests for every type, 6 <= d <= 12')
    for d in range(6, 13):
        V = V_exact(d)
        K = kruskal_bound(d, True)
        K0 = kruskal_bound(d, False)
        ok(V >= K, ('maximum-weight forests are optimal', d))
        print('d=%2d: V_d = %.6f; maximum-weight forests with halving %.6f%s; without halving %.6f'
              % (d, float(V), float(K), ' (equal)' if V == K else '', float(K0)))
        if d == 10:
            V10 = V
    ok(V10 < Fr(6603, 10000), 'V_10 < 0.6603')
    print('V_10 = %d/%d-bit fraction < 0.6603 (exact)' % (V10.numerator.bit_length(), V10.denominator.bit_length()))

    section('S3: the transfer lemma at the level of sub-cells, every type, 3 <= d <= %d' % dmax_tr)
    tot = 0
    for d in range(3, dmax_tr + 1):
        tot += transfer_check(d)
    print('%d chooser edges transferred: sub-cell maps injective, on the image edge, sizes not decreasing;' % tot)
    print('targets strictly more central; extra chooser not an image, extra edge weight = N*(t\') in every case;')
    print('B(t\') <= B(t) beta(N*(t\')) checked exactly for d <= 13')

    section('S4: rho_d <= 1/2 for every d >= 10')
    worst = Fr(0)
    for d in range(10, 41):
        r, t = rho_upper(d)
        ok(r <= Fr(1, 2), ('rho', d))
        worst = max(worst, r)
        if d <= 16 or d % 8 == 0:
            print('d=%2d -> %2d: rho_d <= %.5f (largest load on %s)' % (d, d + 1, float(r), t))
    print('rho_d <= %.5f for 10 <= d <= 40 (exact rational upper bounds)' % float(worst))
    for d in (6, 7, 8, 9):
        r, t = rho_upper(d)
        print('  (for comparison, d=%d: rho_d <= %.4f)' % (d, float(r)))
    # crude bound for d >= 41 (proved in the note); here: its ingredients and its value at d = 41
    for m in range(0, 301):
        ok(comb(m, m // 2) * (m + 1) >= 2 ** m, ('central binomial', m))
        if m >= 1:
            ok(comb(m, (m + 1) // 2 - 1) * 2 * (m + 1) >= 2 ** m, ('near-central binomial', m))
    fP = lambda d: 2 * (d + 1) ** 1.5 * sqrt(2 / pi) * 2 ** (-d / 2)
    fX = lambda d: 4 * (d + 1) ** 2 / (d - 1) * sqrt(2 / pi) * 2 ** (-d / 2)
    ok(max(fP(41), fX(41)) < 0.001, 'crude bound at 41')
    for d in range(41, 65):
        r, t = rho_upper(d)
        ok(float(r) <= max(fP(d), fX(d)) * 1.000001, ('crude bound dominates', d))
    print('crude bound max(2(d+1)^(3/2), 4(d+1)^2/(d-1)) sqrt(2/pi) 2^(-d/2): %.2e at d = 41, decreasing;' % max(fP(41), fX(41)))
    print('it dominates the exact rho_d bound for 41 <= d <= 64; binomial lower bounds checked for m <= 300')

    section('S5: antipodal terms and the conclusion')
    SA = Fr(0)
    for d in range(11, 61):
        A = A_upper(d)
        SA += A
        if d <= 14:
            print('A_%d <= %.3e' % (d, float(A)))
    # tail: A_d <= d 2^(d-2) (3/8) (2 pi (d-1))^(-m/2), m = ceil((d-1)/2) - 1 >= (d-3)/2, and that is <= 2^(-d/2), d >= 61
    for d in range(61, 400):
        lg = log2(d) + (d - 2) + log2(3 / 8) - ((d - 3) / 4) * log2(2 * pi * (d - 1))
        ok(lg <= -d / 2, ('antipodal tail', d))
    tail = 2 ** (-30.5) / (1 - 2 ** -0.5)
    total = float(SA) + tail
    print('sum of A_d over 11 <= d <= 60 <= %.6f; tail d >= 61 <= %.1e (A_d <= 2^(-d/2): proved in the note,' % (float(SA), tail))
    print('  checked numerically for 61 <= d < 400)')
    bound = float(V10) + total
    ok(bound < 0.677, 'final')
    print('V_(d+1) <= rho_d V_d + A_(d+1) with rho_d <= 1/2 gives, for every d >= 10,')
    print('  V_d <= V_10 + sum_{k >= 11} A_k <= %.4f < 1' % bound)
    V11 = V_exact(11)
    r10, _ = rho_upper(10)
    ok(V11 <= r10 * V10 + A_upper(11), 'recursion at d = 10')
    print('consistency: V_11 = %.5f <= rho_10 V_10 + A_11 = %.5f' % (float(V11), float(r10 * V10 + A_upper(11))))

    print()
    print('%d checks' % CHECKS)
    print('ALL CHECKS PASSED')
    sys.stderr.write('[%8.1fs] done\n' % (time.time() - T0))


if __name__ == '__main__':
    main()
