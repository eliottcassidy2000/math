#!/usr/bin/env python3
"""Literal native charges versus a shell decoder; integer arithmetic only."""
from collections import defaultdict
from itertools import product
import json

W = ((1, 0), (0, 1), (-1, -1))
U = ((0, 0), (1, 0), (1, 1))
CHECKS = 0


def need(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def add(*zs):
    return tuple(map(sum, zip(*zs)))


def scale(z, c):
    return tuple(c*x for x in z)


def mul(z, w):
    a, b = z; c, d = w
    return a*c-b*d, a*d+b*c-b*d


def norm(z):
    a, b = z
    return a*a-a*b+b*b


def primitive(s):
    return s > 0 and s % 2 == 1 and s % 8 != 5


def native(s):
    a = 2 if s % 8 == 1 else 1
    y, t = (3*s+1) >> a, a
    while y % 8 == 5:
        y = (y-1)//4
        t += 2
    return y, t, a


def ray(u, a):
    D = 2**a+3
    tau = next(t for t in range(a, 7, 2) if (2**t*(3*u+1)-D) % 9 == 0)
    return 2**tau*(3*u+1), D, tau


def source(C, D, k):
    return (C*64**k-D)//9


def literal_charge(s, N):
    # Native B-chain, without the shell formula or cancellation table.
    lower = ((1, 2), (-2, -1), (1, -1))
    upper = ((0, 1), (-1, -1), (1, 0))
    z = (0, 0)
    while s < 4*N:
        z = add(z, (lower if s < 2*N else upper)[s % 3])
        s = 4*s+1
    return z


def geometric_sum(h, eta):
    if h < 0:
        return 0, 0
    seq = ((1, 0), (0, -1), (0, 0)) if eta == -1 else ((1, 0), (1, 1), (0, 0))
    return seq[h % 3]


def shell_point(C, E, i, k):
    return (C*4**i-E*4**(i-3*k)-3)//9


def prefix(X, C, D, K=None):
    """Sum omega^m over this ray's B-chain terms below X, with k <= K."""
    if (C-D)//9 >= X or K is not None and K < 0:
        return 0, 0
    E = D-3
    eta = -1 if D == 5 else 1
    i = 0
    while shell_point(C, E, i+1, 0) < X:
        i += 1
    ell, r = divmod(i, 3)
    cap = ell if K is None else min(ell, K)
    h = next(k for k in range(cap, -1, -1) if shell_point(C, E, i, k) < X)
    P = mul(geometric_sum(cap, eta), U[r])
    partial = mul(W[r], geometric_sum(h, eta))
    return mul(W[((C-D)//9) % 3], add(P, partial))


def decoded(N, C, D, K=None):
    return add(scale(mul(W[2], prefix(2*N, C, D, K)), -1),
               mul(W[1], prefix(4*N, C, D, K)))


def literal(N, C, D, K=None):
    z = (0, 0); k = 0
    while (K is None or k <= K) and source(C, D, k) < 4*N:
        z = add(z, literal_charge(source(C, D, k), N)); k += 1
    return z


def v2(x):
    return (x & -x).bit_length()-1


def energy_debt(zs):
    tail = (0, 0); debt = 0
    for z in reversed(zs):
        combined = add(z, tail)
        debt += max(0, norm(combined)-norm(z)-norm(tail))
        tail = combined
    return debt, sum(map(norm, zs)), norm(tail)


def block_word(eta, front, back, m, n):
    return [mul(W[(eta*k) % 3], front if k < m else back) for k in range(m+n)]


def block_proof():
    rows = []; slopes = set()
    for eta in (-1, 1):
        for mode in ('upper', 'lower'):
            for r in range(3):
                c = mul((1, 2), U[r])
                front, back = (W[1], c) if mode == 'upper' else (add(W[1], scale(W[(r+2) % 3], -1)), W[1])
                for m, n in product(range(3), repeat=2):
                    calc = lambda x, y: energy_debt(block_word(eta, front, back, x, y))
                    d, S, Q = calc(m, n)
                    delta = calc(m+3, n)[0]-d
                    # Three trailing back terms sum to zero and have no
                    # positive debt. Three leading front terms return to the
                    # same tail, so their nonnegative debt increment repeats.
                    need(calc(m, n+3)[0] == d and calc(m, n+6)[0] == d,
                         'back triple recurrence control')
                    need(calc(m+6, n)[0] == d+2*delta, 'front triple recurrence control')
                    need(3*d <= 2*S and 0 <= delta <= 2*norm(front), 'two-thirds induction table')
                    need(Q <= min(4, S), 'nonexpansion and sharp amplitude table')
                    slopes.add((norm(front), delta))
                    rows.append((eta, mode, r, m, n, d, S, Q, delta))
    # A formula, not a large-number search: it disproves a half-debt bound
    # for an isolated ray while retaining bounded total charge.
    family = []
    C, D, _ = ray(1, 1)
    for h in (0, 1, 4, 5, 8, 20, 100):
        N = (8*4**(9*h+10)-8186)//9
        zs = [literal_charge(source(C, D, k), N) for k in range(3*h+4)]
        expected = [W[(1-k) % 3] for k in range(3*h+2)]+[(-2, -1), (1, 2)]
        need(zs == expected and source(C, D, 3*h+4) >= 4*N, 'explicit sharp debt family word')
        need(energy_debt(zs) == (2*h+2, 3*h+8, 4), 'explicit debt family formula')
        C2, D2, _ = ray(1, 2)
        even = [literal_charge(source(C2, D2, k), N) for k in range(3*h+4)]
        expected_even = [mul(W[k % 3], (-1, -2)) for k in range(3*h+2)]+[(0, 1), (-1, -1)]
        need(even == expected_even and source(C2, D2, 3*h+4) >= 4*N, 'paired even-clock word')
        def increments(word):
            out = []; tail = (0, 0)
            for z in reversed(word):
                new = add(z, tail); out.append(norm(new)-norm(z)-norm(tail)); tail = new
            return list(reversed(out)), tail
        odd_d, odd_t = increments([scale(z, -1) for z in zs])
        even_d, even_t = increments(even)
        need(odd_d == [2, -3, -2]*h+[2, -3, -3, 0], 'odd-stage increments')
        need(even_d == [-3]*(3*h+2)+[-1, 0], 'even-stage increments')
        need(all(x+y <= 0 for x, y in zip(odd_d, even_d)), 'adjacent refinement-stage cancellation')
        cross = norm(add(odd_t, even_t))-norm(odd_t)-norm(even_t)
        need(cross == 4 and norm(add(odd_t, even_t)) == 9, 'paired initial split and endpoint')
        family.append(dict(h=h, N=str(N), positive_debt=2*h+2, source_norm=3*h+8, total_norm=4,
                           paired_source_norm=12*h+16, paired_positive_debt=4,
                           unpaired_two_ray_positive_debt=2*h+6, paired_endpoint_norm=9))
    # Merging different rays needs a separate bound: four is not inherited.
    merged = add(*(scale(literal(1864022, *ray(1, a)[:2]), (-1)**ray(1, a)[2]) for a in (1, 2)))
    need(merged == (-1, -4) and norm(merged) == 13, 'cross-ray amplitude hostile')
    target_seven = []
    for h in (0, 1, 4, 7, 20, 100):
        N = (22*4**(9*h+11)-65530)//9
        data = []
        for a in (1, 2):
            C, D, tau = ray(7, a)
            zs = [scale(literal_charge(source(C, D, k), N), (-1)**tau) for k in range(3*h+4)]
            need(source(C, D, 3*h+4) >= 4*N, 'target-seven full source interval')
            ds = []; tail = (0, 0)
            for z in reversed(zs):
                v = add(z, tail); ds.append(norm(v)-norm(z)-norm(tail)); tail = v
            data.append((zs, list(reversed(ds)), tail))
        odd = [(1, 2), (1, -1), (-2, -1)]*h+[(1, 2), (0, -1), (-1, 0), (1, 1)]
        even = [(0, 1), (-1, -1), (1, 0)]*h+[(0, 1), (-1, -1), (2, 1), (-1, 1)]
        need(data[0][0] == odd and data[1][0] == even, 'target-seven explicit phase words')
        paired = [x+y for x, y in zip(data[0][1], data[1][1])]
        need(paired == [2, -9, -5]*h+[2, -5, -4, 0], 'target-seven paired increments')
        cross = norm(add(data[0][2], data[1][2]))-norm(data[0][2])-norm(data[1][2])
        need(cross == 6 and sum(norm(z) for zs, _, _ in data for z in zs) == 12*h+14,
             'target-seven source energy and cross split')
        target_seven.append(dict(h=h, N=str(N), source_norm=12*h+14, paired_positive_debt=2*h+8,
                                 endpoint_norm=norm(add(data[0][2], data[1][2]))))
    return dict(finite_induction_rows=len(rows), front_triple_slopes=sorted(slopes),
                sharp_debt_family=family, target_seven_unbounded_pair_debt=target_seven,
                cross_ray_hostile=dict(target=1, N=1864022, charge=merged, norm=13))


def main():
    block_certificate = block_proof()
    table = {}
    for eta in (-1, 1):
        for mode in ('upper', 'lower'):
            values = set()
            for r, ell, h in product(range(3), repeat=3):
                P = mul(geometric_sum(ell, eta), U[r])
                H = mul(W[r], geometric_sum(ell, eta))
                part = mul(W[r], geometric_sum(h, eta))
                z = add(mul((1, 2), P), mul(W[1], part)) if mode == 'upper' else add(
                    mul((1, 2), P), mul(W[1], H), scale(mul(W[2], part), -1))
                values.add(norm(z))
            need(values == {0, 1, 3, 4}, ('finite proof table', eta, mode))
            table[str((eta, mode))] = sorted(values)
    targets = [u for u in range(1, 258) if primitive(u)]
    sharp = {}; comparisons = 0
    for u in targets:
        for a in (1, 2):
            C, D, tau = ray(u, a)
            for k in range(12):
                s = source(C, D, k)
                need(primitive(s) and native(s) == (u, tau+6*k, a), ('native ray', u, a, k))
                for j in range(5):
                    m = s
                    for _ in range(j): m = 4*m+1
                    need(m == shell_point(C, D-3, 3*k+j, k), 'shell regrouping')
            # Dense small scales, then both sides of every event in 24 shells.
            Ns = set(range(1, 2049))
            for i in range(24):
                for k in range(i//3+1):
                    m = shell_point(C, D-3, i, k)
                    for c in (2, 4):
                        Ns.update(n for n in (m//c, m//c+1) if n > 0)
            for N in sorted(Ns):
                z = literal(N, C, D)
                need(z == decoded(N, C, D), ('literal versus shells', u, a, N))
                need(norm(z) <= 4, ('uniform bound', u, a, N))
                comparisons += 1
                if norm(z) == 4 and str(a) not in sharp:
                    sharp[str(a)] = dict(target=u, a=a, N=N, unsigned_charge=z)
    # Arbitrarily large integer clocks: 241 shells, all event cuts, selected
    # contiguous index intervals. No floating log, real clock, or residue weight.
    intervals = 0; maxbits = 0; binchecks = 0
    for u in (1, 3, 7, 11, 17, 31, 127, 1009):
        for a in (1, 2):
            C0, D, tau = ray(u, a)
            for start in (0, 1, 5, 31):
                C = C0*64**start
                for i in (0, 1, 2, 3, 7, 19, 63, 120, 240):
                    for k in sorted({0, i//6, i//3}):
                        m = shell_point(C, D-3, i, k)
                        for N in sorted({max(1, m//2), m//2+1, max(1, m//4), m//4+1}):
                            for K in (None, 0, 1, 7, 40):
                                z = literal(N, C, D, K)
                                need(z == decoded(N, C, D, K) and norm(z) <= 4,
                                     ('clock interval', u, a, start, i, N, K))
                                intervals += 1; maxbits = max(maxbits, N.bit_length())
            vals = [source(C0, D, k) for k in range(20)]
            for k in range(20):
                for l in range(k+1, 20):
                    need(v2((vals[l]-vals[k])//2) == v2(C0)-1+6*k, 'exact source LCA')
            for q in range(130):
                bins = defaultdict(list)
                for k, s in enumerate(vals): bins[((s-1)//2) % 2**q].append(k)
                tails = [ks for ks in bins.values() if len(ks) > 1]
                need(len(tails) <= 1 and all(ks == list(range(ks[0], 20)) for ks in tails),
                     'source refinement is a singleton-and-tail comb')
                binchecks += 1
    # Hostile to singleton/prefix substitution: the general fibre need not
    # cancel in one common phase. Record, do not infer a full ROOT inequality.
    return dict(status='FINITE-EXACT checks supporting proved ray lemma',
                targets=len(targets), dense_scales=[1, 2048], event_shells=24,
                literal_comparisons=comparisons, interval_comparisons=intervals,
                largest_N_bits=maxbits, dyadic_bin_checks=binchecks,
                finite_table=table, norm_four_witnesses=sharp, checks=CHECKS,
                block_certificate=block_certificate,
                scope='one target and one native a-ray; arbitrary contiguous source index interval')


if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
