#!/usr/bin/env python3
"""Exact q=1 Taylor jets of odd-cycle Burnside series, by connected diagrams.

Standard library only. All series calculations use Fraction; no floating point.
See README.md for proofs, limits, inherited sources, and reproducibility.
"""
from argparse import ArgumentParser
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import comb, factorial, gcd, lcm
from time import perf_counter


def multiply(a, b):
    n = len(a) - 1
    out = [F(0)] * (n + 1)
    aa = [(i, x) for i, x in enumerate(a) if x]
    bb = [(i, x) for i, x in enumerate(b) if x]
    for i, x in aa:
        for j, y in bb:
            if i + j <= n:
                out[i + j] += x * y
    return out


def compositions(total, slots):
    if slots == 0:
        if total == 0:
            yield ()
    elif slots == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for rest in compositions(total - first, slots - 1):
                yield (first,) + rest


def connected(k, edges):
    seen = {0}
    changed = True
    while changed:
        changed = False
        for u, v in edges:
            if (u in seen) != (v in seen):
                seen.update((u, v))
                changed = True
    return len(seen) == k


@lru_cache(None)
def primitive_subsets(d, b):
    """Psi_b(d): b-subsets of {1,...,d} with gcd(d, all elements)=1."""
    return comb(d, b) - sum(primitive_subsets(c, b)
                            for c in range(1, d) if d % c == 0)


class ClusterCompiler:
    def __init__(self, n, core=False, cycle_weight=1):
        self.n = n
        self.core = core
        self.cycle_weight = F(cycle_weight)
        self.divisors = tuple(range(1, n + 1, 2))
        self.vertex_cache = {}
        self.divisor_cache = {}

    def vertex(self, a, d):
        """V[a,d]=sum_{p odd,d|p} binom((p-1)/2,a) x^p/p."""
        key = a, d
        if key not in self.vertex_cache:
            arr = [F(0)] * (self.n + 1)
            for p in range(d, self.n + 1, 2 * d):
                if not (self.core and p == 1):
                    arr[p] = F(comb((p - 1) // 2, a), p)
            self.vertex_cache[key] = arr
        return self.vertex_cache[key]

    def divisor_profiles(self, k, edges, powers):
        """Collapse edge-divisor assignments by their incident lcm vector."""
        key = k, edges, powers
        if key in self.divisor_cache:
            return self.divisor_cache[key]
        options = [[(d, primitive_subsets(d, b)) for d in self.divisors
                    if primitive_subsets(d, b)] for b in powers]
        profiles = defaultdict(int)
        for choices in product(*options):
            ls = [1] * k
            weight = 1
            for (u, v), (d, w) in zip(edges, choices):
                ls[u] = lcm(ls[u], d)
                ls[v] = lcm(ls[v], d)
                weight *= w
            minimum_mass = sum(max(3, d) if self.core else d for d in ls)
            if minimum_mass <= self.n:
                profiles[tuple(ls)] += weight
        self.divisor_cache[key] = dict(profiles)
        return self.divisor_cache[key]

    def log_layer(self, r):
        """[s^r] log A(x,1+s), truncated at x^n."""
        if r == 0:
            return [self.cycle_weight * x for x in self.vertex(0, 1)]
        out = [self.cycle_weight * x for x in self.vertex(r, 1)]  # one-vertex graph, no edges
        for k in range(2, min(r + 1, self.n) + 1):
            possible = tuple(combinations(range(k), 2))
            for e in range(k - 1, min(r, len(possible)) + 1):
                for edges in combinations(possible, e):
                    if not connected(k, edges):
                        continue
                    # e mandatory edge marks, then distribute the remaining marks.
                    for extra in compositions(r - e, e + k):
                        powers = tuple(1 + t for t in extra[:e])
                        marks = extra[e:]
                        for ds, weight in self.divisor_profiles(k, edges, powers).items():
                            arr = [F(1)] + [F(0)] * self.n
                            for a, d in zip(marks, ds):
                                arr = multiply(arr, self.vertex(a, d))
                            scale = F(weight, factorial(k)) * self.cycle_weight ** k
                            for n, value in enumerate(arr):
                                out[n] += scale * value
        return out

    def explicit_first_two(self):
        """Independent closed Lambert/lcm formulas for log layers 1 and 2."""
        u = self.cycle_weight
        first = [u * x for x in self.vertex(1, 1)]
        second = [u * x for x in self.vertex(2, 1)]
        for d in self.divisors:
            phi = primitive_subsets(d, 1)
            psi2 = primitive_subsets(d, 2)
            vv = multiply(self.vertex(0, d), self.vertex(0, d))
            hv = multiply(self.vertex(1, d), self.vertex(0, d))
            for n in range(self.n + 1):
                first[n] += u ** 2 * F(phi, 2) * vv[n]
                second[n] += u ** 2 * (F(psi2, 2) * vv[n] + phi * hv[n])
            for e in self.divisors:
                if d + e + lcm(d, e) > self.n:
                    continue
                triple = multiply(multiply(self.vertex(0, lcm(d, e)),
                                            self.vertex(0, d)), self.vertex(0, e))
                scale = u ** 3 * F(phi * primitive_subsets(e, 1), 2)
                for n in range(self.n + 1):
                    second[n] += scale * triple[n]
        return first, second


def reconstruct(log_layers):
    """n B[n,r]=sum_{m=1}^n sum_{j=0}^r m L[m,j] B[n-m,r-j]."""
    rmax, nmax = len(log_layers) - 1, len(log_layers[0]) - 1
    b = [[F(0)] * (nmax + 1) for _ in range(rmax + 1)]
    b[0][0] = F(1)
    for n in range(1, nmax + 1):
        for r in range(rmax + 1):
            b[r][n] = sum((m * log_layers[j][m] * b[r - j][n - m]
                           for m in range(1, n + 1) for j in range(r + 1)), F(0)) / n
    return b


def odd_partitions(n, minimum=1):
    if n == 0:
        yield ()
    for p in range(minimum, n + 1, 2):
        for tail in odd_partitions(n - p, p):
            yield (p,) + tail


def burnside_control(nmax, rmax, core=False, cycle_weight=1):
    """Independent odd-partition control, not used by the compiler."""
    b = [[F(0)] * (nmax + 1) for _ in range(rmax + 1)]
    a = []
    for n in range(nmax + 1):
        full = F(0)
        for parts in odd_partitions(n, 3 if core else 1):
            multiplicities = Counter(parts)
            z = 1
            for p, count in multiplicities.items():
                z *= p ** count * factorial(count)
            energy = sum((p - 1) // 2 for p in parts)
            energy += sum(gcd(p, q) for p, q in combinations(parts, 2))
            for r in range(rmax + 1):
                b[r][n] += F(comb(energy, r), z) * cycle_weight ** len(parts)
            full += F(2 ** energy, z) * cycle_weight ** len(parts)
        a.append(full)
    return b, a


def format_row(values):
    return ', '.join(str(v) for v in values)


def run(n, order):
    start = perf_counter()
    compiler = ClusterCompiler(n)
    layers = [compiler.log_layer(r) for r in range(order + 1)]
    compiled = reconstruct(layers)
    elapsed = perf_counter() - start
    direct, a = burnside_control(n, order)
    assert compiled == direct, 'connected diagram compiler vs partition sum'
    assert all(c >= 0 for layer in layers for c in layer), 'log positivity'
    assert all((factorial(i) * c).denominator == 1
               for row in compiled for i, c in enumerate(row)), 'integer marked permutations'
    for j, layer in enumerate(compiler.explicit_first_two(), 1):
        if j <= order:
            assert layer == layers[j], 'closed first/second formula vs graph compiler'
    for i, value in enumerate(compiled[0]):
        assert value == F(comb(2 * (i // 2), i // 2), 4 ** (i // 2))
    expected = [1, 1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056,
                903753248, 154108311168, 48542114686912]
    assert a[:min(len(a), len(expected))] == expected[:len(a)]
    # Hostile control: a finite Taylor jet is a lower bound, not the full census.
    lower = [sum(compiled[r][i] for r in range(order + 1)) for i in range(n + 1)]
    assert all(x <= y for x, y in zip(lower, a))
    strict = [i for i in range(n + 1) if lower[i] < a[i]]
    assert strict, 'test horizon must exceed the exact low-degree region'
    # Independent degree-formula check: first jet's even coefficients are divisor/harmonic sums.
    for m in range(1, n // 2 + 1):
        harmonic = lambda k: sum((F(1, 2 * j - 1) for j in range(1, k + 1)), F(0))
        formula = sum((F(primitive_subsets(d, 1), d) * harmonic(m // d)
                       for d in range(1, m + 1, 2) if m % d == 0), F(0)) / (2 * m)
        assert layers[1][2 * m] == formula
    core_n = min(n, 12)
    core = ClusterCompiler(core_n, core=True)
    core_layers = [core.log_layer(r) for r in range(min(order, 2) + 1)]
    assert reconstruct(core_layers) == burnside_control(core_n, min(order, 2), True)[0]
    for u, is_core in [(2, False), (3, True)]:
        marked = ClusterCompiler(12, core=is_core, cycle_weight=u)
        marked_logs = [marked.log_layer(r) for r in range(3)]
        assert reconstruct(marked_logs) == burnside_control(12, 2, is_core, u)[0]
        assert marked_logs[1:] == list(marked.explicit_first_two())
    print('PASS: cycle-count variable specializations u=2 (all), u=3 (core), through n=12/order=2.')
    print(f'FINITE-EXACT: all layers 0..{order}, degrees 0..{n}, agree with independent odd partitions.')
    print('PASS: closed first/second Lambert formulas; base q=1; positivity; factorial integrality.')
    print(f'PASS: odd-core-only specialization through n={core_n}, order={min(order, 2)}.')
    print(f'Compiler wall seconds (environment-specific): {elapsed:.3f}')
    print('A000568 through tested horizon: ' + format_row(a))
    upto = min(n, 12)
    for r in range(order + 1):
        print(f'n! [x^n s^{r}] A(x,1+s), n=0..{upto}: ' +
              format_row([factorial(i) * compiled[r][i] for i in range(upto + 1)]))
    for r in range(1, order + 1):
        print(f'[x^n s^{r}] log A(x,1+s), n=0..{min(n, 9)}: ' +
              format_row(layers[r][:min(n, 9) + 1]))
    witness = strict[0]
    print(f'HOSTILE: order-{order} Taylor lower bound first misses full A000568 at n={witness}: '
          f'{lower[witness]} < {a[witness]}.')
    # Every coefficient [x^n] has q-degree at most C(n,2). A complete jet recovers q=2.
    exact_n = min(n, 4)
    exact_order = comb(exact_n, 2)
    small = ClusterCompiler(exact_n)
    exact_b = reconstruct([small.log_layer(r) for r in range(exact_order + 1)])
    assert [sum(row[i] for row in exact_b) for i in range(exact_n + 1)] == a[:exact_n + 1]
    print(f'PASS: complete Taylor reconstruction of q=2 through n={exact_n} (order {exact_order}).')


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=18)
    parser.add_argument('--order', type=int, default=3)
    args = parser.parse_args()
    if args.n < 4 or args.order < 1:
        parser.error('use n >= 4 and order >= 1 for the bundled controls')
    run(args.n, args.order)
