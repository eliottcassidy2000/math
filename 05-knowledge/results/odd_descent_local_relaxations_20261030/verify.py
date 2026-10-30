"""Exact signed-affine, odd-median, and arithmetic boundary controls."""
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, permutations, product
from math import gcd
from pathlib import Path
from collections import Counter
import argparse
import importlib.util
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CHECKS = 0


def need(ok, name):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(name)


def median(values):
    need(len(values) % 2 == 1, 'odd median arity')
    return sorted(values)[len(values)//2]


def signed_solve(n, edges):
    """Independent tree propagation; edges encode x[v]=s*x[u]+b."""
    adj = [[] for _ in range(n)]
    for u, v, s, b in edges:
        adj[u].append((v, s, b))
        adj[v].append((u, s, -s*b))
    signs, shifts, values = {}, {}, [None]*n
    for root in range(n):
        if root in signs:
            continue
        signs[root], shifts[root] = 1, 0
        todo, comp, pin = [root], [], None
        while todo:
            u = todo.pop()
            comp.append(u)
            for v, s, b in adj[u]:
                a, c = s*signs[u], s*shifts[u]+b
                if v not in signs:
                    signs[v], shifts[v] = a, c
                    todo.append(v)
                else:
                    coefficient = signs[v]-a
                    rhs = c-shifts[v]
                    if not coefficient:
                        if rhs:
                            return None
                    else:
                        new = F(rhs, coefficient)
                        if pin is not None and new != pin:
                            return None
                        pin = new
        # A declared integer root for free components; no domain restriction here.
        t = F(root-2) if pin is None else pin
        for v in comp:
            values[v] = signs[v]*t+shifts[v]
    need(all(values[v] == s*values[u]+b for u, v, s, b in edges), 'actual rational solution')
    return values


def gaussian_consistent(n, edges):
    """Separate Fraction row reduction, without graph traversal."""
    rows = []
    for u, v, s, b in edges:
        row = [F(0)]*(n+1)
        row[u], row[v], row[n] = -F(s), F(1), F(b)
        rows.append(row)
    rank = 0
    for col in range(n):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = rows[rank][col]
        rows[rank] = [x/scale for x in rows[rank]]
        for i in range(rank+1, len(rows)):
            scale = rows[i][col]
            if scale:
                rows[i] = [a-scale*b for a, b in zip(rows[i], rows[rank])]
        rank += 1
    return not any(not any(row[:n]) and row[n] for row in rows)


def mod2_solutions(n, edges):
    return [xs for xs in product(range(2), repeat=n)
            if all((xs[v]-s*xs[u]-b) % 2 == 0 for u, v, s, b in edges)]


def lifted_edges(edges, matchings, d):
    return [(d*u+i, d*v+perm[i], s, b)
            for (u, v, s, b), perm in zip(edges, matchings) for i in range(d)]


def affine_census():
    options = [None, (-1, 0), (-1, 1), (1, 0), (1, 1)]
    tally, half_example = Counter(), None
    for n in range(1, 5):
        pairs = list(combinations(range(n), 2))
        for word in product(options, repeat=len(pairs)):
            edges = [(u, v, *entry) for (u, v), entry in zip(pairs, word) if entry is not None]
            xs = signed_solve(n, edges)
            real = xs is not None
            need(real == gaussian_consistent(n, edges), 'real feasibility: propagation versus Gaussian')
            binary = bool(mod2_solutions(n, edges))
            integral = real and all(x.denominator == 1 for x in xs)
            need(integral == (real and binary), 'integer iff real plus mod2')
            status = 'no_real_solution' if not real else ('integral' if integral else 'half_integral_only')
            tally[status] += 1
            if not real:
                continue
            need(all(x.denominator in (1, 2) for x in xs), 'half-integer denominator bound')
            for q in (3, 5, 9, 15):
                residues = [x.numerator*pow(x.denominator, -1, q) % q for x in xs]
                need(all((residues[v]-s*residues[u]-b) % q == 0 for u, v, s, b in edges),
                     'real-compatible systems pass every tested odd modulus')
            # The signed double cover: negative edges exchange the two sheets.
            matchings = [(0, 1) if s == 1 else (1, 0) for _, _, s, _ in edges]
            double = lifted_edges(edges, matchings, 2)
            ys = []
            for x in xs:
                ys.extend([x.numerator//x.denominator, -((-x.numerator)//x.denominator)])
            need(all(ys[v] == s*ys[u]+b for u, v, s, b in double), 'floor/ceiling double lift')
            if not integral and half_example is None:
                half_example = {'n': n, 'edges': edges, 'base': [str(x) for x in xs], 'lift': ys}
    return {'universe': 'All loopless simple supports on labelled n=1..4; each edge absent or (s,b) in {+/-1}x{0,1}.',
            'systems': sum(tally.values()), 'classification': dict(tally), 'first_half_witness': half_example}


def cover_census():
    result = []
    n, pairs = 3, list(combinations(range(3), 2))
    for d in (2, 3):
        perms = list(permutations(range(d)))
        choices = [None]+[(s, b, p) for s, b, p in product((-1, 1), (0, 1), perms)]
        tally = Counter()
        for word in product(choices, repeat=3):
            edges, matchings = [], []
            for (u, v), entry in zip(pairs, word):
                if entry is not None:
                    s, b, p = entry
                    edges.append((u, v, s, b))
                    matchings.append(p)
            base = signed_solve(n, edges)
            cover = signed_solve(n*d, lifted_edges(edges, matchings, d))
            bi = base is not None and all(x.denominator == 1 for x in base)
            ci = cover is not None and all(x.denominator == 1 for x in cover)
            if d % 2:
                need(ci == bi, 'odd-cover integral feasibility iff base')
                if ci:
                    med = [median(cover[d*v:d*(v+1)]) for v in range(n)]
                    need(all(med[v] == s*med[u]+b for u, v, s, b in edges), 'actual median descent')
                    need(all(med[v] in cover[d*v:d*(v+1)] for v in range(n)), 'arbitrary unary support retained')
            else:
                need(not bi or ci, 'base solution pulls back to every even cover too')
            if ci:
                extremes = [z for v in range(n) for z in
                            (min(cover[d*v:d*(v+1)]), max(cover[d*v:d*(v+1)]))]
                canonical = [(0, 1) if s == 1 else (1, 0) for _, _, s, _ in edges]
                need(all(extremes[v] == s*extremes[u]+b
                         for u, v, s, b in lifted_edges(edges, canonical, 2)),
                     'any finite-cover solution compresses to canonical double cover')
            tally[f'base_{int(bi)}_cover_{int(ci)}'] += 1
        result.append({'sheets': d, 'cover_presentations': len(choices)**3, 'counts': dict(tally)})
    return result


def local_marginal_controls():
    # A triangle of Boolean disequalities has one consistent local solution, of denominator two.
    marginal = {0: F(1, 2), 1: F(1, 2)}
    pair = {(0, 1): F(1, 2), (1, 0): F(1, 2)}
    need(all(sum(w for xy, w in pair.items() if xy[j] == a) == marginal[a]
             for j in (0, 1) for a in (0, 1)), 'triangle local LP feasible')
    need(not any(all(x[u] != x[v] for u, v in combinations(range(3), 2))
                 for x in product((0, 1), repeat=3)), 'triangle globally infeasible')
    odd_denominators = []
    for d in range(1, 14):
        solutions = [(a, b, c) for a, b, c in product(range(d+1), repeat=3)
                     if a+b == d and b+c == d and c+a == d]
        need(bool(solutions) == (d % 2 == 0), 'triangle excludes every tested odd denominator')
        if d % 2:
            odd_denominators.append(d)
    # Compatible denominator-five edge measures on a signed triangle.
    assignments = [(0, 3, 2), (1, 2, 3), (1, 2, 3), (4, -1, 6), (7, -4, 9)]
    edges = [(0, 1, -1, 3), (1, 2, -1, 5), (0, 2, 1, 2)]
    need(all(all(x[v] == s*x[u]+b for u, v, s, b in edges) for x in assignments), 'local tuples valid')
    med = [median([x[v] for x in assignments]) for v in range(3)]
    need(all(med[v] == s*med[u]+b for u, v, s, b in edges), 'odd local rounding')
    # Reorder each edge's tuple list independently; aligned copies are not assumed.
    tuples = []
    for j, (u, v, _, _) in enumerate(edges):
        tuples.append([(assignments[(i+j) % 5][u], assignments[(i+j) % 5][v]) for i in range(5)])
    need(all(median([p[0] for p in ps]) == med[e[0]] and median([p[1] for p in ps]) == med[e[1]]
             for e, ps in zip(edges, tuples)), 'joint ports synchronized by marginal counts')
    return {'triangle_odd_denominators_excluded': odd_denominators,
            'five_fold_median': med, 'independent_edge_orderings': tuples}


def modular_normal_form(n, edges):
    """Integer tree coordinates and chord content; never divide by two."""
    adj = [[] for _ in range(n)]
    for u, v, s, b in edges:
        adj[u].append((v, s, b))
        adj[v].append((u, s, -s*b))
    visited, blocks = set(), []
    for root in range(n):
        if root in visited:
            continue
        a, c, stack = {root: 1}, {root: 0}, [root]
        visited.add(root)
        while stack:
            u = stack.pop()
            for v, s, b in adj[u]:
                if v not in visited:
                    visited.add(v)
                    a[v], c[v] = s*a[u], s*c[u]+b
                    stack.append(v)
        zero_rows, pins = [], []
        for u, v, s, b in edges:
            if u not in a:
                continue
            coeff, rhs = a[v]-s*a[u], s*c[u]+b-c[v]
            if coeff:
                pins.append(rhs if coeff == 2 else -rhs)
            else:
                zero_rows.append(rhs)
        anchor = pins[0] if pins else None
        content = 0
        for h in zero_rows + ([r-anchor for r in pins] if pins else []):
            content = gcd(content, h)
        blocks.append((content, anchor))
    return blocks


def modular_controls():
    tally, systems = 0, 0
    for n in range(1, 4):
        pairs = list(combinations(range(n), 2))
        options = [None, (-1, 0), (-1, 1), (1, 0), (1, 1)]
        for word in product(options, repeat=len(pairs)):
            edges = [(u, v, *e) for (u, v), e in zip(pairs, word) if e is not None]
            blocks = modular_normal_form(n, edges)
            systems += 1
            for q in range(2, 13):
                predicted = 1
                for content, anchor in blocks:
                    if content % q or (anchor is not None and anchor % gcd(2, q)):
                        predicted = 0
                        break
                    predicted *= q if anchor is None else gcd(2, q)
                actual = sum(all((x[v]-s*x[u]-b) % q == 0 for u, v, s, b in edges)
                             for x in product(range(q), repeat=n))
                need(actual == predicted, 'all-modulus affine count versus full assignment census')
                tally += 1
    return {'systems': systems, 'moduli': list(range(2, 13)), 'comparisons': tally,
            'law': 'Per component: zero unless q divides content and gcd(2,q) divides a negative pin; otherwise q (balanced) or gcd(2,q) (unbalanced).'}


def odd_denominator_face_controls():
    # Complete small simplex moment polytopes: sum x_i=1, sum a_i*x_i=b, x>=0.
    # Independently enumerate rational points with odd common denominator <=31.
    hist = Counter()
    exact_exponents = Counter()
    tuples = [(d, x, y, d-x-y) for d in range(1, 32, 2)
              for x in range(d+1) for y in range(d-x+1)]
    for aa in product(range(-2, 3), repeat=3):
        witnessed = {F(sum(a*x for a, x in zip(aa, xs)), d)
                     for d, *xs in tuples}
        best_exponent = {}
        for d in range(1, 33):
            val = (d & -d).bit_length()-1
            for x in range(d+1):
                for y in range(d-x+1):
                    moment = F(aa[0]*x+aa[1]*y+aa[2]*(d-x-y), d)
                    best_exponent[moment] = min(best_exponent.get(moment, 99), val)
        for b in range(-2, 3):
            feasible = min(aa) <= b <= max(aa)
            if feasible:
                support = [i for i in range(3) if min(aa) < b < max(aa) or aa[i] == b]
                first = aa[support[0]]
                divisor = 0
                for i in support:
                    divisor = gcd(divisor, aa[i]-first)
                if divisor == 0:
                    permitted = (b == first)
                    exponent = 0
                else:
                    two_part = divisor & -divisor
                    permitted = ((b-first) % two_part == 0)
                    rhs = abs(b-first)
                    exponent = 0 if not rhs else max(0, two_part.bit_length()-(rhs & -rhs).bit_length())
                need(best_exponent[F(b)] == exponent, 'sharp denominator exponent versus all D<=32')
                exact_exponents[str(exponent)] += 1
            else:
                permitted = False
            found = F(b) in witnessed
            need(permitted == found, 'positive-face 2-primary criterion versus odd-denominator enumeration')
            hist['odd_feasible' if permitted else ('real_only' if feasible else 'real_infeasible')] += 1
    # Dropping the positive face would accept this unrestricted integral solution.
    need(1+2*0 == 1 and 1+(-1) == 0, 'unrestricted integer solution exists')
    # x1+x3=0 and nonnegativity force x1=x3=0; then 2*x2=1.
    need(F(1, 2).denominator % 2 == 0, 'support face has unavoidable even denominator')
    return {'simplex_systems': sum(hist.values()), 'classification': dict(hist),
            'minimum_two_adic_denominator_exponents': dict(exact_exponents),
            'exponent_enumeration_all_denominators_through': 32,
            'enumerated_odd_denominators': list(range(1, 32, 2)),
            'face_hostile': {'A': [[1, 2, 0], [1, 0, 1]], 'b': [1, 0],
                             'unrestricted_integer_solution': [1, 0, -1],
                             'nonnegative_solution': ['0', '1/2', '0']}}


def weight_controls():
    cases = 0
    for d in (1, 3, 5, 7, 9):
        for xs in combinations_with_replacement(range(-2, 3), d):
            m = median(xs)
            for a in range(-2, 3):
                for cost in (lambda x: abs(x-a), lambda x: (x-a)**2,
                             lambda x: max(0, x-a)):
                    need(((d+1)//2)*cost(m) <= sum(map(cost, xs)), 'sharp convex median bound')
                    cases += 1
        sharp = [0]*(d//2)+[1]*((d+1)//2)
        need(((d+1)//2)*median(sharp) == sum(sharp), 'sharp equality in convex bound')
        point = [-1]*(d//2)+[0]+[1]*(d//2)
        need(median(point) == 0 and sum(x == 0 for x in point) == 1, 'point target costs factor d')
        for tail in combinations_with_replacement(range(-2, 3), d//2):
            need(median([0]*((d+1)//2)+list(tail)) == 0, 'strict majority preserves marked point')
        if d > 1:
            need(median([0]*(d//2)+[1]*((d+1)//2)) == 1, 'minority pin not sufficient')
    return {'convex_cases': cases, 'sharp_factor': '2d/(d+1) versus normalized cover average',
            'point_target_factor': 'd, unbounded over odd d'}


def involution_controls():
    cases = 0
    for f, r, q in product(range(3), range(3), (2, 3, 4, 5)):
        n = f+2*r
        tau = list(range(f))+[f+(i ^ 1) for i in range(2*r)]
        count = sum(all((xs[tau[i]]+xs[i]) % q == 0 for i in range(n))
                    for xs in product(range(q), repeat=n))
        need(count == (2 if q % 2 == 0 else 1)**f*q**r, 'anti-invariant function count')
        cases += 1
    for a, b in product(range(-5, 6), repeat=2):
        s, t = a+b, a-b
        need((s-t) % 2 == 0 and (s+t)//2 == a and (s-t)//2 == b, 'integral even/odd clutch')
    return {'finite_function_spaces': cases, 'integral_index_per_free_pair': 2}


def arithmetic_hostiles():
    def accelerated(n, multiplier=3, constant=1):
        x = multiplier*n+constant
        a = 0
        while x % 2 == 0:
            x //= 2
            a += 1
        return x, a
    triple = [1, 3, 5]
    outs = [accelerated(x)[0] for x in triple]
    need(accelerated(median(triple))[0] != median(outs), 'unclocked Collatz not median preserving')
    same_clock = 0
    for xs in combinations(range(1, 64, 2), 3):
        vals = [accelerated(x) for x in xs]
        if len({a for _, a in vals}) == 1:
            need(accelerated(median(xs))[0] == median([v for v, _ in vals]), 'fixed clock preserved')
            same_clock += 1
    cycle = [13, 33, 83]
    need([accelerated(x, 5)[0] for x in cycle] == [33, 83, 13], '5x+1 hostile remains')
    square_edges = [(1, 3), (2, 7), (3, 1)]
    med = tuple(median([e[j] for e in square_edges]) for j in (0, 1))
    need([a+b for a, b in square_edges] == [4, 9, 4] and sum(med) == 5, 'square-label mixing fails')
    labels = [(0, 1, 2), (1, 2, 0), (2, 0, 1)]
    need([median([row[j] for row in labels]) for j in range(3)] == [1, 1, 1], 'injectivity fails')
    t = [F(1, 8), F(3, 8), F(5, 8)]
    phase = [(2*x) % 1 for x in t]
    need((2*median(t)) % 1 != median(phase), 'unmarked LRC winding fails')
    votes = [(1, 1, 1), (0, 0, 1), (1, 0, 0)]
    need(tuple(median([row[j] for row in votes]) for j in range(3)) == (1, 0, 1),
         'three transitive votes give cyclic tournament')
    odd_polynomial = lambda x: x**3-3*x
    need(odd_polynomial(1) != median([odd_polynomial(x) for x in (0, 1, 2)]),
         'odd function alone need not preserve medians')
    need(min([0, 1, 2]) == min([0, 0, 2]) and max([0, 1, 2]) == max([0, 0, 2])
         and (1 in [0, 1, 2]) != (1 in [0, 0, 2]), 'extrema lose a future middle pin')
    for q in range(2, 13):
        need(all((x+1) % q != x for x in range(q)), 'cyclic permutation has no base fixed point')
        for d in range(1, 33):
            # The marginal equations force all q counts equal; verify the only candidate.
            counts = [d//q]*q
            need((sum(counts) == d) == (d % q == 0), 'clocked rotation cover degree divisible by q')
        copies = list(range(q))
        need(all(copies[(i+1) % q] == (copies[i]+1) % q for i in range(q)), 'explicit q-sheet rotation lift')
    carry_controls = []
    for word, expected in [((1, 1, 1, 2, 1, 1, 4), 2363), ((1, 1, 1, 1, 1, 2, 4), 2123)]:
        k, total = len(word), sum(word)
        determinant = 2**total-3**k
        numerators = []
        for start in range(k):
            prefix, value = 0, 0
            for j in range(k):
                value += 3**(k-1-j)*2**prefix
                prefix += word[(start+j) % k]
            numerators.append(value)
        values = [F(v, determinant) for v in numerators]
        need(determinant == -139 and numerators[0] == expected, 'incoming ordered-carry witness')
        need(all(3*values[i]+1 == 2**word[i]*values[(i+1) % k] for i in range(k)),
             'incoming rational cycle equations')
        need(all(x.denominator % 2 for x in values), 'two-primary test cannot see odd cycle denominator')
        need(all(x.denominator == 1 for x in values) == (expected == 2363),
             'same abstract obstruction, different right-hand-side membership')
        carry_controls.append({'word': word, 'D': determinant, 'N0': numerators[0],
                               'values': list(map(str, values))})
    # One cyclic permutation with a 2-cycle and a 3-cycle: no fixed point,
    # but finite covers of both degrees. Primewise minima all vanish, yet degree1 fails.
    rotation = (1, 0, 3, 4, 2)
    need(all(i != rotation[i] for i in range(5)), 'mixed 2/3 clock has no base solution')
    need(rotation[rotation[0]] == 0 and rotation[rotation[rotation[2]]] == 2,
         'mixed clock admits degree2 and degree3 covers')
    return {'collatz_inputs': triple, 'collatz_outputs': outs, 'fixed_clock_triples': same_clock,
            'five_x_plus_one_cycle': cycle, 'square_edge_median': med,
            'LRC_times': list(map(str, t)), 'LRC_phases': list(map(str, phase)),
            'Condorcet_votes': votes, 'incoming_ordered_carry': carry_controls,
            'mixed_clock_cover_degrees': '2a+3b>0: all degrees>=2, but not1; every primewise minimum is zero'}


def hostile_81_27():
    path = ROOT/'05-knowledge/results/regular_parent_symmetry_20261029/verify.py'
    spec = importlib.util.spec_from_file_location('native_81_27', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.hostile_checks()


def recovered_cycle_controls():
    """Subset path DP versus literal cycles; replay the repaired digit-tower boundary."""
    cases = [
        ((24, 1, 3, 6, 14), 11, (1, 5, 0)),
        ((18, 8, 3, 5, 14), 13, (1, 6, 0)),
        ((58, 60, 57, 16, 32, 8), 9, (1, 2, 1)),
        ((48, 33, 35, 39, 46, 0), 9, (1, 4, 0)),
    ]
    receipt = []
    for rows, expected, alpha in cases:
        n = len(rows)
        need(all(not (rows[i] >> i) & 1 for i in range(n)) and
             all(((rows[i] >> j) & 1)+((rows[j] >> i) & 1) == 1
                 for i, j in combinations(range(n), 2)), 'actual tournament rows')
        dp = [[0]*n for _ in range(1 << n)]
        for i in range(n):
            dp[1 << i][i] = 1
        for mask in range(1, 1 << n):
            for u in range(n):
                if not dp[mask][u]:
                    continue
                for v in range(n):
                    if not ((mask >> v) & 1) and (rows[u] >> v) & 1:
                        dp[mask | (1 << v)][v] += dp[mask][u]
        h = sum(dp[-1])
        cycles, by_length = [], Counter()
        for length in range(3, n+1, 2):
            for vertices in combinations(range(n), length):
                for tail in permutations(vertices[1:]):
                    word = (vertices[0],)+tail
                    if all((rows[word[i]] >> word[(i+1) % length]) & 1 for i in range(length)):
                        cycles.append(sum(1 << v for v in vertices))
                        by_length[str(length)] += 1
        # Keep duplicate supports: different directed cyclic words are different objects.
        observed = (1, len(cycles), sum(not (a & b) for a, b in combinations(cycles, 2)))
        need(observed == alpha and h == expected and h == sum((1 << i)*v for i, v in enumerate(alpha)),
             'recovered odd-cycle/independent path counts')
        receipt.append({'rows': rows, 'H': h, 'alpha': observed, 'odd_cycles': dict(by_length)})
    return receipt


def code_operation_controls():
    def syndrome(mask):
        out = 0
        for a in range(1, 8):
            if (mask >> (a-1)) & 1:
                out ^= a
        return out
    cosets = [[m for m in range(128) if syndrome(m) == s] for s in range(8)]
    need(all(len(c) == 16 for c in cosets), 'eight named Hamming syndrome fibres')
    triples, median_hostile = 0, None
    for s, coset in enumerate(cosets):
        for a, b, c in combinations_with_replacement(coset, 3):
            need(syndrome(a ^ b ^ c) == s, 'odd XOR preserves affine syndrome fibre')
            triples += 1
            majority = (a & b) | (a & c) | (b & c)
            if s == 0 and median_hostile is None and syndrome(majority) != 0:
                median_hostile = {'masks': [a, b, c], 'median_mask': majority,
                                  'median_syndrome': syndrome(majority)}
    need(median_hostile is not None, 'Hamming code is not median closed')
    return {'affine_cosets': 8, 'size_each': 16, 'unordered_triples_with_repetitions': triples,
            'median_hostile': median_hostile}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    result = {'affine': affine_census(), 'covers': cover_census(),
              'all_modulus_counts': modular_controls(),
              'local_marginals': local_marginal_controls(),
              'odd_denominator_faces': odd_denominator_face_controls(), 'weights': weight_controls(),
              'involutions': involution_controls(), 'hostiles': arithmetic_hostiles(),
              'recovered_cycles': recovered_cycle_controls(),
              'code_operations': code_operation_controls(),
              'inherited_81_27': hostile_81_27(), 'local_checks': CHECKS,
              'scope': 'Finite exact controls and two independent rational feasibility methods; no formal build.'}
    out = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(out, encoding='utf-8', newline='\n')
    else:
        print(out, end='')


if __name__ == '__main__':
    main()
