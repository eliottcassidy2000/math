"""Finite quotient loss, certified Collatz components, and selector boundaries.

All claims of undecidability are mathematical reductions in the companion
note, not inferred from this finite experiment. Run normally and with -O.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from pathlib import Path
import json
import subprocess
import sys


def need(ok, message):
    if not ok:
        raise ValueError(message)


def tournament(n, bits):
    need(type(n) is int and n >= 1 and type(bits) is int
         and 0 <= bits < 1 << (n*(n-1)//2), 'finite labelled tournament')
    a = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        a[i][j] = (bits >> k) & 1
        a[j][i] = 1-a[i][j]
    return a


def labelled_lambda(a):
    n = len(a)
    result = {(i, j): 0 for i, j in combinations(range(n), 2)}
    for i, j, k in combinations(range(n), 3):
        if a[i][j]*a[j][k]*a[k][i] or a[j][i]*a[i][k]*a[k][j]:
            for pair in ((i, j), (i, k), (j, k)):
                result[pair] += 1
    return tuple(result.values())


def path_count(a):
    n = len(a)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            for w in range(n):
                if a[v][w] and not (mask >> w) & 1:
                    dp[mask | 1 << w][w] += dp[mask][v]
    return sum(dp[-1])


def paths_by_permutations(a):
    return sum(all(a[u][v] for u, v in zip(order, order[1:]))
               for order in permutations(range(len(a))))


def trace_power(a, k):
    n = len(a)
    out = [[int(i == j) for j in range(n)] for i in range(n)]
    for _ in range(k):
        out = [[sum(out[i][v]*a[v][j] for v in range(n))
                for j in range(n)] for i in range(n)]
    return sum(out[i][i] for i in range(n))


def simple_cycles(a, length):
    total = 0
    for vertices in combinations(range(len(a)), length):
        first, *rest = vertices
        for tail in permutations(rest):
            order = (first,)+tail+(first,)
            total += all(a[u][v] for u, v in zip(order, order[1:]))
    return total


def circle_distance(value):
    remainder = value-value.numerator//value.denominator
    return min(remainder, 1-remainder)


def exact_loneliness(speeds):
    """All affine pieces and pair intersections; no sampled-grid maximum."""
    knots = sorted({F(k, 2*v) for v in speeds for k in range(2*v+1)})
    candidates = set(knots)
    for left, right in zip(knots, knots[1:]):
        middle = (left+right)/2
        lines = []
        for v in speeds:
            x = v*middle
            floor = x.numerator//x.denominator
            lines.append((v, -floor) if x-floor < F(1, 2) else (-v, floor+1))
        for (a, b), (c, d) in combinations(lines, 2):
            if a != c:
                point = F(d-b, a-c)
                if left <= point <= right:
                    candidates.add(point)
    return max(min(circle_distance(v*t) for v in speeds) for t in candidates)


class Components:
    """Finite checked-edge quotient; the graph, not just its labels, is kept."""
    def __init__(self, vertices):
        self.parent = {v: v for v in vertices}

    def find(self, v):
        self.parent.setdefault(v, v)
        if self.parent[v] != v:
            self.parent[v] = self.find(self.parent[v])
        return self.parent[v]

    def add(self, x, y):
        a, b = self.find(x), self.find(y)
        self.parent[max(a, b)] = min(a, b)


def step(n):
    need(type(n) is int and n > 0 and n % 2, 'exact positive odd source')
    out = 3*n+1
    a = (out & -out).bit_length()-1
    return out >> a, a


def replay(n, word, require_root=False):
    need(type(n) is int and n > 0 and n % 2, 'exact positive odd source')
    need(type(word) is tuple and all(type(a) is int and a > 0 for a in word),
         'finite positive exact valuation tuple')
    states = [n]
    for a in word:
        need(n != 1, 'first-hit convention excludes root padding')
        n, actual = step(n)
        need(actual == a, 'literal source guard')
        states.append(n)
    if require_root:
        need(n == 1, 'completed ROOT certificate')
    return tuple(states)


def finite_component_bank(starts, depth):
    components, edges = Components(starts), {}
    for source in starts:
        current = source
        for _ in range(depth):
            if current == 1:
                break
            target, a = step(current)
            need(current not in edges or edges[current] == (target, a), 'deterministic edge')
            edges[current] = (target, a)
            components.add(current, target)
            current = target
    return components, edges


def recorded_path(n, edges):
    states, word, seen = [n], [], {n}
    while n in edges:
        n, a = edges[n]
        word.append(a)
        states.append(n)
        if n in seen:
            break
        seen.add(n)
    return tuple(states), tuple(word)


def component_receipt(n, representative, edges):
    left, lw = recorded_path(n, edges)
    right, rw = recorded_path(representative, edges)
    position = {v: j for j, v in enumerate(right)}
    for i, target in enumerate(left):
        if target in position:
            j = position[target]
            u, v = lw[:i], rw[:j]
            need(replay(n, u)[-1] == replay(representative, v)[-1], 'checked common future')
            return u, v
    raise ValueError('component label lacks its promised finite witness')


def word_code(word):
    return ''.join('1'*a+'0' for a in word)+'0'


def decode_word(code):
    need(type(code) is str and code and set(code) <= {'0', '1'}, 'binary string')
    out, index = [], 0
    while index < len(code):
        if code[index] == '0':
            need(index+1 == len(code), 'terminal marker consumes the whole code')
            return tuple(out)
        start = index
        while index < len(code) and code[index] == '1':
            index += 1
        need(index < len(code), 'complete unary valuation')
        out.append(index-start)
        index += 1
    raise ValueError('missing terminal marker')


def words_of_edge_cost(total):
    if total == 0:
        yield ()
    for first_cost in range(2, total+1):
        for tail in words_of_edge_cost(total-first_cost):
            yield (first_cost-1,)+tail


def valid_root(n, word):
    try:
        replay(n, word, require_root=True)
        return True
    except ValueError:
        return False


def checks():
    report = {}
    for source in (True, 1.0, 0, 2):
        need(not valid_root(source, ()), 'empty ROOT word retains exact source typing')
    report['empty_root_source_type_hostiles'] = 4
    small = {}
    for n in (3, 4, 5):
        fibres = {}
        for bits in range(1 << (n*(n-1)//2)):
            a = tournament(n, bits)
            key, h = labelled_lambda(a), path_count(a)
            need(h == paths_by_permutations(a), 'independent exact Hamiltonian paths')
            need(key not in fibres or fibres[key] == h, 'finite small-order factorization')
            fibres[key] = h
        small[n] = len(fibres)
    report['all_labelled_tournaments_orders_3_4_5_lambda_fibres'] = small
    left, right = (tournament(7, bits) for bits in (1354195, 1351956))
    need(labelled_lambda(left) == labelled_lambda(right), 'identical labelled lambda')
    need((path_count(left), path_count(right)) == (109, 111), 'quotient-loss hostile')
    need((paths_by_permutations(left), paths_by_permutations(right)) == (109, 111),
         'independent permutation witness')
    for i, j in combinations(range(7), 2):
        need((left[i][j] != right[i][j]) == (j < 4), 'precisely four-vertex reversal')
    report['lambda_hostile_masks_H'] = [[1354195, 109], [1351956, 111]]
    report['lambda_hostile_simple_c7'] = [simple_cycles(a, 7) for a in (left, right)]
    strong4 = [[0, 1, 1, 0], [0, 0, 1, 1], [0, 0, 0, 1], [1, 0, 0, 0]]
    need(trace_power(strong4, 7) == 14 and simple_cycles(strong4, 7) == 0,
         'closed-walk statistic is not a simple-cycle count')
    padded = [[0]*7 for _ in range(7)]
    for i in range(7):
        for j in range(i+1, 7):
            padded[i][j] = strong4[i][j] if j < 4 else 1
            padded[j][i] = 1-padded[i][j]
    need(trace_power(padded, 7) == 14 and simple_cycles(padded, 7) == 0,
         'same hostile at historical sample order7')
    report['walk7_hostile_trace_div7_and_simple_count'] = [2, 0]
    groups = ((0, 2), (100,))
    values = [F(v) for group in groups for v in group]
    mean = sum(values)/len(values)
    total = sum((v-mean)**2 for v in values)/len(values)
    within = sum(sum((F(v)-sum(group)/F(len(group)))**2 for v in group)
                 for group in groups)/len(values)
    between = sum(len(group)*(sum(group)/F(len(group))-mean)**2
                  for group in groups)/len(values)
    need(within == F(2, 3) and within+between == total, 'singleton-correct total variance')
    report['exact_variance_singleton_control_within'] = str(within)

    need(exact_loneliness((3, 6)) == F(1, 3), 'unrestricted C-prime is false')
    need(exact_loneliness((2,)) == F(1, 2), 'one-runner scaling boundary')
    scaled = 0
    for size in (1, 2, 3):
        for speeds in combinations(range(1, 7), size):
            base = exact_loneliness(speeds)
            for multiplier in (2, 3):
                need(exact_loneliness(tuple(multiplier*v for v in speeds)) == base,
                     'exact gcd/scale invariance control')
                scaled += 1
    need(circle_distance(F(0)) <= F(1, 3) and circle_distance(F(1)) <= F(1, 3)
         and circle_distance(F(1, 2)) > F(1, 3)-F(1, 2),
         'separate endpoint circle distances do not retain a common lifted centre')
    report['LRC_nonprimitive_counterexamples'] = {'n3_S3_6': '1/3', 'n2_S2': '1/2'}
    report['exact_LRC_scale_controls_subsets_1_to_6_size_1_to_3_scales_2_3'] = scaled

    starts = tuple(range(1, 128, 2))
    previous = {n: n for n in starts}
    stage_rows, witness_count = [], 0
    for depth in (0, 1, 2, 4, 8, 16, 32, 64):
        components, edges = finite_component_bank(starts, depth)
        labels = {n: components.find(n) for n in starts}
        need(all(labels[n] <= previous[n] for n in starts), 'monotone certified minima')
        for n, label in labels.items():
            u, v = component_receipt(n, label, edges)
            if label == 1:
                need(not v, 'ROOT component supplies an actual first-hit source word')
                replay(n, u, require_root=True)
            witness_count += 1
        stage_rows.append([depth, len(edges), sum(label == 1 for label in labels.values())])
        previous = labels
    report['component_stages_depth_edges_rooted_starts'] = stage_rows
    report['independently_replayed_component_receipts'] = witness_count
    # A local sibling automaton BREAK is not a global actual-edge obstruction.
    root_words = {21: (6,), 3: (1, 4), 7: (1, 1, 2, 3, 4)}
    boundary = {}
    for source, sign in ((21, 1), (3, -1), (7, -1)):
        states = replay(source, root_words[source], require_root=True)
        for a, b in zip(states, states[1:]):
            boundary[a] = boundary.get(a, 0)+sign
            boundary[b] = boundary.get(b, 0)-sign
    boundary = {n: coefficient for n, coefficient in boundary.items() if coefficient}
    need(boundary == {21: 1, 1: 1, 3: -1, 7: -1},
         'actual chains repair the fusion relation despite local BREAK')
    need(step(7)[0] == 11 and step(11)[0] == 17 and 3*11-1 == 32
         and 17 == 2**4*1+1, 'exact x7 sibling state plus-sign exponent4 BREAK')
    report['two_sheet_local_BREAK_actual_repair'] = {
        'x': 7, 'state': [17, 1], 'sign': 1, 'exponent': 4,
        'actual_boundary_Q21_minus_Q3_minus_Q7': boundary}
    # Finite late-change control for the all-length halting-pair proof in the note.
    late = []
    for stage in (0, 1, 5, 19, 20, 100):
        dsu = Components(range(1, 9))
        for e, halt_stage in ((0, 1), (1, 5), (3, 20)):
            if halt_stage <= stage:
                dsu.add(2*e+1, 2*e+2)
        late.append([stage, dsu.find(8)])
    need(late[-3:] == [[19, 8], [20, 7], [100, 7]], 'long stability is not a finality test')
    report['finite_delayed_pair_selector'] = late

    codebook = []
    for length in range(1, 19):
        codebook.extend((word_code(word), word) for word in words_of_edge_cost(length-1))
    codebook.sort(key=lambda item: (len(item[0]), item[0]))
    codes = {code for code, _ in codebook}
    for code, word in codebook:
        need(decode_word(code) == word, 'lossless finite certificate syntax')
        need(all(code[:i] not in codes for i in range(1, len(code))), 'prefix-free syntax')
    found = {}
    for n in range(1, 32, 2):
        for code, word in codebook:
            if valid_root(n, word):
                found[n] = len(code)
                break
    need(found[1] == 1 and found[3] == 8 and found[5] == 6 and 27 not in found,
         'decidable bounded certificate search, not a nontermination claim')
    report['prefix_codes_length_at_most_18'] = len(codebook)
    report['bounded_shortlex_root_certificate_lengths_odds_below_32'] = found
    print(json.dumps(report, indent=2, sort_keys=True))
    # Preserve the scope repair's reproduction without treating floating sample
    # statistics as exact theorems or overwriting the historical output file.
    historical = Path(__file__).resolve().parents[1]/'vitali_nonmeasurable.py'
    output = subprocess.run([sys.executable, str(historical)], check=True,
                            capture_output=True, text=True).stdout
    print('HISTORICAL SAMPLE: seed42, 20000 draws; walk7 labels and full-sample variance denominator.')
    for line in output.splitlines():
        if any(key in line for key in ('Distinct lambda graphs', 'Within-lambda variance:',
                                       'H(walk7 | lambda)', 'Fibers with walk7 variation:')):
            print(line.strip())
    print('PASS: finite exact controls; no claim that Collatz equivalence is undecidable.')


if __name__ == '__main__':
    checks()
