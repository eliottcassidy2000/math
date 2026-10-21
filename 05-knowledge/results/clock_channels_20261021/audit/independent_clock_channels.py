"""Independent literal/source-clock audit; standard library, no producer imports.

Usage: python3 independent_clock_channels.py --packet /path/to/packet
The default packet is the sibling session18-clock-channels directory.
"""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
from math import gcd, isqrt, lcm, log2
from pathlib import Path
import json
import re

parser = argparse.ArgumentParser()
parser.add_argument('--packet', type=Path, default=Path(__file__).resolve().parent.parent / 'session18-clock-channels')
args = parser.parse_args()
packet = args.packet.resolve()
proof_path = packet / 'FINITE_OBSERVERS_INTEGER_CLOCKS.md'
proof = proof_path.read_text()
checks = Counter()
EXPECTED_PROOF_SHA = 'a20ec35bcd408192091c93181d8d38c5ce8c7672982a162ac50bef93facc6304'

def need(condition, label):
    if not condition:
        raise RuntimeError(label)
    checks[label] += 1

def digest(path):
    return sha256(path.read_bytes()).hexdigest()

need(digest(proof_path) == EXPECTED_PROOF_SHA, 'frozen_proof_hash')
manifest = []
for line in (packet / 'SHA256SUMS').read_text().splitlines():
    h, relative = line.split(None, 1)
    need(digest(packet / relative) == h, 'producer_manifest_entry')
    manifest.append({'path': relative, 'sha256': h})
need(len(manifest) == 6, 'producer_manifest_count')

# Generate sources by a separate affine recurrence, not by the producer formula.
source_cache = [3]
def source(j):
    while len(source_cache) <= j:
        source_cache.append(64 * source_cache[-1] + 35)
    return source_cache[j]

def v2(n):
    need(n > 0, 'positive_v2_input')
    return (n & -n).bit_length() - 1

def primitive(n):
    if n <= 0 or not n & 1 or n % 8 == 5:
        raise ValueError('Not an inherited primitive source')
    initial_cost = 2 if n % 8 == 1 else 1
    numerator = 3 * n + 1
    need(numerator % (1 << initial_cost) == 0, 'primitive_initial_integrality')
    state = numerator >> initial_cost
    chain = [state]
    while state % 8 == 5:
        quotient, remainder = divmod(state - 1, 4)
        need(remainder == 0 and 0 < quotient < state and quotient & 1, 'primitive_ancestor')
        state = quotient
        chain.append(state)
    return state, initial_cost + 2 * (len(chain) - 1), chain

def ordinary_path(n):
    path = [n]
    for _ in range(100000):
        if path[-1] == 1:
            return path
        z = path[-1]
        path.append(3 * z + 1 if z & 1 else z // 2)
    raise RuntimeError('Literal path did not terminate within explicit audit cap')

def shortcut_path(n):
    path = [n]
    for _ in range(100000):
        if path[-1] == 1:
            return path
        z = path[-1]
        path.append((3 * z + 1) // 2 if z & 1 else z // 2)
    raise RuntimeError('Shortcut path did not terminate within explicit audit cap')

def compress_ordinary(path):
    # Retain the shortcut path by deleting precisely the value after each
    # odd multiplication; infer the fully odd word from ordinary transition runs.
    shortcut = [path[0]]
    odds = [path[0]]
    valuations = []
    i = 0
    while i < len(path) - 1:
        need(path[i] & 1 == 1, 'odd_checkpoint')
        i += 1
        need(path[i] & 1 == 0, 'ordinary_odd_multiplication')
        a = 0
        while path[i] & 1 == 0:
            i += 1
            a += 1
            shortcut.append(path[i])
        valuations.append(a)
        odds.append(path[i])
    return shortcut, odds, valuations

rows = []
for j in range(97):
    n = source(j)
    need(9 * n + 5 == 32 * (1 << (6 * j)), 'recurrence_matches_closed_formula')
    need(n > 0 and n % 8 == 3, 'source_integrality_and_primitive_guard')
    terminal, clock, chain = primitive(n)
    need((terminal, clock, len(chain) - 1) == (1, 6 * j + 3, 3 * j + 1), 'independent_P_clock')
    need(chain[0] == (3 * n + 1) // 2, 'actual_P_initial_value')
    need(all(upper == 4 * lower + 1 for upper, lower in zip(chain, chain[1:])), 'exact_B_ancestry')
    ordinary = ordinary_path(n)
    shortcut = shortcut_path(n)
    grouped, odds, word = compress_ordinary(ordinary)
    need(grouped == shortcut, 'two_independent_shortcut_paths')
    need(odds == [n, chain[0], 1] and word == [1, 6 * j + 4], 'actual_two_letter_odd_word')
    need(len(shortcut) - 1 == 6 * j + 5, 'shortcut_ROOT_clock')
    need(len(ordinary) - 1 == 6 * j + 7, 'ordinary_ROOT_clock')
    need(sum(word) == len(shortcut) - 1 and sum(word) + len(word) == len(ordinary) - 1, 'typed_cost_conversion')
    first_shortcut = next(i for i, z in enumerate(shortcut) if z < n)
    first_ordinary = next(i for i, z in enumerate(ordinary) if z < n)
    need((first_shortcut, first_ordinary) == (4, 6), 'first_literal_descents')
    need(shortcut[1:5] == [(3*n+1)//2, (9*n+5)//4, (9*n+5)//8, (9*n+5)//16], 'printed_first_four_values')
    need(all(z > n for z in shortcut[1:4]) and shortcut[4] < n, 'strict_first_descent_inequalities')
    need(odds[1] > n and odds[2] < n, 'first_odd_descent_at_ROOT')
    later = source(j + 6)
    _, later_clock, _ = primitive(later)
    need(all(n % m == later % m for m in (15, 35, 8, 840)), 'same_printed_residue_registers')
    need(later_clock - clock == 36 and later_clock % 12 == clock % 12, 'full_versus_modular_P_clock')
    need(len(shortcut_path(later)) - len(shortcut) == 36, 'same_register_shortcut_clock_gap')
    need(len(ordinary_path(later)) - len(ordinary) == 36, 'same_register_ordinary_clock_gap')
    rows.append({'j': j, 'source': n, 'P_clock': clock, 'odd_word': word, 'shortcut_ROOT': len(shortcut)-1, 'ordinary_ROOT': len(ordinary)-1, 'first_shortcut_exit': first_shortcut, 'first_ordinary_exit': first_ordinary})

# Parse the actual printed example, rather than hard-coding its decimal source.
match = re.search(r'For example s_0=(\d+) and s_6=(\d+) share the registers, with unstopped P clocks(\d+) and(\d+) and physical shortcut ROOT clocks(\d+) and(\d+)\.', proof)
need(match is not None, 'printed_example_parsed')
a, b, pa, pb, ta, tb = map(int, match.groups())
need((a, b) == (source(0), source(6)), 'printed_decimal_sources')
need((primitive(a)[1], primitive(b)[1]) == (pa, pb), 'printed_P_clocks')
need((len(shortcut_path(a))-1, len(shortcut_path(b))-1) == (ta, tb), 'printed_shortcut_clocks')

# Explicit stopped frontiers exercise the important source-already-in-frontier case.
frontiers = [{1}, {1, source(0)}, {1, source(0), source(2), source(5), 7, 999}]
frontier_rows = []
for frontier in frontiers:
    for j in range(18):
        n = source(j)
        if n in frontier:
            terminal, depth, clock = n, 0, 0
        else:
            terminal, clock, _ = primitive(n)
            depth = 1
            need(terminal in frontier, 'frontier_actual_single_edge_arrival')
        need((depth, clock) == ((0, 0) if n in frontier else (1, 6*j+3)), 'stopped_frontier_membership_control')
        if n > max(frontier):
            need((terminal, depth, clock) == (1, 1, 6*j+3), 'eventual_frontier_statement')
        frontier_rows.append({'frontier': sorted(frontier), 'j': j, 'terminal': terminal, 'depth': depth, 'clock': clock})

# Compute multiplicative order directly and test both sides of the dyadic threshold.
def order64(modulus):
    need(gcd(64, modulus) == 1, 'order_unit_domain')
    residue = 1
    for exponent in range(1, modulus + 1):
        residue = residue * 64 % modulus
        if residue == 1:
            return exponent
    raise RuntimeError('Multiplicative order not found')

period_rows = []
for odd in (1, 3, 5, 7, 9, 11, 35, 105, 315):
    period = order64(9 * odd)
    phi = sum(gcd(k, 9 * odd) == 1 for k in range(1, 9 * odd + 1))
    need(phi % period == 0 and pow(64, phi, 9 * odd) == 1, 'Euler_period_suffices')
    for a2 in (0, 1, 3, 5, 6, 11, 12, 17, 18, 31, 64, 200):
        modulus = (1 << a2) * odd
        threshold = max(0, (a2 - 5 + 5) // 6)
        js = sorted({0, threshold, threshold + 1, max(0, threshold - 1)})
        for j in js:
            difference = source(j + period) - source(j)
            need(9 * difference == 32 * (1 << (6*j)) * ((1 << (6*period)) - 1), 'source_difference_identity')
            need(difference % odd == 0, 'general_odd_divisibility')
            need(v2(difference) == 6*j+5, 'exact_two_adic_valuation')
            need((difference % modulus == 0) == (6*j+5 >= a2), 'sharp_dyadic_threshold_positive_and_negative')
        for clock_modulus in (1, 2, 3, 4, 5, 7, 12, 25):
            combined_period = lcm(period, clock_modulus // gcd(6, clock_modulus))
            j = threshold
            need((source(j+combined_period)-source(j)) % modulus == 0, 'combined_source_register')
            need(6*combined_period > 0 and (6*combined_period) % clock_modulus == 0, 'arbitrary_finite_clock_modulus')
        period_rows.append({'modulus': modulus, 'odd_part': odd, 'v2_modulus': a2, 'base_period': period, 'least_eligible_j': threshold})
need([order64(m) for m in (27, 5, 7)] == [3, 2, 1], 'printed_modular_orders')
need(lcm(15, 35, 8) == 840 and pow(64, 6, 945) == 1, 'printed_840_period')

# Independently recheck only the cited small companion examples, not completeness.
need((63 % 15, 63 % 35) == (273 % 15, 273 % 35) == (3, 28), 'companion_same_registers')
need(primitive(63)[:2] == (95, 1) and primitive(273)[:2] == (51, 4), 'companion_typed_P_examples')
need(95 % 105 != 51 % 105, 'companion_observed_successors_differ')

def reduced_forms(d):
    # D=-4d. Enumerate a,c independently, imposing the exact discriminant.
    forms = []
    bound = isqrt((4*d)//3)
    for a in range(1, bound+1):
        for b in range(-a, a+1):
            if (b*b+4*d) % (4*a):
                continue
            c = (b*b+4*d)//(4*a)
            if a > c or gcd(gcd(a, abs(b)), c) != 1:
                continue
            if b < 0 and (abs(b) == a or a == c):
                continue
            forms.append((a, b, c))
    return forms

def order_exponent_two(d):
    return all(b == 0 or abs(b) == a or a == c for a,b,c in reduced_forms(d))

need(all(order_exponent_two(d) for d in (3,5,7,273)) and not order_exponent_two(11) and not order_exponent_two(63), 'explicit_small_idoneal_flags_only')
need((3*3+1) >> v2(3*3+1) == 5 and (3*7+1) >> v2(3*7+1) == 11, 'actual_flag_factor_counterexample')

# Exhaustive finite models of the exact factor criterion, independently solve f.
for size in range(1,5):
    for observer in product(range(2), repeat=size):
        labels = sorted(set(observer))
        for transition in product(range(size), repeat=size):
            fibre_ok = all(observer[x] != observer[y] or observer[transition[x]] == observer[transition[y]] for x in range(size) for y in range(size))
            existence = any(all(candidate[labels.index(observer[x])] == observer[transition[x]] for x in range(size)) for candidate in product(labels, repeat=len(labels)))
            need(fibre_ok == existence, 'finite_factor_criterion_exhaustive')
            costs = tuple(x % 3 for x in range(size))
            timed_outputs = [(observer[transition[x]], costs[x]) for x in range(size)]
            fibre_timed = all(observer[x] != observer[y] or timed_outputs[x] == timed_outputs[y] for x in range(size) for y in range(size))
            image = sorted(set(timed_outputs))
            timed_existence = any(all(candidate[labels.index(observer[x])] == timed_outputs[x] for x in range(size)) for candidate in product(image, repeat=len(labels)))
            need(fibre_timed == timed_existence, 'finite_timed_factor_criterion_exhaustive')

# Entropy control for static observers of the declared finite uniform experiment.
for J in (1,2,3,5,17,31,97):
    for q in (1,2,3,7,100):
        observer = [j % q for j in range(J)]
        fibre_sizes = Counter(observer)
        conditional_entropy = sum((size/J)*log2(size) for size in fibre_sizes.values())
        need(conditional_entropy + 1e-12 >= log2(J)-log2(q), 'finite_static_entropy_bound_control')
        need(len({6*j+3 for j in range(J)}) == J, 'clock_injectivity_finite_control')

# Guardrails are literal parts of the frozen statement, not inferred stronger claims.
for text in (
    'not a universal Collatz proof',
    'Small s_j that already lie in V instead have stopped first-arrival depth and clock zero.',
    'finite controller with an **unbounded accumulated integer counter**',
    'It also does not show that a sufficient convergence proof must store every exact clock',
    'Neither universal idoneal completeness nor a uniform annulus inequality is used.',
):
    need(text in proof, 'explicit_scope_guardrail')

producer = json.loads((packet/'verification.json').read_text())
optimized = json.loads((packet/'verification.optimized.json').read_text())
need(producer == optimized and producer['total_checks'] == 263, 'producer_receipt_modes_and_count')
for supplied in producer['rows']:
    independently = rows[supplied['j']]
    need(all(independently[key] == value for key, value in supplied.items()), 'producer_rows_match_independent_replay')

print(json.dumps({
    'status':'PASS',
    'proof_sha256':digest(proof_path),
    'method':'Recurrence-generated literal sources; ordinary trajectories and separately replayed shortcuts; B-ancestry replay; parsed printed source/clocks; exact two-sided dyadic controls; exhaustive finite factor controls. No producer import.',
    'checks':dict(sorted(checks.items())),
    'total_checks':sum(checks.values()),
    'rows':rows,
    'parsed_printed_example':{'sources':[a,b], 'P_clocks':[pa,pb], 'shortcut_clocks':[ta,tb]},
    'frontier_controls':frontier_rows,
    'modular_controls':period_rows,
    'producer_manifest':manifest,
    'scope':{'all_j_proofs':'Reviewed algebraically in AUDIT.md; finite checks are controls only.', 'actual_627_frontier_membership':'Not needed or claimed; theorem applies to any fixed finite frontier containing 1.', 'global_idoneal_completeness_audited':False, 'finite_automata_counters_or_parameterized_rules_impossible':False}
},indent=2,sort_keys=True))
