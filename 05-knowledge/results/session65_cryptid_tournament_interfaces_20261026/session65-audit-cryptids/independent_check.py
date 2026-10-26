#!/usr/bin/env python3
"""Independent Session65 audit. Only Python stdlib; no producer-code import.
Native execution stops on control/read events, not on the claimed macro clock.
All checks raise explicitly, so -O cannot remove them.
"""
import hashlib
import itertools as it
import json
from pathlib import Path
import random
import re

ROOT = Path(__file__).resolve().parent
PRODUCER = ROOT.parent / 'session65-cryptids'
REPO = ROOT.parent / 'session65-repo'
CODES = {
    'sandwich': '1RB1LC_1LA1RD_1RA0LC_1LB0RD',
    'bouncer_bell': '1RB0LA_0RC1LA_1RD0RA_0LB1RB',
    'irregular': '1RB0RC_1LC1RA_1RA1LD_0LC0LA',
}
SECTIONS = {'sandwich': 'AB', 'bouncer_bell': 'A', 'irregular': 'AC'}

def must(condition, detail):
    if not condition:
        raise RuntimeError(detail)

def decode_table(code):
    tokens = re.findall(r'[01][LR][ABCD]', code)
    must('_'.join(''.join(tokens[2*i:2*i+2]) for i in range(4)) == code,
         'not an exact total four-state binary table')
    return {(s, b): (int(tokens[2*i+b][0]), tokens[2*i+b][1], tokens[2*i+b][2])
            for i, s in enumerate('ABCD') for b in (0, 1)}
TABLES = {name: decode_table(code) for name, code in CODES.items()}

# Dense-value dictionary, including visited zero cells; unlike set mutation in producer.
def tick(name, tape, head, state):
    old = tape.get(head, 0)
    bit, direction, state = TABLES[name][state, old]
    tape[head] = bit
    return head + {'L': -1, 'R': 1}[direction], state, int(old != bit)

def native_macro(name, support, head, state):
    tape = {x: 1 for x in support}
    initial_state = state
    erase_run = name == 'bouncer_bell' and tape.get(head, 0) == 1
    low = high = head
    changes = 0
    for time in range(1, 1000001):
        head, state, changed = tick(name, tape, head, state)
        changes += changed
        low, high = min(low, head), max(high, head)
        if name == 'sandwich':
            finished = state in 'AB'
        elif name == 'bouncer_bell':
            finished = state == 'A' and (not erase_run or tape.get(head, 0) == 0)
        else:
            finished = state in 'AC' if initial_state == 'A' else state == 'A'
        if finished:
            return (frozenset(x for x, b in tape.items() if b), head, state, time), (changes, low, high)
    raise RuntimeError('native event stop not reached')

# Macro specification reads the ORIGINAL immutable tape, then applies one patch.
def specification(name, support, h, q):
    must(q in SECTIONS[name], 'invalid section state')
    bit = lambda x: int(x in support)
    patch = {}
    if name == 'sandwich':
        d = {'A': -1, 'B': 1}[q]
        if bit(h) == 0:
            patch[h] = 1
            result = (h-d, {'A': 'B', 'B': 'A'}[q], 1)
            branch, k, edits, delta = q+'0', 0, 1, 1
        else:
            k = 0
            while bit(h+(k+1)*d):
                k += 1
            patch.update({h+i*d: 0 for i in range(1, k+1)})
            patch[h+(k+1)*d] = 1
            result = (h+k*d, q, k+2)
            branch, edits, delta = q+'1', k+1, 1-k
    elif name == 'bouncer_bell':
        if bit(h):
            k = 1
            while bit(h-k):
                k += 1
            patch.update({h-i: 0 for i in range(k)})
            result = (h-k, 'A', k)
            branch, edits, delta = 'erase', k, -k
        else:
            k = 0
            while tuple(bit(h+1+3*k+i) for i in range(3)) == (0, 0, 1):
                k += 1
            j = h+1+3*k
            patch = {h: 1, **{h+2+3*i: 1 for i in range(k)}}
            if bit(j):
                result = (j-1, 'A', 3*k+2)
                branch, edits, delta = 'pump1', k+1, k+1
            elif bit(j+1):
                patch[j+1] = 0
                result = (j+2, 'A', 3*k+3)
                branch, edits, delta = 'pump01', k+2, k
            else:
                must(not bit(j+2), 'nonmaximal pump')
                patch[j+1] = 1
                result = (j, 'A', 3*k+5)
                branch, edits, delta = 'pump000', k+2, k+2
    else:
        if q == 'A':
            k = 0
            if bit(h):
                patch[h] = 0
                result = (h+1, 'C', 1)
                branch, edits, delta = 'A1', 1, -1
            elif bit(h+1):
                patch[h] = 1
                result = (h+2, 'A', 2)
                branch, edits, delta = 'A01', 1, 1
            else:
                patch = {h: 1, h+1: 1}
                result = (h, 'C', 2)
                branch, edits, delta = 'A00', 2, 2
        else:
            k = 0
            while (bit(h-2*k-1), bit(h-2*k)) == (0, 1):
                k += 1
            j = h-2*k
            if not bit(j):
                patch[j] = 1
                result = (j+1, 'A', 2*k+1)
                branch, edits, delta = 'sweep0', 1, 1
            else:
                must(bit(j-1) == 1, 'nonmaximal sweep')
                patch[j-1] = 0
                result = (j-2, 'A', 2*k+2)
                branch, edits, delta = 'sweep11', 1, -1
    final = frozenset((set(support) - {x for x, b in patch.items() if not b}) |
                      {x for x, b in patch.items() if b})
    return (final, *result), (branch, k, edits, delta)

COUNTS = {}
CASE_TOTAL = 0
NATIVE_TOTAL = 0

def case(name, support, head, state, group):
    global CASE_TOTAL, NATIVE_TOTAL
    support = frozenset(support)
    expected, (branch, k, edits, delta) = specification(name, support, head, state)
    actual, (changes, low, high) = native_macro(name, support, head, state)
    must(actual == expected, ('native mismatch', name, sorted(support), head, state, actual, expected))
    must((len(support ^ actual[0]), len(actual[0])-len(support)) == (edits, delta), 'edit/mass')
    must(changes == edits, 'unexpected overwritten cell in tested macro')
    must(actual[3] > 0, 'zero macro clock')
    COUNTS.setdefault(group, {}).setdefault(name, {})[branch] = COUNTS.get(group, {}).get(name, {}).get(branch, 0)+1
    CASE_TOTAL += 1
    NATIVE_TOTAL += actual[3]
    return actual, (branch, k, edits, delta), (low, high)

def support(word, origin=0):
    return frozenset(origin+i for i, b in enumerate(word) if int(b))

def macro_checks():
    for n in range(11):
        for word in it.product((0, 1), repeat=n):
            s = support(word)
            for h in range(-1, n+1):
                for name, states in SECTIONS.items():
                    for q in states:
                        case(name, s, h, q, 'producer_universe')
    for k in range(201):
        for tail in ('1', '01', '000'):
            _, info, _ = case('bouncer_bell', support('0'+'001'*k+tail), 0, 'A', 'producer_long')
            must(info[1] == k, 'pump length')
        for prefix in ('0', '11'):
            w = prefix+'01'*k
            _, info, _ = case('irregular', support(w), len(w)-1, 'C', 'producer_long')
            must(info[1] == k and info[2] == 1, 'sweep length/edits')
    must(sum(sum(x.values()) for g in ('producer_universe', 'producer_long') for x in COUNTS[g].values()) == 113645,
         'producer advertised universe cardinality')
    for n in range(8):
        for word in it.product((0, 1), repeat=n):
            for origin in (-17, 23):
                s = support(word, origin) | {-1000000, 1000000}
                for h in range(origin-2, origin+n+2):
                    for name, states in SECTIONS.items():
                        for q in states:
                            case(name, s, h, q, 'signed_heads_and_remote_sentinels')
    rng = random.Random(65065)
    for _ in range(5000):
        s = {x for x in range(-40, 41) if rng.randrange(4) == 0}
        h = rng.randrange(-45, 46)
        name = rng.choice(tuple(CODES))
        q = rng.choice(SECTIONS[name])
        case(name, s, h, q, 'seeded_sparse_signed')
    for k in (0, 1, 2, 17, 201, 512, 1024):
        for q, d in (('A', -1), ('B', 1)):
            s = {-4000, 4000} | {i*d-101 for i in range(k+1)}
            _, info, _ = case('sandwich', s, -101, q, 'extended_long')
            must(info[1] == k, 'sandwich run length')
        for tail in ('1', '01', '000'):
            _, info, _ = case('bouncer_bell', support('0'+'001'*k+tail, -101)|{-4000, 4000}, -101, 'A', 'extended_long')
            must(info[1] == k, 'extended pump length')
        for prefix in ('0', '11'):
            w = prefix+'01'*k
            actual, info, _ = case('irregular', support(w, -101)|{-4000, 4000}, -101+len(w)-1, 'C', 'extended_long')
            must(info[1] == k and info[2] == 1, 'extended sweep')
            if prefix == '0':
                must(actual[3] == 2*k+1, 'clock hostile')
    blank = {}
    for name in CODES:
        s, h, q, t, calls = frozenset(), 0, 'A', 0, 0
        while t < 10000:
            (s, h, q, elapsed), _, _ = case(name, s, h, q, 'blank_orbit')
            t += elapsed
            calls += 1
        blank[name] = dict(native_steps=t, macro_calls=calls, head=h, state=q, ones=len(s), support_min=min(s), support_max=max(s))
    receipt = json.loads((PRODUCER/'normal.json').read_text())
    must(blank == receipt['blank_orbit'], 'producer orbit receipt mismatch')
    must(COUNTS['producer_universe'] == receipt['macros']['branch_counts'], 'producer branch receipt mismatch')
    return dict(total_cases=CASE_TOTAL, native_instructions=NATIVE_TOTAL, branch_counts=COUNTS,
                independent_event_stop=True, blank_orbit=blank, random_seed=65065)

# Adjacency matrices, not the producer's edge closures.
def tournament(word):
    n = 3*len(word)
    a = [[int(i < j) for j in range(n)] for i in range(n)]
    for k, b in enumerate(word):
        if b:
            i, j = 3*k, 3*k+2
            a[i][j], a[j][i] = 0, 1
    return a

def relabel(a, p):
    return [[a[p[i]][p[j]] for j in range(len(a))] for i in range(len(a))]

def distance(a, b):
    n = len(a)
    must(len(b) == n, 'distance order')
    costs = [sum(a[i][j] != b[p[i]][p[j]] for i, j in it.combinations(range(n), 2))
             for p in it.permutations(range(n))]
    return min(costs), len(costs)

def decode_scc(a):
    # Kosaraju two-pass DFS, unlike producer bitset transitive closure.
    n = len(a)
    seen, finish = set(), []
    def dfs(v):
        seen.add(v)
        for u in range(n):
            if a[v][u] and u not in seen:
                dfs(u)
        finish.append(v)
    for v in range(n):
        if v not in seen:
            dfs(v)
    seen, components = set(), []
    def reverse_dfs(v, c):
        seen.add(v); c.append(v)
        for u in range(n):
            if a[u][v] and u not in seen:
                reverse_dfs(u, c)
    for v in reversed(finish):
        if v not in seen:
            c = []; reverse_dfs(v, c); components.append(c)
    # Kosaraju order is source-first; validate every cross-component arc.
    for i, c in enumerate(components):
        for d in components[i+1:]:
            must(all(a[x][y] for x in c for y in d), 'SCC condensation not total forward order')
    answer, singletons = [], 0
    for c in components+[[]]:
        if len(c) == 1:
            singletons += 1
        else:
            must(singletons % 3 == 0, 'undecodable singleton run')
            answer += [0]*(singletons//3)
            singletons = 0
            if c:
                must(len(c) == 3 and all(sum(a[x][y] for y in c) == 1 for x in c), 'non-C3 SCC')
                answer.append(1)
    return tuple(answer)

def twin_pairs(a):
    return [(i, j) for i, j in it.combinations(range(len(a)), 2)
            if all(a[i][x] == a[j][x] for x in range(len(a)) if x not in (i, j))]

def tournament_checks():
    pairs = perms = 0
    invariant_rows = []
    for m in range(3):
        words = list(it.product((0, 1), repeat=m))
        for u in words:
            a = tournament(u); n = len(a)
            triangles = sum(all(sum(a[x][y] for y in triple) == 1 for x in triple)
                            for triple in it.combinations(range(n), 3))
            paths = sum(all(a[p[i]][p[i+1]] for i in range(n-1)) for p in it.permutations(range(n)))
            must((triangles, paths) == (sum(u), 3**sum(u)), 'invariants')
            invariant_rows.append(dict(word=''.join(map(str, u)), c3=triangles, H=paths))
            for v in words:
                d, count = distance(a, tournament(v))
                must(d == sum(x != y for x, y in zip(u, v)), 'unlabelled distance')
                pairs += 1; perms += count
    rng = random.Random(2183)
    decoded_words = decoded_labelings = writes = 0
    for m in range(9):
        for w in it.product((0, 1), repeat=m):
            a = tournament(w)
            labels = list(range(len(a))); rng.shuffle(labels)
            for p in (tuple(reversed(range(len(a)))), tuple(labels)):
                must(decode_scc(relabel(a, p)) == w, 'SCC decode')
                decoded_labelings += 1
            decoded_words += 1
            for k in range(m):
                v = list(w); v[k] ^= 1
                b = tournament(v)
                changed = [(i, j) for i, j in it.combinations(range(len(a)), 2) if a[i][j] != b[i][j]]
                must(changed == [(3*k, 3*k+2)], 'native write not one distinguished arc')
                writes += 1
    a, b = tournament((0, 0)), tournament((1, 0))
    must((2, 3) in twin_pairs(a) and (2, 3) not in twin_pairs(b), 'clone split')
    must(twin_pairs(a) == [(0,1),(1,2),(2,3),(3,4),(4,5)], 'transitive six twin ledger')
    must(twin_pairs(b) == [(3,4),(4,5)], 'cycle plus transitive triple twin ledger')
    return dict(word_pairs=pairs, all_permutation_comparisons=perms, max_order=6,
                invariants=invariant_rows, decoded_words=decoded_words, decoded_labelings=decoded_labelings,
                one_bit_write_checks=writes, twins_00=twin_pairs(a), twins_10=twin_pairs(b))

def hostiles():
    a, _, _ = case('bouncer_bell', {1}, 0, 'A', 'hostiles')
    b, _, _ = case('bouncer_bell', {0}, 0, 'A', 'hostiles')
    must(a != b and a[3] == 2 and b[3] == 1, 'position hostile')
    must(distance(tournament((0, 1)), tournament((1, 0)))[0] == 2, 'position distance')
    # The two-step 01 case includes A0,B1, not the three-step tail-01 branch.
    must(a == (frozenset((0,1)), 0, 'A', 2), '01 exact output')
    must(b == (frozenset(), -1, 'A', 1), '10 exact output')
    sa, _ = native_macro('sandwich', (), 0, 'A')
    sb, _ = native_macro('sandwich', (), 0, 'B')
    must(sa != sb, 'state forgotten')
    sh, _ = native_macro('sandwich', (), 1, 'A')
    must(sa != sh, 'head forgotten')
    t = {3: 1}; h, q = 0, 'A'
    for _ in range(5): h, q, _ = tick('bouncer_bell', t, h, q)
    wrong = (frozenset((0,2,3)), 1, 'A')
    actual = (frozenset(x for x, bit in t.items() if bit), h, q)
    must(actual != wrong, '00 improperly substituted for 000')
    # Fixed singleton window [0,0] cannot even retain the next blank head.
    h, q, _ = tick('sandwich', {}, 0, 'A')
    must(h == 1, 'boundary exit is not rejected')
    return dict(same_population_position_loss=True, missing_head=True, missing_control=True,
                terminal_000_guard_required=True, fixed_window_exit=True,
                all_k_clock_proof='C on 0(01)^k: only cell zero changes, exactly 2k+1 instructions')

def marked_native_steps():
    checked = accepted = exits = 0
    for m in range(1, 6):
        for w in it.product((0, 1), repeat=m):
            a = tournament(w)
            for h in range(m):
                for name in CODES:
                    for q in 'ABCD':
                        t = {i: bit for i, bit in enumerate(w)}
                        h1, q1, changed = tick(name, t, h, q)
                        w1 = tuple(t[i] for i in range(m))
                        b = tournament(w1)
                        differences = [(i,j) for i,j in it.combinations(range(3*m), 2) if a[i][j] != b[i][j]]
                        must(differences == ([(3*h,3*h+2)] if changed else []), 'native marked arc rewrite')
                        must(decode_scc(b) == w1 and q1 in 'ABCD', 'native marked tape/control')
                        checked += 1
                        if 0 <= h1 < m:
                            accepted += 1
                        else:
                            exits += 1
    for invalid in ([[0]], [[0,1],[0,0]],
                    [[int(i != j and (j-i)%5 in (1,2)) for j in range(5)] for i in range(5)]):
        rejected = False
        try:
            decode_scc(invalid)
        except RuntimeError:
            rejected = True
        must(rejected, 'decoder accepted non-image tournament')
    return dict(native_transitions_checked=checked, valid_in_window=accepted,
                explicit_head_exits=exits, nonimage_tournaments_rejected=3)

def input_pins():
    expected = {
        'REPORT.md':'5895328ca829b7937983947d08caffd57d3b85f0ef102cf18b49a7843915a08d',
        'SOURCES.json':'f47ed0f64a976d8f42e1f25f78a9ca873857dbcaf7100dbcb6862e5af1ea3280',
        'check.py':'ad88b82106aa9d517f7f2f4a56d50abf00b0086eacff80c4d0e64a6c3b6fd4b1',
        'normal.json':'a3ae717a869a8f115aaa8347969308a738ad8b67d1c51a9e1caf6f070299f128',
        'optimized.json':'a3ae717a869a8f115aaa8347969308a738ad8b67d1c51a9e1caf6f070299f128',
    }
    for path, wanted in expected.items():
        must(hashlib.sha256((PRODUCER/path).read_bytes()).hexdigest() == wanted, ('producer pin changed', path))
    deps = {
        'THM-2183-order-join-is-an-exact-tournament-metric-product.md':'00c34a935dc73191edbd7088e1bd9ce2029af23c9d89c9889e1af5b5425499bc',
        'THM-2195-transitive-quotients-exactly-control-universal-substitution-products.md':'60ddaa62c9d6e594b994e2e9c7ab2600edf9c763146c950d5ebf009ab025a006',
    }
    for path, wanted in deps.items():
        must(hashlib.sha256((REPO/'01-canon/theorems'/path).read_bytes()).hexdigest() == wanted, ('dependency pin changed', path))
    return dict(producer=expected, dependencies=deps)

if __name__ == '__main__':
    pins = input_pins()
    macro = macro_checks()
    graph = tournament_checks()
    hostile = hostiles()
    marked = marked_native_steps()
    macro['total_cases'] = CASE_TOTAL
    macro['native_instructions'] = NATIVE_TOTAL
    print(json.dumps(dict(status='PASS', pins=pins, native=macro, tournaments=graph,
                         hostiles=hostile, marked_native=marked, total_cases_including_hostiles=CASE_TOTAL,
                         total_native_instructions_including_hostiles=NATIVE_TOTAL), sort_keys=True, indent=2))
