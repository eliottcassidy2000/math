#!/usr/bin/env python3
"""Order-eleven residual certificate: only proved exceptional cell families.

The only external data file read is a byte-pinned regular-nine certificate;
the lower strong-error histograms are separately declared inherited mathematical
inputs.
No full order-eleven enumeration, class count, or regular-eleven count is used.
Explicit checks remain live with -O. Python standard library only.
"""
from collections import Counter
from itertools import combinations, permutations, product
from math import factorial
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
CHECKS = Counter()

def require(test, label):
    CHECKS[label] += 1
    if not test:
        raise RuntimeError(label)

def edge(T, a, b):
    return (T[a] >> b) & 1

def validate(T):
    n = len(T)
    require(all(0 <= r < 1 << n for r in T), 'bit bounds')
    require(all(not edge(T, i, i) for i in range(n)), 'no loops')
    require(all(edge(T, i, j) + edge(T, j, i) == 1 for i in range(n) for j in range(i)), 'tournament axioms')

def cells(T):
    color = [0] * len(T)
    while True:
        labels = sorted(set(color))
        keys = [(color[v],) + tuple(sum(edge(T, v, w) for w in range(len(T)) if color[w] == c) for c in labels) for v in range(len(T))]
        ranks = {k: i for i, k in enumerate(sorted(set(keys)))}
        new = [ranks[k] for k in keys]
        if len(set(new)) == len(set(color)):
            return tuple(tuple(v for v in range(len(T)) if new[v] == c) for c in sorted(set(new)))
        color = new

def strong(T):
    def accessible(reverse):
        reached = {0}
        todo = [0]
        while todo:
            u = todo.pop()
            for v in range(len(T)):
                if v not in reached and (edge(T, v, u) if reverse else edge(T, u, v)):
                    reached.add(v)
                    todo.append(v)
        return len(reached) == len(T)
    return accessible(False) and accessible(True)

def relabel(T, p):
    R = [0] * len(T)
    for i in range(len(T)):
        for j in range(len(T)):
            R[p[i]] |= edge(T, i, j) << p[j]
    return tuple(R)

def isomorphisms(A, B, first=False):
    """Full arc-preserving search; no presentation blocks or refinement assumed."""
    if len(A) != len(B):
        return ()
    n = len(A)
    da, db = [r.bit_count() for r in A], [r.bit_count() for r in B]
    targets = [[j for j in range(n) if da[i] == db[j]] for i in range(n)]
    order = sorted(range(n), key=lambda i: (len(targets[i]), i))
    p, answers = [-1] * n, []
    def go(k, used):
        if k == n:
            answers.append(tuple(p))
            return first
        i = order[k]
        for j in targets[i]:
            if (used >> j) & 1:
                continue
            if any(edge(A, i, h) != edge(B, j, p[h]) for h in order[:k]):
                continue
            p[i] = j
            if go(k + 1, used | (1 << j)):
                return True
            p[i] = -1
        return False
    go(0, 0)
    return tuple(answers)

def action_error(T):
    actions = isomorphisms(T, T)
    require(len(actions) >= 1, 'full automorphism identity')
    orbits = {tuple(sorted({p[v] for p in actions})) for v in range(len(T))}
    sizes = {v: len(C) for C in cells(T) for v in C}
    return len(actions), sum(sizes[O[0]] - len(O) for O in orbits), tuple(sorted(map(len, orbits)))

def tournaments(n):
    pairs = list(combinations(range(n), 2))
    for code in range(1 << len(pairs)):
        T = [0] * n
        for k, (i, j) in enumerate(pairs):
            a, b = (i, j) if (code >> k) & 1 else (j, i)
            T[a] |= 1 << b
        yield tuple(T)

def carousel(n):
    return tuple(sum(1 << ((i + j) % n) for j in range(1, (n + 1) // 2)) for i in range(n))

def substitute(Q, parts):
    starts, T, k = [], [], 0
    for P in parts:
        starts.append(k)
        T += [r << k for r in P]
        k += len(P)
    for i, P in enumerate(parts):
        for j, R in enumerate(parts):
            if edge(Q, i, j):
                for v in range(starts[i], starts[i] + len(P)):
                    T[v] |= ((1 << len(R)) - 1) << starts[j]
    return tuple(T)

def odd_partitions(n, low=1):
    if n == 0:
        yield ()
    for k in range(low, n + 1, 2):
        for rest in odd_partitions(n - k, k):
            yield (k,) + rest

def digest(x):
    return hashlib.sha256(json.dumps(x, separators=(',', ':'), sort_keys=True).encode()).hexdigest()

def regular_sevens():
    # Decide all arcs touching vertex i at step i; row-score pruning is exact.
    found = set()
    def visit(i, T):
        if i == 7:
            require(all(r.bit_count() == 3 for r in T), 'regular seven terminal')
            found.add(tuple(T))
            return
        undecided = tuple(range(i + 1, 7))
        need = 3 - T[i].bit_count()
        if not 0 <= need <= len(undecided):
            return
        for chosen in combinations(undecided, need):
            chosen = set(chosen)
            R = T[:]
            for j in undecided:
                if j in chosen:
                    R[i] |= 1 << j
                else:
                    R[j] |= 1 << i
            if all(R[j].bit_count() <= 3 <= R[j].bit_count() + 6 - i for j in undecided):
                visit(i + 1, R)
    visit(0, [0] * 7)
    require(len(found) == 2640, 'regular seven labelled universe')
    reps, pending, group = [], set(found), tuple(permutations(range(7)))
    while pending:
        T = min(pending)
        orbit = {relabel(T, p) for p in group}
        require(orbit <= pending, 'regular seven orbit partition')
        pending -= orbit
        a, e, sizes = action_error(T)
        require(a * len(orbit) == factorial(7), 'regular seven independent orbit stabilizer')
        reps.append({'rows': T, 'automorphisms': a, 'error': e, 'orbit_sizes': sizes})
    require(sorted((r['automorphisms'], r['error']) for r in reps) == [(3, 14), (7, 0), (21, 0)], 'regular seven classification')
    return reps, {'labelled_count': len(found), 'classes': reps, 'universe_sha256': digest(sorted(found))}

def regular_nines():
    path = HERE / 'sources/session50-regular-nine-receipt.json'
    require(hashlib.sha256(path.read_bytes()).hexdigest() == '74ec9ba505971a5e6741a4d52b83031b2185ca6c1ab1d6f9b94e4c75a68eb9d9', 'regular nine input pin')
    data = json.loads(path.read_text())
    reps, hist = [], Counter()
    for component in data['components']:
        state = component['states'][0]
        P = data['states'][state]['parent']
        require(len(P) == 8 and sorted(r.bit_count() for r in P) == [3] * 4 + [4] * 4, 'almost regular eight source')
        T = list(P) + [0]
        for i, r in enumerate(P):
            if r.bit_count() == 3:
                T[i] |= 1 << 8
            else:
                T[8] |= 1 << i
        T = tuple(T)
        validate(T)
        require(all(r.bit_count() == 4 for r in T), 'regular nine forced extension')
        require(strong(T), 'regular nine strong')
        a, e, sizes = action_error(T)
        require(e == component['error'] and list(sizes) == component['root_orbit_sizes'], 'regular nine full action matches reroot certificate')
        require(not any(isomorphisms(T, R['rows'], first=True) for R in reps), 'regular nine representatives pairwise distinct')
        reps.append({'rows': T, 'automorphisms': a, 'error': e, 'orbit_sizes': sizes})
        hist[e] += 1
    require(dict(hist) == {0: 3, 18: 3, 36: 2, 72: 7}, 'regular nine error distribution')
    require(sum(k * v for k, v in hist.items()) == 630, 'regular nine inherited total independently checked')
    return reps, {'classes': reps, 'histogram': sorted(hist.items()), 'completeness': 'inherited authenticated Session50 regular-nine certificate, not established by this checker'}

def check_uniform_sectors(sevens, nines):
    hist, details = Counter(), {}
    # 9+1+1 has one rooted cycle; 7+3+1 has two differently weighted cycles.
    Qs = [Q for Q in tournaments(3) if strong(Q)]
    records = []
    Q = Qs[0]
    for R in nines:
        T = substitute(Q, [R['rows'], (0,), (0,)])
        validate(T)
        a, e, _ = action_error(T)
        require(sorted(map(len, cells(T))) == [1, 1, 9] and strong(T), 'nine plus two actual strong')
        require(a == R['automorphisms'] and e == R['error'], 'nine plus two full action')
        hist[e] += 1
        records.append({'rows': T, 'error': e, 'automorphisms': a})
    details['nine_plus_two_singletons'] = {'classes': records, 'total_error': sum(r['error'] for r in records)}
    records = []
    for Q in Qs:
        for R in sevens:
            T = substitute(Q, [R['rows'], carousel(3), (0,)])
            validate(T)
            a, e, _ = action_error(T)
            require(sorted(map(len, cells(T))) == [1, 3, 7] and strong(T), 'seven triangle singleton actual strong')
            require(a == 3 * R['automorphisms'] and e == R['error'], 'seven triangle singleton full action')
            hist[e] += 1
            records.append({'rows': T, 'error': e, 'automorphisms': a})
    details['seven_triangle_singleton'] = {'classes': records, 'total_error': sum(r['error'] for r in records)}
    G = tuple((0,) + p for p in permutations(range(1, 5)))
    actual, collapsed, classifications = set(), set(), Counter()
    for Q in tournaments(5):
        T = substitute(Q, [carousel(7)] + [(0,)] * 4)
        pattern = tuple(sorted(map(len, cells(T))))
        require(strong(Q) == strong(T), 'weighted quotient strong equivalence')
        classifications[(strong(Q), pattern)] += 1
        if strong(Q):
            if pattern == (1, 1, 1, 1, 7):
                actual.add(Q)
            else:
                require(pattern == (1, 3, 7), 'weighted five quotient only collapse')
                collapsed.add(Q)
    require(len(actual) == 528 and len(collapsed) == 16, 'rooted five actual and collapsed counts')
    pending, records, quotient_reps = set(actual), [], []
    while pending:
        Q = min(pending)
        orbit = {relabel(Q, p) for p in G}
        require(orbit <= pending and len(orbit) == 24, 'weighted rooted five free quotient orbits')
        pending -= orbit
        quotient_reps.append(Q)
        for R in sevens:
            T = substitute(Q, [R['rows']] + [(0,)] * 4)
            a, e, _ = action_error(T)
            require(a == R['automorphisms'] and e == R['error'], 'seven four singleton full action')
            require(sorted(map(len, cells(T))) == [1, 1, 1, 1, 7], 'seven four singleton actual cells')
            hist[e] += 1
            records.append({'rows': T, 'error': e, 'automorphisms': a, 'quotient': Q})
    require(len(quotient_reps) == 22, 'twenty two rooted weighted quotients')
    pending, collapse_orbits = set(collapsed), []
    while pending:
        Q = min(pending)
        orbit = {relabel(Q, p) for p in G}
        require(orbit <= pending and len(orbit) == 8, 'collapsed rooted five orbit size')
        pending -= orbit
        collapse_orbits.append(Q)
    require(len(collapse_orbits) == 2, 'collapsed quotients belong to two seven triangle cycles')
    details['seven_plus_four_singletons'] = {'classes': records, 'quotient_classes': quotient_reps,
       'classification': [{'strong': s, 'stable_sizes': p, 'presentations': k} for (s, p), k in sorted(classifications.items())],
       'collapsed_quotient_representatives': collapse_orbits, 'total_error': sum(r['error'] for r in records)}
    return hist, details

def matrices(d):
    options = [sum(1 << i for i in S) for S in combinations(range(5), d)]
    out, columns = set(), [0] * 5
    def visit(R):
        if len(R) == 5:
            require(columns == [d] * 5, 'regular matrix terminal')
            out.add(tuple(R))
            return
        for r in options:
            for j in range(5):
                columns[j] += (r >> j) & 1
            if all(c <= d <= c + 4 - len(R) for c in columns):
                visit(R + [r])
            for j in range(5):
                columns[j] -= (r >> j) & 1
    visit([])
    return out

def matrix_shift(M, a, b):
    R = [0] * 5
    for i in range(5):
        for j in range(5):
            R[(i + a) % 5] |= edge(M, i, j) << ((j + b) % 5)
    return tuple(R)

def five_presentation(M, t):
    T = list(carousel(5)) + [r << 5 for r in carousel(5)] + [0]
    for i in range(5):
        for j in range(5):
            a, b = (i, 5 + j) if edge(M, i, j) else (5 + j, i)
            T[a] |= 1 << b
    for i in range(10):
        a, b = (i, 10) if t[i // 5] else (10, i)
        T[a] |= 1 << b
    return tuple(T)

def check_five_sector():
    hist, report = Counter(), []
    require(action_error(carousel(5))[0] == 5, 'C5 full group')
    for d in (1, 2):
        universe = matrices(d)
        require(len(universe) == {1: 120, 2: 2040}[d], 'five crossing labelled matrices')
        actual_orientations = [(0, 1), (1, 0)] if d == 1 else [(0, 1)]
        collisions = Counter()
        for M in universe:
            for t in product(range(2), repeat=2):
                T = five_presentation(M, t)
                pattern = tuple(sorted(map(len, cells(T))))
                expected = (11,) if d == 2 and t == (1, 0) else (1, 5, 5)
                require(pattern == expected, 'five decorated actual cell collision rule')
                require(strong(T) == (t[0] != t[1]), 'five decorated strong iff mixed')
                collisions[(t, pattern, strong(T))] += 1
        pending, classes = set(universe), []
        while pending:
            M = min(pending)
            orbit = {matrix_shift(M, a, b) for a, b in product(range(5), repeat=2)}
            require(orbit <= pending, 'five matrix translation orbit partition')
            pending -= orbit
            stabilizer = 25 // len(orbit)
            require(stabilizer in (1, 5), 'five crossing stabilizer types')
            for t in actual_orientations:
                T = five_presentation(M, t)
                validate(T)
                a, e, _ = action_error(T)
                require(a == stabilizer, 'five decorated unrestricted full automorphism')
                require(e == (40 if a == 1 else 0), 'five decorated orbit error')
                hist[e] += 1
                classes.append({'matrix': M, 'attachment': t, 'rows': T, 'automorphisms': a, 'error': e})
        expected = {1: Counter({0: 8, 40: 8}), 2: Counter({0: 8, 40: 80})}[d]
        require(Counter(r['error'] for r in classes) == expected, 'five decorated class histogram')
        report.append({'degree': d, 'matrix_count': len(universe), 'matrix_sha256': digest(sorted(universe)), 'classes': classes,
            'presentation_classification': [{'attachment': t, 'stable_sizes': p, 'strong': s, 'count': k} for (t, p, s), k in sorted(collisions.items())]})
    # Uniform crossing cycle is exact, and is not part of the positive-error sector.
    T = five_presentation((0,) * 5, (1, 0))
    a, e, _ = action_error(T)
    require(strong(T) and sorted(map(len, cells(T))) == [1, 5, 5] and (a, e) == (25, 0), 'uniform two five cycle zero control')
    hist[0] += 1
    require(hist[40] == 88, 'all positive five decorated classes')
    return hist, {'nonuniform_sectors': report, 'uniform_zero_control': {'rows': T, 'automorphisms': a}, 'total_error': 3520}

def triangle_presentation(Q, attachment, singleton_edge, signs=(1, 1, -1), phases=(0, 0, 0)):
    T = [0] * 11
    for i in range(3):
        for x in range(3):
            T[3 * i + x] |= 1 << (3 * i + (x + 1) % 3)
    for k, (i, j) in enumerate(combinations(range(3), 2)):
        for x in range(3):
            target = (signs[k] * x + phases[k]) % 3
            for y in range(3):
                forward = (y != target) if edge(Q, i, j) else (y == target)
                a, b = (3 * i + x, 3 * j + y) if forward else (3 * j + y, 3 * i + x)
                T[a] |= 1 << b
    for h in range(2):
        for i in range(3):
            for x in range(3):
                a, b = (3 * i + x, 9 + h) if attachment[3 * h + i] else (9 + h, 3 * i + x)
                T[a] |= 1 << b
    a, b = (9, 10) if singleton_edge else (10, 9)
    T[a] |= 1 << b
    return tuple(T)

def compressed_refinement(Q, attachment, se):
    # Outcounts into each intended block, including internal outdegree on diagonal.
    W = [[0] * 5 for _ in range(5)]
    for i in range(3):
        W[i][i] = 1
        for j in range(3):
            if i != j:
                W[i][j] = 1 + edge(Q, i, j)
        for h in range(2):
            W[i][3 + h] = attachment[3 * h + i]
            W[3 + h][i] = 3 * (1 - attachment[3 * h + i])
    W[3][4], W[4][3] = se, 1 - se
    colors = [0] * 5
    while True:
        keys = [(colors[i],) + tuple(sum(W[i][j] for j in range(5) if colors[j] == c) for c in sorted(set(colors))) for i in range(5)]
        ranks = {k: i for i, k in enumerate(sorted(set(keys)))}
        new = [ranks[k] for k in keys]
        if len(set(new)) == len(set(colors)):
            weights = (3, 3, 3, 1, 1)
            return tuple(sorted(sum(weights[i] for i in range(5) if new[i] == c) for c in set(new)))
        colors = new

def skeleton_key(T):
    Q = tuple(sum((sum(edge(T, 3 * i, 3 * j + y) for y in range(3)) == 2) << j for j in range(3) if j != i) for i in range(3))
    attachment = tuple(edge(T, 3 * i, 9 + h) for h in range(2) for i in range(3))
    return Q, attachment, edge(T, 9, 10)

def check_triangles():
    skeletons, collapsed, classification = {}, [], Counter()
    for Q in tournaments(3):
        count_strong = 0
        for attachment in product(range(2), repeat=6):
            for se in range(2):
                key = (Q, attachment, se)
                T = triangle_presentation(*key)
                actual = tuple(sorted(map(len, cells(T))))
                require(actual == compressed_refinement(*key), 'compressed weighted refinement matches literal eleven')
                # Core is strongly connected. Each singleton has uniform source,
                # uniform sink, or mixed attachment to it.
                kinds = [sum(attachment[3*h:3*h+3]) for h in range(2)]
                mixed = [0 < k < 3 for k in kinds]
                expected_strong = (all(mixed) or
                    (mixed[0] and not mixed[1] and (se == (kinds[1] == 0))) or
                    (mixed[1] and not mixed[0] and (se == (kinds[0] == 3))) or
                    (not any(mixed) and kinds[0] != kinds[1] and (se == (kinds[0] == 3))))
                require(strong(T) == expected_strong, 'triangle strong attachment formula')
                count_strong += strong(T)
                cyclic = strong(Q)
                require(actual == ((1, 1, 9) if cyclic and kinds[0] in (0,3) and kinds[1] in (0,3) else (1,1,3,3,3)), 'triangle exact collapse characterization')
                classification[(cyclic, strong(T), actual)] += 1
                if not strong(T):
                    continue
                if actual == (1, 1, 9):
                    collapsed.append({'skeleton': key, 'rows': T})
                else:
                    skeletons[key] = T
        require(count_strong == 98, 'ninety eight strong attachments per majority')
    require(len(skeletons) == 780 and len(collapsed) == 4, 'triangle actual and collapsed strong skeleton counts')
    block_group = tuple(tuple(3 * p[i] + x for i in range(3) for x in range(3)) + (9 + s[0], 9 + s[1]) for p in permutations(range(3)) for s in permutations(range(2)))
    rotations = tuple(tuple(3 * i + (x + shifts[i]) % 3 for i in range(3) for x in range(3)) + (9, 10) for shifts in product(range(3), repeat=3))
    pending, reps = set(skeletons), []
    while pending:
        key = min(pending)
        orbit = {skeleton_key(relabel(skeletons[key], p)) for p in block_group}
        require(orbit <= pending and len(orbit) == 12, 'triangle skeleton free block orbit partition')
        pending -= orbit
        reps.append(key)
    require(len(reps) == 65, 'triangle skeleton classes')
    hist, split, witnesses, classes = Counter(), Counter(), {}, []
    for key in reps:
        universe = {}
        for signs in product((-1, 1), repeat=3):
            negative = signs[0] * signs[1] * signs[2] == -1
            for phases in product(range(3), repeat=3):
                T = triangle_presentation(*key, signs, phases)
                require(T not in universe, 'triangle sign phase presentation injection')
                universe[T] = negative
        require(len(universe) == 216, 'all sign phase presentations per skeleton')
        pending_T = set(universe)
        local = Counter()
        while pending_T:
            T = min(pending_T)
            orbit = {relabel(T, p) for p in rotations}
            require(orbit <= pending_T, 'triangle literal internal rotation quotient')
            pending_T -= orbit
            a, e, sizes = action_error(T)
            negative = universe[T]
            require(a == (1 if negative else 3), 'triangle unrestricted full automorphism')
            require(a * len(orbit) == 27, 'triangle internal orbit stabilizer')
            require(e == (18 if negative else 0), 'triangle full action error')
            validate(T)
            local[e] += 1
            hist[e] += 1
            split[(strong(key[0]), e)] += 1
            classes.append({'skeleton': key, 'rows': T, 'automorphisms': a, 'error': e, 'vertex_orbit_sizes': sizes})
            witnesses.setdefault((strong(key[0]), e), T)
        require(local == Counter({0: 12, 18: 4}), 'four rigid twelve balanced per skeleton')
    require(hist == Counter({0: 780, 18: 260}), 'triangle complete class error distribution')
    require(split[(False,18)] == 196 and split[(True,18)] == 64, 'triangle majority split')
    return hist, {'classification': [{'cyclic_majority': c, 'strong': s, 'stable_sizes': p, 'skeletons': k} for (c,s,p),k in sorted(classification.items())],
       'skeleton_representatives': reps, 'collapsed_strong_controls': collapsed, 'classes': classes,
       'class_histogram': [{'cyclic_majority': c, 'error': e, 'classes': k} for (c,e),k in sorted(split.items())],
       'total_error': 4680, 'raw_actual_presentations': 780 * 216, 'negative_raw_presentations': 780 * 108}

def run():
    parts = sorted(odd_partitions(11))
    expected = [(1,)*11, (1,)*8+(3,), (1,)*6+(5,), (1,)*5+(3,3), (1,)*4+(7,), (1,)*3+(3,5),
                (1,1,3,3,3), (1,1,9), (1,3,7), (1,5,5), (3,3,5), (11,)]
    require(parts == sorted(expected), 'exhaustive twelve odd partitions')
    sevens, seven_report = regular_sevens()
    nines, nine_report = regular_nines()
    uniform_hist, uniform_report = check_uniform_sectors(sevens, nines)
    five_hist, five_report = check_five_sector()
    triangle_hist, triangle_report = check_triangles()
    residual = uniform_hist + five_hist + triangle_hist
    residual.pop(0, None)
    require(residual == Counter({14:24, 18:263, 36:2, 40:88, 72:7}), 'nonregular strong eleven full positive distribution')
    total = sum(e * count for e, count in residual.items())
    require(total == 9166 and sum(residual.values()) == 384, 'nonregular strong eleven count and first moment')
    low_a = (1,1,1,2,4)
    a2 = [sum(low_a[i]*low_a[k-i] for i in range(k+1)) for k in range(5)]
    require(a2 == [1,2,3,6,13], 'low SCC padding coefficients')
    # These two already certified lower-order histograms are explicit inputs.
    inherited_s9 = Counter({14:1,18:7,36:2,72:7})
    inherited_s10 = Counter({14:4,18:32,40:84})
    padded = residual.copy()
    for e,count in inherited_s10.items(): padded[e] += a2[1]*count
    for e,count in inherited_s9.items(): padded[e] += a2[2]*count
    padded[14] += a2[4]
    require(padded == Counter({14:48,18:348,36:8,40:256,72:28}), 'known total eleven positive distribution')
    require(sum(padded.values()) == 688 and sum(e*c for e,c in padded.items()) == 19480, 'known total eleven support and error')
    require(19480 % 11 == 10, 'prospective upper trace residue')
    require(19480 % 22 == 10, 'total error residue modulo twenty two')
    return {'status':'PASS', 'scope':'complete nonregular strong order-eleven correction; regular eleven not computed by this packet',
      'odd_partitions':parts, 'regular_seven':seven_report, 'inherited_regular_nine':nine_report,
      'uniform_substitutions':uniform_report, 'five_cells':five_report, 'triangle_cells':triangle_report,
      'nonregular_strong_positive_histogram':sorted(residual.items()), 'nonregular_strong_error':total,
      'known_total_positive_histogram':sorted(padded.items()), 'known_total_error':19480,
      'A_squared_coefficients_through_four':a2, 'upper_trace_modulo_eleven':10,
      'unknown':{'regular_eleven_error':'11*(P_almost_regular_10 - a_regular_11)',
       'total_error_order_eleven':'regular_eleven_error + 19480',
       'positive_regular_eleven_error_support':[22,44,66,88,110],
       'total_error_modulo_twenty_two':10,
       'order_eleven_count':'(U_11 - regular_eleven_error - 19480)/11'},
      'boundary':{'full_order_eleven_enumeration':False, 'order_eleven_atlas_or_count_input':False,
       'regular_eleven_classification':False, 'quotient_and_exceptional_family_enumeration':True,
       'eleven_vertex_automorphism_tests':'only explicitly generated exceptional-family representatives'},
      'check_counts':dict(sorted(CHECKS.items()))}

if __name__ == '__main__':
    print(json.dumps(run(), indent=2, sort_keys=True))
