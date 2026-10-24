#!/usr/bin/env python3
"""Native fixed-path relocation/holonomy census, Python standard library only."""
from argparse import ArgumentParser
from collections import Counter, defaultdict, deque
from fractions import Fraction
from itertools import permutations
from math import factorial, gcd
from pathlib import Path
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p):
    q = [0] * len(p)
    for i, x in enumerate(p):
        q[x] = i
    return tuple(q)


def cycle_type(p):
    seen = set()
    sizes = []
    for i in range(len(p)):
        if i not in seen:
            j, size = i, 0
            while j not in seen:
                seen.add(j)
                size += 1
                j = p[j]
            sizes.append(size)
    return tuple(sorted(sizes))


def group_closure(generators, n):
    identity = tuple(range(n))
    group = {identity}
    todo = [identity]
    while todo:
        x = todo.pop()
        for g in generators:
            z = compose(x, g)
            if z not in group:
                group.add(z)
                todo.append(z)
    return group


def rows_from_word(n, word):
    rows = [0] * n
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            if j == i + 1 or (word >> bit) & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
            if j != i + 1:
                bit += 1
    return tuple(rows)


def rows_from_label_code(n, code):
    rows = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            u, v = (i, j) if (code >> k) & 1 else (j, i)
            rows[u] |= 1 << v
            k += 1
    return tuple(rows)


def label_code(rows, p):
    k, code = 0, 0
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            code |= ((rows[p[i]] >> p[j]) & 1) << k
            k += 1
    return code


def word_in_order(rows, p):
    code, k = 0, 0
    for i in range(len(p)):
        for j in range(i + 2, len(p)):
            code |= ((rows[p[i]] >> p[j]) & 1) << k
            k += 1
    return code


def is_path(rows, p):
    return all((rows[u] >> v) & 1 for u, v in zip(p, p[1:]))


def is_auto(rows, p):
    return all(((rows[i] >> j) & 1) == ((rows[p[i]] >> p[j]) & 1)
               for i in range(len(p)) for j in range(i + 1, len(p)))


def hamiltonian_count(rows):
    n = len(rows)
    d = [defaultdict(int) for _ in range(1 << n)]
    for v in range(n):
        d[1 << v][v] = 1
    for mask in range(1 << n):
        for v, count in tuple(d[mask].items()):
            remaining = rows[v] & ~mask
            while remaining:
                bit = remaining & -remaining
                remaining ^= bit
                w = bit.bit_length() - 1
                d[mask | bit][w] += count
    return sum(d[-1].values())


def relocations(n, max_span=None):
    identity = tuple(range(n))
    result = []
    for i in range(n):
        q = identity[:i] + identity[i + 1:]
        for j in range(n):
            if i == j or (max_span is not None and abs(i-j)+1 > max_span):
                continue
            p = q[:j] + (i,) + q[j:]
            result.append(p)
    # The two descriptions of an adjacent swap coincide; retain one arrow.
    return sorted(set(result))


def scan_components(n, max_span=None, independent=True):
    """No canonicalization or precomputed isomorphism data enters generation."""
    nwords = 1 << ((n - 1) * (n - 2) // 2)
    moves = relocations(n, max_span)
    component_of = [-1] * nwords
    components = []
    identity = tuple(range(n))
    permutations_n = list(permutations(range(n))) if independent else None
    total_edges = 0
    z = defaultdict(Fraction)
    for root in range(nwords):
        if component_of[root] >= 0:
            continue
        index = len(components)
        root_rows = rows_from_word(n, root)
        transport = {root: identity}
        transport_inv = {root: identity}
        component_of[root] = index
        todo = deque([root])
        group = {identity}
        generators = []
        edge_count = 0
        while todo:
            t = todo.popleft()
            rows = rows_from_word(n, t)
            for p in moves:
                if not is_path(rows, p):
                    continue
                u = word_in_order(rows, p)
                edge_count += 1
                candidate = compose(transport[t], p)
                if component_of[u] < 0:
                    component_of[u] = index
                    transport[u] = candidate
                    transport_inv[u] = inverse(candidate)
                    todo.append(u)
                else:
                    require(component_of[u] == index, 'reverse-closed edge crossed finished component')
                    g = compose(candidate, transport_inv[u])
                    if g not in group:
                        require(is_auto(root_rows, g), 'holonomy is not a source automorphism')
                        generators.append(g)
                        group = group_closure(generators, n)
        H = hamiltonian_count(root_rows)
        b, g = len(transport), len(group)
        require(b * g <= H, 'holonomy saturation inequality failed')
        if independent:
            autos = {p for p in permutations_n if is_auto(root_rows, p)}
            require(group <= autos, 'holonomy subgroup exceeds automorphisms')
            paths = [p for p in permutations_n if is_path(root_rows, p)]
            require(len(paths) == H, 'subset DP and literal path count disagree')
            all_words = {word_in_order(root_rows, p) for p in paths}
            require(H == len(autos) * len(all_words), 'independent fibre theorem failed')
            covered = {compose(gamma, r) for gamma in group for r in transport.values()}
            require(len(covered) == b*g and covered <= set(paths), 'actual lift certificate failed')
            require(set(transport) <= all_words, 'component contains wrong source class')
            require((b * g == H) == (set(transport) == all_words and group == autos),
                    'saturation iff failed')
            if b * g == H:
                require(all(component_of[x] == index for x in all_words), 'certified class split')
        else:
            autos = None
        cycles = Counter(cycle_type(p) for p in group)
        for typ, multiplicity in cycles.items():
            z[typ] += Fraction(multiplicity, g)
        components.append({'root_word': root, 'rows': list(root_rows),
                           'word_count': b, 'holonomy_order': g, 'H': H,
                           'directed_relocation_edges': edge_count,
                           'saturated': b*g == H,
                           'holonomy_generators': [list(x) for x in generators],
                           'cycle_histogram': {','.join(map(str,k)): v for k,v in sorted(cycles.items())}})
        total_edges += edge_count
    return components, component_of, z, total_edges


def partitions(n, minimum=1):
    if n == 0:
        yield ()
    for k in range(minimum, n + 1):
        for rest in partitions(n-k, k):
            yield (k,) + rest


def classical_cycle_index(n):
    result = {}
    for sizes in partitions(n):
        if any(x % 2 == 0 for x in sizes):
            continue
        exponent = sum((x-1)//2 for x in sizes)
        exponent += sum(gcd(x,y) for i,x in enumerate(sizes) for y in sizes[i+1:])
        denom = 1
        for length, count in Counter(sizes).items():
            denom *= length ** count * factorial(count)
        result[sizes] = Fraction(1 << exponent, denom)
    return result


def frac(x):
    x = Fraction(x)
    return str(x.numerator) if x.denominator == 1 else f'{x.numerator}/{x.denominator}'


def independent_labelled_control(n):
    """Full label cube checks cycle-index coefficients, not an input to generation."""
    perms = list(permutations(range(n)))
    terms = defaultdict(Fraction)
    total_paths = 0
    label_masses = set()
    canonical_classes = set()
    for code in range(1 << (n*(n-1)//2)):
        rows = rows_from_label_code(n, code)
        auto = [p for p in perms if is_auto(rows,p)]
        H = sum(is_path(rows,p) for p in perms)
        total_paths += H
        canonical_classes.add(min(label_code(rows,p) for p in perms))
        for p in auto:
            terms[cycle_type(p)] += Fraction(1, factorial(n))
    require(dict(terms) == classical_cycle_index(n), 'full labelled fixed-point coefficients mismatch')
    require(total_paths == factorial(n) * (1 << ((n-1)*(n-2)//2)), 'labelled marked mass mismatch')
    return {'hosts': 1 << (n*(n-1)//2), 'marked_records': total_paths,
            'canonical_classes': len(canonical_classes),
            'cycle_index': {','.join(map(str,k)): frac(v) for k,v in sorted(terms.items())}}


def native_false_controls():
    # The C3 source has one normalized word, but three distinct observers.
    cyclic_rows = rows_from_word(3, 0)
    require(cyclic_rows == (2,4,1), 'C3 convention changed')
    adjacent, _, _, _ = scan_components(3, max_span=2)
    cycle = next(c for c in adjacent if c['root_word'] == 0)
    require(cycle['word_count'] == 1 and cycle['holonomy_order'] == 1 and cycle['H'] == 3,
            'adjacent relocation hostile failed')
    # Suppressing owner labels and silently declaring every loop identity is wrong.
    full, _, _, _ = scan_components(3)
    cycle_full = next(c for c in full if c['root_word'] == 0)
    require(cycle_full['holonomy_order'] == 3, 'C3 owner-holonomy missing')
    return {'adjacent_only_C3': cycle, 'full_relocation_C3': cycle_full,
            'marked_observer_stabilizer_order': 1,
            'ordinary_source_automorphism_order': 3}


def blowup_c3(k):
    n = 3*k
    rows = [0]*n
    for u in range(n):
        for v in range(u+1,n):
            bu,bv = u//k,v//k
            forward = (bu == bv) or ((bv-bu)%3 == 1)
            a,b = (u,v) if forward else (v,u)
            rows[a] |= 1 << b
    return tuple(rows)


def window_orders(n, width):
    identity = tuple(range(n))
    orders = set()
    for a in range(n):
        for b in range(a+1,min(n,a+width)+1):
            for middle in permutations(range(a,b)):
                q = identity[:a] + middle + identity[b:]
                if q != identity:
                    orders.add(q)
    return sorted(orders)


def block_word_count(k):
    # Independent inherited run-multinomial formula for C3[TT_k].
    total = 0
    words = 0
    numerator = factorial(k)**3
    def visit(prefix, counts):
        nonlocal total, words
        if len(prefix) == 3*k:
            denominator = 1
            run = 1
            for i in range(1,len(prefix)):
                if prefix[i] == prefix[i-1]:
                    run += 1
                else:
                    denominator *= factorial(run)
                    run = 1
            denominator *= factorial(run)
            require(numerator % denominator == 0, 'run weight is not integral')
            total += numerator//denominator
            words += 1
            return
        for v in range(3):
            if counts[v] < k and (not prefix or v == prefix[-1] or v == (prefix[-1]+1)%3):
                updated = list(counts)
                updated[v] += 1
                visit(prefix+(v,),tuple(updated))
    visit((),(0,0,0))
    return total, words


def locality_controls():
    result = {'backward_nonpath_sources':[], 'cyclic_chain_inflations':[]}
    for n in range(3,9):
        rows = rows_from_word(n,0)
        orders = window_orders(n,n-1)
        require(not any(is_path(rows,q) for q in orders), 'proper-window rigidity failed')
        rotations = [tuple(range(i,n))+tuple(range(i)) for i in range(n)]
        require(all(is_path(rows,q) for q in rotations), 'global rotations missing')
        autos = [p for p in permutations(range(n)) if is_auto(rows,p)]
        require(len(autos) == (3 if n == 3 else 1), 'backward-nonpath automorphism order')
        result['backward_nonpath_sources'].append({'n':n, 'rows':list(rows),
            'proper_window_orders_tested':len(orders), 'accepted_nonidentity_orders':0,
            'actual_global_rotations':n, 'H':hamiltonian_count(rows), 'aut_order':len(autos)})
    for k in range(1,5):
        n,width = 3*k,k+1
        rows = blowup_c3(k)
        q = tuple((i+k)%n for i in range(n))
        require(is_auto(rows,q), 'block rotation is not an automorphism')
        require(compose(compose(q,q),q) == tuple(range(n)), 'block rotation order')
        orders = window_orders(n,width)
        require(not any(is_path(rows,p) for p in orders), 'inflation window frozen test failed')
        H = hamiltonian_count(rows)
        run_count, block_words = block_word_count(k)
        require(run_count == H, 'block-word sum and Hamiltonian DP disagree')
        result['cyclic_chain_inflations'].append({'k':k,'n':n,'window_width':width,
            'rows':list(rows),'nonidentity_window_orders_tested':len(orders),
            'accepted_nonidentity_orders':0,'certified_source_auto_order_at_least':3,
            'block_rotation':list(q),'H':H,'valid_block_words':block_words,
            'independent_run_multinomial_H':run_count,
            'missing_inertia_mass_at_frozen_word_at_least':frac(Fraction(2,H))})
    return result


def rescue_certificate(n, root_word=0):
    rows = rows_from_word(n,root_word)
    H = hamiltonian_count(rows)
    identity = tuple(range(n))
    transport = {root_word:identity}
    group = {identity}
    generators = []
    steps = [{'word_count':1,'holonomy_order':1,'covered_paths':1}]
    chosen = []
    for q in permutations(range(n)):
        if len(transport)*len(group) == H:
            break
        if not is_path(rows,q):
            continue
        y = word_in_order(rows,q)
        changed = False
        if y not in transport:
            transport[y] = q
            changed = True
        else:
            g = compose(q,inverse(transport[y]))
            if g not in group:
                require(is_auto(rows,g),'rescue loop has invalid source action')
                generators.append(g)
                group = group_closure(generators,n)
                changed = True
        if changed:
            chosen.append(list(q))
            steps.append({'word_count':len(transport),'holonomy_order':len(group),
                          'covered_paths':len(transport)*len(group)})
            require(steps[-1]['covered_paths'] > steps[-2]['covered_paths'], 'rescue made no progress')
            require(steps[-1]['covered_paths'] <= H, 'rescue exceeds true paths')
    require(len(transport)*len(group) == H,'full reobservation rescue did not saturate')
    return {'n':n,'root_word':root_word,'H':H,'chosen_global_observers':chosen,'steps':steps}


def partial_graph_bounds():
    output = []
    true_counts = {3:2,4:4,5:12,6:56}
    for n in range(3,7):
        comps, _, _, edges = scan_components(n,max_span=3,independent=False)
        lower = sum(Fraction(c['word_count']*c['holonomy_order'],c['H']) for c in comps)
        require(lower <= true_counts[n] <= len(comps),'partial-graph count sandwich failed')
        output.append({'n':n,'max_relocation_span':3,'components':len(comps),
                       'inertia_lower_bound':frac(lower),'true_count_control':true_counts[n],
                       'unsaturated_components':sum(not c['saturated'] for c in comps),
                       'directed_edges':edges})
    return output


def main():
    parser = ArgumentParser()
    parser.add_argument('--output', default='normal.json')
    parser.add_argument('--max-n', type=int, default=7)
    args = parser.parse_args()
    require(1 <= args.max_n <= 7, 'this reproducibility packet is bounded at n=7')
    report = {'schema': 'session50-native-observer-holonomy-v1', 'orders': [],
              'full_labelled_controls': {}, 'hostile_controls': native_false_controls(),
              'locality_controls':locality_controls(),
              'global_rescue_certificates':[rescue_certificate(n) for n in range(3,7)],
              'partial_graph_bounds':partial_graph_bounds()}
    known = {1:1,2:1,3:2,4:4,5:12,6:56,7:456}
    for n in range(1,args.max_n+1):
        comps, lookup, cycle, edges = scan_components(n)
        require(all(c['saturated'] for c in comps), f'local relocation failed saturation at n={n}')
        require(len(comps) == known[n], f'known small count failed at n={n}')
        if n >= 3:
            require(edges == Fraction(len(lookup)) * Fraction(n*n+n-4,8),
                    'all-size relocation-arrow formula mismatch')
        expected_cycle = classical_cycle_index(n)
        require(dict(cycle) == expected_cycle, f'classical cycle-index coefficient mismatch at n={n}')
        require(sum(cycle.values()) == len(comps), 'cycle-index count mismatch')
        base = Fraction(1 << (n*(n-1)//2), factorial(n))
        marked_mass = sum(c['word_count'] for c in comps)
        h_sum = sum(c['H'] for c in comps)
        identity_tiling_mass = sum(Fraction(c['word_count'], c['H']) for c in comps)
        inertia_tiling_mass = sum(Fraction(c['word_count']*c['holonomy_order'], c['H']) for c in comps)
        require(identity_tiling_mass == base, 'identity-loop mass does not match labelled mass')
        require(inertia_tiling_mass == len(comps), 'inertia tiling mass mismatch')
        order = {'n':n, 'fixed_path_words':marked_mass, 'classes':len(comps),
                 'all_components_saturated':True, 'directed_relocation_edges':edges,
                 'labelled_groupoid_mass':frac(base),
                 'symmetry_correction':frac(len(comps)-base),
                 'inertia_loop_weighted_count':frac(inertia_tiling_mass),
                 'raw_closed_reobservation_count':h_sum,
                 'aut_order_histogram':dict(sorted(Counter(c['holonomy_order'] for c in comps).items())),
                 'cycle_index':{','.join(map(str,k)):frac(v) for k,v in sorted(cycle.items())},
                 'components':comps}
        report['orders'].append(order)
        if n <= 5:
            report['full_labelled_controls'][str(n)] = independent_labelled_control(n)
        print(f'n={n}: words={marked_mass}, classes={len(comps)}, relocation arrows={edges}, saturated={len(comps)}/{len(comps)}', flush=True)
    Path(args.output).write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print('PASS: native holonomy saturation, independent paths/automorphisms, cycle index, hostile controls')

if __name__ == '__main__':
    main()
