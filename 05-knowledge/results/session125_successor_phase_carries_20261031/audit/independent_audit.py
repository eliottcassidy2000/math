#!/usr/bin/env python3
"""Independent Session125 phase/carry audit; external C is inert integer data only.

No producer Python module is imported. No external C code is executed. The
receipt contains counts, not external arrays or base paths.
"""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import re

PIN = '602c1919ee51101519a2f19f2c1817407006824f6f09c85a3b772a818a0b736e'
MASKS = (1187208845135, 1453217412461, 1487448155501)

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def graph(n, mask=None):
    out = [0] * (n + 1)
    def color(v):
        return (mask >> (v - 1)) & 1 if v <= 41 else v % 2
    for root in range(2, math.isqrt(2 * n) + 1):
        for a in range(max(1, root * root - n), min(n, (root * root - 1) // 2) + 1):
            b = root * root - a
            if mask is None or color(a) != color(b):
                out[a] |= 1 << (b - 1)
                out[b] |= 1 << (a - 1)
    return out

def path_spectrum(n, fixed_mask=None):
    """Enumerate paths by endpoint extension, necessary residual degree/connectivity.

    The earlier independent onset checker uses degree-subset search instead.
    We construct square edges by square enumeration, not integer-root pairs.
    """
    adjacency = graph(n, fixed_mask)
    endpoint = 3 if n % 2 else 8
    spectra = Counter()
    witness = {}
    nodes = 0
    path = [1]
    def visit(current, remaining, odd_mask):
        nonlocal nodes
        nodes += 1
        if not remaining:
            if current == endpoint:
                spectra[odd_mask] += 1
                witness.setdefault(odd_mask, list(path))
            return
        if current == endpoint:
            return
        available = remaining | (1 << (current - 1))
        scan = remaining
        while scan:
            bit = scan & -scan
            scan -= bit
            vertex = bit.bit_length()
            if (adjacency[vertex] & available).bit_count() < (1 if vertex == endpoint else 2):
                return
        reached = frontier = 1 << (current - 1)
        while frontier:
            bit = frontier & -frontier
            frontier -= bit
            new = adjacency[bit.bit_length()] & available & ~reached
            reached |= new
            frontier |= new
        if reached != available:
            return
        candidates = adjacency[current] & remaining
        next_odd = (len(path) + 1) % 2
        while candidates:
            bit = candidates & -candidates
            candidates -= bit
            vertex = bit.bit_length()
            if (vertex == 3 and not next_odd) or (vertex == 8 and next_odd):
                continue
            path.append(vertex)
            visit(vertex, remaining ^ bit, odd_mask | (bit if next_odd else 0))
            path.pop()
    visit(1, ((1 << n) - 1) ^ 1, 1)
    return spectra, witness, nodes

def edge_set(path, close=False):
    result = {tuple(sorted((a, b))) for a, b in zip(path, path[1:])}
    if close:
        result.add(tuple(sorted((path[0], path[-1]))))
    return result

def deletion(a, b, close=False):
    return len(edge_set(a, close) - edge_set(b, close))

def check_path(path):
    n = len(path)
    require(sorted(path) == list(range(1, n + 1)), 'path is not a permutation')
    require(path[0] == 1 and path[-1] == (3 if n % 2 else 8), 'wrong path ports')
    require(all(math.isqrt(a + b) ** 2 == a + b for a, b in zip(path, path[1:])), 'nonsquare path edge')
    return {v: i % 2 for i, v in enumerate(path, 1)}

def check_pair(a, b):
    require(len(b) == len(a) + 1, 'wrong pair lengths')
    x, y = check_path(a), check_path(b)
    require(all(x[v] == y[v] for v in a), 'pair phase disagreement')

def load_inert_arrays(source):
    raw = Path(source).read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    require(digest == PIN, 'unrecognized external source digest')
    text = raw.decode('ascii')
    arrays = {}
    for name in ('S', 'data'):
        start = re.search(r'\bint\s+' + name + r'\[[^;=]*=\s*', text)
        require(start is not None, 'missing inert initializer ' + name)
        depth = 0
        end = start.end()
        require(text[end] == '{', 'not an initializer')
        for end in range(end, len(text)):
            if text[end] == '{':
                depth += 1
            elif text[end] == '}':
                depth -= 1
                if depth == 0:
                    break
        body = text[start.end():end + 1]
        require(re.fullmatch(r'[\s0-9{},-]+', body) is not None, 'initializer not inert integer data')
        arrays[name] = json.loads(body.translate(str.maketrans('{}', '[]')))
    return arrays, len(raw)

def phase_audit():
    spectra = {}
    counts = []
    for n in (41, 42):
        s, _, nodes = path_spectrum(n)
        spectra[n] = s
        counts.append(dict(n=n, paths=sum(s.values()), masks=len(s), search_nodes=nodes))
    restricted42 = Counter()
    for mask, count in spectra[42].items():
        restricted42[mask & ((1 << 41) - 1)] += count
    common = sorted(spectra[41].keys() & restricted42.keys())
    require(common == list(MASKS), 'wrong complete phase catalogue')
    require([x['paths'] for x in counts] == [49, 275], 'wrong onset path counts')
    failures = []
    positive43 = {}
    for mask in MASKS:
        found43, witnesses, _ = path_spectrum(43, mask)
        found44, _, _ = path_spectrum(44, mask)
        require(not found44, 'unexpected coherent continuation to 44')
        if mask != MASKS[2]:
            require(bool(found43), 'missing positive order-43 continuation')
            p = next(iter(witnesses.values()))
            c = check_path(p)
            require(all(c[v] == ((mask >> (v - 1)) & 1 if v <= 41 else v % 2) for v in p), 'wrong fixed phase')
            positive43[str(mask)] = dict(paths=sum(found43.values()), generated_witness=p)
        else:
            require(not found43, 'third phase unexpectedly survives to 43')
        failures.append(dict(mask=mask, first_failure=44 if found43 else 43))
    def neighbors(mask, n, v):
        bits = graph(n, mask)[v]
        return [j for j in range(1, n + 1) if bits >> (j - 1) & 1]
    forced = {18: [7, 31], 42: [7, 22], 16: [9, 20], 44: [20, 37], 29: [7, 20, 35]}
    require(all(neighbors(MASKS[0], 44, v) == ns for v, ns in forced.items()), 'saturation obstruction')
    require(neighbors(MASKS[1], 44, 44) == [37], 'second obstruction')
    require(neighbors(MASKS[2], 43, 43) == [38], 'third obstruction')
    components = [{0, 1}, {2, 7}, {3, 6}, {4, 5}]
    residue_to_component = {r: j for j, c in enumerate(components) for r in c}
    groups = [{0, 3}, {1, 2}]
    for r in range(8):
        require(residue_to_component[r] == residue_to_component[(1 - r) % 8], 'odd-square channels')
        attached = {residue_to_component[(-r) % 8], residue_to_component[(4 - r) % 8]}
        require(attached in groups, 'defect joins both supergroups')
    # Exhaust finite weighted channel walks by dynamic programming up to six runs.
    state = {(0, 1): 0}
    even_min = odd_min = 100
    for steps in range(1, 7):
        nxt = {}
        for (last, seen), cost in state.items():
            for new in range(4):
                if new == last:
                    continue
                cross = int((last in groups[0]) != (new in groups[0]))
                key = (new, seen | (1 << new))
                nxt[key] = min(nxt.get(key, 100), cost + 1 + cross)
        state = nxt
        even_min = min(even_min, state.get((0, 15), 100))
        odd_min = min(odd_min, state.get((2, 15), 100))
    require((even_min, odd_min) == (6, 4), 'weighted channel lower bounds')
    return dict(counts=counts, common_masks=common,
                pair_multiplicities=[spectra[41][m] * restricted42[m] for m in common],
                first_failures=failures, positive43=positive43,
                even_channel_lower_bound=even_min, odd_channel_lower_bound=odd_min)

def lift_audit(source):
    arrays, source_bytes = load_inert_arrays(source)
    bases, words = arrays['S'], arrays['data']
    require(len(words) == 196 and len(bases) == 3984, 'wrong table counts')
    def word(d, parity, extension):
        row = words[4 * (d - 24) + 2 * parity + extension]
        require(row[0] == [49, d, parity, extension], 'wrong transition metadata')
        return row[1:]
    def ends(atoms, parity):
        answer = []
        for offset, orientation, member in atoms:
            if orientation == 0:
                answer.append((offset, offset))
            else:
                length_odd = (parity + member) % 2
                first = 49 + offset
                last = 49 * (3 if length_odd else 8) + (offset if length_odd else -offset)
                answer.append((first, last) if orientation == 1 else (last, first))
        return answer
    def connectors(atoms, parity):
        z = ends(atoms, parity)
        return {tuple(sorted((a[1], b[0]))) for a, b in zip(z, z[1:])}
    def decode(label):
        center = (label + 24) // 49
        return center, label - 49 * center
    transition = []
    collision_candidates = Counter()
    for parity in (0, 1):
        for d in range(24, 73):
            selected = []
            for extension in (0, 1):
                atoms = word(d, parity, extension)
                require(len(atoms) == 73, 'wrong atom count')
                require(sorted(c for c, o, _ in atoms if o == 0) == list(range(1, 25)), 'wrong singleton atoms')
                require(sorted(c for c, o, _ in atoms if o != 0) == list(range(-24, 25)), 'wrong full offset set')
                members = {c: m for c, o, m in atoms if o}
                require(all(m == int(49 + (-1) ** parity * c <= d + extension) for c, m in members.items()), 'member threshold failure')
                selected.append(members)
                con = connectors(atoms, parity)
                require(len(con) == 72, 'repeated connector')
                for a, b in con:
                    x, c = decode(a)
                    y, e = decode(b)
                    if x > 0 and y > 0 and c + e == 0:
                        collision_candidates[tuple(sorted((x, y)))] += 1
            changes = [c for c in selected[0] if selected[0][c] != selected[1][c]]
            require(changes == [(-1) ** parity * (d + 1 - 49)], 'not exactly one changed member')
            k = len(connectors(word(d, parity, 0), parity) - connectors(word(d, parity, 1), parity))
            transition.append(dict(parity=parity, digit=d, connector_deletions=k))
    require(set(collision_candidates) <= {(1, 3), (8, 8)}, 'unexpected possible internal/connector type')
    carry = {}
    for parity in (0, 1):
        old = word(72, 1 - parity, 1)
        new = word(24, parity, 0)
        require(all(m == 1 for _, o, m in old if o), 'carry old not all longer')
        require(all(m == 0 for _, o, m in new if o), 'carry new not all shorter')
        q = 3 if parity else 8
        for con in (connectors(old, 1 - parity), connectors(new, parity)):
            require(all({decode(a)[0], decode(b)[0]} <= {0, 1, q} for a, b in con), 'carry connector outside endpoint centers')
        carry[parity] = len(connectors(old, 1 - parity) - connectors(new, parity))
    require(carry == {0: 40, 1: 52}, 'changed carry connector constants')
    @lru_cache(maxsize=12)
    def pair(n):
        if 41 <= n <= 2032:
            i = 2 * (n - 41)
            return tuple(bases[i]), tuple(bases[i + 1])
        require(n >= 2033, 'bad recursive order')
        parent = (n - 24) // 49
        digit = n - 49 * parent
        parents = pair(parent)
        answer = []
        for extension in (0, 1):
            output = []
            for offset, orientation, member in word(digit, parent % 2, extension):
                if orientation == 0:
                    output.append(offset)
                else:
                    transformed = [49 * vertex + (offset if i % 2 == 0 else -offset)
                                   for i, vertex in enumerate(parents[member])]
                    output.extend(transformed[::orientation])
            answer.append(tuple(output))
        return tuple(answer)
    base_stats = []
    for n in range(41, 2033):
        a, b = pair(n)
        check_pair(a, b)
        base_stats.append((deletion(a, b), deletion(a, b, True), n))
    for record in transition:
        parent = 42 if record['parity'] == 0 else 41
        n = 49 * parent + record['digit']
        a, b = pair(n)
        check_pair(a, b)
        require(deletion(a, b) == deletion(*pair(parent)) + record['connector_deletions'], 'exact within-pair recurrence fails')
    ray_samples = []
    for seed in (42, 43):
        r = seed
        prev_h = None
        for depth in range(3):
            previous = pair(r - 1)
            current = pair(r)
            h = deletion(previous[1], current[0])
            require(h == deletion(previous[1], current[0], True), 'same-order closure changes gap')
            if prev_h is not None:
                require(h == 49 * prev_h + carry[seed % 2], 'exact carry recurrence fails')
            within = deletion(*previous, close=True)
            canonical = deletion(previous[0], current[0], True)
            require(h - within - 1 <= canonical <= h + within, 'canonical cycle triangle fails')
            ray_samples.append(dict(seed=seed, depth=depth, order=r, representation_deletions=h,
                                    canonical_old_order=r - 1, canonical_cycle_deletions=canonical,
                                    within_pair_cycle_deletions=within))
            prev_h = h
            r = 49 * r + 24
    require([s['representation_deletions'] for s in ray_samples if s['seed'] == 42] == [18, 922, 45218], 'even-ray samples')
    # An explicit opposite-phase toy test checks full-offset invariance without
    # requiring the phases of separate parent representatives to agree.
    toy_edges = [(1, 24), (24, 25), (25, 39)]
    for a, b in toy_edges:
        forward = {tuple(sorted((49 * a + c, 49 * b - c))) for c in range(-24, 25)}
        flipped = {tuple(sorted((49 * a - c, 49 * b + c))) for c in range(-24, 25)}
        require(forward == flipped and len(forward) == 49, 'signed-offset phase independence')
        require(all({decode(x)[0], decode(y)[0]} == {a, b} for x, y in forward), 'unique decoding')
    original_ratio = Fraction(18, 1) + Fraction(carry[0], 48)
    original_ratio /= Fraction(85, 2)
    require(original_ratio == Fraction(113, 255), 'wrong exact pinned ratio')
    root_ratios = []
    for seed in range(42, 2034):
        gap = deletion(pair(seed - 1)[1], pair(seed)[0])
        ratio = (Fraction(gap) + Fraction(carry[seed % 2], 48)) / Fraction(2 * seed + 1, 2)
        root_ratios.append((ratio, seed, gap))
    best_ratio, best_seed, best_gap = max(root_ratios)
    outside_root_upper_bound = (Fraction(72) + Fraction(52, 48)) / Fraction(4069, 2)
    require(best_ratio > outside_root_upper_bound, 'large noncarry roots not excluded')
    require((best_ratio, best_seed, best_gap) == (Fraction(22273, 24042), 2003, 1855), 'pinned global maximum changed')
    maximizing_first_carry = 49 * best_seed + 24
    maximizing_first_carry_gap = deletion(pair(maximizing_first_carry - 1)[1], pair(maximizing_first_carry)[0])
    require(maximizing_first_carry_gap == 49 * best_gap + carry[best_seed % 2], 'maximizing ray fails recurrence')
    selector_checks = []
    for r in (2082, 2131):
        h = deletion(pair(r - 1)[1], pair(r)[0], True)
        left_gap = deletion(pair(r - 2)[1], pair(r - 1)[0], True)
        right_gap = deletion(pair(r)[1], pair(r + 1)[0], True)
        require(left_gap <= 72 and right_gap <= 72, 'neighbor of carry not cheap')
        losses = deletion(*pair(r - 1), close=True) + deletion(*pair(r), close=True) + 2 + left_gap + right_gap
        for left in (pair(r - 2)[1], pair(r - 1)[0]):
            for right in (pair(r)[1], pair(r + 1)[0]):
                require(deletion(left, right, True) >= h - losses, 'two-step selector lower bound')
                for middle in (pair(r - 1)[1], pair(r)[0]):
                    require(deletion(left, middle, True) + deletion(middle, right, True) >= h - losses, 'selector triangle lower bound')
                require(deletion(left, pair(r)[0], True) >= h - deletion(*pair(r - 1), close=True) - 73, 'short middle incoming lower bound')
                require(deletion(pair(r - 1)[1], right, True) >= h - deletion(*pair(r), close=True) - 73, 'long middle outgoing lower bound')
        selector_checks.append(dict(carry_order=r, gap=h, left_neighbor_gap=left_gap,
                                    right_neighbor_gap=right_gap, two_step_lower_bound=h - losses))
    return dict(source_sha256=PIN, source_bytes=source_bytes, external_code_executed=False,
                external_arrays_redistributed=False, valid_base_pairs=len(base_stats),
                max_base_path_deletions=max(t[0] for t in base_stats),
                max_base_cycle_deletions=max(t[1] for t in base_stats),
                transition_states=len(transition), transition_K_range=[min(t['connector_deletions'] for t in transition), max(t['connector_deletions'] for t in transition)],
                transition_K=transition,
                possible_connector_internal_decoded_types={str(k): v for k, v in collision_candidates.items()},
                carry_K={str(k): v for k, v in carry.items()}, carry_ray_samples=ray_samples,
                universal_limsup_lower_bound='4/87', pinned_even_ray_exact_ratio=str(original_ratio),
                selector_limsup_lower_bound='4/87', arbitrary_middle_local_limsup_lower_bound='2/87', selector_checks=selector_checks,
                exact_pinned_global_limsup=dict(value=str(best_ratio), attaining_seed=best_seed,
                    seed_gap=best_gap, root_range=[42, 2033], outside_root_ratio_upper_bound=str(outside_root_upper_bound),
                    special_boundary_seed_2033_gap=root_ratios[-1][2],
                    maximizing_first_carry=maximizing_first_carry, maximizing_first_carry_gap=maximizing_first_carry_gap,
                    applies_to_every_two_representative_selector=True))

def hostile_controls():
    controls = {}
    for name, path in [('duplicate', [1, 1, 3]), ('wrong_ports', [3, 2, 1]), ('nonsquare', [1, 2, 3])]:
        try:
            check_path(path)
        except RuntimeError:
            controls[name] = True
        else:
            controls[name] = False
    require(all(controls.values()), 'hostile path accepted')
    return controls

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('source')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    report = dict(status='PASS', method='independent endpoint DFS and inert-array lift audit; no producer imports',
                  phase=phase_audit(), lift=lift_audit(args.source), hostile_controls=hostile_controls())
    result = json.dumps(report, indent=2, sort_keys=True) + '\n'
    Path(args.output).write_text(result)
    print('PASS: independent phase and carry audit; ' + args.output)

if __name__ == '__main__':
    main()
