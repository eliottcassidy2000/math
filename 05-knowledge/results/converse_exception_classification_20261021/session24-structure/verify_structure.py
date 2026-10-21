#!/usr/bin/env python3
"""Independent exact controls for the Session24 structural proof.

This script imports no producer code. Finite controls are not all-order proofs.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def tournament(n, mask):
    a = [[False] * n for _ in range(n)]
    for bit, (i, j) in enumerate(itertools.combinations(range(n), 2)):
        a[i][j] = bool((mask >> bit) & 1)
        a[j][i] = not a[i][j]
    return a


def cyclic(n, steps):
    return [[i != j and (j-i) % n in steps for j in range(n)] for i in range(n)]


def substitution(a, blocks):
    coords = [(i, j) for i, b in enumerate(blocks) for j in range(len(b))]
    return [[(blocks[i][x][y] if i == j else a[i][j])
             for j, y in coords] for i, x in coords]


def signed_groups(a):
    n = len(a)
    degrees = [sum(row) for row in a]
    pairs = list(itertools.combinations(range(n), 2))
    g, anti = [], []
    for p in itertools.permutations(range(n)):
        if all(degrees[p[i]] == degrees[i] for i in range(n)):
            if all(a[p[i]][p[j]] == a[i][j] for i, j in pairs):
                g.append(p)
        if all(degrees[p[i]] == n-1-degrees[i] for i in range(n)):
            if all(a[p[i]][p[j]] != a[i][j] for i, j in pairs):
                anti.append(p)
    return g, anti


def image(mask, p):
    answer = 0
    for i, j in enumerate(p):
        if mask >> i & 1:
            answer |= 1 << j
    return answer


def regular(mask, group):
    return all(p == tuple(range(len(p))) or image(mask, p) != mask for p in group)


def vertex_orbits(g, n):
    remaining = set(range(n))
    result = []
    while remaining:
        i = min(remaining)
        orbit = {p[i] for p in g}
        result.append(orbit)
        remaining -= orbit
    return result


def palettes(g, anti, n):
    reg = [mask for mask in range(1 << n) if regular(mask, g)]
    representatives = sorted({min(image(s, p) for p in g) for s in reg})
    number = {s: i for i, s in enumerate(representatives)}
    palette = {s: number[min(image(s, p) for p in g)] for s in reg}
    alpha = None
    if anti:
        alpha = [palette[image(s, anti[0])] for s in representatives]
    return reg, representatives, palette, alpha


def orbit_toggle_certificate(g, anti, n):
    orbits = vertex_orbits(g, n)
    s = next(s for s in range(1 << n) if regular(s, g))
    if not anti:
        return s
    tau = next(p for p in anti if all(p[p[i]] == i for i in range(n)))
    pair = next((o, {tau[i] for i in o}) for o in orbits
                if {tau[i] for i in o} != o)
    o, other = pair
    left = sum((s >> i) & 1 for i in o)
    right = sum((s >> i) & 1 for i in other)
    if left == right:
        s ^= sum(1 << i for i in o)
    require(regular(s, g + anti), 'Orbit toggle failed to distinguish signed action')
    return s


def census_small():
    rows = []
    for n in range(2, 6):
        row = dict(n=n, labelled_tournaments=0, intransitive=0,
                   self_converse_intransitive=0, exceptional=0, q_two=0)
        for mask in range(1 << (n*(n-1)//2)):
            a = tournament(n, mask)
            g, anti = signed_groups(a)
            reg, reps, palette, alpha = palettes(g, anti, n)
            k = len(vertex_orbits(g, n))
            require(len(reps) >= 2**k, 'Palette lower bound failed')
            signed_regular = [s for s in reg if regular(s, g + anti)]
            require(not anti or len(anti) == len(g), 'Wrong anti coset size')
            require(len(g) % 2 == 1, 'Even ordinary automorphism group')
            require(bool(reg), 'Ordinary regular subset missing')
            if k > 1:
                row['intransitive'] += 1
                row['self_converse_intransitive'] += bool(anti)
                orbit_toggle_certificate(g, anti, n)
                require(bool(signed_regular), 'Intransitive exception')
            if not signed_regular:
                row['exceptional'] += 1
                require(k == 1 and n % 2 == 1, 'Exceptional not odd/transitive')
            if len(reps) == 2:
                row['q_two'] += 1
                require(k == 1, 'q=2 not transitive')
                require(not anti or not signed_regular, 'Self-converse q=2 not exceptional')
            row['labelled_tournaments'] += 1
        rows.append(row)
    return rows


def details(a):
    n = len(a)
    g, anti = signed_groups(a)
    reg, reps, palette, alpha = palettes(g, anti, n)
    signed_reg = [s for s in reg if regular(s, g+anti)]
    result = dict(n=n, ordinary_order=len(g), anti_count=len(anti),
                  vertex_orbit_sizes=[len(o) for o in vertex_orbits(g, n)],
                  q=len(reps), f=sum(i == j for i, j in enumerate(alpha)) if alpha else None,
                  ordinary_regular_subsets=len(reg), signed_regular_subsets=len(signed_reg))
    if len(result['vertex_orbit_sizes']) > 1:
        result['toggle_certificate_mask'] = orbit_toggle_certificate(g, anti, n)
    return result


def palette_recursion_controls():
    rows = []
    for q in range(1, 11):
        for f in range(q % 2, q+1, 2):
            alpha = list(range(q))
            for i in range(f, q, 2):
                alpha[i], alpha[i+1] = i+1, i
            def rep(t):
                return min(t, t[1:]+t[:1], t[2:]+t[:2])
            reps = {rep(t) for t in itertools.product(range(q), repeat=3)
                    if len(set(t)) > 1}
            fixed = sum(rep((alpha[t[0]], alpha[t[2]], alpha[t[1]])) == t for t in reps)
            require(len(reps) == (q**3-q)//3, 'q recursion fails')
            require(fixed == f*(q-1), 'f recursion fails')
            rows.append(dict(q=q, f=f, q_next=len(reps), f_next=fixed))
    return rows


def actual_triangle_substitution_controls(b):
    """Enumerate actual local subset triples, not just palette color triples."""
    n = len(b)
    g, anti = signed_groups(b)
    reg, reps, palette, alpha = palettes(g, anti, n)
    negative_images = {s: {image(s, p) for p in anti} for s in reg}
    positive_images = {s: {image(s, p) for p in g} for s in reg}
    reflections = [(k, (k-1) % 3, (k-2) % 3) for k in range(3)]
    ordinary_count = signed_count = 0
    for t in itertools.product(reg, repeat=3):
        # Every local action is regular. A top rotation can be compensated
        # exactly when the three local positive orbits all coincide.
        rotation = all(t[(i+1) % 3] in positive_images[t[i]] for i in range(3))
        if rotation:
            continue
        ordinary_count += 1
        has_negative = any(all(t[p[i]] in negative_images[t[i]] for i in range(3))
                           for p in reflections)
        signed_count += not has_negative
    q = len(reps)
    f = sum(i == j for i, j in enumerate(alpha))
    group_order = 3 * len(g)**3
    require(ordinary_count == group_order*(q**3-q)//3, 'Actual ordinary subset count mismatch')
    require(signed_count == group_order*((q**3-q)//3-f*(q-1)), 'Actual signed subset count mismatch')
    return dict(base_order=n, q=q, f=f, product_ordinary_order=group_order,
                actual_local_triples_enumerated=len(reg)**3,
                ordinary_regular_subsets=ordinary_count,
                signed_regular_subsets=signed_count,
                signed_regular_orbits=signed_count//(2*group_order))


def proper_kernel_controls():
    """Four nested cyclic/reflection-invariant kernels on three triangles."""
    vectors = list(itertools.product(range(3), repeat=3))
    kernels = {
        'zero': [(0, 0, 0)],
        'diagonal': [v for v in vectors if v[0] == v[1] == v[2]],
        'sum_zero': [v for v in vectors if sum(v) % 3 == 0],
        'full': vectors,
    }
    rows = []
    for name, kernel in kernels.items():
        g, anti = [], []
        for sign in (1, -1):
            for shift in range(3):
                for k in kernel:
                    p = tuple(3*((sign*i+shift) % 3) + (sign*j+k[(sign*i+shift) % 3]) % 3
                              for i in range(3) for j in range(3))
                    (g if sign == 1 else anti).append(p)
        count = sum(regular(s, g+anti) for s in range(1 << 9))
        require((count == 0) == (name == 'full'), 'Proper base kernel unexpectedly exceptional')
        # Fix c=(singleton,doubleton,doubleton), whose unique anti quotient
        # is reflection fixing block zero, and inspect marker coverage directly.
        all_covered = True
        for x in vectors:
            s = (1 << x[0]) | sum(1 << (3*i+j) for i in (1, 2) for j in range(3) if j != x[i])
            require(regular(s, g), 'Marker subset not ordinary regular')
            all_covered &= any(image(s, p) == s for p in anti)
        require(all_covered == (name == 'full'), 'Marker coverage kernel forcing failed')
        rows.append(dict(kernel=name, kernel_order=len(kernel), signed_group_order=len(g+anti),
                         signed_regular_subsets=count, all_27_markers_covered=all_covered))
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='receipt.json')
    args = parser.parse_args()
    c3 = cyclic(3, {1})
    r5 = cyclic(5, {1, 2})
    p7 = cyclic(7, {1, 2, 4})
    named = {
        'C3': details(c3),
        'R5': details(r5),
        'P7': details(p7),
        'C7_123_nonexceptional': details(cyclic(7, {1, 2, 3})),
        'strong_unequal_C3_C3_1_C3': details(substitution(c3, [c3, [[False]], c3])),
        'stacked_two_C3': details(substitution(tournament(2, 1), [c3, c3])),
        'C3_C3': details(substitution(c3, [c3, c3, c3])),
    }
    require(named['strong_unequal_C3_C3_1_C3']['vertex_orbit_sizes'] == [3, 1, 3],
            'Wrong unequal-module control orbits')
    require(named['C3_C3']['ordinary_order'] == 81, 'Wrong full wreath ordinary group')
    require(named['C3_C3']['signed_regular_subsets'] == 0, 'Triangle tower is not exceptional')
    result = {
        'status': 'FINITE-EXACT controls; not a substitute for all-order proofs',
        'universe': 'All labelled tournaments n=2..5; named controls through n=9; all palette involutions represented by (q,f), q<=10; actual local subset triples for C3,R5,P7; four nested base kernels in F3^3.',
        'census': census_small(),
        'named': named,
        'palette_recursions': palette_recursion_controls(),
        'actual_substitution_counts': [actual_triangle_substitution_controls(b) for b in (c3, r5, p7)],
        'proper_kernel_controls': proper_kernel_controls(),
        'all_order_claim': 'Not inferred from finite controls; see STRUCTURE.md Sections3-7.',
    }
    result['checks'] = CHECKS
    result['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    Path(args.output).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'output': args.output, 'checks': CHECKS,
                      'named': named, 'proper_kernel_controls': result['proper_kernel_controls'],
                      'actual_substitution_counts': result['actual_substitution_counts']}, indent=2))


if __name__ == '__main__':
    main()
