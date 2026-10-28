#!/usr/bin/env python3
"""Exact stdlib audit: marked F64/F256, alternating forms and Frobenius.
Integers encode polynomial-basis vectors, least-significant coefficient first.
No floating point, probabilistic tests, external packages, or removable asserts.
"""
import itertools
import json


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def rank(vectors):
    piv = {}
    for v in vectors:
        while v:
            k = v.bit_length() - 1
            if k in piv:
                v ^= piv[k]
            else:
                piv[k] = v
                break
    return len(piv)


def nullspace(rows, n):
    rows = list(rows)
    pivots = []
    r = 0
    for c in range(n):
        p = next((i for i in range(r, len(rows)) if rows[i] >> c & 1), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i] >> c & 1:
                rows[i] ^= rows[r]
        pivots.append(c)
        r += 1
    result = []
    for f in range(n):
        if f not in pivots:
            v = 1 << f
            for i, c in enumerate(pivots):
                if rows[i] >> f & 1:
                    v |= 1 << c
            result.append(v)
    need(all((a & b).bit_count() % 2 == 0 for a in rows for b in result), 'nullspace')
    return result


def span(basis):
    out = [0]
    for b in basis:
        out += [v ^ b for v in out]
    return sorted(out)


class Field:
    def __init__(self, n, polynomial):
        self.n = n
        self.polynomial = polynomial
        self.size = 1 << n

    def mul(self, a, b):
        value = 0
        while b:
            if b & 1:
                value ^= a
            b >>= 1
            a <<= 1
            if a & self.size:
                a ^= self.polynomial
        return value

    def power(self, a, exponent):
        value = 1
        while exponent:
            if exponent & 1:
                value = self.mul(value, a)
            a = self.mul(a, a)
            exponent >>= 1
        return value

    def trace(self, a):
        value = 0
        for _ in range(self.n):
            value ^= a
            a = self.mul(a, a)
        need(value in (0, 1), 'trace is binary')
        return value


def poly_rem(a, b):
    while a and a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def poly_gcd(a, b):
    while b:
        a, b = b, poly_rem(a, b)
    return a


def prime_divisors(n):
    result = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            result.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        result.append(n)
    return result


def audit(n, polynomial):
    f = Field(n, polynomial)
    m = n // 2
    q = 1 << m
    # Rabin irreducibility criterion and complete prime-factor order criterion.
    irred = {str(p): poly_gcd(f.power(2, 1 << (n // p)) ^ 2, polynomial)
             for p in prime_divisors(n)}
    need(f.power(2, 1 << n) == 2 and set(irred.values()) == {1}, 'irreducible')
    order = {str(p): f.power(2, (f.size - 1) // p) for p in prime_divisors(f.size - 1)}
    need(f.power(2, f.size - 1) == 1 and all(v != 1 for v in order.values()), 'primitive')
    pairs = list(itertools.combinations(range(n), 2))

    def wedge(x, y):
        return sum((((x >> i & 1) & (y >> j & 1)) ^
                    ((x >> j & 1) & (y >> i & 1))) << k
                   for k, (i, j) in enumerate(pairs))

    def B(mask, x, y):
        return (mask & wedge(x, y)).bit_count() % 2

    def formrank(mask):
        return rank(sum(B(mask, 1 << i, 1 << j) << j for j in range(n)) for i in range(n))

    frob = [f.mul(1 << i, 1 << i) for i in range(n)]
    constraints = [wedge(frob[i], frob[j]) ^ (1 << k) for k, (i, j) in enumerate(pairs)]
    invariant_basis = nullspace(constraints, len(pairs))
    invariant = span(invariant_basis)
    need(len(invariant_basis) == m, 'invariant dimension')
    nond = [mask for mask in invariant if formrank(mask) == n]
    need(len(nond) == {6: 2, 8: 8}[n], 'nondegenerate invariant count')
    K = [x for x in range(f.size) if f.power(x, q) == x]
    need(len(K) == q, 'subfield size')
    line_sets = sorted({tuple(sorted(f.mul(x, k) for k in K if k)) for x in range(1, f.size)})
    need(len(line_sets) == q + 1 and len(set().union(*map(set, line_sets))) == f.size - 1, 'spread')
    line_bases = []
    line_constraints = []
    for line in line_sets:
        basis = []
        for x in line:
            if rank(basis + [x]) > len(basis):
                basis.append(x)
        need(len(basis) == m, 'line dimension')
        line_bases.append(basis)
        line_constraints.append([wedge(x, y) for x, y in itertools.combinations(basis, 2)])

    def isotropic(mask, index):
        return all((mask & c).bit_count() % 2 == 0 for c in line_constraints[index])

    trace_family = {c: sum(f.trace(f.mul(c, f.mul(1 << i, f.power(1 << j, q)))) << k
                           for k, (i, j) in enumerate(pairs)) for c in K}
    need(len(set(trace_family.values())) == q, 'trace family distinct')
    all_spread = span(nullspace(list(itertools.chain.from_iterable(line_constraints)), len(pairs)))
    need(all_spread == sorted(trace_family.values()), 'full spread classification')
    canonical = trace_family[1]
    need([a for a in nond if all(isotropic(a, i) for i in range(q + 1))] == [canonical], 'unique intersection')
    # Exhaust every subset until first distinguishing cardinality; all minima retained.
    minimum_sets = []
    for count in range(q + 2):
        for inds in itertools.combinations(range(q + 1), count):
            if [a for a in nond if all(isotropic(a, i) for i in inds)] == [canonical]:
                minimum_sets.append(list(inds))
        if minimum_sets:
            break
    need(minimum_sets == [[i] for i in range(1, q + 1)], 'all non-K singleton minima')
    normal = next([f.power(a, 1 << i) for i in range(n)] for a in range(1, f.size)
                  if rank([f.power(a, 1 << i) for i in range(n)]) == n)
    invrows = []
    for mask in invariant:
        coeffs = [B(mask, normal[0], normal[d]) for d in range(1, m + 1)]
        for i in range(n):
            for j in range(n):
                expected = 0 if i == j else coeffs[min((j - i) % n, (i - j) % n) - 1]
                need(B(mask, normal[i], normal[j]) == expected, 'circulant normal form')
        invrows.append({'mask': mask, 'rank': formrank(mask), 'normal_coefficients': coeffs,
                        'isotropic_line_indices': [i for i in range(q + 1) if isotropic(mask, i)]})
    # Explicit hostile pair with identical restrictions on every spread line.
    c = next(c for c in K if c not in (0, 1))
    hostile = trace_family[c]
    witness = next((x, y) for x in range(1, f.size) for y in range(1, f.size)
                   if B(canonical, x, y) != B(hostile, x, y))
    pair_histograms = []
    for i in range(q + 1):
        for j in range(i + 1, q + 1):
            A, C = (0,) + line_sets[i], (0,) + line_sets[j]
            h0 = sum(B(canonical, x, y) for x in A for y in C)
            h1 = sum(B(hostile, x, y) for x in A for y in C)
            need(h0 == h1 == q * (q - 1) // 2, 'cross-line marginal histogram')
            pair_histograms.append(h0)
    # Faithfully record the whole paired table: difference form is nondegenerate.
    all_joint_disagreements = sum(B(canonical, x, y) != B(hostile, x, y)
                                  for x in range(f.size) for y in range(f.size))
    need(all_joint_disagreements == f.size * (f.size - 1) // 2, 'joint disagreement count')
    # Enumerate the ENTIRE centralizer of Frobenius in GL(n,2).
    # In a normal basis F is a cyclic shift. Every commuting map is determined
    # by u=g(normal[0]), and g(F^i normal[0])=F^i u. Thus only 2^n candidates.
    normal_coordinates = {}
    for coordinates in range(f.size):
        value = 0
        for i in range(n):
            if coordinates >> i & 1:
                value ^= normal[i]
        normal_coordinates[value] = coordinates
    centralizer = []
    centralizer_spread = []
    centralizer_form = []
    centralizer_joint = []
    form_only_line_stabilizers = [0] * (q + 1)
    frobenius_columns = {tuple(f.power(1 << i, 1 << r) for i in range(n)) for r in range(n)}
    for u in range(f.size):
        images = [f.power(u, 1 << i) for i in range(n)]
        if rank(images) != n:
            continue
        g = []
        for x in range(f.size):
            value = 0
            for i in range(n):
                if normal_coordinates[x] >> i & 1:
                    value ^= images[i]
            g.append(value)
        columns = tuple(g[1 << i] for i in range(n))
        centralizer.append(columns)
        preserves_form = all(B(canonical, g[1 << i], g[1 << j]) ==
                             B(canonical, 1 << i, 1 << j) for i, j in pairs)
        image_lines = [tuple(sorted(g[x] for x in L)) for L in line_sets]
        preserves_spread = all(L in line_sets for L in image_lines)
        if preserves_form:
            centralizer_form.append(columns)
            for i, L in enumerate(line_sets):
                form_only_line_stabilizers[i] += (image_lines[i] == L)
        if preserves_spread:
            centralizer_spread.append(columns)
        if preserves_form and preserves_spread:
            centralizer_joint.append({'normal_generator_image': u, 'columns': columns,
                                      'line_permutation': [line_sets.index(L) for L in image_lines]})
    need({tuple(g['columns']) for g in centralizer_joint} == frobenius_columns,
         'joint centralizer is exactly Frobenius')
    # A coordinate-free second nondegenerate invariant form in every even n.
    trace_product_mask = sum((f.trace(f.mul(1 << i, 1 << j)) ^
                              (f.trace(1 << i) & f.trace(1 << j))) << k
                             for k, (i, j) in enumerate(pairs))
    need(trace_product_mask in nond, 'trace-product form is invariant nondegenerate')
    need([i for i in range(q + 1) if isotropic(trace_product_mask, i)] == [0],
         'trace-product form vanishes only on K')
    one_port_phase_witnesses = []
    if n == 6:
        for line_index, y in [(2, 3), (4, 5)]:
            x, r = 14, 2
            rotated_y = f.power(y, 1 << r)
            need(rotated_y in line_sets[line_index], 'phase preserves marked line')
            need(B(canonical, x, y) != B(canonical, x, rotated_y), 'one-port phase witness')
            need(B(canonical, f.power(x, 1 << r), rotated_y) == B(canonical, x, y),
                 'simultaneous phase is an isometry')
            one_port_phase_witnesses.append({'line_index': line_index, 'x': x, 'y': y,
                                             'frobenius_power': r, 'rotated_y': rotated_y,
                                             'pairing_before': B(canonical, x, y),
                                             'pairing_after_one_port_rotation': B(canonical, x, rotated_y)})
    # Attach an identification between two K-lines. Scalars preserve lines,
    # Frobenius preserves the full field; neither supplies an absolute phase.
    line_frob = [line_sets.index(tuple(sorted(f.mul(x, x) for x in L))) for L in line_sets]
    orbit_sizes = []
    unseen = set(range(1, f.size))
    while unseen:
        a = min(unseen)
        orb, x = [], a
        while x not in orb:
            orb.append(x)
            x = f.mul(x, x)
        unseen.difference_update(orb)
        orbit_sizes.append(len(orb))
    return {'degree': n, 'q': q, 'polynomial_mask': polynomial,
            'irreducibility_gcd_checks': irred, 'primitive_order_checks': order,
            'upper_triangle_pair_order': pairs, 'normal_basis': normal,
            'frobenius_constraint_rank': rank(constraints), 'invariant_basis_masks': invariant_basis,
            'invariant_form_count': len(invariant), 'nondegenerate_invariant_count': len(nond),
            'all_invariant_forms': invrows, 'trace_family_masks_by_c': trace_family,
            'canonical_mask': canonical, 'spread_line_nonzero_vectors': line_sets,
            'spread_line_bases': line_bases, 'line_frobenius_permutation': line_frob,
            'minimum_marked_lines': len(minimum_sets[0]), 'all_minimum_line_index_sets': minimum_sets,
            'per_line_nondegenerate_survivors': [[a for a in nond if isotropic(a, i)] for i in range(q + 1)],
            'per_line_invariant_solution_dimensions': [len(nullspace(constraints + line_constraints[i], len(pairs))) for i in range(q + 1)],
            'hostile': {'coefficient': c, 'mask': hostile, 'cross_pair_witness': witness,
                        'witness_values': [B(canonical, *witness), B(hostile, *witness)],
                        'each_cross_line_ones': pair_histograms[0], 'each_cross_line_pairs': q * q,
                        'all_ordered_pair_disagreements': all_joint_disagreements},
            'nonzero_frobenius_orbit_sizes': sorted(orbit_sizes),
            'trace_product_form_mask': trace_product_mask,
            'one_port_phase_witnesses': one_port_phase_witnesses,
            'frobenius_centralizer': {
                'invertible_count': len(centralizer),
                'form_preserving_count': len(centralizer_form),
                'spread_preserving_count': len(centralizer_spread),
                'joint_count': len(centralizer_joint),
                'single_line_form_only_stabilizer_orders': form_only_line_stabilizers,
                'joint_maps': centralizer_joint,
                'single_line_stabilizer_orders': [sum(g['line_permutation'][i] == i
                                                    for g in centralizer_joint)
                                                 for i in range(q + 1)]}}


def main():
    result = {'source_pin': '5ed1835235b36a57b8cd11eb705aa259b6bc0f08',
              'convention': 'Polynomial vectors as binary integers; pair masks ordered lexicographically.',
              'cases': [audit(6, 0x43), audit(8, 0x11d)]}
    print(json.dumps(result, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
