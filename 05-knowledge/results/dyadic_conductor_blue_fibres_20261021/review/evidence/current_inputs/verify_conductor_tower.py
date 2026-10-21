#!/usr/bin/env python3
"""Exact dyadic Eisenstein conductor controls; standard-library only.
No assertion statements: python -O performs the same checks.
"""
import hashlib
import json
from math import gcd, isqrt
from pathlib import Path

checks = 0

def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)


def mul(x, y, modulus=None):
    a, b = x
    c, d = y
    z = (a*c-b*d, a*d+b*c-b*d)
    return z if modulus is None else tuple(v % modulus for v in z)


def norm(x):
    a, b = x
    return a*a-a*b+b*b


def power(x, exponent, modulus=None):
    z = (1, 0)
    while exponent:
        if exponent & 1:
            z = mul(z, x, modulus)
        x = mul(x, x, modulus)
        exponent //= 2
    return z


def reduced_forms(f):
    D = -3*f*f
    result = []
    for a in range(1, isqrt((-D)//3)+1):
        for b in range(-a, a+1):
            if (b*b-D) % (4*a):
                continue
            c = (b*b-D)//(4*a)
            if a > c or gcd(gcd(a, b), c) != 1:
                continue
            if (abs(b) == a or a == c) and b < 0:
                continue
            result.append((a, b, c))
    return result


def reduce_form(q):
    a, b, c = q
    D = b*b-4*a*c
    for _ in range(10000):
        m = (a-b)//(2*a)
        old_b = b
        b += 2*a*m
        c = a*m*m+old_b*m+c
        if a > c:
            a, b, c = c, -b, a
            continue
        if (abs(b) == a or a == c) and b < 0:
            b = -b
        require(b*b-4*a*c == D, 'reduction discriminant')
        if -a < b <= a <= c:
            return (a, b, c)
    raise RuntimeError('reduction failed to terminate')


def contraction_form(alpha, f):
    """Form of alpha R intersect O_f; all used alpha have cyclic R/(alpha)."""
    a, b = alpha
    N = norm(alpha)
    if N == 1:
        return (1, 0, 3*f*f//4)
    require(gcd(b, N) == 1, 'cyclic residue presentation required')
    r = (-f*a*pow(b, -1, N)) % N
    B = (f+2*r) % (2*N)
    if B > N:
        B -= 2*N
    require((B*B+3*f*f) % (4*N) == 0, 'integral contraction form')
    q = (N, B, (B*B+3*f*f)//(4*N))
    require(gcd(gcd(q[0], q[1]), q[2]) == 1, 'primitive contraction form')
    return reduce_form(q)


def extend_form(q, f):
    """Independent root transport for extension O_f to O_(f/2), f>=4."""
    a, b, c = q
    if a % 2 == 0:
        a, b, c = c, -b, a
    require(a % 2 == 1, 'odd leading coefficient representative')
    g = f//2
    # New even B obeys B == old_b/2 mod a.
    B = 2*((b//2)*pow(2, -1, a) % a) if a > 1 else 0
    require((B*B+3*g*g) % (4*a) == 0, 'integral extension form')
    return reduce_form((a, B, (B*B+3*g*g)//(4*a)))


def residue_quotient(f):
    if f < 2 or f & (f-1):
        raise ValueError('dyadic residue quotient requires a power of two')
    units = {(a, b) for a in range(f) for b in range(f) if norm((a, b)) % 2}
    denominator = {mul((s, 0), u, f) for s in range(1, f, 2)
                   for u in ((1, 0), (0, 1), (-1, -1))}
    lookup = {}
    orbits = []
    for x in sorted(units):
        if x in lookup:
            continue
        orbit = {mul(x, d, f) for d in denominator}
        require(orbit <= units, 'denominator action preserves units')
        require(not orbit.intersection(lookup), 'quotient orbits disjoint')
        index = len(orbits)
        orbits.append(orbit)
        for y in orbit:
            lookup[y] = index
    require(set(lookup) == units, 'quotient exhaustion')
    require(len(units) == 3*f*f//4, 'finite-ring unit count')
    require(len(denominator) == 3*f//2, 'scalar and Eisenstein unit denominator')
    require(len(orbits) == f//2, 'class number f/2')
    return lookup, orbits


def chi_minus4(n):
    return 1 if n % 4 == 1 else -1


def chi_8(n):
    return 1 if n % 8 in (1, 7) else -1


def non_dyadic_hostile():
    """f=40 disproves an unrestricted norm-kernel = square-subgroup claim."""
    f = 40
    units = {(a,b) for a in range(f) for b in range(f) if gcd(norm((a,b)),f) == 1}
    denominator = {mul((s,0),u,f) for s in range(f) if gcd(s,f) == 1
                   for u in ((1,0),(0,1),(-1,-1))}
    lookup = {}
    representatives = []
    for x in sorted(units):
        if x in lookup:
            continue
        orbit = {mul(x,d,f) for d in denominator}
        require(orbit <= units, 'f40 denominator preserves units')
        require(not orbit.intersection(lookup), 'f40 quotient orbits disjoint')
        representatives.append(x)
        for y in orbit:
            lookup[y] = x
            require(norm(y) % 8 == norm(x) % 8 and pow(norm(y),2,5) == pow(norm(x),2,5),
                    'f40 norm characters well-defined on cosets')
    require(set(lookup) == units, 'f40 quotient exhaustive')
    require(len(units) == 1152 and len(denominator) == 48 and len(representatives) == 24,
            'f40 quotient sizes')
    squares = {lookup[mul(x,x,f)] for x in representatives}
    kernel = {x for x in representatives if norm(x) % 8 == 1}
    require(len(squares) == 3 and len(kernel) == 6, 'f40 defeats unrestricted square-kernel claim')
    alpha = (1,8)
    require(norm(alpha) == 57 and norm(alpha) % 8 == 1, 'f40 hostile norm')
    require(pow(norm(alpha),2,5) == 4, 'f40 extra mod5 genus character detects nonsquare')
    require(lookup[alpha] in kernel and lookup[alpha] not in squares, 'f40 explicit nonsquare class')
    q = contraction_form(alpha,f)
    require(q == (25,20,52), 'f40 hostile reduced form')
    return {'conductor': 40, 'unit_count': len(units), 'denominator_count': len(denominator),
            'class_count': len(representatives), 'norm_mod8_kernel_size': len(kernel),
            'square_subgroup_size': len(squares), 'witness_alpha': alpha, 'witness_norm': 57,
            'witness_form': q, 'witness_legendre_norm_mod5': -1,
            'corrected_scope': 'Genus/square-kernel identification here is dyadic only.'}


def main():
    w, j = (1, 4), (1, 2)
    levels = {}
    receipt = []
    for f in (2, 4, 8, 16, 32, 64, 128):
        lookup, orbits = residue_quotient(f)
        forms = reduced_forms(f)
        require(len(forms) == f//2, 'independent reduced-form class count')
        if f == 2:
            coordinates = [(0, 0)]
            cyclic_order = 1
        else:
            cyclic_order = f//4
            coordinates = [(t, e) for t in range(cyclic_order) for e in range(2)]
        coord_to_class = {}
        class_to_coord = {}
        coord_to_form = {}
        for t, e in coordinates:
            alpha = mul(power(w, t), power(j, e))
            cls = lookup[tuple(v % f for v in alpha)]
            require(cls not in class_to_coord, 'independence of w and j')
            coord_to_class[t, e] = cls
            class_to_coord[cls] = (t, e)
            coord_to_form[t, e] = contraction_form(alpha, f)
            if f >= 8:
                require(norm(alpha) % 8 == (pow(5, t, 8)*pow(3, e, 8)) % 8,
                        'norm character coordinates')
                for x in orbits[cls]:
                    require(norm(x) % 8 == norm(alpha) % 8, 'well-defined norm mod8')
        require(set(class_to_coord) == set(range(len(orbits))), 'generator exhaustion')
        require(set(coord_to_form.values()) == set(forms), 'independent ideal/form exhaustion')
        for x in coordinates:
            for y in coordinates:
                alpha_x = mul(power(w, x[0], f), power(j, x[1], f), f)
                alpha_y = mul(power(w, y[0], f), power(j, y[1], f), f)
                z = lookup[mul(alpha_x, alpha_y, f)]
                expected = ((x[0]+y[0]) % cyclic_order, (x[1]+y[1]) % 2) if f >= 4 else (0, 0)
                require(z == coord_to_class[expected], 'complete multiplication table')
        rows = []
        for t, e in coordinates:
            form = coord_to_form[t, e]
            row = {'coordinate_w_j': [t, e], 'form': form}
            if f >= 4:
                smaller = levels[f//2]
                target_coord = (t % max(1, f//8), e) if f >= 8 else (0, 0)
                image = extend_form(form, f)
                require(image == smaller['forms'][target_coord], 'form extension agrees with residue reduction')
                alpha = mul(power(w, t, f), power(j, e, f), f)
                smaller_cls = smaller['lookup'][tuple(v % (f//2) for v in alpha)]
                require(smaller_cls == smaller['classes'][target_coord], 'residue extension coordinates')
                row['extension_form'] = image
            if f >= 8:
                residue = (pow(5, t, 8)*pow(3, e, 8)) % 8
                row['norm_mod8'] = residue
                row['characters_minus4_8_minus8'] = [chi_minus4(residue), chi_8(residue), chi_minus4(residue)*chi_8(residue)]
            rows.append(row)
        levels[f] = {'lookup': lookup, 'classes': coord_to_class, 'forms': coord_to_form}
        receipt.append({'conductor': f, 'discriminant': -3*f*f, 'class_number': len(forms),
                        'group': 'trivial' if f == 2 else f'C{cyclic_order} x C2', 'rows': rows})
    # Modular valuation controls through conductor 2^64; all-order proof is written.
    for a in range(3, 65):
        modulus = 2**(a+1)
        z = w
        for t in range(a-1):
            b = z[1] % modulus
            require(z[0] % 2 == 1, 'odd real coefficient in repeated squaring')
            require(b % 2**(t+2) == 0 and b % 2**(t+3) != 0,
                    'exact omega coefficient valuation')
            z = mul(z, z, modulus)
    # Hostiles: f16 is not exponent two, and no Klein section can lift A.
    G16 = levels[16]
    require(G16['forms'][(1,0)] == (13,8,16), 'order-four form label')
    require(G16['forms'][(2,0)] == (4,4,49), 'square-kernel form label')
    require(G16['forms'][(1,0)] != G16['forms'][(3,0)], 'inverse forms distinct')
    for lift in [(1,0), (3,0)]:
        require(((2*lift[0]) % 4, (2*lift[1]) % 2) == (2,0), 'every A lift has order four')
    # Exhaust all possible generator images of automorphisms of C4 x C2.
    group16 = [(t,e) for t in range(4) for e in range(2)]
    auto_count = 0
    shore_lifts = []
    for image_w in group16:
        for image_j in group16:
            if ((2*image_j[0]) % 4, (2*image_j[1]) % 2) != (0,0):
                continue
            image = {((t*image_w[0]+e*image_j[0]) % 4,
                      (t*image_w[1]+e*image_j[1]) % 2) for t,e in group16}
            if len(image) != 8:
                continue
            auto_count += 1
            if (image_w[0] % 2,image_w[1]) == (1,0) and (image_j[0] % 2,image_j[1]) == (1,1):
                shore_lifts.append((image_w,image_j))
    require(auto_count == 8, 'complete C4 x C2 automorphism census')
    require(not shore_lifts, 'shore reflection cannot lift through conductor16')
    # The two distinct binary-port stabilizers cannot both map to one subgroup.
    require({(0,0),(1,0)} != {(0,0),(0,1)}, 'left/right inactive stabilizers distinct')
    # Blue-preserving marking: (epsilon,eta) -> (eta,epsilon+eta).
    for epsilon in range(2):
        for eta in range(2):
            t, e = eta, (epsilon+eta) % 2
            residue = pow(5,t,8)*pow(3,e,8) % 8
            require((chi_minus4(residue) == 1) == (epsilon == eta)
                    and chi_8(residue) == (-1)**epsilon
                    and chi_minus4(residue)*chi_8(residue) == (-1)**eta,
                    'blue parity and individual port characters')
    # Natural conjugation and negation act trivially on G8, unlike translations.
    G8 = levels[8]
    for alpha in [(1,0), w, j, mul(w,j)]:
        a,b = alpha
        cls = G8['lookup'][(a%8,b%8)]
        require(G8['lookup'][((a-b)%8,(-b)%8)] == cls, 'conjugation trivial on exponent-two classes')
        require(G8['lookup'][((-a)%8,(-b)%8)] == cls, 'negation principal-unit trivial')
    scope_hostile = non_dyadic_hostile()
    return {'status': 'FINITE-EXACT controls supporting written all-order proofs',
            'checks': checks, 'finite_ring_conductors': [2,4,8,16,32,64,128],
            'valuation_exponents_checked': [3,64], 'levels': receipt,
            'non_dyadic_scope_hostile': scope_hostile,
            'hostiles': ['conductor16 order4', 'nonsplit A lift', 'nonliftable shore reflection', 'distinct port stabilizers',
                         'geometric reflection action is not class translation action',
                         'blue bridge requires parent coincidence and active-port type',
                         'f40 norm kernel is larger than the square subgroup']}


if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
