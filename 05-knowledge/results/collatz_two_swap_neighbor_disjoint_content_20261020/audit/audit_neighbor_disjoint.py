#!/usr/bin/env python3
"""Independent standard-library audit of adjacent-source, disjoint two-swaps."""
import argparse
import itertools
import json
import math
from collections import Counter


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def bezout(a, b):
    old_r, r, old_s, s, old_t, t = a, b, 1, 0, 0, 1
    while r:
        k = old_r // r
        old_r, r = r, old_r-k*r
        old_s, s = s, old_s-k*s
        old_t, t = t, old_t-k*t
    require(old_r == 1, ('not coprime', a, b))
    return old_s, old_t


def root(p, n, modulus):
    if modulus == 1:
        return 0
    a, b = bezout(p, n)
    t = 4 * pow(3, -1, modulus) % modulus
    z = pow(2, a, modulus) * pow(t, b, modulus) % modulus
    require(pow(z, p, modulus) == 2 % modulus, ('p clock', p, n))
    require(pow(z, n, modulus) == t, ('n clock', p, n))
    require(math.gcd(z*(z-1), modulus) == 1, ('nonunit root', p, n))
    return z


def carry(word):
    b = 0
    power = 1
    for bit in word:
        if bit:
            b = 3*b + power
        power *= 2
    return b


def carry_by_positions(word):
    positions = [i for i, bit in enumerate(word) if bit]
    return sum((1 << i) * 3**(len(positions)-1-j) for j, i in enumerate(positions))


def cross(v, w):
    return (v[1]*w[2]-v[2]*w[1], v[2]*w[0]-v[0]*w[2], v[0]*w[1]-v[1]*w[0])


def minima(p, n, u):
    first = None
    for h in range(1, p+1):
        for b, c, m in itertools.product(range(-h, h+1), range(h+1), range(-h, h+1)):
            if max(abs(b), c, abs(m)) != h or (m-b*n-c*u) % p:
                continue
            v = (b, c, m)
            if first is None:
                first = (v, h)
            elif any(cross(first[0], v)):
                require((first[1]*h)**3 <= p*p, ('minima product', p, n, u, first, v, h))
                return first, (v, h)
    raise RuntimeError(('minima not found', p, n, u))


def polynomial(p, n, u, v, direction):
    b, c, m = v
    a = (m-b*n-c*u)//p
    h = max(abs(b), c, abs(m))
    require(abs(a) <= h, ('a bound', p, n, u, v, a))
    delta = int(direction == -1)
    epsilon = -direction
    e = a+2*b+delta*c
    f = c-b
    U = 2**max(-e, 0) * 3**max(-f, 0) * 7**c
    V = epsilon**c * 2**max(e, 0) * 3**max(f, 0)
    coefficients = [0] * (max(max(m, 0)+c, max(-m, 0))+1)
    for j in range(c+1):
        coefficients[max(m, 0)+j] += U * math.comb(c, j) * (-1)**(c-j)
    coefficients[max(-m, 0)] -= V
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    require(any(coefficients), ('zero polynomial', p, n, u, v))
    require(coefficients[0] != 0 and sum(coefficients) != 0, ('zero or one root', p, n, u, v))
    require(len(coefficients)-1 <= 2*h, ('degree', p, n, u, v))
    base = 121 if direction == 1 else 65
    require(sum(map(abs, coefficients)) <= base**h, ('coefficient bound', p, n, u, v))
    return coefficients


def determinant(matrix):
    a = [row[:] for row in matrix]
    size = len(a)
    if not size:
        return 1
    sign, previous = 1, 1
    for k in range(size-1):
        if a[k][k] == 0:
            pivot = next((i for i in range(k+1, size) if a[i][k]), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign *= -1
        pivot = a[k][k]
        for i in range(k+1, size):
            for j in range(k+1, size):
                numerator = a[i][j]*pivot-a[i][k]*a[k][j]
                require(numerator % previous == 0, 'nonexact Bareiss division')
                a[i][j] = numerator//previous
            a[i][k] = 0
        previous = pivot
    return sign*a[-1][-1]


def resultant(f, g):
    d, e = len(f)-1, len(g)-1
    require(d or e, 'two constant polynomials')
    if d == 0:
        return f[0]**e
    if e == 0:
        return g[0]**d
    fd, gd = f[::-1], g[::-1]
    matrix = [[0]*i+fd+[0]*(e-1-i) for i in range(e)]
    matrix += [[0]*i+gd+[0]*(d-1-i) for i in range(d)]
    return determinant(matrix)


def evaluate(f, z, modulus):
    value = 0
    for coefficient in reversed(f):
        value = (value*z+coefficient) % modulus
    return value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-p', type=int, default=120)
    parser.add_argument('--lattice-max-p', type=int, default=36)
    args = parser.parse_args()
    counts = Counter()
    integers, depth_controls, largest = [], [], None
    lattice_cache = {}
    for p in range(2, args.max_p+1):
        for r in range(1, p):
            if math.gcd(p, r) != 1:
                continue
            counts['source_shapes'] += 1
            n, q = p-r, p+r
            D = 2**q-3**p
            modulus = abs(D)
            word = [((i+1)*p)//q-(i*p)//q for i in range(q)]
            positions = [i for i, bit in enumerate(word) if bit]
            labels = {((j+(j*r+p-1)//p-1) % q): (-j*r) % p for j in range(1, p+1)}
            require(set(labels) == set(positions), ('positions', p, r))
            z = root(p, n, modulus)
            actual = set()
            for idx, left in enumerate(positions):
                right = positions[(idx+1) % p]
                if (right-left) % q != 2:
                    continue
                u = labels[left]
                require(0 <= u < r and labels[right] == u+n, ('residue gap', p, r, left, right))
                for direction in (1, -1):
                    edges = [{left, (left+direction) % q}, {right, (right+direction) % q}]
                    initially_legal = not word[(left+direction) % q] and not word[(right+direction) % q]
                    if not initially_legal or not edges[0].isdisjoint(edges[1]):
                        continue
                    target = None
                    for order in ((left, right), (right, left)):
                        changed = word[:]
                        for old in order:
                            new = (old+direction) % q
                            require(changed[old] == 1 and changed[new] == 0, ('illegal second swap', p, r, u, direction, order))
                            changed[old], changed[new] = 0, 1
                        if target is None:
                            target = changed
                        else:
                            require(target == changed, ('noncommutation', p, r, u, direction))
                        actual.add((direction, u, order[0] == left))
                        counts['ordered_paths'] += 1
                        counts['plus_paths' if direction == 1 else 'minus_paths'] += 1
                        counts['positive_clock_paths' if D > 0 else 'negative_clock_paths'] += 1
                        counts['original_cut_crossing_paths'] += int(any(edge == {0, q-1} for edge in edges))
                        bvalue = carry(changed)
                        require(bvalue == carry_by_positions(changed), ('carry disagreement', p, r, u, direction))
                        g = math.gcd(modulus, bvalue)
                        t = pow(z, u, modulus)*(z-1) % modulus
                        gate = 3+7*t if direction == 1 else 6-7*t
                        require(g == math.gcd(modulus, gate), ('whole gcd gate', p, r, u, direction, g))
                        require(math.gcd(g, 42) == 1, ('bad prime content', p, r, u, direction, g))
                        counts['whole_gcd_checks'] += 1
                        if largest is None or g > largest['g']:
                            largest = {'p': p, 'r': r, 'u': u, 'direction': direction, 'g': g}
                        if order[0] == left and len(depth_controls) < 12:
                            square_prime = next((ell for ell in (5,11,13,17,19,23,29,31,37,41,43,47) if g % (ell*ell) == 0), None)
                            if square_prime:
                                depth_controls.append({'p':p,'r':r,'u':u,'direction':direction,'g':g,'square_prime':square_prime})
                        if g == modulus:
                            x, trajectory = bvalue//D, []
                            value = x
                            for bit in changed:
                                require(value % 2 == bit, ('integer parity', p, r, u, direction, value, bit))
                                trajectory.append(value)
                                value = (3*value+1)//2 if bit else value//2
                            require(value == x, ('not a cycle', p,r,u,direction))
                            integers.append({'p':p,'r':r,'u':u,'direction':direction,'first_source_position':order[0], 'word': ''.join(map(str, changed)), 'B':bvalue,'D':D,'fixed_point':x,'trajectory':trajectory})
                    if p <= args.lattice_max_p:
                        key = (p, n, u)
                        if key not in lattice_cache:
                            lattice_cache[key] = minima(p, n, u)
                            counts['lattices'] += 1
                        (v1,h1),(v2,h2) = lattice_cache[key]
                        f1, f2 = polynomial(p,n,u,v1,direction), polynomial(p,n,u,v2,direction)
                        R = resultant(f1,f2)
                        require(R != 0, ('zero resultant', p,r,u,direction,v1,v2))
                        require(evaluate(f1,z,g) == 0 and evaluate(f2,z,g) == 0, ('polynomial gate',p,r,u,direction,g))
                        require(R % g == 0, ('full content resultant',p,r,u,direction,g,R))
                        base = 121 if direction == 1 else 65
                        require(abs(R) <= base**(4*h1*h2), ('resultant bound',p,r,u,direction))
                        counts['resultants'] += 1
            expected = {(direction,u,first) for direction,us in ((1,range(max(0,2*r-p))),(-1,range(n,r))) for u in us for first in (False,True)}
            require(actual == expected, ('legality mismatch',p,r,actual^expected))
            counts['legality_comparisons'] += 1
    # Full-depth and degenerate-resultant controls beyond the small lattice census.
    controls = []
    for p,r,u,direction,expected_g in ((19,14,6,-1,25),(37,29,25,-1,121),(39,34,26,-1,125),(119,72,20,1,104233)):
        n,q = p-r,p+r
        word = [((i+1)*p)//q-(i*p)//q for i in range(q)]
        residue_position = {(-j*r)%p:(j+(j*r+p-1)//p-1)%q for j in range(1,p+1)}
        for residue in (u,u+n):
            old=residue_position[residue]
            new=(old+direction)%q
            require(word[old] == 1 and word[new] == 0, ('depth control legality',p,r,u,direction))
            word[old],word[new]=0,1
        D=2**q-3**p
        g=math.gcd(abs(D),carry(word))
        require(g == expected_g, ('depth control content',p,r,u,direction,g))
        (v1,h1),(v2,h2)=minima(p,n,u)
        f1,f2=polynomial(p,n,u,v1,direction),polynomial(p,n,u,v2,direction)
        R=resultant(f1,f2)
        require(R != 0 and R%g == 0, ('depth control resultant',p,r,u,direction,g,R))
        controls.append({'p':p,'r':r,'u':u,'direction':direction,'g':g,'vectors':[v1,v2],'resultant':R})
    f1=polynomial(3,1,0,(0,1,0),1)
    f2=polynomial(3,1,0,(3,1,0),1)
    constant=polynomial(3,1,0,(3,0,0),1)
    require(f1 == [-4,7] and f2 == [-31,63] and constant == [-5], 'dependent projection polynomials')
    require(abs(resultant(f1,f2)) == 35 and abs(resultant(constant,f1)) == 5, 'dependent projection resultants')
    counts['extra_resultant_controls'] += 6
    # Joint coefficient-envelope hostile: the unrestricted H-box worst corner.
    for h in range(1,21):
        for direction in (1,-1):
            delta = int(direction == -1)
            maximum = 0
            for a,b,c in itertools.product((-h,h),(-h,h),(0,h)):
                e, f = a+2*b+delta*c, c-b
                bound = 2**(max(-e,0)+c)*3**max(-f,0)*7**c + 2**max(e,0)*3**max(f,0)
                maximum=max(maximum,bound)
            require(maximum == (112 if direction == 1 else 56)**h+9**h, ('joint coefficient envelope',h,direction,maximum))
            counts['joint_envelope_checks'] += 1
    print(json.dumps({'status':'FINITE-EXACT','parameters':vars(args),'counts':dict(sorted(counts.items())),'integer_examples':integers,'prime_square_controls':depth_controls,'extra_resultant_controls':controls,'largest_content':largest},indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
