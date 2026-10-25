#!/usr/bin/env python3
"""Independent projective-coordinate census and fixed-family zeta control.

No producer imports. Finite checks supplement, not replace, geometric proofs.
"""
from collections import Counter
import argparse
from fractions import Fraction
from itertools import product
import json
from pathlib import Path


def require(test, message):
    if not test:
        raise RuntimeError(message)


def projective_points(p, dimension):
    for pivot in range(dimension + 1):
        for tail in product(range(p), repeat=dimension-pivot):
            yield (0,)*pivot + (1,) + tail


def rank(matrix, p):
    a = [list(row) for row in matrix]
    r = 0
    for c in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][c] % p), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][c] % p, -1, p)
        a[r] = [v*inv % p for v in a[r]]
        for i in range(len(a)):
            if i != r:
                v = a[i][c] % p
                a[i] = [(x-v*y) % p for x,y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


COEFFICIENTS = [(-2,1,1,0,0,0,0), (-2,0,0,1,1,0,0),
                (-2,0,0,0,0,1,1), (1,-1,0,-1,0,1,0)]


def projective_census(p, recorded):
    nu = next(a for a in range(1,p) if pow(a,(p-1)//2,p) == p-1)
    inverses = [(1,) + tuple(pow(nu,-1,p) if mask >> i & 1 else 1 for i in range(6))
                for mask in range(64)]
    compact, free, nodes = [0]*64, [0]*64, [0]*64
    total = 0
    for point in projective_points(p,6):
        total += 1
        squares = [a*a % p for a in point]
        for mask, inv in enumerate(inverses):
            t = [s*d % p for s,d in zip(squares,inv)]
            if ((t[1]+t[2]-2*t[0]) % p or
                (t[3]+t[4]-2*t[0]) % p or
                (t[5]+t[6]-2*t[0]) % p or
                (t[5]-t[1]-t[3]+t[0]) % p):
                continue
            compact[mask] += 1
            free[mask] += all(point)
            jacobian = [[2*c*y*d % p for c,y,d in zip(row,point,inv)]
                        for row in COEFFICIENTS]
            nodes[mask] += rank(jacobian,p) < 4
    require(total == sum(p**i for i in range(7)), 'ambient projective census')
    t3 = 5 if p == 3 else 3
    require(sum(compact) == 64*(p*p+p+1), 'compact mean')
    require(sum(free) == 64*(p*p-6*p+15-t3), 'open mean')
    require(sum(nodes) == 64*t3, 'singular node mean')
    record = next(c for c in recorded['cases'] if c['q'] == p)
    require(compact == record['twist_projective_counts'], 'all64 compact bins versus root-table implementation')
    require(free == record['twist_open_counts'], 'all64 open bins versus root-table implementation')
    require(nodes == record['twist_node_counts'], 'all64 Jacobian nodes versus triple incidence')
    return {'p':p, 'ambient_projective_points':total,
            'compact_counts':compact, 'free_counts':free, 'jacobian_node_counts':nodes,
            'compact_sum':sum(compact), 'free_sum':sum(free), 'node_sum':sum(nodes)}


def zeta_control():
    # An explicit finite G-set with commuting Frobenius. It is a finite etale
    # Fq-scheme with nonfree strata: stabilizer orders1,4,64.
    G = range(64)
    # A block has quotient G/H encoded by low bits, and an F-cycle factor.
    blocks = [(6,3,19), (4,4,9), (0,2,0)]
    degree = 16
    fixed = [[0]*(degree+1) for _ in G]
    W, B = [0]*(degree+1), [0]*(degree+1)
    for m in range(1,degree+1):
        for bits, cycle, shift in blocks:
            if m % cycle:
                continue
            B[m] += cycle
            for g in G:
                displacement = ((g ^ shift) & ((1 << bits)-1)) if m % 2 else 0
                if displacement == 0:
                    fixed[g][m] += cycle * (1 << bits)
        W[m] = fixed[0][m]
        require(sum(row[m] for row in fixed) == 64*(B[m] if m%2 else W[m]),
                ('fixed-family odd/even',m))
    # Check the claimed zeta identity in formal logarithms, including t^(2r)/(2r).
    lhs = [Fraction(sum(row[m] for row in fixed),m) for m in range(1,degree+1)]
    rhs = []
    for m in range(1,degree+1):
        coefficient = Fraction(64*B[m],m)
        if m % 2 == 0:
            coefficient += Fraction(32*(W[m]-B[m]),m//2)
        rhs.append(coefficient)
    require(lhs == rhs,'zeta exponent32 formal-log identity')
    require(any(W[m] != B[m] for m in range(2,degree+1,2)), 'nontrivial fixed-family hostile')
    return {'degree':degree, 'group_order':64, 'stabilizers':[1,4,64],
            'W_extension_counts':W[1:], 'B_extension_counts':B[1:],
            'log_coefficients':[str(c) for c in lhs]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True,
                        help='Recovery lane twist_finite_normal.json for all64 labeled comparisons')
    args = parser.parse_args()
    recorded = json.loads(args.reference.read_text())
    good_euler = 64*7 + 32*(-19) + 16*12 + 8*3
    require(good_euler == 56 and good_euler+24 == 80, 'Euler carrier distinction')
    print(json.dumps({'status':'PASS', 'direct_projective_checks':[projective_census(p, recorded) for p in (3,5)],
                      'zeta_control':zeta_control(), 'euler_controls':{'V':448,'X':56,'Y':80}},
                     indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
