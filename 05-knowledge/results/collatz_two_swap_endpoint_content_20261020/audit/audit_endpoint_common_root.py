#!/usr/bin/env python3
"""Independent direct bit-word/carry audit. No repository math imports or asserts."""
from collections import Counter
from hashlib import sha256
from math import gcd, isqrt
import json


def require(test, message):
    if not test:
        raise RuntimeError(message)


def bezout(a, b):
    old_r, r, old_s, s, old_t, t = a, b, 1, 0, 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    require(old_r == 1, 'noncoprime clock')
    return old_s, old_t


def bit_carry(bits):
    # Compose each chronological affine branch, independently of the closed sum.
    multiplier, numerator, denominator = 1, 0, 1
    for bit in bits:
        if bit:
            multiplier *= 3
            numerator = 3 * numerator + denominator
        denominator *= 2
    return numerator


def sum_carry(bits):
    positions = [i for i, bit in enumerate(bits) if bit]
    p = len(positions)
    return sum(2 ** i * 3 ** (p - 1 - j) for j, i in enumerate(positions))


def actual_word(p, r, family):
    q = p + r
    source = [((i + 1) * p) // q - (i * p) // q for i in range(q)]
    ordinary = [i for i, bit in enumerate(source) if bit]
    positions = [ordinary[-1]] + ordinary[:-1]
    residues = [(-j * r) % p for j in range(p)]
    move_residues = [p - 1 - r, p - 1] if family == 'plus' else [r, 0]
    delta = 1 if family == 'plus' else -1
    labels = [residues.index(a) for a in move_residues]
    forbidden_first = (positions[labels[1]] + delta) % q
    require(source[forbidden_first] == 1, 'reversed chain order must be blocked')
    current = source.copy()
    words = [''.join(map(str, current))]
    edges = []
    changes = []
    for label in labels:
        at = positions[label]
        to = (at + delta) % q
        require(current[at] == 1 and current[to] == 0, 'non-physical swap')
        edges.append(at if delta == 1 else to)
        changes.append([at, to])
        current[at], current[to] = 0, 1
        positions[label] = to
        require(sum(current) == p, 'weight lost')
        words.append(''.join(map(str, current)))
    return source, current, words, edges, changes


def short_relation(p, r, family, content):
    N = isqrt(p - 1) + 1
    require(N < p, 'short-relation weight scope')
    bins = {}
    for i in range(N+1):
        residue = (i*r) % p
        bucket = N*residue // p
        if bucket in bins:
            j = bins[bucket]
            b = i-j
            a = -((i*r)//p - (j*r)//p)
            m = a*p + b*r
            break
        bins[bucket] = i
    else:
        raise RuntimeError('pigeonhole collision absent')
    require(1 <= b <= N and -N <= a <= 0, 'short-relation coefficient bound')
    require(0 < abs(m)*N < p, 'short-relation strict residue bound')
    require(abs(m) <= N, 'short-relation square-root bound')
    if family == 'plus':
        left, right = (10**m*2**(b-a),13**m*3**b) if m>=0 else (
            13**(-m)*2**(b-a),10**(-m)*3**b)
        base = 52
    else:
        left, right = (9**m*2**(b-a),5**m*3**b) if m>=0 else (
            5**(-m)*2**(b-a),9**(-m)*3**b)
        base = 36
    certificate = left-right
    require(certificate != 0, 'zero short-relation certificate')
    require(certificate % content == 0, 'content does not divide certificate')
    require(content <= abs(certificate) < base**N, 'subexponential content bound')
    if p >= 25:
        require(25*N <= 6*p, 'ceil-square-root linear comparison')
    return {'N':N, 'i':i, 'j':j, 'a':a, 'b':b, 'm':m,
            'certificate':certificate, 'bound_base':base}


def audit(p, r, family, every_cut=True):
    require(p >= 3 and p < 2*r < 2*p and gcd(p, r) == 1, 'scope')
    q = p + r
    source, final, words, edges, changes = actual_word(p, r, family)
    B = bit_carry(final)
    require(B == sum_carry(final), 'two independent carry computations')
    D = 2**q - 3**p
    actual = gcd(abs(D), B)
    if family == 'plus':
        first = 2*13**p - 10**p
        second = 3*13**r - 2*10**r
        numer, denom, excluded = 10, 13, 6*5*13
    else:
        first = 9**p - 2*5**p
        second = 2*9**r - 3*5**r
        numer, denom, excluded = 9, 5, 6*5
    predicted = gcd(first, second)
    require(actual == predicted, f'content mismatch at {p},{r},{family}')
    require(gcd(predicted, excluded) == 1, 'denominator-prime exclusion')
    u, v = bezout(p, r)
    if predicted > 1:
        m = predicted
        inv2 = pow(2, -1, m)
        z = numer * pow(denom, -1, m) % m
        clock = pow(2, u, m) * pow(3*inv2 % m, v, m) % m
        require(clock == z, 'canonical common root not identified')
        require(pow(z, p, m) == 2 % m, 'p clock')
        require(pow(z, r, m) == 3*inv2 % m, 'r clock')
        require(gcd(z-1, m) == 1, 'z-1 not a unit')
        if family == 'plus':
            sparse = (3-pow(z,p-1,m)+pow(z,p-r,m)-pow(z,p-r-1,m)) % m
            require((3*z*sparse - (13*z-10)) % m == 0, 'plus linear reduction')
        else:
            sparse = (3-z-pow(z,r+1,m)+pow(z,r,m)) % m
            require((2*sparse - (9-5*z)) % m == 0, 'minus linear reduction')
        require(sparse == 0, 'sparse common root')
        # Verify whole-modulus carry transport at the final moved label 0.
        ordinary = [i for i, bit in enumerate(source) if bit]
        positions = [ordinary[-1]] + ordinary[:-1]
        residues = [(-j*r) % p for j in range(p)]
        d = [0]*p
        for a in ([p-1-r, p-1] if family == 'plus' else [r, 0]):
            d[residues.index(a)] = 1 if family == 'plus' else -1
        anchor = (positions[0]+d[0]) % q
        anchored = final[anchor:] + final[:anchor]
        S = sum(pow(z,a,m) for a in range(p)) % m
        perturbation = sum((pow(2,move,m)-1)*pow(z,a,m) for a,move in zip(residues,d)) % m
        transported = pow(2,-d[0],m)*pow(3,p-1,m)*(S+perturbation) % m
        require(bit_carry(anchored) % m == transported, 'moved-anchor carry identity')
        require(transported == 0, 'sufficiency at moved anchor')
    cut_checks = 0
    if every_cut:
        prior = B
        for cut in range(q):
            rotated = final[cut:] + final[:cut]
            value = bit_carry(rotated)
            require(gcd(abs(D), value) == actual, 'rotation content loss')
            if cut:
                require(2*value == 3**final[cut-1]*prior + final[cut-1]*D,
                        'exact signed-D one-step cut transport')
            prior = value
            cut_checks += 1
    return {'p':p, 'r':r, 'q':q, 'family':family, 'B':B, 'D':D, 'g':actual,
            'gcd_formula':predicted, 'words':words, 'ordered_edges':edges,
            'moves':changes, 'cut_checks':cut_checks,
            'first_gcd_operand':first, 'second_gcd_operand':second,
            'short_relation':short_relation(p,r,family,actual)}


def main():
    counts = Counter()
    nontrivial = {'plus':[], 'minus':[]}
    examples = []
    trace = sha256()
    shape_count = 0
    for p in range(3,101):
        for r in range(p//2+1,p):
            if gcd(p,r) != 1:
                continue
            shape_count += 1
            for family in ['plus','minus']:
                row = audit(p,r,family)
                counts[family+'_instances'] += 1
                counts['exact_rotated_carry_checks'] += row['cut_checks']
                counts['actual_swaps'] += 2
                counts['blocked_reverse_orders'] += 1
                counts['short_relation_certificates'] += 1
                counts['cut_crossing_swaps'] += sum(edge == p+r-1 for edge in row['ordered_edges'])
                counts['positive_D' if row['D']>0 else 'negative_D'] += 1
                trace.update(json.dumps(row,sort_keys=True,separators=(',',':')).encode()+b'\n')
                if row['g'] > 1:
                    nontrivial[family].append({k:row[k] for k in ['p','r','q','g']})
                if (p,r) in [(3,2),(11,9),(17,12),(28,17),(54,29)]:
                    examples.append(row)
    # These controls deliberately retain nontrivial prime-power depth.
    controls = [audit(418,237,'plus'), audit(203,165,'minus')]
    require(controls[0]['g'] == 17**2*41, 'plus mixed-prime / prime-square control')
    require(controls[1]['g'] == 13**2, 'minus prime-square control')
    for row in controls:
        row['factorization'] = [[17,2],[41,1]] if row['family']=='plus' else [[13,2]]
    receipt = {
        'status':'FINITE-EXACT',
        'universe':'Every 3<=p<=100, p/2<r<p, gcd(p,r)=1; q=p+r; both ordered endpoint chains.',
        'shape_count':shape_count,
        'counts':dict(sorted(counts.items())),
        'full_case_trace_sha256':trace.hexdigest(),
        'nontrivial_content':nontrivial,
        'nontrivial_counts':{family:len(rows) for family,rows in nontrivial.items()},
        'examples':examples,
        'prime_power_controls':controls,
        'scope_notes':[
            'No repository mathematical implementation imported.',
            'Words use floor differences; carries use both affine composition and the positional sum.',
            'All original-cut rotations checked for every enumerated word and both prime-power controls.',
            'Only specified endpoint chains are tested; no cycle exclusion or bound for all two-edit words claimed.',
            'The minus polynomial is H=2Qraw, and its linear identity is 2H=9-5z.',
            'Every case includes an exact nonzero square-root-length relation certificate divisible by its full content.'
        ]
    }
    print(json.dumps(receipt,indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
