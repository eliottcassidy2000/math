"""Exact finite CRT/geometry controls; all checks remain active under -O.

The infinite density statements have separate proofs in the companion note.
No finite sieve is confused with global squarefreeness or an orbit law.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt, prod
from pathlib import Path
import argparse
import json


def require(test, message):
    if not test:
        raise ArithmeticError(message)


def frac(a, b):
    f = Fraction(a, b)
    return {'numerator': f.numerator, 'denominator': f.denominator}


def crt_weights(moduli):
    modulus = prod(moduli)
    return modulus, [(modulus // m) * pow(modulus // m, -1, m) for m in moduli]


def digit_split(n, p, depth):
    a = b = 0
    place = 1
    for _ in range(depth):
        a += (n % p) * place
        n //= p
        b += (n % p) * place
        n //= p
        place *= p
    return a, b


def digit_merge(a, b, p, depth):
    n = 0
    place = 1
    for _ in range(depth):
        n += (a % p) * place
        a //= p
        place *= p
        n += (b % p) * place
        b //= p
        place *= p
    return n


def transport(primes, depth):
    mod, weights = crt_weights([p**depth for p in primes])
    source_mod, inverse_weights = crt_weights([p**(2*depth) for p in primes])

    def forward(n):
        pieces = [digit_split(n, p, depth) for p in primes]
        return (sum(a*w for (a, _), w in zip(pieces, weights)) % mod,
                sum(b*w for (_, b), w in zip(pieces, weights)) % mod)

    def inverse(a, b):
        return sum(digit_merge(a % p**depth, b % p**depth, p, depth)*w
                   for p, w in zip(primes, inverse_weights)) % source_mod
    return mod, source_mod, forward, inverse


def rotate(v):
    a, b = v
    return (-b, a+b)


def quarter_turn(v):
    a, b = v
    return (-b, a)


def orbit(v, operation, modulus=None):
    result = []
    n = v
    while n not in result:
        result.append(n)
        n = operation(n)
        if modulus is not None:
            n = (n[0] % modulus, n[1] % modulus)
    require(n == v, 'noncyclic orbit input')
    return result


def squarefree(n):
    return all(n % (d*d) != 0 for d in range(2, isqrt(abs(n))+1))


def main():
    crt_rows = []
    universes = [([2], 1), ([3], 1), ([5], 1), ([2,3], 1),
                 ([2,3,5], 1), ([2,3,5,7], 1), ([2,3], 2), ([2,3], 3)]
    for primes, depth in universes:
        mod, source_mod, forward, inverse = transport(primes, depth)
        seen = set()
        good = 0
        M = prod(primes)
        for n in range(source_mod):
            a, b = forward(n)
            require(inverse(a, b) == n, ('inverse transport', primes, depth, n))
            require((a,b) not in seen, ('bijection collision', primes, depth, n))
            seen.add((a,b))
            left = all(n % (p*p) != 0 for p in primes)
            right = gcd(gcd(a,b), M) == 1
            require(left == right, ('local predicate', primes, depth, n))
            good += left
            if depth > 1:
                _, _, previous, _ = transport(primes, depth-1)
                require((a % (M**(depth-1)), b % (M**(depth-1))) ==
                        previous(n % (M**(2*depth-2))), 'inverse-limit compatibility')
        expected = prod(p**(2*depth)-p**(2*depth-2) for p in primes)
        require(good == expected and len(seen) == mod*mod, 'CRT population count')
        crt_rows.append({'primes': primes, 'depth': depth, 'source_modulus': source_mod,
                         'target_modulus_each_coordinate': mod, 'good_count': good,
                         'exact_density': frac(good, source_mod)})

    _, _, f2, _ = transport([2], 1)
    require(f2(2) != tuple((x+y) % 2 for x,y in zip(f2(1), f2(1))), 'additive hostile')
    require(f2(0) != tuple(x*x % 2 for x in f2(2)), 'multiplicative hostile')
    _, _, f6, _ = transport([2,3], 1)
    require(f6(25) == (1,2) and not squarefree(25) and gcd(*f6(25)) == 1,
            'omitted-square-prime hostile')
    require(f6(29) == (5,0) and squarefree(29) and gcd(*f6(29)) != 1,
            'omitted-gcd-prime hostile')

    phi = [0] + [sum(gcd(k,n) == 1 for k in range(1,n+1)) for n in range(1,101)]
    geom_rows = []
    for N in [1,2,3,5,10,25,50,100]:
        hexagon = {(a,b) for a in range(-N,N+1) for b in range(-N,N+1)
                   if (a,b) != (0,0) and max(abs(a),abs(b),abs(a+b)) <= N}
        square = {(a,b) for a in range(-N,N+1) for b in range(-N,N+1) if (a,b) != (0,0)}
        primitive_hex = {v for v in hexagon if gcd(*v) == 1}
        primitive_square = {v for v in square if gcd(*v) == 1}
        for population, operation, size in [(hexagon,rotate,6), (square,quarter_turn,4)]:
            for v in population:
                o = orbit(v,operation)
                require(len(o) == size and all(w in population for w in o), 'free geometric action')
                require(all(gcd(*w) == gcd(*v) for w in o), 'content preservation')
                if size == 6:
                    require(o[3] == (-v[0],-v[1]), 'central sign')
                    require({(-a,-b) for a,b in o[::2]} == set(o[1::2]), 'two triangles')
        totient_sum = sum(phi[1:N+1])
        require(len(hexagon) == 3*N*(N+1) and len(square) == 4*N*(N+1), 'ambient counts')
        require(len(primitive_hex) == 6*totient_sum and len(primitive_square) == 8*totient_sum,
                'primitive counts')
        visible = sum(gcd(a,b) == 1 for a in range(1,N+1) for b in range(1,N+1))
        require(visible == 2*totient_sum-1, 'positive square relation')
        for radius in range(1,N+1):
            shell = {v for v in hexagon if max(abs(v[0]),abs(v[1]),abs(sum(v))) == radius}
            require(len(shell) == 6*radius and sum(gcd(*v)==1 for v in shell) == 6*phi[radius],
                    'hexagon shell count')
        require(Fraction(len(primitive_hex),len(hexagon)) ==
                Fraction(len(primitive_square),len(square)), 'equal finite proportions')
        geom_rows.append({'N': N, 'hexagon_count': len(hexagon), 'primitive_hexagon_count': len(primitive_hex),
                          'square_count': len(square), 'primitive_square_count': len(primitive_square),
                          'positive_square_primitive_count': visible,
                          'common_exact_density': frac(len(primitive_hex),len(hexagon)),
                          'hexagon_rotation_orbits': len(hexagon)//6,
                          'square_rotation_orbits': len(square)//4})

    residue_orbits = []
    for M in (2,3,5,6,30):
        left = {(a,b) for a in range(M) for b in range(M) if gcd(gcd(a,b),M) == 1}
        lengths = Counter()
        while left:
            seed = min(left)
            o = orbit(seed,rotate,M)
            require(all(v in left for v in o), 'modular orbit partition')
            left.difference_update(o)
            lengths[len(o)] += 1
        residue_orbits.append({'modulus': M, 'orbit_length_counts': dict(sorted(lengths.items()))})
    require(residue_orbits[0]['orbit_length_counts'] == {3:1}, 'mod2 collapse')
    require(residue_orbits[1]['orbit_length_counts'] == {2:1,6:1}, 'mod3 collapse')
    require(residue_orbits[3]['orbit_length_counts'] == {6:4}, 'mod6 control')

    fixed = {}
    for b in (1,5):
        candidates = []
        for k in range(1, (3+abs(b)).bit_length()):
            gap = 2**k-3
            if b % gap == 0:
                n = b//gap
                numerator, actual = 3*n+b, 0
                while numerator % 2 == 0:
                    numerator //= 2
                    actual += 1
                require(numerator == n and actual == k, 'fixed-point guard')
                candidates.append({'n':n, 'k':k})
        fixed[str(b)] = candidates
    require([x['n'] for x in fixed['1']] == [-1,1], 'parameter1 fixed points')
    require([x['n'] for x in fixed['5']] == [-5,5,1], 'parameter5 fixed points')

    N = 600_000
    sf = bytearray([1])*(N+1)
    sf[0] = 0
    for d in range(2,isqrt(N)+1):
        square = d*d
        for n in range(square,N+1,square):
            sf[n] = 0
    rows = []
    for residue in (1,3,5):
        count = len(range(residue,N+1,6))
        good = sum(sf[n] for n in range(residue,N+1,6))
        rows.append({'residue_mod6':residue,'population':count,'squarefree_count':good,
                     'exact_density':frac(good,count),
                     'limiting_density':'6/pi^2' if residue==3 else '9/pi^2'})
    Q = sum(sf)
    require(Fraction(2*Q,2*N)==Fraction(Q,N), 'sign duplication normalization')
    require(sum(squarefree(n) for n in range(1,10)) == 6, 'squarefree height hostile')
    require(sum(gcd(a,b)==1 for a in range(1,4) for b in range(1,4)) == 7, 'lattice height hostile')
    return {'status':'FINITE-EXACT controls; all infinite statements proved separately',
            'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'CRT_universes':crt_rows, 'CRT_inputs_checked':sum(x['source_modulus'] for x in crt_rows),
            'geometry':geom_rows, 'modular_rotation_orbits':residue_orbits,
            'first_integral_six_orbit':orbit((1,0),rotate), 'fixed_points_by_parameter':fixed,
            'squarefree_prefix':{'X':N,'count':Q,'exact_density':frac(Q,N),'odd_rows':rows},
            'hostiles':{'addition':'f_2(1+1)=(0,1), f_2(1)+f_2(1)=(0,0)',
                        'multiplication':'f_2(2*2 mod4)=(0,0), f_2(2)^2=(0,1)',
                        'finite_sieve_25_to_1_2':True,'finite_sieve_29_to_5_0':True,
                        'global_integer_lift_obstruction':'n=2 forces b=0 at every odd prime, but b=1 modulo2',
                        'height_box':'Q(9)=6, visible positive pairs through3=7'},
            'infinite_Collatz_orbit_claim':False}


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=Path(__file__).with_suffix('.json'))
    args=parser.parse_args()
    result=main()
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':'PASS','output':str(args.output),'CRT_inputs_checked':result['CRT_inputs_checked']}))
