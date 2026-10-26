#!/usr/bin/env python3
"""Independent exact controls. Does not import or execute producer functions."""
from fractions import Fraction as F
from itertools import product
from math import lcm
from hashlib import sha256
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent / 'session64-reverse-valuation'


def check(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def vp(x, p=2):
    x = F(x)
    if not x:
        raise ValueError('infinite valuation')
    def z(n):
        n = abs(n)
        s = 0
        while n % p == 0:
            n //= p
            s += 1
        return s
    return z(x.numerator) - z(x.denominator)


def step(x, kind):
    x = F(x)
    check(vp(x) == 0, 'non-odd input')
    t = 3*x + (1 if kind == 'C' else -1)
    k = vp(t)
    check(k >= 1, 'invalid letter')
    return (t / 2**k) * (1 if kind == 'C' else 3**(k-1)), k


def closed_prefix(w, tail, kind):
    """Use printed series, not repeated inverse branches."""
    A, s = 0, F(0)
    for j, k in enumerate(w):
        s += -F(2**A, 3**(j+1)) if kind == 'C' else F(2, 3)**A / 3
        A += k
    m = F(2**A, 3**len(w)) if kind == 'C' else F(2, 3)**A
    return s + m*tail


def cycle(w, kind):
    s = closed_prefix(w, 0, kind)
    m = F(2**sum(w), 3**len(w)) if kind == 'C' else F(2, 3)**sum(w)
    return s/(1-m)


def modq(x, bits):
    x = F(x)
    return x.numerator * pow(x.denominator, -1, 2**bits) % 2**bits


def periodic_distance(u, v):
    A = 0
    for j in range(lcm(len(u), len(v))):
        p, q = u[j % len(u)], v[j % len(v)]
        if p != q:
            return A + min(p, q)
        A += p
    return None


def exact_cylinders_and_codes():
    B = 12
    words = [w for n in range(5) for w in product(range(1, 5), repeat=n)
             if sum(w) <= B-1]
    counters = {'residue_sources': 0, 'cylinder_words': 0, 'singular_prefixes': 0,
                'periodic_closures': 0, 'periodic_pair_isometries': 0}
    for kind in ['C', 'H']:
        census = {}
        for seed in range(1, 2**B, 2):
            counters['residue_sources'] += 1
            w, A, x = (), 0, F(seed)
            census.setdefault(w, []).append(seed)
            for _ in range(4):
                x, k = step(x, kind)
                A += k
                if A >= B:
                    break
                w += (k,)
                census.setdefault(w, []).append(seed)
        for w in words:
            A = sum(w)
            representative = modq(closed_prefix(w, 1, kind), A+1)
            expected = list(range(representative, 2**B, 2**(A+1)))
            check(census.get(w) == expected, 'full residue cylinder census')
            counters['cylinder_words'] += 1
            x = closed_prefix(w, F(-1 if kind == 'C' else 1, 3), kind)
            for k in w:
                x, actual = step(x, kind)
                check(k == actual, 'singular preimage prefix')
            try:
                step(x, kind)
            except ValueError as exc:
                check(str(exc) == 'infinite valuation', 'wrong zero error')
            else:
                raise RuntimeError('singular numerator admitted')
            counters['singular_prefixes'] += 1
    periodic_words = [w for n in range(1, 4) for w in product(range(1, 5), repeat=n)]
    points = {kind: [cycle(w, kind) for w in periodic_words] for kind in ['C', 'H']}
    for kind in ['C', 'H']:
        for w, x in zip(periodic_words, points[kind]):
            initial = x
            for k in w:
                x, actual = step(x, kind)
                check(actual == k, 'periodic first valuation')
            check(x == initial, 'periodic closure')
            if kind == 'H':
                check(F(1, 3) < x <= 1, 'periodic real interval')
                for pref in [(), (1,), (4,), (1, 4, 2)]:
                    check(F(1, 3) < closed_prefix(pref, x, kind) <= 1,
                          'eventually periodic real interval')
            counters['periodic_closures'] += 1
    for i, u in enumerate(periodic_words):
        for j in range(i+1, len(periodic_words)):
            d = periodic_distance(u, periodic_words[j])
            for kind in ['C', 'H']:
                difference = points[kind][i]-points[kind][j]
                check((difference == 0) if d is None else (vp(difference) == d),
                      'same-code metric sharpness')
            counters['periodic_pair_isometries'] += 1
    check(cycle((2,), 'C') == 1 and cycle((2,), 'H') == F(3, 5), 'integer locus loss')
    check(cycle((1,), 'C') == -1 and cycle((1,), 'H') == 1, 'sign locus loss')
    check(step(3, 'C') == (5, 1) and step(13, 'C') == (5, 3), 'sibling merge')
    check(step(1, 'C') == (1, 2) and step(F(-5, 3), 'C') == (-1, 2), 'prefix is not future')
    check(step(-1, 'C') == (-1, 1), 'negative fixed point')
    return counters


def inverse_integer_controls():
    # A branch with integer x requires v3(y)=k-1; |y|<=2047 bounds k<=7.
    # Enumerate extra k through 12 independently, including impossible branches.
    counts = {'target_inverses': 0, 'inverse_candidates': 0, 'source_images': 0,
              'backward_ray_steps': 0}
    for y in range(-2047, 2048, 2):
        candidates = []
        for k in range(1, 13):
            x = F(1, 3) + F(2, 3)**k*y
            counts['inverse_candidates'] += 1
            if x.denominator == 1:
                check(step(x, 'H') == (y, k), 'integer inverse roundtrip')
                candidates.append((int(x), k))
        r = vp(y, 3)
        m = y // 3**r
        predicted = [] if m % 3 != (-1)**r % 3 else [((2**(r+1)*m+1)//3, r+1)]
        check(candidates == predicted, 'forced v3 integer inverse')
        counts['target_inverses'] += 1
        x = y
        seen = set()
        while x != 1:
            check(x not in seen, 'backward integer cycle')
            seen.add(x)
            rr = vp(x, 3)
            predecessor = F(2**(rr+1)*(x//3**rr)+1, 3)
            if predecessor.denominator != 1:
                break
            predecessor = int(predecessor)
            check((1 < predecessor < x) if x > 1 else (x < predecessor < 0),
                  'backward ray must strictly approach bounded endpoint')
            x = predecessor
            counts['backward_ray_steps'] += 1
    images = set()
    for x in range(-4999, 5000, 2):
        y, _ = step(x, 'H')
        check(y.denominator == 1 and y not in images, 'injective signed integer image')
        check(y > x if x > 1 else y == x if x == 1 else y < x, 'signed integer drift')
        images.add(y)
        counts['source_images'] += 1
    return counts


def macro(a, b):
    check(a >= 0 and b >= 0, 'parameter domain')
    if b % 2:
        return None if a == 0 else (a-1, (3*b+3)//2)
    return a+2, 3*b//2+2


def debts():
    counts = {'words': 0, 'guard_comparisons': 0, 'cuts': 0, 'native_guard_replays': 0}
    for n in range(6):
        for w in product(range(1, 4), repeat=n):
            counts['words'] += 1
            levels = [0]
            for k in w:
                levels.append(levels[-1]+2*k-3)
            required = 1-min(levels)
            r = 1
            for k in reversed(w):
                r = max(1, r+3-2*k)
            check(r == required, 'reverse versus complete prefixes')
            for cut in range(n+1):
                left = (levels[cut], min(levels[:cut+1]))
                right_levels = [z-levels[cut] for z in levels[cut:]]
                right = (right_levels[-1], min(right_levels))
                check((left[0]+right[0], min(left[1], left[0]+right[1])) ==
                      (levels[-1], min(levels)), 'ordered summary composition')
                counts['cuts'] += 1
            source = modq(closed_prefix(w, 1, 'H'), sum(w)+1)
            while source < 5:
                source += 2**(sum(w)+1)
            for a0 in range(11):
                check(all(a0+d > 0 for d in levels) == (a0 >= required), 'inclusive guard iff')
                counts['guard_comparisons'] += 1
                a, b, executed = a0, source-4, True
                for k in w:
                    for j in range(k):
                        check(b % 2 == (1 if j == 0 else 0), 'native word realization parity')
                        nxt = macro(a, b)
                        if nxt is None:
                            executed = False
                            break
                        a, b = nxt
                    if not executed:
                        break
                    check(b % 2 == 1, 'returned odd state')
                execution_test = True if not w else all(a0+d > 0 for d in levels[:-1])
                check(executed == execution_test, 'execution-only final guard distinction')
                check((executed and a > 0) == (a0 >= required), 'macro inclusive guard')
                counts['native_guard_replays'] += 1
    # a=1, k=1 executes one odd macro, returns resource zero, then halts.
    check(macro(1, 1) == (0, 3) and macro(0, 3) is None, 'terminal zero hostile')
    for a in range(1, 25):
        x = 1+2**(a+1)
        for _ in range(a):
            x, k = step(x, 'H')
            check(k == 1, 'integer bad prefix')
        check(x == 1+2*3**a, 'integer bad endpoint')
    return counts


def ruin():
    # Backwards horizon recursion: probabilities from all a, not forward mass/safe sinks.
    f = [F(1)]
    out = {}
    for t in range(1, 49):
        g = [F(1)] + [F(0)]*t
        for a in range(1, t+1):
            k = 1
            while a+2*k-3 < len(f):
                g[a] += F(1, 2**k)*f[a+2*k-3]
                k += 1
        f = g
        if t in [1, 6, 12, 24, 48]:
            out[str(t)] = str(f[6] if t >= 6 else F(0))
    # Pair (a,b) means a+b*q with q^2=1-q; exact quadratic-field algebra.
    def mul(x, y):
        a,b = x; c,d = y
        return a*c+b*d, a*d+b*c-b*d
    q = (0, 1)
    q2 = mul(q, q)
    check(mul(q, (2-q2[0], -q2[1])) == (1, 0), 'exact harmonic kernel identity')
    power = (1, 0)
    for _ in range(6):
        power = mul(power, q)
    check(power == (5, -8), 'q^6 = 9-4sqrt5')
    vals = [F(out[str(t)]) for t in [1, 6, 12, 24, 48]]
    check(vals == sorted(vals) and vals[1] == F(1, 64), 'ruin monotonic and first boundary')
    for p in vals:
        check(9-p > 0 and (9-p)**2 > 80, 'exact finite-horizon lower bound')
    return out


def affine_ledgers():
    count = 0
    for kind in ['C', 'H']:
        for source in range(1, 128, 2):
            x, A, B, S, previous = F(source), 0, 0, F(0), F(source)
            for n in range(1, 13):
                old_A = A
                x, k = step(x, kind)
                A += k
                if kind == 'C':
                    B = 3*B+2**old_A
                    check(x == F(3**n*source+B, 2**A), 'forward affine carry')
                    check((x < source) == ((2**A-3**n)*source > B), 'exact descent inequality')
                else:
                    S += F(2, 3)**old_A/3
                    L = F(2, 3)**A*x
                    check(L == source-S, 'twisted exact normalized ledger')
                    check(source-1 < L < previous, 'twisted normalized real range')
                    previous = L
                count += 1
    check(2**2 > 3 and step(1, 'C')[0] == 1, 'coefficient alone does not imply descent')
    return count


def replay():
    # The producer drives k macros from a precomputed valuation. This replay
    # instead discovers k by uncompressed parameter parity/return detection.
    a, b = 0, 4
    initial = [[a,b]]
    while b % 2 == 0:
        a,b = macro(a,b)
        initial.append([a,b])
    check((a,b) == (6,23), 'blank to first odd source')
    letters = []
    minimum_odd, minimum_macro = a, a
    debt_min, debt = 0, 0
    for n in range(10000):
        old_a, old_x, k = a,b+4,0
        check(a > 0, 'unexpected finite barrier')
        while True:
            nxt = macro(a,b)
            check(nxt is not None, 'native parameter barrier')
            a,b = nxt
            k += 1
            minimum_macro = min(minimum_macro,a)
            if b % 2:
                break
        check(vp(3*old_x-1) == k, 'detected return length equals valuation')
        check(b+4 == 3**(k-1)*(3*old_x-1)//2**k, 'compressed odd return')
        check(a == old_a+2*k-3, 'compressed resource')
        debt += 2*k-3
        debt_min = min(debt_min, debt)
        letters.append(k)
        minimum_odd = min(minimum_odd,a)
        check(a == 6+2*sum(letters)-3*(n+1), 'global counter ledger')
    return {'status':'FINITE_PREFIX_ONLY','odd_returns':10000,'valuation_sum':sum(letters),
            'a':a,'min_odd_section_a':minimum_odd,'min_macro_a_from_first_odd':minimum_macro,
            'x_bit_length':(b+4).bit_length(),'initial_macro_states':initial,
            'first_35_valuations':letters[:35], 'reverse_requirement_for_prefix':1-debt_min,
            'valuation_csv_sha256':sha256(''.join(str(k)+',' for k in letters).encode()).hexdigest(),
            'x_sha256_hex_encoding':sha256(hex(b+4).encode()).hexdigest()}


def main():
    global PACKET
    parser = argparse.ArgumentParser()
    parser.add_argument('--packet', type=Path, default=PACKET)
    PACKET = parser.parse_args().packet.resolve()
    input_hashes = {}
    for line in (PACKET/'SHA256SUMS').read_text().splitlines():
        expected, name = line.split()
        actual = sha256((PACKET/name).read_bytes()).hexdigest()
        check(actual == expected, 'sealed input hash '+name)
        input_hashes[name] = actual
    observed = json.loads((PACKET/'report.json').read_text())
    answer = {'status':'PASS','producer_inputs_sha256':input_hashes,
              'codes':exact_cylinders_and_codes(),'integers':inverse_integer_controls(),
              'debts':debts(),'affine_ledger_steps':affine_ledgers(),
              'finite_horizon_ruin_mass_a6':ruin(),'macro_prefix':replay()}
    check(answer['finite_horizon_ruin_mass_a6'] == observed['finite_horizon_ruin_mass_a6'], 'producer ruin values')
    reference = observed['blank_antihydra_macro_prefix']
    for key,value in answer['macro_prefix'].items():
        if key not in ['min_odd_section_a','min_macro_a_from_first_odd']:
            check(reference[key] == value, 'producer replay field '+key)
    check(answer['macro_prefix']['min_odd_section_a'] == reference['min_a'], 'odd-only minimum scope')
    print(json.dumps(answer, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
