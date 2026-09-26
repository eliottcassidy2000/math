"""Exact fixed-integer/residue tests and actual high-multiplicity segments.

No external packages, no probabilistic decisions, no floating-point filters.
"""
from collections import Counter
from fractions import Fraction
import json


def require(value, message):
    if not value:
        raise AssertionError(message)


def step(n, a=3, b=1):
    return n//2 if n % 2 == 0 else (a*n+b)//2


def word_residue(word):
    """Lift the start address one bit at a time and retain exact affine data."""
    r, A, C, modulus = 0, 1, 0, 1
    representatives = [0]
    for bit in word:
        if ((A*r+C)//modulus) % 2 != bit:
            r += modulus
        require((A*r+C) % modulus == 0, 'prefix arithmetic compatibility')
        require(((A*r+C)//modulus) % 2 == bit, 'next parity')
        if bit:
            A, C = 3*A, 3*C+modulus
        modulus *= 2
        require((A*r+C) % modulus == 0, 'extended affine compatibility')
        representatives.append(r)
    return r, A, C, representatives


def actual_word(n, k, a=3, b=1):
    word, orbit = [], [n]
    for _ in range(k):
        word.append(n % 2)
        n = step(n, a, b)
        orbit.append(n)
    return word, orbit


def audit_addresses():
    cases = 0
    for n in range(1, 513):
        word, orbit = actual_word(n, 24)
        r, A, C, representatives = word_residue(word)
        for k, rk in enumerate(representatives):
            require(rk == n % (1 << k), 'least residues equal fixed start reductions')
            if k:
                require(rk-representatives[k-1] in (0, 1 << (k-1)), 'monotone binary lift')
            cases += 1
        require(r == n and (A*n+C)//(1 << 24) == orbit[-1], 'fixed-start stabilization and affine endpoint')
    negative_word, _ = actual_word(-1, 64)
    _, _, _, negative_reps = word_residue(negative_word)
    require(all(r == (1 << k)-1 for k, r in enumerate(negative_reps)), 'negative fixed point escapes positive height')

    # Computable sparse 2-adic input; its parity prefix is extracted exactly.
    positions = [1, 2, 6, 24, 120]
    sparse = sum(1 << p for p in positions)
    word, _ = actual_word(sparse, 121)
    _, _, _, reps = word_residue(word)
    require(all(reps[k] == sparse % (1 << k) for k in range(122)), 'sparse address reconstruction')
    samples = [{'k': k, 'representative': str(reps[k]),
                'normalized_exact': str(Fraction(reps[k], 1 << k))}
               for k in [7, 23, 25, 119, 121]]
    print(json.dumps({'fixed_start_address_checks': cases, 'negative_one_representatives': '2^k-1',
                      'sparse_factorial_bit_positions_through_depth121': positions,
                      'sparse_samples': samples}))


def hover_word(m):
    word, e, pow3 = [], 0, 1
    for j in range(1, m+1):
        previous = e
        while pow3 <= (1 << j):
            e += 1
            pow3 *= 3
        require(e-previous in (0, 1), 'critical ceiling word is binary')
        require((1 << j) < pow3 < 3*(1 << j), 'prefix gains strictly between1 and3')
        word.append(e-previous)
    return word


def audit_height_carrier():
    cases, carry_witness = 0, None
    for k in range(1, 11):
        modulus = 1 << k
        for residue in range(modulus):
            word, _ = actual_word(residue, k)
            A, C, two, cap = 1, 0, 1, None
            for bit in word:
                if bit:
                    A, C = 3*A, 3*C+two
                two *= 2
                if A < two:
                    bound = C//(two-A)
                    cap = bound if cap is None else min(cap, bound)
            for lift in [0, 1, 2, 17]:
                n = residue+lift*modulus
                if not n:
                    continue
                actual, orbit = actual_word(n, k)
                require(actual == word, 'carrier residue class preserves word')
                no_descent = all(x >= n for x in orbit[1:])
                require(no_descent == (cap is None or n <= cap), 'exact affine height cutoff iff')
                if n > 1 and cap is not None and no_descent and carry_witness is None:
                    carry_witness = {'n': n, 'word': ''.join(map(str, word)), 'cap': cap, 'orbit': orbit}
                cases += 1
    print(json.dumps({'height_carrier_exact_cases': cases,
                      'contracting_prefix_with_positive_carry_no_descent_witness': carry_witness}))


def cylinder_cap(word):
    A, C, two, cap = 1, 0, 1, None
    for bit in word:
        if bit:
            A, C = 3*A, 3*C+two
        two *= 2
        if A < two:
            bound = C//(two-A)
            cap = bound if cap is None else min(cap, bound)
    return cap


def exact_bad_mass(k, z=Fraction(1, 2)):
    modulus, total = 1 << k, Fraction(0)
    for residue in range(modulus):
        word, _ = actual_word(residue, k)
        cap = cylinder_cap(word)
        first = residue if residue >= 2 else residue+modulus
        if cap is not None and first > cap:
            continue
        factor = 1 if cap is None else 1-z**(((cap-first)//modulus+1)*modulus)
        total += z**first*factor/(1-z**modulus)
    return total


def audit_weighted_mass():
    z, H = Fraction(1, 2), 64
    values = [z**2/(1-z)]
    for k in range(1, 10):
        value = exact_bad_mass(k, z)
        low = sum(z**n for n in range(2, H+1)
                  if all(x >= n for x in actual_word(n, k)[1][1:]))
        require(low <= value <= low+z**(H+1)/(1-z), 'uniform fixed-height tail certificate')
        require(value <= values[-1], 'monotone exact weighted masses')
        values.append(value)
    require(values[2] == values[3] == Fraction(2, 15), 'one-step strict contraction is false')
    stopping = {}
    for n in range(2, 28):
        x = n
        for k in range(1, 100):
            x = step(x)
            if x < n:
                stopping[n] = k
                break
        require(n in stopping, 'bounded first descent control')
    start = max(stopping[n] for n in range(2, 27))
    require(stopping[27] == 59, '27 first descent')
    for k in range(start, 59):
        mass = sum(z**n for n in range(2, 28) if stopping[n] > k)
        require(mass == z**27, 'fixed-integer mass plateau')
    print(json.dumps({'full_support_weight': 'z=1/2', 'Z2_exact': str(values[2]),
                      'Z3_exact': str(values[3]), 'one_step_uniform_strict_contraction': False,
                      'exact_mass_checks_through_depth': 9, 'uniform_tail_height': H,
                      'tail_bound_exact': str(z**(H+1)/(1-z)),
                      'height27_plateau_depths_inclusive': [start, 58],
                      'plateau_mass_exact': str(z**27), '27_first_descent': 59}))


def multiplicity_segment(m):
    require(m >= 6, 'hover length')
    D = (m-1).bit_length()+6
    K = m+D
    word = hover_word(m)+[0]*D
    residue, A, C, _ = word_residue(word)
    n = (1 << (4*K))+residue
    actual, orbit = actual_word(n, K)
    X, L, depth = 1 << (4*K+3), 4*K+3, D-3
    require(actual == word, 'actual positive orbit realizes entire specified word')
    require(all(n <= x < 4*n for x in orbit[:m+1]), 'physical hover bounds include carry')
    require(all(0 < x < X for x in orbit), 'one positive segment below declared cutoff')
    require(min(orbit) > 3**K, 'uniform height excludes repeats of period at mostK')
    require(len(set(orbit)) == K+1, 'distinct actual orbit terms')
    require((A*n+C)//(1 << K) == orbit[-1], 'exact affine endpoint')
    landings = Counter()
    source_landing = {}
    for i, x in enumerate(orbit):
        for j in range(i+1, min(len(orbit), i+L+1)):
            if orbit[j]*(1 << depth) < x:
                landings[j] += 1
                source_landing[i] = j
                break
    require(all(i in source_landing for i in range(m)), 'every hover source dips')
    require(all(m+D-4 <= source_landing[i] <= m+D for i in range(m)), 'hover landings share final five positions')
    require(set(i for i in source_landing if i >= m) == {m, m+1, m+2}, 'exactly first three drop sources dip')
    require(sum(landings.values()) == m+3 and len(landings) <= 5, 'linear average multiplicity')
    return {'m': m, 'D': D, 'segment_steps': K, 'L': L,
            'theta_L_exact': depth, 'theta_exact': str(Fraction(depth, L)),
            'source_binary_length': n.bit_length(), 'dippers': m+3,
            'landing_positions_from_drop_start': [[j-m, v] for j, v in sorted(landings.items())],
            'number_of_landings': len(landings),
            'mean_multiplicity_exact': str(Fraction(m+3, len(landings))),
            'mean_divided_by_L_exact': str(Fraction(m+3, len(landings)*L)),
            'endpoint_allowance_absorbs_whole_segment': K < L}


def audit_cycles():
    _, negative_band = actual_word(7, 2, b=-3)
    _, negative_landing = actual_word(6, 2, b=-7)
    require(negative_band == [7,9,12] and all(6<x<=12 for x in negative_band),
            'negative-b one-bit-band premise needs a higher cutoff')
    require(negative_landing == [6,3,1] and 2<6 and 2<3 and 6>4,
            'low odd landing need not give claimed halving band')
    print(json.dumps({'negative_b_band_counterexample':negative_band,
                      'negative_b_odd_landing_counterexample':negative_landing,
                      'safe_scope':'b>0 or windows above5|b|; O_b(k) low-core windows on distinct orbits'}))
    _, minus = actual_word(5, 3, b=-1)
    _, five = actual_word(13, 7, a=5)
    require(minus == [5, 7, 10, 5], 'sign-hostile half-map cycle')
    require(five == [13, 33, 83, 208, 104, 52, 26, 13], 'multiplier-hostile half-map cycle')
    print(json.dumps({'computable_positive_hostile_3n_minus1': minus,
                      'computable_positive_hostile_5n_plus1': five,
                      'computability_implies_termination': False}))


def census_first_descent(limit=1000000, cap=10000):
    """Direct exhaustive census; no unproved stopping-time memoization."""
    records, maximum = [], 0
    for n in range(2, limit+1):
        x = n
        for k in range(1, cap+1):
            x = step(x)
            if x < n:
                break
        else:
            raise AssertionError('unresolved first descent at start '+str(n))
        if k > maximum:
            maximum = k
            y, full, visits27, visits41, peak = n, 0, n == 27, n == 41, n
            while y != 1 and full < cap:
                y = step(y)
                full += 1
                visits27 = visits27 or y == 27
                visits41 = visits41 or y == 41
                peak = max(peak, y)
            require(y == 1, 'record full stopping time resolved')
            records.append({'start': n, 'first_descent_steps': k,
                            'first_descent_value': x, 'steps_to_first1': full,
                            'visits27_including_start': visits27,
                            'visits41_including_start': visits41,
                            'shortcut_peak': peak})
    for i, record in enumerate(records):
        record['minimum_survivor_depth_interval'] = [
            0 if i == 0 else records[i-1]['first_descent_steps'],
            record['first_descent_steps']-1]
    require(records[0]['start'] == 2 and records[0]['first_descent_steps'] == 1,
            'census initial control')
    require(any(row['start'] == 27 and row['first_descent_steps'] == 59 for row in records),
            'census27 control')
    print(json.dumps({'first_descent_census_universe': [2, limit],
                      'map': 'even n/2; odd (3n+1)/2', 'per_start_step_cap': cap,
                      'unresolved_starts': 0, 'strict_record_count': len(records),
                      'records': records,
                      'minimum_survivor_above_census_after_depth': maximum,
                      'asymptotic_record_claimed': False}))


def main():
    audit_addresses()
    audit_height_carrier()
    audit_weighted_mass()
    audit_cycles()
    for m in [6, 8, 12, 16, 24, 40, 64, 96]:
        print(json.dumps(multiplicity_segment(m)))
    census_first_descent()
    print('ALL CHECKS PASSED; NO COLLATZ PROOF CLAIMED')


if __name__ == '__main__':
    main()
