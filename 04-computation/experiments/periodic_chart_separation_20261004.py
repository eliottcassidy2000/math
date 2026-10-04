"""Universal repeated-chart separation; bounded counts and density controls.

Enumerates all positive compositions with p<=12 and 2^sum(w)<3^p.
Keeps all marked primitive words whose every prefix has slope >1.
Mass cutoff: 2<=m<=30. Literal controls: m=2,3, least source per chart.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb, gcd, lcm
from pathlib import Path
import json
import runpy

ROOT = Path(__file__).resolve().parents[2]
MAX_PERIOD = 12
MAX_REPEAT = 30


def need(ok, message):
    if not ok:
        raise ValueError(message)


def growth_bound(p):
    """Largest integer A with 2^A<3^p, without logarithms."""
    return (3**p).bit_length()-1


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for first in range(1, total-length+2):
            for tail in compositions(total-first, length-1):
                yield (first,)+tail


def primitive(word):
    p = len(word)
    return not any(word == word[:d]*(p//d) for d in range(1, p) if p % d == 0)


def summary(word):
    P, Q, S = 1, 1, 0
    for a in word:
        P, S, Q = 3*P, 3*S+Q, Q*2**a
    return P, Q, S


def favorable(word):
    A = 0
    for j, a in enumerate(word, 1):
        A += a
        if 3**j <= 2**A:
            return False
    return True


def mobius(n):
    out, p = 1, 2
    while p*p <= n:
        if n % p == 0:
            n //= p
            if n % p == 0:
                return 0
            out = -out
        p += 1
    return -out if n > 1 else out


def enumerate_charts():
    rows, charts = [], []
    for p in range(1, MAX_PERIOD+1):
        B = growth_bound(p)
        raw = aperiodic = favorable_count = 0
        necklaces = {}
        for A in range(p, B+1):
            for word in compositions(A, p):
                raw += 1
                if not primitive(word):
                    continue
                aperiodic += 1
                rotations = [word[j:]+word[:j] for j in range(p)]
                key = min(rotations)
                P, Q, S = summary(word)
                center = F(S, P-Q)
                need(center.numerator < Q*Q, "all growth words satisfy h<Q squared")
                previous = necklaces.get(key)
                if previous is None or center < previous[0]:
                    necklaces[key] = (center, word)
                if favorable(word):
                    favorable_count += 1
                    charts.append((word, P, Q, center.numerator, center.denominator))
        formula = sum(mobius(d)*comb(growth_bound(p//d), p//d)
                      for d in range(1, p+1) if p % d == 0)
        need(raw == comb(B, p), "independent composition count by hockey-stick identity")
        need(aperiodic == formula == p*len(necklaces), "Mobius and rotation-count cross-check")
        need(raw <= 2**(B-1) < F(3**p, 2), "uniform marked-word count bound")
        need(all(favorable(word) for _, word in necklaces.values()),
             "closest-to-zero canonical phase has every prefix expanding")
        rows.append(dict(period=p, maximum_total_valuation=B, all_growth_marked_words=raw,
                         primitive_growth_marked_words=aperiodic,
                         primitive_growth_necklaces=len(necklaces),
                         favorable_primitive_marked_words=favorable_count))
    need(sum(r['primitive_growth_necklaces'] for r in rows if r['period'] <= 8) == 130,
         "reproduce corrected inherited p<=8 count")
    need(sum(r['primitive_growth_necklaces'] for r in rows) == 5966,
         "reproduce inherited p<=12 count")
    need(len(charts) == 12378, "complete marked favorable atlas through period twelve")
    return rows, charts


def trim(poly):
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def divmod_monic(numerator, denominator):
    need(denominator[-1] == 1, "monic polynomial divisor")
    rem = list(numerator)
    quotient = [0]*max(1, len(rem)-len(denominator)+1)
    while len(rem) >= len(denominator) and any(rem):
        coefficient, shift = rem[-1], len(rem)-len(denominator)
        quotient[shift] = coefficient
        for i, a in enumerate(denominator):
            rem[i+shift] -= coefficient*a
        trim(rem)
    return trim(quotient), trim(rem)


@lru_cache(None)
def cyclotomic(n):
    poly = [-1]+[0]*(n-1)+[1]
    for d in range(1, n):
        if n % d == 0:
            poly, remainder = divmod_monic(poly, cyclotomic(d))
            need(not any(remainder), "exact cyclotomic factor division")
    return tuple(poly)


def Fourier_controls():
    # Exact polynomial remainders replace numerical complex evaluation.
    count = 0
    for T in range(2, 121):
        phi = cyclotomic(T)
        monomial = [0]*(T-1)+[1]
        need(any(divmod_monic(monomial, phi)[1]), "one final-position impulse survives a primitive mode")
        for p in range(1, T):
            if T % p:
                continue
            geometric = [int(j % p == 0) for j in range(T-p+1)]
            need(not any(divmod_monic(geometric, phi)[1]), "nontrivial repetition annihilates primitive mode")
            count += 1
    # Direct finite hostile probe of the repeated-word conclusion itself.
    words = [w for p in range(1, 6) for w in product((1, 2, 3), repeat=p) if primitive(w)]
    pairs = 0
    for i, word in enumerate(words):
        for other in words[i+1:]:
            T = lcm(len(word), len(other))
            if min(T//len(word), T//len(other)) == 1:
                T *= 2
            different = sum(word[j % len(word)] != other[j % len(other)] for j in range(T))
            need(different >= 2, "different repeated primitive words cannot differ only at the last place")
            pairs += 1
    cosine, sine = (1, 0, -1, 0), (0, 1, 0, -1)
    repeated, terminal_defect = (1, 1, 1, 1), (1, 1, 1, 2)
    dot = lambda a, b: sum(x*y for x, y in zip(a, b))
    need(dot(repeated, cosine) == dot(terminal_defect, cosine) == 0 and
         dot(repeated, sine) == 0 and dot(terminal_defect, sine) == -1,
         "fixed cosine component misses a quarter-turn final-position defect")
    return dict(cyclotomic_periods='2..120', exact_geometric_remainders=count,
                primitive_test_words=len(words), distinct_word_pairs=pairs)


def boundary_controls(compiler):
    first = compiler['cylinder']((1, 2), 1)
    second = compiler['cylinder']((1,), 2)
    need(first[:2] == second[:2] == (3, 4), "m=1 overlap hostile")
    need(compiler['verify']((1, 2), 1, 3) == compiler['verify']((1,), 2, 3) == 1,
         "same actual first-descent route 3->5->1")
    # Primitive is also essential: same marked primitive root counted twice.
    repeated = compiler['cylinder']((1, 1), 2)
    simple = compiler['cylinder']((1,), 4)
    need(repeated == simple, "nonprimitive duplicate chart")
    n, route, valuations = 7, [7], []
    while n >= 7:
        n, a = compiler['step'](n)
        route.append(n)
        valuations.append(a)
    need(route == [7, 11, 17, 13, 5] and valuations == [1, 1, 2, 3],
         "small convergent source outside every m>=2 periodic chart")
    for p in (1, 2):
        word = tuple(valuations[:p])
        need(any(word[j % p] != valuations[j] for j in range(3)),
             "all proper period divisors of first descent four fail before the final letter")
    # A countable family of tiny cylinders can contain every positive integer.
    tiny = [(n, n+10) for n in range(1, 101)]
    need(all(n % (1 << k) == r for n, (r, k) in enumerate(tiny, 1)), "every tested n lies in its tiny named cell")
    need(sum((F(1, 1 << k) for _, k in tiny), F()) < F(1, 1024),
         "small cylinder mass does not bound natural density of arbitrary countable unions")
    return dict(single_copy_overlap=dict(first_word=[1, 2], first_m=1, second_word=[1], second_m=2,
                                        cylinder=[3, 4], route=[3, 5, 1]),
                nonprimitive_overlap=dict(first_word=[1, 1], first_m=2, second_word=[1], second_m=4,
                                          cylinder=list(simple[:2])),
                all_period_missed_source=dict(source=7, first_descent_route=route,
                                              valuations=valuations, impossible_periods=[1, 2]),
                density_hostile='C(n,n+10), n>=1: Haar measure <=1/1024, positive-integer union all of N')


def non_descent_density(k):
    B = growth_bound(k)
    return sum((F(comb(A-1, k-1), 1 << (A+1)) for A in range(k, B+1)), F())


def natural_density_controls():
    rows = []
    for k in range(1, 61):
        B = growth_bound(k)
        carry = F(3**k, 2**k)-1
        gap = 1-F(3**k, 2**(B+1))
        threshold = carry/gap
        mass = non_descent_density(k)
        generating_bound = F(1, 2)*F(3, 5)**k*F(4, 3)**B
        need(mass <= generating_bound, "independent finite lower-tail sum vs generating-function bound")
        if k % 5 == 0:
            need(generating_bound < F(1, 2)*F(65536, 84375)**(k//5),
                 "explicit decaying five-step block bound")
        if k in (1, 2, 5, 10, 20, 24, 30, 40, 60):
            rows.append(dict(odd_steps=k, maximum_total_valuation=B,
                             uniform_source_threshold=str(threshold),
                             non_descent_upper_density=str(mass),
                             generating_function_upper=str(generating_bound)))
    # Independent finite exact affine-word and positive-source probes.
    controls = 0
    for k in range(1, 7):
        B = growth_bound(k)
        carry_bound = F(3**k, 2**k)-1
        threshold = carry_bound/(1-F(3**k, 2**(B+1)))
        start = threshold.numerator//threshold.denominator+1
        if start % 2 == 0:
            start += 1
        for n in range(start, start+512, 2):
            current, A, S = n, 0, 0
            for _ in range(k):
                z = 3*current+1
                a = (z & -z).bit_length()-1
                current = z >> a
                S = 3*S+2**A
                A += a
            need(F(S, 2**A) <= carry_bound, "uniform carry bound independent of the word")
            need(current*2**A == 3**k*n+S, "independent affine composition")
            need(A <= B or current < n, "uniform large-source contraction outside the finite valuation tail")
            controls += 1
    return rows, controls


def atlas_mass(charts, compiler):
    histogram, row_counts = Counter(), Counter()
    checks = 0
    max_bits = 0
    for word, P, Q, h, d in charts:
        need(Q*Q > h, "every m>=2 is admissible")
        Pm, Qm = P, Q
        for m in range(2, MAX_REPEAT+1):
            Pm *= P
            Qm *= Q
            t = 1
            while 2**t*(Qm-h) <= Pm-h:
                t += 1
            K = m*sum(word)+t
            histogram[K] += 1
            row_counts[len(word)] += 1
            if m in (2, 3):
                r, KK, tt = compiler['cylinder'](word, m)
                need((K, t) == (KK, tt), "independent mass budget agrees with compiler")
                compiler['verify'](word, m, r)
                max_bits = max(max_bits, r.bit_length())
                checks += 1
    lower = sum((F(count, 1 << K) for K, count in histogram.items()), F())
    repeat_tail = sum((F(1, (P-1)*P**MAX_REPEAT) for _, P, _, _, _ in charts), F())
    period_tail = F(3, 4*(3**(MAX_PERIOD+1)-1))
    need(F('0.11862865') < lower and lower+repeat_tail+period_tail < F('0.11862913'),
         "outward-rounded decimal interval checked against exact rational endpoints")
    need(sum(histogram.values()) == len(charts)*(MAX_REPEAT-1), "complete finite density universe")
    return dict(finite_density_lower=str(lower), repetition_tail_upper=str(repeat_tail),
                period_tail_upper=str(period_tail),
                full_atlas_density_upper=str(lower+repeat_tail+period_tail),
                cylinder_count=sum(histogram.values()),
                exponent_histogram=dict(sorted(histogram.items())),
                literal_first_descent_controls=checks, maximum_replayed_source_bits=max_bits)


def main():
    compiler = runpy.run_path(str(ROOT/'04-computation/experiments/rational_anchor_returns_20261004.py'))
    rows, charts = enumerate_charts()
    Fourier = Fourier_controls()
    hostile = boundary_controls(compiler)
    density_rows, carry_controls = natural_density_controls()
    measured = atlas_mass(charts, compiler)
    report = dict(status='PROVED universal m>=2 separation and natural-density existence; FINITE-EXACT enumerated atlas',
                  universe=dict(maximum_odd_period=MAX_PERIOD, minimum_repetition=2,
                                maximum_mass_repetition=MAX_REPEAT,
                                literal_repetitions=[2, 3], literal_lifts=[0]),
                  per_period_counts=rows, total_favorable_marked_charts=len(charts),
                  total_primitive_growth_necklaces=sum(r['primitive_growth_necklaces'] for r in rows),
                  Fourier_exact_controls=Fourier, boundary_hostiles=hostile,
                  uniform_non_descent_controls=carry_controls,
                  natural_density_tail_controls=density_rows, universal_atlas_mass=measured)
    print(json.dumps(report, indent=2))
    print('Approximate displays of the exact lower/upper endpoints:',
          format(float(F(measured['finite_density_lower'])), '.17g'),
          format(float(F(measured['full_atlas_density_upper'])), '.17g'))
    print('Certified outward-rounded interval: 0.11862865 < delta_atlas < 0.11862913')
    print('PASS: all checks remain active under optimized Python; no universal Collatz-coverage assertion.')


if __name__ == '__main__':
    main()
