"""Independent direct-trajectory audit of the critical-tube price argument.

Does not import another experiment. Finite verification supplements the
all-height interval-counting and cyclic-rotation proofs, not vice versa.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from math import comb


def require(ok, label):
    if not ok:
        raise AssertionError(label)


def step(n):
    return (3*n+1)//2 if n & 1 else n//2


def main():
    for length in (6, 8, 10, 12):
        for cap in (3, 9, 27):
            bins = defaultdict(set)
            served = defaultdict(set)
            count = 0
            for n in range(1, 4*(1 << length)):
                x = n
                e = 0
                slope = Fraction(1)
                carry = Fraction(0)
                rows = []
                valid = True
                for k in range(length):
                    rows.append((k, e, x, slope, carry))
                    if x & 1:
                        carry += Fraction(1, 3)/slope
                        e += 1
                        slope *= Fraction(3, 2)
                    else:
                        slope /= 2
                    x = step(x)
                    require(x == slope*(n+carry), 'actual affine source identity')
                    if slope < 1 or slope > cap:
                        valid = False
                        break
                if not valid:
                    continue
                count += 1
                for k, e, v, slope, carry in rows:
                    require(Fraction(v)/slope-Fraction(k, 3) <= n <= Fraction(v)/slope,
                            'actual-source interval')
                    bins[k, e, v].add(n)
                    served[v].add(n)
            for (k, e, v), sources in bins.items():
                require(len(sources) <= k//3+1, 'integer interval occupancy')
            ell = 0
            while 3**(ell+1) <= cap:
                ell += 1
            bound = (ell+1)*sum(k//3+1 for k in range(length))
            require(all(len(sources) <= bound for sources in served.values()), 'hub capacity')
            pairs = defaultdict(set)
            for v, sources in served.items():
                pairs[(v+1)//2].update(sources)
            require(all(len(sources) <= 2*bound for sources in pairs.values()), 'both pair members')
            print('tube', length, cap, 'sources', count,
                  'max_actual_hub', max(map(len, served.values()), default=0), 'M', bound)
    for block_length in range(2, 17):
        ones = 0
        while 3**ones < 2**block_length:
            ones += 1
        good = set()
        for positions in combinations(range(block_length), ones):
            word = [int(j in positions) for j in range(block_length)]
            slopes = [Fraction(1)]
            slope = Fraction(1)
            for bit in word:
                slope *= Fraction(3 if bit else 1, 2)
                slopes.append(slope)
            cut = min(range(block_length), key=lambda j: slopes[j])
            rotated = word[cut:]+word[:cut]
            slope = Fraction(1)
            for bit in rotated:
                slope *= Fraction(3 if bit else 1, 2)
                require(slope >= 1, 'cyclic-minimum prefix')
            require(slope < 3, 'block excess')
            good.add(tuple(rotated))
        require(len(good)*block_length >= comb(block_length, ones), 'rotation-fibre count')
        print('rotation', block_length, ones, 'distinct_good', len(good),
              'binomial', comb(block_length, ones))
    print('PASS independent actual-source capacity and critical-block rotation audit')


if __name__ == '__main__':
    main()
