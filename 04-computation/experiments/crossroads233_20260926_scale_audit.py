"""Independent controls for the two-cutoff P2 obstruction. Stdlib only."""
from fractions import Fraction
from array import array
from functools import lru_cache
from itertools import product
import json


def require(ok, message):
    if not ok:
        raise AssertionError(message)


def lawful(bits):
    end = len(bits)-1
    for i in range(2, end+1):
        if not i % 2:
            j = 3*i//2
            if j <= end and bits[i]+bits[j] != 1:
                return False
        else:
            j, k = (3*i-1)//2, (3*i+1)//2
            if j <= end and bits[j] > bits[i]:
                return False
            if k <= end and bits[i] > bits[k]:
                return False
    return True


def finite_cost(end, small):
    # Direct costs; no difference or local-toll recurrence is used here.
    costs = [(0, 0)]*(end+1)
    for i in range(end, 1, -1):
        candidates = []
        for bit in (0, 1):
            total = bit*(1+int(i <= small))
            if i % 2 == 0:
                j = 3*i//2
                if j <= end:
                    total += costs[j][1-bit]
            else:
                j, k = (3*i-1)//2, (3*i+1)//2
                if j <= end:
                    total += min(costs[j][:bit+1])
                if k <= end:
                    total += min(costs[k][bit:])
            candidates.append(total)
        costs[i] = tuple(candidates)
    return min(costs[2]) if end >= 2 else 0


@lru_cache(None)
def full_cost(depth, address):
    weight = 1 if depth < 2 else 2
    if not depth:
        return (0, 1)
    if address % 2 == 0:
        a, b = full_cost(depth-1, 3*address//2)
        return b, weight+a
    a, b = full_cost(depth-1, (3*address-1)//2)
    c, d = full_cost(depth-1, (3*address+1)//2)
    return a+min(c, d), weight+min(a, b)+d


def full_toll(depth, address):
    if not depth:
        return 0
    children = [3*address//2] if address % 2 == 0 else [(3*address-1)//2, (3*address+1)//2]
    return min(full_cost(depth, address))-sum(min(full_cost(depth-1, j)) for j in children)


def main():
    brute = 0
    for end in range(2, 16):
        small = 4*end//9
        best = 2*end
        for word in product((0, 1), repeat=end-1):
            bits = [0, 0, *word]
            brute += 1
            if lawful(bits):
                best = min(best, sum(bits)+sum(bits[:small+1]))
        require(best == finite_cost(end, small), 'weighted optimum versus all assignments')
    require(finite_cost(2, 0) == 0 and finite_cost(5, 0) == 1 and finite_cost(5, 2) == 2,
            'minimal two-cutoff incompatible optimum')
    print(json.dumps({'boolean_assignments': brute, 'minimal_cutoffs': [2, 5],
                      'separate_optima': [0, 1], 'joint_optimum': 2}))
    sums = []
    for depth in range(1, 10):
        total = sum(full_toll(depth, i) for i in range(2**depth))
        sums.append(total)
        for i in range(2**depth):
            a, b = full_cost(depth, i)
            require(abs(b-a) <= 2*(depth+1), 'depth bound')
            require(full_toll(depth, i) >= 0, 'nonnegative toll')
            require(full_toll(depth, i) == full_toll(depth, i+2**depth), 'residue period')
    print(json.dumps({'independent_full_tree_coefficients': sums, 'max_depth': 9}))
    for end in (233, 1000, 10000, 100000):
        small = 4*end//9
        value = finite_cost(end, small)
        print(json.dumps({'large_cutoff': end, 'small_cutoff': small,
                          'joint_optimum': value, 'normalized': float(Fraction(value, end+small))}))
    # Exact certificate supplied as integers; its coefficient generation is
    # independently checked above and compared with the lane's output by root.
    alpha_upper = Fraction(3089623223, 10460353203)
    beta20 = Fraction(4538724347, 15109399071)
    oscillating_upper_floor = (13*beta20-9*alpha_upper)/4
    require(beta20 > alpha_upper, 'strict separation')
    require(oscillating_upper_floor == Fraction(13417603, 43046721), 'oscillation floor')
    require(oscillating_upper_floor-alpha_upper == Fraction(170854306, 10460353203),
            'oscillation gap')
    print(json.dumps({'beta20': str(beta20), 'alpha_upper': str(alpha_upper),
                      'forced_limsup': str(oscillating_upper_floor),
                      'forced_oscillation': str(oscillating_upper_floor-alpha_upper)}))
    # Independent full-cost arrays, never using the lane's difference/toll
    # recurrence. Their total has the same telescoping fringe certificate.
    kernel = [128, 192, 288, 432, 648, 972, 1458, 2187]
    normalizer = sum((Fraction(a)*Fraction(2, 3)**j for j, a in enumerate(kernel)), Fraction(0))
    require(normalizer == 1024, 'eight-cutoff normalization')
    zeros, ones = array('q', [0]), array('q', [kernel[0]])
    for depth in range(1, 21):
        modulus = len(zeros)
        weight = sum(kernel[:depth+1])
        next_zero, next_one = array('q'), array('q')
        for address in range(2*modulus):
            if address % 2 == 0:
                child = (3*address//2) % modulus
                zero, one = ones[child], weight+zeros[child]
            else:
                left = ((3*address-1)//2) % modulus
                right = ((3*address+1)//2) % modulus
                zero = zeros[left]+min(zeros[right], ones[right])
                one = weight+min(zeros[left], ones[left])+ones[right]
            next_zero.append(zero); next_one.append(one)
        zeros, ones = next_zero, next_one
    aggregate = sum(min(a, b) for a, b in zip(zeros, ones))
    best = Fraction(aggregate, 3**21*1024)
    require(aggregate == 3286041555236, 'full-cost aggregate certificate')
    require(best == Fraction(821510388809, 2677850419968), 'eight-cutoff lower bound')
    require(best > beta20, 'eight-cutoff strict improvement')
    print(json.dumps({'independent_full_cost_depth': 20, 'residue_types': len(zeros),
                      'aggregate': aggregate, 'kernel_normalizer': str(normalizer),
                      'upper_density_lower': str(best), 'decimal': float(best)}))
    print('ALL INDEPENDENT SCALE CHECKS PASSED')


if __name__ == '__main__':
    main()
