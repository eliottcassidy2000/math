"""Small independent exact controls; does not repeat the incoming 2^32 scan."""
from fractions import Fraction as Q
from itertools import product


def require(test, label):
    if not test:
        raise RuntimeError(label)


cases = empty = 0
for k in range(1, 11):
    for W in (Q(3, 2), Q(2), Q(3), Q(7), Q(100)):
        hits = []
        residual = Q(0)
        for word in product((0, 1), repeat=k):
            value = Q(1)
            first = None
            for bit in word:
                value *= Q(3, 2) if bit else Q(1, 2)
                if first is None and value >= W:
                    first = value
            if first is None:
                residual += value
            else:
                hits.append(first)
        eps = residual / (2 ** k)
        require(sum(hits, Q(0)) / (2 ** k) == 1 - eps, 'undivided stopping identity')
        if hits:
            mean = sum(hits, Q(0)) / len(hits)
            probability = Q(len(hits), 2 ** k)
            require(W <= mean < 3 * W / 2, 'overshoot bound')
            require(probability == (1 - eps) / mean, 'conditional identity on nonempty event')
        else:
            require(eps == 1, 'empty hit event')
            empty += 1
        if k == 1 and W == 2:
            print('Hostile k=1 W=2: probability=0 eps=1 conditional_mean=undefined')
        cases += 1
print('Exact finite stopping cases', cases, 'including empty hit cases', empty)

value = Q(1)
ones = 0
for j in range(1, 2001):
    bit = value < 2
    ones += bit
    value *= Q(3, 2) if bit else Q(1, 2)
    require(1 <= value < 3, 'greedy threshold3 invariant')
    require(value == Q(3 ** ones, 2 ** j), 'actual parity multiplier')
print('Greedy rotation invariant checked for2000 prefixes with no threshold3 hit')
print('PASS: finite controls only; infinite nonattainment uses the written irrational-rotation proof')
