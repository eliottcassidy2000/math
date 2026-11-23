"""Joint Pareto frontier of the two E_e projected modular sieves, all odd exponents."""
from importlib import import_module
from itertools import product
from fractions import Fraction

phase = import_module("square-directions")


def partly_fixed(n, private_prime, target):
    primes = (2,5,29,private_prime)
    tables = [phase.local_counts(p,phase.chart(n),target%p) for p in primes]
    sizes = phase.all_sizes(tables)
    private_signatures = set(tables[-1].values())
    out = {}
    for shared in product(*(t.keys() for t in tables[:3])):
        signature = tuple(t[k] for t,k in zip(tables[:3],shared))
        out[shared] = min(sizes[signature+(last,)] for last in private_signatures)
    return out


def frontier(target):
    left = partly_fixed(4,17,target)
    right = partly_fixed(8,257,target)
    pairs = sorted({(left[s],right[s]) for s in left})
    result = []
    for a,b in pairs:
        if not result or b < result[-1][1]:
            result.append((a,b))
    return result


if __name__ == "__main__":
    counts = {1:0,2:0}
    groups = {}
    for e in range(3,115,2):
        result=frontier(2**(e-1)*(2**e-1)-1)
        counts[len(result)] += 1
        if len(result)==2:
            (a,v),(u,b) = result
            assert Fraction(u,a)==Fraction(v,b)==Fraction(15,14)
            assert e%28 in (7,11,13,23)
            groups.setdefault(tuple(result),[]).append(e%112)
        else:
            assert e%28 not in (7,11,13,23)
    print("Numbers of one- and two-point frontiers:", counts)
    for points,exponents in sorted(groups.items()):
        print(points,sorted(exponents))
