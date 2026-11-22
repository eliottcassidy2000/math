"""Sharp one-sided Pell approximation under the two E61 projection minima."""
from importlib import import_module
from itertools import product
from math import gcd, lcm, prod

phase = import_module("square-directions")
PRIMES = (2,5,17,29,257)
MODULUS = prod(PRIMES)


def local_winners(n, primes):
    tables = [phase.local_counts(p,phase.chart(n),phase.TARGET%p) for p in primes]
    sizes = phase.all_sizes(tables)
    optimum = min(sizes.values())
    winning_tuples = [ds for ds in product(*(t.keys() for t in tables))
                      if sizes[tuple(t[d] for t,d in zip(tables,ds))]==optimum]
    factors = {p:{ds[i] for ds in winning_tuples} for i,p in enumerate(primes)}
    # Equality of cardinalities proves the full set is the Cartesian product.
    assert prod(map(len,factors.values())) == len(winning_tuples)
    return factors


def advance(a,b,modulus):
    return (3*a+2*b)%modulus,(4*a+3*b)%modulus


def direction_period(p):
    slopes=[(1,b) for b in range(p)]+[(0,1)]
    current = slopes
    for k in range(1,p+2):
        current = [phase.normalized(*advance(a,b,p),p) for a,b in current]
        if current == slopes:
            return k
    raise AssertionError(p)


def main():
    four = local_winners(4,(2,5,17,29))
    eight = local_winners(8,(2,5,29,257))
    allowed = {p:set.intersection(*(t[p] for t in (four,eight) if p in t)) for p in PRIMES}
    assert allowed[2] == allowed[5] == {(1,0)}
    assert allowed[17] == {(1,s) for s in (1,2,4,8,9,13,15,16)}
    periods=[direction_period(p) for p in PRIMES]
    assert periods == [1,3,4,5,32]
    period=lcm(*periods)
    for seed in [(1,2),(1,0),(7,6),(9,10)]:
        a,b=seed
        hits=[]
        for index in range(period):
            if all(phase.normalized(a,b,p) in allowed[p] for p in PRIMES):
                hits.append(index)
            a,b=advance(a,b,MODULUS)
        print(seed,"winning indices per",period,len(hits),"first",hits[:5])
        if seed == (1,2):
            assert len(hits)==96 and hits[0]==1
        if seed == (1,0):
            assert hits==[]
        if seed == (7,6):
            assert hits[0]==14
    a,b=7,6
    for _ in range(14):
        a,b=3*a+2*b,4*a+3*b
    assert (a,b)==(293465818991,415023341310) and 2*a*a-b*b==62
    assert gcd(a,b)==1 and all(phase.normalized(a,b,p) in allowed[p] for p in PRIMES)
    print("upper Pell certificate",a,b,"defect",2*a*a-b*b)


if __name__ == "__main__":
    main()
