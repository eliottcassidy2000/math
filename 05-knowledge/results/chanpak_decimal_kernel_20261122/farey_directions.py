"""Invert the lower chart on every positive ray with bounded primitive product."""
from argparse import ArgumentParser
from itertools import combinations
from math import gcd, isqrt

SUPERSINGULAR = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71)

from fermat_receiver import TARGET, address_height
from hyperbola_controls import chart


def direction_scale(data, target, a, b):
    """Return the unique positive scale that could attain the height."""
    h, q, c, t, e, r, k = data
    assert min(a, b) >= 1 and q > h and q+c+t > e+r
    if target < k:
        return 0
    quadratic = q*a*b
    linear = t*a+c*b
    return (isqrt(linear*linear+4*quadratic*(target-k))-linear)//(2*quadratic)


def search(data, target, product_limit):
    h, q, c, t, e, r, k = data
    proofs = near = 0
    hits = []
    differences = []
    controls = {}
    for a in range(1, product_limit+1):
        for b in range(1, 1+product_limit//a):
            if gcd(a, b) != 1:
                continue
            proofs += 1
            scale = direction_scale(data, target, a, b)
            if not scale:
                continue
            x, y = a*scale, b*scale
            residual = target-k-q*x*y-t*x-c*y-r*(x&1)-e*(y&1)
            if 0 <= residual <= h*(x+y):
                difference = address_height(data, x, y)-target
                differences.append(difference)
                for prime in (17, 47):
                    if difference % prime == 0:
                        controls.setdefault(prime, (a,b,scale,difference))
                if residual % (2*h) == 0:
                    near += 1
                if difference == 0:
                    hits.append((x, y))
    return proofs, near, hits, differences, controls


def least_observers(differences):
    for size in range(len(SUPERSINGULAR)+1):
        good = []
        for primes in combinations(SUPERSINGULAR, size):
            modulus = 1
            for prime in primes:
                modulus *= prime
            if all(value % modulus for value in differences):
                good.append(primes)
        if good:
            return good
    return []


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("product_limit", type=int)
    parser.add_argument("spines", type=int, nargs="+")
    args = parser.parse_args()
    combined = []
    for n in args.spines:
        data = chart(n)
        for a, b, scale in ((1, 3, 11), (7, 10, 641), (8, 5, 2**40+101)):
            value = address_height(data, scale*a, scale*b)
            assert direction_scale(data, value, a, b) == scale
        total, near, hits, differences, controls = search(data, TARGET, args.product_limit)
        gbound = isqrt((TARGET-data[-1])//(data[1]*(args.product_limit+1)))
        print("n", n, "primitive directions", total, "residual candidates", near,
              "hits", hits, "gcd bound on other pairs", gbound, flush=True)
        observers = least_observers(differences)
        print("in absolute envelope", len(differences), "least observer sets", len(observers),
              "first", observers[:1], flush=True)
        print("separate local witnesses", controls, flush=True)
        combined.extend(differences)
    observers = least_observers(combined)
    print("combined", len(combined), "least observer sets", len(observers), "observers", observers)
