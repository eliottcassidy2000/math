"""Optimize the projected fixed-parity phase on all primitive rational normals."""
from argparse import ArgumentParser
from functools import lru_cache
from itertools import product
from math import gcd, prod
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "chanpak_decimal_kernel_20261122"))
from hyperbola_controls import chart
from fermat_receiver import TARGET


@lru_cache(None)
def local_counts(p, data, target):
    """For each projective normal store every nonempty intersection of parity masks."""
    _, q, c, t, e, r, k = data
    directions = [(1, b) for b in range(p)] + [(0, 1)]
    points = []
    for i, j in product(range(2), repeat=2):
        if p == 2:
            points.append([(i, j)])
        else:
            points.append([(x, y) for x in range(p) for y in range(p)
                           if (k+t*x+c*y+q*x*y+r*i+e*j-target) % p == 0])
    summaries = {}
    for a, b in directions:
        masks = []
        for fibre in points:
            mask = 0
            for x, y in fibre:
                mask |= 1 << ((a*x+b*y) % p)
            masks.append(mask)
        counts = []
        for subset in range(1, 16):
            common = (1 << p) - 1
            for index in range(4):
                if subset & (1 << index):
                    common &= masks[index]
            counts.append(common.bit_count())
        summaries[a, b] = tuple(counts)
    return summaries


def normalized(a, b, p):
    return (1, b*pow(a, -1, p) % p) if a%p else (0, 1)


def all_sizes(tables):
    signs = [1 if subset.bit_count() % 2 else -1 for subset in range(1, 16)]
    signatures = [set(table.values()) for table in tables]
    return {choice: sum(sign*prod(values) for sign, values in zip(signs, zip(*choice)))
            for choice in product(*signatures)}


def euclid_exponents():
    common = set(range(1,112,2))
    for n, primes in [(4, (2, 5, 17, 29)), (8, (2, 5, 29, 257))]:
        data = chart(n)
        winners = set()
        for e in range(1, 112, 2):
            exponent = e if e > 1 else 113
            target = (1 << (exponent-1))*((1 << exponent)-1)-1
            tables = [local_counts(p, data, target % p) for p in primes]
            sizes = all_sizes(tables)
            sig = tuple(t[normalized(7,10,p)] for t,p in zip(tables,primes))
            if sizes[sig] == min(sizes.values()):
                winners.add(e)
        print(n, "at", len(winners), "odd exponent classes")
        common &= winners
    print("simultaneously optimal exponent residues mod 112:", sorted(common))


def main():
    for n, primes in [(4, (2, 5, 17, 29)), (8, (2, 5, 29, 257))]:
        tables = [local_counts(p, chart(n), TARGET % p) for p in primes]
        # Normals giving equal intersection counts have the same global objective.
        sizes = all_sizes(tables)
        best = min(sizes.values())
        min_dirs = 0
        for keys in product(*(table.keys() for table in tables)):
            sig = tuple(table[key] for table, key in zip(tables, keys))
            min_dirs += sizes[sig] == best
        witness = None
        for largest in range(1, 301):
            for a, b in (pair for d in range(1,largest+1)
                         for pair in ((largest,d),(d,largest))):
                if gcd(a, b) != 1:
                    continue
                sig = tuple(t[normalized(a,b,p)] for p,t in zip(primes,tables))
                if sizes[sig] == best:
                    witness = a, b
                    break
            if witness:
                break
        assert witness is not None
        print(n, "classes", prod(p+1 for p in primes), "minimum", best,
              "of", prod(primes), "attained by", min_dirs, "classes; weights", witness)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--exponents", action="store_true")
    if parser.parse_args().exponents:
        euclid_exponents()
    else:
        main()
