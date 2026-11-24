"""Self-converse presentation mass via the affine events of reversing involutions."""
from argparse import ArgumentParser
from itertools import combinations
from math import comb, factorial


def involution_count(n):
    half = n//2
    return factorial(n)//(2**half * factorial(half))


def marked_from_runs(n):
    """Place disjoint paired even runs and, optionally, a centered even run."""
    total = 0
    for pairs in range(n//6 + 1):
        for half_sum in range(pairs, n//4 + 1):
            if pairs == 0 and half_sum:
                continue
            ways = comb(half_sum-1, pairs-1) if pairs else 1
            for center_half in range(1 + (n//2 if n%2 else 0)):
                centered = int(center_half > 0)
                runs = 2*pairs + centered
                vertices = 4*half_sum + 2*pairs + 2*center_half + centered
                if vertices > n:
                    continue
                placements = factorial(runs)*comb(n-vertices+runs, runs)//factorial(pairs)
                total += ways*placements*involution_count(n-vertices)
    orbit_count = (comb(n,2) + n//2)//2
    return total << (orbit_count-(n-1))


def matchings(vertices, fixed):
    if not vertices:
        if fixed == 0:
            yield []
        return
    a, *others = vertices
    if fixed:
        for pairs in matchings(others, 0):
            yield [(a,a), *pairs]
    for j,b in enumerate(others):
        for pairs in matchings(others[:j]+others[j+1:], fixed):
            yield [(a,b), *pairs]


def involutions(n):
    for pairs in matchings(list(range(n)), n%2):
        pi = [0]*n
        for a,b in pairs:
            pi[a],pi[b] = b,a
        yield pi


def affine_event(n, pi, tiles):
    if any(pi[i+1] == pi[i]+1 for i in range(n-1)):
        return None
    seed, generators = 0, []
    for i,j in combinations(range(n), 2):
        a,b = sorted((pi[i],pi[j]))
        if (a,b) < (i,j):
            continue
        # If images are in ascending order, this anti-map reverses the bit.
        flip = int(pi[i] < pi[j])
        first = tiles.get((i,j), 0)
        other = tiles.get((a,b), 0)
        if first and other:
            generators.append(first | other)
            if flip:
                seed |= other
        elif flip:
            seed |= first | other
    return seed, generators


def small_inclusion_exclusion(events):
    # Intersections as sets for the very small independent IE check.
    total = 0
    for size in range(1,len(events)+1):
        sign = 1 if size%2 else -1
        total += sign*sum(len(set.intersection(*parts)) for parts in combinations(events, size))
    return total


def census(n):
    tiles = {edge:1<<i for i,edge in enumerate(
        (e for e in combinations(range(n),2) if e[1] > e[0]+1))}
    space = 1<<len(tiles)
    seen = bytearray((space+7)//8)
    twice = bytearray(len(seen))
    marked, count, multiple, compatible = 0, 0, 0, 0
    small = []
    for pi in involutions(n):
        event = affine_event(n, pi, tiles)
        if event is None:
            continue
        compatible += 1
        value, generators = event
        size = 1<<len(generators)
        marked += size
        board_set = set()
        for index in range(size):
            if index:
                value ^= generators[(index & -index).bit_length()-1]
            block, bit = divmod(value,8)
            mask = 1<<bit
            if not seen[block]&mask:
                count += 1
                seen[block] |= mask
            elif not twice[block]&mask:
                multiple += 1
                twice[block] |= mask
            if n <= 5:
                board_set.add(value)
        if n<=5:
            small.append(board_set)
    if n<=5:
        assert small_inclusion_exclusion(small) == count
    assert marked == marked_from_runs(n)
    reflection = affine_event(n, list(range(n-1,-1,-1)), tiles)
    blue = 1<<len(reflection[1])
    return space, compatible, marked, count, multiple, blue, space-count


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--through", type=int, choices=range(2,10), default=9)
    args = parser.parse_args()
    print("n total compatible marked self-converse multiple blue doubled")
    for n in range(2,args.through+1):
        print(n,*census(n),flush=True)
