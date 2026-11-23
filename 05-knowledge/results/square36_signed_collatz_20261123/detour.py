"""Build four intermediate tournaments in the two shortest E61 square-sum detours."""
from functools import cache
from itertools import combinations
from math import isqrt, prod
from pathlib import Path
import sys

RESULTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RESULTS/"chanpak_prime_frontier_selector_20261118"))
sys.path.insert(0, str(RESULTS/"chanpak_decimal_kernel_20261122"))
from join_loads import analyze
from e15_pair import literal_count, star_rows

JOINS = {
    3:(1,2), 5:(2,2), 31:(2,11), 37:(2,13), 47:(1,25), 49:(5,5),
    379:(7,35), 5051:(58,74), 7243:(32,221), 17107:(78,169),
    27127:(140,182), 48119:(42,891), 77551:(162,329), 197573:(204,837),
    411883:(74,4085), 493657:(181,1919), 732541:(539,983), 1391029:(124,10211),
    2129389:(801,1445), 3085921:(756,3481), 5288203:(54,72577),
    10663217:(2935,3536), 10803413:(302,27585), 47311571:(1720,25713),
    65145583:(1307,31145),
}
STARS = {83:(10,5), 109:(1,10), 1597:(11,83), 84121:(4,7593)}

# The rows give a divisor of E-1, factors of u, s-u and s+u, respectively.
PATHS = (
    (3187801395226305, (5,83,7243,2129389,65145583),
     (49,47,48119,10663217,493657), (109,197573,1391029,47311571)),
    (73039056357735, (3,5,37,5051,27127,77551,3085921),
     (31,31,379,17107,411883,5288203), (47,1597,84121,732541,10803413)),
)


@cache
def block(label):
    if label in JOINS:
        proof = analyze(*JOINS[label])
        assert proof['H'] == label
        return tuple(proof['rows'])
    mask, address = STARS[label]
    host = star_rows((2,8,3,5), mask, 0, address)
    assert literal_count(host, [4]+[1]*address.bit_length()) == label
    return tuple(host)


def tournament(factors):
    rows = []
    for f in factors:
        source = block(f)
        n = len(rows)
        new_vertices = ((1<<len(source))-1) << n
        rows = [r|new_vertices for r in rows] + [r<<n for r in source]
    for i,r in enumerate(rows):
        assert r>>len(rows) == 0 and ((r>>i)&1) == 0
        assert all((((r>>j)^(rows[j]>>i))&1) == 1 for j in range(i))
    return rows


def main():
    e = (1<<60)*((1<<61)-1)
    a, b = e-1, e+1
    paths = []
    for divisor, uf, left, right in PATHS:
        u = prod(uf)
        root1 = a//divisor-u
        x, y = u*u, prod(left)*prod(right)
        s = prod(left)+u
        assert prod(right) == s+u and root1*root1 == a+x and s*s == x+y
        root3 = isqrt(y+b)
        assert root3**2 == y+b
        order1 = len(tournament(uf+uf))
        order2 = len(tournament(left+right))
        print("intermediate heights", x,y, "orders",order1,order2, "edge roots",root1,s,root3)
        paths.append((x,y))
    labels = [a, *paths[0], b, *reversed(paths[1])]
    assert len(set(labels)) == 6
    for i,j in combinations(range(6),2):
        total = labels[i]+labels[j]
        assert (isqrt(total)**2 == total) == ((j-i) in (1,5))
    print("the six labels induce exactly the six cycle edges")


if __name__ == "__main__":
    main()
