"""Integer certificates for the shifted-square intersection of the rigid block bank."""
from functools import lru_cache
from itertools import permutations, product
from math import gcd, isqrt
from pathlib import Path
import subprocess
import sys
import tempfile

PREVIOUS = Path(__file__).resolve().parent.parent / "chanpak_tiling_viewpoint_20261123"
sys.path.insert(0, str(PREVIOUS))
from precedence_bank import SEEDS


def first_norm_two(discriminant):
    root = isqrt(discriminant)
    term = root
    shift, divisor = 0, 1
    p_previous, p = 0, 1
    q_previous, q = 1, 0
    norms = []
    for _ in range(26):
        p_previous, p = p, term*p + p_previous
        q_previous, q = q, term*q + q_previous
        norm = p*p - discriminant*q*q
        norms.append(norm)
        if norm == 2:
            return p, q, norms
        shift = divisor*term - shift
        divisor = (discriminant-shift*shift)//divisor
        term = (root+shift)//divisor
    raise AssertionError("The stated certificate window did not contain norm two")


def ordinal(left, right):
    split = len(left)
    suffix = ((1 << len(right)) - 1) << split
    return tuple(row | suffix for row in left) + tuple(row << split for row in right)


def height(rows):
    n = len(rows)
    full = (1 << n) - 1
    counts = [[0]*n for _ in range(1 << n)]
    for vertex in range(n):
        counts[1 << vertex][vertex] = 1
    for used in range(1, full):
        absent = full ^ used
        for last, multiplicity in enumerate(counts[used]):
            if not multiplicity:
                continue
            choices = absent & rows[last]
            while choices:
                bit = choices & -choices
                choices ^= bit
                counts[used | bit][bit.bit_length()-1] += multiplicity
    return sum(counts[full])


def unit_order(d, x, p):
    # The norm-two base has y=1 in the three applications.
    g, h = (x*x-1) % p, x % p
    u, v = 1, 0
    for n in range(1, p*p):
        u, v = (u*g+d*v*h) % p, (u*h+v*g) % p
        if (u,v) == (1,0):
            return n
    raise AssertionError("Norm-one element did not return")


def pell_denominator(d, x, n):
    assert n > 0 and n % 2 == 1 and x*x-d == 2
    u, v = x, 1
    for _ in range((n-1)//2):
        u, v = u*(x*x-1)+d*v*x, u*x+v*(x*x-1)
    return v


def augmented_inventory():
    allowed_small = [h for h in range(1,26,2) if h not in (7,21)]
    inert = (3,5,11,13,19,37)
    assert [h for h in allowed_small if all(h % p for p in inert)] == [1,17,23]
    assert first_norm_two(391)[:2] == (2709,137)
    allowed_product = 17*23*751
    assert gcd(allowed_product, 271*137) == 1
    for d, x, p, order, r, forbidden in (
        (23,5,17,9,3,7),
        (23,5,751,125,5,2351),
        (12767,113,23,11,11,857),
    ):
        assert unit_order(d,x,p) == order and order % r == 0
        assert pell_denominator(d,x,r) % forbidden == 0
        assert forbidden > 1 and gcd(forbidden,allowed_product) == 1


def raw_small_inventory():
    for d, x, p, r, valuation in ((7,3,17,3,1), (7,3,23,12,None),
                                  (119,11,23,11,1), (23,5,7,3,2),
                                  (23,5,17,9,1)):
        assert unit_order(d,x,p) == r
        if valuation is not None:
            y = pell_denominator(d,x,r)
            assert y % p**valuation == 0 and y % p**(valuation+1) != 0
    assert pell_denominator(23,5,5) == 2351
    for s, h in ((3,7),(5,23),(11,7*17),(45,7*17**2),(235,23*7**4)):
        assert s*s-2 == h


def physical_seams():
    for root, rows, pair in (
        (5, [2,16,11,35,13,23], (0,4)),
        (11, [246,124,120,241,224,192,128,6], (4,7)),
        (45, [476,113,66,310,420,389,56,78,198], (6,8)),
        (235, [1766,500,952,691,1601,1872,1164,304,601,194,910], (0,10)),
        (113, [594,360,539,209,898,349,916,807,525,42], (4,9))
    ):
        assert height(rows) == root*root-2
        u,v = pair
        rows[u] ^= 1 << v
        rows[v] ^= 1 << u
        assert height(rows) == root*root
        print("minimum-order seam",root*root-2,root*root,"order",len(rows))


@lru_cache(None)
def canonical(rows, marked=()):
    # Isomorphism-invariant refinement. Exhaust the permutations within every tie.
    colors = [(v in marked, rows[v].bit_count()) for v in range(len(rows))]
    while True:
        palette = sorted(set(colors))
        cells = [[v for v, color in enumerate(colors) if color == key] for key in palette]
        refined = [(palette.index(colors[v]),
                    *(sum((row >> w) & 1 for w in cell) for cell in cells))
                   for v, row in enumerate(rows)]
        if len(set(refined)) == len(cells):
            break
        colors = refined
    keys = []
    for pieces in product(*(permutations(cell) for cell in cells)):
        order = sum(pieces, ())
        matrix = tuple(sum(((rows[v] >> w) & 1) << i for i,w in enumerate(order))
                       for v in order)
        keys.append((tuple(v in marked for v in order), matrix))
    return min(keys)


def converse(rows):
    full = (1 << len(rows))-1
    return tuple(full ^ row ^ (1 << v) for v,row in enumerate(rows))


def census():
    here = Path(__file__).resolve().parent
    inherited = here.parent / "chanpak_arc_fiber_grammar_20261122" / "birth_cut_cover.cpp"
    with tempfile.TemporaryDirectory() as temp:
        obj, exe = Path(temp)/"heights.o", Path(temp)/"scan"
        subprocess.run(["c++","-std=c++17","-O3","-Dmain=birth_cut_cover_main","-c",
                        str(inherited),"-o",str(obj)],check=True)
        subprocess.run(["c++","-std=c++17","-O3",str(here/"bank_scan.cpp"),
                        str(obj),"-o",str(exe)],check=True)
        for n, h, expected in ((5,23,(12,4,2,4,4,2)), (9,12767,(48,24,22,24,24,12))):
            with subprocess.Popen(["gentourng","-q",str(n)],stdout=subprocess.PIPE) as gen:
                records = subprocess.check_output([str(exe),str(n),str(h)],stdin=gen.stdout,
                                                   text=True).splitlines()
                assert gen.wait() == 0
            lowers, uppers, pairs, marked, merged = (set() for _ in range(5))
            for line in records:
                *r, i, j = map(int,line.split())
                t = r.copy()
                t[i] ^= 1 << j
                t[j] ^= 1 << i
                r, t = tuple(r), tuple(t)
                assert (height(r),height(t)) == (h,h+2)
                low, high = canonical(r), canonical(t)
                assert low != canonical(converse(r)) and high != canonical(converse(t))
                m = canonical(r,(i,j))
                lowers.add(low)
                uppers.add(high)
                pairs.add((low,high))
                marked.add(m)
                merged.add(min(m,canonical(converse(r),(i,j))))
            observed = tuple(map(len,(records,lowers,uppers,pairs,marked,merged)))
            assert observed == expected
            print("minimum-order marked census",h,observed)


def main():
    for factor, prime in ((5,5), (37,37), (77,11), (129,3)):
        assert factor % prime == 0 and 2 not in {x*x % prime for x in range(prime)}
    x, y, norms = first_norm_two(751)
    assert (x,y) == (2700614460969,98546821297)
    assert len(norms) == 26 and y % 271 == 0 and gcd(271,17*751) == 1
    assert first_norm_two(12767)[:2] == (113,1)
    for a,b in ((17,751),(751,17)):
        rows = ordinal(SEEDS[a], SEEDS[b])
        assert len(rows) == 15 and height(rows) == 113**2-2 == a*b
    augmented_inventory()
    raw_small_inventory()
    physical_seams()
    print("D=751 convergent norms:", *norms)
    print("Exactly one positive shifted square in the six-block bank: 113^2-2=17*751=12767.")


if __name__ == "__main__":
    main()
    if sys.argv[1:] == ["--census"]:
        census()
