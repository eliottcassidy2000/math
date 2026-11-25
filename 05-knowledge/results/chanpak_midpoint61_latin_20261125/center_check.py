"""Replay the exact midpoint, packing window, and target-directed ternary word."""
from fractions import Fraction
from pathlib import Path
import sys

base = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(base / "chanpak_gap_star_ports_20261122"))
sys.path.insert(0, str(base / "chanpak_square25_dilation_20261123"))
from replay import flipped, strong, subpaths
from selector import skeleton

LOW = (246,252,152,241,96,196,132,16)
HIGH = flipped(LOW,1,5)


def inflate(outer, m, k):
    sites = [(i, rank) for i in range(8)
             for rank in range(m if i == 0 else k if i == 7 else 1)]
    return tuple(sum(1 << y for y,(j,s) in enumerate(sites)
                     if (i == j and r < s) or (i != j and outer[i] >> j & 1))
                 for i,r in sites)


def check():
    for source, weights in [(LOW,(26,-1,-5,-1)),(HIGH,(38,-1,-1,1))]:
        a,b,c,d = weights
        words = skeleton(source,(0,7))
        assert words == {(1,1):4*a+2*b+2*c+d, (2,1):2*a+b, (1,2):2*a+c, (2,2):a}
        assert strong(source)
        for m,k in [(1,1),(2,3),(3,2),(4,4)]:
            graph = inflate(source,m,k)
            assert strong(graph)
            assert subpaths(graph)[-1] == (a << (m+k)) + (b << m) + (c << k) + d
    for p in range(4,128):
        x = 1 << (p-3)
        mid = (1 << (p-1))*((1 << p)-1)
        assert (26*x*x-6*x-1 + 38*x*x-2*x+1)//2 == mid
        assert 38*x*x-2*x+1-mid == 6*x*x+2*x+1
    total,s = 3250,857
    mean = Fraction(57*s,total)
    variance = Fraction(64*s*(total-s),total*total)
    assert total*variance//(variance+(20-mean)**2) == 1087
    perfect = 2**60*(2**61-1)
    n = perfect-1
    digits = []
    while n > 1:
        j = (2*n)%3
        digits.append(j)
        n = (n-2*j)//3
    for j in reversed(digits):
        n = 3*n+2*j
    assert n == perfect-1 and len(digits) == 76 and sum(digits) == 71
    assert 9*5**51 > perfect+1 >= 9*5**50
    print("E61 centered counts", perfect-(6*2**116+2**59+1),
          perfect+(6*2**116+2**59+1))
    print("ternary digits", ''.join(map(str,reversed(digits))))


if __name__ == "__main__":
    check()
