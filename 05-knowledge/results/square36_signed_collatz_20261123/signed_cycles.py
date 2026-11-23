"""Integral bounds and certificates for two choices of signed odd acceleration."""
import argparse
from math import gcd


def v2(n):
    return (n & -n).bit_length()-1


def accelerated(n, sign):
    raw = 3*n + sign
    return raw >> v2(raw)


def child(y, a):
    assert y > 0 and y % 2 and y % 3 and a >= 2
    sign_y = 1 if y % 3 == 1 else -1
    return ((y << a) - (-1)**a * sign_y)//3


def sterile_residue(y, k):
    """Digit lift for exponents whose child is a multiple of 3**k."""
    assert k >= 0 and y % 3
    sign_y = 1 if y % 3 == 1 else -1
    a, period = 0, 1
    for level in range(1, k+1):
        modulus = 3**(level+1)
        hits = [a+j*period for j in range(3)
                if (pow(-2, a+j*period, modulus)*y-sign_y) % modulus == 0]
        assert len(hits) == 1
        a = hits[0]
        period *= 3
    return a


def bracket_excludes(m, period):
    """True if the cycle minimum cannot reach m, by integer powers alone."""
    lower = (9*m-3)**period
    upper = (9*m+3)**period
    scale = m**period
    a = (3**period).bit_length()-1
    return lower > (scale << (2*a)) and (scale << (2*a+2)) > upper


def prove_period_bound(max_period, cutoff):
    for period in range(2, max_period+1, 2):
        assert bracket_excludes(cutoff, period), (cutoff, period)


def census(max_period, cutoff):
    cycles = set()
    for m in range(1, cutoff, 2):
        if m%3 == 0:
            continue
        for initial_sign in (-1, 1):
            n, sign = m, initial_sign
            for step in range(1, max_period+1):
                n, sign = accelerated(n, sign), -sign
                if n < m:
                    break
                if (n,sign) == (m,initial_sign):
                    start = (n, sign)
                    orbit = [start]
                    for _ in range(step-1):
                        n, sign = accelerated(n,sign), -sign
                        orbit.append((n,sign))
                    i = orbit.index(min(orbit))
                    cycles.add(tuple(orbit[i:] + orbit[:i]))
                    break
    return sorted(cycles, key=len)


def residue_controls():
    # y mod 3*36*37 fixes the child residues at both moduli.
    for y in range(1, 3*36*37, 2):
        if gcd(y,3*37)>1:
            continue
        eps=1 if y%3==1 else -1
        children=[child(y,a) for a in range(2,38)]
        assert {x%36 for x in children[:18]}==set(range(1,36,2))
        for a in range(2,20):
            x,z=children[a-2],children[a+16]
            assert x%36==z%36 and x%37!=z%37
            assert (3*(x+z)+2*(-1)**a*eps)%37==0
    assert pow(-2,18,37)==36 and pow(-2,12,37)!=1
    assert child(1,2)%37==1 and child(1,20)%37==23


def controls():
    residue_controls()
    minus_cycles = ((1,), (5,7), (17,25,37,55,41,61,91))
    for orbit in minus_cycles:
        assert all(accelerated(n,-1)==orbit[(i+1)%len(orbit)] for i,n in enumerate(orbit))
    for y in (1,5,7,11,17, 41,61):
        residue = sterile_residue(y, 5)
        for a in range(2, 2+3**5):
            x = child(y, a)
            sign = 1 if x%4 == 1 else -1
            assert v2(3*x+sign)==a and accelerated(x,sign)==y
            assert (x%3**5==0) == ((a-residue)%3**5==0)
    for m in range(1, 1024, 2):
        for sign in (-1,1):
            n,s = m,sign
            for _ in range(v2(5*m+sign)-1):
                assert v2(3*n+s)==1
                n,s=accelerated(n,s),-s
            assert v2(3*n+s)>=2


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extended', action='store_true', help='classify periods through 4096')
    args=parser.parse_args()
    controls()
    period, minimum_cutoff = (4096,2538899) if args.extended else (1000,56135)
    prove_period_bound(period,minimum_cutoff)
    cycles=census(period,minimum_cutoff)
    expected=(13,19,29,43,65,97,73,109,41,61,23,17)
    assert len(cycles)==2 and [x for x,s in cycles[0]]==[1,1]
    assert tuple(x for x,s in cycles[1])==expected and cycles[1][0][1]==-1
    print('all even periods through',period,'have minimum below',minimum_cutoff)
    for orb in cycles:
        print('period', len(orb), 'cycle', orb, 'valuation word',
              [v2(3*x+s) for x,s in orb])
