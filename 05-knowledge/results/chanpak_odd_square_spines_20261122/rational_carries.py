"""Invert two binary suffixes by clearing a common rational phase denominator."""
from argparse import ArgumentParser
from collections import Counter
from itertools import product
from math import gcd, isqrt, lcm
from pathlib import Path
from random import Random
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "chanpak_decimal_kernel_20261122"))
from fermat_receiver import TARGET, address_height
from hyperbola_controls import chart


def reflected(word, width):
    return int(f"{word:0{width}b}"[::-1], 2)


def joint_approximants(width, i, j):
    """Use an exact LLL basis in rank three to choose one shared denominator."""
    from fractions import Fraction
    from sympy import Matrix
    base = 1 << width
    if not (0 <= i < base and 0 <= j < base):
        raise ValueError("suffix exceeds the selected width")
    z, w = (base-1-reflected(a, width) for a in (i, j))
    if width < 2:
        return (i,z,base,z), (j,w,base,w)
    lattice = Matrix([[1,z,w], [0,base,0], [0,0,base]])
    reduced = lattice.lll(delta=Fraction(3,4))
    d, ei, ej = map(int, reduced.row(0))
    if d < 0:
        d, ei, ej = -d, -ei, -ej
    # Test the LLL norm bound and the lattice certificate in integers.
    assert (d*d+ei*ei+ej*ej)**3 <= 8*base**4 and d > 0
    u, rem_i = divmod(d*z-ei, base)
    v, rem_j = divmod(d*w-ej, base)
    assert rem_i == rem_j == 0 and 0 <= u <= d and 0 <= v <= d
    return (i,u,d,z), (j,v,d,w)


def phase_bank(width, limit, radius):
    """Deduplicate suffixes. Each phase cell carries a usable rational approximant."""
    base = 1 << width
    bank = {}
    for denominator in range(1, limit + 1):
        for numerator in range(denominator + 1):
            if gcd(numerator, denominator) != 1:
                continue
            center = (2 * base * numerator + denominator) // (2 * denominator)
            for z in range(max(0, center-radius), min(base-1, center+radius)+1):
                word = reflected(base-1-z, width)
                score = abs(denominator*z-base*numerator) + denominator
                if word not in bank or score < bank[word][0]:
                    bank[word] = (score, numerator, denominator, z)
    return [(word, *data[1:]) for word, data in sorted(bank.items())]


def one_side(data, target, width, left, right, stats):
    h, q, c, t, e, r, k = data
    base = 1 << width
    i, u, d1, z = left
    j, v, d2, w = right
    if target < k:
        return
    den = lcm(d1, d2)
    u *= den // d1
    v *= den // d2
    ci = den*z-base*u
    cj = den*w-base*v
    constant = address_height(data, i, j)
    if target < constant:
        return
    delta = target-constant
    # Cases X=0 and Y=0 have a direct one-dimensional interval inverse.
    for swapped in (False, True):
        low_slope = t+q*j if swapped else c+q*i
        low_phase = z if swapped else w
        slope = base*low_slope + h*low_phase
        floor_high = max(int(swapped), (delta+slope+h-1)//(slope+h))
        ceil_high = delta//slope
        g1 = gcd(slope,h)
        if delta%g1:
            continue
        stride = h//g1
        residue = (delta//g1)*pow(slope//g1,-1,stride)%stride
        for high in range(floor_high+(residue-floor_high)%stride,ceil_high+1,stride):
            x,y = (base*high+i,j) if swapped else (i,base*high+j)
            stats["face_reads"] += 1
            if address_height(data,x,y) == target:
                yield x,y
    minimum_x, minimum_y = base+i, base+j
    cap = (target-k)//q
    if minimum_x*minimum_y > cap:
        return
    sum_bound = max(minimum_x + cap//minimum_x, minimum_y + cap//minimum_y)
    high_total = (sum_bound-i-j)//base
    # Exact deletion of the low words leaves a carry depending on the two high words.
    lo = high_total*min(0, ci, cj)
    hi = high_total*max(0, ci+den, cj+den)
    P = den*q*base
    A = den*(t+q*j)+h*u
    C = den*(c+q*i)+h*v
    remainder = den*(target-constant)
    common = base*gcd(P, gcd(A, C))
    g = gcd(common, h)
    if remainder%g:
        return
    step = common//g
    residue = (remainder//g)*pow(h//g, -1, step)%step
    short_high = isqrt(cap)//base
    for carry in range(lo+(residue-lo)%step, hi+1, step):
        stats["carries"] += 1
        value = (remainder-h*carry)//base
        if value < 0:
            continue
        # The smaller high word is bounded by the geometric mean, including zero.
        for swapped in (False, True):
            slope, other = (C, A) if swapped else (A, C)
            g2 = gcd(other, P)
            g3 = gcd(slope, g2)
            if value%g3:
                continue
            stride = g2//g3
            start = (value//g3)*pow(slope//g3,-1,stride)%stride
            first = 1+(start-1)%stride
            for small in range(first, min(short_high, value//slope)+1, stride):
                big, rem = divmod(value-slope*small, P*small+other)
                if rem or big < small + int(swapped):
                    continue
                X, Y = (big, small) if swapped else (small, big)
                x, y = base*X+i, base*Y+j
                stats["divisor_pairs"] += 1
                if address_height(data,x,y) == target:
                    yield x,y


def census(data, target, width, limit, radius):
    bank = phase_bank(width, limit, radius)
    stats = Counter(words=len(bank), pairs=len(bank)**2, carries=0, divisor_pairs=0, face_reads=0)
    hits = []
    for left,right in product(bank, repeat=2):
        hits.extend(one_side(data, target, width, left, right, stats))
    return hits, stats


def checks():
    rng = Random(2311)
    for n in (4, 8):
        data = chart(n)
        for case in range(30):
            width = 1+case%4
            base = 1<<width
            x,y = rng.randrange(1,80),rng.randrange(1,75)
            i,j = x%base,y%base
            # Deliberately nonoptimal approximants and distinct denominators.
            def entry(word, denominator):
                return word, rng.randrange(denominator+1), denominator, base-1-reflected(word,width)
            left,right = entry(i,11),entry(j,3)
            target = address_height(data,x,y)+2*(case%2)
            expected = {(a,b) for a in range(i,(target-data[-1])//data[3]+1,base)
                        for b in range(j, (target-data[-1]-data[3]*a)//(data[2]+data[1]*a)+1,base)
                        if address_height(data,a,b)==target}
            actual = set(one_side(data,target,width,left,right,Counter()))
            assert actual == expected, (case, actual, expected)
    print("60 complete small fibers passed")


if __name__ == "__main__":
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    parser.add_argument("--width",type=int,default=42)
    parser.add_argument("--limit",type=int,default=11)
    parser.add_argument("--radius",type=int,default=0)
    parser.add_argument("--spine",type=int,choices=(4,8))
    parser.add_argument("--suffix",nargs=2,type=int,metavar=("I","J"),
                        help="choose both rational phases from arbitrary suffixes with SymPy LLL")
    args=parser.parse_args()
    if args.check:
        checks()
    else:
        if args.width < 1 or args.limit < 1 or args.radius < 0:
            parser.error("width and limit must be positive and radius nonnegative")
        for n in ((args.spine,) if args.spine else (4,8)):
            if args.suffix is None:
                hits, stats = census(chart(n),TARGET,args.width,args.limit,args.radius)
                print(n,args.width,args.limit,args.radius,dict(stats),"hits",hits,flush=True)
            else:
                left, right = joint_approximants(args.width, *args.suffix)
                stats = Counter()
                hits = list(one_side(chart(n),TARGET,args.width,left,right,stats))
                print(n,"numerators",left[1],right[1],"denominator",left[2],
                      dict(stats),"hits",hits,flush=True)
