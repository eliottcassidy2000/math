"""Invert two unrestricted address words on affine lines with a sparse determinant."""
from argparse import ArgumentParser
from itertools import combinations
from math import isqrt

from fermat_receiver import TARGET, address_height
from hyperbola_controls import chart


def offsets(limit):
    powers = [1 << i for i in range(limit.bit_length()+1)]
    yield 0
    for lo, hi in combinations(powers, 2):
        if hi-lo <= limit:
            yield hi-lo
    for lo, hi, lo2, hi2 in combinations(powers, 4):
        value = hi2-lo2+hi-lo
        if value <= limit:
            yield value


def lines(data, target, a, b, negative):
    """All nonnegative solutions of b*x-a*y=+/-v, v having <=2 binary one runs."""
    h,q,c,t,e,r,k = data
    if negative:
        a,b = b,a
        c,t,e,r = t,c,r,e
    inv = pow(a, -1, b)
    limit = b*((target-k)//t)
    quadratic = q*a*b
    linear = t*a+c*b
    total = envelope = local = 0
    hits = []
    for v in offsets(limit):
        if negative and v == 0:
            continue
        total += 1
        y0 = -v*inv % b
        x0 = (v+a*y0)//b
        rem = target-(k+q*x0*y0+t*x0+c*y0)
        if rem < 0:
            continue
        slope = linear+q*(a*y0+b*x0)
        step = (isqrt(slope*slope+4*quadratic*rem)-slope)//(2*quadratic)
        x,y = x0+a*step, y0+b*step
        residual = rem-step*(slope+quadratic*step)-r*(x&1)-e*(y&1)
        if 0 <= residual <= h*(x+y):
            envelope += 1
            if residual % (2*h) == 0:
                local += 1
                physical = (y,x) if negative else (x,y)
                if address_height(data,*physical) == target:
                    hits.append(physical)
    return total,envelope,local,hits


def directions(n):
    h = 2**n+1
    r = (3*h-1)//5
    return [(1,1)] + [(a,b) for a in (3,4,5) for b in (3,4,5) if a!=b] + [(h,r),(r,h)]


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("n", type=int, choices=(4,8))
    parser.add_argument("--direction", type=int, choices=range(9))
    args = parser.parse_args()
    data = chart(args.n)
    for index,(a,b) in enumerate(directions(args.n)):
        if args.direction is not None and index!=args.direction:
            continue
        # A multiple of the direction plus the line b*x-a*y=+/-1.
        for negative in (False,True):
            y0 = (-pow(a,-1,b))%b
            x0 = (1+a*y0)//b
            x,y = x0+2*a, y0+2*b
            if negative:
                # Swap an actual solution from the transposed direction.
                x0 = (-pow(b,-1,a))%a
                x,y = x0+2*a, (1+b*x0)//a+2*b
            planted = address_height(data,x,y)
            assert (x,y) in lines(data,planted,a,b,negative)[3]
        for negative in (False,True):
            print(args.n, a, b, '-' if negative else '+',
                  *lines(data,TARGET,a,b,negative), flush=True)
