"""Group all controlled histories by slope, leaving the other flank free."""
from argparse import ArgumentParser
from fractions import Fraction
from math import isqrt

from cycle_inverse import CHARTS, TARGET, controlled, height


def ceil(a,b):
    return -(-a//b)


def slopes(depth):
    rates = [Fraction(3**m,2**(d+m)) for d in range(depth+1) for m in range(d+1)]
    radii = {}
    for left in rates:
        for right in rates:
            if left < right:
                continue
            slope = left/right
            displacement = 2/left-1-1/slope
            radius = ceil(displacement.numerator,displacement.denominator)
            radii[slope] = max(radii.get(slope,0),radius)
    return radii


def intervals(chart, target, radii, transposed):
    _,q,c,t,_,_,k=chart
    if transposed:
        t,c = c,t
    ranges = []
    for slope, radius in radii.items():
        a,b = slope.as_integer_ratio()
        A,B = q*a,t*b+c*a
        anchor = (isqrt(B*B+4*A*b*(target-k))-B)//(2*A)
        ranges.append((max(1,anchor-radius),anchor+radius))
    merged = []
    for lo,hi in sorted(ranges):
        if hi < lo:
            continue
        if merged and lo <= merged[-1][1]+1:
            merged[-1] = (merged[-1][0],max(hi,merged[-1][1]))
        else:
            merged.append((lo,hi))
    return merged


def search(chart, target, radii, transposed):
    h,q,c,t,e,r,k=chart
    if transposed:
        t,c,e,r=c,t,r,e
    ranges=intervals(chart,target,radii,transposed)
    row_count=sum(b-a+1 for a,b in ranges)
    band=congruent=0
    hits=set()
    for lo,hi in ranges:
        # Restriction on the physical left address follows from H mod 5.
        # For target changes recompute its permitted residues modulo 10.
        left_residues = {j for j in range(10)
                         if (chart[-1]+chart[3]*j+chart[5]*(j%2)-target)%5==0}
        if transposed:
            progressions=(range(lo,hi+1),)
        else:
            progressions=(range(lo+(j-lo)%10,hi+1,10) for j in left_residues)
        for progression in progressions:
            for x in progression:
                rem=target-k-t*x-r*(x%2)
                denom=c+q*x
                ymin=max(1,ceil(rem-e-h*x,denom+h))
                ymax=rem//denom
                for y in range(ymin,ymax+1):
                    if transposed and y%10 not in left_residues:
                        continue
                    gap=rem-denom*y-e*(y%2)
                    if not 0<=gap<=h*(x+y):
                        continue
                    band+=1
                    if gap%(2*h):
                        continue
                    congruent+=1
                    pair=(y,x) if transposed else (x,y)
                    if height(chart,*pair)==target:
                        hits.add(pair)
    return row_count, len(ranges),band,congruent,sorted(hits)


def controls():
    radii=slopes(3)
    for chart in CHARTS.values():
        for x,y in ((7,10),(10,7),(1,2),(17,13),(53,71)):
            assert {controlled(controlled(controlled(x))),controlled(controlled(x)),controlled(x),x} & {
                controlled(controlled(controlled(y))),controlled(controlled(y)),controlled(y),y}
            n=height(chart,x,y)
            found=set()
            for sw in (False,True):
                found.update(search(chart,n,radii,sw)[-1])
            assert (x,y) in found,(x,y)
    print("slope strip planted controls passed")


if __name__=="__main__":
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("--depth",type=int,default=18)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    if args.check:
        controls()
    else:
        radii=slopes(args.depth)
        print("depth",args.depth,"distinct oriented slopes >= 1",len(radii),flush=True)
        for n, chart in CHARTS.items():
            for sw in (False,True):
                print("spine",n,"transposed",sw,*search(chart,TARGET,radii,sw),flush=True)
