"""Select nonaxial, nonconstant-response addresses in equal-deletion-height charts."""
import argparse
from collections import Counter
from functools import lru_cache
from math import gcd, isqrt


def read(x):
    total = 0
    while x:
        v, bit = divmod(x, 2)
        if bit == 0:
            total += v + v % 2
        x = v
    return total


def height(z, x, y):
    h, q, c, t, e, r, k = z
    return k + e*(y%2) + r*(x%2) + q*x*y + c*y + t*x + h*(read(x)+read(y))


@lru_cache(maxsize=4096)
def signed_divisors(n):
    small = [i for i in range(1, isqrt(n)+1) if n%i == 0]
    return [a for i in small for d in ({i, n//i}) for a in (d, -d)]


def ceil_div(a,b):
    return -(-a//b)


def branches(d, px, py, gap):
    dh, a, b, c, ey, ex, k = d
    assert dh == 0
    rhs = gap - k - ey*py - ex*px
    if a:
        product = a*rhs+b*c
        if product == 0:
            return
        for factor in signed_divisors(abs(product)):
            x0, y0 = factor-b, product//factor-c
            if x0%a or y0%a:
                continue
            x, y = x0//a, y0//a
            if x>0 and y>0 and x%2==px and y%2==py:
                yield ('point', x,y,0,0)
        return
    if not b or not c:
        return
    # c*x+b*y=rhs. Start at the smallest admissible x.
    g = gcd(b,c)
    if rhs%g:
        return
    if (rhs-c*px-b*py)%(2*g):
        return
    B, C, R = b//g, c//g, (rhs-c*px-b*py)//(2*g)
    X = (R*pow(C, -1, abs(B)))%abs(B)
    x = 2*X+px
    y = (rhs-c*x)//b
    dx = 2*abs(B)
    dy = -(c*dx)//b
    lower = ceil_div(1-x,dx)
    if dy>0:
        lower = max(lower,ceil_div(1-y,dy))
        yield ('ray', x+dx*lower,y+dy*lower,dx,dy)
    else:
        upper = (y-1)//(-dy)
        for s in range(lower,upper+1):
            yield ('point',x+dx*s,y+dy*s,0,0)


def ray_candidate(z, branch, target):
    _, x, y, dx, dy = branch
    h,q,c,t,e,r,k = z
    A=q*dx*dy
    B=q*(x*dy+y*dx)+c*dy+t*dx
    C=k+e*(y%2)+r*(x%2)+q*x*y+c*y+t*x
    if target < C:
        return None
    s=(isqrt(B*B+4*A*(target-C))-B)//(2*A)
    return x+s*dx, y+s*dy


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bank', help='nine_switch.cpp output')
    parser.add_argument('centers', nargs='*', type=int)
    args=parser.parse_args()
    centers=args.centers or [33506, *(b**3+744 for b in (32,96,960,5280,640320)), 2**60*(2**61-1)]
    counts=Counter()
    hits=Counter()
    nearest={}
    largest_point=0
    with open(args.bank) as bank:
        for line in bank:
            values=list(map(int,line.split()))
            label,z,d=values[:4],values[4:11],values[11:]
            for px in (0,1):
                for py in (0,1):
                    for gap in (-2,2):
                        for branch in branches(d,px,py,gap):
                            kind=branch[0]
                            counts[kind]+=1
                            counts[(kind,'spoke' if label[1] in label[2:] else 'internal')]+=1
                            if kind=='point':
                                largest_point=max(largest_point,height(z,*branch[1:3])+gap//2)
                            for center in centers:
                                pos=branch[1:3] if kind=='point' else ray_candidate(z,branch,center-gap//2)
                                if pos is None:
                                    pos=branch[1:3]
                                x,y=pos
                                actual=height(z,x,y)+gap//2
                                tested=[(actual,x,y)]
                                if kind=='ray':
                                    dx,dy=branch[3:]
                                    u,v=(x+dx,y+dy) if actual<center else (x-dx,y-dy)
                                    if u>0 and v>0: tested.append((height(z,u,v)+gap//2,u,v))
                                for near,u,v in tested:
                                    if center not in nearest or abs(near-center)<abs(nearest[center][0]):
                                        nearest[center]=(near-center,tuple(label),u,v,gap)
                                if actual==center:
                                    assert height(d,x,y)==gap
                                    if hits[center]==0:
                                        print('hit',center,label,x,y,gap,z,d,flush=True)
                                    hits[center]+=1
    print('branches',dict(counts),'largest point midpoint',largest_point)
    for c in centers: print('nearest',c,nearest.get(c))
    print('oriented hits',[(c,hits[c]) for c in centers])


if __name__=='__main__':
    main()
