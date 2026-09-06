"""Independent no-import replay of the sharp LRC full-cap finite base.

Solve for the middle carrier coordinate, and compare integer-scaled roofs.
Uses only the standard library and does not import the primary companion.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import gcd


def check(ok, message):
    if not ok:
        raise RuntimeError(message)


def carriers(a,b,c):
    lx,ly,lz=[(3*t-1)//14 for t in (b+c,a+c,a+b)]
    points=[]
    for x in range(-lx,lx+1):
        if x%3 == 0:
            continue
        for z in range(-lz,lz+1):
            if z%3 == 0:
                continue
            v=-a*x-c*z
            if v%b:
                continue
            y=v//b
            if abs(y)<=ly and y%3:
                points.append((x,y,z))
    return points


def cap(points):
    # Projection onto coordinates 1 and 3 is injective on this kernel plane.
    for p,q,r in combinations(points,3):
        if (q[0]-p[0])*(r[2]-p[2]) == (r[0]-p[0])*(q[2]-p[2]):
            return False
    return True


def main():
    count=Counter()
    sharp=Fraction(0)
    equal=[]
    choices=[v for v in range(1,101) if v%2 and v%3]
    for a,b,c in combinations(choices,3):
        if gcd(gcd(a,b),c)!=1:
            continue
        count['universe']+=1
        points=carriers(a,b,c)
        if not cap(points):
            continue
        count['cap']+=1
        count['cap_size_'+str(len(points))]+=1
        check(len(points)<=8,'midpoint cap bound')
        if len(points)<3:
            continue
        count['multidirection_cap']+=1
        # Common denominator 14abc clears both clipped roofs exactly.
        totals=[]
        for i,w in enumerate((a,b,c)):
            totals.append(sum(min(6*a*b,w*(3*(a+b+c-w)-14*abs(p[i]))) for p in points))
        value=Fraction(min(totals),14*a*b*c)
        if value>sharp:
            sharp,equal=value,[(a,b,c)]
        elif value==sharp:
            equal.append((a,b,c))
    check(count['universe']==5409,'complete primitive finite base')
    check(count['cap']==1723 and count['multidirection_cap']==1230,'cap classes')
    check(sharp==Fraction(204,5957) and equal==[(23,29,37)],'unique exact extremum')
    check(Fraction(24,707)<sharp<Fraction(6,77),'strict all-height tail')
    print('Independent middle-coordinate carrier enumeration; integer-scaled projection roofs')
    print(dict(sorted(count.items())))
    print('sharp=',sharp,'equality=',equal)
    print('tail_at_c101=',Fraction(24,707),'target_gap=',Fraction(6,77)-sharp)
    print('PASS; all-height consequence additionally uses the analytic full-dictionary midpoint lemma.')


if __name__=='__main__':
    main()
