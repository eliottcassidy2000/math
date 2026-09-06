"""Independent audit: direct two-half-lift Boolean wall decomposition."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import gcd, isqrt
import sys

sys.stdout.reconfigure(newline="\n")
GATES = 0

def need(ok, text):
    global GATES
    GATES += 1
    if not ok:
        raise ArithmeticError(text)


def inert_sums(limit):
    prime = [True] * (limit + 1)
    prime[0] = prime[1] = False
    for i in range(2, isqrt(limit) + 1):
        if prime[i]:
            for j in range(i*i, limit+1, i):
                prime[j] = False
    sums = {1}
    for p in range(2, limit+1):
        if prime[p] and p % 3 == 2:
            sums = {n * p**e for n in sums for e in range(3) if n * p**e <= limit}
    return sums


SUMS = inert_sums(600)

def admissible(p, q):
    return p < q and p % 2 == q % 2 == 1 and p % 3 and q % 3 and gcd(p, q) == 1 and p + q in SUMS


def distance(x):
    y = x % 1
    return min(y, 1-y)


def carrier(p, q, y):
    return all(any(distance(s * (y+e) / 2) < F(1,14) for s in (p,q)) for e in (0,1))


@lru_cache(None)
def geometry(p, q):
    # Independent walls are pulled from each physical half-lift, not from
    # intersections of parity-labelled teeth.
    walls = {F(0), F(1)}
    for s in (p, q):
        for k in range(s+1):
            for e in (0,1):
                for sign in (-1,1):
                    y = (2*F(k) + F(sign,7)) / s - e
                    if 0 <= y <= 1:
                        walls.add(y)
    walls = sorted(walls)
    intervals = []
    mass = F(0)
    for left, right in zip(walls, walls[1:]):
        if carrier(p, q, (left+right)/2):
            mass += right-left
            if intervals and intervals[-1][1] == left and carrier(p, q, left):
                intervals[-1] = (intervals[-1][0], right)
            else:
                intervals.append((left,right))
    need(not carrier(p,q,F(0)) and not carrier(p,q,F(1)), "no circular joining at chart cut")
    for a,b in intervals:
        need(not carrier(p,q,a) and not carrier(p,q,b), "literal strict excluded component endpoints")
    width = max((b-a for a,b in intervals), default=F(0))
    need(mass == sum((b-a for a,b in intervals),F(0)), "component measure agrees with cell measure")
    return mass, width, tuple(intervals)


def bernoulli_mass(p,q):
    def B(x):
        u=x%1
        return u*u-u+F(1,6)
    return F(2,49) + 2*(B(F(1,2)+F(q-p,14))-B(F(1,2)+F(q+p,14)))/(p*q)


def main():
    cap, widthcap = F(20,469), F(1,49)
    need(cap-F(2,49)==F(6,3283) and F(1,548)<F(6,3283), "strict infinite mass tail")
    head = [(p,q) for p in range(1,274) for q in range(p+1,274) if p*q<=273 and admissible(p,q)]
    need(len(head)==21, "independent full product head")
    masses = {pair:geometry(*pair)[0] for pair in head}
    need(max(masses.values())==cap and [k for k,v in masses.items() if v==cap]==[(1,67)], "unique sharp mass")
    lowq = [(p,q) for p in range(1,15) for q in range(p+1,15) if admissible(p,q)]
    need(lowq==[(7,13)] and F(2,105)<widthcap, "strict infinite width tail and complete head")
    need(geometry(7,13)==(F(2,49),widthcap,((F(1,7),F(8,49)),(F(41,49),F(6,7)))), "exact sharp-width components")
    head67 = [(p,q) for p in range(1,68) for q in range(p+1,68) if admissible(p,q)]
    need(len(head67)==46, "independent complete Pareto head")
    profiles = {pair:geometry(*pair)[:2] for pair in head67}
    actual = {(p,q,*value) for (p,q),value in profiles.items() if not any(other!=value and other[0]>=value[0] and other[1]>=value[1] for other in profiles.values())}
    expected = {(7,13,F(2,49),F(1,49)),(19,25,F(138,3325),F(37,3325)),(5,41,F(12,287),F(2,287)),(5,53,F(78,1855),F(2,371)),(1,67,F(20,469),F(2,469))}
    need(actual==expected, "independently reconstructed complete five-profile Pareto frontier")
    for value in profiles.values():
        need(any(value[0]<=m and value[1]<=b for _,_,m,b in expected), "every head profile dominated")
    controls=0
    for q in range(3,76,2):
        for p in range(1,q,2):
            if gcd(p,q)!=1:
                continue
            m,b,_=geometry(p,q)
            need(m==bernoulli_mass(p,q), "physical Boolean walls versus analytic Bernoulli formula")
            need(b<=F(2,7*q), "width bound on unrestricted odd control bank")
            controls+=1
    atlas=[(p,q) for p in range(1,356) for q in range(p+1,356) if p+q<=356 and admissible(p,q)]
    need(len(atlas)==548, "independent inherited-atlas subclass")
    need(all((p,q) in atlas for p,q,_,_ in expected), "all five extremals belong to actual ratio atlas")
    for p,q in atlas:
        m=bernoulli_mass(p,q)
        need(m<=cap, "full atlas mass consumer")
    for pair,value in (((1,9),F(4,63)),((5,11),F(18,385)),((1,11),F(4,77))):
        need(geometry(*pair)[0]==value>cap, "arithmetic-hypothesis hostile")
    need(geometry(1,3)[:2]==(F(0),F(0)), "empty-carrier control")
    M, scaled_L=F(12,287),F(2,371)
    need(all(M>=m or scaled_L>=b for _,_,m,b in expected) and M<cap and scaled_L<widthcap, "genuine gain for joint sufficient sidecar inequality")
    need(2*cap==F(40,469)<F(8,91), "inclusive q4 absorption arithmetic")
    need(cap+F(8,63)==F(716,4221)<F(20,117), "inclusive q2 absorption arithmetic")
    need((F(8,91)-F(40,469))/F(8,91)==F(2,67), "relative q4 improvement")
    print("AUDIT PASS: all-height five-profile inert odd-3-unit envelope")
    print("GEOMETRY: literal half-lift Boolean evaluation on rational wall cells; excluded walls kept")
    print("FINITE UNIVERSES: product head",len(head),"; q<=67 head",len(head67),"; unrestricted primitive odd controls",controls,"; inert atlas subclass",len(atlas))
    for p,q,m,b in sorted(actual):
        print("PROFILE:",p,q,"mass",m,"component",b)
    print("CONSUMERS: inclusive mass/component disjunction and inherited absorption arithmetic PASS")
    print("SCOPE: no universal body lower bound or LRC14 closure")
    print("ACTIVE GATES:",GATES)


if __name__=="__main__":
    main()
