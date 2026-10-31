#!/usr/bin/env python3
"""Exact native lifts along a strictly growing lcm flag; hostile controls."""
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd, lcm
import json


def norm(x):
    y = x % 1
    return min(y, 1-y)


def lift(d, speeds, s):
    j, period = 0, 1
    stages = []
    for v in speeds:
        h = d // gcd(d, v)
        H = h // gcd(h, period)
        if H == 1:
            return None
        alpha = v * (s+j) / d
        a = v * period * H // d
        # Nearest grid point to one half, with exact rational tie handling.
        x = H * (F(1, 2)-alpha)
        m = (x + F(1, 2)).numerator // (x + F(1, 2)).denominator
        t = (m * pow(a, -1, H)) % H
        j += period*t
        period *= H
        stages.append((h, H, period))
    return j, stages


def check(d, speeds, s):
    result = lift(d, speeds, s)
    if result is None:
        return False
    j, stages = result
    if not 0 <= j < d or any(norm(v*(s+j)/d) < F(1,4) for v in speeds):
        raise ValueError((d, speeds, s, result))
    if stages[-1][2] != lcm(*(d//gcd(d,v) for v in speeds)):
        raise ValueError("wrong retained period")
    return True


def main():
    pairs, triples = 0, 0
    phases = [F(0),F(1,10),F(1,4),F(1,2),F(7,13),F(1)]
    for d in range(2,65):
        for u,v in combinations(range(1,2*d+1),2):
            hu,hv=d//gcd(d,u),d//gcd(d,v)
            if min(hu,hv)<=1 or hu==hv:
                continue
            # Select the direction in which the second order adds a factor.
            row=(u,v) if hv//gcd(hu,hv)>1 else (v,u)
            for s in phases:
                if not check(d,row,s):
                    raise ValueError("unequal orders lost a lift")
                pairs+=1
    for d in range(2,13):
        for row in product(range(1,d+1),repeat=3):
            for s in phases:
                if check(d,row,s):
                    triples+=1
    # Equal order is a genuine pointwise obstruction, including at1/14.
    if any(min(norm(v*(F(7,50)+j)/2) for v in (1,7)) >= F(1,14) for j in range(2)):
        raise ValueError("lost equal-order hostile")
    # No lcm flag exists for orders2,3,6. This fibre really fails at1/4.
    if any(min(norm(v*(F(1,10)+j)/6) for v in (3,2,1)) >= F(1,4) for j in range(6)):
        raise ValueError("lost redundant-lcm hostile")
    d=2**150*3**100*5**80
    orders=(2**50,2**100,3**100,d)
    row=tuple(d//h for h in orders)
    s=F(17,31)
    if not check(d,row,s):
        raise ValueError("large native lift failed")
    j,stages=lift(d,row,s)
    print(json.dumps({"unequal_order_pair_checks":pairs,"strict_flag_triple_checks":triples,
        "radius":"1/4","phases":[str(x) for x in phases],"pair_universe":"2<=d<=64,1<=u<v<=2d, distinct nontrivial orders",
        "triple_universe":"2<=d<=12, ordered rows in[1,d]^3 passing the strict flag condition",
        "equal_order_hostile":[2,1,7,"7/50","1/14"],
        "redundant_lcm_hostile":{"d":6,"speeds":[3,2,1],"s":"1/10","radius":"1/4"},
        "large_d":str(d),"large_j":str(j),"large_phases":[str(norm(v*(s+j)/d)) for v in row],"pass":True},indent=2))


if __name__=="__main__":
    main()
