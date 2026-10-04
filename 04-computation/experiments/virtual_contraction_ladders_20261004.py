"""Exact sibling contractions before coefficient stopping, with 19-adic lifts.

Run from the repository root. No third-party packages or orbit-search oracle.
"""
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache

from inverse_ray_ternary_addresses_20261004 import (
    ROOT, bit_bounds, expand, extend, kappa, mod3, ranks,
)
from paley_geometry_observations_20261004 import residue as certificate_residue

CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def v2(n):
    if n <= 0:
        raise ValueError("positive integer required")
    return (n & -n).bit_length()-1


def odd_step(n):
    a = v2(3*n+1)
    return (3*n+1) >> a, a


@dataclass(frozen=True)
class Chart:
    s: int
    j: int
    P: int
    Q: int
    R: int
    C: int
    r: int
    modulus: int


@lru_cache(None)
def chart(s):
    if type(s) is not int or s < 1:
        raise ValueError("s must be a positive integer")
    j = 1
    while 3**(j+1) < 2**(j+1+2*s):
        j += 1
    P, Q = 3**j, 2**j
    R = Q*4**s
    C = P-Q*(1+(4**s-1)//3)
    modulus = 4*R
    r = ((3*R-C)*pow(P,-1,modulus)) % modulus
    need(P < R and 3*P > 2*R and 0 < C < P, "chart inequalities")
    need(v2(r+1) == j+1, "source leading-one guard")
    return Chart(s,j,P,Q,R,C,r,modulus)


def recognize(n):
    """Return a verified smaller dependency, never a completed-home flag."""
    if type(n) is not int or n < 3 or n % 2 == 0:
        return None
    j = v2(n+1)-1
    if j < 3:
        return None
    P = 3**j
    s = (P.bit_length()-j+1)//2
    if s < 1:
        return None
    row = chart(s)
    if row.j != j or n % row.modulus != row.r:
        return None
    y = (row.P*n+row.C)//row.R
    if not (0 < y < n and y % 4 == 3):
        return None
    return row,y


def suffix_residue(b, modulus):
    """The known two-step source (32*64**b-5)/9, without expansion."""
    return ((32*pow(64,b,9*modulus)-5) % (9*modulus))//9


@lru_cache(None)
def suffix_address(s):
    row = chart(s)
    target = (row.P*row.r+row.C)//row.R
    b, modulus = 0, 1
    for _ in range(row.j):
        nxt = 3*modulus
        choices = [b+d*modulus for d in range(3)
                   if suffix_residue(b+d*modulus,nxt) == target % nxt]
        need(len(choices) == 1, "one ternary suffix address")
        b, modulus = choices[0], nxt
    return b


def completed_residue(s, h, modulus):
    row = chart(s)
    b = suffix_address(s)+row.P*h
    # Retain P extra precision before dividing; valid even if 3 divides modulus.
    numerator = (row.R*suffix_residue(b,row.P*modulus)-row.C) % (row.P*modulus)
    need(numerator % row.P == 0, "completed source integer guard")
    return numerator//row.P


def inverse_exponent(parent, a):
    residue = mod3(parent,2)
    source_mod3 = ((pow(2,a,9)*residue-1) % 9)//3
    first = kappa(residue,source_mod3)
    need(a >= first and (a-first)%6 == 0, "inverse exponent channel")
    return extend(parent,source_mod3,(a-first)//6)


def completed_certificate(s, h):
    if type(h) is not int or h < 0:
        raise ValueError("h must be a nonnegative integer")
    row = chart(s)
    b = suffix_address(s)+row.P*h
    node = inverse_exponent(ROOT,4+6*b)
    child = inverse_exponent(node,1)
    node = inverse_exponent(node,2*s+1)
    for _ in range(row.j):
        node = inverse_exponent(node,1)
    return node,child,b


def route_digit(s, target, depth):
    """Select h mod19**(depth-1) inside this chart's one first-digit fibre."""
    if type(depth) is not int or depth < 1:
        raise ValueError("positive depth required")
    target %= 19**depth
    if completed_residue(s,0,19) != target % 19:
        return None
    row = chart(s)
    gamma = (32*row.R*pow(64,suffix_address(s),19)*pow(9,-1,19)) % 19
    h, parameter_modulus, source_modulus = 0, 1, 19
    for _ in range(1,depth):
        nxt = 19*source_modulus
        error = (target-completed_residue(s,h,nxt)) % nxt
        need(error % source_modulus == 0, "retained source digits")
        digit = ((error//source_modulus)*pow(gamma,-1,19)) % 19
        h += digit*parameter_modulus
        parameter_modulus *= 19
        source_modulus = nxt
        need(completed_residue(s,h,source_modulus) == target % source_modulus,
             "selected next 19-adic digit")
    return h,parameter_modulus


def main():
    rows = [chart(s) for s in range(1,31)]
    literal = 0
    for row in rows:
        need(3*row.s <= row.j < 4*row.s, "linear clock bounds")
        for k in range(1,9):
            n = row.r+row.modulus*k
            found = recognize(n)
            need(found is not None, "source recognizer")
            got,y = found
            need(got == row and row.P*n+row.C == row.R*y, "dependency formula")
            x = n
            total = 0
            for i,a in enumerate([1]*row.j+[2*row.s+1],1):
                x,actual = odd_step(x)
                total += actual
                need(actual == a and x > n, "actual growing prefix")
                need(3**i > 2**total, "no coefficient crossing yet")
            need(x == odd_step(y)[0], "real common future")
            for t in range(3):
                m = 19**(t+1)
                for target in range(min(m,19)):
                    lift = ((target-row.r)*pow(row.modulus,-1,m)) % m
                    n2 = row.r+row.modulus*(lift+m)
                    need(n2 % m == target and recognize(n2) is not None,
                         "every residue has a smaller-dependency source")
        for h in range(3):
            cert,child,b = completed_certificate(row.s,h)
            need(ranks(cert)[0] == row.j+2 and ranks(child)[0] == 2, "first-hit odd ranks")
            need(ranks(cert)[1] == 2*row.j+2*row.s+6*b+7, "ordinary first-hit rank")
            for m in [3,19,361,64,3**5*19**2]:
                need(completed_residue(row.s,h,m) == certificate_residue(cert,m),
                     "independent AST source residue")
            if bit_bounds(cert)[1] <= 100000:
                n,y = expand(cert,100000),expand(child,100000)
                need(y == (32*64**b-5)//9, "literal suffix formula")
                need(n == (row.R*y-row.C)//row.P and recognize(n) == (row,y),
                     "literal completed family")
                x = n
                for i,a in enumerate([1]*row.j+[2*row.s+1,4+6*b]):
                    x,actual = odd_step(x)
                    need(actual == a, "literal route exact valuations")
                    need(x > n if i <= row.j else x == 1, "first return at final step")
                literal += 1
    for s in range(1,9):
        row = chart(s)
        need(pow(64,row.P,361) == (1+19*row.P)%361, "19-adic multiplier derivative")
        for depth in range(1,4):
            base = completed_residue(s,0,19)
            seen = {completed_residue(s,h,19**depth) for h in range(19**(depth-1))}
            need(seen == {base+19*d for d in range(19**(depth-1))}, "full first-digit fibre")
            for target in seen:
                answer = route_digit(s,target,depth)
                need(answer is not None and completed_residue(s,answer[0],19**depth)==target,
                     "recursive digit selector")
    huge_s,huge_depth = 20,24
    base = completed_residue(huge_s,0,19)
    target = base+19*(19**(huge_depth-1)-1)
    h,period = route_digit(huge_s,target,huge_depth)
    # h+period gives the same residue, and h>=1 places it in the proved tail.
    huge,_,b = completed_certificate(huge_s,h+period)
    need(completed_residue(huge_s,h+period,19**huge_depth)==target, "huge selected residue")
    need(certificate_residue(huge,19**huge_depth)==target, "huge AST residue independently")
    density = sum((Fraction(1,chart(s).modulus) for s in range(5,31)),Fraction())
    need(all(chart(s).j+1 > 16 for s in range(5,31)), "outside depth16 crossing bank")
    print("Virtual charts:s=1..30;240 literal growing prefixes; all firstj+1 coefficients expand.")
    print("First family:n=79+128k; smaller dependency67+108k; actual join after4 steps.")
    print("Completed families:90 symbolic certificates;",literal,"literal replays under100000bits.")
    print("First suffix addresses:",[(s,chart(s).j,suffix_address(s)) for s in range(1,7)])
    print("First19-adic fibres:",[(s,completed_residue(s,0,19)) for s in range(1,13)])
    print("Exact19-digit census:8 charts,depth1..3; one complete first-digit fibre per chart.")
    print("Huge selected certificate:odd rank",ranks(huge)[0],"source bit lower bound",bit_bounds(huge)[0])
    print("Smaller-dependency cylinder density,s5..30 (not root coverage):",density)
    print("Remaining density afters30 <= 1/(124*2^150), fromj>=3s.")
    print("Boundary: dependencies need certificates;19-adic digit support is not integer coverage.")
    print("PASS:",CHECKS,"explicit exact checks, active under-O.")


if __name__ == '__main__':
    main()
