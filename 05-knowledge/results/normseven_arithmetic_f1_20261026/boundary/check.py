"""Independent exact controls for the central trichotomy and exceptional conics."""
import sys as _output_sys
_output_sys.stdout.reconfigure(newline="\n")
from fractions import Fraction as F
from itertools import product
from math import gcd, isqrt
import json

gates = 0


def require(condition, message):
    global gates
    gates += 1
    if not condition:
        raise RuntimeError(message)


def square(x):
    x = F(x)
    return x >= 0 and isqrt(x.numerator)**2 == x.numerator and isqrt(x.denominator)**2 == x.denominator


def values(d, x, z):
    return list(map(F,[0, d[0]*x, -d[1]*x, d[2]*z, -d[3]*z, d[4]*(x+z), -d[5]*(x+z)]))


def direct_lift_count(v):
    nonzero = [t for t in v if t]
    return 2**(len(nonzero)-1) if all(square(F(t)/F(nonzero[0])) for t in nonzero) else 0


def table(d):
    q = [-F(d[1], d[0]), -F(d[3], d[2]), -F(d[5], d[4])]
    ok = [square(t) for t in q]
    triples = [ok[1] and ok[2] and square(F(d[4], d[2])),
               ok[0] and ok[2] and square(F(d[4], d[0])),
               ok[0] and ok[1] and square(-F(d[2], d[0]))]
    return ok, triples


def factor(n):
    out = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0)+1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0)+1
    return out


def norm_criterion(x):
    x = F(x)
    return x > 0 and all(p % 4 != 3 or e % 2 == 0 for p,e in factor(x.numerator*x.denominator).items())


def two_square(n):
    for a in range(isqrt(n)+1):
        b = isqrt(n-a*a)
        if a*a+b*b == n:
            return (a,b)
    return None


base = [(x,z) for x in range(-5,6) for z in range(-5,6)
        if gcd(x,z) == 1 and (x > 0 or (x == 0 and z > 0))]
special = [(0,1),(1,0),(1,-1)]
finite_histogram = {"zero":0, "eight":0, "all_pair_classes_square_not_decided_here":0}
point_checks = 0
for d in product((-2,-1,1,2), repeat=6):
    ok, triples = table(d)
    for y, passes in zip(special, triples):
        require(direct_lift_count(values(d,*y)) == 8*passes, ("triple",d,y))
    if all(ok):
        finite_histogram["all_pair_classes_square_not_decided_here"] += 1
        continue
    require(sum(triples) <= 1, ("two special fibres without all square classes",d))
    finite_histogram["eight" if any(triples) else "zero"] += 1
    observed = 0
    for y in base:
        count = direct_lift_count(values(d,*y))
        expected = 8 if y in special and triples[special.index(y)] else 0
        require(count == expected, ("finite central locus",d,y,count,expected))
        observed += count
        point_checks += 1
    require(observed == 8*sum(triples), ("finite total",d))

# Direct surface equations, not just squareclass tests.
node_checks = 0
for delta in (F(-1),F(1),F(2),F(3),F(5),F(7),F(1,2),F(2,3)):
    d = [2/delta,F(1),F(1),F(1),F(-1),F(1),F(-1)]
    for eps in product((-1,1),repeat=3):
        s = [F(0),F(0),F(0),F(1),*map(F,eps)]
        E,x,z = F(0),F(0),F(1)
        forms = [E,E+x,E-x,E+z,E-z,E+x+z,E-x-z]
        require(all(s[i]**2 == d[i]*forms[i] for i in range(7)), ("node",delta,eps))
        node_checks += 1
    # Coefficient identities before choosing a direction or solving its conic.
    for a,c in product(map(F,(-2,-1,0,1,2)),repeat=2):
        E=delta*c*c/2; x=a*a-E; z=1-E
        require(d[4]*(E-z) == 1-delta*c*c, "s4 chart")
        require(d[5]*(E+x+z) == 1+a*a-delta*c*c/2, "s5 chart")
        require(d[6]*(E-x-z) == 1+a*a-3*delta*c*c/2, "s6 chart")
        require(d[2]*(E-x) == delta*c*c-a*a, "tangent cone")

# Primitive modulo-nine obstruction for delta=3. Modulo three alone leaves
# c nonzero; modulo nine plus the equation forces all three coordinates zero mod3.
primitive_mod9 = [(a,b,c) for a,b,c in product(range(9),repeat=3)
                  if (a*a+b*b-3*c*c)%9 == 0 and any(t%3 for t in (a,b,c))]
require(primitive_mod9 == [], "delta=3 primitive mod9 hostile")

integer_norm_checks = 0
for n in range(1,501):
    require(norm_criterion(n) == (two_square(n) is not None), ("norm",n))
    integer_norm_checks += 1
rational_norm_checks = 0
for m in range(1,51):
    for n in range(1,26):
        # Direct integral norm mn yields rational norm (a/n,b/n).
        witness = two_square(m*n)
        require(norm_criterion(F(m,n)) == (witness is not None), ("rational norm",m,n))
        if witness:
            a,b = witness
            require(F(a,n)**2+F(b,n)**2 == F(m,n), "rational witness map")
        rational_norm_checks += 1
for delta in (-1,-2,-3,-5):
    require(not norm_criterion(delta), "negative norm")

parameter_checks=0
for u,v in ((F(1),F(1)),(F(1),F(2)),(F(3,5),F(4,5)),(F(1,2),F(1,2))):
    delta=u*u+v*v
    for num,den in product(range(-6,7),range(1,7)):
        t=F(num,den)
        a=(u*(1-t*t)-2*v*t)/(1+t*t)
        b=(v*(1-t*t)+2*u*t)/(1+t*t)
        require(a*a+b*b==delta, "norm-one parametrization")
        parameter_checks += 1

print(json.dumps({
    "status":"FINITE-EXACT",
    "universe":{
        "outer_twists":4096,"twist_values":[-2,-1,1,2],
        "central_base_height":5,"central_projective_base_points":len(base),
        "all_square_pair_twists_deliberately_not_decided_by_bounded_search":True,
        "integer_norm_range":[1,500],"rational_norm_grid":"m=1..50,n=1..25",
    },
    "finite_twist_histogram":finite_histogram,
    "central_point_checks":point_checks,
    "node_equation_checks":node_checks,
    "primitive_delta3_solutions_mod9":primitive_mod9,
    "integer_norm_checks":integer_norm_checks,
    "rational_norm_checks":rational_norm_checks,
    "parametrization_checks":parameter_checks,
    "gates":gates,
    "claims_not_made":["A bounded search decides arbitrary conic solubility",
                       "The central divisor determines rational resolution lifts",
                       "The finite controls prove the unbounded trichotomy"],
},indent=2,sort_keys=True))
