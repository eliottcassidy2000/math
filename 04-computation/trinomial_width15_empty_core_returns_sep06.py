#!/usr/bin/env python3
"""Exact unbounded certificate for (-15,2g-15,3g-15), g>=8, gcd(g,15)=1.

Symbolic trace/norm positivity is the proof; seven named exact Laurent
controls corroborate it. No repository producer is imported. Requires SymPy.
"""
import hashlib
import json
from math import factorial, gcd

import sympy as S

T, G, Z = S.symbols("tau g s")
CHECKS = 0


def require(test, label):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(label)


def falling(value, degree):
    return S.prod(value-k for k in range(degree))


def raw_channels(a, b, c, mass):
    result = []
    for x in range(mass+1):
        numerator = a*x-b*(mass-x)
        if numerator % (c-b):
            continue
        z = numerator//(c-b)
        y = mass-x-z
        if min(y, z) >= 0:
            result.append(((x,y,z), factorial(mass)//
                           (factorial(x)*factorial(y)*factorial(z))))
    return sorted(result, key=lambda item: item[0][2])


def literal_rows(a, b, c, maximum):
    """Direct multiplication with gamma degree retained independently."""
    state = {(0,0): 1}
    result = []
    for _ in range(maximum):
        following = {}
        for (charge,z), coefficient in state.items():
            for delta,dz in ((-a,0),(b,0),(c,1)):
                key = charge+delta,z+dz
                following[key] = following.get(key,0)+coefficient
        state = following
        result.append({z:v for (charge,z),v in state.items() if charge == 0})
    return result


def main():
    p = S.Poly(6*T**2+20*(G-5)*T+(G-5)*(G-6), T)
    content = falling(G,5)/720
    first = S.Poly(sum(falling(G,7-j)*T**j/
                       (factorial(6-3*j)*factorial(1+2*j)) for j in range(3)), T)
    q = S.Poly(sum(falling(2*G,15-j)*T**j/
                   (factorial(15-3*j)*factorial(2*j)) for j in range(6)), T)
    require(S.cancel(first.as_expr()-content*p.as_expr()) == 0, "first positive content")
    inverse = -(6*T+20*(G-5))/((G-5)*(G-6))
    require(S.cancel(S.rem(T*inverse-1,p.as_expr(),T)) == 0, "exact inverse carry")
    remainder = S.Poly(S.rem(q.as_expr()*inverse,p.as_expr(),T), T)

    k = G*(G-4)*(G-3)*(G-2)*(G-1)*(2*G-9)*(2*G-7)*(2*G-5)*(2*G-3)*(2*G-1)
    cpoly = 1106459*G**3-17466623*G**2+91449662*G-158894736
    dpoly = 27352253*G**3-404127357*G**2+1990079478*G-3266244982
    c = k*(G-5)*cpoly/12259447200
    d = k*dpoly/15324309000
    require(remainder.degree() == 1 and
            S.cancel(remainder.nth(0)-c) == 0 and S.cancel(remainder.nth(1)-d) == 0,
            "complete anchored symbolic remainder")
    j = 106089635*G**3-1564109559*G**2+7685968926*G-12588295720
    h = (8051838961249*G**7-295063248863031*G**6+4630220682540559*G**5
         -40333487193583421*G**4+210638390555236696*G**3
         -659512619763975524*G**2+1146316900083491136*G-853264618035013184)
    trace = 2*c-S.Rational(10,3)*(G-5)*d
    norm = c*c-S.Rational(10,3)*(G-5)*c*d+(G-5)*(G-6)*d*d/6
    expected_trace = -k*(G-5)*j/18389170800
    expected_norm = k*k*(G-5)*h/3757351141239696000000
    require(S.cancel(trace-expected_trace) == 0, "exact quotient trace factorization")
    require(S.cancel(norm-expected_norm) == 0, "exact quotient norm factorization")
    require(S.cancel(S.discriminant(p.as_expr(),T)-8*(G-5)*(47*G-232)) == 0,
            "distinct negative roots at g>=8")
    shifts = {}
    for name, polynomial in (("J",j),("H",h)):
        coefficients = S.Poly(S.expand(polynomial.subs(G,Z+8)),Z).all_coeffs()
        require(all(value > 0 for value in coefficients), "shift positivity "+name)
        shifts[name] = [str(value) for value in coefficients]
    # Confirm all other factors in both signs are strictly positive for s>=0.
    for factor in (G-5,G-4,G-3,G-2,G-1,G,2*G-9,2*G-7,2*G-5,2*G-3,2*G-1):
        require(all(value > 0 for value in S.Poly(factor.subs(G,Z+8),Z).all_coeffs()),
                "positive prefactor at g>=8")
    print("P_normalized",p.as_expr())
    print("first_content",S.factor(content))
    print("anchored_R_constant",S.factor(c))
    print("anchored_R_slope",S.factor(d))
    print("trace",S.factor(expected_trace))
    print("norm",S.factor(expected_norm))
    print("strict_positive_shifts",json.dumps(shifts,sort_keys=True))

    manifest = []
    for g in (8,11,13,16,19,23,31):
        a,b,c_charge = 15,2*g-15,3*g-15
        require(gcd(g,15) == 1 and gcd(a,gcd(b,c_charge)) == 1, "primitive control")
        native = literal_rows(a,b,c_charge,2*g)
        require(all(not row for row in native[:g-1]), "first support return mass")
        raw = [raw_channels(a,b,c_charge,mass) for mass in (g,2*g)]
        for mass,channels in zip((g,2*g),raw):
            require(native[mass-1] == {v[2]:weight for v,weight in channels},
                    "literal multiplication versus original charge enumeration")
        require([v for v,_ in raw[0]] == [(g-7+j,6-3*j,1+2*j) for j in range(3)],
                "complete three-channel first row")
        require([v for v,_ in raw[1]] == [(2*g-15+j,15-3*j,2*j) for j in range(6)],
                "complete six-channel second row with lower carry")
        pf,qf = [S.Poly(sum(weight*T**i for i,(_,weight) in enumerate(row)),T)
                 for row in raw]
        require(S.cancel(pf.as_expr()-first.as_expr().subs(G,g)) == 0 and
                S.cancel(qf.as_expr()-q.as_expr().subs(G,g)) == 0, "symbolic rows specialize exactly")
        rf = S.Poly(S.rem(qf.as_expr()*S.invert(T,pf.as_expr(),T),pf.as_expr(),T),T)
        require(S.cancel(rf.as_expr()-remainder.as_expr().subs(G,g)) == 0,
                "independent specialized inverse/remainder")
        # Reconstruct the multiplication matrix directly, not using c,d traces.
        multiplication = S.zeros(2)
        for column in range(2):
            entry = S.Poly(S.rem(T**(column+1),pf.as_expr(),T),T)
            for row in range(2):
                multiplication[row,column] = entry.nth(row)
        response = rf.nth(0)*S.eye(2)+rf.nth(1)*multiplication
        matrix = S.Matrix(2,2,lambda i,j:S.trace(response*multiplication**(i+j)))
        require(S.trace(response) == expected_trace.subs(G,g), "literal quotient trace")
        require(response.det() == expected_norm.subs(G,g), "literal quotient norm")
        require(matrix[0,0] < 0 and matrix.det() > 0, "negative definite Hermite form")
        require(S.gcd(pf,qf).degree() == 0, "coprime complete rows")
        require(pf.count_roots(-S.oo,0) == 2 and S.gcd(pf,pf.diff()).degree() == 0,
                "two simple attainable negative scalar roots")
        if g == 8:
            require(rf.as_expr() == 4743488+47087024*T, "named carry correction")
            require(matrix[0,0] == -461383264 and matrix.det() == 587632760217600,
                    "named Hermite minors")
        print("named_control",(-a,b,c_charge),"first_mass",g,"channels",[len(row) for row in raw],
              "negative_trace",S.trace(response)<0,"positive_norm",response.det()>0)
        manifest.append([g,[str(weight) for _,weight in raw[0]],
                         [str(weight) for _,weight in raw[1]],str(rf.as_expr()),
                         str(S.trace(response)),str(response.det())])
    # Hostile to dropping the primitive hypothesis, not to the mass-row identity.
    hostile = literal_rows(15,3,12,9)
    first_mass = next(i for i,row in enumerate(hostile,1) if row)
    require(first_mass == 3 and first_mass != 9, "excluded nonprimitive first-return hostile")
    print("nonprimitive_hostile (-15,3,12), formal g=9, actual first support return=3")
    print("PROVED: g>=8 integral, gcd(g,15)=1; first nonzero return g or 2g, both attainable")
    print("SCOPE: this fixed (a,A,B)=(15,2,3) family; no endpoint-fifteen strip or general closure")
    encoded = json.dumps({"shifts":shifts,"controls":manifest},sort_keys=True,
                         separators=(",",":")).encode()
    print("semantic_sha256",hashlib.sha256(encoded).hexdigest())
    print("explicit_checks",CHECKS)


if __name__ == "__main__":
    main()
