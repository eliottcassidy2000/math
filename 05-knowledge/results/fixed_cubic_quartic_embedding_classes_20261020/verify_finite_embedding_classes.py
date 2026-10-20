#!/usr/bin/env python3
"""Exact controls for the fixed <3,4> finite-embedding-class theorem."""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import sympy as sp

y, t, x, u, v = sp.symbols("y t x u v")
a, b, c, beta, gamma = sp.symbols("a b c beta gamma")
checks = {}


def check(name, condition):
    if not bool(condition):
        raise RuntimeError(name)
    checks[name] = True


def zero(expr):
    return sp.cancel(sp.expand(expr)) == 0


def rem(expr, divisor, var=y):
    return sp.rem(sp.Poly(expr,var),sp.Poly(divisor,var)).as_expr()


def degree_gaps_from_basis(degrees, n):
    residues = [d % n for d in degrees]
    check(f"basis:{degrees}:distinct_residues", sorted(residues) == list(range(n)))
    gaps = []
    for d in degrees:
        gaps.extend(range(d % n,d,n))
    return sorted(gaps)


def characteristic_patterns():
    allowed = [3,4,6,7]
    result = []
    for length in range(2,5):
        for seq in itertools.product(allowed,repeat=length):
            if seq[0] <= seq[1] or seq[0] % seq[1] == 0:
                continue
            if any(z > 5 for z in seq[2:]) or not {3,4}.issubset(seq):
                continue
            gcds = [seq[0]]
            for z in seq[1:]:
                gcds.append(math.gcd(gcds[-1],z))
            ns = [gcds[i-1]//gcds[i] for i in range(1,length)]
            if gcds[-1] != 1 or any(z <= 1 for z in ns):
                continue
            if any(seq[i] >= seq[i-1]*ns[i-2] for i in range(2,length)):
                continue
            conductor = 1-seq[0]+sum((ns[i-1]-1)*seq[i] for i in range(1,length))
            if conductor != 6:
                continue
            result.append(list(seq))
    return sorted(result)


def monogenic_example(p,q,name,expected_delta):
    n = sp.degree(p,y)
    relation = p-u
    def column(expr):
        poly = sp.Poly(expr,y)
        return sp.Matrix([poly.nth(i) for i in range(n)])
    M = sp.Matrix.hstack(*(column(rem(q*y**j,relation)) for j in range(n)))
    H = (v*sp.eye(n)-M).det().expand()
    gcd = sp.gcd(sp.Poly(H,v,domain=sp.QQ.frac_field(u)),
                 sp.Poly(sp.diff(H,v),v,domain=sp.QQ.frac_field(u)))
    check(name+":primitive_characteristic", gcd.degree() == 0)
    conductor = sp.cancel(sp.diff(H,v).subs({u:p,v:q})/sp.diff(p,y))
    check(name+":polynomial_conductor", sp.denom(conductor) == 1)
    check(name+":delta_from_conductor", sp.degree(conductor,y) == 2*expected_delta)
    U = sp.cancel((p.subs(y,t)-p)/(t-y)).expand()
    V = sp.cancel((q.subs(y,t)-q)/(t-y)).expand()
    V = rem(V,U,t)
    check(name+":divided_difference_conductor",
          zero(sp.resultant(U,V,t)-conductor))
    return {"p":str(p),"q":str(q),"conductor":str(sp.factor(conductor)),
            "delta":expected_delta}


def run():
    h = y**3+a*y
    k = y**4+b*y**2+c*y
    p = k+beta*h
    q = h**2+gamma*h
    r = sp.Poly(sp.expand(q*q-p**3),y)
    check("general:degree12_cancels", r.nth(12) == 0)
    check("general:degree11_is_minus3beta", r.nth(11) == -3*beta)
    rho = 4*a-3*b
    r0 = sp.Poly(sp.expand(q*q-k**3-rho*k*q),y)
    for degree in (12,11,10):
        check(f"general:corrected_degree{degree}_zero", r0.nth(degree) == 0)
    check("general:corrected_degree9", r0.nth(9) == 2*gamma-3*c)
    check("general:pq_no_degree9", sp.Poly(k*q,y).nth(9) == 0)

    # Rows are 1,q,q²,q³; columns are 1,q,r,qr.
    P, Rho = sp.symbols("P Rho")
    change = sp.Matrix([[1,0,-P**3,0],[0,1,-Rho*P,-P**3],
                        [0,0,1,-Rho*P],[0,0,0,1]])
    check("basis:unimodular_change", change.det() == 1)
    Q = sp.symbols("Q")
    old = sp.Matrix([[1,Q,Q**2,Q**3]])
    new = old*change
    rformal = Q**2-Rho*P*Q-P**3
    check("basis:third_column", zero(new[2]-rformal))
    check("basis:fourth_column", zero(new[3]-Q*rformal))

    qstar = h*h+sp.Rational(3,2)*c*h
    eta = -2*a*(a-b)
    rstar = sp.Poly(sp.expand(qstar*qstar-k**3-rho*k*qstar-eta*k*k),y)
    for d in range(8,13):
        check(f"candidate:degree{d}_zero", rstar.nth(d) == 0)
    check("candidate:degree7", sp.factor(rstar.nth(7)) == -c*(2*a-3*b)/2)
    kstar = k.subs(b,sp.Rational(2,3)*a)
    recovery = qstar*qstar-kstar**3-2*a*kstar*qstar+sp.Rational(2,3)*a*a*kstar*kstar\
               +(2*a**3/sp.Integer(27)+3*c*c/sp.Integer(4))*qstar\
               -(a**4/sp.Integer(9)+a*c*c)*kstar-c**3*h/sp.Integer(8)
    check("candidate:universal_recovery_identity", zero(recovery))
    gaps5 = degree_gaps_from_basis([0,6,7,13],4)
    check("candidate:delta5", len(gaps5) == 5)
    check("candidate:gaps5", gaps5 == [1,2,3,5,9])

    gaps7 = degree_gaps_from_basis([0,6,11,17],4)
    gaps6 = degree_gaps_from_basis([0,6,9,15],4)
    check("basis:delta7", len(gaps7) == 7)
    check("basis:delta6", len(gaps6) == 6)
    check("basis:gaps7", gaps7 == [1,2,3,5,7,9,13])
    check("basis:gaps6", gaps6 == [1,2,3,5,7,11])
    check("basis:delta7_formula", (sum([0,6,11,17])-sum(range(4)))//4 == 7)
    check("basis:delta6_formula", (sum([0,6,9,15])-sum(range(4)))//4 == 6)
    patterns = characteristic_patterns()
    check("semigroup:reduced_patterns", patterns == [[4,3],[6,4,3]])

    H,K = y**3,y**4+y
    Qsharp = H**2+sp.Rational(3,2)*H
    check("sharp:generates", zero(8*(Qsharp**2-K**3)+6*Qsharp-H))
    Z = 8*Qsharp+3
    check("sharp:previous_second_embedding", zero(Z-(8*y**6+12*y**3+3)))
    check("candidate:need_not_generate",
          sp.expand(((y**3+y)**2).subs(y,-y)-(y**3+y)**2) == 0
          and (y**4).subs(y,-y) == y**4
          and (y**3+y).subs(y,-y) != y**3+y)

    examples = []
    for av,bv,cv,be,ga,delta in [
        (0,0,1,1,0,7),
        (1,2,3,2,1,7),
        (0,0,1,0,0,6),
        (1,2,3,0,2,6),
        (0,0,1,0,sp.Rational(3,2),3),
        (0,0,2,0,3,3),
        (1,2,3,0,sp.Rational(9,2),5),
    ]:
        HH = h.subs(a,av)
        KK = k.subs({b:bv,c:cv})
        pp = KK+be*HH
        qq = HH**2+ga*HH
        examples.append(monogenic_example(pp,qq,f"example:{len(examples)}",delta))

    # The ramified cusp has no distinct normalization collision, but the
    # chain-rule obstruction forces the Jacobian to vanish along y=0.
    cuspP = y**3+x*y**4+x**2
    cuspQ = y**4+x*y**6+x
    check("ramification:cusp_py_zero", sp.diff(cuspP,y).subs(y,0) == 0)
    check("ramification:cusp_qy_zero", sp.diff(cuspQ,y).subs(y,0) == 0)
    cuspJ = sp.diff(cuspP,x)*sp.diff(cuspQ,y)-sp.diff(cuspP,y)*sp.diff(cuspQ,x)
    check("ramification:cusp_jacobian_zero", cuspJ.subs(y,0) == 0)

    # Entire generating families obtained by ambient R-automorphisms.
    # Each integral has the same endpoint difference as its reference pair.
    period_records = []
    for index,(P0,Q0,reference) in enumerate([
        (H,K,sp.Rational(9,28)),
        (H+x,K,sp.Rational(9,28)),
        (H+x*K**2,K,sp.Rational(9,28)),
        (H,K+x*H**2,sp.Rational(9,28)),
        (K,Z,-sp.Rational(9,5)),
        (K+x,Z,-sp.Rational(9,5)),
        (K+x*Z**2,Z,-sp.Rational(9,5)),
        (K,Z+x*K**2,-sp.Rational(9,5)),
    ]):
        L = sp.integrate(sp.expand(P0*sp.diff(Q0,y)),y)
        reduced = sp.expand(rem(L,y**3+1))
        check(f"period:{index}:linear_endpoint",
              zero(sp.diff(reduced,y)-reference))
        check(f"period:{index}:x_independent_difference",
              zero(sp.diff(reduced,x,y)))
        period_records.append(str(reduced))
    return {
        "status":"PASS: exact controls, not proof of the cited semigroup theorem",
        "sympy_version":sp.__version__,
        "check_count":len(checks),
        "checks":checks,
        "basis_gap_sets":{"delta7":gaps7,"delta6":gaps6,"delta5":gaps5},
        "reduced_patterns":patterns,
        "examples":examples,
        "fixed_triple_period_remainders":period_records,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output",type=Path)
    args = parser.parse_args()
    output = json.dumps(run(),sort_keys=True,indent=2)+"\n"
    if args.output:
        args.output.write_text(output)
    print(output,end="")
