#!/usr/bin/env python3
"""Independent audit of the smaller-endpoint<=8 trinomial certificates.

Rebuilds coefficient polynomials from literal falling products and verifies
each stored resultant up to a nonzero constant by independent fraction-free
Sylvester determinants at enough integer points to force polynomial identity.
The interpolation bound is proved from the Sylvester degree; it is not a
bounded-parameter inference. Does not call resultant or import producer code.
"""
import ast
from fractions import Fraction
from functools import reduce
from math import comb, factorial, gcd, lcm
from pathlib import Path
import json
import sys
import sympy as sp

sys.stdout.reconfigure(newline="\n")
G,V = sp.symbols("g v")
CHECKS = 0


def require(condition,label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def expression(source):
    """Parse only arithmetic in g,v and rational constants; never eval text."""
    def rec(node):
        if isinstance(node,ast.Constant) and type(node.value) is int:
            return sp.Integer(node.value)
        if isinstance(node,ast.Name) and node.id in ("g","v"):
            return {"g":G,"v":V}[node.id]
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):
            return -rec(node.operand)
        if isinstance(node,ast.BinOp):
            a,b = rec(node.left),rec(node.right)
            if isinstance(node.op,ast.Add): return a+b
            if isinstance(node.op,ast.Sub): return a-b
            if isinstance(node.op,ast.Mult): return a*b
            if isinstance(node.op,ast.Div): return a/b
            if isinstance(node.op,ast.Pow): return a**b
        raise ValueError("non-arithmetic certificate text")
    return rec(ast.parse(source,mode="eval").body)


def member(n,A,B):
    return n >= 0 and any(A*y+B*z == n for y in range(n//A+1)
                          for z in range(n//B+1))


def literal_compressed(a,A,B,k):
    reps = [(y,z) for y in range(k*a//A+1) for z in range(k*a//B+1)
            if A*y+B*z == k*a]
    low = min(y for y,z in reps)
    polynomial = sp.Integer(0)
    for y,z in reps:
        require((y-low) % B == 0,"multiplicity exponent step")
        numerator = sp.Integer(1)
        for i in range(y+z):
            numerator *= k*G-i
        polynomial += numerator*V**((y-low)//B)/(factorial(y)*factorial(z))
    return low,sp.expand(polynomial)


def integer_coefficients(poly):
    terms = sp.Poly(poly,G,V,domain=sp.QQ).terms()
    denominator = lcm(*(int(coefficient.q) for powers,coefficient in terms))
    degree_v = max(powers[1] for powers,coefficient in terms)
    degree_g = max(powers[0] for powers,coefficient in terms)
    rows = [[0]*(degree_g+1) for _ in range(degree_v+1)]
    for (gpower,vpower),coefficient in terms:
        rows[vpower][gpower] = int(coefficient*denominator)
    return rows,degree_g


def evaluate(coefficients,x):
    answer = 0
    for coefficient in reversed(coefficients):
        answer = answer*x+coefficient
    return answer


def bareiss(entries):
    a = [list(row) for row in entries]
    n = len(a)
    if not n:
        return 1
    old,sign = 1,1
    for k in range(n-1):
        i = next((i for i in range(k,n) if a[i][k]),None)
        if i is None:
            return 0
        if i != k:
            a[k],a[i] = a[i],a[k]
            sign *= -1
        pivot = a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                value = pivot*a[i][j]-a[i][k]*a[k][j]
                require(value % old == 0,"Sylvester fraction-free division")
                a[i][j] = value//old
            a[i][k] = 0
        old = pivot
    return sign*a[-1][-1]


def sylvester(f,g):
    m,n = len(f)-1,len(g)-1
    rows = []
    for shift in range(n):
        rows.append([0]*shift+list(reversed(f))+[0]*(n-1-shift))
    for shift in range(m):
        rows.append([0]*shift+list(reversed(g))+[0]*(m-1-shift))
    return bareiss(rows)


def direct_channels(a,b,c,m):
    result = []
    for y in range(m+1):
        for z in range(m-y+1):
            x = m-y-z
            if -a*x+b*y+c*z == 0:
                result.append((x,y,z))
    return result


def main():
    path = Path("05-knowledge/results/overnight_20260906_moments_width8_certificates.json")
    data = json.loads(path.read_text(encoding="utf-8"))
    expected = {(a,A,B) for a in range(1,9) for A in range(1,9) for B in range(A+1,9)
                if gcd(A,B) == 1 and member(a-A*B,A,B)}
    found = [(r["a"],r["A"],r["B"]) for r in data["families"]]
    require(len(found) == len(set(found)) == 30,"thirty distinct family certificates")
    require(set(found) == expected,"complete independent normal-form universe")
    determinant_points = positive_coefficients = 0
    for row in data["families"]:
        a,A,B = row["a"],row["A"],row["B"]
        require(row["gmin"] == a//A+1,"strict positive-middle-charge threshold")
        stored = []
        for k,label in ((1,"first"),(2,"second")):
            low,full = literal_compressed(a,A,B,k)
            content = expression(row[label+"_content"])
            primitive = expression(row[label+"_primitive"])
            require(low == row[label+"_ymin"],"removed monomial exponent")
            require(sp.expand(full-content*primitive) == 0,"exact literal content reconstruction")
            sp.Poly(content,G,domain=sp.QQ)
            sp.Poly(primitive,G,V,domain=sp.QQ)
            # Positivity at every admissible real integer g follows from
            # y+z <= ka/A < kg, not from evaluating these checks.
            stored.append(integer_coefficients(primitive))
        (pc,dg_p),(qc,dg_q) = stored
        degree_p,degree_q = len(pc)-1,len(qc)-1
        bound = degree_q*dg_p+degree_p*dg_q
        R = row["resultant_primitive_ascending"]
        require(len(R)-1 <= bound,"Sylvester polynomial degree bound")
        require(R[-1] > 0 and reduce(gcd,(abs(v) for v in R)) == 1,"primitive positive-leading normalization")
        shift = row["gmin"]
        shifted = [sum(R[j]*comb(j,i)*shift**(j-i) for j in range(i,len(R)))
                   for i in range(len(R))]
        require(shifted == row["shifted_ascending"],"independent binomial shift identity")
        require(all(v > 0 for v in shifted),"strict coefficient positivity")
        positive_coefficients += len(shifted)
        ratio = None
        # Both polynomials have degree <=bound. Equality at bound+1
        # distinct points is therefore an exact certificate identity.
        for value in range(bound+1):
            det = sylvester([evaluate(c,value) for c in pc],[evaluate(c,value) for c in qc])
            target = evaluate(R,value)
            if target:
                if ratio is None:
                    ratio = Fraction(det,target)
                    require(ratio != 0,"nonzero resultant scalar")
                require(det == ratio*target,"exact Sylvester certificate value")
            else:
                require(det == 0,"exact zero certificate value")
            determinant_points += 1
        require(ratio is not None,"nonzero polynomial sampled beyond degree")
        print("certified",(a,A,B),"Sylvester_degree_bound",bound,"points",bound+1)
    require(positive_coefficients == data["shifted_coefficient_count"] == 352,"all shifted coefficients")
    exceptions = []
    for c in range(2,9):
        for b in range(1,c):
            for a in range(9,(c-1)*(c-2)-c+1):
                if gcd(a,gcd(b,c)) != 1:
                    continue
                bound = min((a+b)//gcd(a,b),(a+c)//gcd(a,c))
                rows = next(rr for m in range(1,bound+1) if (rr:=direct_channels(a,b,c,m)))
                if len(rows) > 1:
                    exceptions.append((-a,b,c))
    require(exceptions == [tuple(r["support"]) for r in data["exceptional_supports"]],
            "exhaustive opposite-endpoint original-charge universe")
    for record in data["exceptional_supports"]:
        a,b,c = -record["support"][0],record["support"][1],record["support"][2]
        polynomials = []
        for k,label in ((1,"first"),(2,"second")):
            mass = k*record["g"]
            rows = direct_channels(a,b,c,mass)
            direct = sum(sp.Rational(factorial(mass),factorial(x)*factorial(y)*factorial(z))*V**y
                         for x,y,z in rows)
            require(sp.expand(direct-expression(record[label+"_raw_moment"])) == 0,
                    "exception literal multinomial identity")
            polynomials.append(sp.Poly(direct,V,domain=sp.QQ))
        common = sp.gcd(*polynomials)
        require(len(common.terms()) == 1,"exception torus saturation")
    print("PASS: complete smaller-endpoint degree<=8 certificate audit")
    print("families:",len(found),"; positive shifted coefficients:",positive_coefficients)
    print("independent integer Sylvester evaluations:",determinant_points)
    print("opposite-endpoint supports:",exceptions)
    print("exact_checks:",CHECKS)
    print("Parameter g is unbounded; degree-bounded interpolation verifies full polynomial identities.")


if __name__ == "__main__":
    main()
