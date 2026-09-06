#!/usr/bin/env python3
"""Independent bounded exact audit of THM-4443; no producer/referee imports.

Literal Hasse matrices, rational elimination and integer Smith form are
compared with reciprocal coefficients obtained by polynomial convolution
and power-series inversion. All mathematical gates survive python -O.
"""
from fractions import Fraction as Q
from hashlib import sha256
from math import comb, lcm, prod
from pathlib import Path
import subprocess
import sys
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

sys.stdout.reconfigure(newline="\n")
REPO = Path(__file__).resolve().parents[1]
COMMIT = "058a8ded98cfa25fd90f061efe5d54903c9b7379"
GATES = 0
CASES = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def valuation(x, p=2):
    x = Q(x)
    if not x:
        return float("inf")
    a, b, out = abs(x.numerator), x.denominator, 0
    while a % p == 0:
        a //= p
        out += 1
    while b % p == 0:
        b //= p
        out -= 1
    return out


def multiply(a, b):
    out = [Q(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def local_series(nodes, mult):
    result = []
    for i, (x, m) in enumerate(zip(nodes, mult)):
        q = [Q(1)]
        for j, (y, k) in enumerate(zip(nodes, mult)):
            if i != j:
                q = multiply(q, [Q(comb(k, r) * (x-y)**(k-r)) for r in range(k+1)])
        inverse = [1 / q[0]]
        for r in range(1, m):
            inverse.append(-sum(q[s]*inverse[r-s] for s in range(1,min(r,len(q)-1)+1))/q[0])
        for r in range(m):
            need(sum(q[s]*inverse[r-s] for s in range(min(r,len(q)-1)+1)) == (r == 0),
                 ("reciprocal coefficient recurrence", nodes, mult, i, r))
        result.append(inverse)
    return result


def gauss_inverse(H):
    n = len(H)
    rows = [[Q(x) for x in row] + [Q(i == j) for j in range(n)] for i, row in enumerate(H)]
    for k in range(n):
        pivot = next(i for i in range(k,n) if rows[i][k])
        rows[k], rows[pivot] = rows[pivot], rows[k]
        c = rows[k][k]
        rows[k] = [x/c for x in rows[k]]
        for i in range(n):
            if i != k and rows[i][k]:
                c = rows[i][k]
                rows[i] = [x-c*y for x,y in zip(rows[i],rows[k])]
    return [row[n:] for row in rows]


def matrix_case(label, nodes, mult, expected_loss=None):
    global CASES
    n = sum(mult)
    H = [[comb(j,r)*x**(j-r) if j >= r else 0 for j in range(n)]
         for x,m in zip(nodes,mult) for r in range(m)]
    inv = gauss_inverse(H)
    for i in range(n):
        for j in range(n):
            need(sum(H[i][k]*inv[k][j] for k in range(n)) == (i == j), ("literal inverse",label,i,j))
    series = local_series(nodes,mult)
    c = 0
    for m, row in zip(mult,series):
        for r in range(m):
            need(inv[-1][c] == row[m-r-1], ("all proposed denominators attained in top row",label,c))
            c += 1
    denom_inv = lcm(*(x.denominator for row in inv for x in row))
    denom_series = lcm(*(x.denominator for row in series for x in row))
    smith = smith_normal_form(Matrix(H),domain=ZZ)
    factors = [abs(int(smith[i,i])) for i in range(n)]
    need(all(factors[i+1] % factors[i] == 0 for i in range(n-1)), ("Smith factor divisibility",label))
    need(denom_inv == denom_series == factors[-1], ("entire largest integer factor",label))
    determinant = prod(abs(nodes[i]-nodes[j])**(mult[i]*mult[j])
                       for i in range(len(nodes)) for j in range(i))
    need(prod(factors) == determinant, ("Hasse determinant has no ordinary-derivative factorial",label))
    exponents = tuple(valuation(x) for x in factors)
    if expected_loss is not None:
        need(exponents[-1] == expected_loss, ("predicted dyadic sharp loss",label,exponents,expected_loss))
    CASES += 1
    return exponents, series


def gamma(t):
    t %= 16
    if t % 4 == 3:
        return 3
    if t % 8 == 5:
        return 2
    return 1 if t == 1 else 0


def unit_law():
    for t in range(-127,128,2):
        nu = min(valuation(t*t+3*t+4),valuation(t*t-7*t+14))
        need(nu == 4-gamma(t), ("complete polynomial residue control",t,nu))
        need(gamma(t) == gamma(2-t), ("closest endpoint exchange",t))
    for e in (0,2,5):
        for t in range(1,16,2):
            nodes = tuple(2**e*x for x in (0,2,t))
            expected = max(7*e+4,8*e+gamma(t))
            matrix_case(("one-level",e,t),nodes,(3,3,3),expected)
    for e in (0,1,3):
        for d in (2,3,5):
            for t in (1,3):
                nodes = tuple(2**e*x for x in (0,2**d,t))
                matrix_case(("deep-pair",e,d,t),nodes,(3,3,3),8*e+5*d-1)
    for t in (1,3,5,9):
        original = tuple(4*x for x in (0,2,t))
        moved = tuple(7-3*x for x in original)
        matrix_case(("signed-unit-affine",t),moved,(3,3,3),max(18,16+gamma(t)))
    print("DYADIC_UNIT_LAW all8 residues at e=0,2,5;18 deep-pair cases;4 signed unit-affine controls")
    print("SIMULTANEOUS_NUMERATOR_VALUATION minima1,2,3,4; Gamma3,2,1,0; endpoint exchange retained")


def hostiles():
    for e in (0,1,2,4):
        A = tuple(2**e*x for x in (0,1,2))
        B = tuple(2**e*x for x in (1,0,3))
        need(tuple(valuation(A[i]-A[j]) for i in range(3) for j in range(i)) ==
             tuple(valuation(B[i]-B[j]) for i in range(3) for j in range(i)), "labelled uniform isometry")
        sa,_ = matrix_case(("uniform-A",e),A,(3,3,3),max(7*e+4,8*e+1))
        sb,_ = matrix_case(("uniform-B",e),B,(3,3,3),max(7*e+4,8*e+3))
        print("UNIFORM_HOSTILE e",e,"A",sa,"B",sb)
        if e == 2:
            need(sa == (0,0,0,2,6,8,11,18,18) and sb == (0,0,0,2,6,8,11,17,19), "full 9x9 hostile partitions")
    for e in (0,1,2,5):
        A = tuple(2**e*x for x in (0,2,1))
        B = tuple(2**e*x for x in (1,3,0))
        need(tuple(valuation(A[i]-A[j]) for i in range(3) for j in range(i)) ==
             tuple(valuation(B[i]-B[j]) for i in range(3) for j in range(i)), "multiplicity-labelled nonuniform isometry")
        sa,_ = matrix_case(("mixed-A",e),A,(2,2,1),max(3*e+2,4*e+1))
        sb,_ = matrix_case(("mixed-B",e),B,(2,2,1),max(3*e+2,4*e))
        print("NONUNIFORM_HOSTILE e",e,"A",sa,"B",sb)


def boundary_controls():
    for m in (1,2,5):
        smith,_ = matrix_case(("singleton",m),(-7,),(m,),0)
        need(smith == (0,)*m,"singleton integral translated monomials")
    for k in range(1,5):
        nodes = (-3,9)
        smith,series = matrix_case(("two-node",k),nodes,(k,k))
        for p in (2,3):
            expected = max((k+r)*valuation(12,p)-valuation(comb(k+r-1,r),p) for r in range(k))
            observed = max(-valuation(x,p) for row in series for x in row)
            need(observed == expected, ("two-node uniform metric-only formula",k,p))
    for nodes,mult in (((-5,2,11),(3,1,2)),((0,3,8,13),(1,2,3,1))):
        matrix_case(("signed heterogeneous",nodes,mult),nodes,mult)
    # Same three-node uniform pair, now complete two-jet observations.
    for e in (0,2,5):
        A = tuple(2**e*x for x in (0,1,2))
        B = tuple(2**e*x for x in (1,0,3))
        a,_ = matrix_case(("two-jet-metric-A",e),A,(2,2,2))
        b,_ = matrix_case(("two-jet-metric-B",e),B,(2,2,2))
        need(a[-1] == b[-1] == max(4*e+2,5*e+2), "old complete two-jet terminal law survives")
    print("BOUNDARY_CONTROLS singleton; two-node k1..4 p2/3; signed heterogeneous; complete two-jet twin controls")


def main():
    print("SCOPE independent THM-4443 written-proof audit plus bounded full integer Hasse inverse/Smith controls")
    print("INCOMING_COMMIT",COMMIT)
    for path in ("01-canon/theorems/THM-4443-arbitrary-jet-precision-and-dyadic-unit-boundary.md",
                 "05-knowledge/results/hermite-higher-jet-unit-boundary-overnight-hexagon-sep05.md"):
        data = subprocess.check_output(["git","show",COMMIT+":"+path],cwd=REPO)
        print("INPUT_SHA256",sha256(data).hexdigest(),path)
    unit_law()
    hostiles()
    boundary_controls()
    print("EXACT_MATRIX_CASES",CASES)
    print("PASS_OPTIMIZATION_LIVE_GATES",GATES)


if __name__ == "__main__":
    main()
