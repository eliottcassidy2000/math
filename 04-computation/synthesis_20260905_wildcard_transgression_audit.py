#!/usr/bin/env python3
"""Independent literal truncated-ring audit of root's collision curvature.

Uses Fraction arithmetic, compiler substitution in Q[s]/(s^3), and direct
linear solution. Does not import root's script, SymPy, or differential jets.
Four affine-basis coefficient tests recover the complete a,b,c dependence;
additional mixed points are hostile controls. Affineness follows separately
from P having double zeros at the three endpoints and from the compiler's
fixed first jets there. This is a finite exact arithmetic audit supporting
the symbolic proof, not an independent universal deformation theorem.
"""
from fractions import Fraction as F
import sys

sys.stdout.reconfigure(newline="\n")
CHECKS = 0


def require(value, label):
    global CHECKS
    CHECKS += 1
    if not value:
        raise RuntimeError(label)


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def scale(a, b):
    return tuple(b*x for x in a)


def mul(a, b):
    return tuple(sum(a[i]*b[j-i] for i in range(j+1)) for j in range(3))


ONE = (F(1), F(0), F(0))


def power(a, n):
    result = ONE
    for _ in range(n):
        result = mul(result, a)
    return result


def polynomial(coefficients, x):
    result = (F(0), F(0), F(0))
    for degree, coefficient in enumerate(coefficients):
        result = add(result, scale(power(x, degree), coefficient))
    return result


def residual(a, b, c, h, etas):
    result = []
    for z, eta in zip((-1, 0, 1), etas):
        x = (F(z), F(eta), F(0))
        x2 = mul(x, x)
        P = mul(x2, power(add(x2, scale(ONE, -1)), 2))
        Q1 = polynomial((F(-3,4), 1, F(-27,4), -2, F(9,2), 1), x)
        Q = add(Q1, mul(P, polynomial((a, b, c), x)))
        q = add(Q, mul((F(0), F(1), F(0)), polynomial(h, x)))
        D = add(ONE, mul(x2, q))
        C = mul(mul(x, D), add(D, scale(ONE, 2)))
        E = mul(q, add(D, scale(ONE, 3)))
        require(C[0] == 0 and E[0] == -3, "literal base collision")
        result.extend((C[2], E[2]))
    return result


def solve(matrix, rhs):
    n = len(rhs)
    a = [[F(v) for v in row]+[F(rhs[i])] for i, row in enumerate(matrix)]
    determinant = F(1)
    for k in range(n):
        i = next(i for i in range(k, n) if a[i][k])
        if i != k:
            a[k], a[i] = a[i], a[k]
            determinant *= -1
        pivot = a[k][k]
        determinant *= pivot
        a[k] = [v/pivot for v in a[k]]
        for i in range(n):
            if i != k:
                factor = a[i][k]
                a[i] = [x-factor*y for x, y in zip(a[i], a[k])]
    return tuple(row[-1] for row in a), determinant


def main():
    J = ((3,0,0,-1,0,-2),(-9,0,0,0,-1,2),
         (0,3,0,-1,0,0),(0,4,0,0,-1,0),
         (0,0,3,-1,0,-2),(0,0,9,0,-1,-2))
    cases = ((0,0,0),(1,0,0),(0,1,0),(0,0,1),
             (-3,7,11),(17,-23,5),(F(2,3),F(-7,5),F(13,11)))
    for a, b, c in cases:
        for label, h, etas, expected in (
            ("constant", (1,), (F(-2,3),0,F(2,3)),
             -(36*F(a)+16*F(b)+36*F(c)+259)/9),
            ("quadratic", (0,-9,4), (F(-14,3),4,F(2,3)),
             -(3072*F(a)-7424*F(b)+8256*F(c)+102511)/144),
        ):
            values = residual(a, b, c, h, etas)
            solution, determinant = solve(J, tuple(-v for v in values))
            require(determinant == -288, "formal Jacobian determinant")
            require(solution[5] == expected, ("literal curvature formula", label,a,b,c))
            print(label, (str(a),str(b),str(c)), "c2="+str(solution[5]))
    print("PASS: independent literal truncated-ring curvature and determinant")
    print("exact_checks:", CHECKS)


if __name__ == "__main__":
    main()
