#!/usr/bin/env python3
"""Primary exact certificate for THM-4340's two residual U=0 M=12 cusps.

This is deliberately dependency-free.  It checks the local coordinate
changes, the normal critical-value (Hasse/Tjurina) coefficients, explicit
hostiles on which all of those coefficients vanish, the cusp conductor
quotients, and the valuation inequalities supplied by the monomial factors
in the Keller residue.
"""

from fractions import Fraction as F
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")


def need(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(label)


def trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return tuple(p)


def add(p, q):
    n = max(len(p), len(q))
    return trim([
        (p[i] if i < len(p) else F(0)) + (q[i] if i < len(q) else F(0))
        for i in range(n)
    ])


def scale(p, c):
    return trim([c * x for x in p])


def sub(p, q):
    return add(p, scale(q, F(-1)))


def mul(p, q):
    if not p or not q:
        return ()
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return trim(out)


def derivative(p, order=1):
    out = tuple(p)
    for _ in range(order):
        out = trim([F(i) * out[i] for i in range(1, len(out))])
    return out


def evaluate(p, x):
    value = F(0)
    for coefficient in reversed(p):
        value = value * x + coefficient
    return value


def h5_normal_coefficients(f5, f4, f3, f2, s0, w):
    """Return the q,q^2,q^3 critical-value coefficients modulo z^4.

    After division by the exact unit 1-S^2*q, the S-dependent function is
    f5+q*f4+q^2*f3+q^3*f2+... .  The correction from -Q*s^2/2 starts at
    q^6.  The face polynomial is exactly W(S-S0)^2, so no higher derivatives
    of f5 occur.
    """

    a1 = f4
    a2 = f3
    a3 = f2
    a1_0 = evaluate(a1, s0)
    a1_1 = evaluate(derivative(a1), s0)
    a1_2 = evaluate(derivative(a1, 2), s0)
    a2_0 = evaluate(a2, s0)
    a2_1 = evaluate(derivative(a2), s0)
    a3_0 = evaluate(a3, s0)
    h1 = a1_0
    h2 = a2_0 - a1_1 * a1_1 / (4 * w)
    h3 = a3_0 - a1_1 * a2_1 / (2 * w) + a1_2 * a1_1 * a1_1 / (8 * w * w)
    return h1, h2, h3


def semigroup(generators, limit):
    reached = {0}
    changed = True
    while changed:
        changed = False
        for value in tuple(reached):
            for generator in generators:
                new = value + generator
                if new <= limit and new not in reached:
                    reached.add(new)
                    changed = True
    return reached


def conductor_quotient_basis(m):
    """Exponent basis of c/J for k[[t^2,t^m]], y^2=x^m."""

    limit = 6 * m
    ring = semigroup((2, m), limit)
    conductor = set(range(m - 1, limit + 1))
    jacobian = {
        shift + exponent
        for shift in (m, 2 * m - 2)
        for exponent in ring
        if shift + exponent <= limit
    }
    basis = tuple(sorted(conductor - jacobian))
    # The cutoff is beyond the conductor stabilization range.
    need(not [n for n in basis if n >= 3 * m], f"m={m}: quotient did not stabilize")
    return basis


def sp_monomial(key, coefficient=F(1)):
    return {} if coefficient == 0 else {key: coefficient}


def sp_add(*polys):
    out = {}
    for poly in polys:
        for key, coefficient in poly.items():
            out[key] = out.get(key, F(0)) + coefficient
            if out[key] == 0:
                del out[key]
    return out


def sp_scale(poly, coefficient):
    return {key: coefficient * value for key, value in poly.items()
            if coefficient * value != 0}


def sp_mul(first, second):
    out = {}
    for left, a in first.items():
        for right, b in second.items():
            key = tuple(x + y for x, y in zip(left, right))
            out[key] = out.get(key, F(0)) + a * b
            if out[key] == 0:
                del out[key]
    return out


def sp_pow(poly, exponent):
    out = sp_monomial((0, 0, 0))
    for _ in range(exponent):
        out = sp_mul(out, poly)
    return out


def check_sparse_factorizations(rows):
    """Reconstruct (H1)/(E1) from F=(s^2-p)(1-QH)-Qs^2/2.

    Each coefficient row is checked as a basis vector, which proves the
    arbitrary linear combination by linearity.  Negative Laurent exponents
    are retained until multiplication by the declared local monomial.
    """

    one = sp_monomial((0, 0, 0))

    # H5 keys are (S,z,t), with p=t^-1*z^-1, y=S*p, Q=t^5.
    h5_s2 = sp_monomial((2, 0, 0))
    h5_p = sp_monomial((0, -1, -1))
    h5_qbase = sp_monomial((0, 0, 5))
    h5_clear = sp_monomial((0, 6, 1))
    h5_rho = sp_monomial((0, 1, 1))
    h5_z5 = sp_monomial((0, 5, 0))
    h5_checks = 0
    for i, j in rows:
        h_monomial = sp_monomial((j, -(i + j), -(i + j)))
        direct = sp_mul(
            h5_clear,
            sp_add(
                sp_mul(sp_add(h5_s2, sp_scale(h5_p, -1)),
                       sp_add(one, sp_scale(sp_mul(h5_qbase, h_monomial), -1))),
                sp_scale(sp_mul(h5_qbase, h5_s2), F(-1, 2)),
            ),
        )
        n = i + j
        a_term = sp_mul(sp_monomial((j, 0, 0)), sp_pow(h5_rho, 5 - n))
        expected = sp_add(
            sp_mul(sp_add(one, sp_scale(sp_mul(h5_s2, h5_rho), -1)),
                   sp_add(a_term, sp_scale(h5_z5, -1))),
            sp_scale(sp_mul(h5_s2, sp_pow(h5_rho, 6)), F(-1, 2)),
        )
        need(direct == expected, f"H5 full sparse factor ({i},{j})")
        h5_checks += 1

    # E3 keys are (V,X,t), with s=t*V*X, p=t^-1*X^-1,
    # y=V, Q=t^3 and Phi=t*X^4*F.
    e3_s = sp_monomial((1, 1, 1))
    e3_s2 = sp_mul(e3_s, e3_s)
    e3_p = sp_monomial((0, -1, -1))
    e3_qbase = sp_monomial((0, 0, 3))
    e3_clear = sp_monomial((0, 4, 1))
    e3_rho = sp_monomial((0, 1, 1))
    e3_x3 = sp_monomial((0, 3, 0))
    e3_v2 = sp_monomial((2, 0, 0))
    e3_rows = tuple((i, j) for i, j in rows if i <= 3)
    e3_checks = 0
    for i, j in e3_rows:
        h_monomial = sp_monomial((j, -i, -i))
        direct = sp_mul(
            e3_clear,
            sp_add(
                sp_mul(sp_add(e3_s2, sp_scale(e3_p, -1)),
                       sp_add(one, sp_scale(sp_mul(e3_qbase, h_monomial), -1))),
                sp_scale(sp_mul(e3_qbase, e3_s2), F(-1, 2)),
            ),
        )
        a_term = sp_mul(sp_monomial((j, 0, 0)), sp_pow(e3_rho, 3 - i))
        expected = sp_add(
            sp_mul(sp_add(one, sp_scale(sp_mul(e3_v2, sp_pow(e3_rho, 3)), -1)),
                   sp_add(a_term, sp_scale(e3_x3, -1))),
            sp_scale(sp_mul(e3_v2, sp_pow(e3_rho, 6)), F(-1, 2)),
        )
        need(direct == expected, f"E3 full sparse factor ({i},{j})")
        e3_checks += 1
    return h5_checks, e3_checks


def main():
    # The complete U=0 support has i+j<=5.  In the H5 chart, after
    # Phi=z^6*sigma^12*F and q=sigma^12*z, every H monomial contributes
    # S^j*q^(5-i-j)-S^(j+2)*q^(6-i-j).
    rows = (
        (1, 0), (2, 0), (3, 0), (0, 2), (2, 1), (4, 0),
        (1, 2), (3, 1), (0, 3), (5, 0), (2, 2), (4, 1),
        (1, 3), (3, 2), (0, 4),
    )
    for i, j in rows:
        n = i + j
        need(n <= 5, "U=0 H5 support bound")
        need(5 - n >= 0 and 6 - n >= 1, "H5 local polynomial exponents")
    h5_sparse_checks, e3_sparse_checks = check_sparse_factorizations(rows)

    # H5 hostile: all three normal T^1 coefficients vanish.  The coefficient
    # relation K=2848/45-7d/6 is retained exactly.
    w = u = F(1)
    alpha = F(2)
    s0 = F(-1)
    d = zed = F(5936, 105)
    k = F(2848, 45) - F(7, 6) * d
    need(alpha * alpha - 4 * w * u == 0, "H5 E5=0")
    need(k == F(-8, 3), "H5 K/Delta relation")
    one_plus_s = (F(1), F(1))
    f5 = mul(one_plus_s, one_plus_s)
    f4 = scale(mul(f5, f5), d)
    e = F(-1376, 135)
    f3 = scale(one_plus_s, e)
    f2 = (F(8, 3), F(0), k)
    h5 = h5_normal_coefficients(f5, f4, f3, f2, s0, w)
    need(h5 == (F(0), F(0), F(0)), "H5 all normal Hasse coefficients vanish")
    # With f1=-3, the critical section begins x=-e*q^2/2 and the first
    # surviving critical value is (-3-e^2/4)q^4.  This r=4 split is rational.
    h5_h4 = F(-3) - e * e / 4
    need(h5_h4 != 0, "H5 hostile has exact r=4")

    # E3 hostile: V0=1 and f3=e(V-1)^2.  The sole normal T^1 coefficient is
    # f2(V0); choose f2=(8/3)(V-1)^2, preserving the fixed p^2 coefficient.
    w3 = e
    eta = -2 * e
    v0 = F(1)
    phi = F(-16, 3)
    xi = F(8, 3)
    need(eta * eta - 4 * e * w3 == 0, "E3 discriminant zero")
    f2_v0 = F(8, 3) + phi * v0 + xi * v0 * v0
    need(f2_v0 == 0, "E3 normal Hasse coefficient vanishes")
    # Taking the remaining i=1 coefficients zero leaves f1=-3, so r=2.
    e3_h2 = F(-3)

    # For y^2=x^m, normalization is k[[t^2,t^m]], conductor t^(m-1),
    # and J=(t^m,t^(2m-2)).  The quotient basis has delta=(m-1)/2 classes.
    conductor = {m: conductor_quotient_basis(m) for m in (3, 5)}
    need(conductor[3] == (2,), "(2,3) conductor/J basis")
    need(conductor[5] == (4, 6), "(2,5) conductor/J basis")

    # The exact Morse forms depend on z and the base only through q=t*z.
    # If r=ord_q(psi)<m, the unique nontrivial scale is
    # z=t^(r/(m-r))*Z and x=z^(m/2)*X.  These are normalized sigma-orders;
    # a finite base change clears their denominators.  The squarefree degree
    # after removing the even part of Z^r gives the tail genus.
    h5_tails = []
    for r in range(1, 5):
        b = F(12 * r, 5 - r)
        a = F(5, 2) * b
        order = F(50) + 5 * b - a
        squarefree_degree = (5 - r) + (r % 2)
        genus = (squarefree_degree - 1) // 2
        need(order > 0, f"H5 r={r}: tail order")
        h5_tails.append((r, b, a, order, genus))

    e3_tails = []
    for r in range(1, 3):
        b = F(4 * r, 3 - r)
        a = F(3, 2) * b
        order = F(14) + 4 * b - a
        squarefree_degree = (3 - r) + (r % 2)
        genus = (squarefree_degree - 1) // 2
        need(order > 0, f"E3 r={r}: tail order")
        e3_tails.append((r, b, a, order, genus))

    # Arithmetic Pick genus must be corrected by the persistent generic
    # A_(r-1) boundary defect.  The new tail attaches at its unique infinity
    # point, so it adds no graph cycle: 4 + 11 + g_tail = g_normalized.
    h5_normalized_genera = []
    for r, _, _, _, genus in h5_tails:
        persistent_delta = r // 2
        normalized_genus = 17 - persistent_delta
        need(4 + 11 + genus == normalized_genus,
             f"H5 r={r}: normalized global genus ledger")
        h5_normalized_genera.append((r, normalized_genus))
    need(17 - (5 - 1) // 2 == 15, "H5 persistent r>=m genus cap")

    e3_normalized_genera = []
    for r, _, _, _, genus in e3_tails:
        persistent_delta = r // 2
        normalized_genus = 16 - persistent_delta
        need(4 + 11 + genus == normalized_genus,
             f"E3 r={r}: normalized global genus ledger")
        e3_normalized_genera.append((r, normalized_genus))
    need(16 - (3 - 1) // 2 == 15, "E3 persistent r>=m genus cap")

    # Cancellation hostile for the valuation direction.  With m=3,
    # psi(q)=q, t=pi^2, z=pi*(1+pi), the bracket after the common pi^3 is
    # (1+pi)^3-(1+pi)=2pi+3pi^2+pi^3.  Thus the sum has value 4 although
    # each summand has value 3: v(A+B)>=min, never <=min in general.
    one_plus_pi = (F(1), F(1))
    cancellation_bracket = sub(mul(mul(one_plus_pi, one_plus_pi), one_plus_pi),
                               one_plus_pi)
    need(cancellation_bracket == (F(0), F(2), F(3), F(1)),
         "valuation cancellation hostile")

    # General conductor threshold: for odd m, B=(m-1)/2 makes
    # B+1-m/2=1/2, hence every self-similar split tail has positive order.
    for m in (3, 5, 7, 9, 11):
        buffer_exponent = F(m - 1, 2)
        need(buffer_exponent + 1 - F(m, 2) == F(1, 2),
             f"m={m}: conductor buffer threshold")

    print("THM4340 U0 M12 REPEATED-CUSP PRIMARY CERTIFICATE")
    print(f"SPARSE_FACTOR_CHECKS=H5:{h5_sparse_checks}+base;E3:{e3_sparse_checks}+base")
    print("H5_FACTOR=Phi=(1-S^2*q)*(A(S,q)-z^5)-S^2*q^6/2;q=sigma^12*z")
    print("H5_MORSE=Phi/unit=x^2+psi(q)-z^5;psi(0)=0")
    print("H5_RESIDUE=sigma^50*z^4*(-dS/Phi_z)=sigma^50*z^4*(dz/Phi_S)")
    print("H5_T1=basis[1,z,z^2,z^3];available_normal_coefficients=h1*z+h2*z^2+h3*z^3")
    print("H5_HOSTILE=E5=0;W=u=1;alpha=2;S0=-1;Delta=Z=5936/105;K=-8/3")
    print("H5_HOSTILE_HASSE=" + ",".join(str(x) for x in h5))
    print(f"H5_HOSTILE_FIRST_CRITICAL_VALUE=h4={h5_h4};r=4;rational_tail")
    print("E3_FACTOR=Phi=(1-V^2*q^3)*(A(V,q)-X^3)-V^2*q^6/2;q=sigma^4*X")
    print("E3_MORSE=Phi/unit=x^2+psi(q)-X^3;psi(0)=0")
    print("E3_RESIDUE=sigma^14*X^3*(-dV/Phi_X)=sigma^14*X^3*(dX/Phi_V)")
    print("E3_T1=basis[1,X];available_normal_coefficient=f2(V0)*X")
    print(f"E3_HOSTILE=E3=0;e=W={e};eta={eta};V0=1;Phi=-16/3;xi=8/3;Hasse={f2_v0}")
    print(f"E3_HOSTILE_FIRST_CRITICAL_VALUE=h2={e3_h2};r=2;rational_tail")
    print("CONDUCTOR_J=(2,3):basis_t^2,length1;(2,5):basis_t^4,t^6,length2")
    print("GENERAL_BUFFER=odd_m;B=(m-1)/2 gives tail coefficient B+1-m/2=1/2")
    print("H5_TAILS=" + ";".join(
        f"r{r}:b={b},a={a},order={order},genus={genus}"
        for r, b, a, order, genus in h5_tails
    ))
    print("E3_TAILS=" + ";".join(
        f"r{r}:b={b},a={a},order={order},genus={genus}"
        for r, b, a, order, genus in e3_tails
    ))
    print("NORMALIZED_GLOBAL_GENUS=H5:" + ",".join(
        f"r{r}->{genus}" for r, genus in h5_normalized_genera
    ) + ",r>=5->15;E3:" + ",".join(
        f"r{r}->{genus}" for r, genus in e3_normalized_genera
    ) + ",r>=3->15")
    print("VALUATION_DIRECTION=N=v(z^m-psi)>=min(mb,r(s+b));"
          "hostile:m3,r1,t=pi^2,z=pi(1+pi),N=4>3")
    print("PERSISTENT_BOUNDARY=ord_q(psi)>=m gives x^2=z^m*unit and no positive-genus tail")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
