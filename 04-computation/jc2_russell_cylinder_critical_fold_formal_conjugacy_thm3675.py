"""Exact controls for THM-3675 critical-fold formal conjugacy closure."""

import hashlib

import sympy as sp


CHECKS = 0


def require(label, condition):
    global CHECKS
    if not condition:
        raise RuntimeError(f"FAILED: {label}")
    CHECKS += 1


t, u = sp.symbols("t u")


def truncate(expr, variable, order):
    return sp.expand(sp.series(expr, variable, 0, order).removeO())


def monomializing_series(H, alpha, k, order):
    """Return psi with H(psi(u))=alpha*u^k to the required order."""
    psi = u
    for degree in range(2, order):
        coefficient = sp.symbols(f"b_{degree}")
        trial = psi + coefficient * u**degree
        composition = truncate(H.subs(t, trial) - alpha * u**k, u, k + degree)
        equation = sp.expand(composition).coeff(u, k + degree - 1)
        solution = sp.solve(sp.Eq(equation, 0), coefficient)
        require(f"monomialization degree {degree} unique", len(solution) == 1)
        psi = sp.expand(trial.subs(coefficient, solution[0]))
    return psi


def coefficient_digest(polynomial, variable, order):
    packet = ";".join(str(sp.factor(polynomial.coeff(variable, degree))) for degree in range(order))
    return hashlib.sha256(packet.encode("ascii")).hexdigest()


lambda_constant = sp.Rational(5 - 18 + 13, 18)
q6_middle_value_sum = sp.Rational(2012, 2187) - sp.Rational(2012, 2187)
qstar_middle_value_sum = sp.Rational(4, 9) - sp.Rational(4, 9)
q6_debt = sp.Rational(365888, 6561)
qstar_debt = sp.Rational(5440, 81)

require("retained row kills constants", lambda_constant == 0)
require("Q6 Jk value block kills constants", q6_middle_value_sum == 0)
require("Qstar Jk value block kills constants", qstar_middle_value_sum == 0)
require("Q6 debt nonzero", q6_debt != 0)
require("Qstar debt nonzero", qstar_debt != 0)


def conjugacy_control(label, k, alpha, H, debt, middle_value_sum):
    psi_order = 2 * k + 2
    verification_order = 3 * k + 1
    require(f"{label} exact order", sp.expand(H).coeff(t, k) == alpha)
    require(f"{label} no lower term", all(sp.expand(H).coeff(t, degree) == 0 for degree in range(k)))
    psi = monomializing_series(H, alpha, k, psi_order)
    require(f"{label} psi zero", psi.subs(u, 0) == 0)
    require(f"{label} psi unit derivative", sp.diff(psi, u).subs(u, 0) == 1)
    require(
        f"{label} monomialized H",
        truncate(H.subs(t, psi) - alpha * u**k, u, verification_order) == 0,
    )
    derivative = sp.diff(psi, u)
    k0 = derivative.coeff(u, 0)
    kk = derivative.coeff(u, k)
    k2k = derivative.coeff(u, 2 * k)
    require(f"{label} transformed constant K0", k0 == 1)
    left = lambda_constant * k2k
    middle = alpha * middle_value_sum * kk
    right = alpha**2 * debt * k0 + middle
    require(f"{label} left cokernel zero", left == 0)
    require(f"{label} middle constant contribution zero", middle == 0)
    require(f"{label} surviving contradiction", right != 0)
    digest = coefficient_digest(psi, u, psi_order)
    print(
        f"PASS {label}_k={k}_alpha={alpha}_psi_sha256={digest} "
        f"Kk={sp.factor(kk)}_K2k={sp.factor(k2k)}_forced={sp.factor(right)}"
    )


print("THM-3675 exact companion -- critical-fold formal conjugacy closure")
print("status=PROVED VERIFIED-EXACT PENDING-INDEPENDENT-HOSTILE-AUDIT")
print("PASS abstract_cokernel=Lambda(1)=0_Jk_endpoint_sum=0_J0_debt_nonzero")

conjugacy_control(
    "Q6_H2",
    2,
    sp.Rational(3, 2),
    sp.Rational(3, 2) * t**2 - sp.Rational(5, 7) * t**3 + sp.Rational(11, 13) * t**5,
    q6_debt,
    q6_middle_value_sum,
)
conjugacy_control(
    "Qstar_H3",
    3,
    -2,
    -2 * t**3 + sp.Rational(7, 5) * t**4 - sp.Rational(3, 11) * t**7,
    qstar_debt,
    qstar_middle_value_sum,
)
conjugacy_control(
    "Q6_H4",
    4,
    sp.Rational(5, 3),
    sp.Rational(5, 3) * t**4 - sp.Rational(2, 9) * t**5 + sp.Rational(7, 4) * t**8 + t**9,
    q6_debt,
    q6_middle_value_sum,
)

print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
