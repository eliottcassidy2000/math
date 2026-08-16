#!/usr/bin/env python3
"""Exact three-face probe for the fixed Keller norm tower.

For a polynomial P(x,y,z), the single face max(i-k) controls the two
inverse sheets which diverge over L=0, but it does *not* determine the
top face of L^e Norm(P).  The latter also uses the opposite faces

    min(i-k),                 min(i-j-2k).

This companion reconstructs THM-3495's 66,146-term J and checks the full
three-face packet for L, H, and J.  It then evaluates the resulting face
transform exactly.  The first two transforms are calibrated against

    L Norm(L) = H/64,
    L^7 Norm(H) = J/2^35.

The third transform hostilely tests, refutes, and repairs the predicted
weight/factor pair (259,87): the exact pair is (271,99) for
G=L^43 Norm(J), without constructing G.  A separate finite-field cubic
norm at the hostile finite inverse point tests whether that sheet can
cancel the next pole.  ``require`` remains active under ``python -O``.

Scope: one fixed map and one norm orbit.  This is not an all-polynomial
recurrence, an irreducibility result for G, or a Jacobian-conjecture claim.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import math
import pickle
import runpy
from fractions import Fraction
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def digest_fraction(value: Fraction) -> str:
    payload = f"{value.numerator}/{value.denominator}".encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def nonzero_scalar_multiple(actual, model, variables, label: str):
    """Return c when actual=c*model in Q[variables], and fail otherwise."""

    actual_poly = sp.Poly(actual, *variables, domain=sp.QQ)
    model_poly = sp.Poly(model, *variables, domain=sp.QQ)
    require(actual_poly.monoms() == model_poly.monoms(), f"{label}: support changed")
    ratios = {
        actual_poly.coeff_monomial(monomial) / model_poly.coeff_monomial(monomial)
        for monomial in model_poly.monoms()
    }
    require(len(ratios) == 1, f"{label}: coefficients are not a scalar multiple")
    scalar = next(iter(ratios))
    require(scalar != 0, f"{label}: scalar vanished")
    return scalar


ROOT = Path(__file__).resolve().parents[1]
GLOBAL_PROBE = ROOT / "04-computation/keller_level_three_global_norm_probe_20260816.py"
H_PATH = ROOT / "05-knowledge/results/keller_L2_polynomial_opus_20260728.pkl"
J_LEDGER_SHA256 = "9aca78e67d33351b2f2fb4dbe8ab5bdff06373fdbd8ef9ec73d29b15bffedefe"

# Reconstruct J through the frozen exact route.  Suppressing the inherited
# report keeps this companion's output deterministic and focused.
captured = io.StringIO()
with contextlib.redirect_stdout(captured):
    namespace = runpy.run_path(str(GLOBAL_PROBE))
require(
    captured.getvalue().rstrip().endswith("all exact checks passed"),
    "the imported global J reconstruction did not finish cleanly",
)
J = namespace["J"]
require(namespace["coefficient_hash"](J) == J_LEDGER_SHA256, "global J ledger changed")
sp = namespace["sp"]

H = pickle.loads(H_PATH.read_bytes())
H_poly = namespace["sp"].Poly(H)

L_terms = {
    (2, 0, 2): 27,
    (1, 1, 1): -18,
    (1, 0, 0): 16,
    (0, 3, 1): 1,
    (0, 2, 0): -1,
}
H_terms = {tuple(map(int, monomial)): int(coefficient) for monomial, coefficient in H_poly.terms()}
J_terms = {
    tuple(map(int, monomial)): int(coefficient)
    for monomial, coefficient in J.to_dict().items()
}


def binomial_top_face(e: int, m: int, coefficient: int) -> dict[tuple[int, int, int], int]:
    """Coefficient dictionary for coefficient*x^e*(3*x*z-2*y)^m."""

    return {
        (e + r, m - r, r): coefficient
        * math.comb(m, r)
        * 3**r
        * (-2) ** (m - r)
        for r in range(m + 1)
    }


def gaussian_small_value(
    face: dict[tuple[int, int, int], int]
) -> tuple[Fraction, Fraction, Fraction]:
    """Evaluate the small-edge face at A=1 and W^2=-1/4.

    The two small inverse roots have leading coordinates

        (x,y,z)=(W,6W,-26),    W^2=-1/4.

    Return p,q,norm for value p+qW and its conjugate product
    p^2+q^2/4.
    """

    p = Fraction(0)
    q = Fraction(0)
    for (i, j, k), coefficient in face.items():
        power = i + j
        scalar = Fraction(coefficient * 6**j * (-26) ** k)
        half = power // 2
        scalar *= Fraction((-1) ** half, 4**half)
        if power % 2:
            q += scalar
        else:
            p += scalar
    norm = p * p + q * q / 4
    return p, q, norm


def packet(
    name: str,
    terms: dict[tuple[int, int, int], int],
    expected_e: int,
    expected_m: int,
    expected_top_coefficient: int,
) -> dict[str, object]:
    require(expected_m % 3 == 0, f"{name}: packet requires m divisible by 3")
    require(0 <= expected_m <= expected_e, f"{name}: packet exponent range failed")
    lambda_values = {monomial: monomial[0] - monomial[2] for monomial in terms}
    beta_values = {
        monomial: monomial[0] - monomial[1] - 2 * monomial[2]
        for monomial in terms
    }
    gamma_values = {
        monomial: monomial[0] - monomial[1] - 5 * monomial[2]
        for monomial in terms
    }
    lambda_max = max(lambda_values.values())
    lambda_min = min(lambda_values.values())
    beta_min = min(beta_values.values())
    gamma_min = min(gamma_values.values())
    top = {m: c for m, c in terms.items() if lambda_values[m] == lambda_max}
    low = {m: c for m, c in terms.items() if lambda_values[m] == lambda_min}
    small = {m: c for m, c in terms.items() if beta_values[m] == beta_min}
    gamma = {m: c for m, c in terms.items() if gamma_values[m] == gamma_min}
    z_degree = max(monomial[2] for monomial in terms)
    z_top = {m: c for m, c in terms.items() if m[2] == z_degree}

    require(lambda_max == expected_e, f"{name}: top lambda weight changed")
    require(
        top == binomial_top_face(expected_e, expected_m, expected_top_coefficient),
        f"{name}: pinned top face changed",
    )
    require(lambda_min == -expected_e, f"{name}: opposite lambda weight is not -e")
    require(len(low) == 1, f"{name}: opposite lambda face is not monomial")
    expected_low_monomial = next(iter(low))
    renewal_m = 3 * expected_e - 2 * expected_m
    require(
        expected_low_monomial == (0, renewal_m, expected_e),
        f"{name}: opposite lambda face is not y^r*z^e; term={expected_low_monomial}",
    )
    small_p, small_q, small_norm = gaussian_small_value(small)
    require(small_norm != 0, f"{name}: small conjugate resultant vanished")
    xx, yy, zz = sp.symbols("xx yy zz")
    small_factor = sp.factor(
        sum(coefficient * yy**j * zz**k for (i, j, k), coefficient in small.items())
    )
    z_top_factor = sp.factor(
        sum(coefficient * xx**i * yy**j for (i, j, k), coefficient in z_top.items())
    )
    gamma_factor = sp.factor(
        sum(coefficient * xx**i * yy**j * zz**k for (i, j, k), coefficient in gamma.items())
    )
    beta_model = (
        yy**renewal_m
        * zz ** (expected_e - expected_m)
        * (yy**2 + 27 * zz) ** (2 * expected_m // 3)
        * (yy**2 + 108 * zz) ** (expected_m // 3)
    )
    z_top_model = (
        xx ** (2 * expected_e - 4 * expected_m // 3)
        * zz ** (2 * expected_e - 2 * expected_m // 3)
    )
    gamma_model = (
        zz**expected_e
        * (27 * xx**2 * zz + yy**3) ** (expected_e - 2 * expected_m // 3)
    )
    beta_scalar = nonzero_scalar_multiple(
        small_factor, beta_model, (xx, yy, zz), f"{name}: beta face"
    )
    z_top_scalar = nonzero_scalar_multiple(
        z_top_factor * zz**z_degree,
        z_top_model,
        (xx, yy, zz),
        f"{name}: z-top face",
    )
    gamma_scalar = nonzero_scalar_multiple(
        gamma_factor, gamma_model, (xx, yy, zz), f"{name}: gamma face"
    )

    require(beta_min == -5 * expected_e + 2 * expected_m, f"{name}: beta weight changed")
    require(
        z_degree == 2 * expected_e - 2 * expected_m // 3,
        f"{name}: z degree changed",
    )
    require(gamma_min == -8 * expected_e + 2 * expected_m, f"{name}: gamma weight changed")

    next_e = expected_e - lambda_min - beta_min
    next_m = expected_low_monomial[1]
    require(next_e == 7 * expected_e - 2 * expected_m, f"{name}: matrix e-step changed")
    require(next_m == renewal_m, f"{name}: matrix m-step changed")
    low_coefficient = low[expected_low_monomial]
    raw_next_coefficient = (
        Fraction(16**expected_e)
        * low_coefficient
        * Fraction((-1) ** next_m, 2**next_m)
        * small_norm
    )
    return {
        "name": name,
        "e": expected_e,
        "m": expected_m,
        "top_terms": len(top),
        "lambda_min": lambda_min,
        "low_coefficient": low_coefficient,
        "beta_min": beta_min,
        "small_terms": len(small),
        "small_p": small_p,
        "small_q": small_q,
        "small_norm": small_norm,
        "small_factor": small_factor,
        "beta_scalar": beta_scalar,
        "z_degree": z_degree,
        "z_top_terms": len(z_top),
        "z_top_factor": z_top_factor,
        "z_top_scalar": z_top_scalar,
        "gamma_min": gamma_min,
        "gamma_terms": len(gamma),
        "gamma_factor": gamma_factor,
        "gamma_scalar": gamma_scalar,
        "next_e": next_e,
        "next_m": next_m,
        "raw_next_coefficient": raw_next_coefficient,
    }


H_TOP_COEFFICIENT = -63078912
J_TOP_COEFFICIENT = -(2**58) * (3**51) * (13**8) * (79**4) * (313**2)

packets = [
    packet("L", L_terms, 1, 0, 16),
    packet("H", H_terms, 7, 3, H_TOP_COEFFICIENT),
    packet("J", J_terms, 43, 15, J_TOP_COEFFICIENT),
]

EXPECTED_FACE_HASHES = {
    "L": (
        "e43868d45dc09e7bcb8e00230b853462cac5fc867b9ab868c6299051082698b0",
        "776084fe3906c1280e8eae27292fbc3967e204a10b0fb0c0712e2235d631cfab",
    ),
    "H": (
        "6dbe7847c8e8b46b01f44688222712e9732f3eec021e9955b32b73be148a75a1",
        "fa77606cba17ed3cb90053bd355c249d73eb3259ea407adbce34636aa4341b81",
    ),
    "J": (
        "14890d5df0a74de980b2c416e60ec4ebb9f45f24bb1074c94aeb28185b6d58ef",
        "fa3ed2569db09620551a553a22c8392e193ba12ec8783ad66456c3af9fe6f49f",
    ),
}
for row in packets:
    expected_small_hash, expected_raw_hash = EXPECTED_FACE_HASHES[row["name"]]
    require(
        digest_fraction(row["small_norm"]) == expected_small_hash,
        f"{row['name']}: quadratic face norm ledger changed",
    )
    require(
        digest_fraction(row["raw_next_coefficient"]) == expected_raw_hash,
        f"{row['name']}: transformed face coefficient ledger changed",
    )

# The first two rows independently calibrate the tropical product against
# the already proved exact norm normalizations.
require(
    packets[0]["raw_next_coefficient"] == Fraction(H_TOP_COEFFICIENT, 64),
    "L-to-H raw face coefficient misses Norm(L)=H/(64L)",
)
require(
    packets[1]["raw_next_coefficient"] == Fraction(J_TOP_COEFFICIENT, 2**35),
    "H-to-J raw face coefficient misses Norm(H)=J/(2^35L^7)",
)
require(
    (packets[0]["next_e"], packets[0]["next_m"]) == (7, 3),
    "L-to-H face pair changed",
)
require(
    (packets[1]["next_e"], packets[1]["next_m"]) == (43, 15),
    "H-to-J face pair changed",
)
require(
    (packets[2]["lambda_min"], packets[2]["beta_min"], packets[2]["next_e"], packets[2]["next_m"])
    == (-43, -185, 271, 99),
    "J-to-G repaired face packet changed",
)


def matrix_step(state: tuple[int, int]) -> tuple[int, int]:
    e_value, m_value = state
    return 7 * e_value - 2 * m_value, 3 * e_value - 2 * m_value


formal_states = [(1, 0)]
for _ in range(4):
    formal_states.append(matrix_step(formal_states[-1]))
require(
    formal_states == [(1, 0), (7, 3), (43, 15), (271, 99), (1699, 615)],
    "matrix orbit changed",
)
for index in range(4):
    left_e, left_m = formal_states[index]
    right_e, right_m = formal_states[index + 1]
    require(
        left_e * right_m - left_m * right_e == 3 * (-8) ** index,
        f"Cassini determinant changed at index {index}",
    )
for index, (e_value, m_value) in enumerate(formal_states[1:], start=1):
    require(e_value % 6 == 1 and m_value % 6 == 3, f"mod-six orbit changed at {index}")
    require(math.gcd(e_value, m_value) == 1, f"primitive ratio changed at {index}")

pythagorean_triples = [
    ((e_value**2 - m_value**2) // 2, e_value * m_value, (e_value**2 + m_value**2) // 2)
    for e_value, m_value in formal_states[1:4]
]
require(
    pythagorean_triples
    == [(20, 21, 29), (812, 645, 1037), (31820, 26829, 41621)],
    "primitive Pythagorean sidecar changed",
)


# --- A finite-sheet hostile for the *next* norm ---------------------------

def inv_mod(value: int, prime: int) -> int:
    return pow(value % prime, prime - 2, prime)


def fraction_mod(value: Fraction, prime: int) -> int:
    return value.numerator % prime * inv_mod(value.denominator, prime) % prime


def alg_add(left, right, prime: int):
    return tuple((left[i] + right[i]) % prime for i in range(3))


def alg_scale(value, scalar: int, prime: int):
    return tuple(scalar * entry % prime for entry in value)


def alg_mul(left, right, prime: int, cubic_p: int, cubic_q: int):
    raw = [0] * 5
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            raw[i + j] = (raw[i + j] + left_value * right_value) % prime
    # w^3=-cubic_p*w-cubic_q for the monic core.
    for degree in (4, 3):
        coefficient = raw[degree] % prime
        raw[degree] = 0
        raw[degree - 2] = (raw[degree - 2] - coefficient * cubic_p) % prime
        raw[degree - 3] = (raw[degree - 3] - coefficient * cubic_q) % prime
    return tuple(raw[:3])


def alg_power_table(base, maximum: int, prime: int, cubic_p: int, cubic_q: int):
    values = [(1, 0, 0)]
    for _ in range(maximum):
        values.append(alg_mul(values[-1], base, prime, cubic_p, cubic_q))
    return values


def determinant3(columns, prime: int) -> int:
    matrix = [[columns[column][row] for column in range(3)] for row in range(3)]
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    ) % prime


def invariants(a: int, b: int, c: int, prime: int) -> dict[str, int]:
    values = {}
    values["L"] = (27 * a * a * c * c - 18 * a * b * c + 16 * a + b**3 * c - b * b) % prime
    values["T"] = (4 - 3 * b * c) % prime
    values["S"] = (27 * a * c * c - 9 * b * c + 8) % prime
    values["K"] = (9 * a * c - b) % prime
    values["M"] = (27 * a * c * c - 9 * b * c + 26) % prime
    values["Y0"] = (81 * a * b * c * c - 72 * a * c - 15 * b * b * c + 16 * b) % prime
    values["A1"] = (27 * a * b * c * c + 54 * a * c - 9 * b * b * c + 2 * b) % prime
    values["A2"] = (
        27 * a * b * b * c * c
        + 18 * a * b * c
        - 48 * a
        - 9 * b**3 * c
        + 10 * b * b
    ) % prime
    values["Z0"] = (-9 * values["A2"] * values["T"] - 4 * values["M"] * values["L"]) % prime
    return values


def evaluate_j_norm_at_target(target: tuple[Fraction, Fraction, Fraction], prime: int):
    aa, bb, cc = (fraction_mod(value, prime) for value in target)
    inv = invariants(aa, bb, cc, prime)
    require(inv["L"] != 0 and inv["S"] != 0, f"p={prime}: inverse chart denominator vanished")
    cubic_p = inv["T"] * inv_mod(inv["L"], prime) % prime
    cubic_q = -2 * cc * inv_mod(inv["L"], prime) % prime
    discriminant = (-4 * cubic_p**3 - 27 * cubic_q**2) % prime
    require(discriminant != 0, f"p={prime}: cubic algebra is not etale")

    w = (0, 1, 0)
    y = tuple(
        coefficient * inv_mod(2 * inv["S"], prime) % prime
        for coefficient in (
            inv["Y0"],
            6 * inv["L"],
            -3 * inv["K"] * inv["L"],
        )
    )
    z = tuple(
        coefficient * inv_mod(8 * inv["S"], prime) % prime
        for coefficient in (
            inv["Z0"],
            6 * inv["L"] * inv["A1"],
            -9 * inv["L"] * inv["A2"],
        )
    )

    x_powers = alg_power_table(w, 86, prime, cubic_p, cubic_q)
    y_powers = alg_power_table(y, 129, prime, cubic_p, cubic_q)
    z_powers = alg_power_table(z, 76, prime, cubic_p, cubic_q)
    value = (0, 0, 0)
    for (i, j, k), coefficient in J_terms.items():
        term = alg_mul(x_powers[i], y_powers[j], prime, cubic_p, cubic_q)
        term = alg_mul(term, z_powers[k], prime, cubic_p, cubic_q)
        value = alg_add(value, alg_scale(term, coefficient, prime), prime)

    basis = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    columns = [alg_mul(value, vector, prime, cubic_p, cubic_q) for vector in basis]
    norm_value = determinant3(columns, prime)
    return norm_value, inv["L"], discriminant


finite_q = (Fraction(2), Fraction(5, 6), Fraction(-7, 8))
finite_target = (Fraction(2, 27), Fraction(1), Fraction(1))


def exact_F(point: tuple[Fraction, Fraction, Fraction]):
    x, y, z = point
    u = 1 + x * y
    return (
        u**3 * z + y**2 * u * (4 + 3 * x * y),
        y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y),
        2 * x - 3 * x**2 * y - x**3 * z,
    )


require(exact_F(finite_q) == finite_target, "finite hostile point no longer maps to L=0 target")

finite_controls = []
for prime in (101, 103, 107):
    norm_j, l_at_q, cubic_discriminant = evaluate_j_norm_at_target(finite_q, prime)
    g_at_q = pow(l_at_q, 43, prime) * norm_j % prime
    # Nonzero reduction at even one good prime proves the rational value is
    # nonzero.  Requiring all three is a reproducibility/hostility check.
    require(norm_j != 0 and g_at_q != 0, f"p={prime}: finite G sheet vanished")
    finite_controls.append((prime, norm_j, l_at_q, g_at_q, cubic_discriminant))

require(
    finite_controls
    == [(101, 39, 16, 27, 56), (103, 89, 12, 98, 91), (107, 51, 38, 27, 14)],
    "finite-sheet control ledger changed",
)


print("== fixed Keller tropical norm-face transform ==")
print(f"J coefficient-ledger sha256={J_LEDGER_SHA256}")
print("inverse-chart tropical split: one linear edge plus one quadratic edge")
print("  large edge: (x,y,z)~((c/2)t, -(3ac-2b)/2, a/t)")
print("  small edge: (x,y,z)~(W*t^(1/2), 6aW*t^(-1/2), -26a/t), 4aW^2+1=0")
for row in packets:
    small_norm = row["small_norm"]
    raw_coefficient = row["raw_next_coefficient"]
    print(
        f"{row['name']} packet: (e,m)=({row['e']},{row['m']}); "
        f"lambda_min={row['lambda_min']}; beta_min={row['beta_min']}; "
        f"low=y^{row['next_m']}*z^{row['e']}; "
        f"small_terms={row['small_terms']}; small_norm_sha256={digest_fraction(small_norm)}"
    )
    print(f"  beta_face_factor={row['small_factor']}")
    print(
        f"  z_top: degree={row['z_degree']}; terms={row['z_top_terms']}; "
        f"coefficient_factor={row['z_top_factor']}"
    )
    print(
        f"  gamma_face: min={row['gamma_min']}; terms={row['gamma_terms']}; "
        f"factor={row['gamma_factor']}"
    )
    print(
        f"  transform -> (e_next,m_next)=({row['next_e']},{row['next_m']}); "
        f"raw_face_coefficient_sha256={digest_fraction(raw_coefficient)}"
    )
print("calibration: L*Norm(L)=H/64 and L^7*Norm(H)=J/2^35: PASS")
g_coefficient = packets[2]["raw_next_coefficient"]
print(
    f"G=L^43*Norm(J): in_lambda(G)=C*x^{packets[2]['next_e']}"
    f"*(3*x*z-2*y)^{packets[2]['next_m']}, C!=0"
)
print(
    f"G raw C: sign={1 if g_coefficient > 0 else -1}; "
    f"numerator_bits={abs(g_coefficient.numerator).bit_length()}; "
    f"denominator={g_coefficient.denominator}; sha256={digest_fraction(g_coefficient)}"
)
for prime, norm_j, l_at_q, g_at_q, cubic_discriminant in finite_controls:
    print(
        f"finite G control mod {prime}: Norm(J)={norm_j}; L={l_at_q}; "
        f"G={g_at_q}; cubic_disc={cubic_discriminant}"
    )
print(
    f"finite sheet G(2,5/6,-7/8) != 0; "
    f"hence v_L(Norm(G))=-{packets[2]['next_e']}"
)
print("hostile verdict: predicted (259,87) is REFUTED; exact next pair is (271,99)")
print("five-face one-step law: (e_next,m_next)=(7e-2m,3e-2m)")
print("conditional iteration: u_(n+2)=5u_(n+1)+8u_n for u=e and u=m")
print("transported G faces: lambda_min=-271 with y^615*z^271; beta_min=-1157")
print("R5=L^271*Norm(G): exposed pair=(1699,615); next finite-sheet gate remains open")
print("Cassini sidecar: det(v_n,v_(n+1))=3*(-8)^n; reduced ratios are not Farey neighbors")
print(f"primitive Pythagorean triples: {pythagorean_triples}")
print("scope: fixed norm orbit only; no irreducibility, all-level induction, or general JC claim")
print("all exact checks passed")
