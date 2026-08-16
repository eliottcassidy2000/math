#!/usr/bin/env python3
"""Exact global reconstruction of the fixed map's level-three norm.

The computation stays in the cubic inverse algebra and never constructs the
degree-27 eliminant or its discriminant.  FLINT integer multivariate
arithmetic reduces H(q(w)) to the cubic basis, takes an eight-term resultant,
and primitively normalizes J=2^35*L^7*Norm(H).  ``require`` remains active
under ``python -O``.
"""

from __future__ import annotations

import hashlib
import math
import pickle
from pathlib import Path

import sympy as sp
from flint import fmpz_mpoly_ctx, fmpz_poly


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def coefficient_hash(poly) -> str:
    """Implementation-independent lexicographic coefficient ledger."""

    ledger = "\n".join(
        f"{monomial}:{int(coefficient)}"
        for monomial, coefficient in sorted(poly.to_dict().items(), reverse=True)
    )
    return hashlib.sha256(ledger.encode("ascii")).hexdigest()


ROOT = Path(__file__).resolve().parents[1]
H_PATH = ROOT / "05-knowledge/results/keller_L2_polynomial_opus_20260728.pkl"
H_RAW_SHA256 = "5a9459b3149e500c1b00b67bd804aa7e607de06bf4610c7cdf5fa26d41d74ce9"
H_LEDGER_SHA256 = "b146c11f33e895c08f72303d282e2668d955e0a58d9268a1b445d4d5202016c2"
EXPECTED_J_LEDGER_SHA256 = "9aca78e67d33351b2f2fb4dbe8ab5bdff06373fdbd8ef9ec73d29b15bffedefe"
EXPECTED_SLICE_HASHES = {
    (1, 2): "de18f1f38b29b92cbbe46c913a1446fe017a4566615d60e37a96d82307ae84a7",
    (3, 1): "f55ce6cbdeadcb27f30a7aa0ca75094b5ded4a344f6d976cf31021b2eb6e6579",
    (1, 3): "c3c98f832b92316343ee24da5f3898c5e3465fc621edd3f063441686cc62febf",
}

raw = H_PATH.read_bytes()
require(hashlib.sha256(raw).hexdigest() == H_RAW_SHA256, "transported H pickle changed")
H = pickle.loads(raw)
sa, sb, sc = sp.symbols("a b c")
H_terms = sp.Poly(H, sa, sb, sc).terms()
H_ledger = "\n".join(
    f"{monomial}:{coefficient}" for monomial, coefficient in H_terms
)
require(
    hashlib.sha256(H_ledger.encode("ascii")).hexdigest() == H_LEDGER_SHA256,
    "transported H coefficient ledger changed",
)

CTX = fmpz_mpoly_ctx.get(["a", "b", "c"], ordering="degrevlex")
a, b, c = CTX.gens()
ZERO, ONE = CTX.constant(0), CTX.constant(1)

L = 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2
T = 4 - 3 * b * c
S = 27 * a * c**2 - 9 * b * c + 8
K = 9 * a * c - b
M = 27 * a * c**2 - 9 * b * c + 26
Y0 = 81 * a * b * c**2 - 72 * a * c - 15 * b**2 * c + 16 * b
A1 = 27 * a * b * c**2 + 54 * a * c - 9 * b**2 * c + 2 * b
A2 = 27 * a * b**2 * c**2 + 18 * a * b * c - 48 * a - 9 * b**3 * c + 10 * b**2
Z0 = -9 * A2 * T - 4 * M * L

# In the quotient by E=L*w^3+T*w-2c, write Y=2S*q_y and Z=8S*q_z.
Y = [Y0, 6 * L, -3 * K * L]
Z = [Z0, 6 * L * A1, -9 * L * A2]


def add(left, right):
    length = max(len(left), len(right))
    result = [ZERO for _ in range(length)]
    for index in range(length):
        if index < len(left):
            result[index] += left[index]
        if index < len(right):
            result[index] += right[index]
    while len(result) > 1 and not result[-1]:
        result.pop()
    return result


def multiply(left, right):
    result = [ZERO for _ in range(len(left) + len(right) - 1)]
    for i, left_coefficient in enumerate(left):
        if not left_coefficient:
            continue
        for j, right_coefficient in enumerate(right):
            if right_coefficient:
                result[i + j] += left_coefficient * right_coefficient
    return result


print("== exact global level-three Keller norm ==", flush=True)
print(f"transported H raw/ledger sha256={H_RAW_SHA256}/{H_LEDGER_SHA256}", flush=True)

Y_powers = [[ONE]]
for _ in range(21):
    Y_powers.append(multiply(Y_powers[-1], Y))
Z_powers = [[ONE]]
for _ in range(12):
    Z_powers.append(multiply(Z_powers[-1], Z))
S_powers = [ONE]
for _ in range(25):
    S_powers.append(S_powers[-1] * S)
L_powers = [ONE]
for _ in range(74):
    L_powers.append(L_powers[-1] * L)
print("inverse-coordinate power tables: PASS", flush=True)

# R=8^25*S^25*H(w,Y/(2S),Z/(8S)).  Total degree(H)=25 makes every
# displayed scalar exponent nonnegative.
R = [ZERO for _ in range(51)]
for (i, j, k), coefficient in H_terms:
    yz = multiply(Y_powers[j], Z_powers[k])
    scalar = int(coefficient) * 2 ** (75 - j - 3 * k) * S_powers[25 - j - k]
    for degree, value in enumerate(yz):
        if value:
            R[i + degree] += scalar * value
require(len(R) == 51 and R[-1], "transported H produced an unexpected w-degree")
print(
    f"scaled H(q): w-degree=50; coefficient terms={sum(len(value) for value in R)}",
    flush=True,
)

# Natural reductions w^d=W_d/L^e_d in the basis 1,w,w^2.
W = [
    ([ONE, ZERO, ZERO], 0),
    ([ZERO, ONE, ZERO], 0),
    ([ZERO, ZERO, ONE], 0),
]
for degree in range(3, 51):
    left, left_exponent = W[degree - 2]
    right, right_exponent = W[degree - 3]
    common = max(left_exponent, right_exponent)
    triple = [
        -T * left[index] * L_powers[common - left_exponent]
        + 2 * c * right[index] * L_powers[common - right_exponent]
        for index in range(3)
    ]
    W.append((triple, common + 1))
common_L_exponent = max(exponent for _triple, exponent in W)
require(common_L_exponent == 24, "cubic basis reduction changed its L exponent")

A = [ZERO, ZERO, ZERO]
for degree, coefficient in enumerate(R):
    if not coefficient:
        continue
    triple, exponent = W[degree]
    scale = L_powers[common_L_exponent - exponent]
    for index in range(3):
        if triple[index]:
            A[index] += coefficient * triple[index] * scale

coefficient_gcd = A[0].gcd(A[1]).gcd(A[2])
expected_gcd = 2**40 * L**21 * S**24
require(coefficient_gcd == expected_gcd, "reduced cubic coefficients have the wrong gcd")
B0, B1, B2 = [coefficient // coefficient_gcd for coefficient in A]
print(
    "H(q(w))=(B0+B1*w+B2*w^2)/(2^35*S*L^3); "
    f"B term counts={(len(B0), len(B1), len(B2))}",
    flush=True,
)

# Res(L*w^3+T*w+C, B2*w^2+B1*w+B0), verified once abstractly here.
sw, sL, sT, sC, sB0, sB1, sB2 = sp.symbols("w L T C B0 B1 B2")
abstract_formula = (
    sB0**3 * sL**2
    - 2 * sB0**2 * sB2 * sL * sT
    + sB0 * sB1**2 * sL * sT
    + 3 * sB0 * sB1 * sB2 * sC * sL
    + sB0 * sB2**2 * sT**2
    - sB1**3 * sC * sL
    - sB1 * sB2**2 * sC * sT
    + sB2**3 * sC**2
)
require(
    sp.expand(
        sp.resultant(sL * sw**3 + sT * sw + sC, sB2 * sw**2 + sB1 * sw + sB0, sw)
        - abstract_formula
    )
    == 0,
    "eight-term resultant formula changed",
)

C = -2 * c
resultant = (
    B0**3 * L**2
    - 2 * B0**2 * B2 * L * T
    + B0 * B1**2 * L * T
    + 3 * B0 * B1 * B2 * C * L
    + B0 * B2**2 * T**2
    - B1**3 * C * L
    - B1 * B2**2 * C * T
    + B2**3 * C**2
)
require(resultant.content() == 2**70, "raw resultant content changed")
J_scaled, remainder = divmod(resultant, L**4 * S**3)
require(not remainder, "resultant is not divisible by L^4*S^3")
content, J = J_scaled.primitive()
require(int(content) == 2**70, "primitive global numerator has the wrong content")
require(
    (len(J), J.degrees(), J.total_degree()) == (66146, (86, 129, 76), 157),
    "primitive global numerator degree/term ledger changed",
)
require(J.gcd(L) == ONE, "global numerator contains L")
H_flint = CTX.from_dict(
    {monomial: int(coefficient) for monomial, coefficient in H_terms}
)
require(J.gcd(H_flint) == ONE, "global numerator contains H")
squarefree_unit, squarefree_factors = J.factor_squarefree()
require(
    squarefree_unit == 1
    and len(squarefree_factors) == 1
    and squarefree_factors[0][0] == J
    and squarefree_factors[0][1] == 1,
    "global numerator is not squarefree",
)
J_digest = coefficient_hash(J)
require(J_digest == EXPECTED_J_LEDGER_SHA256, "global J coefficient ledger changed")

print(
    "resultant=2^70*L^4*S^3*J; hence Norm(H)=J/(2^35*L^7)",
    flush=True,
)
print("J: primitive, squarefree, gcd(J,LH)=1", flush=True)
print(
    f"J: terms={len(J)}; multidegree={J.degrees()}; total degree={J.total_degree()}",
    flush=True,
)
print(f"J coefficient-ledger sha256={J_digest}", flush=True)


def exact_slice(poly, b_value: int, c_value: int):
    coefficients = [0 for _ in range(poly.degrees()[0] + 1)]
    for (i, j, k), coefficient in poly.to_dict().items():
        coefficients[i] += int(coefficient) * b_value**j * c_value**k
    content_value = 0
    for coefficient in coefficients:
        content_value = math.gcd(content_value, abs(coefficient))
    require(content_value > 0, "global J vanished on a control slice")
    primitive_coefficients = [coefficient // content_value for coefficient in coefficients]
    primitive = fmpz_poly(primitive_coefficients)
    ledger = "\n".join(
        f"({degree},):{primitive_coefficients[degree]}"
        for degree in range(len(primitive_coefficients) - 1, -1, -1)
        if primitive_coefficients[degree]
    )
    return content_value, primitive, hashlib.sha256(ledger.encode("ascii")).hexdigest()


expected_slice_content = {(1, 2): 2**14, (3, 1): 1, (1, 3): 1}
for slice_values in ((1, 2), (3, 1), (1, 3)):
    slice_content, primitive_slice, slice_digest = exact_slice(J, *slice_values)
    require(slice_content == expected_slice_content[slice_values], "slice content changed")
    require(primitive_slice.degree() == 86, "slice degree changed")
    require(slice_digest == EXPECTED_SLICE_HASHES[slice_values], "slice ledger changed")
    if slice_values == (1, 2):
        factor_unit, factors = primitive_slice.factor()
        require(
            abs(int(factor_unit)) == 1
            and len(factors) == 1
            and factors[0][0] == primitive_slice * int(factor_unit)
            and factors[0][1] == 1,
            "the (1,2) slice is not irreducible",
        )
    print(
        f"slice (b,c)={slice_values}: content={slice_content}; degree=86; sha256={slice_digest}",
        flush=True,
    )

print("three old slice artifacts recovered coefficient-for-coefficient", flush=True)
print("scope: one fixed three-dimensional Keller map; no general JC conclusion", flush=True)
print("all exact checks passed", flush=True)
