#!/usr/bin/env python3
"""Finite exact controls for THM-4285.

The theorem itself proves the all-n identity by modular forms.  This audit is
deliberately finite: it compares a direct A2-lattice enumeration and exact
convolution with an independent divisor-sieve evaluation of the claimed
Fourier coefficients.  It also freezes the zero-coordinate strata that
become the order-four/order-two contact strata in THM-4280's V-consumer.
"""

from __future__ import annotations

from hashlib import sha256
from math import isqrt


LIMIT = 4096


def require(condition: bool, message: str) -> None:
    """Optimization-safe audit gate."""
    if not condition:
        raise RuntimeError(message)


def norm(x: int, y: int) -> int:
    """The Eisenstein/A2 norm x^2-xy+y^2."""
    return x * x - x * y + y * y


def digest(values: list[int]) -> str:
    """Stable digest of a labelled integer vector."""
    payload = "\n".join(f"{index},{value}" for index, value in enumerate(values))
    return sha256((payload + "\n").encode("ascii")).hexdigest()


# Since 2Q(x,y)=(x-y)^2+x^2+y^2 >= x^2+y^2, every vector of norm
# at most LIMIT lies in this square.  This makes the direct enumeration
# complete, not heuristic.
coordinate_bound = isqrt(2 * LIMIT) + 1
rho = [0] * (LIMIT + 1)
rho_plus = [0] * (LIMIT + 1)
for x in range(-coordinate_bound, coordinate_bound + 1):
    for y in range(-coordinate_bound, coordinate_bound + 1):
        value = norm(x, y)
        if value <= LIMIT:
            rho[value] += 1
        plus_value = x * x + x * y + y * y
        if plus_value <= LIMIT:
            rho_plus[plus_value] += 1

require(rho == rho_plus, "the y -> -y norm-convention control failed")
require(rho[0] == 1, "the A2 zero vector is not unique")

# Ordered representations by two norms.  For n <= LIMIT, neither summand can
# exceed n, so the truncated norm histogram is complete for this convolution.
representation_count = [0] * (LIMIT + 1)
for n in range(LIMIT + 1):
    representation_count[n] = sum(rho[k] * rho[n - k] for k in range(n + 1))

# Independent divisor-sieve side: sigma[n] is computed without using the
# lattice data.  The n=0 coefficient is a separate boundary case.
sigma = [0] * (LIMIT + 1)
for divisor in range(1, LIMIT + 1):
    for multiple in range(divisor, LIMIT + 1, divisor):
        sigma[multiple] += divisor

formula_count = [1]
for n in range(1, LIMIT + 1):
    correction = 3 * sigma[n // 3] if n % 3 == 0 else 0
    formula_count.append(12 * (sigma[n] - correction))

require(
    representation_count == formula_count,
    "direct A2 convolution disagrees with the divisor formula",
)
require(all(value > 0 for value in representation_count), "universality failed")

# If n=3^a*m with 3 not dividing m, the divisor coefficient collapses to
# 12*sigma(m).  Powers of three are the most cancellation-sensitive rows.
collapsed_count = [1]
for n in range(1, LIMIT + 1):
    prime_to_three = n
    while prime_to_three % 3 == 0:
        prime_to_three //= 3
    collapsed_count.append(12 * sigma[prime_to_three])
require(formula_count == collapsed_count, "3-adic coefficient collapse failed")

powers_three = []
power = 1
while power <= LIMIT:
    powers_three.append(power)
    require(representation_count[power] == 12, "power-of-three hostile changed")
    power *= 3

# In the V-shell dictionary, a=0 gives the order-four stratum and contributes
# exactly rho[n] ordered classes; a!=0 gives the order-two stratum.
order_four_count = rho[:]
order_four_count[0] = 0
order_two_count = [
    representation_count[n] - order_four_count[n] for n in range(LIMIT + 1)
]
order_two_count[0] = 0
require(order_two_count[0] == 0, "the zero shell acquired an order-two row")
for n in range(1, LIMIT + 1):
    if order_four_count[n] == 0:
        require(
            order_two_count[n] == representation_count[n] > 0,
            "a nonnorm shell lacks its required order-two witnesses",
        )

inert_hostiles = [2, 5, 8, 11]
for n in inert_hostiles:
    require(rho[n] == 0, f"expected nonnorm hostile {n} became a norm")
    require(order_two_count[n] > 0, f"nonnorm hostile {n} lost universality")

print("THM4285_TWO_EISENSTEIN_NORM_EXACT_AUDIT_V1")
print("LIMIT", LIMIT, "COORDINATE_BOUND", coordinate_bound)
print("R_ZERO", representation_count[0], "FORMULA_BOUNDARY_SEPARATE")
print("R_1_TO_15", ",".join(map(str, representation_count[1:16])))
print("POWERS_THREE", ",".join(map(str, powers_three)), "COUNT_EACH", 12)
print("INERT_NONNORM_HOSTILES", ",".join(map(str, inert_hostiles)), "PASS")
print("RHO_SHA256", digest(rho))
print("REPRESENTATION_SHA256", digest(representation_count))
print("ORDER_TWO_SHA256", digest(order_two_count))
print("VERDICT PASS FINITE-EXACT; ALL-N PROOF IN THM-4285")
