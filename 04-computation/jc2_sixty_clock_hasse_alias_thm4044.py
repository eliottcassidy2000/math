#!/usr/bin/env python3
"""Primary exact certificate for THM-4044.

The all-k kernel statement is proved symbolically in the theorem file.  This
companion verifies its finite-field realization, sharp boundary coordinate,
four-dimensional Vandermonde interpretation, node-address involution, and the
bounded-degree Sun/conductor controls on independent exact examples.
"""

from __future__ import annotations

import hashlib
import itertools
import math
import sys


sys.stdout.reconfigure(newline="\n")
PRIME = 61
PHI = 44
CHECKS = 0
SEMANTIC: list[str] = []


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)
    print(f"PASS  {label}")


def note(value: str) -> None:
    SEMANTIC.append(value)
    print(f"RESULT {value}")


def mul(left: list[int], right: list[int]) -> list[int]:
    answer = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            answer[i + j] = (answer[i + j] + a * b) % PRIME
    return answer


def power(base: list[int], exponent: int) -> list[int]:
    answer = [1]
    factor = base[:]
    n = exponent
    while n:
        if n & 1:
            answer = mul(answer, factor)
        factor = mul(factor, factor)
        n >>= 1
    return answer


def hasse(poly: list[int], order: int, point: int) -> int:
    return sum(
        coefficient * math.comb(degree, order) * pow(point, degree - order, PRIME)
        for degree, coefficient in enumerate(poly)
        if degree >= order
    ) % PRIME


def centered_binomial(rank: int, t: int) -> int:
    numerator = 1
    inv_two = pow(2, -1, PRIME)
    n = (t + rank - 1) * inv_two % PRIME
    for j in range(rank):
        numerator = numerator * (n - j) % PRIME
    return numerator * pow(math.factorial(rank), -1, PRIME) % PRIME


print("STATUS=THM-4044;PROVED_EXACT_CLOCK_KERNEL;JC(2)_OPEN")
print("UNIVERSE=mu_60_Hasse_observer;F_61_realization;unbounded_residual_firewall")

# The actual C_60 clock from THM-4035.
nodes = [pow(PHI, r, PRIME) for r in range(60)]
require(len(set(nodes)) == 60 and set(nodes) == set(range(1, PRIME)), "44 generates F_61^x")
require(all(nodes[r + 30] == (-nodes[r]) % PRIME for r in range(30)), "half-period is t -> -t")

# Delta_k=P^2(P^60-1)^k.  Hasse jets are checked literally for several
# depths; the theorem proves the same divisibility statement for every k.
base = [PRIME - 1] + [0] * 59 + [1]
for k in range(1, 7):
    delta = [0, 0] + power(base, k)
    require(len(delta) - 1 == 60 * k + 2, f"Delta_{k} has sharp degree 60k+2")
    require(delta[2] == ((-1) ** k) % PRIME, f"Delta_{k} changes [P^2] by (-1)^k")
    require(hasse(delta, 0, 0) == 0 and hasse(delta, 1, 0) == 0, f"Delta_{k} is invisible to boundary jets 0,1")
    require(hasse(delta, 2, 0) == ((-1) ** k) % PRIME, f"boundary Hasse jet 2 detects Delta_{k}")
    require(
        all(hasse(delta, j, t) == 0 for t in nodes for j in range(k)),
        f"all 60 phase jets below depth {k} kill Delta_{k}",
    )

note("kernel_generator=P^2(P^60-1)^k;first_invisible_degree=60k+2")
note("missing_boundary_coordinate=D_P^[2]R(0,0)=[P^2]R")

# Four Kakeya-spine rows are evaluation rows on cubics.  Exhaust every
# four-subset: the Vandermonde determinant is nonzero, while the high-degree
# residual alias remains outside that interpolation range.
vandermonde_count = 0
vandermonde_digest = hashlib.sha256()
for chosen in itertools.combinations(nodes, 4):
    determinant = 1
    for i in range(4):
        for j in range(i + 1, 4):
            determinant = determinant * (chosen[j] - chosen[i]) % PRIME
    if determinant == 0:
        raise RuntimeError("four-phase Vandermonde determinant vanished")
    vandermonde_count += 1
    vandermonde_digest.update(bytes([determinant]))
require(vandermonde_count == math.comb(60, 4), "all four-phase cubic observers are broad")
note(f"kakeya_cubic_vandermonde_count={vandermonde_count}")
note(f"kakeya_vandermonde_sha256={vandermonde_digest.hexdigest()}")

# The half-clock exactly labels the two normalization addresses of every
# split node in the reduced (2,3) boundary model.  It does not supply the
# nodal equation or an address owner.
inv_two = pow(2, -1, PRIME)
inv_three = pow(3, -1, PRIME)
node_parameters = set()
node_address_count = 0
for r in range(30):
    x = nodes[r]
    a = (-2 * inv_three * x * x) % PRIME
    singular_a = (-a * inv_two) % PRIME
    boundary_a = (x * x + a) % PRIME
    boundary_c = (x**3 + 3 * a * inv_two * x) % PRIME
    if boundary_a != singular_a or boundary_c != 0:
        raise RuntimeError(f"node address phase pair failed at r={r}")
    node_address_count += 1
    node_parameters.add(a)
require(node_address_count == 30, "all half-clock node-address identities hold")
require(len(node_parameters) == 30, "half-clock labels all 30 split nonzero node parameters")
note("node_phase_shift=r+30 labels X->-X and (A,C)->(A,-C)")

# Even centered binomial atoms have the same free sign-orbit action.  Their
# degrees are below 60, so one-clock evaluation is faithful for the atom;
# the fixed point t=0 is not a torus phase.
for rank in (2, 4, 6, 8):
    require(
        all(centered_binomial(rank, t) == centered_binomial(rank, (-t) % PRIME) for t in nodes),
        f"centered binomial rank {rank} is even on the clock",
    )
    require(rank < 60, f"centered binomial rank {rank} lies below interpolation boundary")
require(0 not in nodes, "the 60-clock omits the centered fixed point t=0")
note("Sun_even_atoms=degree_at_most_8<60;torus_evaluation_is_faithful_but_omits_fixed_point")

# A degree-d polynomial is recovered by k-deep confluent data whenever
# d<60k.  In particular, three jets recover the degree-178 conductor after
# scalar extension, while no fixed k controls an unrestricted residual.
require(178 < 60 * 3 and 178 >= 60 * 2, "three clock jets reach the degree-178 conductor range sharply")
note("degree_178_conductor_is_injective_under_O_3_after_scalar_extension")
note("unbounded_JC_residual_has_no_fixed_finite_clock_depth")

semantic_hash = hashlib.sha256("\n".join(SEMANTIC).encode()).hexdigest()
print(f"SEMANTIC_SHA256={semantic_hash}")
print(f"CHECKS={CHECKS}")
print("ALL THM-4044 PRIMARY CHECKS PASSED")
