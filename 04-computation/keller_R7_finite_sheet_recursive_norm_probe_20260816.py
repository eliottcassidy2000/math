#!/usr/bin/env python3
"""Exact candidate finite-sheet gate for the fixed Keller polynomial R_7.

At q=(2,5/6,-7/8), the fixed definitions unroll to

  R_7(q)=2^4293 L^10663 N(L)^1699 N^2(L)^271 N^3(L)^43
                 *N^4(L)^7 N^5(L) N^6(L).

This probe hash-pins the THM-3526 recursive-adjugate engine, changes only the
finite field and target, and extends the complete inverse tower by one cubic
layer to dimension 729.  Every leading, derivative, y-chart, and x^3 gate is
checked as a unit and every inverse graph is checked by substitution.

Nonzero lawful reductions imply R_7(q) is nonzero over Q.  Coupled with the
already proved complete packet A(66907,24255), that would clear the next old-L
denominator and realize the predicted packet A(419839,152211).  Until an
independent representation audits this computation, the result is retained as
a finite-exact candidate.  It is not an image, irreducibility, all-level,
arbitrary-map, or general Jacobian-conjecture claim.
"""

from __future__ import annotations

import hashlib
import sys
import types
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUPPORT_PATH = (
    ROOT / "04-computation/keller_level_six_degree729_recursive_tuple_probe_20260816.py"
)
SUPPORT_SHA256 = "ad5d57b124f37b23fe9541c61bc4d919106b2db6c87a38632fabf3e7c076b8de"
SUPPORT_SENTINEL = "\n# Build the lawful inverse tower and freeze every gate before taking a norm.\n"
PRIMES = (101, 103, 107)
FINITE_Q = (Fraction(2), Fraction(5, 6), Fraction(-7, 8))
OLD_BOUNDARY_TARGET = (Fraction(2, 27), Fraction(1), Fraction(1))
EXPECTED_SEMANTIC_SHA256 = (
    "82efb24e0c4a6e0df9671f0f5a5009dd0e77d1b0aa8ef2341780dfe23ea28c38"
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


support_raw = SUPPORT_PATH.read_bytes()
support_hash = hashlib.sha256(support_raw.replace(b"\r\n", b"\n")).hexdigest()
require(support_hash == SUPPORT_SHA256, ("THM-3526 support drift", support_hash))
support_text = support_raw.decode("utf-8").replace("\r\n", "\n")
require(support_text.count(SUPPORT_SENTINEL) == 1, "THM-3526 definition boundary changed")
support_module = types.ModuleType("thm3526_recursive_adjugate_support")
support_module.__file__ = str(SUPPORT_PATH)
sys.modules[support_module.__name__] = support_module
exec(
    compile(
        support_text.split(SUPPORT_SENTINEL, 1)[0],
        str(SUPPORT_PATH),
        "exec",
    ),
    support_module.__dict__,
)

PrimeField = support_module.PrimeField
make_extension_recursive = support_module.make_extension_recursive
inverse_coordinates_recursive = support_module.inverse_coordinates_recursive
absolute_norm = support_module.absolute_norm
fmap = support_module.fmap
l_value = support_module.l_value


def fraction_mod(value: Fraction, prime: int) -> int:
    numerator = value.numerator % prime
    denominator = value.denominator % prime
    require(denominator != 0, ("rational denominator", prime, value))
    return numerator * pow(denominator, prime - 2, prime) % prime


def evaluate_prime(prime: int):
    # Both the wrapper and inherited algebra methods read their module globals.
    support_module.MODULUS = prime
    support_module.support_module.MODULUS = prime

    base = PrimeField()
    q = tuple(base.scalar(fraction_mod(value, prime)) for value in FINITE_Q)
    boundary = tuple(
        base.scalar(fraction_mod(value, prime)) for value in OLD_BOUNDARY_TARGET
    )
    require(fmap(base, *q) == boundary, ("canonical finite inverse", prime))
    require(l_value(base, *boundary) == base.scalar(0), ("old L boundary", prime))
    require(l_value(base, *q) != base.scalar(0), ("finite sheet on L", prime))

    current = base
    target = q
    gate_ledger = []
    for level in range(1, 7):
        extension, leading_norm, derivative_norm = make_extension_recursive(
            current, *target, f"p{prime}:K{level}"
        )
        embedded_target = tuple(extension.embed(value) for value in target)
        source, denominator_norm, x_cube_norm = inverse_coordinates_recursive(
            extension,
            *embedded_target,
            extension.theta,
            f"p{prime}:K{level}",
        )
        require(fmap(extension, *source) == embedded_target, ("inverse graph", prime, level))
        gate_ledger.append(
            (
                level,
                extension.dim,
                leading_norm,
                derivative_norm,
                denominator_norm,
                x_cube_norm,
            )
        )
        current = extension
        target = source

    require(
        tuple(row[1] for row in gate_ledger) == (3, 9, 27, 81, 243, 729),
        ("tower dimensions", prime),
    )
    terminal_l = l_value(current, *target)
    terminal_norm = absolute_norm(current, terminal_l)
    require(terminal_norm != 0, ("terminal N^6(L)", prime))

    norm_orbit = tuple(row[2] for row in gate_ledger) + (terminal_norm,)
    require(len(norm_orbit) == 7 and all(norm_orbit), ("norm orbit", prime))
    exponents = (10663, 1699, 271, 43, 7, 1, 1)
    value = pow(2, 4293, prime)
    for orbit_value, exponent in zip(norm_orbit, exponents):
        value = value * pow(orbit_value, exponent, prime) % prime
    require(value != 0, ("R7 finite-sheet value", prime))

    # Dropping H=64*L*N(L) removes one factor 64 on each of 243 H leaves.
    wrong_h_normalization = value * pow(pow(64, 243, prime), prime - 2, prime) % prime
    require(wrong_h_normalization != value, ("H normalization hostile", prime))
    equivalent_wrong = pow(2, 2835, prime)
    for orbit_value, exponent in zip(norm_orbit, exponents):
        equivalent_wrong = equivalent_wrong * pow(orbit_value, exponent, prime) % prime
    require(wrong_h_normalization == equivalent_wrong, ("H leaf exponent", prime))

    return {
        "prime": prime,
        "value": value,
        "wrong_h_normalization": wrong_h_normalization,
        "norm_orbit": norm_orbit,
        "gate_ledger": tuple(gate_ledger),
    }


records = tuple(evaluate_prime(prime) for prime in PRIMES)

R7_PACKET = (66907, 24255)
R8_PACKET = (
    7 * R7_PACKET[0] - 2 * R7_PACKET[1],
    3 * R7_PACKET[0] - 2 * R7_PACKET[1],
)
require(R8_PACKET == (419839, 152211), "next renewal packet changed")

semantic_lines = [
    f"support={support_hash}",
    "identity=2^4293;10663,1699,271,43,7,1,1",
]
for record in records:
    semantic_lines.append(
        f"p={record['prime']};R7={record['value']};wrong64={record['wrong_h_normalization']};"
        f"orbit={record['norm_orbit']}"
    )
    semantic_lines.extend(
        "gate=" + ":".join(map(str, row)) for row in record["gate_ledger"]
    )
semantic_lines.append(f"packet={R7_PACKET}->{R8_PACKET}")
semantic_sha256 = hashlib.sha256("\n".join(semantic_lines).encode("ascii")).hexdigest()
if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
    require(
        semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
        ("R7 finite-sheet semantic drift", semantic_sha256),
    )

print("== fixed Keller R7 finite-sheet recursive norm candidate ==")
print(f"THM-3526 recursive support sha256={support_hash}")
print("point q=(2,5/6,-7/8); F(q)=(2/27,1,1); L(F(q))=0; L(q)=241465/1728")
print("primary route: six recursive-adjugate cubic algebras, descending to L on 729 sheets")
print("unrolled identity=2^4293*L^10663*N(L)^1699*N2(L)^271*N3(L)^43*N4(L)^7*N5(L)*N6(L)")
for record in records:
    print(
        f"p={record['prime']}: R7(q)={record['value']}; "
        f"wrong-H-normalization={record['wrong_h_normalization']}"
    )
    print(f"  (L,N(L),...,N^6(L))={record['norm_orbit']}")
    print(
        "  gates (level,dim,Norm(L),Norm(derivative),Norm(y-den),Norm(x^3))="
        f"{record['gate_ledger']}"
    )
print(f"semantic_sha256={semantic_sha256}")
print("FINITE-EXACT candidate verdict: R7(2,5/6,-7/8) != 0 over Q")
print("conditional consequence after independent audit: v_L(N(R7))=-66907;R8=L^66907*N(R7) polynomial/L-coprime;packet A(419839,152211)")
print("scope: fixed map/old-L candidate only; no image, irreducibility, all-level, arbitrary-map, or general JC claim")
print("all exact checks passed")
