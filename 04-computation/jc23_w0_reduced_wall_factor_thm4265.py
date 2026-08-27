#!/usr/bin/env python3
"""Exact reduced-wall-factor certificate for the THM-4260 frontier.

The program loads the frozen dependency-free polynomial engine and clean-room
class compiler used to audit THM-4260.  Its reciprocal hook records every
reduced coefficient denominator immediately before the already-audited gcd
test.  It then evaluates the denominator and its derivative at t=+1,-1.

Nonzero derivatives modulo the proof prime 397 show that the unavoidable
characteristic-zero wall factors are simple.  The independent prime 577 is a
hostile specialization control.  No W-adic deformation is constructed here.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from pathlib import Path
import subprocess
import types


REVISION = "bafd69bdd607fe4dfa704edc63bef23cbf2e2ce9"
DEPENDENCY = (
    "04-computation/"
    "jc23_w0_canonical_node_attachment_exclusion_independent_audit_thm4260.py"
)
DEPENDENCY_SHA256 = (
    "238af711b734a0e054a3da6d8a836f5f2996e34550d5e507111718d5eff2314f"
)
EMBEDDINGS = ((397, 157, 161, 27), (577, 57, 224, 25))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def load_cleanroom(repo: Path):
    source = subprocess.run(
        ["git", "-C", str(repo), "show", f"{REVISION}:{DEPENDENCY}"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout
    require(sha256(source).hexdigest() == DEPENDENCY_SHA256, "dependency hash changed")
    module = types.ModuleType("jc4265_cleanroom")
    module.__file__ = f"{repo}@{REVISION}:{DEPENDENCY}"
    exec(compile(source, module.__file__, "exec"), module.__dict__)
    module.REPO = str(repo)
    return module


def polynomial_value(polynomial: tuple[int, ...], value: int, prime: int) -> int:
    answer = 0
    for coefficient in reversed(polynomial):
        answer = (answer * value + coefficient) % prime
    return answer


def polynomial_derivative_value(
    polynomial: tuple[int, ...], value: int, prime: int
) -> int:
    answer = 0
    for exponent, coefficient in enumerate(polynomial[1:], start=1):
        answer = (answer + exponent * coefficient * pow(value, exponent - 1, prime)) % prime
    return answer


def audit_embedding(cleanroom, representatives, embedding):
    prime = embedding[0]
    recorded: list[tuple[int, ...]] = []
    original_reciprocal = cleanroom.reciprocal

    def recording_reciprocal(polynomial):
        recorded.append(polynomial)
        return original_reciprocal(polynomial)

    cleanroom.reciprocal = recording_reciprocal
    try:
        group_digests = cleanroom.audit_prime(cleanroom.load_independent(), representatives, embedding)
    finally:
        cleanroom.reciprocal = original_reciprocal

    require(len(recorded) == 280, "reciprocal hook did not see all 280 classes")
    rows = []
    for denominator in recorded:
        value_plus = polynomial_value(denominator, 1, prime)
        value_minus = polynomial_value(denominator, -1, prime)
        derivative_plus = polynomial_derivative_value(denominator, 1, prime)
        derivative_minus = polynomial_derivative_value(denominator, -1, prime)
        require(value_plus == value_minus == 0, "geometric wall root disappeared")
        require(
            derivative_plus != 0 and derivative_minus != 0,
            "geometric wall factor became nonreduced",
        )
        require(derivative_plus == derivative_minus, "odd denominator lost even derivative")
        rows.append((len(denominator) - 1, derivative_plus, derivative_minus))

    profile = Counter(degree for degree, _, _ in rows)
    require(
        profile
        == Counter({11: 8, 19: 36, 23: 24, 27: 64, 35: 52, 39: 72, 43: 24}),
        "denominator-degree profile changed",
    )
    ledger = b"".join(
        value.to_bytes(8, "little") for row in rows for value in row
    )
    return profile, sha256(ledger).hexdigest(), group_digests


def profile_text(profile: Counter[int]) -> str:
    return ",".join(f"{degree}:{profile[degree]}" for degree in sorted(profile))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("."))
    arguments = parser.parse_args()
    repo = arguments.repo.resolve()
    cleanroom = load_cleanroom(repo)
    independent = cleanroom.load_independent()
    representatives = cleanroom.build_representatives(independent)
    require(
        {degree: len(rows) for degree, rows in representatives.items()} == {34: 176, 42: 104},
        "class universe changed",
    )

    print("THM-4265 W0 REDUCED WALL FACTOR")
    print(f"revision={REVISION}")
    print(f"dependency={DEPENDENCY} sha256={DEPENDENCY_SHA256}")
    print("classes=d34:176 d42:104 total:280")
    for embedding in EMBEDDINGS:
        profile, ledger_hash, group_digests = audit_embedding(
            cleanroom, representatives, embedding
        )
        print(
            f"q={embedding[0]} wall_values_zero=280/280 "
            f"wall_derivatives_nonzero=280/280 derivative_even=280/280 "
            f"degree_profile={profile_text(profile)} "
            f"wall_ledger_sha256={ledger_hash} groups={len(group_digests)}"
        )
    print("proof_prime=397 characteristic_zero_wall_factors=t-1,t+1 both_simple")
    print("hostile_prime=577 same_reduced_wall_result")
    print("local_observer=det(d_t,d_W;dstar_t,dstar_W)_at_(W=0,t=+-1)")
    print("scope=no_W_deformation_no_transverse_derivative_no_off_fibre_exclusion")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
