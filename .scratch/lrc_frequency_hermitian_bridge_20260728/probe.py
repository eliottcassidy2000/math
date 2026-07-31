#!/usr/bin/env python3
"""Exact bridge from the THM-2868 split atlas to the THM-2861 edge.

This is a scratch probe, not canonical evidence.  It works projectively:
the common full-current scale cancels, leaving t_r=U_r/V_r.
"""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


RESULTS = ROOT / "05-knowledge/results"
PINNED = {
    COMP / "lrc14_two_origin_endpoint_projective_kummer_thm2868.py":
        "3282526029d27210a5cadaa10f702a974f926e217e6838913d11d9dbc888d9b5",
    RESULTS / "lrc14_two_origin_endpoint_projective_kummer_thm2868.out":
        "ce4b8961b3dfa468d51b884b91002872440b05bf70d7e5eecfc3318d4100e2f9",
    COMP / "lrc14_endpoint_hermitian_edge_holonomy_thm2861.py":
        "57bad76968ec9c61d2202331e007860f2817d15c606d8ba558ab8b8d3c41f20c",
    RESULTS / "lrc14_endpoint_hermitian_edge_holonomy_thm2861.out":
        "9bc846b6269b6ca967d32b5b4091ec506b3ede632c58a249c78211e1ecc8b43d",
}
for path, expected in PINNED.items():
    actual = sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    require(actual == expected, f"pinned dependency changed: {path.name}")

import lrc14_two_origin_endpoint_projective_kummer_thm2868 as atlas


P = 13
PROJECTIVE_EXPONENTS = (
    955, 1501, 2047, 227, 773, 1319, 1865,
    45, 591, 1137, 1683, 2229, 409,
)
FIELD_ROWS = {
    352341050142921841: {
        "u0": 305272472600558724,
        "v0": 26104004983225775,
        "s0": 331376477583784499,
        "t0": 304174682735739267,
    },
    956354278959359281: {
        "u0": 869313819406323966,
        "v0": 482074204885093064,
        "s0": 395033745332057749,
        "t0": 223267090638481188,
    },
}


def dft(values: tuple[int, ...], omega: int, prime: int) -> tuple[int, ...]:
    return tuple(
        sum(
            value * pow(omega, (-a * r) % P, prime)
            for r, value in enumerate(values)
        )
        % prime
        for a in range(P)
    )


def main() -> None:
    summaries = []
    for prime, root in atlas.endpoint.MODS:
        row = FIELD_ROWS[prime]
        xi = pow(root, atlas.endpoint.NN // 2366, prime)
        omega = pow(xi, 182, prime)
        require(
            pow(xi, 2366, prime) == 1
            and pow(xi, 1183, prime) != 1
            and pow(omega, P, prime) == 1
            and omega != 1,
            "root normalization changed",
        )

        t = tuple(pow(xi, exponent, prime)
                  for exponent in PROJECTIVE_EXPONENTS)
        require(
            t[0] == row["t0"]
            and all(t[(r + 1) % P] == omega**3 * t[r] % prime
                    for r in range(P)),
            "projective atlas rotation changed",
        )
        zeta = pow(xi, 2, prime)
        canonical_a = pow(zeta, 624, prime)
        canonical_b = pow(zeta, 510, prime)
        canonical_c = tuple(
            (
                canonical_a
                - canonical_b * pow(omega, 3 * r, prime)
            )
            % prime
            for r in range(P)
        )
        common_source = atlas.COMMON_SOURCE[prime]
        require(
            t[0]
            == -canonical_b * pow(canonical_a, -1, prime) % prime
            and row["v0"] == common_source * canonical_a % prime
            and row["s0"] == common_source * canonical_c[0] % prime
            and all(
                (1 + t[r]) % prime
                == canonical_c[r] * pow(canonical_a, -1, prime) % prime
                for r in range(P)
            ),
            "THM-2868 ratio did not normalize the THM-2861 scalar",
        )

        u = tuple(row["u0"] * pow(omega, 3 * r, prime) % prime
                  for r in range(P))
        v = (row["v0"],) * P
        s = tuple((u_r + v_r) % prime for u_r, v_r in zip(u, v))
        require(
            s[0] == row["s0"]
            and all(value != 0 for value in s)
            and all(
                u_r * pow(v_r, -1, prime) % prime == t_r
                for u_r, v_r, t_r in zip(u, v, t)
            ),
            "recombined signed-current atlas changed",
        )

        # Remove the common V_0 scale.  Then s'_r=1+t_r and conjugation
        # is t_r -> t_r^{-1}.  This avoids assuming that the common source
        # coefficient itself lies in the smaller cyclotomic field.
        normalized_s = tuple((1 + value) % prime for value in t)
        edge = tuple(
            normalized_s[(r + 1) % P]
            * (1 + pow(t[r], -1, prime))
            % prime
            for r in range(P)
        )
        conjugate_edge = tuple(
            (1 + pow(t[(r + 1) % P], -1, prime))
            * normalized_s[r]
            % prime
            for r in range(P)
        )
        require(
            all(edge[r] == pow(omega, 3, prime) * conjugate_edge[r] % prime
                for r in range(P))
            and len(set(edge)) == P
            and all(edge),
            "THM-2861 Hermitian phase/separation law failed",
        )

        transform = dft(edge, omega, prime)
        support = tuple(a for a, value in enumerate(transform) if value)
        require(support == (0, 3, 10),
                "Hermitian edge Fourier support changed")

        # The known oriented phase line makes the coefficient-field
        # symmetric trace algebraically sufficient:
        # E = omega^3/(1+omega^3) (E+bar(E)).
        phase = pow(omega, 3, prime)
        recovery = phase * pow((1 + phase) % prime, -1, prime) % prime
        symmetric_trace = tuple(
            (edge[r] + conjugate_edge[r]) % prime for r in range(P)
        )
        require(
            all(symmetric_trace)
            and all(
                edge[r] == recovery * symmetric_trace[r] % prime
                for r in range(P)
            ),
            "symmetric-trace recovery failed",
        )

        reversed_edge = tuple(conjugate_edge)
        require(
            all(
                reversed_edge[r]
                == pow(omega, -3, prime) * edge[r] % prime
                for r in range(P)
            )
            and tuple(
                (reversed_edge[r] + edge[r]) % prime for r in range(P)
            ) == symmetric_trace,
            "orientation hostile changed",
        )
        require(row["v0"] != 0, "recombined edge scale vanished")
        summaries.append((
            prime,
            support,
            len(set(edge)),
            recovery,
            symmetric_trace[0],
        ))

    print("THM-2868 -> THM-2861 FREQUENCY HERMITIAN BRIDGE")
    print(
        "recombined=S_r=U_r+V_r; "
        "U_(r+1)=omega^3 U_r; V_(r+1)=V_r"
    )
    for prime, support, distinct, recovery, first_trace in summaries:
        print(
            f"field={prime}; edge_support={support}; "
            f"edge_distinct={distinct}; "
            f"E_from_E_plus_barE_factor={recovery}; "
            f"first_symmetric_trace={first_trace}"
        )
    print(
        "identification=1+t_r=c_r/A for the canonical THM-2861 scalar; "
        "V0=P*A and S_r=P*c_r, so the normalized edge is exactly "
        "c_(r+1)*bar(c_r) and the full signed edge differs by P*bar(P)"
    )
    print(
        "positive=Prony summands reconstructed from the same signed q3 "
        "selector give the full three-channel Hermitian edge; its symmetric "
        "trace determines the edge only after the proved forward frequency "
        "orientation is retained"
    )
    print(
        "hostile=edge reversal has the same symmetric trace; "
        "frequency orientation is external data and is not a physical "
        "ancestry orientation"
    )
    print(
        "scope=frequency-dual coefficient bridge only; variable multiplier "
        "charts, positivity, same-triangle ancestry, q11-to-q7 E3 transport, "
        "a physical polarization measurement, and row closure remain absent"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
