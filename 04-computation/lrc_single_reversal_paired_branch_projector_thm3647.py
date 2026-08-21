#!/usr/bin/env python3
"""Exact companion for THM-3647's sparse branch-orbit projector.

The thirteen THM-3625 operators have point-order-reversal covariance
r -> -r-1.
This companion proves that each one of the seven orbit packets already has
exactly two scalar channels, endpoint and middle, and hence affinely recovers
THM-3636's endpoint projector.  All conclusions remain static over one pinned
finite field.
"""

from __future__ import annotations

import ast
from hashlib import sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMPUTATION = ROOT / "04-computation"
sys.path.insert(0, str(COMPUTATION))

import lrc_arc_reversal_middle_address_quotient_thm3636 as D


A = D.A
M = D.M
MOD = A.MOD
POINTS = A.POINTS
INV2 = pow(2, -1, MOD)
CHECKS = 0

THM3636_FILES = (
    (
        "theorem",
        ROOT / "01-canon/theorems/THM-3636-lrc-arc-reversal-and-middle-address-quotient-intertwiner.md",
        "20afb5600c990f7319a4142dae7eb0ce39012d8c1c58c155fc711c0be0b67f92",
    ),
    (
        "script",
        COMPUTATION / "lrc_arc_reversal_middle_address_quotient_thm3636.py",
        "796a80daae5401c37f280ed97d7c80b40e47c06c23c3d03a844f38f797595001",
    ),
    (
        "output",
        ROOT / "05-knowledge/results/lrc_arc_reversal_middle_address_quotient_thm3636.out",
        "c96d5b7b107c6b67a0500020727f349478908baffcce87ce336642c316ede82a",
    ),
)

# (representative, mate, c0, c1, alpha_on_W, beta_on_M), where
# Pi_W=c0 I+c1 B and B=alpha Pi_W+beta Pi_M.
EXPECTED_RECORDS = (
    (0, 12, 469347427394193908097540, 743694593282835709257775,
     583299261868940099523963, 460700142682906771021811),
    (1, 11, 11866950042139627589926, 480803294042397325270155,
     43106534571494424015806, 107225529790717415656375),
    (2, 10, 89945454758491288419934, 289822667326357863318398,
     353459988660952602113073, 107225529790717415656375),
    (3, 9, 254360568621407416847893, 468406141801133917520805,
     32271302543699805871617, 189845688026952442643144),
    (4, 8, 143267695071213814847927, 434643651478141871468116,
     440916124876065854346257, 745044336390295863842285),
    (5, 7, 637022687995553438312113, 455370956553165074603945,
     327919687265179724083408, 122334475151651243469948),
    (6, 6, 327751551088300035847056, 278406788612936780834335,
     485148529749842658283160, 533745727702934015947346),
)


def require(condition: bool, payload: object) -> None:
    global CHECKS
    if condition is not True:
        raise RuntimeError(payload)
    CHECKS += 1


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def standard(indices: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(int(index in indices) for index in range(POINTS))


def diagonal(values: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(values[row] if row == column else 0 for column in range(POINTS))
        for row in range(POINTS)
    )


def permutation_matrix(permutation: tuple[int, ...]):
    return tuple(
        tuple(int(column == permutation[row]) for column in range(POINTS))
        for row in range(POINTS)
    )


def main() -> None:
    parent_files = D.PARENT_FILES + THM3636_FILES
    parent_hashes = tuple((name, lf_sha256(path))
                          for name, path, _expected in parent_files)
    require(all(observed == expected
                for (_name, _path, expected), (_label, observed)
                in zip(parent_files, parent_hashes)), "parent source drift")

    _zeta, _generators, addresses = A.build_addresses()
    identity = A.identity_matrix()
    zero = A.zero_matrix()
    arc = permutation_matrix((1, 0, 3, 2, 5, 4))
    point_reversal = permutation_matrix((5, 4, 3, 2, 1, 0))
    pi_w = diagonal(standard((0, 1, 4, 5)))
    pi_m = A.matsub(identity, pi_w)
    pi_k = A.scalar(A.matmul(pi_w, A.matadd(identity, arc)), INV2)
    require(D.digest(arc) == D.EXPECTED_DIGESTS["S"], "arc operator")
    require(D.digest(pi_w) == D.EXPECTED_DIGESTS["Pi_W"], "endpoint projector")
    require(D.digest(pi_m) == D.EXPECTED_DIGESTS["Pi_M"], "middle projector")
    require(D.digest(pi_k) == D.EXPECTED_DIGESTS["Pi_K"], "rigidity projector")

    # Check the covariance directly rather than importing the theorem sentence.
    for r in range(13):
        mate = (-r - 1) % 13
        require(A.matmul(A.matmul(point_reversal, addresses[r]), point_reversal)
                == addresses[mate],
                ("branch covariance", r, mate))

    packets = []
    records = []
    for r, mate, c0, c1, alpha, beta in EXPECTED_RECORDS:
        packet = (addresses[r] if r == mate
                  else A.matadd(addresses[r], addresses[mate]))
        packets.append(packet)
        require(A.rank(packet) == POINTS, ("packet invertibility", r))
        require(alpha != beta and alpha != 0 and beta != 0,
                ("distinct nonzero spectra", r))
        spectral = A.matadd(A.scalar(pi_w, alpha), A.scalar(pi_m, beta))
        require(packet == spectral, ("two-channel spectral form", r))
        require(A.rank(A.matsub(packet, A.scalar(identity, alpha))) == 2,
                ("endpoint eigenspace", r))
        require(A.rank(A.matsub(packet, A.scalar(identity, beta))) == 4,
                ("middle eigenspace", r))
        require(A.matmul(A.matsub(packet, A.scalar(identity, alpha)),
                         A.matsub(packet, A.scalar(identity, beta))) == zero,
                ("quadratic polynomial", r))
        recovered = A.matadd(A.scalar(identity, c0), A.scalar(packet, c1))
        require(recovered == pi_w, ("affine endpoint recovery", r))
        require(c1 == pow((alpha - beta) % MOD, -1, MOD)
                and c0 == -beta * c1 % MOD, ("Lagrange coefficients", r))
        recovered_k = A.scalar(A.matmul(recovered, A.matadd(identity, arc)), INV2)
        require(recovered_k == pi_k, ("rigidity recovery", r))
        records.append((r, mate, c0, c1, alpha, beta))

    require(A.rank(tuple(A.flat(packet) for packet in packets), 36) == 2,
            "paired packet span")
    packet_total = zero
    for packet in packets:
        packet_total = A.matadd(packet_total, packet)
    require(packet_total == identity, "branch-orbit partition sum")

    k_rows = (standard((0, 1)), standard((4, 5)))
    require(D.canonical(tuple(pi_k)) == D.canonical(k_rows), "Pi_K image")

    semantic = {
        "field": MOD,
        "orbits": tuple((r, mate) for r, mate, *_rest in records),
        "records": tuple(records),
        "packet_span": 2,
        "packet_ranks": tuple(A.rank(packet) for packet in packets),
        "digests": {
            "packets": digest(tuple(packets)),
            "Pi_W": D.digest(pi_w),
            "Pi_K": D.digest(pi_k),
        },
    }

    source = Path(__file__).resolve()
    source_bytes = source.read_bytes()
    require(b"\r\n" not in source_bytes, "source raw LF")
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source_bytes.decode("utf-8")))),
            "Python assert node present")

    print("== THM-3647 LRC single reversal-paired branch projector ==")
    print(f"field=p:{MOD};point_reversal_branch_involution:r->-r-1;orbits={semantic['orbits']}")
    print(f"parent_sha256_lf={parent_hashes}")
    print("packet_algebra=span_dimension:2;sum:I;ranks:(6,6,6,6,6,6,6)")
    print("spectral_form=B_r=alpha_r*Pi_W+beta_r*Pi_M;alpha_r,beta_r:distinct_nonzero")
    print(f"records={tuple(records)}")
    print("recovery=Pi_W=(B_r-beta_r*I)/(alpha_r-beta_r);Pi_K=Pi_W*(I+S)/2;all_7_orbits")
    print(f"packet_digest={semantic['digests']['packets']}")
    print(f"semantic_sha256={digest(semantic)}")
    print(f"CHECKS={CHECKS};imported_address_checks={A.CHECKS}")
    print("scope=static pinned finite-field branch packet;no chronology/current/characteristic-zero/LRC14")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
