#!/usr/bin/env python3
"""Close the two projected k=2 scalar-atlas rows at z1=1732.

The repaired z1736 exact-descent engine is reused without changing its proof
object: periodic ray maxima, denominator states, crude capacities, canonical
common-status infeasibility instances, and (only when necessary) the lossless
projected-cell residual.  The scalar atlas has exactly two rows at this height
and no rows at 1725..1731, so their closure descends the cap to z1 <= 1724.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT_SOURCE = ROOT / "04-computation/lrc14_j7_k2_z1736_exact_descent_closure_thm2941.py"
PARENT_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k2_z1736_exact_descent_closure_thm2941.out"
OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_j7_k2_z1732_exact_descent_closure_thm2941.out"

EXPECTED_PARENT_SOURCE_SHA256 = (
    "1a77ea867330548751e3346cabaa4dd02c5dcb74c10948e9a899c412cd4d0e65"
)
EXPECTED_PARENT_OUTPUT_SHA256 = (
    "746266780aad0bf1b5979940a1d792aace863f76b1d007f3b435fa7f5efc087c"
)
EXPECTED_PROFILE_SHA256 = (
    "516a45f2ce0490f00c05772b1ca2a633e58a2ea71891c3deecf3321f6394e473"
)
EXPECTED_SEMANTIC_SHA256 = (
    "ce4fda72c1607b3172768e0b21ced0586f777c6ee7dd3f94bbd405e6bb0e5099"
)

FIRST = 1732
CASES = (
    (FIRST, (1, 4, 8, 10, 12, 14)),
    (FIRST, (2, 4, 8, 10, 12, 14)),
)
EXPECTED_HEIGHTS = (
    (1683, 1),
    (1694, 10),
    (1702, 3),
    (1708, 14),
    (1722, 11),
    (1724, 2),
    (1732, 2),
    (1736, 15),
)
# scalar, crude kills, common-status kills, status survivors, literal packets.
EXPECTED_COUNTS = {
    CASES[0]: (188, 33, 155, 0, 0),
    CASES[1]: (24, 1, 3, 20, 20),
}
EMPTY_DIGEST = sha256(b"()").hexdigest()
EXPECTED_EXACT_AUDIT = {
    CASES[0]: (
        (
            "c4ce1fde263280dda337184ded196f4f31fae14ed1d0adb197478847d54d7467",
            "b26ca570b73cfa3c4f867e4895cc4078179e232dd5167d34f592e519be79cd85",
            "abc170518ac0c2ee30d1788a0504cb9f24141195267d1c240b97ac82ef34a15a",
            EMPTY_DIGEST,
        ),
        None,
        0,
        None,
        EMPTY_DIGEST,
    ),
    CASES[1]: (
        (
            "14b8157df7fb8300a5a4d41bafa2bd3786d4995e8e9b2a90f4daa05e49d7d70d",
            "678cd0cd9a69539e98a0618175e0354641c600e6d0ec63135f589974c8f5cf94",
            "c33bb1089cfe5bdae95b956d7d2a8fbae51a002361262734054f60ce96780202",
            "be0363f9e947ad7ab2bf768614433e933b2712e02bfb1276b2d42a2dd46ab03e",
        ),
        F(903353, 7131943),
        2,
        F(1),
        "70462ba71eb90368300dc9c17e75a82eb70fca2ee04b1fa405f163bbc8018d02",
    ),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    if value is None:
        return "NONE"
    return f"{value.numerator}/{value.denominator}"


require(
    file_sha256(PARENT_SOURCE) == EXPECTED_PARENT_SOURCE_SHA256,
    "z1736 exact-descent source changed",
)
require(
    file_sha256(PARENT_OUTPUT) == EXPECTED_PARENT_OUTPUT_SHA256,
    "z1736 exact-descent output changed",
)
SPEC = spec_from_file_location("k2_z1732_parent", PARENT_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1736 parent")
Z = module_from_spec(SPEC)
SPEC.loader.exec_module(Z)
require(Z.ATLAS_HEIGHTS == EXPECTED_HEIGHTS, "scalar-atlas height profile changed")
require(tuple(key for key in Z.ATLAS_KEYS if key[0] == FIRST) == CASES, "z1732 universe changed")
require(not any(1724 < first < FIRST for first, _body in Z.ATLAS_KEYS), "descent gap occupied")
require(max(first for first, _body in Z.ATLAS_KEYS if first < FIRST) == 1724, "next height changed")


def profile(case):
    first, body = case
    carrier = Z.E.U.suffix.A.carrier_for(body)
    require(Z.E.P.A.carrier_for(body) == carrier, (case, "carrier engines disagree"))
    h = F(sum(right - left for left, right in carrier), Z.E.U.suffix.A.RULER)
    lower = h * Z.E.U.suffix.ETAS[2]
    L = 14 * lcm(*body)
    require(L == 11760, (case, "unexpected exact-body ruler", L))
    result = Z.E.ray_and_status(first, body, carrier, h, lower, L)
    if result[9]:
        projected = Z.E.projected_packets(
            first, body, carrier, h, lower, L, result[0], result[4], result[9]
        )
    else:
        projected = (0, 0, 0, None, 0, None, None, EMPTY_DIGEST, ())
    counts = (
        len(result[6]),
        len(result[7]),
        len(result[8]),
        len(result[9]),
        projected[1],
    )
    require(counts == EXPECTED_COUNTS[case], (case, "counts changed", counts))
    require(projected[1] == projected[2], (case, "a projected packet survived"))
    audit = (result[10], projected[3], projected[4], projected[6], projected[7])
    require(audit == EXPECTED_EXACT_AUDIT[case], (case, "exact audit changed", audit))
    return (
        first,
        body,
        h,
        len(carrier),
        L,
        lower,
        result[4],
        result[5],
        result[1],
        result[2],
        result[3],
        counts,
        result[10],
        *projected[:-1],
    )


def render(profiles):
    require(len(profiles) == len(CASES) == 2, "profile universe changed")
    totals = tuple(sum(row[11][i] for row in profiles) for i in range(5))
    projected_kills = sum(row[15] for row in profiles)
    require(totals == (212, 34, 158, 20, 20), "global count ledger changed")
    require(projected_kills == totals[-1] == 20, "a packet survived globally")
    positive_margins = tuple(row[16] for row in profiles if row[16] is not None)
    require(min(positive_margins) == F(903353, 7131943), "minimum margin changed")
    require(max(row[17] for row in profiles) == 2, "maximum prefix changed")
    profile_hash = sha256(repr(tuple(profiles)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_payload = (
        EXPECTED_PARENT_SOURCE_SHA256,
        EXPECTED_PARENT_OUTPUT_SHA256,
        FIRST,
        CASES,
        totals,
        projected_kills,
        F(903353, 7131943),
        1724,
        profile_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=2 exact descent closure at z1=1732",
        f"parent_source_sha256={file_sha256(PARENT_SOURCE)}",
        f"parent_output_sha256={file_sha256(PARENT_OUTPUT)}",
        "scope=both scalar-atlas rows;all distinct later nonaligned labels;no finite label horizon",
        "routes=periodic ray/status closure;lossless projected residual only for the second row",
        "atlas=top:1732;empty:1725..1731;next occupied height:1724",
        "projected_identity=P_(E,Z)=phi_L(C_E minus union_z D_z);two-aligned cap=25/91",
        (
            f"global_counts=scalar:{totals[0]};crude:{totals[1]};"
            f"status:{totals[2]};status_survivors:{totals[3]};"
            f"literal_packets:{totals[4]};projected_kills:{projected_kills};survivors:0"
        ),
        "global_minimum_margin=903353/7131943;maximum_prefix_cells=2",
    ]
    for row in profiles:
        counts = row[11]
        lines.append(
            "CASE;"
            f"z1={row[0]};E={','.join(map(str, row[1]))};h={ftext(row[2])};"
            f"carrier_components={row[3]};L={row[4]};lower={ftext(row[5])};"
            f"delta1={ftext(row[6])};first_d={row[7]};ray_sha256={row[8]};"
            f"divisors={row[9]};trials={row[10]};scalar={counts[0]};"
            f"crude={counts[1]};status={counts[2]};status_survivors={counts[3]};"
            f"stage_sha256={'/'.join(row[12])};cells={row[13]};packets={row[14]};"
            f"projected_kills={row[15]};minimum_margin={ftext(row[16])};"
            f"maximum_prefix_cells={row[17]};minimum_row={row[18]};"
            f"direct_mass={ftext(row[19])};state_sha256={row[20]};survivors=0"
        )
    lines.extend(
        (
            "global_z1732_rows=2;empty=2;survivors=0",
            "consequence=projected k=2 first drift label z1<=1724",
            f"profile_sha256={profile_hash}",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    profiles = tuple(profile(case) for case in CASES)
    text = render(profiles)
    args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
