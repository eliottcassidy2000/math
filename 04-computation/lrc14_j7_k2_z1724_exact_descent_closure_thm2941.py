#!/usr/bin/env python3
"""Close both projected k=2 scalar-atlas rows at z1=1724 exactly."""

from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT_SOURCE = ROOT / "04-computation/lrc14_j7_k2_z1732_exact_descent_closure_thm2941.py"
PARENT_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k2_z1732_exact_descent_closure_thm2941.out"
OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_j7_k2_z1724_exact_descent_closure_thm2941.out"
EXPECTED_PARENT_SOURCE_SHA256 = "40718586f9cc537a4e21d85cc8f5e272d5236c487dc4a77cc3c6aaa43007cf1b"
EXPECTED_PARENT_OUTPUT_SHA256 = "d201c2ef454c7895db4801c6fbd552d4658b1188af6d4813c1c24437f4bab978"
EXPECTED_PROFILE_SHA256 = "e6b8de7255cd0b009b24d317de8fa4e7b921efd1bad5f4f3549cafd2448cb735"
EXPECTED_SEMANTIC_SHA256 = "6222d8f465e39b199781fad71902dab4da042c4598085fe351361b8b6f7e75be"

FIRST = 1724
CASES = (
    (FIRST, (1, 4, 8, 10, 12, 14)),
    (FIRST, (2, 4, 8, 10, 12, 14)),
)
# scalar, crude kills, common-status kills, status survivors, literal packets.
EXPECTED_COUNTS = {
    CASES[0]: (462, 91, 370, 1, 1),
    CASES[1]: (116, 8, 24, 84, 97),
}
EMPTY_DIGEST = sha256(b"()").hexdigest()
EXPECTED_EXACT_AUDIT = {
    CASES[0]: (
        (
            "ca094df3d0778ce40b543f249186d0510b4cf0d34b2c1041b3a49555b62b5fd5",
            "a77304b8945bbab76af5d240c744a40bbaa6759ad9c0fe1834c0d41997d01319",
            "224b2bc0f4ad97ad6e6789c4c74cc5b340cc96ca44034a62afb8bf90daa8fce6",
            "0dfc035bf962c388dc39047659c2064a3f79e47659e4f84be608bdfdb02a0430",
        ),
        F(1588507, 4039763),
        2,
        F(1),
        "51869de5c4e729b2d8a5180cbf8db8f37ee9e0f849ee9b4bea6b229df707c3b5",
    ),
    CASES[1]: (
        (
            "9cf44022620714fdbe76b352bdba48b744d47cc68a8c87230b8e170b51ae5713",
            "16560b82d715ede634928645c49613969d431d7a6d4c00f60b3760b423700fc1",
            "b8482ba4e80f17b610279566bd86c6b796b3b9dac1ffd1163bc8b08ecb2effb0",
            "8a1e482c95812241f8ac9aabd9254aee6df1cc7ce4c9785a552cf43c1b60942e",
        ),
        F(681, 2821),
        2,
        F(1),
        "aaa62b017a0397a0bd743385a683c0ce6f5ccca164856c65f744cfd8d3e3adf8",
    ),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return "NONE" if value is None else f"{value.numerator}/{value.denominator}"


require(digest(PARENT_SOURCE) == EXPECTED_PARENT_SOURCE_SHA256, "z1732 source changed")
require(digest(PARENT_OUTPUT) == EXPECTED_PARENT_OUTPUT_SHA256, "z1732 output changed")
SPEC = spec_from_file_location("k2_z1724_parent", PARENT_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1732 parent")
PARENT = module_from_spec(SPEC)
SPEC.loader.exec_module(PARENT)
ENGINE = PARENT.Z.E
ATLAS_KEYS = PARENT.Z.ATLAS_KEYS
require(tuple(key for key in ATLAS_KEYS if key[0] == FIRST) == CASES, "z1724 universe changed")
require(not any(1722 < first < FIRST for first, _body in ATLAS_KEYS), "descent gap occupied")
require(max(first for first, _body in ATLAS_KEYS if first < FIRST) == 1722, "next height changed")


def profile(case):
    first, body = case
    carrier = ENGINE.U.suffix.A.carrier_for(body)
    require(ENGINE.P.A.carrier_for(body) == carrier, (case, "carrier disagreement"))
    h = F(sum(right - left for left, right in carrier), ENGINE.U.suffix.A.RULER)
    lower = h * ENGINE.U.suffix.ETAS[2]
    L = 14 * lcm(*body)
    result = ENGINE.ray_and_status(first, body, carrier, h, lower, L)
    projected = ENGINE.projected_packets(
        first, body, carrier, h, lower, L, result[0], result[4], result[9]
    )
    counts = (
        len(result[6]), len(result[7]), len(result[8]), len(result[9]), projected[1]
    )
    require(counts == EXPECTED_COUNTS[case], (case, "counts changed", counts))
    require(projected[1] == projected[2], (case, "a projected packet survived"))
    audit = (result[10], projected[3], projected[4], projected[6], projected[7])
    require(audit == EXPECTED_EXACT_AUDIT[case], (case, "exact audit changed", audit))
    return (
        first, body, h, len(carrier), L, lower, result[4], result[5], result[1],
        result[2], result[3], counts, result[10], *projected[:-1],
    )


def render(rows):
    totals = tuple(sum(row[11][i] for row in rows) for i in range(5))
    kills = sum(row[15] for row in rows)
    require(totals == (578, 99, 394, 85, 98), "global ledger changed")
    require(kills == totals[-1] == 98, "a packet survived globally")
    require(min(row[16] for row in rows) == F(681, 2821), "minimum margin changed")
    require(max(row[17] for row in rows) == 2, "maximum prefix changed")
    profile_hash = sha256(repr(tuple(rows)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_hash = sha256(repr((
        EXPECTED_PARENT_SOURCE_SHA256, EXPECTED_PARENT_OUTPUT_SHA256, CASES,
        totals, kills, F(681, 2821), 1722, profile_hash,
    )).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=2 exact descent closure at z1=1724",
        f"parent_source_sha256={digest(PARENT_SOURCE)}",
        f"parent_output_sha256={digest(PARENT_OUTPUT)}",
        "scope=both scalar-atlas rows;all distinct later nonaligned labels;no finite label horizon",
        "routes=periodic rays;canonical common-status instances;lossless projected residual",
        "atlas=top:1724;empty:1723;next occupied height:1722",
        "projected_identity=P_(E,Z)=phi_L(C_E minus union_z D_z);two-aligned cap=25/91",
        (
            f"global_counts=scalar:{totals[0]};crude:{totals[1]};status:{totals[2]};"
            f"status_survivors:{totals[3]};literal_packets:{totals[4]};"
            f"projected_kills:{kills};survivors:0"
        ),
        "global_minimum_margin=681/2821;maximum_prefix_cells=2",
    ]
    for row in rows:
        counts = row[11]
        lines.append(
            "CASE;"
            f"z1={row[0]};E={','.join(map(str, row[1]))};h={ftext(row[2])};"
            f"carrier_components={row[3]};L={row[4]};lower={ftext(row[5])};"
            f"delta1={ftext(row[6])};first_d={row[7]};ray_sha256={row[8]};"
            f"divisors={row[9]};trials={row[10]};scalar={counts[0]};crude={counts[1]};"
            f"status={counts[2]};status_survivors={counts[3]};"
            f"stage_sha256={'/'.join(row[12])};cells={row[13]};packets={row[14]};"
            f"projected_kills={row[15]};minimum_margin={ftext(row[16])};"
            f"maximum_prefix_cells={row[17]};minimum_row={row[18]};"
            f"direct_mass={ftext(row[19])};state_sha256={row[20]};survivors=0"
        )
    lines.extend((
        "global_z1724_rows=2;empty=2;survivors=0",
        "consequence=projected k=2 first drift label z1<=1722",
        f"profile_sha256={profile_hash}",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ))
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    text = render(tuple(profile(case) for case in CASES))
    args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
