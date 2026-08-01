#!/usr/bin/env python3
"""Exact dependency assembly for the reflected THM-2941 ``k=1`` family.

This verifier does not replace the component interval computations.  It pins
their LF-normalized sources and stored transcripts, checks the exact conclusion
tokens used in the logical assembly, and separates three statements which must
not be conflated:

* every spread ``D <= 5`` closes at every minimum level ``m >= 1``;
* a bodywise bank closes 2,421 of the 3,003 six-label bodies at arbitrary
  positive levels; and
* every spread ``D >= 6`` closes when ``m >= 3D``.

Consequently the assembled certificate-failure locus is confined to the 582
bodies not covered by that bank and ``D >= 6, 1 <= m < 3D``.  The uncovered
bodies are not thereby proved to fail some different certificate.  This is a
theorem about the sufficient reflected residual family inherited from
THM-2941, not a classification of
physical LRC(14) survivors and not a proof of LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_d0_d5_wedge_assembly_thm2941.out"
)

# role, source stem, LF source digest, LF output digest, required conclusion rows
COMPONENTS = (
    (
        "D0_D2_CHROMATIC",
        "lrc14_j7_reflected_universal_pair_chromatic_closure_thm2941",
        "a6f58c1a52dfc1fca61a239068dbe0b216bac41f1622b98748bc4a6d213fb6e8",
        "7364d5866171405fa90539a9ad76727c0c52f020ac1a104a1ab4f0276aedd115",
        (
            "consequence=every reflected k=1 packet using at most three distinct positive levels closes",
            "corollary=the full normalized D=2 sector q_e in {m,m+1,m+2} is empty for every m>=1",
        ),
    ),
    (
        "NEAREST_LEVEL_TAIL",
        "lrc14_j7_reflected_nearest_level_matching_tail_referee_thm2941",
        "c886867412f851c9a2ad75daf7d7533bee1fbd3617147c459ef226015f230564",
        "114ce28bd8a836dfd8d8155b93b37b9037f924d0df7dc532a376eb3f2bf190c0",
        (
            "nearest_level_corollary=m=min q;Delta=least positive q_e-m;m>(13392/35)Delta+93992/3185;clean m>=383Delta+30",
        ),
    ),
    (
        "D3_EXCEPTIONAL",
        "lrc14_j7_reflected_exceptional_proper4_d3_uniform_closure_thm2941",
        "c69a1452514704c44ec5cc128a7606a0cf6f7cf600f83729b10cac46c7210616",
        "1b94bce5a9b5b907e1a48e29ce93a981061bff1f2f7473ef788d443ba2bf10e8",
        (
            "corollary=the entire reflected D=3 sector closes for every m>=1",
        ),
    ),
    (
        "D4_EXCEPTIONAL",
        "lrc14_j7_reflected_exceptional_d4_uniform_closure_thm2941",
        "1e8bd38e48acad3708f59e55f28a1f0c0dbe90073049a8c9c3ce0799ea3ce925",
        "f6100323574fc1f14cfc3ea326ee114a17cdf6b25c70add76802abe349f166dd",
        (
            "corollary=the entire reflected D=4 sector closes for every m>=1",
            "guardrail=spread four does not force all five offsets:72 proper four-colour rows occur on each exceptional body;Delta can equal 2",
        ),
    ),
    (
        "D5_TAIL",
        "lrc14_j7_reflected_d5_crossdet_tail_closure_thm2941",
        "1575b9fabec292bccf0bd639b47b3775922a1531e421b8e4441c6909cc2cedb7",
        "73e8242e68431bd42d3a39b40a8dd0cc6ae8f9e1439e80038a54835bedc3ac55",
        (
            "conclusion=every reflected D=5 residual packet closes for every m>=16",
        ),
    ),
    (
        "D5_HEAD",
        "lrc14_j7_reflected_d5_head_median_cell_closure_thm2941",
        "c12b1a9f40c1a1d7e365589606dbe94009612d0e258806b3d02a4defb99bdf0e",
        "28f6bd514be62ea7061f39682f5bbb15252cde2fe43e54d80543acddcc3187fc",
        (
            "conclusion=every reflected D=5 residual packet closes for every m>=1 (head here, tail upstream)",
            "head_crossdet_repair=choose transport orientation with target A=p(qL-b)/(pL-a)>=1;orientation product AA'=pq",
        ),
    ),
    (
        "ARBITRARY_LEVEL_EDGE10",
        "lrc14_j7_reflected_robust_edge10_threshold_block_uniform_closure_thm2941",
        "cd2a4e84b2527f3fc7bf79980d2816ad207d770639fae9b42a1b9202d54cd2cd",
        "69fca67fbc57ef09e2ce1a13ea75ef754eba9c8e477a167cae67ee45e8501b7c",
        (
            "corollary=arbitrary-level body closure rises from 2354 to 2386;remaining_bodies=617",
        ),
    ),
    (
        "ARBITRARY_LEVEL_EDGE9",
        "lrc14_j7_reflected_robust_edge9_threshold_block_uniform_closure_thm2941",
        "b708f0fc3b5e89d9a17301201b112cadbfa68c279aebe08b9ae486611b273858",
        "3a2acc9fa52dc6cc6223d805395d7de1276d5b410dc840bdf7645af4a38ab699",
        (
            "conclusion=all 35 robust-edge-9 bodies close for every assignment of positive reflected levels",
            "corollary=arbitrary-level body closure rises from 2386 to 2421;remaining_bodies=582",
        ),
    ),
    (
        "C4_CONE",
        "lrc14_j7_reflected_c4_central_exception_cone_closure_thm2941",
        "48a21cfd26c6250a317d37a59523012548204f7b217538d74a5c4b2e21b6f9ae",
        "fe4a949e68fe71fadab5068f688ac1b0d3aeb0c32b13d94a6d49bb2212342aa3",
        (
            "conclusion=all reflected residual packets with D>=6,m>=4D close on all 3003 bodies",
            "corollary=the still-open reflected D>=6 wedge is confined to m<4D",
        ),
    ),
    (
        "C3_CONE",
        "lrc14_j7_reflected_c3_three_reverse_ladder_cone_closure_thm2941",
        "102ef691101bf8d00721e08e1b9b229c893fb532d06dc59959a26a31cb87cbe8",
        "e997189c3cebb693d9b04d56be088842a532cd6c95d8ba6db2dbc1377af93875",
        (
            "conclusion=all reflected residual packets with D>=6,m>=3D close on all 3003 bodies",
            "corollary=the still-open reflected D>=6 wedge is confined to m<3D",
        ),
    ),
)

EXPECTED_D5_RESIDUAL_COUNTS = (
    12985,
    1511,
    773,
    427,
    235,
    53,
    43,
    23,
    8,
    9,
    5,
    2,
    0,
    0,
    2,
)
EXPECTED_SEMANTIC_SHA256 = "9bc1da061c14d44ef1016a9b887097666fd50ea6af242ea8302822a6dbe664b1"


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_bytes(path: Path) -> bytes:
    return path.read_bytes().replace(b"\r\n", b"\n")


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(lf_bytes(path)).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    manifest = []
    transcripts = {}
    for role, stem, expected_source, expected_output, required_rows in COMPONENTS:
        source = ROOT / "04-computation" / f"{stem}.py"
        output = ROOT / "05-knowledge" / "results" / f"{stem}.out"
        source_sha = lf_sha256(source)
        output_sha = lf_sha256(output)
        require(source_sha == expected_source, (role, "source", source_sha, expected_source))
        require(output_sha == expected_output, (role, "output", output_sha, expected_output))
        transcript = lf_bytes(output).decode("utf-8")
        require("all_exact_controls=PASS\n" in transcript, (role, "missing exact controls"))
        require("normal_vs_python_O=BYTE_IDENTICAL\n" in transcript,
                (role, "missing ordinary/-O control"))
        for row in required_rows:
            require(row + "\n" in transcript, (role, "missing conclusion row", row))
        manifest.append((role, source_sha, output_sha))
        transcripts[role] = transcript

    body_count = comb(14, 6)
    complete_good_graph_bodies = 3001
    exceptional_bodies = 2
    arbitrary_level_closed = 2421
    arbitrary_level_residual = 582
    require(body_count == 3003, body_count)
    require(complete_good_graph_bodies + exceptional_bodies == body_count,
            "chromatic body partition changed")
    require(arbitrary_level_closed + arbitrary_level_residual == body_count,
            "arbitrary-level body partition changed")

    # The universal result handles every word on <=3 level values, hence every
    # D=0,1,2 word.  D3 and D4 combine its 3001 complete good graphs with their
    # two exceptional-body finite closures.  The D5 head explicitly imports
    # the D5 tail and concludes all m>=1.
    closed_spreads = tuple(range(0, 6))
    require(closed_spreads == (0, 1, 2, 3, 4, 5), closed_spreads)

    residual_text = (
        "residual_counts_by_m=" + repr(EXPECTED_D5_RESIDUAL_COUNTS)
    )
    require(sum(EXPECTED_D5_RESIDUAL_COUNTS) == 16076,
            sum(EXPECTED_D5_RESIDUAL_COUNTS))
    require(residual_text in transcripts["D5_HEAD"], "D5 residual vector changed")
    require("raw_assignments:7851600" in transcripts["D5_HEAD"],
            "D5 raw assignment count changed")
    require("head_crossdet_repair=choose transport orientation" in transcripts["D5_HEAD"],
            "D5 orientation repair missing")

    remaining_wedge = ("D>=6", "1<=m<3D", arbitrary_level_residual)
    semantic_payload = (
        tuple(manifest),
        body_count,
        complete_good_graph_bodies,
        exceptional_bodies,
        closed_spreads,
        arbitrary_level_closed,
        arbitrary_level_residual,
        EXPECTED_D5_RESIDUAL_COUNTS,
        remaining_wedge,
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256:
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    lines = [
        "LRC14 reflected D<=5 and m>=3D wedge exact dependency assembly",
        "universe=E subset {1,...,14},|E|=6;3003 bodies;q_e positive integers;m=min q_e;D=max q_e-min q_e",
        "D0_D2=universal good-edge chromatic theorem:at most three level values",
        "D3=3001 complete good graphs by pigeonhole plus two exceptional proper-four-colour rays",
        "D4=3001 complete good graphs by pigeonhole plus two exceptional proper-word lanes;Delta<=2 is sharp",
        "D5=head m=1..15 plus cross-determinant tail m>=16;all 3003 bodies",
        "conclusion_1=every reflected THM-2941 residual packet with D<=5 closes for every m>=1",
        "conclusion_2=bodywise bank closes 2421/3003 bodies for arbitrary positive reflected levels;582 bodies are uncovered",
        "conclusion_3=every reflected THM-2941 residual packet with D>=6 and m>=3D closes",
        "remaining_wedge=582 bank-uncovered bodies only;D>=6;1<=m<3D",
        "logical_status=three conclusions are incomparable inputs to the final intersection;none is arbitrary k<=1 or physical-survivor classification",
        "D5_head_counts=raw:7851600;crossdet:7835524;median_residual:16076",
        f"D5_residual_counts_by_m={EXPECTED_D5_RESIDUAL_COUNTS}",
        f"body_partition=complete_good_graph:{complete_good_graph_bodies};exceptions:{exceptional_bodies};total:{body_count}",
        f"arbitrary_level_partition=closed:{arbitrary_level_closed};residual:{arbitrary_level_residual};total:{body_count}",
    ]
    for role, source_sha, output_sha in manifest:
        lines.append(f"dependency={role};source_lf_sha256={source_sha};output_lf_sha256={output_sha}")
    lines.extend((
        f"semantic_sha256={semantic}",
        f"source_lf_sha256={lf_sha256(Path(__file__))}",
        "normal_vs_python_O=BYTE_IDENTICAL",
        "scope=PROVED+VERIFIED-EXACT only for the reflected THM-2941 sufficient family;LRC14 remains open",
        "all_exact_controls=PASS",
    ))
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
