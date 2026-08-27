#!/usr/bin/env python3
"""Exact proof-graph consequences through the THM-4266 CEGAR candidate.

Each carrier deletion is the exhaustively audited semantic high-endpoint set
in the current post-THM-4256 residual, less its literal retained boundary
witness.  The raw round-three closure is retained as a component, then its
overlaps with proved THM-4261 and THM-4262 are removed to obtain THM-4266's
proof-graph-new contribution.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

from reconstruct_post_thm4254_residual import edge_fnv, edge_sha, reconstruct, require


EXPECTED_RAY_SCALES = (
    *range(146, 207),
    *range(208, 213),
    *range(214, 218),
    220,
    221,
    230,
)

EXPECTED_RATIO34_SCALES = (
    *range(97, 167),
    172,
    180,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("repo", type=Path)
    parser.add_argument("--final-edges", action="store_true")
    parser.add_argument("--round2-final-edges", action="store_true")
    parser.add_argument("--post-ray-edges", action="store_true")
    parser.add_argument("--round3-final-edges", action="store_true")
    parser.add_argument("--post-thm4261-4262-edges", action="store_true")
    parser.add_argument("--thm4266-novel-edges", action="store_true")
    parser.add_argument("--thm4266-final-edges", action="store_true")
    args = parser.parse_args()

    _, _, post_4254 = reconstruct(args.repo.resolve())
    ray = [
        edge for edge in post_4254
        if edge[0] * 3 == edge[1] * 2 and edge[0] % 2 == 0
        and edge[0] // 2 >= 146
    ]
    ray_scales = tuple(q // 2 for q, _ in ray)
    require(ray_scales == EXPECTED_RAY_SCALES,
            "semantic THM-4256 ray scales changed")
    require(len(ray) == 73 and edge_fnv(ray) == 0x6BE23222A3A20764,
            "THM-4256 ray ledger changed")
    require(edge_sha(ray) ==
            "9d43ad3311533711e21c496141b05052c45ed68b8b03bcbce59737df3c7391ea",
            "THM-4256 ray SHA changed")
    ray_set = set(ray)
    post_ray = [edge for edge in post_4254 if edge not in ray_set]
    require(len(post_ray) == 180_991 and
            edge_fnv(post_ray) == 0x021BF0ED1581657F,
            "post-THM-4256 residual ledger changed")
    require(edge_sha(post_ray) ==
            "9192c5d73aa5f123ddd10f0115dcaf7231fa518980610042e4cd3f8e73afd44f",
            "post-THM-4256 residual SHA changed")

    high = [edge for edge in post_ray if edge[1] >= 701]
    retained_high = [(542, 732)]
    require(all(edge in high for edge in retained_high),
            "retained high compression witness vanished")
    carrier_closed = [edge for edge in high if edge not in retained_high]
    require(len(high) == 1_935 and len(carrier_closed) == 1_934,
            "carrier high-layer universe changed")
    require(edge_fnv(carrier_closed) == 0xA4C27BE33B8AC001,
            "carrier deletion FNV changed")
    require(edge_sha(carrier_closed) ==
            "ff1031c96058f6165cdcc51cdeff54c0f58ba5cae9fd40acaed7759087789570",
            "carrier deletion SHA changed")
    carrier_set = set(carrier_closed)
    require(not (ray_set & carrier_set), "ray/carrier deletion overlap")

    final = [edge for edge in post_ray if edge not in carrier_set]
    require(len(final) == 179_057, "final residual count changed")
    require(edge_fnv(final) == 0x8101FE1407E16253,
            "final residual FNV changed")
    require(edge_sha(final) ==
            "35905f011237b651918fa90f41df4afb67a44a8543ec0d7a87c281a62dcab81b",
            "final residual SHA changed")
    maximum = max(r for _, r in final)
    top = [edge for edge in final if edge[1] == maximum]
    require(maximum == 732 and top == [(542, 732)],
            "final residual top layer changed")

    round2_high = [edge for edge in post_ray if edge[1] >= 694]
    round2_retained = [(384, 694)]
    require(all(edge in round2_high for edge in round2_retained),
            "round-two retained witness vanished")
    round2_closed = [
        edge for edge in round2_high if edge not in round2_retained
    ]
    require(len(round2_high) == 2_638 and len(round2_closed) == 2_637,
            "round-two high-layer universe changed")
    require(edge_fnv(round2_closed) == 0xE55D3F21BA767347,
            "round-two carrier deletion FNV changed")
    require(edge_sha(round2_closed) ==
            "c456d5bf01ec9b1f9e997b5dc8306dc065db71cbffe61bb8da8bd2f6b16e8d53",
            "round-two carrier deletion SHA changed")
    round2_set = set(round2_closed)
    require(not (ray_set & round2_set),
            "ray/round-two carrier deletion overlap")
    round2_final = [edge for edge in post_ray if edge not in round2_set]
    require(len(round2_final) == 178_354,
            "round-two final residual count changed")
    require(edge_fnv(round2_final) == 0xB1DBDE11F5DA1635,
            "round-two final residual FNV changed")
    require(edge_sha(round2_final) ==
            "b974b2a1e96d68f081a3e8e322b366a89c0e19cb8655bbe139605757f42b3d3a",
            "round-two final residual SHA changed")
    round2_maximum = max(r for _, r in round2_final)
    round2_top = [
        edge for edge in round2_final if edge[1] == round2_maximum
    ]
    require(round2_maximum == 694 and round2_top == [(384, 694)],
            "round-two final top layer changed")

    round3_high = [edge for edge in post_ray if edge[1] >= 688]
    round3_retained = [(520, 688)]
    require(all(edge in round3_high for edge in round3_retained),
            "round-three retained witness vanished")
    round3_closed = [
        edge for edge in round3_high if edge not in round3_retained
    ]
    require(len(round3_high) == 3_337 and len(round3_closed) == 3_336,
            "round-three high-layer universe changed")
    require(edge_fnv(round3_closed) == 0xC1BA162D8A364FB3,
            "round-three carrier deletion FNV changed")
    require(edge_sha(round3_closed) ==
            "95a9c0eb185847f1d64d949cef2b1343e85701dcf416c27f3790c31a91c40854",
            "round-three carrier deletion SHA changed")
    round3_set = set(round3_closed)
    require(not (ray_set & round3_set),
            "ray/round-three carrier deletion overlap")
    round3_final = [edge for edge in post_ray if edge not in round3_set]
    require(len(round3_final) == 177_655,
            "round-three final residual count changed")
    require(edge_fnv(round3_final) == 0x05D884E33AFD6C65,
            "round-three final residual FNV changed")
    require(edge_sha(round3_final) ==
            "50b78270a3d16f7c78249f19b3a39ab4b66aa95be1d9820d6a473e85fef12c06",
            "round-three final residual SHA changed")
    round3_maximum = max(r for _, r in round3_final)
    round3_top = [
        edge for edge in round3_final if edge[1] == round3_maximum
    ]
    require(round3_maximum == 688 and round3_top == [(520, 688)],
            "round-three final top layer changed")

    thm4261_band = [
        edge for edge in post_ray if 733 <= edge[1] <= 754
    ]
    require(len(thm4261_band) == 297 and
            edge_fnv(thm4261_band) == 0xE923D1494185B820,
            "THM-4261 band ledger changed")
    require(edge_sha(thm4261_band) ==
            "745ef7c8809335e6d6e9623314beff917edc71cfaaaa88e7210ede9dcd97d11b",
            "THM-4261 band SHA changed")
    thm4261_set = set(thm4261_band)

    thm4262_ray = [
        edge for edge in post_ray
        if edge[0] * 4 == edge[1] * 3 and edge[0] % 3 == 0
        and edge[0] // 3 >= 97
    ]
    thm4262_scales = tuple(q // 3 for q, _ in thm4262_ray)
    require(thm4262_scales == EXPECTED_RATIO34_SCALES,
            "semantic THM-4262 ray scales changed")
    require(len(thm4262_ray) == 72 and
            edge_fnv(thm4262_ray) == 0x512CBBA28E2235FD,
            "THM-4262 ray ledger changed")
    require(edge_sha(thm4262_ray) ==
            "b1c89073d9b82351b663e97c18c807f03f3fd2d40ddcfafe038d8cad0535cb2c",
            "THM-4262 ray SHA changed")
    thm4262_set = set(thm4262_ray)
    require(not (thm4261_set & thm4262_set),
            "THM-4261/THM-4262 deletion overlap")

    proved_prior = thm4261_set | thm4262_set
    proved_prior_ordered = [edge for edge in post_ray if edge in proved_prior]
    require(len(proved_prior_ordered) == 369 and
            edge_fnv(proved_prior_ordered) == 0xCB13865C00A1A670 and
            edge_sha(proved_prior_ordered) ==
            "61a1c14cf58eeae4010ec0e6b8384e38d24aefa14bac4b70a4aeb7cf5f59c34c",
            "THM-4261+4262 union ledger changed")
    post_4261_4262 = [edge for edge in post_ray if edge not in proved_prior]
    require(len(post_4261_4262) == 180_622 and
            edge_fnv(post_4261_4262) == 0x0CEF4E2887C8F24E,
            "post-THM-4261+4262 residual ledger changed")
    require(edge_sha(post_4261_4262) ==
            "fa1c5672b0f2cd2490413e9b69a4720bf1dc4eef8aee694c1c73d390aba58e11",
            "post-THM-4261+4262 residual SHA changed")

    overlap_4261 = [
        edge for edge in round3_closed if edge in thm4261_set
    ]
    require(overlap_4261 == thm4261_band,
            "raw carrier/THM-4261 overlap changed")
    overlap_4262 = [
        edge for edge in round3_closed if edge in thm4262_set
    ]
    require(overlap_4262 == [(516, 688), (540, 720)] and
            edge_fnv(overlap_4262) == 0xC3A79960D665F57D,
            "raw carrier/THM-4262 overlap changed")
    require(edge_sha(overlap_4262) ==
            "3110bd0e1e067bdb0fb9b5ba16f74afd62dc1e6e2388f96afff4970063ef5427",
            "raw carrier/THM-4262 overlap SHA changed")

    thm4266_novel = [
        edge for edge in round3_closed if edge not in proved_prior
    ]
    require(len(thm4266_novel) == 3_037 and
            edge_fnv(thm4266_novel) == 0x24B36D7047589076,
            "THM-4266 novel deletion ledger changed")
    require(edge_sha(thm4266_novel) ==
            "fcfec867819898ec1a0e1072f27747aec29b6785794328a153bc2b85956ba112",
            "THM-4266 novel deletion SHA changed")
    require(len(thm4266_novel) + len(overlap_4261) +
                len(overlap_4262) == len(round3_closed),
            "raw round-three decomposition changed")
    thm4266_set = set(thm4266_novel)

    current_high = [edge for edge in post_4261_4262 if edge[1] >= 688]
    require(len(current_high) == 3_038 and
            edge_fnv(current_high) == 0x4ED8CEB63E3EDC9E and
            edge_sha(current_high) ==
            "4333d7306bb1f8df464b3bd3261559e9d91e6dc88c835e37518206b8a9d0e643" and
            [edge for edge in current_high if edge not in thm4266_set] ==
                [(520, 688)],
            "current THM-4266 high universe changed")
    thm4266_final = [
        edge for edge in post_4261_4262 if edge not in thm4266_set
    ]
    require(len(thm4266_final) == 177_585 and
            edge_fnv(thm4266_final) == 0x6CE05D05EB01DAED,
            "THM-4266 final residual ledger changed")
    require(edge_sha(thm4266_final) ==
            "009614651bb81e9763b2a9ff4b580497bfb6978a6c69d18cf986346e369374d9",
            "THM-4266 final residual SHA changed")
    thm4266_maximum = max(r for _, r in thm4266_final)
    thm4266_top = [
        edge for edge in thm4266_final if edge[1] == thm4266_maximum
    ]
    require(thm4266_maximum == 688 and thm4266_top == [(520, 688)],
            "THM-4266 final top layer changed")

    if args.final_edges:
        for q, r in final:
            print(f"{q},{r}")
        return
    if args.post_ray_edges:
        for q, r in post_ray:
            print(f"{q},{r}")
        return
    if args.round3_final_edges:
        for q, r in round3_final:
            print(f"{q},{r}")
        return
    if args.round2_final_edges:
        for q, r in round2_final:
            print(f"{q},{r}")
        return
    if args.post_thm4261_4262_edges:
        for q, r in post_4261_4262:
            print(f"{q},{r}")
        return
    if args.thm4266_novel_edges:
        for q, r in thm4266_novel:
            print(f"{q},{r}")
        return
    if args.thm4266_final_edges:
        for q, r in thm4266_final:
            print(f"{q},{r}")
        return

    print("LRC14_CARRIER_CEGAR_PROOF_GRAPH_V1")
    print(
        f"POST_THM4254 COUNT {len(post_4254)} "
        f"FNV {edge_fnv(post_4254):016x} SHA256 {edge_sha(post_4254)}"
    )
    print(
        f"THM4256_RAY_REMOVED COUNT {len(ray)} FNV {edge_fnv(ray):016x} "
        f"SHA256 {edge_sha(ray)} MAX_ENDPOINT {max(r for _, r in ray)}"
    )
    print(
        f"POST_THM4256 COUNT {len(post_ray)} FNV {edge_fnv(post_ray):016x} "
        f"SHA256 {edge_sha(post_ray)}"
    )
    print(
        f"CARRIER_HIGH_UNIVERSE COUNT {len(high)} RANGE_R_GE_701 "
        "RETAINED 542,732"
    )
    print(
        f"CARRIER_REMOVED COUNT {len(carrier_closed)} "
        f"FNV {edge_fnv(carrier_closed):016x} "
        f"SHA256 {edge_sha(carrier_closed)}"
    )
    print(f"RAY_CARRIER_OVERLAP {len(ray_set & carrier_set)}")
    print(
        f"FINAL_RESIDUAL COUNT {len(final)} FNV {edge_fnv(final):016x} "
        f"SHA256 {edge_sha(final)}"
    )
    print(
        f"TOP_LAYER ENDPOINT {maximum} COUNT {len(top)} EDGES "
        + " ".join(f"{q},{r}" for q, r in top)
    )
    print(
        f"ROUND2_CARRIER_HIGH_UNIVERSE COUNT {len(round2_high)} "
        "RANGE_R_GE_694 RETAINED 384,694"
    )
    print(
        f"ROUND2_CARRIER_REMOVED COUNT {len(round2_closed)} "
        f"FNV {edge_fnv(round2_closed):016x} "
        f"SHA256 {edge_sha(round2_closed)}"
    )
    print(f"RAY_ROUND2_CARRIER_OVERLAP {len(ray_set & round2_set)}")
    print(
        f"ROUND2_FINAL_RESIDUAL COUNT {len(round2_final)} "
        f"FNV {edge_fnv(round2_final):016x} "
        f"SHA256 {edge_sha(round2_final)}"
    )
    print(
        f"ROUND2_TOP_LAYER ENDPOINT {round2_maximum} "
        f"COUNT {len(round2_top)} EDGES "
        + " ".join(f"{q},{r}" for q, r in round2_top)
    )
    print(
        f"ROUND3_CARRIER_HIGH_UNIVERSE COUNT {len(round3_high)} "
        "RANGE_R_GE_688 RETAINED 520,688"
    )
    print(
        f"ROUND3_CARRIER_REMOVED COUNT {len(round3_closed)} "
        f"FNV {edge_fnv(round3_closed):016x} "
        f"SHA256 {edge_sha(round3_closed)}"
    )
    print(f"RAY_ROUND3_CARRIER_OVERLAP {len(ray_set & round3_set)}")
    print(
        f"ROUND3_FINAL_RESIDUAL COUNT {len(round3_final)} "
        f"FNV {edge_fnv(round3_final):016x} "
        f"SHA256 {edge_sha(round3_final)}"
    )
    print(
        f"ROUND3_TOP_LAYER ENDPOINT {round3_maximum} "
        f"COUNT {len(round3_top)} EDGES "
        + " ".join(f"{q},{r}" for q, r in round3_top)
    )
    print(
        f"THM4261_INHERITED COUNT {len(thm4261_band)} "
        f"FNV {edge_fnv(thm4261_band):016x} "
        f"SHA256 {edge_sha(thm4261_band)}"
    )
    print(
        f"THM4262_INHERITED COUNT {len(thm4262_ray)} "
        f"FNV {edge_fnv(thm4262_ray):016x} "
        f"SHA256 {edge_sha(thm4262_ray)}"
    )
    print(f"THM4261_THM4262_OVERLAP {len(thm4261_set & thm4262_set)}")
    print(
        f"THM4261_THM4262_UNION COUNT {len(proved_prior_ordered)} "
        f"FNV {edge_fnv(proved_prior_ordered):016x} "
        f"SHA256 {edge_sha(proved_prior_ordered)}"
    )
    print(
        f"POST_THM4261_4262 COUNT {len(post_4261_4262)} "
        f"FNV {edge_fnv(post_4261_4262):016x} "
        f"SHA256 {edge_sha(post_4261_4262)}"
    )
    print(
        f"RAW_ROUND3_OVERLAP_THM4261 COUNT {len(overlap_4261)} "
        f"FNV {edge_fnv(overlap_4261):016x} "
        f"SHA256 {edge_sha(overlap_4261)}"
    )
    print(
        f"RAW_ROUND3_OVERLAP_THM4262 COUNT {len(overlap_4262)} "
        f"EDGES {' '.join(f'{q},{r}' for q, r in overlap_4262)} "
        f"FNV {edge_fnv(overlap_4262):016x} "
        f"SHA256 {edge_sha(overlap_4262)}"
    )
    print(
        f"THM4266_CURRENT_HIGH_UNIVERSE COUNT {len(current_high)} "
        f"FNV {edge_fnv(current_high):016x} SHA256 {edge_sha(current_high)} "
        "RANGE_R_GE_688 RETAINED 520,688"
    )
    print(
        f"THM4266_NOVEL_REMOVED COUNT {len(thm4266_novel)} "
        f"FNV {edge_fnv(thm4266_novel):016x} "
        f"SHA256 {edge_sha(thm4266_novel)}"
    )
    print(
        f"THM4266_FINAL_RESIDUAL COUNT {len(thm4266_final)} "
        f"FNV {edge_fnv(thm4266_final):016x} "
        f"SHA256 {edge_sha(thm4266_final)}"
    )
    print(
        f"THM4266_TOP_LAYER ENDPOINT {thm4266_maximum} "
        f"COUNT {len(thm4266_top)} EDGES "
        + " ".join(f"{q},{r}" for q, r in thm4266_top)
    )
    print("VERDICT EXACT_THM4266_PROOF_GRAPH_CONSEQUENCE LRC14_OPEN")


if __name__ == "__main__":
    main()
