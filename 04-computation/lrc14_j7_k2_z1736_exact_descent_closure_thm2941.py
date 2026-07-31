#!/usr/bin/env python3
"""Uniform closure of all fifteen projected k=2 rows at z1=1736.

The exact all-body atlas on 1680..1742 has no scalar row on 1737..1742
and leaves exactly fifteen rows at its top occupied height 1736.  Nine have
exact suffix maximizers.  Of the six HIGH-TAIL rows, five die under the exact
forced-high decreasing-ray maximum.  The last is treated under the stronger
unrestricted all-label relaxation.

For the ten all-label rows this verifier runs the complete lossless chain

  residue rays -> denominator multisets -> all-divisor capacity
  -> common K5 Hunter status -> finite scalar-slack literal packets
  -> projected residual.

The chain closes every packet, including all unrestricted packets on the
hard HIGH-TAIL body (1,8,10,12,13,14).  There is no finite label horizon in
the ray/status/packet stages; the horizon 7000 belongs only to the rigorous
parent scalar envelope and carries its separate omitted-slot bound.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXACT_ENGINE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k2_exact_descent_1800_1824_closure_thm2941.py"
)
EXACT_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_exact_descent_1800_1824_closure_thm2941.out"
)
HIGH_ENGINE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k2_high_wall_descent_1800_1810_closure_thm2941.py"
)
HIGH_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_high_wall_descent_1800_1810_closure_thm2941.out"
)
BAND_SOURCE = (
    ROOT / "04-computation" / "lrc14_j7_k2_scalar_band_1680_1742_thm2941.py"
)
BAND_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_scalar_band_1680_1742_thm2941.out"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_z1736_exact_descent_closure_thm2941.out"
)

EXPECTED_EXACT_ENGINE_SHA256 = (
    "1bc6674fbb9b6f4c8979c229c164d267e60911ed582fb3184813d45c21da2adf"
)
EXPECTED_EXACT_OUTPUT_SHA256 = (
    "8d65c8758897bc97c9c463280c1d833e0570ef4763553fdb0c187d1a8996d42e"
)
EXPECTED_HIGH_ENGINE_SHA256 = (
    "c12c17297aa8a96cbdd7d9d529838c776b160ace86b92d0030f5df447fe6877b"
)
EXPECTED_HIGH_OUTPUT_SHA256 = (
    "8038d1c15e157998abd1c8995a92b0ec5b28436a7d04cc2406468a079eecca41"
)
EXPECTED_BAND_SOURCE_SHA256 = (
    "1224d5594571f21c91c55fe3ab165c4fc34ba7968719862d12660d24efac919d"
)
EXPECTED_BAND_OUTPUT_SHA256 = (
    "c20607cb478ed491d7000f2b8a49213f57d1606a5152700ac3b50c836e2dc66c"
)
EXPECTED_PROFILE_SHA256 = (
    "c1047c004b372c2fb742969b28b5a2bef3059ba8e14b5c13c670fe299dff6fbe"
)
EXPECTED_SEMANTIC_SHA256 = (
    "0e65288aa44704bf5a49075ae561d8166aa852f9a1ff3261c0096b1023025d2c"
)

FIRST = 1736
QUANTIFIER = "distinct later nonaligned labels"
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

ORDINARY_CASES = tuple(
    (FIRST, body)
    for body in (
        (1, 2, 8, 10, 12, 14),
        (1, 4, 8, 10, 12, 14),
        (1, 6, 8, 10, 12, 14),
        (2, 3, 8, 10, 12, 14),
        (2, 4, 6, 8, 10, 14),
        (2, 4, 8, 10, 12, 14),
        (2, 6, 8, 10, 12, 14),
        (3, 4, 8, 10, 12, 14),
        (4, 6, 8, 10, 12, 14),
    )
)
HIGH_ALL_LABEL_CASE = (FIRST, (1, 8, 10, 12, 13, 14))
ALL_LABEL_CASES = (*ORDINARY_CASES, HIGH_ALL_LABEL_CASE)
HIGH_SCALAR_CASES = tuple(
    (FIRST, body)
    for body in (
        (1, 2, 10, 11, 12, 14),
        (1, 2, 10, 12, 13, 14),
        (1, 4, 10, 12, 13, 14),
        (2, 3, 10, 11, 12, 14),
        (2, 8, 9, 10, 12, 14),
    )
)
ALL_CASES = tuple(sorted((*ALL_LABEL_CASES, *HIGH_SCALAR_CASES)))
ALL_HIGH_CASES = tuple(sorted((*HIGH_SCALAR_CASES, HIGH_ALL_LABEL_CASE)))

EXPECTED_HIGH_GAPS = {
    HIGH_SCALAR_CASES[0]: F(-17704278379, 151099294846500),
    HIGH_SCALAR_CASES[1]: F(-2122333151, 17733467524500),
    HIGH_SCALAR_CASES[2]: F(-7353136579, 118459563063660),
    HIGH_SCALAR_CASES[3]: F(-47509753883, 178751453380500),
    HIGH_SCALAR_CASES[4]: F(-7891025249, 164964861300785),
}

# scalar, crude kills, common-status kills, status survivors, literal packets.
EXPECTED_COUNTS = {
    ORDINARY_CASES[0]: (46, 0, 16, 30, 43),
    ORDINARY_CASES[1]: (617, 0, 611, 6, 8),
    ORDINARY_CASES[2]: (30, 0, 23, 7, 10),
    ORDINARY_CASES[3]: (3, 0, 1, 2, 2),
    ORDINARY_CASES[4]: (1, 0, 1, 0, 0),
    ORDINARY_CASES[5]: (198, 0, 109, 89, 194),
    ORDINARY_CASES[6]: (30, 0, 9, 21, 29),
    ORDINARY_CASES[7]: (3, 0, 3, 0, 0),
    ORDINARY_CASES[8]: (1, 0, 1, 0, 0),
    HIGH_ALL_LABEL_CASE: (827, 2, 607, 218, 749),
}

EMPTY_DIGEST = sha256(b"()").hexdigest()
EXPECTED_EXACT_AUDIT = {
    ORDINARY_CASES[0]: (
        (
            "d7f2e2b289c739bb667af041a6617fc915ce88f0ddd47d6e56837057cea99cb7",
            EMPTY_DIGEST,
            "f63ad05ba13638760bab101f00e442afadd564c388f9b31a4ba8e42507a87f4b",
            "ac04d1fbda030f66542a66e06989e6abbe92b2593534ef3fc1bd4b0d465a1131",
        ),
        F(23897, 1434069),
        4,
        "77491cb5de95588fee59c489e0f8427b9365799aa40d02350f8c68765e35e5e3",
    ),
    ORDINARY_CASES[1]: (
        (
            "45d8ffc36a4488ae7e95a7dffd7f7f5840c4d00c198696c11f5fbc8f72b58ebe",
            EMPTY_DIGEST,
            "098f806b570e9c46c9dd886a9cd4b9c3071d0d7a51bb7b8231f3db0d7bf24b9f",
            "97997123d9471015e1971c89dbf2a4e74eaf2c0e21af90603532affe64207b6e",
        ),
        F(161107, 6585579),
        3,
        "9cf36665aa9b9aa3b5d1845bb8abc3d8b23d1c917610ba951cb2c1cd651cc490",
    ),
    ORDINARY_CASES[2]: (
        (
            "60fbbf3a845c112400e2d50b49879b4edfb7dd09ffccc9d6444daccdee212a3b",
            EMPTY_DIGEST,
            "9ebd1120f093db27c7d7fa7a6817489c9d521051c1d7934ff65404033c37d190",
            "c6a2125986b0135103f443cd8b8cd5b3873d15477a735baad86b02e4ec02781a",
        ),
        F(590, 2821),
        2,
        "ae15f835044caddbfe3a81d637ce360b8d7d3a331c7e9906d2bbbbf86405a2d1",
    ),
    ORDINARY_CASES[3]: (
        (
            "cebb57972815261050413c464f8e1a83754632091322901fcc57f738d3551692",
            EMPTY_DIGEST,
            "c712d71c153772e7eb342ffd7375756ca23de4ea849f59bad04bd31915ad1958",
            "650d42500026f5fb381e6a6fd204a0555b014dabd4b094e90e434f145e7e94df",
        ),
        F(681, 2821),
        1,
        "bb6d1f3147315ebda37641915eb8a385283bfad504e12c5e6fd4e27162958439",
    ),
    ORDINARY_CASES[4]: (
        (
            "d6891c20082b2dd5793b175241f569fc3227a8b44e7626b1ae136018c3707d8a",
            EMPTY_DIGEST,
            "69a8c539f9c849aad0470f170d4cc72aa14f9369bb239d658cc6ab657d6141f7",
            EMPTY_DIGEST,
        ),
        None,
        0,
        EMPTY_DIGEST,
    ),
    ORDINARY_CASES[5]: (
        (
            "31dfd0f4e4ed72f597d4906f4267d649bd5095833b23ead5352175dbc5b4bf9f",
            EMPTY_DIGEST,
            "809e5388ba976146aed82124a982900d4ef4b53b1e74c802df460a9be87842c5",
            "7b2555bbbd300cd6b33fa4505587bfeb2991837a2e10625575acd9a7d48bd454",
        ),
        F(681, 2821),
        2,
        "d3484bf4398240d2a64b78dd5f4075037bfde530a599d8f9fea77f9123140748",
    ),
    ORDINARY_CASES[6]: (
        (
            "1261a50f3b553e4b6b96d573fb5c3207ec5679687d47247db37dbea669d8c825",
            EMPTY_DIGEST,
            "cd6f775c81d5516783df19d9514ee20a2b2f5b3a52bf72e4a56d7b5159187d69",
            "de7635dde176a1f80604e781a7102d5869880883e668773c78c1b4511821233b",
        ),
        F(49915, 510601),
        2,
        "34c8a735ae7198acdef49ff96ceff9b151226e97cf0880c98abbba038f4b6e1c",
    ),
    ORDINARY_CASES[7]: (
        (
            "c659774c4dbcdce1977a3dd711a88e955c7d10772d170bca6485ca8039de9a7d",
            EMPTY_DIGEST,
            "7464e0b63efdbf03d659bc90412f19a9edc130bfdde91b35c0ff16f034768cb4",
            EMPTY_DIGEST,
        ),
        None,
        0,
        EMPTY_DIGEST,
    ),
    ORDINARY_CASES[8]: (
        (
            "e8f6ee0c0fdfe775f96b8041771f97cfc7d2d615177d10e9f62adbe99d367f1c",
            EMPTY_DIGEST,
            "1e39ca5ff073368d3dd4cd5177efb89916889dab47e9134d1626f2c6484b7eb7",
            EMPTY_DIGEST,
        ),
        None,
        0,
        EMPTY_DIGEST,
    ),
    HIGH_ALL_LABEL_CASE: (
        (
            "c7767065928c6fe39ac1ca28ee3ea74a5181d8ee4ebd8ea3d8049b7c3e39d3a7",
            "40ea404638775bddc1cd84d1a37e92d19bcd1c139d0bf652b6fc99039e7183f2",
            "15ef06189a2497cfb5cf7b52c51707840b6c5bec3d1498e0914efda97ca683bd",
            "2f484f27da3fb5f049044c3b2233ea185b27bbd0f52cbc7ce32db8983b310597",
        ),
        F(121, 18109),
        31,
        "aa644ef827e84923f548a30c13d415d02f8633b949787b115f9d7b5925774ae5",
    ),
}
CHECK_EXPECTED_EXACT_AUDIT = True


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    if value is None:
        return "NONE"
    return f"{value.numerator}/{value.denominator}"


for path, expected, label in (
    (EXACT_ENGINE, EXPECTED_EXACT_ENGINE_SHA256, "exact engine"),
    (EXACT_OUTPUT, EXPECTED_EXACT_OUTPUT_SHA256, "exact output"),
    (HIGH_ENGINE, EXPECTED_HIGH_ENGINE_SHA256, "high engine"),
    (HIGH_OUTPUT, EXPECTED_HIGH_OUTPUT_SHA256, "high output"),
    (BAND_SOURCE, EXPECTED_BAND_SOURCE_SHA256, "scalar band"),
    (BAND_OUTPUT, EXPECTED_BAND_OUTPUT_SHA256, "scalar output"),
):
    require(file_sha256(path) == expected, f"{label} changed")

ESPEC = spec_from_file_location("k2_z1736_exact_engine", EXACT_ENGINE)
require(ESPEC is not None and ESPEC.loader is not None, "cannot load exact engine")
E = module_from_spec(ESPEC)
ESPEC.loader.exec_module(E)
HSPEC = spec_from_file_location("k2_z1736_high_engine", HIGH_ENGINE)
require(HSPEC is not None and HSPEC.loader is not None, "cannot load high engine")
H = module_from_spec(HSPEC)
HSPEC.loader.exec_module(H)
H.EXPECTED_GAPS.update(EXPECTED_HIGH_GAPS)


def atlas_partition():
    keys = []
    high_keys = []
    heights = Counter()
    for line in BAND_OUTPUT.read_text(encoding="utf-8").splitlines():
        if not line.startswith("SURVIVOR;"):
            continue
        fields = dict(
            field.split("=", 1) for field in line.split(";")[1:] if "=" in field
        )
        key = (int(fields["z1"]), tuple(map(int, fields["E"].split(","))))
        keys.append(key)
        heights[key[0]] += 1
        if "HIGH-TAIL" in line:
            high_keys.append(key)
    return tuple(sorted(keys)), tuple(sorted(high_keys)), tuple(sorted(heights.items()))


ATLAS_KEYS, ATLAS_HIGH_KEYS, ATLAS_HEIGHTS = atlas_partition()
require(ATLAS_HEIGHTS == EXPECTED_HEIGHTS, "band height profile changed")
require(tuple(key for key in ATLAS_KEYS if key[0] == FIRST) == ALL_CASES, "z1736 universe changed")
require(tuple(key for key in ATLAS_HIGH_KEYS if key[0] == FIRST) == ALL_HIGH_CASES, "HIGH partition changed")
require(not any(FIRST < first <= 1742 for first, _body in ATLAS_KEYS), "upper gap occupied")
require(max(first for first, _body in ATLAS_KEYS if first < FIRST) == 1732, "next height changed")


def all_label_profile(case):
    first, body = case
    carrier = E.U.suffix.A.carrier_for(body)
    require(E.P.A.carrier_for(body) == carrier, (case, "carrier engines disagree"))
    h = F(sum(right - left for left, right in carrier), E.U.suffix.A.RULER)
    lower = h * E.U.suffix.ETAS[2]
    L = 14 * lcm(*body)
    (
        amplitudes,
        ray_digest,
        divisor_count,
        trials,
        first_delta,
        first_d,
        scalar,
        crude_kills,
        status_kills,
        states,
        stage_digests,
    ) = E.ray_and_status(first, body, carrier, h, lower, L)
    if states:
        projected = E.projected_packets(
            first, body, carrier, h, lower, L, amplitudes, first_delta, states
        )
    else:
        projected = (0, 0, 0, None, 0, None, None, EMPTY_DIGEST, ())
    counts = (
        len(scalar),
        len(crude_kills),
        len(status_kills),
        len(states),
        projected[1],
    )
    require(counts == EXPECTED_COUNTS[case], (case, "counts changed", counts))
    require(projected[1] == projected[2], (case, "a projected packet survived"))
    audit = (stage_digests, projected[3], projected[4], projected[7])
    if CHECK_EXPECTED_EXACT_AUDIT:
        require(audit == EXPECTED_EXACT_AUDIT[case], (case, "audit changed", audit))
    return (
        "ALL-LABEL",
        first,
        body,
        h,
        len(carrier),
        L,
        lower,
        first_delta,
        first_d,
        ray_digest,
        divisor_count,
        trials,
        counts,
        stage_digests,
        *projected[:-1],
    )


def profile(task):
    route, first, body = task
    if route == "ALL-LABEL":
        return all_label_profile((first, body))
    require(route == "HIGH-SCALAR", (task, "unknown route"))
    return (route, *H.profile((first, body)))


TASKS = tuple(
    [("ALL-LABEL", first, body) for first, body in ALL_LABEL_CASES]
    + [("HIGH-SCALAR", first, body) for first, body in HIGH_SCALAR_CASES]
)


def render(profiles):
    require(len(profiles) == len(TASKS) == 15, "profile universe changed")
    exact_rows = tuple(row for row in profiles if row[0] == "ALL-LABEL")
    high_rows = tuple(row for row in profiles if row[0] == "HIGH-SCALAR")
    require(len(exact_rows) == 10 and len(high_rows) == 5, "route partition changed")
    require(all(row[-1] < 0 for row in high_rows), "a forced-high row survived")
    totals = tuple(sum(row[12][i] for row in exact_rows) for i in range(5))
    projected_kills = sum(row[16] for row in exact_rows)
    require(totals == (1756, 2, 1381, 373, 1035), "global exact ledger changed")
    require(projected_kills == totals[-1] == 1035, "a packet survived globally")
    positive_margins = tuple(row[17] for row in exact_rows if row[17] is not None)
    require(min(positive_margins) == F(121, 18109), "minimum margin changed")
    require(max(row[18] for row in exact_rows) == 31, "maximum prefix changed")
    profile_hash = sha256(repr(tuple(profiles)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_payload = (
        EXPECTED_EXACT_ENGINE_SHA256,
        EXPECTED_HIGH_ENGINE_SHA256,
        EXPECTED_BAND_SOURCE_SHA256,
        FIRST,
        QUANTIFIER,
        ALL_LABEL_CASES,
        HIGH_SCALAR_CASES,
        totals,
        projected_kills,
        tuple(row[-1] for row in high_rows),
        profile_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=2 exact descent closure at z1=1736",
        f"exact_engine_sha256={file_sha256(EXACT_ENGINE)}",
        f"exact_output_sha256={file_sha256(EXACT_OUTPUT)}",
        f"high_engine_sha256={file_sha256(HIGH_ENGINE)}",
        f"high_output_sha256={file_sha256(HIGH_OUTPUT)}",
        f"scalar_band_source_sha256={file_sha256(BAND_SOURCE)}",
        f"scalar_band_output_sha256={file_sha256(BAND_OUTPUT)}",
        (
            f"scope=ten rows over all {QUANTIFIER};five rows under one forced-high "
            "wall obligation;no finite label horizon"
        ),
        (
            "routes=ten unrestricted all-label ray/status/projected closures;"
            "five exact forced-high scalar maxima"
        ),
        "atlas=top:1736;empty:1737..1742;next occupied height:1732",
        "projected_identity=P_(E,Z)=phi_L(C_E minus union_z D_z);two-aligned cap=25/91",
        (
            f"all_label_counts=scalar:{totals[0]};crude:{totals[1]};"
            f"status:{totals[2]};status_survivors:{totals[3]};"
            f"literal_packets:{totals[4]};projected_kills:{projected_kills};survivors:0"
        ),
        "global_minimum_margin=121/18109;maximum_prefix_cells=31",
    ]
    for row in profiles:
        if row[0] == "ALL-LABEL":
            (
                _route,
                first,
                body,
                h,
                components,
                L,
                lower,
                first_delta,
                first_d,
                ray_digest,
                divisor_count,
                trials,
                counts,
                stage_digests,
                body_cells,
                packet_count,
                killed_count,
                minimum_margin,
                maximum_cells,
                minimum_row,
                direct_mass,
                state_digest,
            ) = row
            lines.append(
                f"CASE;route=ALL-LABEL;z1={first};E={','.join(map(str, body))};"
                f"h={ftext(h)};r={components};L={L};lower={ftext(lower)};"
                f"delta1={ftext(first_delta)};first_d={first_d};"
                f"ray_sha256={ray_digest};denominators={divisor_count};trials={trials};"
                f"counts={counts};stage_sha256={stage_digests};body_cells={body_cells};"
                f"packets={packet_count};kills={killed_count};"
                f"min_margin={ftext(minimum_margin)};max_prefix_cells={maximum_cells};"
                f"direct_control_mass={ftext(direct_mass)};state_sha256={state_digest};"
                f"minimum_row={minimum_row};conclusion=EMPTY"
            )
        else:
            (
                _route,
                first,
                body,
                h,
                components,
                L,
                high_floor,
                first_delta,
                lower,
                positive,
                negative,
                zero,
                ray_digest,
                rank4,
                omitted_max,
                best_high,
                branch,
                constrained,
                upper,
                gap,
            ) = row
            lines.append(
                f"CASE;route=HIGH-SCALAR;z1={first};E={','.join(map(str, body))};"
                f"h={ftext(h)};r={components};L={L};high_floor={high_floor};"
                f"delta1={ftext(first_delta)};lower={ftext(lower)};"
                f"ray_signs=+{positive}/-{negative}/0:{zero};ray_sha256={ray_digest};"
                f"unrestricted_top4={rank4};first_omitted={omitted_max};"
                f"best_high={best_high};branch={branch};constrained={constrained};"
                f"upper={ftext(upper)};gap={ftext(gap)};conclusion=SCALAR-EMPTY"
            )
    lines.extend(
        (
            "global_z1736_rows=15;empty=15;survivors=0",
            "consequence=projected k=2 first drift label z1<=1732",
            f"profile_sha256={profile_hash}",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(4, mp.cpu_count() or 1))
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    require(args.hash_seed >= 0, "hash seed must be nonnegative")
    os.environ["PYTHONHASHSEED"] = str(args.hash_seed)
    if args.workers == 1:
        profiles = [profile(task) for task in TASKS]
    else:
        with mp.get_context("spawn").Pool(min(args.workers, len(TASKS))) as pool:
            profiles = list(pool.imap(profile, TASKS))
    order = {task: index for index, task in enumerate(TASKS)}
    profiles.sort(key=lambda row: order[(row[0], row[1], row[2])])
    output = render(tuple(profiles))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
