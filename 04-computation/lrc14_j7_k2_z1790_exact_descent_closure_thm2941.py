#!/usr/bin/env python3
"""Uniform closure of all eight projected k=2 scalar rows at z1=1790.

The exact all-body atlas on 1750..1799 has no row on 1791..1799 and leaves
exactly eight rows at its top occupied height 1790.  This verifier closes the
five exact-suffix rows, and also the third HIGH-TAIL row, by the all-label
pipeline

  residue rays -> denominator multisets -> crude all-divisor capacities
  -> common K5 Hunter status -> exact scalar-slack label packets
  -> lossless projected residual.

This stronger route quantifies over every choice of four distinct later
nonaligned labels, so it does not need the HIGH-TAIL obligation.  The other
two HIGH-TAIL rows die already under the exact forced-high scalar maximum:
the unrestricted top three singleton excesses plus the best wall-eligible
singleton.  No finite label horizon is imposed in either route.
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
    ROOT / "04-computation" / "lrc14_j7_k2_scalar_band_1750_1799_thm2941.py"
)
BAND_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_scalar_band_1750_1799_thm2941.out"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_z1790_exact_descent_closure_thm2941.out"
)

EXPECTED_EXACT_ENGINE_SHA256 = (
    "f864262ef94e97827b6f54b8db3bee38e1420f7af4c31a549157cbbd40d3df86"
)
EXPECTED_EXACT_OUTPUT_SHA256 = (
    "b2470d0acc71cdb731ca443ae6683391e582df9a90aba5bd4192a6c383fcd33b"
)
EXPECTED_HIGH_ENGINE_SHA256 = (
    "a12b60f8ff9c711313f57fcf0993049d455f73ebcc5ed74c06e87ebdaf8fe12d"
)
EXPECTED_HIGH_OUTPUT_SHA256 = (
    "51697a317d0618cba183f5e0ea82875f186f5a6c57ebff9a9a5396a3cefdf34d"
)
EXPECTED_BAND_SOURCE_SHA256 = (
    "9ad7c58575b25c79b41bcf4226710c934720bf3e1394abf573a647eec5af87ff"
)
EXPECTED_BAND_OUTPUT_SHA256 = (
    "2ce806d361d7eb97d9ae2d23e438898c8e1f895a89501c9a1847e51f61ca8009"
)
EXPECTED_PROFILE_SHA256 = (
    "80c30db678b446602ea48ed1f338c64c6e27eea04add82e0a6bb2822e9c9120d"
)
EXPECTED_SEMANTIC_SHA256 = (
    "b5ea0a00c850ec19a4d8c3955dc85d5e933609bee62271baa5881b91a51bc32b"
)

QUANTIFIER = "distinct later nonaligned labels"
HIGH_OBLIGATION = "one later label at or beyond the projected wall"
FIRST = 1790
EXPECTED_HEIGHTS = (
    (1750, 12),
    (1758, 1),
    (1768, 1),
    (1776, 1),
    (1780, 5),
    (1784, 6),
    (1788, 1),
    (1790, 8),
)
EXACT_CASES = (
    (FIRST, (1, 4, 8, 10, 12, 14)),
    (FIRST, (1, 6, 8, 10, 12, 14)),
    (FIRST, (2, 4, 8, 10, 12, 14)),
    (FIRST, (2, 6, 8, 10, 12, 14)),
    (FIRST, (4, 6, 8, 10, 12, 14)),
)
HIGH_SCALAR_CASES = (
    (FIRST, (1, 8, 10, 12, 13, 14)),
    (FIRST, (2, 6, 8, 9, 10, 14)),
)
HIGH_ALL_LABEL_CASE = (FIRST, (2, 8, 9, 10, 12, 14))
ALL_LABEL_CASES = (*EXACT_CASES, HIGH_ALL_LABEL_CASE)
ALL_CASES = tuple(sorted((*ALL_LABEL_CASES, *HIGH_SCALAR_CASES)))
EXPECTED_HIGH_GAPS = {
    HIGH_SCALAR_CASES[0]: F(-87183907, 706985472648),
    HIGH_SCALAR_CASES[1]: F(-57185331157, 1085575216301940),
}

# (scalar, crude kills, common-status kills, status survivors, literal packets)
EXPECTED_COUNTS = {
    EXACT_CASES[0]: (122, 35, 87, 0, 0),
    EXACT_CASES[1]: (10, 4, 5, 1, 1),
    EXACT_CASES[2]: (12, 2, 0, 10, 10),
    EXACT_CASES[3]: (21, 2, 1, 18, 36),
    EXACT_CASES[4]: (1, 0, 1, 0, 0),
    HIGH_ALL_LABEL_CASE: (64, 1, 16, 47, 68),
}

# Filled with the independent exact run below.  Each value is
# (four stage digests, minimum projected margin, maximum prefix cells,
#  projected-state digest).  Keeping the structure explicit makes a changed
# intermediate certificate fail before the aggregate digest.
EXPECTED_EXACT_AUDIT = {
    EXACT_CASES[0]: (
        (
            "7ebb8d9ecdde60d26cfdfd05c7f9c822efec3f2dee693bc4f5289447c3dab79b",
            "30b86634de2c6b61f264a7ccbc12e55c7b54fa202fb10128aa0515ac6060b3e0",
            "db370a5116159096d4c6229fd6da9e35c4c41aa8bc6a761f3ac3638f8dde94e5",
            "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
        ),
        None,
        0,
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
    ),
    EXACT_CASES[1]: (
        (
            "24e26b85aa9cddda1b667b024a56264a1c2d5259a8ae190b46162caed988d77a",
            "3c4afc13c2f7afd0ca2a1841dd47b9650ebbbba144dc32b8d101eb2c8928a13e",
            "9e791e2d81733fc04f8844615a7835de6792ddc93640dbd224d41c44c5f91322",
            "7f530705f83e058fae0ce97693deabd5efaf4140ce5d0934705955b0d08324af",
        ),
        F(66, 91),
        2,
        "5ce619224167b19c4e8b31a16f3302c89a2bc99002f2a6e0e81c843d7bb3c934",
    ),
    EXACT_CASES[2]: (
        (
            "3fff19f7e5d545f466afde922387d57cad871bbd94a0c62077bbdf0346d481dc",
            "1826251bb98c76a39a83cdfd87374ad15efb70ec8f974708768b59475813dedf",
            "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
            "7c3a276e4d85adcaf5aa9395304f53559f76d33dc22193d1a689e393600bc2f8",
        ),
        F(1896001, 2948309),
        2,
        "2cc6504b4ae5df6a7d609615dcf4fdda10f6568addd11550d13996c2ebf44fad",
    ),
    EXACT_CASES[3]: (
        (
            "026e71beb6d15a19077986c49e075fb7f520e81ca8ac9ee31fcc14d2b1a0cf38",
            "cc6dd7effaf0d792808c080699eff6137dd78709712eb19b9d653c5cea16a43f",
            "8bcd8cebf8fcf908ef40d5b6a5e96a1d0dc5201752e2ea5bec9418a9a0e3d5e8",
            "a50aa091fef8852f27e29dcbc6d6240f2d2a7a95896eb763fbdea3244e830a78",
        ),
        F(45959, 477659),
        4,
        "7b42d4bed4da115ec917e20efc88329207427cf81b2b3b60ce6584e75a85f887",
    ),
    EXACT_CASES[4]: (
        (
            "e9c03f269420efa4a66e29d1c6e1e8b22eebe9abb646caf7c8a834446ec76c85",
            "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
            "8f1640967eeeacec91b1d18b46bb3a5134900d0ddead79117934a227d224f06f",
            "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
        ),
        None,
        0,
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
    ),
    HIGH_ALL_LABEL_CASE: (
        (
            "7bda59f70f20fdb8f9f139143db6c80d9917a9aac538b292c190d07229240bd9",
            "96e5052ef61a29079b92b63dde0cdcd8041c333334bbd6fc4268c55162a8eb68",
            "4a1354ed04812bd8f57c738036bf32461e3d5a9450f0a04887ae933ec1ed792c",
            "64e3e9e03165824b06003368343779f06f748a3b745a3c010acf41ae960ad36a",
        ),
        F(891, 5369),
        11,
        "fe4b810ce6c833704d4d96889a7fe7c58f029f9ebf15a7552bf5a45c002a8478",
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


require(
    file_sha256(EXACT_ENGINE) == EXPECTED_EXACT_ENGINE_SHA256,
    "exact descent engine changed",
)
require(
    file_sha256(EXACT_OUTPUT) == EXPECTED_EXACT_OUTPUT_SHA256,
    "exact descent output changed",
)
require(
    file_sha256(HIGH_ENGINE) == EXPECTED_HIGH_ENGINE_SHA256,
    "high-wall engine changed",
)
require(
    file_sha256(HIGH_OUTPUT) == EXPECTED_HIGH_OUTPUT_SHA256,
    "high-wall output changed",
)
require(
    file_sha256(BAND_SOURCE) == EXPECTED_BAND_SOURCE_SHA256,
    "1750..1799 band source changed",
)
require(
    file_sha256(BAND_OUTPUT) == EXPECTED_BAND_OUTPUT_SHA256,
    "1750..1799 band output changed",
)

ESPEC = spec_from_file_location("k2_z1790_exact_engine", EXACT_ENGINE)
require(ESPEC is not None and ESPEC.loader is not None, "cannot load exact engine")
E = module_from_spec(ESPEC)
ESPEC.loader.exec_module(E)
HSPEC = spec_from_file_location("k2_z1790_high_engine", HIGH_ENGINE)
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
require(max(first for first, _body in ATLAS_KEYS) == FIRST, "band top changed")
require(
    tuple(key for key in ATLAS_KEYS if key[0] == FIRST) == ALL_CASES,
    "z1790 case universe changed",
)
require(
    tuple(key for key in ATLAS_HIGH_KEYS if key[0] == FIRST)
    == tuple(sorted((*HIGH_SCALAR_CASES, HIGH_ALL_LABEL_CASE))),
    "z1790 HIGH-TAIL partition changed",
)
require(not any(first == 1789 for first, _body in ATLAS_KEYS), "z1789 became occupied")


def all_label_profile(case):
    first, body = case
    carrier = E.U.suffix.A.carrier_for(body)
    require(E.P.A.carrier_for(body) == carrier, (case, "carrier engines disagree"))
    h = F(sum(right - left for left, right in carrier), E.U.suffix.A.RULER)
    lower = h * E.U.suffix.ETAS[2]
    L = 14 * lcm(*body)
    require(E.P.A.RULER % L == 0, (case, "body ruler left master ruler", L))
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
        projected = (0, 0, 0, None, 0, None, None, sha256(b"()").hexdigest(), ())
    counts = (
        len(scalar),
        len(crude_kills),
        len(status_kills),
        len(states),
        projected[1],
    )
    require(counts == EXPECTED_COUNTS[case], (case, "stage counts changed", counts))
    require(
        projected[1] == projected[2],
        (case, "a scalar-eligible projected packet survived"),
    )
    if projected[1]:
        require(projected[3] is not None and projected[3] > 0, (case, "bad margin"))
    audit = (stage_digests, projected[3], projected[4], projected[7])
    expected_audit = EXPECTED_EXACT_AUDIT[case]
    if CHECK_EXPECTED_EXACT_AUDIT and expected_audit is not None:
        require(audit == expected_audit, (case, "exact audit changed", audit))
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
    case = (first, body)
    if route == "ALL-LABEL":
        return all_label_profile(case)
    require(route == "HIGH-SCALAR", (task, "unknown route"))
    return (route, *H.profile(case))


TASKS = tuple(
    [("ALL-LABEL", first, body) for first, body in ALL_LABEL_CASES]
    + [("HIGH-SCALAR", first, body) for first, body in HIGH_SCALAR_CASES]
)


def render(profiles):
    require(len(profiles) == len(TASKS) == 8, "profile universe changed")
    exact_rows = tuple(row for row in profiles if row[0] == "ALL-LABEL")
    high_rows = tuple(row for row in profiles if row[0] == "HIGH-SCALAR")
    require(len(exact_rows) == 6 and len(high_rows) == 2, "route partition changed")
    require(all(row[-1] < 0 for row in high_rows), "a high scalar row survived")
    scalar_states = sum(row[12][0] for row in exact_rows)
    crude_kills = sum(row[12][1] for row in exact_rows)
    status_kills = sum(row[12][2] for row in exact_rows)
    status_survivors = sum(row[12][3] for row in exact_rows)
    packets = sum(row[15] for row in exact_rows)
    projected_kills = sum(row[16] for row in exact_rows)
    require(
        (scalar_states, crude_kills, status_kills, status_survivors, packets, projected_kills)
        == (230, 44, 110, 76, 115, 115),
        "global exact ledger changed",
    )
    profile_hash = sha256(repr(tuple(profiles)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_payload = (
        EXPECTED_EXACT_ENGINE_SHA256,
        EXPECTED_HIGH_ENGINE_SHA256,
        EXPECTED_BAND_SOURCE_SHA256,
        QUANTIFIER,
        HIGH_OBLIGATION,
        EXPECTED_HEIGHTS,
        ALL_LABEL_CASES,
        HIGH_SCALAR_CASES,
        scalar_states,
        crude_kills,
        status_kills,
        status_survivors,
        packets,
        projected_kills,
        tuple(row[-1] for row in high_rows),
        profile_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=2 exact descent closure at z1=1790",
        f"exact_engine_sha256={file_sha256(EXACT_ENGINE)}",
        f"exact_output_sha256={file_sha256(EXACT_OUTPUT)}",
        f"high_engine_sha256={file_sha256(HIGH_ENGINE)}",
        f"high_output_sha256={file_sha256(HIGH_OUTPUT)}",
        f"scalar_band_source_sha256={file_sha256(BAND_SOURCE)}",
        f"scalar_band_output_sha256={file_sha256(BAND_OUTPUT)}",
        (
            f"scope=six rows over all {QUANTIFIER};two rows under {HIGH_OBLIGATION};"
            "no finite label horizon"
        ),
        (
            "routes=six all-label ray/status/projected closures;"
            "two exact forced-high scalar maxima"
        ),
        "atlas=top:1790;empty:1791..1799 and 1789;next occupied height:1788",
        "projected_identity=P_(E,Z)=phi_L(C_E minus union_z D_z);two-aligned cap=25/91",
        (
            f"all_label_counts=scalar:{scalar_states};crude:{crude_kills};"
            f"status:{status_kills};status_survivors:{status_survivors};"
            f"literal_packets:{packets};projected_kills:{projected_kills};survivors:0"
        ),
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
            "global_z1790_rows=8;empty=8;survivors=0",
            "consequence=projected k=2 first drift label z1<=1788",
            f"profile_sha256={profile_hash}",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, mp.cpu_count() or 1))
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
