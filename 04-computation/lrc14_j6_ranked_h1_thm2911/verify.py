#!/usr/bin/env python3
"""Locked verifier for the THM-2911 finite ranked-H1 closure package.

The package has four layers.

1. Thirty-two hash-pinned shard outputs certify the exact finite r4/H1
   census and its 132 route roots (127 hard-terminal plus five
   scalar-terminal).
2. The hash-pinned Hunter output certifies the 24 nonstar repairs left by
   the ordered-pivot test.
3. The hardened THM-2905 parser reconstructs all 14,806 hard branch keys,
   retaining the excluded prefix, and reconstructs the G5-positive keys.
4. A hash-pinned THM-2904 hostile-centre ledger recovers its 279
   pivot-closed keys and composes all three branch routes exactly.

This verifier does not rerun the interval scans in the shards.  The scout
and Hunter sources are included so those computations can be replayed;
the THM-2904 source reproduces its promoted hostile-centre ledger.
"""

from __future__ import annotations

import ast
import hashlib
import importlib.util
import math
import re
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[2]
RESULTS = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j6_ranked_h1_thm2911"
)
SCOUT = HERE / "scout.py"
HUNTER_BASE = HERE / "hunter_two_star_exceptions.py"
HUNTER_ALL = HERE / "hunter_all_24_star_survivors.py"
HUNTER_OUTPUT = RESULTS / "hunter_all_24_star_survivors.out"
REPLAY_OUTPUT = RESULTS / "ordinary_optimized_replay.out"
G5_SOURCE = (
    ROOT
    / "04-computation"
    / "lrc14_j6_all_hard_hunter_star_envelope_census_codex_20260729.py"
)
PIVOT_SOURCE = (
    ROOT
    / "04-computation"
    / "lrc14_j6_all_hard_ranked_h1_hunter_pivot_census_codex_20260729.py"
)
PIVOT_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j6_all_hard_ranked_h1_hunter_pivot_census_codex_20260729.out"
)
PIVOT_LEDGER = RESULTS / "thm2904_hostile_centre.ledger.out"

SCOUT_SHA256 = (
    "6abc06972cda64bb4f53db26cca62bbea5362ca178bdbf5bf7398e8b0f28317a"
)
HUNTER_BASE_SHA256 = (
    "c452781c7cf6a2ada6be8984b9c9cbe7aab8369c5a9333721ef8ecc8e0207393"
)
HUNTER_ALL_SHA256 = (
    "cd401a6172218c42898adc5349ee28ace718aefb0b428b08c1c0b64ef7c14692"
)
HUNTER_OUTPUT_SHA256 = (
    "c142f6389a38549f3b2096b1d4b785f2c19b4062c8879eeeebd2ae04d7be5c7f"
)
REPLAY_OUTPUT_SHA256 = (
    "b70c6852da6e2edc31cc65f933f37a407f1cfd4f41c19572069239da8452cf2b"
)
G5_SOURCE_SHA256 = (
    "794111b992e912ec8471c8334a867d7b2db1d248f4b08f744f52faf7f50b86c3"
)
PIVOT_SOURCE_SHA256 = (
    "644104b0de90654466e75c6531109736b0445aadb357eee2413e8787ac3a53fa"
)
PIVOT_OUTPUT_SHA256 = (
    "0933c67a108b6d588e36737fb2b17b325ca36146976cfb035bebe036a6234036"
)
PIVOT_LEDGER_SHA256 = (
    "bec35518329b5d9e6ba2c9a8c87bfb20234a0c07dc1a5c5f2babec21888d452a"
)
PIVOT_LEDGER_SEMANTIC_SHA256 = (
    "ec878244b922ba5f48633614a86a1f9706c1fbdd0ebd6c61f020291cfd737bab"
)

SHARD_HASHES = (
    "3db3b6cd23f153d372cb0ab50c67c9e36241ee26a05adab798ad6c121e604144",
    "cf6497d63081795e1061815561a96c493e7d9da89b54b7a31dcc2e7b70cd4f1a",
    "420517a4296082b7b0ed631a2bb58d65e8e106d469811ebe3f5aae4c46c45b23",
    "05a849a174cb21c9b535917d0aee2de6fe608bdad9a2ccc10a780ba964e871e9",
    "ed098a7dca3277aa218d7551aa0c8b61cea7694e1e4a0837efaa57e50b518415",
    "d566eb5cd153bfd285d45a7930eed6b44cf2bbfa284f12bbc466331509951db6",
    "ace1736219c9c982006060c9ff3795b0445ffde3b583d97cc97c0a5a5280d46d",
    "fd80456ca425c6415e15ea9e2070e07ed342d68aeb958c954a4a37e10900d15c",
    "4847f7307e44666c75aa4f38a2556bb1b4de34bdaafd8558d8cbee2b54dc52b9",
    "5c46837c492818e10b2fc3a74051195f6c5c141c5c1b980fb9118cddc31cd540",
    "2139f43f41a2dc00204819789d692c808ef4a8469f11de30d1bf9a225dad6285",
    "39d39b4db24dd4fc0262a697739934595c80249b3828ac1cb027dc1a1a7ad5f7",
    "6ded1ef928a485f0ea9be2b42cacb0ebc20afd3146794b08b9934a3d5a439ef5",
    "5c84303d4dc9823816521252af72f9f4988e19a98d777d5d969144affc7b811b",
    "6a39f3fa31d62b1500f80b5363332fde315e65a69f93d7942bc74c69eaddabd4",
    "3f7b7bd52c6f100063d4bbc843aa95e2bd4c4298c0a2e7367c7072dda732c4c1",
    "bb028f81d973ad9ecc602575c1e79dd05c2ac38bf43604364adebdbdf423e870",
    "21cf2246efd8a7550c9f7e00303f8890e3c7d375cfae2845d76f3097437ac9d0",
    "cac225d5168b14281999a37ec54ed83e67b5644d23670ba8a6ab964be959893a",
    "ca9b63778c0c5ee444a3b00d546e59aa14d96b079502ea1f03d4469815bf343c",
    "455e6818abb3dc6e290ab7bdfff55b99091d801e5f05e5fdbf78b9d95504573b",
    "92bdf6728d2300c332224f629fe98acfed54332a47ee5f5d339ce2bf7d9302b2",
    "9db59ff1e03707a6b250b7df743b349845528f24961c4f6ab3842bf4f799b989",
    "7f965b38500c5f9a0e26b9ed146d498f68e03e0721170121b62ee2aec6e9c4d4",
    "13f3cbff4c03fbb00d23b56014d9b24fac8fed513e780433d15cd98c294d4041",
    "be0023225c9cf9e21590902fc31fa039b3be9a202405b0316e16550f53b202b9",
    "90891f23ec139baab5af42d36d09394ff8d3db8826c2d43d31d364b180a0a22b",
    "928f2515020608495e4ab213a1d845a961ad0cf71055ed86b102870fad129309",
    "d0d585157bd5220aacef1f3d00998279a3fb5f876862794b1500ba0cfceefd6e",
    "9ff29bdff1c1433f5ec79cebe2d65c275277801538a450fff250c17ed7975d03",
    "9487b8fab4c4f81f65290b43b4e6f0c8216671daeaa66c7d810f19df19772a42",
    "3177dcd74f11e6743d81660967b8059dc3f0e3f92b09255aaa825872c78aa0d5",
)
EXPECTED_SHARD_FILE_DIGEST = (
    "d686d1afd01b494ca6aa9a4733816c774f7e8e4ed398081f556c44e4af7869a0"
)
EXPECTED_SHARD_LEDGER_DIGEST = (
    "f02aea1d70f89b2c634b84b313da9bba78b26718e85fe1731dbb9ef496559a9b"
)
EXPECTED_REPLAY_CORPUS_DIGEST = (
    "7a30bf625c8b0f8f3d1789d051582be2272686312e4df4e30c49c0acc398bfa0"
)
EXPECTED_H1_ROOT_DIGEST = (
    "c5eb5a86bbe82c2dc46113527e67f6a1f50add87b100cad916ee6367d7177a13"
)
EXPECTED_HG_ROOT_DIGEST = (
    "01c7f66d905d432c45f1fa3a7700ffd0ba8b4127fec3dcfd0f9fe71e42ee71dc"
)
EXPECTED_HGP_ROOT_DIGEST = (
    "4c9a5f8f79b3a42849fb4e5edd16e3958f6eebc40f3a06f4312cd135dd787e0e"
)
EXPECTED_PROVED82_DIGEST = (
    "b1d4f1a63637f59f45da770c7b8dab830c6b455198a2e59e9913b10486c0158e"
)
EXPECTED_PROVED88_DIGEST = (
    "49dce548bc8d9aaf85408852bf5a4bfcce3ba3824e824c9884caa3d9aaa309f7"
)
EXPECTED_HG_ADDITIVE_DIGEST = (
    "55e5a9b8561cd64ee8ec946a07e6fb7977aaf0a8a2630c98f62069f6c760f404"
)
EXPECTED_PROVED180_DIGEST = (
    "75629a1862ca65c522679e17b3c35ce7bba3fd0ca412d099f4efcb82de84d7df"
)
EXPECTED_HGP_ADDITIVE_DIGEST = (
    "525e65b26807e1ea0309914c22ae7d853da0a6d443f449382f81291f77224da8"
)
EXPECTED_PROVED181_DIGEST = (
    "1921d07cab2653f6b610d9a6044d6bb4836f71d663e66d8bc9bd24b116c7b4f6"
)

EXPECTED_H1_BRANCH_DIGEST = (
    "4a5c0d887a71f5b4896fc6a3453c1861498616fb89902a273b890decbf48790a"
)
EXPECTED_G5_BRANCH_DIGEST = (
    "385b7c4638dfa3a4568c6649f3d7e784f3267f198754af45985acbc9946ab241"
)
EXPECTED_HG_BRANCH_UNION_DIGEST = (
    "8113e676f58539ca0f3d0ca615f97030c2cc3053d0b98c1a250ec4aae810cf30"
)
EXPECTED_PIVOT_BRANCH_DIGEST = (
    "2a47a0860d6c78f19a55118f677991e8ff9b0891aa08eadd2da2d38691e82d51"
)
EXPECTED_HGP_BRANCH_UNION_DIGEST = (
    "b2cff827a29c1257aa7f769fd51b8710abeb419947a5e97f40fe905f23e73e03"
)

S2 = F(99, 70)
MAX_CUTOFF = 15_000
RAW_PATH_READ_BYTES = Path.read_bytes


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    """Hash repository text independently of LF/CRLF checkout policy."""

    payload = RAW_PATH_READ_BYTES(path)
    require(
        b"\r" not in payload.replace(b"\r\n", b""),
        f"{path.name}: unexpected lone carriage return",
    )
    return hashlib.sha256(payload.replace(b"\r\n", b"\n")).hexdigest()


def lf_read_bytes(path: Path) -> bytes:
    """Read repository text on the canonical LF byte basis."""

    payload = RAW_PATH_READ_BYTES(path)
    require(
        b"\r" not in payload.replace(b"\r\n", b""),
        f"{path.name}: unexpected lone carriage return",
    )
    return payload.replace(b"\r\n", b"\n")


def line_after(text: str, prefix: str) -> str:
    rows = [
        line[len(prefix) :]
        for line in text.splitlines()
        if line.startswith(prefix)
    ]
    require(len(rows) == 1, f"expected one {prefix!r} line")
    return rows[0]


def body_digest(bodies: set[tuple[int, ...]]) -> str:
    payload = "".join(
        ",".join(map(str, body)) + "\n" for body in sorted(bodies)
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def branch_digest(
    keys: set[
        tuple[tuple[int, ...], int, int, int, tuple[int, ...]]
    ],
    header: str,
) -> str:
    digest = hashlib.sha256((header + "\n").encode())
    for body, gate_size, rank, apex, prefix in sorted(keys):
        digest.update(
            (
                f"{','.join(map(str, body))};K={gate_size};"
                f"rank={rank};a={apex};P={','.join(map(str, prefix))}\n"
            ).encode()
        )
    return digest.hexdigest()


def nearest_quantiles(
    histogram: Counter[int],
) -> tuple[tuple[int, int], ...]:
    values = []
    for value, count in histogram.items():
        values.extend([value] * count)
    values.sort()
    require(values, "empty H1 histogram")
    size = len(values)
    return tuple(
        (
            pct,
            values[0 if pct == 0 else math.ceil(pct * size / 100) - 1],
        )
        for pct in (0, 25, 50, 75, 90, 95, 99, 100)
    )


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    # Nested Hunter/G5 audit modules predate checkout-independent hashes.
    # Canonicalize their import-time reads, including transitive imports.
    original_read_bytes = Path.read_bytes
    Path.read_bytes = lf_read_bytes
    try:
        spec.loader.exec_module(module)
    finally:
        Path.read_bytes = original_read_bytes
    if hasattr(module, "file_sha256"):
        module.file_sha256 = file_sha256
    return module


def verify_replay_attestation() -> None:
    require(
        file_sha256(REPLAY_OUTPUT) == REPLAY_OUTPUT_SHA256,
        "ordinary/optimized replay attestation changed",
    )
    text = REPLAY_OUTPUT.read_text()
    require(
        line_after(text, "shards=") == "32;byte_matches=32;failures=0",
        "replay count changed",
    )
    require(
        line_after(text, "worker_schedule=")
        == "(shard0:4,shards1-31:1)",
        "worker schedule changed",
    )
    require(
        line_after(text, "shard_file_digest=")
        == EXPECTED_SHARD_FILE_DIGEST,
        "replay shard digest changed",
    )
    require(
        line_after(text, "ordinary_replay_corpus_digest=")
        == EXPECTED_REPLAY_CORPUS_DIGEST,
        "ordinary replay corpus changed",
    )
    require(text.endswith("all_exact_controls=PASS\n"), "replay did not pass")


def verify_hunter() -> tuple[tuple[tuple[int, ...], int, int, int], ...]:
    require(file_sha256(HUNTER_BASE) == HUNTER_BASE_SHA256, "Hunter base changed")
    require(file_sha256(HUNTER_ALL) == HUNTER_ALL_SHA256, "Hunter audit changed")
    require(
        file_sha256(HUNTER_OUTPUT) == HUNTER_OUTPUT_SHA256,
        "Hunter output changed",
    )
    hunter = load_module("thm2911_hunter_targets", HUNTER_ALL)
    source_targets = tuple(hunter.TARGETS)
    require(
        len(source_targets) == len(set(source_targets)) == 24,
        "Hunter source target battery changed",
    )
    text = HUNTER_OUTPUT.read_text()
    require(
        line_after(text, "targets=")
        == "24;statuses=(('HUNTER_CLOSED', 24),)",
        "Hunter target status changed",
    )
    rows = [line for line in text.splitlines() if line.startswith("E=")]
    require(len(rows) == 24, "Hunter row count changed")
    hostile = 0
    repairs = 0
    hard = 0
    margins = []
    output_records = []
    for line in rows:
        body_match = re.search(r"^E=(\([^;]+\))", line)
        rank_match = re.search(r";rank=(\d+)", line)
        core_match = re.search(r";H=(\d+)", line)
        pivot_match = re.search(r";pivot=(\d+)", line)
        hostile_match = re.search(r";star_hostile=(\d+)", line)
        repairs_match = re.search(r";hunter_repairs=(\d+)", line)
        hard_match = re.search(r";hunter_hard=(\d+)", line)
        margin_match = re.search(r";margin=(-?\d+/\d+)", line)
        require(
            all(
                match is not None
                for match in (
                    body_match,
                    rank_match,
                    core_match,
                    pivot_match,
                    hostile_match,
                    repairs_match,
                    hard_match,
                    margin_match,
                )
            ),
            "Hunter row schema changed",
        )
        hostile += int(hostile_match.group(1))
        repairs += int(repairs_match.group(1))
        hard += int(hard_match.group(1))
        margins.append(F(margin_match.group(1)))
        output_records.append(
            (
                ast.literal_eval(body_match.group(1)),
                int(rank_match.group(1)),
                int(pivot_match.group(1)),
                int(core_match.group(1)),
            )
        )
    require(
        tuple(record[:3] for record in output_records) == source_targets,
        "Hunter output identities differ from pinned source targets",
    )
    require(
        (hostile, repairs, hard) == (54, 54, 0),
        "Hunter repair counts changed",
    )
    require(min(margins) == F(1_799_771, 75_716_368), "Hunter margin changed")
    require(
        line_after(text, "aggregate_sha256=")
        == "030c5af7cdb8fe5c2c8188b3debb78e8e0aabb79a7e553308d0f11a411ded74b",
        "Hunter semantic digest changed",
    )
    require(text.endswith("all_exact_controls=PASS\n"), "Hunter did not pass")
    return tuple(output_records)


def verify_shards() -> tuple[
    set[tuple[int, ...]],
    tuple[tuple[tuple[int, ...], int, int, int], ...],
]:
    require(file_sha256(SCOUT) == SCOUT_SHA256, "H1 scout changed")
    require(len(SHARD_HASHES) == 32, "shard hash battery changed")
    counts = [0] * 9
    statuses: Counter[str] = Counter()
    h_histogram: Counter[int] = Counter()
    closures: list[tuple[int, ...]] = []
    exceptions = []
    ledger_hashes = []
    margins = []
    file_digest = hashlib.sha256(b"LRC14/j6/H1-shards/v1\n")

    for index, expected_hash in enumerate(SHARD_HASHES):
        path = RESULTS / f"all_s{index:02d}_of32.out"
        data = lf_read_bytes(path)
        actual_hash = hashlib.sha256(data).hexdigest()
        require(actual_hash == expected_hash, f"shard {index} hash changed")
        text = data.decode()
        require(text.endswith("all_exact_controls=PASS\n"), f"shard {index} open")
        expected_workers = 4 if index == 0 else 1
        require(
            (
                f"parameters=workers:{expected_workers},"
                "max_cutoff:15000,max_combinations:50000000000,"
                f"scope:all,shard:{index}/32"
            )
            in text,
            f"shard {index} parameters changed",
        )
        file_digest.update(f"{index:02d};{actual_hash}\n".encode())
        row_counts = ast.literal_eval(line_after(text, "counts="))
        require(len(row_counts) == len(counts), f"shard {index} counts changed")
        counts = [left + right for left, right in zip(counts, row_counts)]
        statuses.update(dict(ast.literal_eval(line_after(text, "statuses="))))
        h_histogram.update(
            dict(ast.literal_eval(line_after(text, "H1_size_histogram=")))
        )
        closures.extend(ast.literal_eval(line_after(text, "root_closures=")))
        ledger_hashes.append(line_after(text, "ledger_sha256="))
        margins.append(F(line_after(text, "depth1_min_margin=")))

        active = False
        for line in text.splitlines():
            if line == "DEPTH1_OPEN":
                active = True
                continue
            if line == "ROOT_CLOSURES":
                active = False
            if not active or not line.startswith("E="):
                continue
            body_match = re.search(r"^E=(\([^;]+\))", line)
            rank_match = re.search(r";rank=(\d+)", line)
            core_match = re.search(r";H=(\d+)", line)
            open_match = re.search(r";d1open=(.*?);d1digest=", line)
            require(
                all(
                    match is not None
                    for match in (
                        body_match,
                        rank_match,
                        core_match,
                        open_match,
                    )
                ),
                f"shard {index} exception row changed",
            )
            open_text = open_match.group(1)
            open_indices = tuple(
                map(int, re.findall(r"\((\d+), \d+, Fraction", open_text))
            )
            open_pivots = tuple(
                map(int, re.findall(r"\(\d+, (\d+), Fraction", open_text))
            )
            exceptions.append(
                (
                    ast.literal_eval(body_match.group(1)),
                    int(rank_match.group(1)),
                    int(core_match.group(1)),
                    open_indices,
                    open_pivots,
                )
            )

    ledger_digest = hashlib.sha256(b"LRC14/j6/H1-shard-ledgers/v1\n")
    for index, value in enumerate(ledger_hashes):
        ledger_digest.update(f"{index:02d};{value}\n".encode())
    closure_set = set(closures)

    require(
        tuple(counts) == (3432, 41_415, 14_806, 6180, 5999, 5999, 5999, 0, 132),
        "aggregate shard counts changed",
    )
    require(
        tuple(sorted(statuses.items()))
        == (("CLOSED", 24), ("CUTOFF_SKIP", 181), ("DEPTH1_CLOSED", 5975)),
        "aggregate shard statuses changed",
    )
    require(
        nearest_quantiles(h_histogram)
        == (
            (0, 5),
            (25, 8),
            (50, 16),
            (75, 36),
            (90, 83),
            (95, 141),
            (99, 283),
            (100, 444),
        ),
        "H1 quantiles changed",
    )
    require(
        len(exceptions) == 24
        and all(row[3] == (0,) and len(row[4]) == 1 for row in exceptions),
        "star-exception anatomy changed",
    )
    require(min(margins) == F(1, 3_166_020), "ordered-pivot margin changed")
    require(
        file_digest.hexdigest() == EXPECTED_SHARD_FILE_DIGEST,
        "aggregate shard-file digest changed",
    )
    require(
        ledger_digest.hexdigest() == EXPECTED_SHARD_LEDGER_DIGEST,
        "aggregate shard-ledger digest changed",
    )
    require(len(closures) == len(closure_set) == 132, "closure roots changed")
    require(body_digest(closure_set) == EXPECTED_H1_ROOT_DIGEST, "H1 root digest")
    exception_targets = tuple(
        (body, rank, pivots[0], core_size)
        for body, rank, core_size, _indices, pivots in exceptions
    )
    return closure_set, exception_targets


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def h1_cutoff(row: dict[str, object]) -> int | None:
    r4 = sum(row["qs"][:4], F(0))
    epsilon = row["mass"] - r4 - row["mass"] / 7
    if epsilon <= 0:
        return None
    return ceiling(S2 * row["components"] / (7 * epsilon)) - 1


def branch_key(
    row: dict[str, object],
) -> tuple[tuple[int, ...], int, int, int, tuple[int, ...]]:
    return (
        row["body"],
        row["gate_size"],
        row["rank"],
        row["apex"],
        row["prefix"],
    )


def verify_pivot_ledger(
    rows: list[dict[str, object]],
) -> tuple[
    set[tuple[tuple[int, ...], int, int, int, tuple[int, ...]]],
    set[tuple[int, ...]],
    set[tuple[int, ...]],
]:
    require(
        file_sha256(PIVOT_SOURCE) == PIVOT_SOURCE_SHA256,
        "THM-2904 source changed",
    )
    require(
        file_sha256(PIVOT_OUTPUT) == PIVOT_OUTPUT_SHA256,
        "THM-2904 output changed",
    )
    require(
        file_sha256(PIVOT_LEDGER) == PIVOT_LEDGER_SHA256,
        "THM-2904 promoted ledger changed",
    )

    output = PIVOT_OUTPUT.read_text()
    require(
        output.endswith("all_exact_controls=PASS\n"),
        "THM-2904 output did not pass",
    )
    counts = ast.literal_eval(line_after(output, "counts="))
    require(
        (
            counts[0],
            counts[1],
            counts[11],
            counts[12],
            counts[13],
            counts[14],
            counts[15],
            counts[16],
            counts[17],
            counts[18],
            counts[19],
            counts[20],
            counts[21],
            counts[22],
            counts[23],
            counts[24],
            counts[25],
        )
        == (
            11_842,
            52,
            55_293,
            4_071,
            51_222,
            4,
            4,
            0,
            279,
            0,
            11_563,
            3_411,
            10,
            4,
            6,
            88,
            3_344,
        ),
        "THM-2904 salient counts changed",
    )
    require(line_after(output, "mode=") == "LOCKED", "THM-2904 mode changed")
    require(
        line_after(output, "ledger_sha256=")
        == PIVOT_LEDGER_SEMANTIC_SHA256,
        "THM-2904 output ledger digest changed",
    )
    output_closed_roots = set(
        ast.literal_eval(line_after(output, "closed_roots="))
    )
    output_additive_roots = set(
        ast.literal_eval(line_after(output, "additive_roots="))
    )

    lookup: dict[
        tuple[tuple[int, ...], int, int, tuple[int, ...]],
        tuple[tuple[int, ...], int, int, int, tuple[int, ...]],
    ] = {}
    for row in rows:
        short_key = (
            row["body"],
            row["rank"],
            row["apex"],
            row["prefix"],
        )
        require(short_key not in lookup, "hard rows collide after dropping K")
        lookup[short_key] = branch_key(row)
    require(len(lookup) == len(rows) == 14_806, "hard lookup changed")

    ledger_lines = PIVOT_LEDGER.read_text().splitlines(keepends=True)
    require(
        ledger_lines[0]
        == "LRC14 j6 all-hard ranked H1 Hunter-pivot ledger\n",
        "THM-2904 ledger header changed",
    )
    semantic = hashlib.sha256(
        b"LRC14/j6/all-hard/ranked-H1-Hunter-pivots/v2\n"
    )
    seen_short_keys = set()
    pivot_keys = set()
    for line in ledger_lines:
        if not line.startswith("H1;"):
            continue
        semantic_line = line.removeprefix("H1;")
        semantic.update(semantic_line.encode())
        fields = {
            item.split("=", 1)[0]: item.split("=", 1)[1]
            for item in semantic_line.rstrip("\n").split(";")
            if "=" in item
        }
        short_key = (
            tuple(map(int, fields["E"].split(","))),
            int(fields["rank"]),
            int(fields["a"]),
            tuple(map(int, fields["P"].split(","))),
        )
        require(
            short_key in lookup and short_key not in seen_short_keys,
            "THM-2904 ledger/hard-row join changed",
        )
        seen_short_keys.add(short_key)
        pivots = fields["pivot"].split(",") if fields["pivot"] else []
        require(pivots, "THM-2904 hostile-centre core became empty")
        closed_flags = []
        for pivot in pivots:
            parts = pivot.split(":")
            require(len(parts) == 5, "THM-2904 pivot schema changed")
            closed, equality_repair = map(int, parts[-2:])
            require(
                closed in (0, 1)
                and equality_repair in (0, 1)
                and equality_repair <= closed,
                "THM-2904 pivot flag changed",
            )
            closed_flags.append(bool(closed))
        if all(closed_flags):
            pivot_keys.add(lookup[short_key])

    require(
        len(seen_short_keys) == 11_842 and len(pivot_keys) == 279,
        "THM-2904 ledger row/closure count changed",
    )
    require(
        semantic.hexdigest() == PIVOT_LEDGER_SEMANTIC_SHA256,
        "THM-2904 semantic ledger changed",
    )
    require(
        line_after(PIVOT_LEDGER.read_text(), "ledger_sha256=")
        == PIVOT_LEDGER_SEMANTIC_SHA256,
        "THM-2904 promoted-ledger footer changed",
    )

    g5_survivor_short_keys = {
        (
            row["body"],
            row["rank"],
            row["apex"],
            row["prefix"],
        )
        for row in rows
        if row["margin"] <= 0
    }
    require(
        seen_short_keys == g5_survivor_short_keys,
        "THM-2904 ledger is not the exact G5-survivor universe",
    )
    survivor_by_body: dict[
        tuple[int, ...],
        list[tuple[tuple[int, ...], int, int, int, tuple[int, ...]]],
    ] = defaultdict(list)
    for short_key in seen_short_keys:
        survivor_by_body[short_key[0]].append(lookup[short_key])
    recomposed_closed_roots = {
        body
        for body, body_keys in survivor_by_body.items()
        if all(key in pivot_keys for key in body_keys)
    }
    require(
        recomposed_closed_roots == output_closed_roots
        and len(output_additive_roots) == 6,
        "THM-2904 whole-root recomposition changed",
    )
    return pivot_keys, output_closed_roots, output_additive_roots


def verify_join(
    shard_roots: set[tuple[int, ...]],
) -> dict[str, object]:
    require(file_sha256(G5_SOURCE) == G5_SOURCE_SHA256, "THM-2905 changed")
    g5 = load_module("thm2911_three_route_join", G5_SOURCE)
    thm2895, thm2898, thm2899, thm2901, thm2902 = (
        g5.canonical_root_sets()
    )
    one_hard_controls = g5.thm2903_controls()
    rows = g5.parse_pair_rows(g5.parse_hard_rows())
    pivot_keys, pivot_closed_roots, pivot_additive_roots = (
        verify_pivot_ledger(rows)
    )
    by_body: dict[tuple[int, ...], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_body[row["body"]].append(row)

    h1_keys = {
        branch_key(row)
        for row in rows
        if (
            (cutoff := h1_cutoff(row)) is not None
            and cutoff <= MAX_CUTOFF
        )
    }
    g5_keys = {
        branch_key(row)
        for row in rows
        if row["margin"] > 0
    }
    exception_keys = {
        branch_key(row) for row in rows if row["pair_margin"] <= 0
    }

    def hard_roots(
        keys: set[
            tuple[tuple[int, ...], int, int, int, tuple[int, ...]]
        ],
    ) -> set[tuple[int, ...]]:
        return {
            body
            for body, body_rows in by_body.items()
            if all(branch_key(row) in keys for row in body_rows)
        }

    scalar_roots = set(thm2899)
    h1_hard_roots = hard_roots(h1_keys)
    g5_hard_roots = hard_roots(g5_keys)
    hg_hard_roots = hard_roots(h1_keys | g5_keys)
    hgp_hard_roots = hard_roots(h1_keys | g5_keys | pivot_keys)
    h1_roots = h1_hard_roots | scalar_roots
    g5_roots = g5_hard_roots | scalar_roots
    hg_roots = hg_hard_roots | scalar_roots
    hgp_roots = hgp_hard_roots | scalar_roots
    require(h1_roots == shard_roots, "ledger H1 roots differ from shard roots")

    prior_fifteen = thm2895 | thm2898 | thm2899 | thm2901
    one_hard = {
        body
        for body, body_rows in by_body.items()
        if len(body_rows) == 1 and body_rows[0]["direct_margin"] <= 0
    }
    through_2903 = prior_fifteen | thm2902 | one_hard
    through_2905 = through_2903 | g5_hard_roots
    require(
        len(one_hard) == one_hard_controls["closed_roots"],
        "THM-2903 root bank changed",
    )
    require(
        pivot_additive_roots == pivot_closed_roots - through_2905,
        "THM-2904 additive-root difference changed",
    )
    through_2904 = through_2905 | pivot_additive_roots

    atoms = Counter()
    for key in exception_keys | h1_keys | g5_keys | pivot_keys:
        atom = "".join(
            label
            for label, key_set in (
                ("E", exception_keys),
                ("G", g5_keys),
                ("H", h1_keys),
                ("P", pivot_keys),
            )
            if key in key_set
        )
        atoms[atom] += 1
    require(
        tuple(sorted(atoms.items()))
        == (
            ("E", 52),
            ("G", 55),
            ("GH", 2_909),
            ("H", 2_875),
            ("HP", 215),
            ("P", 64),
        ),
        "four-way branch atoms changed",
    )
    require(
        (
            len(rows),
            len(h1_keys),
            len(g5_keys),
            len(pivot_keys),
            len(h1_keys & g5_keys),
            len(h1_keys & pivot_keys),
            len(g5_keys & pivot_keys),
            len(h1_keys & g5_keys & pivot_keys),
            len(h1_keys | g5_keys),
            len(h1_keys | g5_keys | pivot_keys),
        )
        == (14_806, 5_999, 2_964, 279, 2_909, 215, 0, 0, 6_054, 6_118),
        "three-route branchwise join counts changed",
    )
    require(
        (
            len(h1_hard_roots),
            len(h1_roots),
            len(g5_hard_roots),
            len(g5_roots),
            len(hg_hard_roots),
            len(hg_roots),
            len(hgp_hard_roots),
            len(hgp_roots),
        )
        == (127, 132, 16, 21, 129, 134, 133, 138),
        "three-route whole-root counts changed",
    )
    require(
        hg_roots - h1_roots
        == {
            (1, 2, 3, 4, 5, 8, 9),
            (5, 7, 8, 10, 11, 13, 14),
        },
        "two cross-route roots changed",
    )
    require(
        (
            len(through_2903),
            len(through_2905),
            len(pivot_closed_roots),
            len(pivot_additive_roots),
            len(through_2904),
        )
        == (76, 82, 10, 6, 88),
        "proved baselines changed",
    )
    require(
        (
            len(h1_roots & through_2904),
            len(h1_roots - through_2904),
            len(h1_roots | through_2904),
            len(hg_roots & through_2904),
            len(hg_roots - through_2904),
            len(hg_roots | through_2904),
            len(hgp_roots & through_2904),
            len(hgp_roots - through_2904),
            len(hgp_roots | through_2904),
            3_432 - len(hgp_roots | through_2904),
        )
        == (40, 92, 180, 42, 92, 180, 45, 93, 181, 3_251),
        "current-root composition changed",
    )
    require(
        h1_roots - through_2904 == hg_roots - through_2904
        and h1_roots | through_2904 == hg_roots | through_2904,
        "G5 changed the current finite-H1 root union",
    )
    new_beyond_hg_current = (
        hgp_roots | through_2904
    ) - (hg_roots | through_2904)
    require(
        new_beyond_hg_current == {(1, 2, 3, 5, 6, 8, 13)},
        "unique H1/pivot composition root changed",
    )

    mixed = []
    mixed_details = []
    for body in sorted(hg_hard_roots):
        classes = []
        for row in by_body[body]:
            row_key = branch_key(row)
            h1 = row_key in h1_keys
            star = row_key in g5_keys
            classes.append(
                (
                    "both" if h1 and star else "H1" if h1 else "G5",
                    row["rank"],
                    row["apex"],
                    row["prefix"],
                )
            )
            if body == (1, 2, 3, 4, 5, 8, 9):
                epsilon = (
                    row["mass"]
                    - sum(row["qs"][:4], F(0))
                    - row["mass"] / 7
                )
                mixed_details.append(
                    (
                        row["rank"],
                        row["apex"],
                        row["prefix"],
                        epsilon,
                        h1_cutoff(row),
                        row["margin"],
                    )
                )
        if any(item[0] == "H1" for item in classes) and any(
            item[0] == "G5" for item in classes
        ):
            mixed.append((body, tuple(classes)))
    require(
        mixed
        == [
            (
                (1, 2, 3, 4, 5, 8, 9),
                (
                    ("G5", 1, 22, (22,)),
                    ("H1", 2, 33, (22, 33)),
                ),
            )
        ],
        "genuinely mixed root changed",
    )
    require(
        tuple(mixed_details)
        == (
            (1, 22, (22,), F(-193, 2_522_520), None, F(43, 84_084)),
            (
                2,
                33,
                (22, 33),
                F(26_111, 2_522_520),
                585,
                F(-1, 10_010),
            ),
        ),
        "mixed branch margins changed",
    )

    pivot_only_roots = hard_roots(pivot_keys) | scalar_roots
    hp_roots = hard_roots(h1_keys | pivot_keys) | scalar_roots
    gp_roots = hard_roots(g5_keys | pivot_keys) | scalar_roots
    pair_synergies = {
        "HG": hg_roots - (h1_roots | g5_roots),
        "HP": hp_roots - (h1_roots | pivot_only_roots),
        "GP": gp_roots - (g5_roots | pivot_only_roots),
    }
    require(
        (
            len(pivot_only_roots),
            len(hp_roots),
            len(gp_roots),
            tuple(map(len, pair_synergies.values())),
            hgp_roots - (hg_roots | hp_roots | gp_roots),
        )
        == (9, 136, 31, (1, 2, 6), set()),
        "pair-route synergy anatomy changed",
    )
    require(
        pair_synergies["GP"] == pivot_additive_roots,
        "THM-2904 six-root gain is no longer the G5/pivot synergy",
    )

    new_root_details = []
    for row in by_body[(1, 2, 3, 5, 6, 8, 13)]:
        row_key = branch_key(row)
        route = "".join(
            label
            for label, keys in (
                ("H", h1_keys),
                ("G", g5_keys),
                ("P", pivot_keys),
            )
            if row_key in keys
        )
        new_root_details.append(
            (
                route,
                row["gate_size"],
                row["rank"],
                row["apex"],
                row["prefix"],
                h1_cutoff(row),
                row["margin"],
            )
        )
    require(
        tuple(new_root_details)
        == (
            (
                "H",
                8,
                1,
                22,
                (22,),
                6_930,
                F(-174_401, 12_612_600),
            ),
            (
                "P",
                8,
                2,
                21,
                (22, 21),
                None,
                F(-383, 254_800),
            ),
        ),
        "new H1/pivot root anatomy changed",
    )

    h1_branch_digest = branch_digest(
        h1_keys,
        "LRC14/j6/r4-H1-finite-prefix-branch-keys/v1",
    )
    g5_branch_digest = branch_digest(
        g5_keys,
        "LRC14/j6/G5-positive-prefix-branch-keys/v1",
    )
    hg_branch_digest = branch_digest(
        h1_keys | g5_keys,
        "LRC14/j6/G5-or-finite-H1-prefix-branch-keys/v1",
    )
    pivot_branch_digest = branch_digest(
        pivot_keys,
        "LRC14/j6/hostile-centre-pivot-prefix-branch-keys/v1",
    )
    hgp_branch_digest = branch_digest(
        h1_keys | g5_keys | pivot_keys,
        "LRC14/j6/G5-or-finite-H1-or-hostile-centre-pivot-prefix-branch-keys/v1",
    )
    require(
        h1_branch_digest == EXPECTED_H1_BRANCH_DIGEST,
        "H1 prefix-branch digest changed",
    )
    require(
        g5_branch_digest == EXPECTED_G5_BRANCH_DIGEST,
        "G5 prefix-branch digest changed",
    )
    require(
        hg_branch_digest == EXPECTED_HG_BRANCH_UNION_DIGEST,
        "H1/G5 prefix-branch digest changed",
    )
    require(
        pivot_branch_digest == EXPECTED_PIVOT_BRANCH_DIGEST,
        "hostile-centre pivot prefix-branch digest changed",
    )
    require(
        hgp_branch_digest == EXPECTED_HGP_BRANCH_UNION_DIGEST,
        "three-route prefix-branch digest changed",
    )

    require(body_digest(h1_roots) == EXPECTED_H1_ROOT_DIGEST, "H1 root set")
    require(body_digest(hg_roots) == EXPECTED_HG_ROOT_DIGEST, "H1/G5 roots")
    require(body_digest(hgp_roots) == EXPECTED_HGP_ROOT_DIGEST, "H1/G5/P roots")
    require(
        body_digest(through_2905) == EXPECTED_PROVED82_DIGEST,
        "proved-82 baseline",
    )
    require(
        body_digest(through_2904) == EXPECTED_PROVED88_DIGEST,
        "proved-88 baseline",
    )
    require(
        body_digest(hg_roots - through_2904)
        == EXPECTED_HG_ADDITIVE_DIGEST,
        "H1/G5 additions against proved 88",
    )
    require(
        body_digest(hg_roots | through_2904)
        == EXPECTED_PROVED180_DIGEST,
        "H1/G5 proved union",
    )
    require(
        body_digest(hgp_roots - through_2904)
        == EXPECTED_HGP_ADDITIVE_DIGEST,
        "three-route additions against proved 88",
    )
    require(
        body_digest(hgp_roots | through_2904)
        == EXPECTED_PROVED181_DIGEST,
        "final proved-181 union",
    )
    return {
        "h1_branch_digest": h1_branch_digest,
        "g5_branch_digest": g5_branch_digest,
        "hg_branch_digest": hg_branch_digest,
        "pivot_branch_digest": pivot_branch_digest,
        "hgp_branch_digest": hgp_branch_digest,
        "atoms": tuple(sorted(atoms.items())),
        "mixed": tuple(mixed),
        "mixed_details": tuple(mixed_details),
        "new_root_details": tuple(new_root_details),
        "pair_synergy_counts": tuple(map(len, pair_synergies.values())),
        "h1_roots": tuple(sorted(h1_roots)),
        "hgp_roots": tuple(sorted(hgp_roots)),
        "additive_roots": tuple(sorted(hgp_roots - through_2904)),
    }


def main() -> None:
    verify_replay_attestation()
    hunter_targets = verify_hunter()
    shard_roots, shard_exception_targets = verify_shards()
    require(
        shard_exception_targets == hunter_targets,
        "shard exception identities differ from Hunter target battery",
    )
    joined = verify_join(shard_roots)

    print("LRC14 j6 finite ranked-H1/G5/hostile-pivot composition")
    print(
        "finite_H1=(eligible=6180,cutoff_le_15000=5999,"
        "ordered_pivot=5975,Hunter=24,deferred=181,witnesses=0)"
    )
    print(
        "Hunter=(targets=24,hostile_five_sets=54,repairs=54,"
        "hard=0,target_identity_join=24/24,"
        "min_margin=1799771/75716368)"
    )
    print(
        "branches=(hard=14806,H1=5999,G5=2964,pivot=279,"
        "HG_union=6054,HGP_union=6118)"
    )
    print(f"branch_atoms={joined['atoms']}")
    print(
        "roots=(H1=132,G5_hard=16,G5_with_scalar=21,"
        "HG=134,HGP=138,pure_three_way_synergy=0)"
    )
    print(
        "current_composition=(baseline=88,H1_intersection=40,"
        "H1_additive=92,H1_union=180,HG_intersection=42,"
        "HG_additive=92,HG_union=180,HGP_intersection=45,"
        "HGP_additive=93,proved_union=181,residual=3251)"
    )
    print(f"genuinely_mixed={joined['mixed']}")
    print(f"mixed_branch_details={joined['mixed_details']}")
    print(f"new_root_details={joined['new_root_details']}")
    print(f"pair_synergy_counts_HG_HP_GP={joined['pair_synergy_counts']}")
    print(f"H1_root_digest={EXPECTED_H1_ROOT_DIGEST}")
    print(f"HG_root_digest={EXPECTED_HG_ROOT_DIGEST}")
    print(f"HGP_root_digest={EXPECTED_HGP_ROOT_DIGEST}")
    print(f"proved82_digest={EXPECTED_PROVED82_DIGEST}")
    print(f"proved88_digest={EXPECTED_PROVED88_DIGEST}")
    print(f"HG_additive_digest={EXPECTED_HG_ADDITIVE_DIGEST}")
    print(f"proved180_digest={EXPECTED_PROVED180_DIGEST}")
    print(f"HGP_additive_digest={EXPECTED_HGP_ADDITIVE_DIGEST}")
    print(f"proved181_digest={EXPECTED_PROVED181_DIGEST}")
    print(f"H1_prefix_branch_digest={joined['h1_branch_digest']}")
    print(f"G5_prefix_branch_digest={joined['g5_branch_digest']}")
    print(f"HG_prefix_branch_digest={joined['hg_branch_digest']}")
    print(f"pivot_prefix_branch_digest={joined['pivot_branch_digest']}")
    print(f"HGP_prefix_branch_digest={joined['hgp_branch_digest']}")
    print(f"THM2904_ledger_file_digest={PIVOT_LEDGER_SHA256}")
    print(f"THM2904_ledger_semantic_digest={PIVOT_LEDGER_SEMANTIC_SHA256}")
    print(f"shard_file_digest={EXPECTED_SHARD_FILE_DIGEST}")
    print(f"shard_ledger_digest={EXPECTED_SHARD_LEDGER_DIGEST}")
    print(f"ordinary_replay_corpus_digest={EXPECTED_REPLAY_CORPUS_DIGEST}")
    print("EXPLICIT_ROOT_LISTS")
    for body in joined["h1_roots"]:
        print("H1_ROOT=" + ",".join(map(str, body)))
    for body in joined["hgp_roots"]:
        print("THREE_ROUTE_ROOT=" + ",".join(map(str, body)))
    for body in joined["additive_roots"]:
        print("ADDITIVE_ROOT=" + ",".join(map(str, body)))
    print("mode=LOCKED_FINITE_EXACT")
    print(
        "scope=5999 finite r4/H1,2964 G5,279 hostile-centre pivot "
        "hard-branch certificates;181 H1 cutoff tails deferred;not LRC14"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
