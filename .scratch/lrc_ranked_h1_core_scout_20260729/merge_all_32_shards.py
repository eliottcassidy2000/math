#!/usr/bin/env python3
"""Locked merge of the 32 all-root ranked-H1 finite-core shards."""

from __future__ import annotations

import ast
import hashlib
import math
import re
from collections import Counter
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
SCOUT = HERE / "scout.py"
SCOUT_SHA = "6abc06972cda64bb4f53db26cca62bbea5362ca178bdbf5bf7398e8b0f28317a"
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
EXPECTED_COUNTS = (3432, 41_415, 14_806, 6180, 5999, 5999, 5999, 0, 132)
EXPECTED_STATUSES = (
    ("CLOSED", 24),
    ("CUTOFF_SKIP", 181),
    ("DEPTH1_CLOSED", 5975),
)
EXPECTED_H_QUANTILES = (
    (0, 5),
    (25, 8),
    (50, 16),
    (75, 36),
    (90, 83),
    (95, 141),
    (99, 283),
    (100, 444),
)
EXPECTED_SHARD_FILE_DIGEST = (
    "d686d1afd01b494ca6aa9a4733816c774f7e8e4ed398081f556c44e4af7869a0"
)
EXPECTED_SHARD_LEDGER_DIGEST = (
    "f02aea1d70f89b2c634b84b313da9bba78b26718e85fe1731dbb9ef496559a9b"
)
EXPECTED_CLOSURE_DIGEST = (
    "c5eb5a86bbe82c2dc46113527e67f6a1f50add87b100cad916ee6367d7177a13"
)
EXPECTED_NEW_DIGEST = (
    "3b133f35fa91191819b35a2760dcdf822cf8523818776584b02de36939658b47"
)
EXPECTED_UNION_DIGEST = (
    "8536e577993aa395a6d2e676ea2c6ea0b29b99d1d15787dee4c7480fa0eebffa"
)

PROVED_17 = frozenset(
    {
        # THM-2899 scalar terminal.
        (1, 2, 3, 4, 5, 6, 13),
        (1, 2, 3, 4, 6, 7, 14),
        (1, 2, 3, 4, 6, 11, 13),
        (1, 2, 3, 4, 6, 12, 13),
        (7, 8, 9, 11, 12, 13, 14),
        # THM-2895 parity descent.
        (2, 8, 9, 10, 11, 13, 14),
        (1, 3, 9, 10, 11, 12, 14),
        (2, 5, 9, 11, 12, 13, 14),
        (2, 3, 4, 5, 6, 7, 8),
        # THM-2898 unique K=25 closure.
        (1, 8, 10, 11, 12, 13, 14),
        # THM-2901 exact pair-cap terminal.
        (1, 2, 3, 4, 6, 11, 12),
        (1, 2, 3, 5, 6, 10, 13),
        (1, 2, 4, 5, 6, 12, 13),
        (1, 3, 4, 5, 6, 7, 14),
        (5, 7, 8, 10, 11, 13, 14),
        # THM-2902 omission-six H1 terminal.
        (1, 2, 3, 4, 5, 10, 12),
        (1, 2, 3, 4, 5, 12, 13),
    }
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def line_after(text: str, prefix: str) -> str:
    rows = [line[len(prefix) :] for line in text.splitlines()
            if line.startswith(prefix)]
    require(len(rows) == 1, f"expected one {prefix!r} line")
    return rows[0]


def body_digest(bodies: tuple[tuple[int, ...], ...]) -> str:
    payload = "\n".join(",".join(map(str, body)) for body in bodies) + "\n"
    return hashlib.sha256(payload.encode()).hexdigest()


def quantiles(histogram: Counter[int]) -> tuple[tuple[int, int], ...]:
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


def main() -> None:
    require(
        hashlib.sha256(SCOUT.read_bytes()).hexdigest() == SCOUT_SHA,
        "H1 scout source changed",
    )
    require(len(SHARD_HASHES) == 32, "shard hash battery changed")
    counts = [0] * len(EXPECTED_COUNTS)
    statuses: Counter[str] = Counter()
    h_histogram: Counter[int] = Counter()
    closures: list[tuple[int, ...]] = []
    exceptions = []
    ledger_hashes = []
    minimum_margins = []
    file_digest = hashlib.sha256(b"LRC14/j6/H1-shards/v1\n")

    for index, expected_hash in enumerate(SHARD_HASHES):
        path = HERE / f"all_s{index:02d}_of32.out"
        data = path.read_bytes()
        actual_hash = hashlib.sha256(data).hexdigest()
        require(actual_hash == expected_hash, f"shard {index} hash changed")
        text = data.decode()
        require(
            text.endswith("all_exact_controls=PASS\n"),
            f"shard {index} did not finish",
        )
        require(
            f"scope:all,shard:{index}/32" in text,
            f"shard {index} parameter line changed",
        )
        file_digest.update(
            f"{index:02d};{actual_hash}\n".encode()
        )
        row_counts = ast.literal_eval(line_after(text, "counts="))
        require(
            isinstance(row_counts, tuple)
            and len(row_counts) == len(counts),
            f"shard {index} count row changed",
        )
        counts = [x + y for x, y in zip(counts, row_counts)]
        statuses.update(dict(ast.literal_eval(line_after(text, "statuses="))))
        h_histogram.update(
            dict(ast.literal_eval(line_after(text, "H1_size_histogram=")))
        )
        shard_closures = ast.literal_eval(
            line_after(text, "root_closures=")
        )
        closures.extend(shard_closures)
        ledger_hashes.append(line_after(text, "ledger_sha256="))
        minimum_margins.append(
            F(line_after(text, "depth1_min_margin="))
        )

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
                map(
                    int,
                    re.findall(r"\((\d+), \d+, Fraction", open_text),
                )
            )
            open_pivots = tuple(
                map(
                    int,
                    re.findall(r"\(\d+, (\d+), Fraction", open_text),
                )
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

    file_digest_value = file_digest.hexdigest()
    ledger_digest = hashlib.sha256(
        b"LRC14/j6/H1-shard-ledgers/v1\n"
    )
    for index, value in enumerate(ledger_hashes):
        ledger_digest.update(f"{index:02d};{value}\n".encode())
    ledger_digest_value = ledger_digest.hexdigest()
    closure_tuple = tuple(sorted(closures))
    require(len(closure_tuple) == len(set(closure_tuple)), "duplicate root")
    new_tuple = tuple(sorted(set(closure_tuple) - PROVED_17))
    union_tuple = tuple(sorted(set(closure_tuple) | PROVED_17))

    require(tuple(counts) == EXPECTED_COUNTS, "aggregate counts changed")
    require(
        tuple(sorted(statuses.items())) == EXPECTED_STATUSES,
        "aggregate statuses changed",
    )
    require(sum(h_histogram.values()) == 5999, "H1 histogram mass changed")
    require(quantiles(h_histogram) == EXPECTED_H_QUANTILES, "H1 quantiles")
    require(
        len(exceptions) == 24
        and all(row[3] == (0,) and len(row[4]) == 1 for row in exceptions),
        "star-exception anatomy changed",
    )
    require(min(minimum_margins) == F(1, 3_166_020), "minimum margin")
    require(
        file_digest_value == EXPECTED_SHARD_FILE_DIGEST,
        "shard file digest changed",
    )
    require(
        ledger_digest_value == EXPECTED_SHARD_LEDGER_DIGEST,
        "shard semantic digest changed",
    )
    require(
        len(closure_tuple) == 132
        and body_digest(closure_tuple) == EXPECTED_CLOSURE_DIGEST,
        "closure root set changed",
    )
    require(
        len(PROVED_17) == 17
        and len(new_tuple) == 121
        and body_digest(new_tuple) == EXPECTED_NEW_DIGEST,
        "new-root set changed",
    )
    require(
        len(union_tuple) == 138
        and body_digest(union_tuple) == EXPECTED_UNION_DIGEST,
        "proved-union set changed",
    )

    print("LRC14 j6 ranked-H1 32-shard locked merge")
    print(f"counts={tuple(counts)}")
    print(f"statuses={tuple(sorted(statuses.items()))}")
    print(f"H1_quantiles={EXPECTED_H_QUANTILES}")
    print(f"star_exceptions={len(exceptions)};anatomy=one_index0_pivot")
    print(f"minimum_positive_depth1_margin={min(minimum_margins)}")
    print(f"closure_roots={len(closure_tuple)}")
    for body in closure_tuple:
        print("CLOSED_ROOT=" + ",".join(map(str, body)))
    print(f"previously_proved_roots={len(PROVED_17)}")
    print(f"new_roots={len(new_tuple)}")
    for body in new_tuple:
        print("NEW_ROOT=" + ",".join(map(str, body)))
    print(f"proved_union_after_merge={len(union_tuple)}")
    print(f"closure_root_digest={body_digest(closure_tuple)}")
    print(f"new_root_digest={body_digest(new_tuple)}")
    print(f"proved_union_digest={body_digest(union_tuple)}")
    print(f"shard_file_digest={file_digest_value}")
    print(f"shard_ledger_digest={ledger_digest_value}")
    print("scope=cutoff<=15000 finite r4/H1 layer;181 tails deferred")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
