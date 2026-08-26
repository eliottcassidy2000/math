#!/usr/bin/env python3
"""Aggregate exact rooted-pair census shard summaries.

The tensor engine emits one summary per gentourng residue shard.  This merger
checks completeness and exact invariants, composes the two multiset semantic
accumulators, selects extrema by integer cross-products, and freezes a SHA-256
of the numerically ordered raw shard-summary files.
"""

from __future__ import annotations

import argparse
import hashlib
import re
from dataclasses import dataclass
from pathlib import Path


MASK64 = (1 << 64) - 1
EXPECTED_QUOTIENTS = 9_355_949
EXPECTED_ROOTED = 93_559_490
EXPECTED_EMPTY_SHARDS = 568
EXPECTED_INNER_LATTICE_NON2 = 1_454
EXPECTED_SCHEDULE_SHA256 = "79c212eda5863fffb201a064398e9d425f6c21f0f14579f2d4483e160bfb2013"
EXPECTED_RAW_SUMMARY_SHA256 = "ec09112ccf6b7c527dcf194490628c184402f5674909a7b7f03e9530c1d68c39"
EXPECTED_BOUNDARY_MANIFEST_SHA256 = "dd1380ab4f6204f27c54a94a6402b8dd3c44a7f9f604da9556544e62e598fc3b"
EXPECTED_SEMANTIC_SUM64 = 0xB93EC47D6452758C
EXPECTED_SEMANTIC_XOR64 = 0x7EF3CE694E8D49FE
HEADER_KEYS = frozenset({
    "quotients",
    "rooted_presentations",
    "rational_failures",
    "coset_failures",
    "inner_lattice_non2",
})
RECORD_KEYS = frozenset({
    "label", "H", "W2", "D4x4", "Chdx4", "exact_rho", "coset_margin"
})
QUOTIENT_KEYS = frozenset({"q", "root"})
BAD_QUOTIENT_KEYS = frozenset({"q", "root", "bad_mask", "max_gcd"})
SEMANTIC_KEYS = frozenset({"semantic_sum64", "semantic_xor64"})


def strict_fields(
    line: str,
    expected_keys: frozenset[str],
    tag: str | None = None,
) -> dict[str, str]:
    """Parse a whole summary line, rejecting duplicate or unlabelled tokens."""
    tokens = line.split()
    if tag is not None:
        if not tokens or tokens[0] != tag:
            raise RuntimeError(f"expected {tag!r}, got {line!r}")
        tokens = tokens[1:]
    parsed: dict[str, str] = {}
    for token in tokens:
        if token.count("=") != 1:
            raise RuntimeError(f"malformed summary token {token!r} in {line!r}")
        key, value = token.split("=", 1)
        if not re.fullmatch(r"[A-Za-z0-9_]+", key) or not value:
            raise RuntimeError(f"malformed summary field {token!r} in {line!r}")
        if key in parsed:
            raise RuntimeError(f"duplicate summary field {key!r} in {line!r}")
        parsed[key] = value
    if parsed.keys() != expected_keys:
        raise RuntimeError(
            f"summary keys {sorted(parsed)} != expected {sorted(expected_keys)} "
            f"in {line!r}"
        )
    return parsed


@dataclass(frozen=True)
class Record:
    tag: str
    label: str
    h: int
    w2: int
    d4x4: int
    chdx4: int
    margin: int
    q: str
    root: int
    bad_mask: int = 0
    max_gcd: int = 0

    @property
    def tie_key(self) -> tuple[str, str, int]:
        return self.label, self.q, self.root


def expand_label(q: str, root: int) -> str:
    """Expand quotient vertex ``root`` into the ordered arc root->10."""
    adjacency = [[False] * 10 for _ in range(10)]
    bit = 0
    for left in range(10):
        for right in range(left + 1, 10):
            adjacency[left][right] = q[bit] == "1"
            adjacency[right][left] = not adjacency[left][right]
            bit += 1
    child = [[False] * 11 for _ in range(11)]
    for left in range(10):
        for right in range(10):
            child[left][right] = adjacency[left][right]
    for vertex in range(10):
        if vertex == root:
            continue
        child[10][vertex] = adjacency[root][vertex]
        child[vertex][10] = adjacency[vertex][root]
    child[root][10] = True
    return "".join(
        "1" if child[left][right] else "0"
        for left in range(11)
        for right in range(left + 1, 11)
    )


def record(lines: list[str], tag: str) -> Record | None:
    for index, line in enumerate(lines):
        if not line.startswith(tag + " "):
            continue
        if index + 1 >= len(lines) or not lines[index + 1].startswith("q="):
            raise RuntimeError(f"{tag}: missing quotient row")
        a = strict_fields(line, RECORD_KEYS, tag)
        q_keys = BAD_QUOTIENT_KEYS if tag == "FIRST_INNER_LATTICE_NON2" else QUOTIENT_KEYS
        b = strict_fields(lines[index + 1], q_keys)
        item = Record(
            tag=tag,
            label=a["label"],
            h=int(a["H"]),
            w2=int(a["W2"]),
            d4x4=int(a["D4x4"]),
            chdx4=int(a["Chdx4"]),
            margin=int(a["coset_margin"]),
            q=b["q"],
            root=int(b["root"]),
            bad_mask=int(b.get("bad_mask", "0")),
            max_gcd=int(b.get("max_gcd", "0")),
        )
        expected_rho = f"{2 * abs(item.chdx4)}/{item.d4x4}"
        if a["exact_rho"] != expected_rho:
            raise RuntimeError(
                f"{tag}: exact_rho={a['exact_rho']} != {expected_rho}"
            )
        if not re.fullmatch(r"[01]{55}", item.label):
            raise RuntimeError(f"{tag}: malformed rooted label {item.label!r}")
        if not re.fullmatch(r"[01]{45}", item.q):
            raise RuntimeError(f"{tag}: malformed quotient label {item.q!r}")
        if not (-1 <= item.root < 10):
            raise RuntimeError(f"{tag}: root {item.root} outside sentinel/quotient range")
        if item.root >= 0 and expand_label(item.q, item.root) != item.label:
            raise RuntimeError(f"{tag}: rooted label is not the stated pair expansion")
        return item
    return None


def max_rho(left: Record | None, right: Record) -> Record:
    if left is None:
        return right
    lhs = abs(left.chdx4) * right.d4x4
    rhs = abs(right.chdx4) * left.d4x4
    if rhs > lhs or (rhs == lhs and right.tie_key < left.tie_key):
        return right
    return left


def min_margin(left: Record | None, right: Record) -> Record:
    if left is None:
        return right
    if right.margin < left.margin or (
        right.margin == left.margin and right.tie_key < left.tie_key
    ):
        return right
    return left


def show(item: Record) -> str:
    return (
        f"{item.tag} label={item.label} H={item.h} W2={item.w2} "
        f"D4x4={item.d4x4} Chdx4={item.chdx4} "
        f"exact_rho={2 * abs(item.chdx4)}/{item.d4x4} "
        f"coset_margin={item.margin}\nq={item.q} root={item.root}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("--shards", type=int, default=1024)
    args = parser.parse_args()
    if args.shards != 1024:
        raise RuntimeError("the frozen census contract requires exactly 1024 shards")

    expected_names = {f"shard_{shard}.out" for shard in range(args.shards)}
    required_names = expected_names | {"completed.log"}
    actual_names = {path.name for path in args.directory.iterdir()}
    missing_names = sorted(required_names - actual_names)
    extra_names = sorted(actual_names - required_names)
    if missing_names:
        raise RuntimeError(f"missing shard entries: {missing_names[:8]}")
    if extra_names:
        raise RuntimeError(f"unexpected census entries: {extra_names[:8]}")
    completion_path = args.directory / "completed.log"
    if not completion_path.is_file():
        raise RuntimeError("completed.log is not a regular file")
    completion_lines = completion_path.read_text(encoding="ascii").splitlines()
    if any(re.fullmatch(r"[0-9]+", line) is None for line in completion_lines):
        raise RuntimeError("completed.log contains a malformed shard identifier")
    completion_ids = [int(line) for line in completion_lines]
    if len(completion_ids) != args.shards or sorted(completion_ids) != list(range(args.shards)):
        raise RuntimeError("completed.log is not a one-to-one certificate for shards 0..1023")

    totals = {
        "quotients": 0,
        "rooted_presentations": 0,
        "rational_failures": 0,
        "coset_failures": 0,
        "inner_lattice_non2": 0,
    }
    semantic_sum = 0
    semantic_xor = 0
    raw_summary_digest = hashlib.sha256()
    boundary_manifest_digest = hashlib.sha256()
    schedule_digest = hashlib.sha256()
    global_max = None
    global_min = None
    first_bad = None
    empty_shards = 0

    for shard in range(args.shards):
        path = args.directory / f"shard_{shard}.out"
        if not path.is_file():
            raise RuntimeError(f"missing shard {shard}: {path}")
        raw = path.read_bytes()
        raw_summary_digest.update(raw)
        boundary_manifest_digest.update(
            f"{shard} {len(raw)} {hashlib.sha256(raw).hexdigest()}\n".encode("ascii")
        )
        lines = raw.decode("ascii").splitlines()
        if not lines:
            raise RuntimeError(f"empty shard {shard}")
        head = strict_fields(lines[0], HEADER_KEYS)
        if any(int(head[key]) < 0 for key in HEADER_KEYS):
            raise RuntimeError(f"shard {shard}: negative summary counter")
        shard_quotients = int(head["quotients"])
        shard_rooted = int(head["rooted_presentations"])
        shard_bad_count = int(head["inner_lattice_non2"])
        expected_line_count = 8 if shard_bad_count > 0 else 6
        if len(lines) != expected_line_count:
            raise RuntimeError(
                f"shard {shard}: {len(lines)} summary lines != {expected_line_count}"
            )
        strict_fields(lines[1], RECORD_KEYS, "ROOTED_MAX_RHO")
        strict_fields(lines[2], QUOTIENT_KEYS)
        strict_fields(lines[3], RECORD_KEYS, "ROOTED_MIN_COSET")
        strict_fields(lines[4], QUOTIENT_KEYS)
        if shard_bad_count > 0:
            strict_fields(lines[5], RECORD_KEYS, "FIRST_INNER_LATTICE_NON2")
            strict_fields(lines[6], BAD_QUOTIENT_KEYS)
        semantic_index = expected_line_count - 1
        sem = strict_fields(lines[semantic_index], SEMANTIC_KEYS)
        schedule_digest.update(f"{shard} {shard_quotients}\n".encode("ascii"))
        if shard_rooted != 10 * shard_quotients:
            raise RuntimeError(
                f"shard {shard}: rooted={shard_rooted} != 10*{shard_quotients}"
            )
        for key in totals:
            totals[key] += int(head[key])
        shard_semantic_sum = int(sem["semantic_sum64"], 16)
        shard_semantic_xor = int(sem["semantic_xor64"], 16)
        semantic_sum = (semantic_sum + shard_semantic_sum) & MASK64
        semantic_xor ^= shard_semantic_xor

        shard_max = record(lines, "ROOTED_MAX_RHO")
        shard_min = record(lines, "ROOTED_MIN_COSET")
        if shard_max is None or shard_min is None:
            raise RuntimeError(f"shard {shard}: missing extrema")
        shard_bad = record(lines, "FIRST_INNER_LATTICE_NON2")
        if (shard_bad_count > 0) != (shard_bad is not None):
            raise RuntimeError(f"shard {shard}: inconsistent inner-lattice witness")

        # gentourng's modular split legitimately has many empty residues.
        # Their engine extrema are sentinel 0/0 rows and must not enter the
        # global comparison; all actual counters and semantics must be zero.
        if shard_quotients == 0:
            empty_shards += 1
            if any(int(head[key]) for key in (
                "rooted_presentations", "rational_failures",
                "coset_failures", "inner_lattice_non2"
            )):
                raise RuntimeError(f"shard {shard}: nonzero empty-shard counter")
            if shard_semantic_sum or shard_semantic_xor:
                raise RuntimeError(f"shard {shard}: nonzero empty-shard semantics")
            for item in (shard_max, shard_min):
                if (
                    item.label != "0" * 55
                    or item.q != "0" * 45
                    or item.root != -1
                    or any((item.h, item.w2, item.d4x4, item.chdx4))
                    or item.margin != (1 << 63) - 1
                ):
                    raise RuntimeError(f"shard {shard}: malformed empty-shard sentinel")
            continue

        if shard_max.root < 0 or shard_min.root < 0:
            raise RuntimeError(f"shard {shard}: sentinel extremum on a nonempty shard")
        if shard_max.d4x4 <= 0 or 2 * abs(shard_max.chdx4) >= shard_max.d4x4:
            raise RuntimeError(f"shard {shard}: noncentral rational maximum")
        if shard_min.margin <= 0:
            raise RuntimeError(f"shard {shard}: nonpositive exact-coset margin")
        global_max = max_rho(global_max, shard_max)
        global_min = min_margin(global_min, shard_min)
        if first_bad is None and shard_bad is not None:
            first_bad = shard_bad

    if totals["quotients"] != EXPECTED_QUOTIENTS:
        raise RuntimeError(f"quotient total {totals['quotients']} != {EXPECTED_QUOTIENTS}")
    if totals["rooted_presentations"] != EXPECTED_ROOTED:
        raise RuntimeError(f"rooted total {totals['rooted_presentations']} != {EXPECTED_ROOTED}")
    if empty_shards != EXPECTED_EMPTY_SHARDS:
        raise RuntimeError(f"empty shards {empty_shards} != {EXPECTED_EMPTY_SHARDS}")
    if totals["inner_lattice_non2"] != EXPECTED_INNER_LATTICE_NON2:
        raise RuntimeError(
            f"inner-lattice total {totals['inner_lattice_non2']} "
            f"!= {EXPECTED_INNER_LATTICE_NON2}"
        )
    schedule_sha256 = schedule_digest.hexdigest()
    if schedule_sha256 != EXPECTED_SCHEDULE_SHA256:
        raise RuntimeError(
            f"quotient schedule {schedule_sha256} != {EXPECTED_SCHEDULE_SHA256}"
        )
    if totals["rational_failures"] or totals["coset_failures"]:
        raise RuntimeError(f"failure totals are nonzero: {totals}")
    assert global_max is not None and global_min is not None
    raw_summary_sha256 = raw_summary_digest.hexdigest()
    boundary_manifest_sha256 = boundary_manifest_digest.hexdigest()
    if raw_summary_sha256 != EXPECTED_RAW_SUMMARY_SHA256:
        raise RuntimeError(
            f"ordered raw-summary digest {raw_summary_sha256} "
            f"!= {EXPECTED_RAW_SUMMARY_SHA256}"
        )
    if boundary_manifest_sha256 != EXPECTED_BOUNDARY_MANIFEST_SHA256:
        raise RuntimeError(
            f"boundary manifest {boundary_manifest_sha256} "
            f"!= {EXPECTED_BOUNDARY_MANIFEST_SHA256}"
        )
    if semantic_sum != EXPECTED_SEMANTIC_SUM64 or semantic_xor != EXPECTED_SEMANTIC_XOR64:
        raise RuntimeError("semantic accumulators do not match the frozen census")

    expected_max = Record(
        "ROOTED_MAX_RHO",
        "1011111101111111111111111111111111110011111101110111110",
        2727, 110646, 3736874076, 1052690016, 818,
        "101111110111111111111111111111110011111111111", 4,
    )
    expected_min = Record(
        "ROOTED_MIN_COSET",
        "1011011111111111111111111111111111111110110101110110000",
        243, 23742, 169212908, 3515272, 380,
        "101101111111111111111111111111111111101111110", 3,
    )
    expected_first_bad = Record(
        "FIRST_INNER_LATTICE_NON2",
        "1111111100111111100111111111111111110011110111111111111",
        21135, 463230, 68646931020, 4844970720, 7074,
        "111111110111111101111111111111110011101111111", 9, 1020, 6,
    )
    if global_max != expected_max or global_min != expected_min or first_bad != expected_first_bad:
        raise RuntimeError("extrema or first non-parity witness differ from the frozen census")

    print("ORDER11 ROOTED PAIR TENSOR CENSUS")
    print(f"shards={args.shards}")
    print(f"nonempty_shards={args.shards-empty_shards} empty_shards={empty_shards}")
    print(
        " ".join(f"{key}={value}" for key, value in totals.items())
    )
    print("one_exact_maximum_witness=below")
    print(show(global_max))
    print("one_exact_minimum_margin_witness=below")
    print(show(global_min))
    if first_bad is not None:
        print(show(first_bad))
        print(f"bad_mask={first_bad.bad_mask} max_gcd={first_bad.max_gcd}")
    print(f"semantic_sum64={semantic_sum:016x} semantic_xor64={semantic_xor:016x}")
    print(f"quotient_schedule_sha256={schedule_sha256}")
    print(f"ordered_raw_summary_files_sha256={raw_summary_sha256}")
    print(f"boundary_manifest_sha256={boundary_manifest_sha256}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
