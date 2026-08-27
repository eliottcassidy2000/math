#!/usr/bin/env python3
"""Portable independent certificate for the round-five minimum augmentation.

This checker consumes the frozen full-universe pattern atlas and the compact
carrier transcript.  It does not solve the cover problem by trusting the
chosen six masks: a six-obligation packing proves the lower bound, while the
six displayed atlas patterns prove the matching upper bound.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import re


FAIL_256 = (
    0x06166401, 0x07067001, 0x07106409, 0x07126088, 0x07126401,
    0x07162401, 0x07163400, 0x07166008, 0x0D107401, 0x0D10E401,
    0x0D146401, 0x0D186401, 0x0D246401, 0x0D506401, 0x0D906401,
    0x0F106401, 0x0F142401, 0x0F142408, 0x0F143400, 0x0F146008,
    0x15923400, 0x17162400, 0x17922008, 0x1D106401, 0x1D902401,
    0x1F142400,
)
FAIL_384 = (0x0D186401,)
ALL_OBLIGATIONS = (1 << 27) - 1
PACKING_BITS = (0, 9, 12, 16, 19, 20)
PACKING_MASK = sum(1 << bit for bit in PACKING_BITS)
SELECTED = (
    (0x00289285, 0x026A0080),
    (0x0260812C, 0x05904D00),
    (0x18689040, 0x000000BD),
    (0x20C0C124, 0x02270060),
    (0x302C1006, 0x0000E21C),
    (0x30580888, 0x00001002),
)
POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
    120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252, 264,
    286, 290,
)

OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def fnv_words(words: list[int]) -> str:
    state = OFFSET
    for word in words:
        for shift in range(0, 64, 8):
            state ^= (word >> shift) & 0xFF
            state = state * PRIME & MASK64
    return f"{state:016x}"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def labels(mask: int) -> str:
    return "{" + ",".join(
        str(label) for bit, label in enumerate(POOL) if mask & (1 << bit)
    ) + "}"


def parse_atlas(path: Path) -> list[tuple[int, int, int]]:
    text = path.read_text(encoding="ascii")
    assert text.startswith("ROUND5_FULL_UNIVERSE_AUGMENTATION_V1\n")
    census = re.search(
        r"^ACTIVE_256 (\d+) ACTIVE_384 (\d+) COVERING_CANDIDATES (\d+) "
        r"DISTINCT_PATTERNS (\d+)$", text, re.MULTILINE
    )
    assert census is not None
    assert tuple(map(int, census.groups())) == (1_721_339, 2_444_056,
                                                19_156, 330)
    rows = [
        (int(pattern, 16), int(least, 16), int(cover))
        for pattern, least, cover in re.findall(
            r"^PATTERN ([0-9a-f]+) LEAST ([0-9a-f]+) COVER (\d+)$",
            text,
            re.MULTILINE,
        )
    ]
    assert len(rows) == 330
    assert [pattern for pattern, _, _ in rows] == sorted(
        {pattern for pattern, _, _ in rows}
    )
    assert all(0 < pattern <= ALL_OBLIGATIONS for pattern, _, _ in rows)
    assert all(least.bit_count() == 8 and least < (1 << 30)
               for _, least, _ in rows)
    assert all(pattern.bit_count() == cover for pattern, _, cover in rows)
    assert text.rstrip().endswith("VERDICT PASS FULL_REPAIR_PATTERN_ATLAS")
    return rows


def parse_mask_list(text: str, q: int, r: int) -> tuple[int, ...]:
    found = re.search(
        rf"^ROUND4_FAILURE_LIST PAIR {q},{r} ACTIVE \d+ FAILURES \d+ "
        r"MASKS_HEX ([0-9a-f,]+)$",
        text,
        re.MULTILINE,
    )
    assert found is not None
    return tuple(int(word, 16) for word in found.group(1).split(","))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("atlas", type=Path)
    parser.add_argument("compact_transcript", type=Path)
    parser.add_argument("discovery_256", type=Path)
    parser.add_argument("discovery_384", type=Path)
    args = parser.parse_args()

    rows = parse_atlas(args.atlas)
    assert sha(args.atlas) == (
        "3717bc24ed8cb6167cd0194606a453b443ece91b7677ae9551573d9e864ab274"
    )
    compact = args.compact_transcript.read_text(encoding="ascii")
    assert parse_mask_list(compact, 256, 671) == FAIL_256
    assert parse_mask_list(compact, 384, 671) == FAIL_384

    discovery256 = args.discovery_256.read_text(encoding="ascii")
    discovery384 = args.discovery_384.read_text(encoding="ascii")
    assert "PAIR 256,671 " in discovery256
    assert "REPAIRS 5852925 ACTIVE 1721339 INACTIVE 4131586 EQUALITIES 0 " \
           in discovery256
    assert "PREFIX_CERT SIZE 13086 FNV 6c1f44ec0661d8c " in discovery256
    assert "PAIR 384,671 " in discovery384
    assert "REPAIRS 5852925 ACTIVE 2444056 INACTIVE 3408869 EQUALITIES 0 " \
           in discovery384
    assert "PREFIX_CERT SIZE 6986 FNV 3a02a5774d3641ab " in discovery384

    pattern_to_least = {pattern: least for pattern, least, _ in rows}
    least_to_pattern = {least: pattern for pattern, least, _ in rows}
    assert len(pattern_to_least) == 330
    for least, pattern in SELECTED:
        assert pattern_to_least[pattern] == least
        assert least_to_pattern[least] == pattern

    selected_patterns = [pattern for _, pattern in SELECTED]
    assert len(set(selected_patterns)) == 6
    assert set().union(*(
        {bit for bit in range(27) if pattern & (1 << bit)}
        for pattern in selected_patterns
    )) == set(range(27))
    assert sum(selected_patterns, 0) != ALL_OBLIGATIONS  # OR, not addition.
    union = 0
    for pattern in selected_patterns:
        union |= pattern
    assert union == ALL_OBLIGATIONS

    # Packing lower bound: no realizable repair pattern covers two of these
    # six obligations. Hence every augmentation uses at least six repairs.
    assert all((pattern & PACKING_MASK).bit_count() <= 1
               for pattern, _, _ in rows)
    assert all((pattern & PACKING_MASK).bit_count() == 1
               for pattern in selected_patterns)
    for left_index, left in enumerate(PACKING_BITS):
        for right in PACKING_BITS[left_index + 1:]:
            assert not any(
                pattern & (1 << left) and pattern & (1 << right)
                for pattern, _, _ in rows
            )

    atlas_words: list[int] = []
    for pattern, least, cover in rows:
        atlas_words.extend((pattern, least, cover))
    selected_words: list[int] = []
    for least, pattern in SELECTED:
        selected_words.extend((least, pattern))

    print("ROUND5_MINIMUM_SET_COVER_INDEPENDENT_AUDIT_V1")
    print("OBLIGATIONS 27 PAIR256 26 PAIR384 1")
    print("ATLAS PATTERNS", len(rows), "FNV", fnv_words(atlas_words),
          "SHA256", sha(args.atlas))
    print("PACKING_BITS " + ",".join(map(str, PACKING_BITS)))
    print("PACKING_BODIES " + ",".join(
        f"{FAIL_256[bit]:x}" for bit in PACKING_BITS
    ))
    print("PACKING_LABELS " + ";".join(
        labels(FAIL_256[bit]) for bit in PACKING_BITS
    ))
    print("PAIRWISE_INCOMPATIBLE_PAIRS 15 MAX_PACKING_HITS_PER_PATTERN 1")
    print("LOWER_BOUND_MECHANISM EACH_REPAIR_PAYS_AT_MOST_ONE_PACKING_OBLIGATION")
    print("LOWER_BOUND 6")
    print("SELECTED_MASKS " + ",".join(f"{least:x}" for least, _ in SELECTED))
    print("SELECTED_PATTERNS " + ",".join(
        f"{pattern:x}" for _, pattern in SELECTED
    ))
    print("SELECTED_FNV", fnv_words(selected_words),
          "UNION", f"{union:x}")
    print("UPPER_BOUND 6 MINIMUM 6")
    print("VERDICT PASS FULL_UNIVERSE_ATLAS_AND_MINIMUM_CERTIFICATE")


if __name__ == "__main__":
    main()
