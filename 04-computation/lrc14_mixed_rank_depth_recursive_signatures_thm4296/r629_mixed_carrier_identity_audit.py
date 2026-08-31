#!/usr/bin/env python3
"""Maintained independent reconstruction of the post-r629 mixed carrier."""

from __future__ import annotations

import pathlib
import sys

MASK64 = (1 << 64) - 1
OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3

DELETED = (
    0x00003E1A, 0x000132A3, 0x00017464, 0x00033388,
    0x000A16C2, 0x000F8118, 0x00142A1A, 0x00154348,
    0x00184BA0, 0x001AA260, 0x00202C2B, 0x002066A4,
    0x002B018A, 0x0030C2A2,
)
EXCHANGE = (
    0x18468880, 0x080E8281, 0x22081017, 0x08422A82,
    0x004CAC40, 0x19C04044, 0x00C08EC0, 0x10443016,
    0x01609124, 0x10413209, 0x01611640, 0x00606449,
    0x0128D084, 0x08806449,
)
REPAIR_A = (0x2040C641, 0x00508325, 0x002A8641)
REPAIR_B = (0x00619324, 0x201813A4, 0x21888126)
POST_R630 = (0x0010E125,)
POST_R629 = (0x002AC4C0, 0x3882A082, 0x0041C325, 0x08C28E40)


def fnv_words(words: list[int] | tuple[int, ...]) -> int:
    state = OFFSET
    for word in words:
        for byte in word.to_bytes(8, "little"):
            state ^= byte
            state = (state * PRIME) & MASK64
    return state


def read_masks(path: pathlib.Path, count: int, digest: int, rank: int) -> list[int]:
    masks = [int(token, 16) for token in path.read_text(encoding="utf-8").split()]
    if (len(masks) != count or len(set(masks)) != count or
            any(mask.bit_count() != rank or mask >= 1 << 30 for mask in masks) or
            fnv_words(masks) != digest):
        raise RuntimeError(f"input identity changed: {path}")
    return masks


def append_distinct(carrier: list[int], additions: tuple[int, ...] | list[int]) -> None:
    present = set(carrier)
    for mask in additions:
        if mask in present:
            raise RuntimeError(f"duplicate addition {mask:08x}")
        carrier.append(mask)
        present.add(mask)


def main() -> None:
    if len(sys.argv) != 5:
        raise RuntimeError("usage: carrier BASE8951 ADD45 SUFFIX9 JOINT421")
    base = read_masks(pathlib.Path(sys.argv[1]), 8951, 0x188F82AB9DD1695A, 8)
    add45 = read_masks(pathlib.Path(sys.argv[2]), 45, 0xEC083B65CC8C34E3, 8)
    suffix = read_masks(pathlib.Path(sys.argv[3]), 9, 0x02B936529030E4BC, 8)
    joint = read_masks(pathlib.Path(sys.argv[4]), 421, 0x20D63DD42FE8150E, 8)

    carrier = list(base)
    append_distinct(carrier, add45)
    append_distinct(carrier, [0x014C9084])
    append_distinct(carrier, suffix)
    if len(carrier) != 9006 or not set(DELETED) <= set(carrier):
        raise RuntimeError("pre-exchange carrier changed")
    deleted = set(DELETED)
    carrier = [mask for mask in carrier if mask not in deleted]
    append_distinct(carrier, EXCHANGE)
    if len(carrier) != 9006 or fnv_words(carrier) != 0x8062CE6D5728DA1F:
        raise RuntimeError("endpoint-636 exchange identity changed")
    append_distinct(carrier, REPAIR_A)
    append_distinct(carrier, REPAIR_B)
    append_distinct(carrier, POST_R630)
    append_distinct(carrier, POST_R629)

    rank8 = sum(mask.bit_count() == 8 for mask in carrier)
    rank9 = sum(mask.bit_count() == 9 for mask in carrier)
    digest = fnv_words(carrier)
    if (len(carrier), rank8, rank9, digest) != (
            9017, 9011, 6, 0x07689A1534CE7327):
        raise RuntimeError("post-r629 mixed carrier identity changed")
    if not set(joint) <= set(carrier):
        raise RuntimeError("joint common deck absent from mixed carrier")

    print("R629_MIXED_CARRIER_IDENTITY_AUDIT_V1")
    print("BASE 8951 FNV 188f82ab9dd1695a ADD45 45 FNV ec083b65cc8c34e3")
    print("EXCHANGED 9006 FNV 8062ce6d5728da1f")
    print("R632_REPAIRS 6 R630_REPAIRS 1 R629_REPAIRS 4")
    print(f"FINAL {len(carrier)} FNV {digest:016x} RANK8 {rank8} RANK9 {rank9}")
    print("JOINT421_SUBSET 1 DISTINCT 1")
    print("VERDICT PASS EXACT_MIXED_CARRIER_IDENTITY")


if __name__ == "__main__":
    main()
