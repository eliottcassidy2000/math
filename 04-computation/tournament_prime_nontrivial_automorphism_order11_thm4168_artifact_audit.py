#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import os
import re
import subprocess
import sys
from pathlib import Path


N = 11
ROOT = Path(__file__).resolve().parents[1]
LABELS = ROOT / "05-knowledge/results/tournament_prime_nontrivial_automorphism_order11_thm4168.labels"
DIGRAPH6 = ROOT / "05-knowledge/results/tournament_prime_nontrivial_automorphism_order11_thm4168.d6"
LABELG = os.environ.get("LABELG_BIN", "/opt/homebrew/bin/labelg")
COUNTG = os.environ.get("COUNTG_BIN", "/opt/homebrew/bin/countg")
HOSTILE = "0000001100100001000000000000111111011101111111111000101"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def adjacency_from_label(label: str, order: int) -> list[int]:
    if len(label) != order * (order - 1) // 2 or set(label) > {"0", "1"}:
        raise ValueError("invalid tournament label")
    out = [0] * order
    offset = 0
    for i in range(order):
        for j in range(i + 1, order):
            if label[offset] == "1":
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            offset += 1
    return out


def label_from_adjacency(out: list[int]) -> str:
    return "".join(
        "1" if out[i] & (1 << j) else "0"
        for i in range(len(out))
        for j in range(i + 1, len(out))
    )


def decode_digraph6(code: str) -> list[int]:
    if not code.startswith("&"):
        raise ValueError("not digraph6")
    order = ord(code[1]) - 63
    bits = "".join(f"{ord(char) - 63:06b}" for char in code[2:])
    if len(bits) < order * order or set(bits[order * order :]) > {"0"}:
        raise ValueError("bad digraph6 padding")
    out = [0] * order
    for i in range(order):
        if bits[i * order + i] != "0":
            raise ValueError("digraph6 loop")
        for j in range(order):
            if bits[i * order + j] == "1":
                out[i] |= 1 << j
    for i in range(order):
        for j in range(i + 1, order):
            if ((out[i] >> j) & 1) + ((out[j] >> i) & 1) != 1:
                raise ValueError("digraph6 is not a tournament")
    return out


def encode_digraph6(out: list[int]) -> str:
    order = len(out)
    bits = "".join(
        "1" if out[i] & (1 << j) else "0"
        for i in range(order)
        for j in range(order)
    )
    bits += "0" * (-len(bits) % 6)
    payload = "".join(
        chr(63 + int(bits[start : start + 6], 2))
        for start in range(0, len(bits), 6)
    )
    return "&" + chr(63 + order) + payload


def module_prime(out: list[int]) -> tuple[bool, int | None]:
    order = len(out)
    full = (1 << order) - 1
    for subset in range(1 << order):
        size = subset.bit_count()
        if size < 2 or size == order:
            continue
        outside = full ^ subset
        while outside:
            bit = outside & -outside
            vertex = bit.bit_length() - 1
            outside ^= bit
            toward = (out[vertex] & subset).bit_count()
            if toward not in (0, size):
                break
        else:
            return False, subset
    return True, None


def converse(out: list[int]) -> list[int]:
    result = [0] * len(out)
    for i in range(len(out)):
        for j in range(len(out)):
            if out[j] & (1 << i):
                result[i] |= 1 << j
    return result


def delete_vertex(out: list[int], vertex: int) -> list[int]:
    keep = [old for old in range(len(out)) if old != vertex]
    position = {old: new for new, old in enumerate(keep)}
    card = [0] * len(keep)
    for old_i in keep:
        for old_j in keep:
            if out[old_i] & (1 << old_j):
                card[position[old_i]] |= 1 << position[old_j]
    return card


def canonicalize(codes: list[str]) -> list[str]:
    process = subprocess.run(
        [LABELG, "-q"],
        input="\n".join(codes) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )
    return process.stdout.splitlines()


def countg(codes: list[str], *arguments: str) -> str:
    process = subprocess.run(
        [COUNTG, *arguments, "-q"],
        input="\n".join(codes) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )
    return re.sub(r"; cpu=.*", "", process.stdout.strip())


def main() -> None:
    codes = DIGRAPH6.read_text("ascii").splitlines()
    derived_labels = [label_from_adjacency(decode_digraph6(code)) for code in codes]
    if "--emit-labels" in sys.argv[1:]:
        if len(codes) != 12_155 or codes != sorted(set(codes)):
            raise RuntimeError("unexpected canonical digraph6 stream")
        print("\n".join(derived_labels))
        return

    labels = LABELS.read_text("ascii").splitlines()
    if len(labels) != len(codes) or len(labels) != 12_155:
        raise RuntimeError("unexpected census size")
    if labels != sorted(set(labels)) or codes != sorted(set(codes)):
        raise RuntimeError("census is not sorted and unique")

    nonprime: list[tuple[int, int | None]] = []
    for line_number, code in enumerate(codes, 1):
        out = decode_digraph6(code)
        prime, witness = module_prime(out)
        if not prime:
            nonprime.append((line_number, witness))
    if derived_labels != labels or nonprime:
        raise RuntimeError(f"label/prime failure: {nonprime[:3]}")

    hostile_out = adjacency_from_label(HOSTILE, N)
    hostile_prime, hostile_module = module_prime(hostile_out)
    hostile_pair = canonicalize(
        [encode_digraph6(hostile_out), encode_digraph6(converse(hostile_out))]
    )
    cards = [delete_vertex(hostile_out, vertex) for vertex in range(N)]
    card_codes = [encode_digraph6(card) for card in cards]
    card_canons = canonicalize(card_codes)
    prime_cards = [
        vertex for vertex, card in enumerate(cards) if module_prime(card)[0]
    ]

    print(f"labels_sha256={sha256(LABELS)}")
    print(f"digraph6_sha256={sha256(DIGRAPH6)}")
    print(f"rows={len(labels)} rowwise_label_match=1 independently_prime={len(labels)}")
    print("automorphism_distribution=" + countg(codes, "--a").replace("\n", " | "))
    print(f"hostile_label={HOSTILE}")
    print(f"hostile_prime={int(hostile_prime)} module_witness={hostile_module}")
    print(f"hostile_self_converse={int(hostile_pair[0] == hostile_pair[1])}")
    print("hostile_automorphism=" + countg([encode_digraph6(hostile_out)], "--a").replace("\n", " | "))
    print(f"hostile_deletion_unique={len(set(card_canons))}")
    print("hostile_prime_card_indices=" + ",".join(map(str, prime_cards)))


if __name__ == "__main__":
    main()
