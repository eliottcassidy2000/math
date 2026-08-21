#!/usr/bin/env python3
"""Exact minimal carry-closed observer bank for the LRC exceptional detector."""

from __future__ import annotations

import ast
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT_SCRIPT = ROOT / "04-computation/lrc_exceptional_detector_simple_spectrum_thm3661.py"
PARENT_OUTPUT = ROOT / "05-knowledge/results/lrc_exceptional_detector_simple_spectrum_thm3661.out"
EXPECTED_PARENT_HASHES = (
    "7c8ebacaff5318b58ed1588d2e38edf1bc3c14ea81beec069d781450db986251",
    "61893d6e628232e877c73d6102f042654327fddfe7d51d1cf79c1aa6bc65419f",
)
EXPECTED_SEMANTIC_SHA256 = "953994bb4440423e1569e8c35165a2a60c8d7c3fa438b92508203ab64f3085a9"

P = 13
N = 169
X_PLUS = frozenset(((12, 0), (0, 11), (6, 5), (9, 3)))
X_MINUS = frozenset(((0, 12), (12, 1), (6, 7), (3, 9)))
THREE_COLUMNS = frozenset((12, 0, 1))
FULL_BOUNDARY = (
    (0, 10), (0, 11), (0, 12),
    (3, 8), (3, 9),
    (6, 4), (6, 5), (6, 6), (6, 7),
    (9, 2), (9, 3),
    (12, 0), (12, 1), (12, 12),
)
CARRY_BOUNDARY = FULL_BOUNDARY[:-3]


def require(condition: bool, payload: object) -> None:
    if condition is not True:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def split_add(left, right):
    return ((left[0] + right[0]) % P, (left[1] + right[1]) % P)


def split_subtract(left, right):
    return ((left[0] - right[0]) % P, (left[1] - right[1]) % P)


def split_negate(label):
    return ((-label[0]) % P, (-label[1]) % P)


def assemble(label):
    return label[0] + P * label[1]


def disassemble(value):
    value %= N
    return value % P, value // P


def cyclic_add(left, right):
    return disassemble(assemble(left) + assemble(right))


def cyclic_subtract(left, right):
    return disassemble(assemble(left) - assemble(right))


def cyclic_negate(label):
    return disassemble(-assemble(label))


def reversal(label):
    return (12 - label[0]) % P, (12 - label[1]) % P


def main() -> None:
    parent_hashes = (lf_sha256(PARENT_SCRIPT), lf_sha256(PARENT_OUTPUT))
    require(parent_hashes == EXPECTED_PARENT_HASHES,
            ("THM-3661 parent hashes", parent_hashes))
    labels = tuple((r0, r1) for r0 in range(P) for r1 in range(P))
    detector = {
        label: int(label in X_PLUS) - int(label in X_MINUS)
        for label in labels
    }
    require(all(detector[reversal(label)] == -detector[label]
                for label in labels), "detector reversal")

    laws = {
        "split": (split_add, split_subtract, split_negate),
        "cyclic": (cyclic_add, cyclic_subtract, cyclic_negate),
    }

    def bank(columns):
        return tuple((column, high) for column in sorted(columns)
                     for high in range(P))

    def codes(offsets, subtract):
        return tuple(
            tuple(detector[subtract(label, offset)] for offset in offsets)
            for label in labels
        )

    observer_bank = bank(THREE_COLUMNS)
    require(len(observer_bank) == 39, "observer bank size")
    vertical = (0, 1)
    law_records = {}
    for name, (add, subtract, negate) in laws.items():
        observer_codes = codes(observer_bank, subtract)
        require(len(set(observer_codes)) == N,
                (name, "observer not injective"))
        index = {offset: position for position, offset in enumerate(observer_bank)}
        carry_permutation = tuple(
            index[subtract(offset, vertical)] for offset in observer_bank
        )
        reversal_permutation = tuple(
            index[negate(offset)] for offset in observer_bank
        )
        require(all(
            observer_codes[labels.index(add(label, vertical))][position]
            == observer_codes[label_index][carry_permutation[position]]
            for label_index, label in enumerate(labels)
            for position in range(len(observer_bank))
        ), (name, "carry coordinate action"))
        require(all(
            observer_codes[labels.index(reversal(label))][position]
            == -observer_codes[label_index][reversal_permutation[position]]
            for label_index, label in enumerate(labels)
            for position in range(len(observer_bank))
        ), (name, "reversal coordinate action"))

        successes_by_size = {}
        for size in (1, 2, 3):
            successes = []
            for columns in combinations(range(P), size):
                offsets = bank(columns)
                if len(set(codes(offsets, subtract))) == N:
                    successes.append(columns)
            successes_by_size[size] = tuple(successes)
        expected_three = tuple(
            tuple(sorted(((center - 1) % P, center, (center + 1) % P)))
            for center in range(P)
        )
        expected_three = tuple(sorted(set(expected_three)))
        require(successes_by_size[1] == () and successes_by_size[2] == (),
                (name, "smaller carry-closed bank separates"))
        require(tuple(sorted(successes_by_size[3])) == expected_three
                and len(expected_three) == P,
                (name, "three-column success atlas", successes_by_size[3]))

        carry_hostile = codes(CARRY_BOUNDARY, subtract)
        full_hostile = codes(FULL_BOUNDARY, subtract)
        zero = (0,) * len(CARRY_BOUNDARY)
        full_zero = (0,) * len(FULL_BOUNDARY)
        index00 = labels.index((0, 0))
        index01 = labels.index((0, 1))
        require(carry_hostile[index00] == carry_hostile[index01] == zero,
                (name, "carry-boundary hostile"))
        require(full_hostile[index00] == full_hostile[index01] == full_zero,
                (name, "full-boundary hostile"))

        law_records[name] = (
            digest_json(tuple(zip(labels, observer_codes))),
            carry_permutation,
            reversal_permutation,
            successes_by_size,
            digest_json(tuple(zip(labels, carry_hostile))),
            digest_json(tuple(zip(labels, full_hostile))),
        )

    require(law_records["split"][0] != law_records["cyclic"][0],
            "group-law observers unexpectedly identical")
    support_columns = frozenset(label[0] for label in X_PLUS | X_MINUS)
    require(support_columns == frozenset((0, 3, 6, 9, 12)),
            "detector column support")
    unrestricted_lower = (N - 2 + len(X_PLUS | X_MINUS) - 1) // len(X_PLUS | X_MINUS)
    require(unrestricted_lower == 21, "unrestricted coverage lower bound")

    semantic = digest_json((
        P, N, parent_hashes,
        tuple(sorted(X_PLUS)), tuple(sorted(X_MINUS)),
        observer_bank, tuple(sorted(support_columns)),
        law_records, unrestricted_lower, CARRY_BOUNDARY, FULL_BOUNDARY,
    ))
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    source = Path(__file__).resolve()
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source.read_text(encoding="utf-8")))),
            "Python assert node present")
    print("== THM-3663 LRC minimal carry-closed exceptional observer ==")
    print(f"address_set=169;detector_support=8;columns={tuple(sorted(support_columns))}")
    print(f"parent_sha256_lf={parent_hashes}")
    print(f"observer_columns={tuple(sorted(THREE_COLUMNS))};coordinates=39")
    print("observer_injective=split:169/169,cyclic:169/169")
    print("operation_actions=vertical_carry:coordinate_permutation;reversal:signed_coordinate_permutation")
    print("carry_closed_minimum=3 columns=39 coordinates")
    print("success_counts_by_columns=(size1:0,size2:0,size3:13)")
    print("successful_size3_banks=exactly cyclically consecutive column triples")
    print(f"split_semantic_sha256={law_records['split'][0]}")
    print(f"cyclic_semantic_sha256={law_records['cyclic'][0]}")
    print("boundary_bank_hostile=O(0,0)=O(0,1)=zero for both 11 and 14 offsets, both laws")
    print(f"unrestricted_translate_minimum_interval=[{unrestricted_lower},39]")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("scope=address-level observer congruence;not physical current/first-return/entry/LRC14")
    print("PASS")


if __name__ == "__main__":
    main()
