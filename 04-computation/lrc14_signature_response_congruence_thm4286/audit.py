#!/usr/bin/env python3
"""Audit whether THM-4282 certificate responses factor through its signature.

This is a derived-ledger audit.  It does not recompute wall activity or body
coverage; those inputs are inherited from the independently audited THM-4282
packet.  The theorem-level consequence checked here is purely finite: a target
predicate factors through a quotient exactly when every quotient fibre is
pure for that predicate.
"""

from __future__ import annotations

import csv
import hashlib
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


REPO = Path(__file__).resolve().parents[2]
PACKET = (
    REPO
    / "05-knowledge/results/lrc14_inactive_signature_deck_surgery_thm4282"
)


@dataclass(frozen=True)
class InputSpec:
    name: str
    path: Path
    sha256: str
    count: int


SIGNATURE = InputSpec(
    "signature421",
    PACKET
    / "components/index367/results/post_thm4281_full421_inactive_signatures.csv",
    "4f0e8da3fdab8bd5a0e14f3b4fa30050602025f63486aa35e0cf03374e6e3832",
    24_223,
)

TARGETS = (
    InputSpec(
        "surgery520_367",
        PACKET / "components/surgery520/results/combined_gain_520_367.csv",
        "04aed4c107b3244a5e488266c46d1cae3bfffb20cddd831fc2d474c7b8c16a0e",
        586,
    ),
    InputSpec(
        "surgery256",
        PACKET / "components/surgery256/results/surgery_common.csv",
        "5946ff653c51a74eba09a14430c494074e53b5aba87c3159bd17bafbe9e605d5",
        188,
    ),
    InputSpec(
        "carrier90",
        PACKET / "results/carrier90.csv",
        "222afb7618d887f32847b4531ffedf5616f20c2196e92f52ca2c11b09e1eab1f",
        90,
    ),
    InputSpec(
        "union850",
        PACKET / "results/deletion_union850.csv",
        "7ad581bccd253e1778b972e8a303207da44534e6b995fa3ba15bd34b2801505b",
        850,
    ),
)

EXPECTED_TARGET_SUMMARIES = {
    # all-negative fibres, mixed fibres, all-positive fibres,
    # mixed rows, positive mixed rows, negative mixed rows
    "surgery520_367": (17_593, 6, 5, 770, 555, 215),
    "surgery256": (17_587, 15, 2, 367, 186, 181),
    "carrier90": (17_522, 34, 48, 1_460, 42, 1_418),
    "union850": (17_502, 48, 54, 2_351, 770, 1_581),
}

EXPECTED_JOINT_WORD_HISTOGRAM = {1: 17_556, 2: 42, 3: 5, 4: 1}
TOP_ROWS = (
    (220, 644),
    (256, 644),
    (258, 644),
    (294, 644),
    (366, 644),
    (416, 644),
    (512, 644),
)

Pair = tuple[int, int]
Signature = tuple[int, ...]
Word = tuple[int, ...]


def fail(message: str) -> None:
    raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    # The frozen packet hashes canonical LF text.  Git-for-Windows may expose
    # an older CSV with CRLF when that packet predates the repository-wide CSV
    # attribute, so normalize only the transport newline before hashing.
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def require_input(spec: InputSpec) -> None:
    if not spec.path.is_file():
        fail(f"missing input: {spec.path}")
    actual = file_sha256(spec.path)
    if actual != spec.sha256:
        fail(f"{spec.name} SHA-256 changed: {actual}")


def signature_indices(signature: Signature) -> tuple[int, ...]:
    return tuple(
        64 * word_index + bit
        for word_index, word in enumerate(signature)
        for bit in range(64)
        if (word >> bit) & 1
    )


def read_signatures(
    spec: InputSpec,
) -> tuple[dict[Signature, list[Pair]], dict[Pair, Signature]]:
    require_input(spec)
    expected_header = ["q", "r", "inactive_count"] + [
        f"w{i}" for i in range(7)
    ]
    fibres: dict[Signature, list[Pair]] = defaultdict(list)
    signature_of: dict[Pair, Signature] = {}
    with spec.path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != expected_header:
            fail(f"signature header changed: {reader.fieldnames}")
        for row_number, row in enumerate(reader, 2):
            pair = (int(row["q"]), int(row["r"]))
            if not (0 < pair[0] < pair[1]):
                fail(f"bad pair at signature row {row_number}: {pair}")
            if pair in signature_of:
                fail(f"duplicate signature pair: {pair}")
            signature = tuple(int(row[f"w{i}"], 16) for i in range(7))
            if signature[6] >> 37:
                fail(f"signature uses an index above 420 at pair {pair}")
            inactive_count = int(row["inactive_count"])
            if sum(word.bit_count() for word in signature) != inactive_count:
                fail(f"inactive-count mismatch at pair {pair}")
            signature_of[pair] = signature
            fibres[signature].append(pair)
    if len(signature_of) != spec.count:
        fail(f"signature row count changed: {len(signature_of)}")
    for pairs in fibres.values():
        pairs.sort()
    return dict(fibres), signature_of


def read_pair_set(spec: InputSpec, universe: set[Pair]) -> set[Pair]:
    require_input(spec)
    pairs: set[Pair] = set()
    with spec.path.open(newline="", encoding="utf-8") as handle:
        for row_number, row in enumerate(csv.reader(handle), 1):
            if len(row) != 2:
                fail(f"bad {spec.name} row {row_number}: {row}")
            try:
                pair = (int(row[0]), int(row[1]))
            except ValueError as error:
                fail(f"nonnumeric {spec.name} row {row_number}: {row}")
                raise AssertionError from error
            if pair in pairs:
                fail(f"duplicate {spec.name} pair: {pair}")
            if pair not in universe:
                fail(f"{spec.name} pair outside the signature universe: {pair}")
            pairs.add(pair)
    if len(pairs) != spec.count:
        fail(f"{spec.name} row count changed: {len(pairs)}")
    return pairs


def word_text(word: Iterable[int]) -> str:
    return "".join(str(bit) for bit in word)


def target_summary(
    fibres: dict[Signature, list[Pair]], target: set[Pair]
) -> tuple[tuple[int, int, int, int, int, int], tuple[object, ...]]:
    classes = Counter()
    mixed_rows = mixed_positive = mixed_negative = 0
    split_witnesses: list[tuple[object, ...]] = []
    for signature, pairs in fibres.items():
        positives = tuple(pair for pair in pairs if pair in target)
        negatives = tuple(pair for pair in pairs if pair not in target)
        if not positives:
            classes["negative"] += 1
        elif not negatives:
            classes["positive"] += 1
        else:
            classes["mixed"] += 1
            mixed_rows += len(pairs)
            mixed_positive += len(positives)
            mixed_negative += len(negatives)
            split_witnesses.append(
                (
                    len(pairs),
                    signature_indices(signature),
                    positives,
                    negatives,
                )
            )
    summary = (
        classes["negative"],
        classes["mixed"],
        classes["positive"],
        mixed_rows,
        mixed_positive,
        mixed_negative,
    )
    return summary, min(split_witnesses)


def format_pairs(pairs: Iterable[Pair]) -> str:
    return ";".join(f"{q},{r}" for q, r in pairs) or "-"


def main() -> None:
    fibres, signature_of = read_signatures(SIGNATURE)
    universe = set(signature_of)
    targets = {spec.name: read_pair_set(spec, universe) for spec in TARGETS}

    if targets["union850"] != (
        targets["surgery520_367"]
        | targets["surgery256"]
        | targets["carrier90"]
    ):
        fail("the 850-row union is not the union of its three typed targets")

    nonsingleton = [pairs for pairs in fibres.values() if len(pairs) > 1]
    if (len(fibres), len(nonsingleton), sum(map(len, nonsingleton))) != (
        17_604,
        803,
        7_422,
    ):
        fail("base signature-fibre census changed")

    print("THM-4286 SIGNATURE-RESPONSE FACTORIZATION AUDIT")
    print("STATUS FINITE-EXACT DERIVED FROM AUDITED THM-4282 LEDGERS")
    for spec in (SIGNATURE,) + TARGETS:
        print(
            f"INPUT {spec.name} ROWS {spec.count} SHA256 {spec.sha256} "
            f"PATH {spec.path.relative_to(REPO).as_posix()}"
        )
    print(
        "BASE ROWS 24223 FIBRES 17604 NONSINGLETON_FIBRES 803 "
        "NONSINGLETON_ROWS 7422"
    )

    for spec in TARGETS:
        summary, witness = target_summary(fibres, targets[spec.name])
        if summary != EXPECTED_TARGET_SUMMARIES[spec.name]:
            fail(f"{spec.name} summary changed: {summary}")
        negative, mixed, positive, rows, pos_rows, neg_rows = summary
        _, indices, positives, negatives = witness
        print(
            f"TARGET {spec.name} SIZE {len(targets[spec.name])} "
            f"ALL_NEGATIVE_FIBRES {negative} MIXED_FIBRES {mixed} "
            f"ALL_POSITIVE_FIBRES {positive} MIXED_ROWS {rows} "
            f"MIXED_POSITIVE_ROWS {pos_rows} MIXED_NEGATIVE_ROWS {neg_rows}"
        )
        print(
            f"LEAST_SPLIT {spec.name} SIGNATURE "
            f"{','.join(map(str, indices))} POSITIVE {format_pairs(positives)} "
            f"NEGATIVE {format_pairs(negatives)}"
        )

    word_histogram: Counter[int] = Counter()
    four_word_fibres: list[tuple[Signature, Counter[Word]]] = []
    target_sets = tuple(targets[spec.name] for spec in TARGETS)
    for signature, pairs in fibres.items():
        words = Counter(
            tuple(int(pair in target) for target in target_sets) for pair in pairs
        )
        word_histogram[len(words)] += 1
        if len(words) == 4:
            four_word_fibres.append((signature, words))
    if dict(word_histogram) != EXPECTED_JOINT_WORD_HISTOGRAM:
        fail(f"joint response-word histogram changed: {word_histogram}")
    if len(four_word_fibres) != 1:
        fail(f"four-word fibre count changed: {len(four_word_fibres)}")
    four_signature, four_words = four_word_fibres[0]
    expected_four_words = {
        (0, 0, 0, 0): 23,
        (0, 1, 0, 1): 9,
        (1, 0, 0, 1): 21,
        (1, 1, 0, 1): 9,
    }
    if signature_indices(four_signature) != (107,) or dict(four_words) != expected_four_words:
        fail("the unique four-word hostile fibre changed")
    print(
        "JOINT_RESPONSE_WORD_CARDINALITY "
        + " ".join(
            f"{word_count}:{fibre_count}"
            for word_count, fibre_count in sorted(word_histogram.items())
        )
    )
    print(
        "FOUR_WORD_FIBRE SIGNATURE 107 "
        + " ".join(
            f"{word_text(word)}:{count}" for word, count in sorted(four_words.items())
        )
    )

    for pair in TOP_ROWS:
        if pair not in signature_of:
            fail(f"top row left the base universe: {pair}")
        signature = signature_of[pair]
        own_word = tuple(int(pair in target) for target in target_sets)
        if own_word != (0, 0, 0, 0):
            fail(f"top row unexpectedly belongs to a target: {pair} {own_word}")
        response_classes: dict[Word, list[Pair]] = defaultdict(list)
        for fibre_pair in fibres[signature]:
            word = tuple(int(fibre_pair in target) for target in target_sets)
            response_classes[word].append(fibre_pair)
        nonzero_classes = {
            word: pairs
            for word, pairs in response_classes.items()
            if word != (0, 0, 0, 0)
        }
        response_text = (
            "|".join(
                f"{word_text(word)}:{format_pairs(pairs)}"
                for word, pairs in sorted(nonzero_classes.items())
            )
            or "-"
        )
        print(
            f"TOP {pair[0]},{pair[1]} SIGNATURE "
            f"{','.join(map(str, signature_indices(signature)))} "
            f"FIBRE_SIZE {len(fibres[signature])} NONZERO_RESPONSE_PEERS "
            f"{response_text}"
        )

    print("VERDICT PASS OLD_SIGNATURE_DOES_NOT_FACTOR_CURRENT_RESPONSES")
    print("SCOPE NO_NEW_ROW_REMOVAL_NO_PHYSICAL_ENTRY_NO_LRC14")


if __name__ == "__main__":
    main()
