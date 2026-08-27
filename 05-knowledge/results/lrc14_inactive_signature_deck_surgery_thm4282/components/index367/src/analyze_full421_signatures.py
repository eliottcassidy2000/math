#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
from collections import Counter, defaultdict
from pathlib import Path


OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def fnv_words(rows) -> str:
    state = OFFSET
    for row in rows:
        for word in row:
            for shift in range(0, 64, 8):
                state ^= (word >> shift) & 0xFF
                state = state * PRIME & MASK64
    return f"{state:016x}"


def read_deck(path: Path) -> list[int]:
    deck = [int(line.strip(), 16) for line in path.read_text().splitlines() if line.strip()]
    assert len(deck) == 421 and len(set(deck)) == 421
    assert all(mask.bit_count() == 8 and mask < (1 << 30) for mask in deck)
    return deck


def read_rows(path: Path):
    edges: list[tuple[int, int]] = []
    signatures: list[int] = []
    weights: list[int] = []
    with path.open(newline="") as handle:
        rows = csv.reader(handle)
        assert next(rows) == ["q", "r", "inactive_count", "w0", "w1", "w2", "w3", "w4", "w5", "w6"]
        for fields in rows:
            assert len(fields) == 10
            q, r, weight = map(int, fields[:3])
            words = [int(value, 16) for value in fields[3:]]
            signature = sum(word << (64 * index) for index, word in enumerate(words))
            assert signature >> 421 == 0
            assert signature.bit_count() == weight > 0
            edges.append((q, r))
            signatures.append(signature)
            weights.append(weight)
    assert len(edges) == 24223 and edges == sorted(set(edges))
    return edges, signatures, weights


def mask_label(index: int, deck: list[int]) -> str:
    return f"{index}:{deck[index]:08x}"


def sig_indices(signature: int) -> list[int]:
    out = []
    while signature:
        bit = (signature & -signature).bit_length() - 1
        out.append(bit)
        signature &= signature - 1
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("signatures", type=Path)
    parser.add_argument("deck", type=Path)
    args = parser.parse_args()
    deck = read_deck(args.deck)
    edges, signatures, weights = read_rows(args.signatures)
    rows = len(edges)
    universe = (1 << rows) - 1

    columns = [0] * 421
    for row, signature in enumerate(signatures):
        for index in sig_indices(signature):
            columns[index] |= 1 << row

    print("INPUT_SHA256", hashlib.sha256(args.signatures.read_bytes()).hexdigest())
    print("ROWS", rows, "DISTINCT_SIGNATURES", len(set(signatures)),
          "INACTIVE_INCIDENCES", sum(weights))

    singleton_counts = Counter(
        signature.bit_length() - 1
        for signature in signatures if signature.bit_count() == 1
    )
    forced = sorted(singleton_counts)
    forced_cover = 0
    for index in forced:
        forced_cover |= columns[index]
    forced_ledger = fnv_words(
        (index, deck[index], singleton_counts[index]) for index in forced
    )
    uncovered_rows = [row for row in range(rows)
                      if ((forced_cover >> row) & 1) == 0]
    extension_masks = [index for index, column in enumerate(columns)
                       if all((column >> row) & 1 for row in uncovered_rows)]
    assert len(forced) == 215 and sum(singleton_counts.values()) == 5038
    assert uncovered_rows == [edges.index((70, 302))]
    assert sig_indices(signatures[uncovered_rows[0]]) == [78, 368]
    assert extension_masks == [78, 368]
    assert all((forced_cover | columns[index]) == universe
               for index in extension_masks)
    print("FORCED_SINGLETON_CORE", len(forced), "SINGLETON_ROWS",
          sum(singleton_counts.values()), "LEDGER_FNV", forced_ledger,
          "COVERED_ROWS", forced_cover.bit_count(), "UNCOVERED_ROWS",
          len(uncovered_rows))
    print("FORCED_SINGLETON_MASKS",
          " ".join(mask_label(index, deck) for index in forced))
    print("SOLE_UNCOVERED_ROW 70,302 INACTIVE",
          " ".join(mask_label(index, deck) for index in extension_masks))
    print("EXACT_MINIMUM_HOSTILE_BASES SIZE 216 COUNT 2 EXTENSIONS",
          " ".join(mask_label(index, deck) for index in extension_masks),
          "LOWER_BOUND 215_FORCED_PLUS_ONE_UNCOVERED_ROW")

    weight_distribution = Counter(weights)
    print("WEIGHT_RANGE", min(weights), max(weights), "MEAN_NUM", sum(weights), "MEAN_DEN", rows)
    print("WEIGHT_DISTRIBUTION")
    for weight, count in sorted(weight_distribution.items()):
        print(weight, count)

    fibres = Counter(signatures)
    fibre_sizes = Counter(fibres.values())
    print("FIBRE_SIZE_DISTRIBUTION size,n_fibres,n_rows")
    for size, count in sorted(fibre_sizes.items()):
        print(size, count, size * count)
    print("LARGEST_FIBRES")
    for signature, count in sorted(fibres.items(), key=lambda item: (-item[1], item[0]))[:20]:
        members = [edges[i] for i, value in enumerate(signatures) if value == signature]
        print(count, signature.bit_count(), ";".join(f"{q},{r}" for q, r in members[:12]),
              "..." if len(members) > 12 else "")

    counts = [column.bit_count() for column in columns]
    print("MASK_INACTIVE_RANGE", min(counts), max(counts), "ZERO", counts.count(0), "ALL", counts.count(rows))
    print("MASKS_MOST_INACTIVE")
    for index in sorted(range(421), key=lambda i: (-counts[i], i))[:30]:
        print(mask_label(index, deck), counts[index])
    print("MASKS_LEAST_NONZERO_INACTIVE")
    for index in sorted((i for i in range(421) if counts[i]), key=lambda i: (counts[i], i))[:30]:
        print(mask_label(index, deck), counts[index])
    print("MASKS_ZERO_INACTIVE", " ".join(mask_label(i, deck) for i, count in enumerate(counts) if count == 0))

    equality_groups: dict[int, list[int]] = defaultdict(list)
    for index, column in enumerate(columns):
        equality_groups[column].append(index)
    nontrivial_equalities = [group for group in equality_groups.values() if len(group) > 1]
    print("COLUMN_EQUIVALENCE_CLASSES", len(equality_groups), "NONTRIVIAL", len(nontrivial_equalities))
    for group in sorted(nontrivial_equalities, key=lambda g: (-len(g), g)):
        print("EQUIV", columns[group[0]].bit_count(), " ".join(mask_label(i, deck) for i in group))

    implications = []
    for i, left in enumerate(columns):
        for j, right in enumerate(columns):
            if i != j and left & ~right == 0:
                implications.append((i, j))
    strict_implications = [(i, j) for i, j in implications if columns[i] != columns[j]]
    print("DIRECTED_IMPLICATIONS", len(implications), "STRICT", len(strict_implications))
    outgoing = Counter(i for i, _ in strict_implications)
    incoming = Counter(j for _, j in strict_implications)
    print("TOP_STRICT_IMPLICATION_SOURCES")
    for i in sorted(range(421), key=lambda x: (-outgoing[x], x))[:20]:
        if outgoing[i]:
            targets = [j for a, j in strict_implications if a == i]
            print(mask_label(i, deck), counts[i], outgoing[i], "TARGETS",
                  " ".join(mask_label(j, deck) for j in targets[:20]))
    print("TOP_STRICT_IMPLICATION_TARGETS")
    for j in sorted(range(421), key=lambda x: (-incoming[x], x))[:20]:
        if incoming[j]:
            sources = [i for i, b in strict_implications if b == j]
            print(mask_label(j, deck), counts[j], incoming[j], "SOURCES",
                  " ".join(mask_label(i, deck) for i in sources[:20]))

    pair_stats = []
    for i in range(421):
        for j in range(i + 1, 421):
            intersection = (columns[i] & columns[j]).bit_count()
            union = counts[i] + counts[j] - intersection
            pair_stats.append((intersection, union, i, j))
    print("TOP_COINACTIVE_INTERSECTION")
    for intersection, union, i, j in sorted(pair_stats, key=lambda x: (-x[0], x[2], x[3]))[:30]:
        print(mask_label(i, deck), mask_label(j, deck), "INTER", intersection, "UNION", union,
              "LEFT", counts[i], "RIGHT", counts[j])
    print("TOP_JACCARD_NONIDENTICAL")
    candidates = [x for x in pair_stats if x[1] and columns[x[2]] != columns[x[3]]]
    for intersection, union, i, j in sorted(candidates,
            key=lambda x: (-(x[0] / x[1]), -x[0], x[2], x[3]))[:30]:
        print(mask_label(i, deck), mask_label(j, deck), "INTER", intersection, "UNION", union,
              "JACCARD_NUM_DEN", intersection, union)

    edge_to_row = {edge: i for i, edge in enumerate(edges)}
    top_edges = [(256, 663), (366, 663), (520, 663)]
    print("TOP663_SIGNATURES_AND_NEAREST")
    for edge in top_edges:
        row = edge_to_row[edge]
        signature = signatures[row]
        members = [edges[i] for i, value in enumerate(signatures) if value == signature]
        distances = [(signature ^ value).bit_count() for value in signatures]
        positive = [distance for index, distance in enumerate(distances) if index != row]
        minimum = min(positive)
        nearest = [edges[index] for index, distance in enumerate(distances)
                   if index != row and distance == minimum]
        print("TOP", f"{edge[0]},{edge[1]}", "WEIGHT", signature.bit_count(),
              "FIBRE", len(members), "MIN_DISTANCE", minimum,
              "NEAREST_COUNT", len(nearest), "NEAREST",
              ";".join(f"{q},{r}" for q, r in nearest[:40]))
        print("INACTIVE", " ".join(mask_label(i, deck) for i in sig_indices(signature)))
        for peer in nearest[:10]:
            peer_sig = signatures[edge_to_row[peer]]
            removed = signature & ~peer_sig
            added = peer_sig & ~signature
            print("DELTA", f"{peer[0]},{peer[1]}", "REMOVED",
                  " ".join(mask_label(i, deck) for i in sig_indices(removed)),
                  "ADDED", " ".join(mask_label(i, deck) for i in sig_indices(added)))

if __name__ == "__main__":
    main()
