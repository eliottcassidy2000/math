"""Exact bounded charge decompositions and retained-resource hostile controls.

No third-party packages; all validation remains active under python -O.
--write regenerates the frozen receipt; default mode compares it exactly.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


class Group:
    def __init__(self, moduli):
        self.moduli = tuple(moduli)
        self.zero = (0,) * len(self.moduli)
        self.elements = tuple(product(*(range(m) for m in self.moduli)))
        self.D = 1 + sum(m - 1 for m in self.moduli)

    def add(self, a, b):
        return tuple((x + y) % m for x, y, m in zip(a, b, self.moduli))

    def sub(self, a, b):
        return tuple((x - y) % m for x, y, m in zip(a, b, self.moduli))

    def total(self, sequence):
        answer = self.zero
        for g in sequence:
            answer = self.add(answer, g)
        return answer


def shortest_zero(group, sequence):
    """Minimum-cardinality nonempty zero-sum subset, represented by bitmask."""
    best = {}
    for i, g in enumerate(sequence):
        additions = [(g, 1 << i)]
        additions.extend((group.add(a, g), mask | (1 << i))
                         for a, mask in tuple(best.items()))
        for a, mask in additions:
            if a not in best or (mask.bit_count(), mask) < (best[a].bit_count(), best[a]):
                best[a] = mask
    return best.get(group.zero)


def literal_zero(group, sequence):
    for size in range(1, len(sequence) + 1):
        for indices in combinations(range(len(sequence)), size):
            if group.total(sequence[i] for i in indices) == group.zero:
                return size
    return None


def paired_partition(group, positive, negative):
    if len(positive) != len(negative):
        raise ValueError("Unequal integer margins")
    if group.total(positive) != group.total(negative):
        raise ValueError("Unequal charge")
    remaining = list(range(len(positive)))
    blocks = []
    while remaining:
        window = remaining[:group.D]
        differences = [group.sub(positive[i], negative[i]) for i in window]
        mask = shortest_zero(group, differences)
        require(mask is not None, "Davenport witness missing")
        block = [window[i] for i in range(len(window)) if mask & (1 << i)]
        require(0 < len(block) <= group.D, "Invalid block size")
        require(group.total(positive[i] for i in block) ==
                group.total(negative[i] for i in block), "Block charge debt")
        blocks.append(block)
        used = set(block)
        remaining = [i for i in remaining if i not in used]
    require(sorted(i for block in blocks for i in block) == list(range(len(positive))),
            "Partition lost or recycled an index")
    return blocks


def sharp_sequence(group):
    sequence = []
    for i, modulus in enumerate(group.moduli):
        unit = tuple(int(j == i) for j in range(len(group.moduli)))
        sequence.extend([unit] * (modulus - 1))
    sequence.append(tuple(1 for _ in group.moduli))
    return sequence


def test_small_sequences():
    rows = []
    # Universes: ALL ordered D-term sequences avoiding zero in these groups.
    # Sequences containing zero have the immediate singleton witness.
    for moduli in [(2,), (3,), (4,), (2, 2), (4, 2), (3, 3)]:
        group = Group(moduli)
        nonzero = [x for x in group.elements if x != group.zero]
        count = 0
        max_shortest = 0
        for sequence in product(nonzero, repeat=group.D):
            mask = shortest_zero(group, sequence)
            require(mask is not None, "Exhaustive p-group bound failed")
            require(group.total(sequence[i] for i in range(group.D)
                                if mask & (1 << i)) == group.zero,
                    "Subset witness has nonzero charge")
            # Independent literal subset path, every sequence in stated universe.
            size = literal_zero(group, sequence)
            require(size == mask.bit_count(), "DP/literal disagreement")
            max_shortest = max(max_shortest, size)
            count += 1
        rows.append({"moduli": moduli, "D": group.D,
                     "ordered_nonzero_sequences": count,
                     "maximum_minimum_zero_subset": max_shortest})
    return rows


def test_balanced_supplies():
    rows = []
    # ALL ordered positive and negative label lists of length <=3.
    for moduli in [(2,), (3,), (4,), (2, 2)]:
        group = Group(moduli)
        histogram = Counter()
        count = 0
        for m in range(1, 4):
            labels = list(product(group.elements, repeat=m))
            for positive in labels:
                for negative in labels:
                    if group.total(positive) != group.total(negative):
                        continue
                    blocks = paired_partition(group, positive, negative)
                    changed = set()
                    for block in blocks:
                        changed.update(block)
                        require(group.total(positive[i] for i in changed) ==
                                group.total(negative[i] for i in changed),
                                "Intermediate current changed")
                        # Actual names distinguish 2m slots, each with two
                        # private resources. Charge coincidences merge no labels.
                        resources = [(side, i, int(i in changed))
                                     for side in (0, 1) for i in range(m)]
                        require(len(set(resources)) == 2*m, "Private resource collision")
                        require(2*len(block) <= min(2*m, 2*group.D), "Budget exceeded")
                    histogram[str(max(2*len(b) for b in blocks))] += 1
                    count += 1
        rows.append({"moduli": moduli, "balanced_ordered_endpoint_pairs": count,
                     "max_block_support_histogram": dict(sorted(histogram.items()))})
    return rows


def test_sharp_supplies():
    rows = []
    for moduli in [(), (2,), (3,), (4,), (2, 2), (4, 2),
                   (4, 2, 2, 2, 2), (128, 2, 2, 2, 2)]:
        group = Group(moduli)
        sequence = sharp_sequence(group)
        require(len(sequence) == group.D, "Wrong extremal sequence length")
        require(group.total(sequence) == group.zero, "Extremal sequence not zero")
        mask = shortest_zero(group, sequence)
        require(mask == (1 << group.D) - 1, "Extremal sequence has proper zero subset")
        blocks = paired_partition(group, sequence, [group.zero]*group.D)
        require(len(blocks) == 1 and len(blocks[0]) == group.D,
                "Extremal balanced supply split")
        rows.append({"moduli": moduli, "group_size": len(group.elements),
                     "D": group.D, "forced_changed_slots": 2*group.D,
                     "minimal_zero_sequence": sequence})
    return rows


def test_vacancy_paths():
    rows = []
    for m in range(2, 15):
        q = [0] + [1]*(m-1) + [0]
        old_charge = sum(q[:m]) % 2
        injective = []
        guarded = []
        balanced_guarded = []
        for bits in product((0, 1), repeat=m):
            resources = [i + bits[i] for i in range(m)]
            literal = len(set(resources)) == m
            implications = all(bits[i] <= bits[i+1] for i in range(m-1))
            require(literal == implications, "Resource/implication disagreement")
            if literal:
                injective.append(bits)
                new_charge = sum(q[j] for j in resources) % 2
                delta = sum(bits[i]*(q[i+1]-q[i]) for i in range(m)) % 2
                require((new_charge-old_charge) % 2 == delta, "Non-telescoping current")
                if new_charge == old_charge:
                    guarded.append(bits)
                if m % 2 == 0:
                    target_margin = sum((1 if i % 2 == 0 else -1)*bits[i]
                                        for i in range(m))
                    source_charge = (bits[0]+bits[-1]) % 2
                    require(source_charge == delta,
                            "Balanced source code/resource charge disagreement")
                    if target_margin == 0 and source_charge == 0:
                        balanced_guarded.append(bits)
        require(len(injective) == m+1, "Vacancy path has wrong suffix count")
        require(guarded == [(0,)*m, (1,)*m], "Unexpected guarded intermediate")
        # Independent resources with identical charge changes admit first/last
        # together and each middle alone. These are genuine complete endpoints.
        charge_moves = [(0, m-1)] + [(i,) for i in range(1, m-1)]
        changed = set()
        for move in charge_moves:
            changed.update(move)
            require(sum(q[i+1]-q[i] for i in changed) % 2 == 0,
                    "Independent-resource positive control fails")
        require(len(changed) == m, "Positive control lost a slot")
        if m % 2 == 0:
            require(balanced_guarded == [(0,)*m, (1,)*m],
                    "Unexpected balanced U/V intermediate")
            pair_moves = [(0, m-1)] + [(i, i+1) for i in range(1, m-1, 2)]
            changed = set()
            for move in pair_moves:
                changed.update(move)
                require(sum(1 if i % 2 == 0 else -1 for i in changed) == 0,
                        "Independent balanced control changes target margin")
                require(sum(i in (0, m-1) for i in changed) % 2 == 0,
                        "Independent balanced control changes source charge")
            require(len(changed) == m, "Balanced positive control lost a slot")
        rows.append({"slots": m, "binary_assignments": 2**m,
                     "injective_mixtures": len(injective), "guarded_states": 2,
                     "forced_budget": m, "independent_resource_budget": 2,
                     "balanced_UV_variant_checked": m % 2 == 0,
                     "balanced_UV_guarded_states": len(balanced_guarded) if m % 2 == 0 else None})
    return rows


def dependencies():
    paths = [
        "05-knowledge/results/session45_injective_multisource_moves_20261023/"
        + stem + "/THEOREM.md"
        for stem in ["dyadic-constraints", "injective-completion", "balanced-circuits"]
    ]
    code_path = "05-knowledge/results/injective_native_blocks_20261023/electrical/receipt.json"
    paths.append(code_path)
    result = {}
    for name in paths:
        data = (ROOT/name).read_bytes().replace(b"\r\n", b"\n")
        result[name] = sha256(data).hexdigest()
    receipt = json.loads((ROOT/code_path).read_text(encoding="utf-8"))
    codes = receipt["codes"]
    require(len(codes) == 3528, "Wrong original merged-owner universe")
    require(len(set(tuple(row) for row in codes)) == 64, "Wrong effective mod8 code image")
    require(codes[1965] == codes[2234], "Selected target paired codes differ")
    return {"sha256_lf_normalized": result,
            "native_mod8_target_codes": {str(i): codes[i] for i in (1965, 2234)},
            "code_image_size": 64,
            "scope": "Code and theorem authentication only; no new native supply is inferred."}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = {"status": "FINITE-EXACT supporting PROVED scoped theorem",
              "dependencies": dependencies(),
              "p_group_sequences": test_small_sequences(),
              "balanced_supplies": test_balanced_supplies(),
              "sharp_abstract_supplies": test_sharp_supplies(),
              "vacancy_path_hostiles": test_vacancy_paths(),
              "scope": "Ordered small universes as recorded; abstract full-record controls; no broad native enumeration."}
    result["checks_before_replay"] = CHECKS
    result = json.loads(json.dumps(result))
    path = HERE/"receipt.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    else:
        require(result == json.loads(path.read_text(encoding="utf-8")), "Frozen receipt differs")
    print(f"PASS {CHECKS} checks; {'wrote' if args.write else 'replayed'} {path.name}")


if __name__ == "__main__":
    main()
