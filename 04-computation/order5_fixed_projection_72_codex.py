#!/usr/bin/env python3
"""
Order-5 fixed-projection gate for a putative extremal Type II [72,36,16] code.

If an automorphism of order 5 has cycle shape 5^14 1^2, the fixed subcode
projects to a binary self-dual length-16 code on the fourteen 5-cycles plus the
two fixed coordinates.  Doubly-evenness makes that projection Type II.  There
are only two Type II [16,8] codes: e8+e8 and d16+.

This script marks two coordinates as the original fixed coordinates and assigns
weighted length

    wt_72(w) = 5 * (# marked-cycle coordinates in w) + (# fixed coordinates in w).

An extremal [72,36,16] code cannot have any nonzero fixed word of wt_72 < 16.
Thus a tetrad in the projected Type II [16,8] code containing both fixed marks
is fatal: it would lift to weight 5+5+1+1 = 12.

Tournament Analysis note:
  We do not assume the useful vertices are runners/arcs.  Candidate vertex sets
  considered: order-5 automorphism cases, projected Type II code types, marked
  coordinate-pair classes, tetrads, Fano heptads, low-weight leakage events,
  nonfixed eigenspaces, Fourier modes over F_16, and proof obligations.  The
  quotient used here takes vertices to marked-pair branch classes, preserving
  the fixed-code low-weight obstruction and the fixed-A16 congruence while
  destroying all information about the nonfixed eigenspaces and glue.
  Challenged assumption: "small automorphism" need not mean "runner-like";
  the first hard datum is a support-marking obstruction in a tiny Type II code.
"""

from __future__ import annotations

import itertools
from collections import Counter, defaultdict


N72_A16 = 249_849
N72_A16_MOD5 = N72_A16 % 5


def wt(x: int) -> int:
    return x.bit_count()


def span(gens: list[int]) -> list[int]:
    words = [0]
    for g in gens:
        words += [w ^ g for w in words]
    return sorted(words)


def hamming_e8_words(shift: int = 0) -> list[int]:
    """Extended Hamming [8,4,4] words, shifted into coordinates shift..shift+7."""
    rows_bits = [
        [1, 0, 0, 0, 0, 1, 1, 1],
        [0, 1, 0, 0, 1, 0, 1, 1],
        [0, 0, 1, 0, 1, 1, 0, 1],
        [0, 0, 0, 1, 1, 1, 1, 0],
    ]
    gens = [
        sum(bit << (shift + i) for i, bit in enumerate(row))
        for row in rows_bits
    ]
    return span(gens)


def e8_plus_e8_words() -> list[int]:
    left = hamming_e8_words(0)
    right = hamming_e8_words(8)
    return sorted(a ^ b for a in left for b in right)


def d16_plus_words() -> list[int]:
    """d16+ as even sums of 11-pairs plus the 1010...10 glue word."""
    pair = lambda i: 0b11 << (2 * i)
    gens = [pair(i) ^ pair(i + 1) for i in range(7)]
    glue = sum(1 << (2 * i) for i in range(8))
    return span(gens + [glue])


def coordinates(mask: int) -> list[int]:
    return [i for i in range(16) if (mask >> i) & 1]


def pair_class_e8(pair: tuple[int, int]) -> str:
    a, b = pair
    return "same-e8-block" if (a // 8) == (b // 8) else "split-e8-block"


def pair_class_d16(pair: tuple[int, int]) -> str:
    a, b = pair
    return "same-dplus-pair" if (a // 2) == (b // 2) else "different-dplus-pairs"


def classify(code_name: str, pair: tuple[int, int]) -> str:
    if code_name == "e8+e8":
        return pair_class_e8(pair)
    if code_name == "d16+":
        return pair_class_d16(pair)
    raise ValueError(code_name)


def analyze_code(name: str, words: list[int]) -> dict:
    tetrads = [w for w in words if wt(w) == 4]
    pair_cover = Counter()
    pair_witnesses: dict[tuple[int, int], list[int]] = defaultdict(list)
    for t in tetrads:
        cs = coordinates(t)
        for p in itertools.combinations(cs, 2):
            pair_cover[p] += 1
            pair_witnesses[p].append(t)

    pairs = list(itertools.combinations(range(16), 2))
    rows = []
    class_rows = defaultdict(list)
    for pair in pairs:
        fixed_mask = (1 << pair[0]) | (1 << pair[1])
        lifted = Counter()
        min_nonzero = None
        for w in words:
            fixed_hits = wt(w & fixed_mask)
            cycle_hits = wt(w) - fixed_hits
            lifted_weight = 5 * cycle_hits + fixed_hits
            lifted[lifted_weight] += 1
            if w and (min_nonzero is None or lifted_weight < min_nonzero):
                min_nonzero = lifted_weight

        tetrads_both = pair_cover[pair]
        tetrads_exactly_one = sum(
            1 for t in tetrads if wt(t & fixed_mask) == 1
        )
        valid = min_nonzero is not None and min_nonzero >= 16
        row = {
            "pair": pair,
            "class": classify(name, pair),
            "pair_cover": tetrads_both,
            "min_nonzero": min_nonzero,
            "valid": valid,
            "fixed_A16": lifted[16],
            "fixed_A16_mod5_ok": lifted[16] % 5 == N72_A16_MOD5,
            "lifted_distribution": dict(sorted(lifted.items())),
        }
        rows.append(row)
        class_rows[row["class"]].append(row)

    return {
        "name": name,
        "words": words,
        "weight_distribution": Counter(wt(w) for w in words),
        "tetrads": tetrads,
        "pair_cover": pair_cover,
        "rows": rows,
        "class_rows": dict(class_rows),
    }


def summarize_code(data: dict) -> list[str]:
    out = []
    out.append(f"== {data['name']} ==")
    out.append(f"dimension check: |C|={len(data['words'])} = 2^8")
    out.append(f"weight distribution: {dict(sorted(data['weight_distribution'].items()))}")
    out.append(f"tetrads A4={len(data['tetrads'])}")
    covered_pairs = sum(1 for p in itertools.combinations(range(16), 2) if data["pair_cover"][p])
    out.append(f"covered coordinate pairs by tetrads: {covered_pairs}/120")
    for cls, rows in sorted(data["class_rows"].items()):
        valid = [r for r in rows if r["valid"]]
        minvals = Counter(r["min_nonzero"] for r in rows)
        a16vals = Counter(r["fixed_A16"] for r in rows)
        covervals = Counter(r["pair_cover"] for r in rows)
        out.append(
            f"  class {cls}: pairs={len(rows)}, valid={len(valid)}, "
            f"pair-cover={dict(sorted(covervals.items()))}, "
            f"min-nonzero={dict(sorted(minvals.items()))}, "
            f"fixed-A16={dict(sorted(a16vals.items()))}"
        )
        if valid:
            sample = valid[0]
            out.append(
                f"    sample valid pair {sample['pair']}: lifted distribution "
                f"{sample['lifted_distribution']}"
            )
    fatal_examples = [r for r in data["rows"] if not r["valid"]]
    if fatal_examples:
        sample = fatal_examples[0]
        out.append(
            f"  sample fatal pair {sample['pair']} has min lift {sample['min_nonzero']} "
            f"and tetrad-pair cover {sample['pair_cover']}"
        )
    return out


def branch_rows(e8: dict, d16: dict) -> list[dict]:
    branches = []
    for data in (e8, d16):
        for cls, rows in data["class_rows"].items():
            valid_count = sum(r["valid"] for r in rows)
            min_lift = min(r["min_nonzero"] for r in rows if r["min_nonzero"] is not None)
            fixed_a16_values = sorted({r["fixed_A16"] for r in rows})
            congruent_count = sum(
                r["valid"] and r["fixed_A16_mod5_ok"] for r in rows
            )
            total_weight12 = sum(r["lifted_distribution"].get(12, 0) for r in rows)
            branches.append(
                {
                    "vertex": f"{data['name']}::{cls}",
                    "pairs": len(rows),
                    "valid_pairs": valid_count,
                    "min_lift": min_lift,
                    "fixed_A16_values": tuple(fixed_a16_values),
                    "congruent_valid_pairs": congruent_count,
                    "weight12_total": total_weight12,
                }
            )
    return sorted(branches, key=lambda r: r["vertex"])


def dominates(a: dict, b: dict) -> int:
    """Return +1 if a beats b under the tie Hamiltonian path, else -1."""
    key_a = (
        a["valid_pairs"] > 0,
        a["congruent_valid_pairs"] > 0,
        a["min_lift"],
        -a["weight12_total"],
        a["valid_pairs"],
        -len(a["fixed_A16_values"]),
        a["vertex"],
    )
    key_b = (
        b["valid_pairs"] > 0,
        b["congruent_valid_pairs"] > 0,
        b["min_lift"],
        -b["weight12_total"],
        b["valid_pairs"],
        -len(b["fixed_A16_values"]),
        b["vertex"],
    )
    return 1 if key_a > key_b else -1


def tournament_summary(branches: list[dict]) -> list[str]:
    score = Counter({b["vertex"]: 0 for b in branches})
    edges = {}
    for i, a in enumerate(branches):
        for b in branches[i + 1 :]:
            if dominates(a, b) > 0:
                winner, loser = a["vertex"], b["vertex"]
            else:
                winner, loser = b["vertex"], a["vertex"]
            score[winner] += 1
            edges[(winner, loser)] = 1

    vertices = [b["vertex"] for b in branches]
    c3 = 0
    for triple in itertools.combinations(vertices, 3):
        cyclic = (
            (triple[0], triple[1]) in edges
            and (triple[1], triple[2]) in edges
            and (triple[2], triple[0]) in edges
        ) or (
            (triple[0], triple[2]) in edges
            and (triple[2], triple[1]) in edges
            and (triple[1], triple[0]) in edges
        )
        c3 += int(cyclic)

    ordered = sorted(vertices, key=lambda v: (-score[v], v))
    out = []
    out.append("== Tournament Analysis: marked-pair branch quotient ==")
    out.append(
        "observable: valid weighted-min>=16 branch wins; then A16 mod-5 match; "
        "then larger min lift; then fewer weight-12 leaks; then more valid markings"
    )
    out.append("tie Hamiltonian path: lexicographic vertex name after the numeric keys")
    for b in branches:
        out.append(
            f"  {b['vertex']}: pairs={b['pairs']}, valid={b['valid_pairs']}, "
            f"min_lift={b['min_lift']}, fixed_A16_values={b['fixed_A16_values']}, "
            f"congruent_valid={b['congruent_valid_pairs']}, "
            f"weight12_total={b['weight12_total']}"
        )
    out.append(f"score histogram: {dict(sorted(Counter(score.values()).items()))}")
    out.append(f"directed 3-cycles: {c3}")
    out.append(f"Hamiltonian path by score: {' > '.join(ordered)}")
    out.append(
        "SCCs: singleton under this transitive quotient; the loss of cycles is a warning "
        "that the next TA should add nonfixed-eigenspace/coset-glue observables."
    )
    return out


def main() -> None:
    e8 = analyze_code("e8+e8", e8_plus_e8_words())
    d16 = analyze_code("d16+", d16_plus_words())

    print("Order-5 fixed projection gate for the putative Type II [72,36,16] code")
    print(f"A16(W72)={N72_A16}, so fixed minimum words must be congruent to {N72_A16_MOD5} mod 5")
    print("Assumed cycle shape: 5^14 1^2; fixed projection length = 14+2 = 16")
    print("Extremal obstruction: no nonzero fixed word may lift to weighted length < 16")
    print()

    for line in summarize_code(e8):
        print(line)
    print()
    for line in summarize_code(d16):
        print(line)
    print()

    print("== Derived order-5 gate ==")
    print(
        "e8+e8 survives exactly when the two fixed coordinates lie in different e8 blocks: "
        "64 valid marked pairs, each with fixed_A16=14 congruent to W72 A16 mod 5."
    )
    print(
        "d16+ has every coordinate pair covered by a tetrad, so every marking creates a "
        "fixed word of lifted weight 12; d16+ is excluded from an extremal order-5 fixed projection."
    )
    print(
        "Therefore any order-5 automorphism of an extremal [72,36,16] code with shape 5^14 1^2 "
        "forces a partition of the fourteen 5-cycles into two heptads, one heptad plus one "
        "fixed point in each e8 block."
    )
    print(
        f"The fixed weight-16 words are exactly 14, leaving {N72_A16 - 14} nonfixed minimum "
        f"words = {(N72_A16 - 14) // 5} orbits of size 5."
    )
    print()

    for line in tournament_summary(branch_rows(e8, d16)):
        print(line)


if __name__ == "__main__":
    main()
