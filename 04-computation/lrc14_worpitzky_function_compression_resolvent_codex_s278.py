#!/usr/bin/env python3
"""HYP-3151 scout: Worpitzky/function compression bridge for LRC14.

This executable synthesis merges four nearby lanes:

* HYP-3150: the reserved function-compression factor-through wall;
* HYP-3147: the n=3 two-class edge-flip kernel and Worpitzky descent word;
* HYP-3148/HYP-3149: the n=4 live-core/canary/filler quotient tables;
* HYP-3132/HYP-3142: the k=8 quartic/bounded-core resolvent;
* the prompt's function split: a+b, a*b are unordered; a^b, b^a are ordered.

The guiding rule tested here is:

  a compression is proof-legal for a target function only if the target is
  constant on compression fibers, or a named ordered/canary sidecar reconstructs
  the destroyed coordinate.

This is not an LRC14 proof.  It is a finite exact scout for the next packet
schema.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb
from typing import Callable, Dict, Iterable, List, Sequence, Tuple


Bit3 = Tuple[int, int, int]
Bit2 = Tuple[int, int]

N3_CLASSES = {
    (0, 0, 0): "C",
    (1, 1, 1): "C",
}

N4_CLASS_BY_BITS: Dict[Bit3, str] = {
    (0, 0, 0): "T",
    (1, 0, 0): "+",
    (0, 1, 0): "-",
    (1, 1, 0): "S",
    (0, 0, 1): "S",
    (1, 0, 1): "S",
    (0, 1, 1): "S",
    (1, 1, 1): "S",
}

N4_TARGET_CLASS: Dict[Bit2, str] = {
    (0, 0): "T",
    (1, 0): "+",
    (0, 1): "-",
    (1, 1): "S",
}


def n3_class(bits: Bit3) -> str:
    return N3_CLASSES.get(bits, "T")


def flip(bits: Bit3, idx: int) -> Bit3:
    out = list(bits)
    out[idx] ^= 1
    return tuple(out)  # type: ignore[return-value]


def n3_edge_kernel() -> Dict[str, object]:
    counts: Counter[Tuple[str, str]] = Counter()
    state_rows = []
    for bits in product((0, 1), repeat=3):
        bits3 = tuple(bits)  # type: ignore[assignment]
        cls = n3_class(bits3)
        exits = []
        for idx in range(3):
            nxt = flip(bits3, idx)
            nxt_cls = n3_class(nxt)
            counts[(cls, nxt_cls)] += 1
            exits.append((idx, nxt, nxt_cls))
        state_rows.append((bits3, cls, exits))
    row_C = (Fraction(counts[("C", "C")], 6), Fraction(counts[("C", "T")], 6))
    row_T = (Fraction(counts[("T", "C")], 18), Fraction(counts[("T", "T")], 18))
    return {
        "counts": dict(sorted(counts.items())),
        "matrix_rows_C_T": (row_C, row_T),
        "stationary_C_T": (Fraction(1, 4), Fraction(3, 4)),
        "eigenvalues": (Fraction(1, 1), Fraction(-1, 3)),
        "state_rows": state_rows,
    }


def descents(order: Sequence[int]) -> int:
    return sum(1 for a, b in zip(order, order[1:]) if a > b)


def eulerian_number(n: int, k: int) -> int:
    # A(n,k): permutations of n items with k descents.
    total = 0
    for j in range(k + 1):
        total += (-1) ** j * comb(n + 1, j) * (k + 1 - j) ** n
    return total


def eulerian_row(n: int) -> Tuple[int, ...]:
    return tuple(eulerian_number(n, k) for k in range(n))


def worpitzky_check(n: int, limit: int = 8) -> List[Tuple[int, int, int]]:
    row = eulerian_row(n)
    out = []
    for x in range(limit + 1):
        lhs = x**n
        rhs = sum(row[k] * comb(x + n - 1 - k, n) for k in range(n))
        out.append((x, lhs, rhs))
    return out


def function_swap_audit(max_value: int = 7) -> Dict[str, object]:
    pairs = [(a, b) for a in range(1, max_value + 1) for b in range(1, max_value + 1) if a != b]

    def invariant_count(fn: Callable[[int, int], object]) -> int:
        return sum(1 for a, b in pairs if fn(a, b) == fn(b, a))

    ordered_examples = []
    for a, b in pairs:
        if a == 1 or b == 1:
            continue
        if a**b != b**a:
            ordered_examples.append((a, b, a**b, b**a))
        if len(ordered_examples) == 6:
            break

    return {
        "ordered_pair_count": len(pairs),
        "symmetric_invariant_counts": {
            "a+b": invariant_count(lambda a, b: a + b),
            "a*b": invariant_count(lambda a, b: a * b),
        },
        "ordered_invariant_counts": {
            "a^b": invariant_count(lambda a, b: a**b),
            "ordered_pair(a^b,b^a)": invariant_count(lambda a, b: (a**b, b**a)),
        },
        "ordered_examples": ordered_examples,
    }


def or_compression(bits: Bit3) -> Bit2:
    a, b, c = bits
    return (a | c, b | c)


def n4_compression_audit() -> Dict[str, object]:
    rows = []
    fiber_by_xy: Dict[Bit2, List[Bit3]] = {(0, 0): [], (1, 0): [], (0, 1): [], (1, 1): []}
    for bits in product((0, 1), repeat=3):
        bits3 = tuple(bits)  # type: ignore[assignment]
        xy = or_compression(bits3)
        cls = N4_CLASS_BY_BITS[bits3]
        target = N4_TARGET_CLASS[xy]
        rows.append((bits3, cls, xy, target, cls == target))
        fiber_by_xy[xy].append(bits3)

    affine_count = 0
    affine_examples = []
    for offset in product((0, 1), repeat=2):
        for matrix_bits in product((0, 1), repeat=6):
            rows2 = (matrix_bits[:3], matrix_bits[3:])

            def affine(bits: Bit3) -> Bit2:
                vals = []
                for r in rows2:
                    vals.append(offset[len(vals)] ^ ((r[0] & bits[0]) ^ (r[1] & bits[1]) ^ (r[2] & bits[2])))
                return tuple(vals)  # type: ignore[return-value]

            if all(N4_CLASS_BY_BITS[bits] == N4_TARGET_CLASS.get(affine(bits), "?") for bits in N4_CLASS_BY_BITS):
                affine_count += 1
                if len(affine_examples) < 3:
                    affine_examples.append((offset, rows2))

    return {
        "rows": rows,
        "fiber_by_xy": {str(k): v for k, v in fiber_by_xy.items()},
        "class_preserving": all(row[-1] for row in rows),
        "affine_class_preserving_count": affine_count,
        "affine_examples": affine_examples,
    }


def quartic_resolvent_audit() -> Dict[str, object]:
    # g(t)=(t-1)(t-2)(t-4)(t-5).  Return expanded coefficients in t and in u=t-3.
    roots_t = (1, 2, 4, 5)
    # Coefficients high-to-low.
    coeff = [1]
    for r in roots_t:
        coeff = [coeff[0]] + [coeff[i] - r * coeff[i - 1] for i in range(1, len(coeff))] + [-r * coeff[-1]]

    # In centered coordinate u, roots are -2,-1,1,2.
    roots_u = tuple(r - 3 for r in roots_t)
    centered = [1]
    for r in roots_u:
        centered = [centered[0]] + [centered[i] - r * centered[i - 1] for i in range(1, len(centered))] + [-r * centered[-1]]
    # centered should be u^4 + 0u^3 -5u^2 + 0u +4.
    disc_u2 = 25 - 16
    return {
        "roots_t": roots_t,
        "coefficients_t_high_to_low": tuple(coeff),
        "roots_u": roots_u,
        "coefficients_u_high_to_low": tuple(centered),
        "odd_centered_coefficients": (centered[1], centered[3]),
        "u2_quadratic_discriminant": disc_u2,
        "degree_ceiling": 4,
        "below_abel_ruffini_quintic_wall": True,
        "worpitzky_degree4_row": eulerian_row(4),
    }


@dataclass(frozen=True)
class Carrier:
    name: str
    function_legality: int
    compression_legality: int
    orientation_payload: int
    canary_payload: int
    resolvent_payload: int
    lrc_transfer: int

    def score(self) -> int:
        return (
            3 * self.function_legality
            + 3 * self.compression_legality
            + 2 * self.orientation_payload
            + 2 * self.canary_payload
            + 2 * self.resolvent_payload
            + 3 * self.lrc_transfer
        )


def tournament_fingerprint(carriers: Sequence[Carrier]) -> Dict[str, object]:
    def beats(a: Carrier, b: Carrier) -> bool:
        if a.score() != b.score():
            return a.score() > b.score()
        return a.name < b.name

    directed_3cycles = 0
    for a, b, c in combinations(carriers, 3):
        ab, bc, ca = beats(a, b), beats(b, c), beats(c, a)
        if (ab and bc and ca) or ((not ab) and (not bc) and (not ca)):
            directed_3cycles += 1

    hpaths = 0
    for order in permutations(carriers):
        if all(beats(order[i], order[i + 1]) for i in range(len(order) - 1)):
            hpaths += 1

    ordered = sorted(carriers, key=lambda c: (-c.score(), c.name))
    return {
        "score_hist": dict(sorted(Counter(c.score() for c in carriers).items())),
        "directed_3cycles": directed_3cycles,
        "hamiltonian_path_count": hpaths,
        "selected_path": [c.name for c in ordered],
    }


def proof_carriers() -> Dict[str, object]:
    carriers = (
        Carrier("function_legal_packet_schema", 5, 5, 5, 5, 4, 5),
        Carrier("nonlinear_or_canary_compression", 5, 5, 3, 5, 3, 5),
        Carrier("n3_edge_flip_minus_one_third_kernel", 4, 4, 5, 2, 2, 4),
        Carrier("worpitzky_descent_function_basis", 5, 4, 4, 2, 4, 4),
        Carrier("k8_biquadratic_resolvent_ceiling", 4, 4, 3, 3, 5, 5),
        Carrier("ordered_exponential_sidecar", 4, 3, 5, 2, 1, 4),
        Carrier("symmetric_sum_product_shadow", 2, 2, 1, 1, 1, 2),
        Carrier("raw_score_class_compression", 1, 1, 1, 1, 1, 1),
    )
    return tournament_fingerprint(carriers)


def main() -> None:
    print("HYP-3151 / codex-2026-06-28-S278")
    print("Worpitzky/function compression bridge for LRC14")
    print("namespace=HYP-3151/T1216/LTI-277/LTT-175")
    print()

    print("1. N=3 TWO-CLASS EDGE-FLIP KERNEL")
    kernel = n3_edge_kernel()
    print(f"counts={kernel['counts']}")
    print(f"matrix_rows_C_T={kernel['matrix_rows_C_T']}")
    print(f"stationary_C_T={kernel['stationary_C_T']}")
    print(f"eigenvalues={kernel['eigenvalues']}")
    print()

    print("2. WORPITZKY DESCENT BASIS")
    for n in (3, 4):
        print(f"degree_{n}_eulerian_row={eulerian_row(n)}")
        rows = worpitzky_check(n, 6)
        print(f"degree_{n}_identity_verified={all(lhs == rhs for _, lhs, rhs in rows)}")
        print(f"degree_{n}_sample_rows={rows[:5]}")
    print()

    print("3. FUNCTION ORDER AUDIT")
    fs = function_swap_audit()
    print(f"ordered_pair_count={fs['ordered_pair_count']}")
    print(f"symmetric_invariant_counts={fs['symmetric_invariant_counts']}")
    print(f"ordered_invariant_counts={fs['ordered_invariant_counts']}")
    print(f"ordered_examples={fs['ordered_examples']}")
    print("rule=symmetric functions survive unordered quotient; ordered functions need orientation sidecar")
    print()

    print("4. N=4 CANARY OR-COMPRESSION")
    comp = n4_compression_audit()
    print(f"class_preserving_OR={comp['class_preserving']}")
    print(f"affine_class_preserving_count={comp['affine_class_preserving_count']}")
    print("rows=(a,b,c),class,OR(x,y),target_class,ok")
    for row in comp["rows"]:  # type: ignore[index]
        print(f"  {row}")
    print(f"fiber_by_xy={comp['fiber_by_xy']}")
    print("rule=the useful compression is nonlinear OR, not a GF(2) linear quotient")
    print()

    print("5. K=8 DEGREE-4 RESOLVENT CEILING")
    res = quartic_resolvent_audit()
    for key, value in res.items():
        print(f"{key}={value}")
    print("rule=bounded-core dual reaches degree 4 at k=8, still below the Abel-Ruffini quintic wall")
    print()

    print("6. TOURNAMENT ANALYSIS OVER PROOF CARRIERS")
    fp = proof_carriers()
    for key, value in fp.items():
        print(f"{key}={value}")
    print()

    print("7. PROPOSED PACKET FIELDS")
    fields = [
        "target_function_id",
        "function_swap_parity",
        "symmetric_shadow_status",
        "ordered_sidecar_status",
        "edge_flip_kernel_mode",
        "minority_edge_gate",
        "worpitzky_descent_word",
        "compression_map_kind",
        "compression_fiber_function_constancy",
        "canary_or_restoration_sidecar",
        "resolvent_degree",
        "centered_odd_coefficient_status",
        "abel_ruffini_wall_status",
        "terminal_exit_or_named_debt",
    ]
    for field in fields:
        print(f"  {field}")


if __name__ == "__main__":
    main()
