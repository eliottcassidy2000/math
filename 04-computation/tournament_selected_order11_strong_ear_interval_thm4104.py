#!/usr/bin/env python3
"""Primary exact referee for the provisional THM-4104 selected order-11 image.

The construction reproduces THM-4102's deterministic first labelled order-10
witness for each of its 7,566 selected values, then evaluates all 2^10-2
nonconstant ears above those parents.  It is a selected lower image, not a
complete census of strong tournaments of order eleven.

The fast evaluator is the exact quadratic cut recurrence obtained from the
THM-4097 (Start, End, Q) boundary formula.  It is cross-checked against that
original formula on complete lower-order banks and frozen order-10 controls.
Every executable gate uses ``require`` and therefore survives ``python -O``.
"""

from __future__ import annotations

import importlib.util
import json
from hashlib import sha256
from math import isqrt, prod
from pathlib import Path


HERE = Path(__file__).resolve().parent
THM4097_COMPILER = HERE / "tournament_order9_strong_ear_spectrum_thm4097.py"
THM4102_COMPILER = HERE / "tournament_selected_order10_strong_ear_interval_thm4102.py"

THM4097_COMPILER_SHA256 = (
    "610ca5850b272e0e75c574f2c1a710a0b96c75cc7191b1e1f1a03dfbdd1378d6"
)
THM4102_COMPILER_SHA256 = (
    "e7049a6347100e6dc54c7b6c03b299cde7dcfaca811954797971b8e7552421a8"
)
THM4102_SELECTED_ORDER9_BANK_SHA256 = (
    "c03c203943e734d09bee4b8818227b8f184405ce4c5092dd56d0fdb6107d528c"
)
SELECTED_ORDER10_BANK_SHA256 = (
    "2f3fbd5d7f56de24a1f08ea08585dd029c70344ef444830915b5ea0d203e4b92"
)
SELECTED_ORDER10_SOURCE_SHA256 = (
    "c4b073e7e0ba7d965bb978635a2cd4ec60dc8acfb815d4388b32bab644d71980"
)
SELECTED_ORDER11_VALUES_SHA256 = (
    "3bead7d8a2539f3c217540199e4b71b208fed526f80cd43bf435e17dd1b0c328"
)
SEMANTIC_SHA256 = "11241b1a8d55d9a1b725b2343b0cf8543397a48a2eb9ba102bb3132b8b020067"

PRIMARY_INTERVAL = (429, 80_265, 39_919)
SECONDARY_INTERVAL = (80_875, 84_259, 1_693)
KEY_VALUES = (225, 429, 14_657, 14_777, 80_265, 80_875, 84_259, 93_751)
FULL_FORMULA_CONTROL_PARENTS = (
    125,
    249,
    14_649,
    14_653,
    14_655,
    15_055,
    15_551,
    15_621,
)
SAMPLED_FORMULA_SIGNATURES = (1, 2, 255, 511, 767, 1_022)

EXPECTED_KEY_ROWS: dict[int, tuple[int, int, int]] = {
    225: (125, 32, 34_095_368_048_213_567),
    429: (125, 116, 33_813_343_248_580_159),
    14_657: (773, 898, 4_980_304_472_833_599),
    14_777: (697, 283, 25_105_481_806_249_279),
    80_265: (13_253, 903, 2_181_427_721_976_895),
    80_875: (15_059, 182, 33_617_973_800_943_135),
    84_259: (12_667, 15, 36_024_248_539_118_623),
    93_751: (15_581, 841, 3_645_901_060_717_615),
}


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"FAILED: {label}")


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_thm4097():
    require(file_sha256(THM4097_COMPILER) == THM4097_COMPILER_SHA256,
            "THM-4097 compiler hash")
    require(file_sha256(THM4102_COMPILER) == THM4102_COMPILER_SHA256,
            "THM-4102 compiler hash")
    spec = importlib.util.spec_from_file_location("thm4097_for_thm4104", THM4097_COMPILER)
    require(spec is not None and spec.loader is not None, "THM-4097 import specification")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def all_ear_values(
    state: tuple[list[int], list[int], list[list[int]]]
) -> list[int]:
    """Evaluate every cut by the exact quadratic recurrence.

    For a cut mask M, the original boundary polynomial is

      sum_(b in M) Start(b) + sum_(a notin M) End(a)
        + sum_(a notin M,b in M) Q(a,b).

    Starting from the empty mask, toggling v into M adds a fixed linear term
    and subtracts Q(u,v)+Q(v,u) for every previously toggled u.
    """
    starts, ends, exposed = state
    order = len(starts)
    base = sum(ends)
    linear = [
        starts[vertex] - ends[vertex]
        + sum(exposed[other][vertex] for other in range(order))
        for vertex in range(order)
    ]
    values = [0] * (1 << order)
    values[0] = base
    for mask in range(1, 1 << order):
        bit = mask & -mask
        vertex = bit.bit_length() - 1
        old = mask ^ bit
        penalty = 0
        rest = old
        while rest:
            other_bit = rest & -rest
            other = other_bit.bit_length() - 1
            penalty += exposed[other][vertex] + exposed[vertex][other]
            rest ^= other_bit
        values[mask] = values[old] + linear[vertex] - penalty
    return values


def canonical_digest(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return sha256(payload).hexdigest()


def is_prime(value: int) -> bool:
    return value >= 2 and all(value % divisor for divisor in range(2, isqrt(value) + 1))


def prime_factors(value: int) -> list[int]:
    factors: list[int] = []
    divisor = 3
    while divisor * divisor <= value:
        while value % divisor == 0:
            factors.append(divisor)
            value //= divisor
        divisor += 2
    if value > 1:
        factors.append(value)
    return factors


def first_missing_prime(values: set[int], start: int) -> int:
    candidate = start if start % 2 else start + 1
    while True:
        if is_prime(candidate) and candidate not in values:
            return candidate
        candidate += 2


def first_missing_seven_prime(values: set[int], start_prime: int) -> int:
    candidate = start_prime if start_prime % 2 else start_prime + 1
    while True:
        if is_prime(candidate) and 7 * candidate not in values:
            return candidate
        candidate += 2


def thm4094_factor_atoms(value: int) -> list[int]:
    """Return the exact factor atoms used in THM-4094 Theorem 6.1."""
    if value == 1:
        return []
    seven_order = 0
    residual = value
    while residual % 7 == 0:
        seven_order += 1
        residual //= 7

    if seven_order == 0:
        return prime_factors(residual)
    if seven_order % 2 == 0:
        return [49] * (seven_order // 2) + prime_factors(residual)
    if seven_order >= 3:
        return [343] + [49] * ((seven_order - 3) // 2) + prime_factors(residual)

    residual_factors = prime_factors(residual)
    if all(factor == 3 for factor in residual_factors):
        require(len(residual_factors) >= 2, "excluded 7 or 21 entered prefix factorization")
        return [63] + residual_factors[2:]

    exceptional_prime = next(factor for factor in residual_factors if factor >= 5)
    remainder = list(residual_factors)
    remainder.remove(exceptional_prime)
    return [7 * exceptional_prime] + remainder


def main() -> None:
    parent = load_thm4097()
    engine = parent.load_engine()
    representatives, counts = engine.generate(8)
    require(counts == parent.A000568, "A000568 class counts through order eight")

    # Rebuild the exact deterministic order-nine bank inherited by THM-4102.
    selected9: dict[int, list[int]] = {}
    selected9_codes: dict[int, int] = {}
    strong_order8 = 0
    order8_ears = 0
    order8_formula_checks = 0
    for old in representatives[8]:
        if not engine.is_strong(8, old):
            continue
        strong_order8 += 1
        state = parent.boundary_state(old)
        fast = all_ear_values(state)
        for signature in range(1, (1 << 8) - 1):
            order8_ears += 1
            original = parent.insertion_h(state, signature)
            require(fast[signature] == original, "quadratic/original formula at order eight")
            order8_formula_checks += 1
            value = fast[signature]
            if value not in selected9:
                child = engine.extend(old, 8, signature)
                selected9[value] = child
                selected9_codes[value] = parent.raw_code(child)

    require(strong_order8 == 6_008, "strong order-eight parent count")
    require(order8_ears == 1_526_032, "complete order-eight ear bank")
    historical, all_classes, strong_classes = parent.read_historical_strong_values()
    require(set(selected9) == historical, "selected order-nine bank equals historical spectrum")
    require((all_classes, strong_classes) == (191_536, 178_133),
            "historical order-nine class controls")
    require(len(selected9) == 1_482, "selected order-nine bank size")
    selected9_bank = sorted(selected9_codes.items())
    selected9_digest = canonical_digest(selected9_bank)
    require(selected9_digest == THM4102_SELECTED_ORDER9_BANK_SHA256,
            "THM-4102 selected order-nine parent bank digest")

    # Rebuild THM-4102's first labelled order-ten witness for every selected
    # value.  Full formula equality is checked on every one of its parent ears.
    selected10: dict[int, list[int]] = {}
    selected10_codes: dict[int, int] = {}
    selected10_sources: dict[int, tuple[int, int]] = {}
    order9_ears = 0
    order9_formula_checks = 0
    for parent_h, old in sorted(selected9.items()):
        state = parent.boundary_state(old)
        fast = all_ear_values(state)
        for signature in range(1, (1 << 9) - 1):
            order9_ears += 1
            original = parent.insertion_h(state, signature)
            require(fast[signature] == original, "quadratic/original formula at order nine")
            order9_formula_checks += 1
            value = fast[signature]
            if value not in selected10:
                child = engine.extend(old, 9, signature)
                selected10[value] = child
                selected10_codes[value] = parent.raw_code(child)
                selected10_sources[value] = (parent_h, signature)

    require(order9_ears == 755_820, "complete THM-4102 selected-ear universe")
    require(len(selected10) == 7_566, "THM-4102 selected value count")
    require((min(selected10), max(selected10)) == (125, 15_621),
            "THM-4102 selected extrema")

    selected10_parent_failures = 0
    for value, old in sorted(selected10.items()):
        selected10_parent_failures += (
            engine.Hcount(10, old) != value or not engine.is_strong(10, old)
        )
    require(selected10_parent_failures == 0, "selected order-ten parent DP/strongness")
    selected10_bank = sorted(selected10_codes.items())
    selected10_digest = canonical_digest(selected10_bank)
    require(selected10_digest == SELECTED_ORDER10_BANK_SHA256,
            "selected order-ten parent bank digest")
    selected10_source_digest = canonical_digest(
        [(value, *selected10_sources[value]) for value in sorted(selected10_sources)]
    )
    require(selected10_source_digest == SELECTED_ORDER10_SOURCE_SHA256,
            "selected order-ten source digest")

    # Evaluate all 7,566*(2^10-2) selected order-eleven ears.  The original
    # formula is checked at six spread masks for every parent and at every mask
    # for eight inherited boundary parents.  Strongness of every construction
    # follows from the proved nonconstant-ear lemma.
    values: set[int] = set()
    first_witness: dict[int, tuple[int, int, int]] = {}
    order10_ears = 0
    order10_formula_checks = 0
    full_control_set = set(FULL_FORMULA_CONTROL_PARENTS)
    sample_set = set(SAMPLED_FORMULA_SIGNATURES)
    for parent_h, old in sorted(selected10.items()):
        state = parent.boundary_state(old)
        fast = all_ear_values(state)
        control_signatures = (
            range(1, (1 << 10) - 1)
            if parent_h in full_control_set
            else SAMPLED_FORMULA_SIGNATURES
        )
        for signature in control_signatures:
            original = parent.insertion_h(state, signature)
            require(fast[signature] == original, "quadratic/original formula at order ten")
            order10_formula_checks += 1
        for signature in range(1, (1 << 10) - 1):
            order10_ears += 1
            value = fast[signature]
            values.add(value)
            if value not in first_witness:
                child = engine.extend(old, 10, signature)
                first_witness[value] = (parent_h, signature, parent.raw_code(child))

    expected_order10_formula_checks = (
        len(selected10) * len(sample_set)
        + len(full_control_set) * ((1 << 10) - 2 - len(sample_set))
    )
    require(order10_formula_checks == expected_order10_formula_checks == 53_524,
            "order-ten original-formula control census")
    require(order10_ears == 7_566 * 1_022 == 7_732_452,
            "complete selected order-eleven ear universe")
    require(len(values) == 43_251, "selected order-eleven value count")
    require((min(values), max(values)) == (225, 93_751), "selected image extrema")
    require(all(value % 2 == 1 for value in values), "Redei parity in selected image")
    selected11_values_digest = canonical_digest(sorted(values))
    require(selected11_values_digest == SELECTED_ORDER11_VALUES_SHA256,
            "selected order-eleven sorted-value digest")

    intervals = parent.odd_intervals(values)
    require(PRIMARY_INTERVAL in intervals, "primary solid interval")
    require(SECONDARY_INTERVAL in intervals, "secondary solid interval")
    require(427 not in values and 80_267 not in values,
            "primary interval adjacent selected-image holes")
    require(80_873 not in values and 84_261 not in values,
            "secondary interval adjacent selected-image holes")

    key_rows: dict[str, dict[str, int]] = {}
    key_failures = 0
    for value in KEY_VALUES:
        require(value in first_witness, f"retained key witness H={value}")
        parent_h, signature, code = first_witness[value]
        child = parent.decode_raw(code, 11)
        direct_h = engine.Hcount(11, child)
        strong = engine.is_strong(11, child)
        key_failures += direct_h != value or not strong
        key_rows[str(value)] = {
            "parent_h": parent_h,
            "signature": signature,
            "cut_weight": signature.bit_count(),
            "code": code,
            "dp_H": direct_h,
        }
    require(key_failures == 0, "key order-eleven DP/strongness controls")
    for value, expected in EXPECTED_KEY_ROWS.items():
        actual = first_witness[value]
        require(actual == expected, f"stable labelled row H={value}")

    ordinary_missing = first_missing_prime(values, 14_657)
    seven_missing_prime = first_missing_seven_prime(values, 2_111)
    require(ordinary_missing == 80_407, "next ordinary-prime lane target")
    require(seven_missing_prime == 11_527 and 7 * seven_missing_prime == 80_689,
            "next seven-prime lane target")
    previous_ordinary_prime = max(
        value for value in range(3, ordinary_missing, 2) if is_prime(value)
    )
    previous_seven_prime = max(
        value for value in range(3, seven_missing_prime, 2) if is_prime(value)
    )
    require(previous_ordinary_prime == 80_387,
            "largest supplied ordinary prime before next target")
    require(previous_seven_prime == 11_519,
            "largest supplied exceptional prime before next target")

    def ordinary_atom_available(prime: int) -> bool:
        return prime != 7 and (prime <= 14_653 or prime in values)

    def seven_atom_available(prime: int) -> bool:
        return prime != 3 and (prime <= 2_099 or 7 * prime in values)

    prefix_checks = 0
    prefix_direct = 0
    prefix_composed = 0
    composed_values: list[int] = []
    for value in range(1, 80_406, 2):
        if value in (7, 21):
            continue
        atoms = thm4094_factor_atoms(value)
        require(prod(atoms) == value, "THM-4094 factor product")
        for atom in atoms:
            if atom in (49, 63, 343):
                continue
            if atom % 7 == 0:
                require(seven_atom_available(atom // 7), "available seven-prime atom")
            else:
                require(is_prime(atom) and ordinary_atom_available(atom),
                        "available ordinary-prime atom")
        prefix_direct += value <= 14_655 or value in values
        if value > 14_655 and value not in values:
            prefix_composed += 1
            composed_values.append(value)
        prefix_checks += 1
    require(prefix_checks == 40_201, "global allowed-prefix census")
    require((prefix_direct, prefix_composed) == (40_189, 12),
            "direct/composed global-prefix split")
    require(is_prime(80_407) and 80_407 not in values,
            "first ordinary lane value is unforced by selected image")
    require(is_prime(11_527) and 80_689 not in values,
            "first seven-prime lane value is unforced by selected image")

    longest = sorted(intervals, key=lambda row: row[2], reverse=True)[:15]
    ledger = {
        "input_hashes": [THM4097_COMPILER_SHA256, THM4102_COMPILER_SHA256],
        "class_counts": counts,
        "selected9_bank_digest": selected9_digest,
        "selected10_parent_codes": selected10_bank,
        "selected10_sources": sorted(selected10_sources.items()),
        "selected10_bank_digest": selected10_digest,
        "formula_checks": [
            order8_formula_checks,
            order9_formula_checks,
            order10_formula_checks,
        ],
        "order10_formula_control_parents": list(FULL_FORMULA_CONTROL_PARENTS),
        "order10_formula_sample_signatures": list(SAMPLED_FORMULA_SIGNATURES),
        "ear_universes": [order8_ears, order9_ears, order10_ears],
        "selected11_values": sorted(values),
        "intervals": intervals,
        "key_rows": key_rows,
        "lane_next": [ordinary_missing, seven_missing_prime, 7 * seven_missing_prime],
        "global_prefix": [
            80_405,
            prefix_checks,
            prefix_direct,
            prefix_composed,
            composed_values,
        ],
    }
    semantic = canonical_digest(ledger)
    require(semantic == SEMANTIC_SHA256, "canonical semantic ledger digest")

    print("THM-4104 SELECTED ORDER-ELEVEN STRONG-EAR PRIMARY REFEREE")
    print("status=PROVISIONAL_PRIMARY_PASS")
    print("A000568_counts=", [counts[order] for order in range(1, 9)])
    print("strong_order8_parents=", strong_order8,
          "complete_order8_formula_checks=", order8_formula_checks)
    print("selected_order9_parents=", len(selected9),
          "complete_order9_formula_checks=", order9_formula_checks)
    print("selected_order9_bank_sha256=", selected9_digest)
    print("selected_order10_parents=", len(selected10),
          "direct_DP_strong_failures=", selected10_parent_failures)
    print("selected_order10_bank_sha256=", selected10_digest)
    print("selected_order10_source_sha256=", selected10_source_digest)
    print("order10_original_formula_controls=", order10_formula_checks, "failures= 0")
    print("nonconstant_order11_ears_checked=", order10_ears)
    print("selected_value_count=", len(values), "min=", min(values), "max=", max(values))
    print("selected_value_set_sha256=", selected11_values_digest)
    print("primary_solid_interval= [429,80265] count=39919 adjacent_holes=(427,80267)")
    print("secondary_solid_interval= [80875,84259] count=1693 adjacent_holes=(80873,84261)")
    print("first_intervals=", intervals[:15])
    print("longest_intervals=", longest)
    print("direct_key_witness_checks=", len(key_rows), "failures=", key_failures)
    print("key_witnesses=", key_rows)
    print("global_allowed_prefix_through=80405 checks=", prefix_checks,
          "direct=", prefix_direct, "composed=", prefix_composed)
    print("multiplicatively_composed_prefix_values=", composed_values)
    print("ordinary_prime_lane_through=", previous_ordinary_prime,
          "next_unforced=", ordinary_missing)
    print("seven_prime_lane_p_through=", previous_seven_prime,
          "next_unforced_value=", 7 * seven_missing_prime)
    print("semantic_sha256=", semantic)
    print("semantic_source=THM-4102 deterministic selected labelled order-ten bank")
    print("semantic_map=adjoin every nonconstant labelled cut and evaluate the quadratic boundary")
    print("semantic_preserved=Hamiltonian-path count and strongness by the strong-ear lemma")
    print("semantic_destroyed=unselected equal-H parents, isomorphism multiplicity, and full order-eleven image")
    print("semantic_sidecar=labelled parent code plus Start/End/Q boundary and cut signature")
    print("scope=selected lower image only; full order-eleven census and global H-spectrum remain OPEN")
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()
