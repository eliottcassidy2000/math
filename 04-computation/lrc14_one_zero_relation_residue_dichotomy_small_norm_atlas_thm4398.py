#!/usr/bin/env python3
"""Exact primary verifier for the THM-4398 one-zero relation atlas.

The script imports only the audited exact carrier primitives of THM-4393.
It proves the mod-three relation/owner dichotomy by a complete finite-field
truth table, derives a one-live-class affine roof, closes every coefficient
pattern of l1 norm at most fourteen by an exact integral tail, and checks the
remaining finite windows against the definition-level physical comb.
"""

from fractions import Fraction
from itertools import permutations, product
import importlib.util
from math import gcd
from pathlib import Path
import sys


Q = Fraction
R = Q(3, 14)
TARGET = Q(6, 77)
CHECKS = 0


class VerificationError(RuntimeError):
    pass


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise VerificationError(message)


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "04-computation" / "lrc14_minimal_ternary_unit_norm18_shell_thm4393.py"
SPEC = importlib.util.spec_from_file_location("thm4393_base", BASE_PATH)
BASE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)


EXPECTED = {
    (1, 2, 3): (Q(8, 245), ((1, 5, 7),), {-1: Q(4, 245), 1: Q(4, 245)}),
    (1, 1, 6): (Q(6, 77), ((1, 5, 11),), {-1: Q(3, 77), 1: Q(3, 77)}),
    (1, 3, 4): (Q(6, 77), ((1, 5, 11),), {-1: Q(3, 77), 1: Q(3, 77)}),
    (2, 3, 5): (
        Q(6, 77),
        ((1, 5, 11),),
        {-2: Q(0), -1: Q(3, 77), 1: Q(3, 77), 2: Q(0)},
    ),
    (1, 2, 9): (
        Q(46, 665),
        ((5, 7, 19),),
        {-2: Q(8, 665), -1: Q(3, 133), 1: Q(3, 133), 2: Q(8, 665)},
    ),
    (1, 3, 8): (
        Q(58, 833),
        ((5, 7, 17),),
        {-2: Q(8, 833), -1: Q(3, 119), 1: Q(3, 119), 2: Q(8, 833)},
    ),
    (1, 5, 6): (
        Q(58, 833),
        ((5, 7, 17),),
        {-2: Q(8, 833), -1: Q(3, 119), 1: Q(3, 119), 2: Q(8, 833)},
    ),
    (2, 3, 7): (
        Q(6, 77),
        ((1, 5, 11),),
        {-2: Q(0), -1: Q(3, 77), 1: Q(3, 77), 2: Q(0)},
    ),
    (3, 4, 5): (
        Q(6, 91),
        ((1, 7, 13),),
        {-2: Q(0), -1: Q(3, 91), 1: Q(3, 91), 2: Q(0)},
    ),
    (1, 1, 12): (
        Q(12, 161),
        ((1, 11, 23),),
        {-2: Q(3, 161), -1: Q(3, 161), 1: Q(3, 161), 2: Q(3, 161)},
    ),
    (1, 3, 10): (
        Q(12, 161),
        ((1, 11, 23),),
        {-2: Q(3, 161), -1: Q(3, 161), 1: Q(3, 161), 2: Q(3, 161)},
    ),
    (1, 4, 9): (
        Q(6, 77),
        ((1, 5, 11),),
        {-2: Q(3, 77), -1: Q(0), 1: Q(0), 2: Q(3, 77)},
    ),
    (1, 6, 7): (
        Q(46, 665),
        ((5, 7, 19),),
        {-2: Q(8, 665), -1: Q(3, 133), 1: Q(3, 133), 2: Q(8, 665)},
    ),
    (3, 4, 7): (
        Q(8, 119),
        ((5, 11, 17),),
        {-2: Q(1, 119), -1: Q(3, 119), 1: Q(3, 119), 2: Q(1, 119)},
    ),
}


def gcd_many(values):
    result = 0
    for value in values:
        result = gcd(result, abs(value))
    return result


def fraction_text(value):
    return f"{value.numerator}/{value.denominator}"


def one_zero_patterns(max_norm):
    result = []
    for norm in range(4, max_norm + 1):
        if norm % 2:
            continue
        for first in range(1, norm):
            for second in range(first, norm):
                third = norm - first - second
                if second > third:
                    continue
                pattern = (first, second, third)
                if gcd_many(pattern) != 1:
                    continue
                if sum(value % 3 == 0 for value in pattern) != 1:
                    continue
                result.append(pattern)
    return tuple(result)


def cross(first, second):
    return (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )


def dot(first, second):
    return sum(x * y for x, y in zip(first, second))


def residue_dichotomy_truth_table():
    """Return counts for the complete F_3 relation/owner gate."""
    type_counts = {"all_unit": 0, "one_zero": 0}
    owner_rows = {"all_unit": 0, "one_zero": 0}
    for w in product((1, 2), repeat=3):
        for c in product(range(3), repeat=3):
            if c == (0, 0, 0) or dot(c, w) % 3:
                continue
            zero_count = sum(value == 0 for value in c)
            require(zero_count in (0, 1), ("relation-residue-dichotomy", w, c))
            kind = "all_unit" if zero_count == 0 else "one_zero"
            weighted = tuple(c[i] * w[i] % 3 for i in range(3))
            if kind == "all_unit":
                require(len(set(weighted)) == 1 and weighted[0] != 0,
                        ("all-unit-weighted-residues", w, c, weighted))
            else:
                zero_index = c.index(0)
                other = tuple(i for i in range(3) if i != zero_index)
                require(weighted[zero_index] == 0, ("weighted-zero", w, c))
                require(weighted[other[0]] == -weighted[other[1]] % 3,
                        ("opposite-weighted-residues", w, c, weighted))
            type_counts[kind] += 1

            for owners in permutations((0, 1, 2)):
                n = tuple((-w[i] * owners[i]) % 3 for i in range(3))
                carrier = cross(w, n)
                require(all(value % 3 for value in carrier),
                        ("distinct-owner-carrier-gate", w, owners, carrier))
                delta = dot(c, n) % 3
                live = tuple(
                    k for k in range(3)
                    if all((carrier[i] + k * c[i]) % 3 for i in range(3))
                )
                if kind == "all_unit":
                    require(delta == 0, ("all-unit-zero-defect", w, c, owners, delta))
                    require(len(live) == 2, ("all-unit-two-live-classes", w, c, live))
                else:
                    require(delta != 0, ("one-zero-nonzero-defect", w, c, owners, delta))
                    require(len(live) == 1, ("one-zero-one-live-class", w, c, live))
                    zero_index = c.index(0)
                    first, second = tuple(i for i in range(3) if i != zero_index)
                    q = c[first] * w[first] % 3
                    orientation = (owners[second] - owners[first]) % 3
                    require(delta == q * orientation % 3,
                            ("defect-orientation-bit", w, c, owners, delta))
                owner_rows[kind] += 1
    require(type_counts == {"all_unit": 16, "one_zero": 48},
            ("residue-type-counts", type_counts))
    require(owner_rows == {"all_unit": 96, "one_zero": 288},
            ("owner-row-counts", owner_rows))
    return type_counts, owner_rows


def defect_set(relation):
    """All nonzero-mod-three integer defects in the strict error window."""
    norm = sum(abs(value) for value in relation)
    maximum = (3 * norm - 1) // 14
    return tuple(delta for delta in range(-maximum, maximum + 1) if delta % 3)


def relation_components(w, relation):
    """Raw carriers in the one-zero affine chart, with (delta,k) sidecars."""
    section = BASE.bezout_vector(relation)
    components = {}
    metadata = {}
    for delta in defect_set(relation):
        carrier_zero = BASE.cross(w, tuple(delta * value for value in section))
        lower, upper = BASE.integer_line_support(w, relation, carrier_zero)
        live_residues = tuple(
            k for k in range(3)
            if all((carrier_zero[i] + k * relation[i]) % 3 for i in range(3))
        )
        require(len(live_residues) == 1,
                ("one-live-affine-class", w, relation, delta, live_residues))
        for k in range(lower, upper + 1):
            carrier = tuple(carrier_zero[i] + k * relation[i] for i in range(3))
            if any(value % 3 == 0 for value in carrier):
                continue
            length = BASE.component_length(w, carrier)
            if not length:
                continue
            require(carrier not in components,
                    ("unique-affine-address", w, relation, carrier))
            components[carrier] = length
            metadata[carrier] = (delta, k)
    return components, metadata


def analyze_pattern(pattern):
    maximum, expected_winners, expected_masses = EXPECTED[pattern]
    deltas = defect_set(pattern)
    integrals = {abs(delta): BASE.slice_integral(pattern, delta) for delta in deltas}
    bulk = sum((BASE.slice_integral(pattern, delta) / 3 for delta in deltas), Q(0))
    error_numerator = Q(len(deltas) * 3, 7)
    require(maximum > bulk, ("positive-tail-gap", pattern, maximum, bulk))
    threshold = error_numerator / (maximum - bulk)
    # Every unscanned integral maximum W is at least proof_height+1, which is
    # strictly larger than the exact threshold (or one larger if integral).
    proof_height = threshold.numerator // threshold.denominator
    require(bulk + error_numerator / (proof_height + 1) < maximum,
            ("strict-tail-cutoff", pattern, proof_height))

    rows = BASE.generated_relation_map(pattern, proof_height)
    require(rows, ("nonempty-sector-window", pattern, proof_height))
    values = {}
    carrier_count = 0
    relation_count = 0
    for w, relations in rows.items():
        literal, _ = BASE.literal_physical_components(w)
        reference = None
        for relation in relations:
            relation_count += 1
            require(sum(abs(value) for value in relation) == sum(pattern),
                    ("relation-norm", pattern, relation))
            require(sum(value % 3 == 0 for value in relation) == 1,
                    ("relation-one-zero", pattern, relation))
            raw, metadata = relation_components(w, relation)
            require(raw == literal, ("definition-level-control", pattern, w, relation))
            require(set(delta for delta, _ in metadata.values()).issubset(set(defect_set(relation))),
                    ("defect-support", pattern, w, relation))
            value = sum(raw.values(), Q(0))
            if reference is None:
                reference = value
            else:
                require(value == reference,
                        ("presentation-independent-mass", pattern, w, relation))
            carrier_count += len(raw)
        values[w] = reference

    observed_maximum = max(values.values())
    winners = tuple(sorted(w for w, value in values.items() if value == observed_maximum))
    require(observed_maximum == maximum,
            ("sector-maximum", pattern, observed_maximum, maximum))
    require(winners == expected_winners, ("sector-winners", pattern, winners))

    winner = winners[0]
    winner_relation = next(iter(rows[winner]))
    raw, metadata = relation_components(winner, winner_relation)
    masses = {
        delta: sum((length for carrier, length in raw.items()
                    if metadata[carrier][0] == delta), Q(0))
        for delta in defect_set(winner_relation)
    }
    require(masses == expected_masses, ("winner-defect-masses", pattern, masses))
    return {
        "pattern": pattern,
        "maximum": maximum,
        "winner": winner,
        "masses": masses,
        "integrals": integrals,
        "bulk": bulk,
        "threshold": threshold,
        "proof_height": proof_height,
        "triples": len(rows),
        "relations": relation_count,
        "carriers": carrier_count,
    }


def mapping_text(mapping):
    return "{" + ",".join(f"{key}:{fraction_text(mapping[key])}" for key in sorted(mapping)) + "}"


def main():
    sys.stdout.reconfigure(newline="\n")
    expected_patterns = tuple(EXPECTED)
    patterns = one_zero_patterns(14)
    require(patterns == expected_patterns, ("complete-pattern-list", patterns))
    require(all(sum(pattern) % 2 == 0 for pattern in patterns), "even norms")
    type_counts, owner_rows = residue_dichotomy_truth_table()

    records = tuple(analyze_pattern(pattern) for pattern in patterns)
    equality = tuple(
        (record["pattern"], record["winner"])
        for record in records if record["maximum"] == TARGET
    )
    expected_equality = tuple(
        (pattern, (1, 5, 11))
        for pattern in ((1, 1, 6), (1, 3, 4), (2, 3, 5), (2, 3, 7), (1, 4, 9))
    )
    require(equality == expected_equality, ("five-equality-presentations", equality))
    require(all(record["maximum"] <= TARGET for record in records), "all sectors close")

    print("LRC THM4398 ONE-ZERO RELATION DICHOTOMY AND NORM-AT-MOST-14 ATLAS")
    print("status=PASS VERIFIED_EXACT; scope=one_zero_mod3_relation_presentations_l1<=14; LRC14=OPEN")
    print("residue_types=all_unit_or_exactly_one_zero; weighted_types=all_equal_or_zero_plus_opposites")
    print("owner_gate=all_unit:defect0,two_live_k_classes; one_zero:defect_nonzero,one_live_k_class")
    print(f"truth_table_relation_types={type_counts}; owner_rows={owner_rows}")
    print("patterns=" + str(patterns))
    for record in records:
        print(
            "pattern=" + str(record["pattern"])
            + "; maximum=" + fraction_text(record["maximum"])
            + "; unique_winner=" + str(record["winner"])
            + "; defect_masses=" + mapping_text(record["masses"])
            + "; slice_integrals=" + mapping_text(record["integrals"])
            + "; bulk=" + fraction_text(record["bulk"])
            + "; threshold=" + fraction_text(record["threshold"])
            + "; proof_height=" + str(record["proof_height"])
            + "; triples=" + str(record["triples"])
            + "; relations=" + str(record["relations"])
            + "; relation_chart_carrier_incidences=" + str(record["carriers"])
        )
    print("global_maximum=6/77")
    print("equality_presentations=" + str(tuple(pattern for pattern, _ in equality)))
    print("equality_physical_comb=(1,5,11); presentation_overlap=yes; sector_values_must_not_be_added")
    print("physical_controls=definition_level_for_every_triple_in_every_proof_window")
    print("tail=one_coset_unimodal_sum<=slice_integral/3+3/(7*max_speed)")
    print("inherited_thm4393_checks=" + str(BASE.CHECKS))
    print("local_checks=" + str(CHECKS))
    print("optimization_safe_checks=yes")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
