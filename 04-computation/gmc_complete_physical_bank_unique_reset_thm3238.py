#!/usr/bin/env python3
"""Exact complete physical-bank unique-reset certificate for THM-3238."""

import ast
import hashlib
from collections import Counter
from fractions import Fraction
from functools import lru_cache, reduce
from itertools import combinations_with_replacement
from math import gcd, lcm
from multiprocessing import Pool
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
UPSTREAM = HERE / "gmc_pole_prefix_hasse_current_scout.py"
UPSTREAM_BANK = HERE / "gmc_product_gamma_arbitrary_anchored_schur_thm3110.py"
THM3216_SCRIPT = HERE / "gmc_depth_nine_unique_reset_face_omega_boundary_thm3216.py"
THM3216_OUTPUT = ROOT / "05-knowledge/results/gmc_depth_nine_unique_reset_face_omega_boundary_thm3216.out"
DEPENDENCIES = (
    (UPSTREAM,
     "151edb9b8cee4807d3f08dc17af32e45420021ba30dfd116c04d9fcaf8bbd5b7"),
    (UPSTREAM_BANK,
     "15b94691d53afbcdc6aefda89fc7cd5497534ca70fb780a686892dabb5646d6f"),
    (THM3216_SCRIPT,
     "3951b7e9a5e08199fd03ab2b3827dcbeb8c1039790acd4d218713c78d44ab9cf"),
    (THM3216_OUTPUT,
     "d303a7fb18b07f2091e1fe5c705b15e62ecd2f9049fe6f29a003834288a24636"),
)


def lf_hash(path):
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


for path, expected in DEPENDENCIES:
    require(lf_hash(path) == expected, ("dependency hash drift", str(path)))


def load_upstream_prefix(maximum_degree):
    tree = ast.parse(UPSTREAM.read_text(encoding="utf-8"))
    prefix = []
    for node in tree.body:
        if (isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name)
                        and target.id == "MAXIMUM_DEGREE"
                        for target in node.targets)):
            node.value = ast.Constant(maximum_degree)
        prefix.append(node)
        if (isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name)
                        and target.id == "UNIVERSE"
                        for target in node.targets)):
            break
    module = ast.fix_missing_locations(ast.Module(body=prefix, type_ignores=[]))
    namespace = {"__file__": str(UPSTREAM)}
    exec(compile(module, str(UPSTREAM), "exec"), namespace)
    return namespace


UP = load_upstream_prefix(14)
UP["residual_roots"] = lru_cache(maxsize=None)(UP["residual_roots"])
UP["all_monomial_values"] = lru_cache(maxsize=None)(UP["all_monomial_values"])
coefficient_vectors = UP["coefficient_vectors"]
dominant_row = UP["dominant_row"]
reduced_poles = UP["reduced_poles"]
residual_roots = UP["residual_roots"]
partitions = UP["partitions"]
BANK = UP["BANKS"][1]


@lru_cache(maxsize=None)
def coarsens(fine, coarse):
    if sum(fine) != sum(coarse) or len(fine) < len(coarse):
        return False
    pieces = tuple(sorted(fine, reverse=True))
    bins = tuple(sorted(coarse, reverse=True))

    @lru_cache(maxsize=None)
    def place(index, capacities):
        if index == len(pieces):
            return not any(capacities)
        piece = pieces[index]
        tried = set()
        for slot, capacity in enumerate(capacities):
            if capacity < piece or capacity in tried:
                continue
            tried.add(capacity)
            changed = list(capacities)
            changed[slot] -= piece
            changed.sort(reverse=True)
            if place(index + 1, tuple(changed)):
                return True
        return False

    return place(0, bins)


def ones(count):
    return (1,) * count


NORMALIZERS = {
    5: 1300713120,
    6: 591017995680,
    7: 67528413334656,
    8: 9558239814808320,
    9: 1141773245941342464,
    10: 98177407096144199040,
    11: 6402822693369065077104,
    12: 623312546828577253020720,
    13: 39394150378793250375693600,
    14: 2528571670236939601303479360,
}


# Each entry is degree, upset size, minimal generators, and a positive exact
# factor a_i.  The lawful raw response multiplier is a_i/M_degree.  The first
# eleven rows are inherited from THM-3216, the next three are principal, and
# the final eight are nonprincipal full-bank stitches.
CERTIFICATE = (
    (5, 1, ((5,),), Fraction(7, 916609201)),
    (14, 130, ((3, 2) + ones(9), (2, 2, 2) + ones(8)),
     Fraction(7177018, 925091523)),
    (10, 40, ((3,) + ones(7), (2, 2) + ones(6)),
     Fraction(113613, 103885675)),
    (12, 74, ((2, 2) + ones(8),), Fraction(14821346, 514459311)),
    (14, 128, ((4,) + ones(10), (3, 2) + ones(9),
               (2, 2, 2, 2, 2, 2, 1, 1)),
     Fraction(1800950, 915068561)),
    (14, 132, ((2, 2) + ones(10),), Fraction(234911529, 973260547)),
    (14, 129, ((5,) + ones(9), (3, 3) + ones(8),
               (2, 2, 2) + ones(8)), Fraction(12760257, 713572348)),
    (14, 127, ((5,) + ones(9), (4, 2) + ones(8),
               (3, 3) + ones(8), (2, 2, 2, 2) + ones(6)),
     Fraction(1171444, 999246305)),
    (10, 41, ((2,) + ones(8),), Fraction(1806523, 275215201)),
    (11, 55, ((2,) + ones(9),), Fraction(5625173, 775067384)),
    (11, 54, ((3,) + ones(8), (2, 2) + ones(7)),
     Fraction(3559607, 824858825)),
    (13, 98, ((2, 2) + ones(9),), Fraction(30666900, 415449457)),
    (13, 100, ((2,) + ones(11),), Fraction(7423537, 97767608)),
    (14, 134, ((2,) + ones(12),), Fraction(386487777, 737119657)),
    (5, 5, ((3, 1, 1), (2, 2, 1)), Fraction(2, 599801583)),
    (9, 28, ((3,) + ones(6), (2, 2) + ones(5)),
     Fraction(5467, 645985968)),
    (11, 51, ((3, 2) + ones(6), (2, 2, 2) + ones(5)),
     Fraction(81670, 325987229)),
    (13, 96, ((4,) + ones(9), (3, 2) + ones(8),
              (2, 2, 2, 2) + ones(5)), Fraction(1928501, 649376836)),
    (13, 97, ((4,) + ones(9), (3, 2) + ones(8),
              (2, 2, 2) + ones(7)), Fraction(1131782, 357891057)),
    (14, 122, ((5,) + ones(9), (4, 2) + ones(8),
               (3, 3, 3) + ones(5), (2, 2, 2, 2, 2) + ones(4)),
     Fraction(9291, 416817778)),
    (11, 53, ((3,) + ones(8), (2, 2, 2) + ones(5)),
     Fraction(768511, 977040433)),
    (14, 126, ((5,) + ones(9), (4, 2) + ones(8),
               (3, 3) + ones(8), (3, 2, 2, 2) + ones(5),
               (2, 2, 2, 2, 2) + ones(4)),
     Fraction(452405, 829943856)),
)


# The promoted depth-nine cocircuit is retained only to certify the sharp
# failure of a one-row degree-five repair on the complete bank.
OLD_CERTIFICATE = (
    (5, 1, ((5,),), Fraction(5, 72755559)),
    (14, 130, ((3, 2) + ones(9), (2, 2, 2) + ones(8)),
     Fraction(1577093, 47102691)),
    (8, 21, ((2,) + ones(6),), Fraction(53234, 56508147)),
    (10, 40, ((3,) + ones(7), (2, 2) + ones(6)),
     Fraction(151685, 66318847)),
    (12, 76, ((2,) + ones(10),), Fraction(30203115, 76749488)),
    (12, 74, ((2, 2) + ones(8),), Fraction(57683, 3902927)),
    (14, 128, ((4,) + ones(10), (3, 2) + ones(9),
               (2, 2, 2, 2, 2, 2, 1, 1)), Fraction(283819, 63667394)),
    (14, 121, ((7,) + ones(7), (3, 3, 2) + ones(6),
               (2, 2, 2, 2) + ones(6)), Fraction(2861, 13045678)),
    (14, 132, ((2, 2) + ones(10),), Fraction(10395173, 58298789)),
    (14, 129, ((5,) + ones(9), (3, 3) + ones(8),
               (2, 2, 2) + ones(8)), Fraction(1951835, 98897507)),
    (14, 127, ((5,) + ones(9), (4, 2) + ones(8),
               (3, 3) + ones(8), (2, 2, 2, 2) + ones(6)),
     Fraction(3979, 1773976)),
    (10, 41, ((2,) + ones(8),), Fraction(4381327, 51259840)),
    (11, 55, ((2,) + ones(9),), Fraction(6670053, 32445565)),
    (9, 29, ((2,) + ones(7),), Fraction(1319165, 92555776)),
    (13, 97, ((3,) + ones(10), (2, 2, 2, 2) + ones(5)),
     Fraction(54743, 8093902)),
    (11, 54, ((3,) + ones(8), (2, 2) + ones(7)),
     Fraction(198899, 26117585)),
    (14, 113, ((7,) + ones(7), (6, 2) + ones(6),
               (5, 3) + ones(6), (5, 2, 2) + ones(5),
               (4, 4) + ones(6), (4, 3, 2) + ones(5),
               (3, 3, 3) + ones(5), (3, 3, 2, 2, 2, 1, 1),
               (2, 2, 2, 2, 2, 2, 1, 1)), Fraction(2057, 73733672)),
    (13, 98, ((3,) + ones(10), (2, 2, 2) + ones(7)),
     Fraction(1313225, 43330928)),
)


def build_upsets(certificate):
    answer = []
    for degree, expected_size, generators, _ in certificate:
        bank = partitions(degree)
        upset = frozenset(
            shape for shape in bank
            if any(coarsens(generator, shape) for generator in generators)
        )
        require(len(upset) == expected_size,
                ("upset size drift", degree, generators))
        require((degree,) in upset and ones(degree) not in upset,
                ("upset not nonempty proper", degree, generators))
        require(all(not coarsens(fine, coarse) or coarse in upset
                    for fine in upset for coarse in bank),
                ("upset lost closure", degree, generators))
        answer.append(upset)
    return tuple(answer)


UPSETS = build_upsets(CERTIFICATE)
OLD_UPSETS = build_upsets(OLD_CERTIFICATE)
POLES, _ = reduced_poles(1, BANK, 1, 3)
RESET = residual_roots(1, dominant_row(1), 1, 3)
require(POLES == (8, 7, 6, 5, 5, 4, 4, 3, 3, 2, 2, 2, 1, 1, 1, 1),
        "physical pole bank drift")
require(RESET == (1, 3, 3, 4, 5, 6, 7, 8), "reset drift")
MULTIPLICITY = Counter(POLES)
VALUES = tuple(sorted(MULTIPLICITY))
LAYERS = tuple(
    tuple(
        state for state in combinations_with_replacement(VALUES, depth)
        if all(Counter(state)[value] <= MULTIPLICITY[value]
               for value in set(state))
    )
    for depth in range(1, len(POLES) + 1)
)
EXPECTED_DEPTHS = (8, 33, 93, 200, 348, 507, 631, 678, 631, 507,
                   348, 200, 93, 33, 8, 1)
require(tuple(map(len, LAYERS)) == EXPECTED_DEPTHS, "full-bank census drift")
STATES = tuple(state for layer in LAYERS for state in layer)
require(len(STATES) == 4319 == reduce(
    lambda left, right: left * right,
    (count + 1 for count in MULTIPLICITY.values()), 1) - 1,
    "physical multiset exhaustion drift")


def compute_state(state):
    return state, coefficient_vectors(1, BANK, 1, 3, state)


def clear_multipliers(certificate):
    rational = tuple(weight / NORMALIZERS[degree]
                     for degree, _, _, weight in certificate)
    common = reduce(lcm, (weight.denominator for weight in rational), 1)
    primitive = tuple(weight.numerator * (common // weight.denominator)
                      for weight in rational)
    divisor = reduce(gcd, primitive)
    return tuple(value // divisor for value in primitive)


def main():
    with Pool(processes=4) as pool:
        vectors = dict(pool.imap_unordered(compute_state, STATES, chunksize=2))
    require(len(vectors) == len(STATES), "response reconstruction lost states")
    observed_normalizers = {
        degree: max(abs(vectors[state][degree][shape].numerator)
                    for state in STATES for shape in partitions(degree))
        for degree in range(5, 15)
    }
    require(observed_normalizers == NORMALIZERS, "normalizer drift")
    require(all(vectors[state][degree][shape].denominator == 1
                for state in STATES for degree in range(5, 15)
                for shape in partitions(degree)),
            "response integrality drift")

    def rows_for(certificate, upsets):
        return tuple(
            tuple(sum(vectors[state][degree][shape] for shape in upset).numerator
                  for state in STATES)
            for (degree, _, _, _), upset in zip(certificate, upsets)
        )

    rows = rows_for(CERTIFICATE, UPSETS)
    multipliers = clear_multipliers(CERTIFICATE)
    require(all(value > 0 for value in multipliers)
            and reduce(gcd, multipliers) == 1,
            "full-bank dual lost primitive positivity")
    coordinates = tuple(
        sum(multiplier * row[index]
            for multiplier, row in zip(multipliers, rows))
        for index in range(len(STATES))
    )
    zero_states = tuple(state for state, value in zip(STATES, coordinates)
                        if value == 0)
    require(zero_states == (RESET,), "full-bank dual lost unique reset")
    require(all(value < 0 for state, value in zip(STATES, coordinates)
                if state != RESET), "full-bank dual lost strict negativity")
    sign_census = tuple(
        (sum(value < 0 for state, value in zip(STATES, coordinates)
             if len(state) == depth),
         sum(value == 0 for state, value in zip(STATES, coordinates)
             if len(state) == depth),
         sum(value > 0 for state, value in zip(STATES, coordinates)
             if len(state) == depth))
        for depth in range(1, 17)
    )
    require(sign_census == tuple(
        (count - int(depth == 8), int(depth == 8), 0)
        for depth, count in enumerate(EXPECTED_DEPTHS, 1)),
        "full-bank sign census drift")
    closest = max(((state, value) for state, value in zip(STATES, coordinates)
                   if value), key=lambda item: item[1])
    farthest = min(((state, value) for state, value in zip(STATES, coordinates)
                    if value), key=lambda item: item[1])
    require(closest[0] == (1,) and farthest[0] == (8,),
            "extreme-state identity drift")

    old_rows = rows_for(OLD_CERTIFICATE, OLD_UPSETS)
    old_multipliers = clear_multipliers(OLD_CERTIFICATE)
    old_coordinates = tuple(
        sum(multiplier * row[index]
            for multiplier, row in zip(old_multipliers, old_rows))
        for index in range(len(STATES))
    )
    old_tail_census = tuple(
        (depth,
         sum(value < 0 for state, value in zip(STATES, old_coordinates)
             if len(state) == depth),
         sum(value == 0 for state, value in zip(STATES, old_coordinates)
             if len(state) == depth),
         sum(value > 0 for state, value in zip(STATES, old_coordinates)
             if len(state) == depth))
        for depth in range(10, 17)
    )
    require(old_tail_census == (
        (10, 488, 0, 19), (11, 313, 0, 35), (12, 155, 0, 45),
        (13, 52, 0, 41), (14, 9, 0, 24), (15, 0, 0, 8),
        (16, 0, 0, 1)), "THM-3216 resurrection census drift")
    full_state = tuple(sorted(POLES))
    delete_eight = tuple(value for value in full_state if value != 8)
    full_index, delete_index = STATES.index(full_state), STATES.index(delete_eight)
    r5 = old_rows[0]
    require(all(r5[index] == 1440 * (
        sum(value ** 5 for value in RESET)
        - sum(value ** 5 for value in state))
                for index, state in enumerate(STATES)),
            "global degree-five power identity drift")
    require(old_coordinates[full_index] > 0 and r5[full_index] == -6117120,
            "full-state one-row hostile drift")
    require(old_coordinates[delete_index] > 0
            and r5[delete_index] == 41068800,
            "delete-eight one-row hostile drift")

    # Exhaust every principal upset.  Only twenty-eight point negatively on
    # both endpoint hostiles, and none admits a nonnegative one-row repair of
    # the old cocircuit on all 4,319 states.
    endpoint_candidates = 0
    principal_repairs = 0
    principal_total = 0
    for degree in range(5, 15):
        bank = partitions(degree)
        for generator in bank:
            upset = frozenset(shape for shape in bank
                              if coarsens(generator, shape))
            if len(upset) == len(bank):
                continue
            principal_total += 1
            row = tuple(sum(vectors[state][degree][shape] for shape in upset).numerator
                        for state in STATES)
            if row[full_index] >= 0 or row[delete_index] >= 0:
                continue
            endpoint_candidates += 1
            lower = Fraction(0)
            upper = None
            failed = False
            for index, state in enumerate(STATES):
                if state == RESET:
                    continue
                h, value = old_coordinates[index], row[index]
                if value < 0:
                    lower = max(lower, Fraction(h, -value))
                elif value > 0:
                    bound = Fraction(-h, value)
                    upper = bound if upper is None else min(upper, bound)
                elif h >= 0:
                    failed = True
                    break
            if not failed and (upper is None or lower < upper) \
                    and (upper is None or upper > 0):
                principal_repairs += 1
    require((endpoint_candidates, principal_repairs) == (28, 0),
            "single-principal-row no-go drift")

    inherited = sum(any(degree == old_degree and upset == old_upset
                        for old_degree, old_upset in zip(
                            (item[0] for item in OLD_CERTIFICATE), OLD_UPSETS))
                    for (degree, _, _, _), upset in zip(CERTIFICATE, UPSETS))
    principal = sum(len(generators) == 1
                    for _, _, generators, _ in CERTIFICATE[inherited:])
    require((inherited, principal, len(CERTIFICATE) - inherited - principal)
            == (11, 3, 8), "row provenance census drift")

    state_digest = hashlib.sha256(repr(STATES).encode("ascii")).hexdigest()
    row_digest = hashlib.sha256(repr(CERTIFICATE).encode("ascii")).hexdigest()
    multiplier_text = ",".join(map(str, multipliers)).encode("ascii")
    print("THM-3238 complete physical-bank unique reset exact control")
    print("dependency_hash_checks=" + repr(len(DEPENDENCIES)))
    print("poles_reset=" + repr((POLES, RESET)))
    print("state_counts_depth1_through16_total="
          + repr((EXPECTED_DEPTHS, len(STATES))))
    print("state_bank_sha256=" + state_digest)
    print("certificate_rows=" + repr(CERTIFICATE))
    print("certificate_rows_sha256=" + row_digest)
    print("row_provenance_inherited_principal_nonprincipal="
          + repr((inherited, principal, len(CERTIFICATE) - inherited - principal)))
    print("degree_normalizers=" + repr(NORMALIZERS))
    print("primitive_multiplier_count_gcd_minbits_maxbits=" + repr((
        len(multipliers), reduce(gcd, multipliers),
        min(value.bit_length() for value in multipliers),
        max(value.bit_length() for value in multipliers))))
    print("primitive_multipliers=" + repr(multipliers))
    print("primitive_multipliers_sha256="
          + hashlib.sha256(multiplier_text).hexdigest())
    print("sign_census_by_depth_negative_zero_positive=" + repr(sign_census))
    print("unique_zero_state=" + repr(zero_states))
    print("closest_negative_state_coordinate=" + repr(closest))
    print("farthest_negative_state_coordinate=" + repr(farthest))
    print("old_tail_resurrection_census=" + repr(old_tail_census))
    print("r5_endpoint_hostiles=" + repr((
        (full_state, old_coordinates[full_index], r5[full_index]),
        (delete_eight, old_coordinates[delete_index], r5[delete_index]))))
    print("principal_upsets_total_both_endpoint_negative_repairs="
          + repr((principal_total, endpoint_candidates, principal_repairs)))
    print("scope=support_(1,3)_bank_I2_complete_physical_multiset_bank_only")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
