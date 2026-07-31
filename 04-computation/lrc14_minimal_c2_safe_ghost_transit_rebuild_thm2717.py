#!/usr/bin/env python3
"""Self-contained exact minimal one-gate ghost repair and safe-complement audit.

The minimal deleted-gate word is

  D_H * product_i g_(q_i) * g_(c3) * g_(c1,ell) * I_(2h+kappa),

with the c2 factor omitted.  The safe-repair word inserts g_(c2).  The old
guard-danger word inserts D_(c2), and deletion = safe + old_danger is checked
coefficientwise on every locally used rail.  This is a changed transit
grammar, not an old terminal word/current.
"""

from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
import hashlib
import importlib.util
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = (ROOT / "04-computation" /
        "lrc14_central_half_odometer_full_local_cycle_thm2698.py")
BASE_SHA256 = "45cc393a856c00342fdf84875a0bc5a6d4c3df196ab35bb9ac2aad3cfc966c25"
SPEC = importlib.util.spec_from_file_location("half", BASE)
half = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(half)

private = half.private
core = half.core
old = private.old
P, Q7, T = private.P, private.Q7, private.T
R, S = private.R, private.S


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path):
    data = path.read_bytes().replace(b"\r\n", b"\n")
    return hashlib.sha256(data).hexdigest()


def freeze(value):
    return tuple(freeze(item) for item in value) if isinstance(value, list) else value


def transit_word(module, mode):
    deletion = module.make_comb(
        module.W[module.GUARD], 91, -13, 13
    )
    for index in module.UNIT_IDX:
        deletion = module.subtract_comb(
            deletion, module.W[index], 182, -13, 13
        )
    deletion = module.subtract_comb(
        deletion, module.C3, 182, -13, 13
    )
    if mode == "deletion":
        return deletion
    if mode == "safe":
        return module.subtract_comb(
            deletion, module.C2, 182, -13, 13
        )
    if mode == "old_danger":
        return module.intersect_comb(
            deletion, module.C2, 182, -13, 13
        )
    raise ValueError(f"unknown transit grammar mode {mode}")


def transit_prefixes(module, mode):
    """Return [clock][h][kappa] prefixes for one changed grammar."""
    word = transit_word(module, mode)
    result = []
    masses = []
    for ell in range(Q7):
        qell = module.subtract_comb(
            word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        by_h = []
        by_h_mass = []
        for h in range(P):
            pair = []
            pair_mass = []
            for kappa in range(2):
                digit = 2 * h + kappa
                interval = old.sat.intersect_interval(
                    qell, digit * T // (2 * P),
                    (digit + 1) * T // (2 * P)
                )
                prefix = module.make_prefix(interval)
                pair.append(prefix)
                pair_mass.append(prefix[2][-1])
            by_h.append(tuple(pair))
            by_h_mass.append(tuple(pair_mass))
        result.append(tuple(by_h))
        masses.append(tuple(by_h_mass))
    require(all(mass % P == 0
                for by_h in masses for pair in by_h for mass in pair),
            "changed delayed half-digit mass does not descend")
    return tuple(result), tuple(masses)


def shard_mode(task):
    bounds, mode = task
    start, stop = bounds
    module, _, _, _, rails, present, starts = core.build_carrier_data()
    prefixes, _ = transit_prefixes(module, mode)
    caches = [[{} for _ in range(P)] for _ in range(Q7)]
    content = positive = partition_checks = singleton_checks = 0
    metadata = []
    rows_out = []
    for j in range(start, stop):
        source, owner, theta, pieces = rails[j]
        metadata.append((source, owner, theta))
        raw = [[[[[0] * Q7 for _ in range(P)] for _ in range(2)]
                for _ in range(P)] for _ in range(2)]
        for h in range(P):
            for ell5 in range(Q7):
                overlap = old.intersect_weighted_union(
                    pieces, present[ell5, (-h) % P],
                    starts[ell5, (-h) % P]
                )
                for root in range(1, P):
                    halves = (
                        old.intersect_weighted_comb(
                            overlap, module.C3, 182,
                            14 * root - 13, 14 * root
                        ),
                        old.intersect_weighted_comb(
                            overlap, module.C3, 182,
                            14 * root, 14 * root + 13
                        ),
                    )
                    for edge, half_tooth in enumerate(halves):
                        values = private.delayed_carry_pair(
                            half_tooth, prefixes[ell5][h],
                            caches[ell5][h]
                        )
                        partition_checks += 2
                        for carry in range(P):
                            for kappa in range(2):
                                digit = 2 * h + kappa
                                predicted = (
                                    2 * carry + digit // P + (edge == 0)
                                ) % P
                                value = values[carry][kappa]
                                if root == predicted:
                                    raw[edge][carry][kappa][h][ell5] = value
                                    if value:
                                        positive += 1
                                        content = gcd(content, value)
                                else:
                                    require(value == 0,
                                            "changed row is not private")
                                singleton_checks += 1
        rows_out.append(freeze(raw))
    return (bounds, content, positive, partition_checks, singleton_checks,
            tuple(metadata), tuple(rows_out))


def primitive_data(values, root, content):
    require(root != 0 and content > 0,
            "primitive test needs nonzero root/content")
    require(all(value % content == 0 for value in values),
            "profile escaped common content")
    scale = pow(root, -1, P)
    normalized = tuple((value // content) * scale % P for value in values)
    reduced = tuple((normalized[index] - normalized[-1]) % P
                    for index in range(Q7 - 1))
    determinant = old.sat.multiplication_determinant_7(reduced)
    return normalized, reduced, determinant


def v_p(value, prime):
    valuation = 0
    while value and value % prime == 0:
        value //= prime
        valuation += 1
    return valuation, value


def combine(results):
    content = positives = partition_checks = singleton_checks = 0
    metadata = []
    rows = []
    shard_contents = []
    indices = []
    for bounds, shard_content, positive, pc, sc, meta, raw in results:
        shard_contents.append(shard_content)
        content = gcd(content, shard_content)
        positives += positive
        partition_checks += pc
        singleton_checks += sc
        metadata.extend(meta)
        rows.extend(raw)
        indices.extend(range(*bounds))
    return {
        "content": content,
        "positives": positives,
        "partition_checks": partition_checks,
        "singleton_checks": singleton_checks,
        "metadata": tuple(metadata),
        "rows": tuple(rows),
        "indices": tuple(indices),
        "shard_contents": tuple(shard_contents),
    }


def flatten(value):
    if isinstance(value, tuple):
        for item in value:
            yield from flatten(item)
    else:
        yield value


def slice_stats(row, edge_h_kappa, global_content=None):
    h, kappa = edge_h_kappa
    local_content = 0
    profiles = []
    for edge in range(2):
        for carry in range(P):
            values = row[edge][carry][kappa][h]
            if not any(values):
                continue
            root = (2 * carry + (2 * h + kappa) // P
                    + (edge == 0)) % P
            require(root != 0, "positive local profile has zero root")
            for value in values:
                local_content = gcd(local_content, value)
            profiles.append((edge, carry, root, values))
    local_units = []
    global_units = []
    determinant_hist = Counter()
    for edge, carry, root, values in profiles:
        normalized, reduced, determinant = primitive_data(
            values, root, local_content
        )
        determinant_hist[determinant] += 1
        if determinant:
            local_units.append((edge, carry, root, normalized, determinant))
        if global_content:
            _, _, global_det = primitive_data(
                values, root, global_content
            )
            if global_det:
                global_units.append((edge, carry, root, global_det))
    v13, residual = v_p(local_content, P)
    return {
        "content": local_content,
        "v13": v13,
        "residual13": residual,
        "profiles": tuple(profiles),
        "profile_count": len(profiles),
        "local_units": tuple(local_units),
        "local_unit_count": len(local_units),
        "global_unit_count": len(global_units),
        "determinant_hist": tuple(sorted(determinant_hist.items())),
    }


def row_map(bank):
    return dict(zip(bank["indices"], bank["rows"]))


def bank_profile_stats(rows, content):
    profiles = units = 0
    determinant_hist = Counter()
    nonunit_examples = []
    for j, row in enumerate(rows):
        for edge in range(2):
            for carry in range(P):
                for kappa in range(2):
                    for h in range(P):
                        values = row[edge][carry][kappa][h]
                        if not any(values):
                            continue
                        profiles += 1
                        root = (2 * carry + (2 * h + kappa) // P
                                + (edge == 0)) % P
                        require(root != 0,
                                "positive global profile has zero root")
                        normalized, reduced, determinant = primitive_data(
                            values, root, content
                        )
                        determinant_hist[determinant] += 1
                        if determinant:
                            units += 1
                        elif len(nonunit_examples) < 8:
                            nonunit_examples.append(
                                (j, edge, carry, kappa, h, root,
                                 normalized, reduced)
                            )
    return {
        "profiles": profiles,
        "units": units,
        "determinant_hist": tuple(sorted(determinant_hist.items())),
        "nonunit_examples": tuple(nonunit_examples),
    }


def paired_bank_unit_table(deletion_rows, safe_rows, content):
    """Cross-tab deletion, safe, and derived old-danger units at G=26."""
    table = Counter()
    examples = {}
    for j, (deletion_row, safe_row) in enumerate(
            zip(deletion_rows, safe_rows)):
        for edge in range(2):
            for carry in range(P):
                for kappa in range(2):
                    for h in range(P):
                        deletion_values = deletion_row[edge][carry][kappa][h]
                        safe_values = safe_row[edge][carry][kappa][h]
                        require(all(0 <= safe <= deleted
                                    for deleted, safe
                                    in zip(deletion_values, safe_values)),
                                "safe row escaped deleted-gate row")
                        if not any(deletion_values) and not any(safe_values):
                            continue
                        require(any(deletion_values) and any(safe_values),
                                "global deletion/safe profile supports differ")
                        old_values = tuple(
                            deleted - safe for deleted, safe
                            in zip(deletion_values, safe_values)
                        )
                        require(any(old_values),
                                "global old-danger profile unexpectedly empty")
                        root = (2 * carry + (2 * h + kappa) // P
                                + (edge == 0)) % P
                        deletion_det = primitive_data(
                            deletion_values, root, content
                        )[-1]
                        safe_det = primitive_data(
                            safe_values, root, content
                        )[-1]
                        old_det = primitive_data(
                            old_values, root, content
                        )[-1]
                        key = (bool(deletion_det), bool(safe_det),
                               bool(old_det))
                        table[key] += 1
                        examples.setdefault(
                            key, (j, edge, carry, kappa, h, root,
                                  deletion_det, safe_det, old_det)
                        )
    return tuple(sorted(table.items())), tuple(sorted(examples.items()))


def multiply_phi7(left, right):
    """Multiply six-term representatives in F13[z]/(Phi7)."""
    require(len(left) == len(right) == 6,
            "Phi7 multiplication needs six-term representatives")
    coefficients = [0] * 11
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            coefficients[i + j] = (coefficients[i + j] + x * y) % P
    for degree in range(10, 5, -1):
        leading = coefficients[degree]
        coefficients[degree] = 0
        for shift in range(1, 7):
            coefficients[degree - shift] = (
                coefficients[degree - shift] - leading
            ) % P
    return tuple(coefficients[:6])


def truth_sheet_slice_stats(deletion_row, safe_row, old_row, h, kappa,
                            content):
    """Karoubi truth-sheet units and swap norms on one 22-profile slice."""
    profiles = []
    norm_hist = Counter()
    cancellation = []
    pair_units = norm_units = 0
    for edge in range(2):
        for carry in range(P):
            deletion_values = deletion_row[edge][carry][kappa][h]
            safe_values = safe_row[edge][carry][kappa][h]
            old_values = old_row[edge][carry][kappa][h]
            if not any(deletion_values):
                continue
            root = (2 * carry + (2 * h + kappa) // P
                    + (edge == 0)) % P
            deletion = primitive_data(deletion_values, root, content)
            safe = primitive_data(safe_values, root, content)
            old_part = primitive_data(old_values, root, content)
            deletion_reduced = deletion[1]
            safe_reduced = safe[1]
            old_reduced = old_part[1]
            require(deletion_reduced == tuple(
                        (x + y) % P
                        for x, y in zip(safe_reduced, old_reduced)),
                    "truth-sheet codiagonal failed in Phi7 quotient")
            truth_unit = bool(safe[-1] and old_part[-1])
            pair_units += truth_unit
            norm = multiply_phi7(safe_reduced, old_reduced)
            require(norm == multiply_phi7(old_reduced, safe_reduced),
                    "truth-sheet swap norm is not invariant")
            norm_det = old.sat.multiplication_determinant_7(norm)
            require(bool(norm_det) == truth_unit,
                    "truth-sheet norm/unit criterion failed")
            norm_units += bool(norm_det)
            norm_hist[norm_det] += 1
            difference = tuple((x - y) % P
                               for x, y in zip(safe_reduced, old_reduced))
            discriminant = tuple(
                (x - 4 * y) % P
                for x, y in zip(
                    multiply_phi7(deletion_reduced, deletion_reduced), norm
                )
            )
            require(discriminant == multiply_phi7(difference, difference),
                    "truth-sheet discriminant identity failed")
            inverse_two = pow(2, -1, P)
            require(tuple(inverse_two * (x + y) % P
                          for x, y in zip(deletion_reduced, difference))
                    == safe_reduced
                    and tuple(inverse_two * (x - y) % P
                              for x, y
                              in zip(deletion_reduced, difference))
                    == old_reduced,
                    "truth-sheet quadratic reconstruction failed")
            discriminant_det = old.sat.multiplication_determinant_7(
                discriminant
            )
            if not any(deletion_reduced):
                require(old_reduced == tuple((-x) % P for x in safe_reduced),
                        "zero codiagonal is not an anti-diagonal pair")
                require(discriminant_det != 0,
                        "reverse anti-diagonal discriminant lost etaleness")
                cancellation.append((edge, carry, root))
            profiles.append((edge, carry, root, deletion_reduced,
                             safe_reduced, old_reduced, norm, norm_det,
                             discriminant, discriminant_det))
    return {
        "profile_count": len(profiles),
        "pair_units": pair_units,
        "norm_units": norm_units,
        "norm_hist": tuple(sorted(norm_hist.items())),
        "cancellation": tuple(cancellation),
        "profiles": tuple(profiles),
    }


def compact_slice(stats):
    return (
        stats["content"], stats["v13"], stats["residual13"],
        stats["profile_count"], stats["local_unit_count"],
        stats["determinant_hist"],
    )


def compact_forced(profile):
    root, _, normalized, _, determinant = profile
    return root, normalized, determinant


def reduce_phi7(coefficients):
    require(len(coefficients) == Q7, "Phi7 reduction needs seven clocks")
    return tuple((coefficients[index] - coefficients[-1]) % P
                 for index in range(Q7 - 1))


def canonical_profile(coefficients):
    return reduce_phi7(coefficients) + (0,)


def rotate_profile(coefficients, shift):
    """Multiply a C7 group-ring profile by z^shift."""
    return canonical_profile(tuple(
        coefficients[(index - shift) % Q7] for index in range(Q7)
    ))


def invert_profile(coefficients):
    """Apply z -> z^-1, which is also absolute 13-Frobenius."""
    return canonical_profile(tuple(
        coefficients[(-index) % Q7] for index in range(Q7)
    ))


def profile_orbit(coefficients):
    """Orbit under scalars, cyclic shifts, and optional inversion."""
    coefficients = canonical_profile(coefficients)
    return {
        tuple(scalar * value % P for value in rotate_profile(base, shift))
        for inverse in range(2)
        for base in ((invert_profile(coefficients)
                      if inverse else coefficients),)
        for shift in range(Q7)
        for scalar in range(1, P)
    }


def orbit_relations(source, target):
    source = canonical_profile(source)
    target = canonical_profile(target)
    relations = []
    for inverse in range(2):
        base = invert_profile(source) if inverse else source
        for shift in range(Q7):
            shifted = rotate_profile(base, shift)
            for scalar in range(1, P):
                if tuple(scalar * value % P for value in shifted) == target:
                    relations.append((inverse, shift, scalar))
    return tuple(relations)


def quadratic_remainder(poly, middle):
    """Remainder modulo z^2+middle*z+1 over F13."""
    constant = linear = 0
    for coefficient in reversed(poly):
        constant, linear = (
            (coefficient - linear) % P,
            (constant - middle * linear) % P,
        )
    return constant, linear


def inverse_pair_components(coefficients):
    """Three F169 components for inverse character pairs (1,6),(2,5),(3,4)."""
    reduced = reduce_phi7(coefficients)
    components = tuple(
        (pair, middle, quadratic_remainder(reduced, middle))
        for pair, middle in (((1, 6), 3), ((2, 5), 6), ((3, 4), 5))
    )
    require(all(component != (0, 0) for _, _, component in components),
            "forced primitive unit has a zero inverse-pair component")
    return components


def component_frobenius(components):
    """Conjugate a+b*z by z -> z^-1=-middle-z in each F169 factor."""
    return tuple(
        (pair, middle, ((constant - middle * linear) % P, (-linear) % P))
        for pair, middle, (constant, linear) in components
    )


def component_multiply(left, right, middle):
    """Multiply in F13[z]/(z^2+middle*z+1)."""
    a, b = left
    c, d = right
    return ((a * c - b * d) % P,
            (a * d + b * c - middle * b * d) % P)


def component_power(value, exponent, middle):
    result = (1, 0)
    base = value
    while exponent:
        if exponent & 1:
            result = component_multiply(result, base, middle)
        base = component_multiply(base, base, middle)
        exponent //= 2
    return result


def component_order(value, middle):
    require(value != (0, 0), "multiplicative order needs a nonzero component")
    current = (1, 0)
    for order in range(1, P * P):
        current = component_multiply(current, value, middle)
        if current == (1, 0):
            require((P * P - 1) % order == 0,
                    "component order does not divide 168")
            return order
    raise RuntimeError("nonzero F169 component exceeded order 168")


def component_quotient(numerator, denominator, middle):
    inverse = component_power(denominator, P * P - 2, middle)
    require(component_multiply(denominator, inverse, middle) == (1, 0),
            "F169 inverse failed")
    return component_multiply(numerator, inverse, middle)


def component_order_report(forward, reverse):
    report = []
    quotients = []
    for ((pair_f, middle_f, value_f),
         (pair_r, middle_r, value_r)) in zip(forward, reverse):
        require((pair_f, middle_f) == (pair_r, middle_r),
                "inverse-pair component labels disagree")
        quotient = component_quotient(value_r, value_f, middle_f)
        quotients.append((pair_f, middle_f, quotient))
        report.append((pair_f, middle_f,
                       component_order(value_f, middle_f),
                       component_order(value_r, middle_f),
                       quotient, component_order(quotient, middle_f)))
    common_order4 = tuple(
        exponent for exponent in range(1, P * P)
        if all(component_order(
                   component_power(value, exponent, middle), middle) == 4
               for _, middle, value in quotients)
    )
    common_sqrt_minus_one = tuple(
        (exponent, scalar)
        for exponent in range(1, P * P)
        for scalar in (5, 8)
        if all(component_power(value, exponent, middle) == (scalar, 0)
               for _, middle, value in quotients)
    )
    common_scalar_powers = []
    for exponent in range(1, P * P):
        powers = tuple(component_power(value, exponent, middle)
                       for _, middle, value in quotients)
        if len(set(powers)) == 1 and powers[0][1] == 0:
            common_scalar_powers.append((exponent, powers[0][0]))
    quotient_order = lcm(*(entry[-1] for entry in report))
    require(quotient_order == 168
            and tuple(common_scalar_powers)
            == ((42, 5), (84, 12), (126, 8), (168, 1)),
            "order-168 lift / central scalar C4 changed")
    return (tuple(report), quotient_order, common_order4,
            common_sqrt_minus_one, tuple(common_scalar_powers))


def verify_sum(deletion, safe, old_danger):
    left = tuple(flatten(deletion))
    right_safe = tuple(flatten(safe))
    right_old = tuple(flatten(old_danger))
    require(len(left) == len(right_safe) == len(right_old),
            "repair row shapes differ")
    require(all(a == b + c for a, b, c in zip(left, right_safe, right_old)),
            "deleted c2 row does not split as safe plus old danger")


def point_profile(row, h, kappa, edge, carry, content):
    root = (2 * carry + (2 * h + kappa) // P + (edge == 0)) % P
    values = row[edge][carry][kappa][h]
    normalized, reduced, determinant = primitive_data(
        values, root, content
    )
    return root, values, normalized, reduced, determinant


def point_geometry(module, rails, present, starts, point, desired_clock):
    x, y = point["x"], point["y"]
    h, kappa = point["h"], point["kappa"]
    require(half.frac(R * x) == y
            and half.floor_fraction(R * x) == point["N"]
            and point["N"] % P == point["carry"],
            "point address/carry changed")
    require(half.shallow(x) == desired_clock
            and half.owner(x) == desired_clock,
            "point clock typing changed")
    coordinate = x * half.T
    require(half.strict_interval_member(
                coordinate, present[desired_clock, (-h) % P],
                starts[desired_clock, (-h) % P])
            and half.strict_weighted_member(
                coordinate, rails[point["j"]][3],
                [a for a, _, _ in rails[point["j"]][3]])
            and half.half_membership(
                module, x, point["edge"], point["root"]),
            "point rail/present/private-root geometry changed")


def split_point(source, middle, target, macro_lift):
    a = (middle["N"] - P * source["N"]
         - half.floor_fraction(P * source["y"])) % R
    b = (target["N"] - P * middle["N"]
         - half.floor_fraction(P * middle["y"])) % R
    require((P * a + b) % R == macro_lift
            and half.frac(P * source["x"] + F(a, R)) == middle["x"]
            and half.frac(P * middle["x"] + F(b, R)) == target["x"],
            "forced point affine split changed")
    source_roots = half.half_roots(
        MODULE, half.frac(P * source["x"]), middle["edge"]
    )
    middle_roots = half.half_roots(
        MODULE, half.frac(P * middle["x"]), target["edge"]
    )
    require(any((root + 2 * a) % P == middle["root"]
                for root in source_roots)
            and any((root + 2 * b) % P == target["root"]
                    for root in middle_roots),
            "forced point root split changed")
    return a, b, source_roots, middle_roots


MODULE = None


def main():
    global MODULE
    require(lf_sha256(BASE) == BASE_SHA256,
            "load-bearing THM2698 companion hash changed")
    full_tasks = tuple(
        (bounds, mode)
        for mode in ("deletion", "safe")
        for bounds in core.SHARDS
    )
    with ProcessPoolExecutor(max_workers=4) as pool:
        full_results = list(pool.map(shard_mode, full_tasks))
    deletion = combine(full_results[:4])
    safe = combine(full_results[4:])
    v13_deletion, residual_deletion = v_p(deletion["content"], P)
    v13_safe, residual_safe = v_p(safe["content"], P)
    deletion_profile_stats = bank_profile_stats(
        deletion["rows"], deletion["content"]
    )
    safe_profile_stats = bank_profile_stats(safe["rows"], safe["content"])
    global_unit_table, global_unit_examples = paired_bank_unit_table(
        deletion["rows"], safe["rows"], deletion["content"]
    )
    require(global_unit_table
            == (((False, False, False), 2618),
                ((False, True, True), 1565),
                ((True, False, True), 120),
                ((True, True, False), 155),
                ((True, True, True), 19414)),
            "global deletion/safe/old unit cross-tab changed")

    local_bounds = ((2, 4), (8, 10))
    tasks = tuple((bounds, "old_danger") for bounds in local_bounds)
    with ProcessPoolExecutor(max_workers=4) as pool:
        local_results = list(pool.map(shard_mode, tasks))
    old_danger = combine(local_results)

    drows = row_map(deletion)
    srows = row_map(safe)
    orows = row_map(old_danger)
    for j in (2, 3, 8, 9):
        verify_sum(drows[j], srows[j], orows[j])

    forward_stats = {
        mode: slice_stats(rows[9], (0, 1),
                          deletion["content"]
                          if mode in ("deletion", "safe") else None)
        for mode, rows in (("deletion", drows), ("safe", srows),
                           ("old_danger", orows))
    }
    reverse_stats = {
        mode: slice_stats(rows[2], (12, 0),
                          deletion["content"]
                          if mode in ("deletion", "safe") else None)
        for mode, rows in (("deletion", drows), ("safe", srows),
                           ("old_danger", orows))
    }
    require(forward_stats["deletion"]["profile_count"] == 22
            and forward_stats["deletion"]["local_unit_count"] == 22,
            "forward minimal deletion local census changed")
    require(reverse_stats["deletion"]["profile_count"] == 22
            and reverse_stats["deletion"]["local_unit_count"] == 2,
            "reverse minimal deletion local census changed")
    require(compact_slice(forward_stats["deletion"])
            == (2225652, 1, 171204, 22, 22,
                ((1, 11), (3, 1), (10, 1), (12, 9)))
            and compact_slice(forward_stats["safe"])
            == (2225652, 1, 171204, 22, 22,
                ((2, 10), (4, 1), (6, 1), (7, 1), (9, 1), (11, 8)))
            and compact_slice(forward_stats["old_danger"])
            == (472491719940240, 1, 36345516918480, 22, 22,
                ((1, 1), (2, 3), (5, 6), (8, 6), (11, 5), (12, 1))),
            "forward deletion/safe/old-danger exact census changed")
    require(compact_slice(reverse_stats["deletion"])
            == (1092, 1, 84, 22, 2, ((0, 20), (1, 1), (12, 1)))
            and compact_slice(reverse_stats["safe"])
            == (1092, 1, 84, 22, 22, ((1, 10), (12, 12)))
            and compact_slice(reverse_stats["old_danger"])
            == (16558902360, 1, 1273761720, 22, 22,
                ((1, 12), (12, 10))),
            "reverse deletion/safe/old-danger exact census changed")
    require(deletion["content"] == 26
            and deletion["shard_contents"] == (26, 26, 26, 26)
            and deletion["positives"] == 59424
            and deletion_profile_stats["profiles"] == 23872
            and deletion_profile_stats["units"] == 19689
            and deletion_profile_stats["determinant_hist"]
            == ((0, 4183), (1, 5429), (2, 337), (3, 1196),
                (4, 972), (5, 1106), (6, 889), (7, 887),
                (8, 1039), (9, 947), (10, 1294), (11, 400),
                (12, 5193)),
            "global deleted-gate control census changed")
    require(safe["content"] == 26
            and safe["shard_contents"] == (26, 26, 26, 26)
            and safe["positives"] == 59424
            and safe_profile_stats["profiles"] == 23872
            and safe_profile_stats["units"] == 21134
            and safe_profile_stats["determinant_hist"]
            == ((0, 2738), (1, 4553), (2, 605), (3, 1830),
                (4, 1276), (5, 1061), (6, 1137), (7, 1100),
                (8, 1058), (9, 1418), (10, 1769), (11, 593),
                (12, 4734)),
            "global c2-safe repair census changed")
    require(forward_stats["safe"]["profile_count"] == 22
            and forward_stats["safe"]["local_unit_count"] == 22
            and reverse_stats["safe"]["profile_count"] == 22
            and reverse_stats["safe"]["local_unit_count"] == 22,
            "forced local safe-repair census changed")
    reverse_deletion_support = tuple(
        (edge, carry, root)
        for edge, carry, root, _, _
        in reverse_stats["deletion"]["local_units"]
    )
    require(reverse_deletion_support == ((0, 0, 2), (1, 0, 1)),
            "reverse deleted-gate unit support changed")
    forward_truth = truth_sheet_slice_stats(
        drows[9], srows[9], orows[9], 0, 1, deletion["content"]
    )
    reverse_truth = truth_sheet_slice_stats(
        drows[2], srows[2], orows[2], 12, 0, deletion["content"]
    )
    require(forward_truth["profile_count"] == 22
            and forward_truth["pair_units"] == 22
            and forward_truth["norm_units"] == 22
            and forward_truth["norm_hist"]
            == ((3, 10), (4, 8), (9, 4)),
            "forward truth-sheet unit/norm census changed")
    require(reverse_truth["profile_count"] == 22
            and reverse_truth["pair_units"] == 22
            and reverse_truth["norm_units"] == 22
            and reverse_truth["norm_hist"] == ((1, 20), (12, 2))
            and len(reverse_truth["cancellation"]) == 20
            and all(carry != 0
                    for _, carry, _ in reverse_truth["cancellation"]),
            "reverse anti-diagonal truth-sheet census changed")
    require(all(profile[9] == 1
                for profile in reverse_truth["profiles"]
                if not any(profile[3])),
            "reverse anti-diagonal discriminant determinant changed")
    require(all(
                all(value == 0 for index, value in enumerate(profile[4])
                    if index != 2)
                and all(value == 0 for index, value in enumerate(profile[5])
                        if index != 2)
                for profile in reverse_truth["profiles"]
            ),
            "reverse truth-sheet arms escaped the z^2 line")

    MODULE, _, _, _, rails, present, starts = core.build_carrier_data()
    prefixes = {
        mode: transit_prefixes(MODULE, mode)[0]
        for mode in ("deletion", "safe", "old_danger")
    }
    for y, ell, h, kappa in ((F(1, 17), 4, 0, 1),
                              (F(16, 17), 1, 12, 0)):
        require(half.strict_interval_member(
                    y, half.prefix_intervals(prefixes["deletion"][ell][h][kappa]))
                and half.strict_interval_member(
                    y, half.prefix_intervals(prefixes["safe"][ell][h][kappa]))
                and not half.strict_interval_member(
                    y, half.prefix_intervals(prefixes["old_danger"][ell][h][kappa])),
                "ghost point deletion/safe/old-danger typing changed")

    source0 = {"x": F(39123022, 82055753), "N": 2301354,
               "y": F(4, 17), "carry": 3, "edge": 0, "root": 7}
    source1 = {"x": F(41305372, 82055753), "N": 2429727,
               "y": F(13, 17), "carry": 1, "edge": 0, "root": 4}
    forward_point = {"x": F(41513423, 82055753), "N": 2441966,
                     "y": F(1, 17), "j": 9, "carry": 7,
                     "h": 0, "kappa": 1, "edge": 0, "root": 2}
    reverse_point = {"x": F(38400313, 82055753), "N": 2258841,
                     "y": F(16, 17), "j": 2, "carry": 0,
                     "h": 12, "kappa": 0, "edge": 0, "root": 2}
    point_geometry(MODULE, rails, present, starts, forward_point, 4)
    point_geometry(MODULE, rails, present, starts, reverse_point, 1)
    forward_split = split_point(source0, forward_point, source1, 4472391)
    reverse_split = split_point(source1, reverse_point, source0, 1956127)

    forced_profiles = {}
    for name, point, j in (("forward", forward_point, 9),
                           ("reverse", reverse_point, 2)):
        forced_profiles[name] = {}
        for mode, rows, stats in (
                ("deletion", drows, forward_stats if name == "forward" else reverse_stats),
                ("safe", srows, forward_stats if name == "forward" else reverse_stats),
                ("old_danger", orows, forward_stats if name == "forward" else reverse_stats)):
            # Use the one common global lattice for the displayed profiles.
            # The directly integrated old-danger rows are divisible by 26
            # because deletion=safe+old_danger coefficientwise.
            content = deletion["content"]
            if content:
                forced_profiles[name][mode] = point_profile(
                    rows[j], point["h"], point["kappa"], point["edge"],
                    point["carry"], content
                )
            else:
                forced_profiles[name][mode] = None

    require(forced_profiles["forward"]["safe"][-1] != 0
            and forced_profiles["reverse"]["safe"][-1] != 0,
            "forced safe-repair ghost profile lost primitivity")
    forward_safe = forced_profiles["forward"]["safe"][2]
    reverse_safe = forced_profiles["reverse"]["safe"][2]
    forward_components = inverse_pair_components(forward_safe)
    reverse_components = inverse_pair_components(reverse_safe)
    forward_frobenius = invert_profile(forward_safe)
    reverse_frobenius = invert_profile(reverse_safe)
    require(inverse_pair_components(forward_frobenius)
            == component_frobenius(forward_components)
            and inverse_pair_components(reverse_frobenius)
            == component_frobenius(reverse_components),
            "profile inversion disagrees with F169 component Frobenius")
    forward_orbit = profile_orbit(forward_safe)
    reverse_orbit = profile_orbit(reverse_safe)
    relations = orbit_relations(forward_safe, reverse_safe)
    (component_orders, quotient_order, common_order4_exponents,
     common_sqrt_minus_one,
     common_scalar_powers) = component_order_report(
         forward_components, reverse_components
    )
    require(len(forward_orbit) in (84, 168)
            and len(reverse_orbit) in (84, 168),
            "unexpected scalar/shift/inversion stabilizer")
    require((reverse_safe in forward_orbit) == bool(relations),
            "orbit relation test disagrees with exhaustive orbit")

    print("LRC14 MINIMAL GHOST TRANSIT DELETION / SAFE-REPAIR EXACT SCOUT")
    print(f"base_thm2698_lf_sha256={BASE_SHA256}")
    print("global_deletion="
          f"(content={deletion['content']},v13={v13_deletion},"
          f"residual13={residual_deletion},positives={deletion['positives']},"
          f"profiles={deletion_profile_stats['profiles']},"
          f"units={deletion_profile_stats['units']},"
          f"shards={deletion['shard_contents']},"
          f"det_hist={deletion_profile_stats['determinant_hist']})")
    print("global_safe_repair="
          f"(content={safe['content']},v13={v13_safe},"
          f"residual13={residual_safe},positives={safe['positives']},"
          f"profiles={safe_profile_stats['profiles']},"
          f"units={safe_profile_stats['units']},"
          f"shards={safe['shard_contents']},"
          f"det_hist={safe_profile_stats['determinant_hist']})")
    print("exact_checks="
          f"(partition={safe['partition_checks']},"
          f"singleton={safe['singleton_checks']})")
    print("global_deletion_safe_unit_crosstab="
          f"(counts={global_unit_table},examples={global_unit_examples})")
    print("forward_slices="
          f"(deletion={compact_slice(forward_stats['deletion'])},"
          f"safe={compact_slice(forward_stats['safe'])},"
          f"old_danger={compact_slice(forward_stats['old_danger'])})")
    print("reverse_slices="
          f"(deletion={compact_slice(reverse_stats['deletion'])},"
          f"safe={compact_slice(reverse_stats['safe'])},"
          f"old_danger={compact_slice(reverse_stats['old_danger'])},"
          f"deletion_unit_support={reverse_deletion_support})")
    print("forced_forward="
          f"(deletion={compact_forced(forced_profiles['forward']['deletion'])},"
          f"safe={compact_forced(forced_profiles['forward']['safe'])},"
          f"old_danger={compact_forced(forced_profiles['forward']['old_danger'])},"
          f"split={forward_split})")
    print("forced_reverse="
          f"(deletion={compact_forced(forced_profiles['reverse']['deletion'])},"
          f"safe={compact_forced(forced_profiles['reverse']['safe'])},"
          f"old_danger={compact_forced(forced_profiles['reverse']['old_danger'])},"
          f"split={reverse_split})")
    print("forced_safe_polynomials="
          f"(forward7={forward_safe},forward6={reduce_phi7(forward_safe)},"
          f"reverse7={reverse_safe},reverse6={reduce_phi7(reverse_safe)})")
    print("forced_safe_symmetry="
          f"(forward_frobenius={forward_frobenius},"
          f"reverse_frobenius={reverse_frobenius},"
          f"forward_inverse_pairs={forward_components},"
          f"reverse_inverse_pairs={reverse_components},"
          f"forward_orbit_size={len(forward_orbit)},"
          f"reverse_orbit_size={len(reverse_orbit)},"
          f"forward_to_reverse_relations={relations})")
    print("forced_safe_component_orders="
          f"(pairs={component_orders},"
          f"quotient_order={quotient_order},"
          f"common_order4_exponents={common_order4_exponents},"
          f"common_sqrt_minus_one={common_sqrt_minus_one},"
          f"common_scalar_powers={common_scalar_powers})")
    print("forward_truth_sheet="
          f"(pair_units={forward_truth['pair_units']},"
          f"norm_units={forward_truth['norm_units']},"
          f"norm_hist={forward_truth['norm_hist']},"
          f"cancellation={forward_truth['cancellation']},"
          f"norm_discriminant_profiles={tuple((p[0], p[1], p[2], p[6], p[7], p[8], p[9]) for p in forward_truth['profiles'])})")
    print("reverse_truth_sheet="
          f"(pair_units={reverse_truth['pair_units']},"
          f"norm_units={reverse_truth['norm_units']},"
          f"norm_hist={reverse_truth['norm_hist']},"
          f"cancellation={reverse_truth['cancellation']},"
          f"norm_discriminant_profiles={tuple((p[0], p[1], p[2], p[6], p[7], p[8], p[9]) for p in reverse_truth['profiles'])})")
    print("row_identity=deletion=safe_repair+old_c2_danger on rails 2,3,8,9: PASS")
    print("scope=global-content-26 c2-safe changed transit universe with two forced units; deletion-only unit inference is invalid; no terminal word/current, live-row transfer, ledger decrement, or LRC conclusion")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
