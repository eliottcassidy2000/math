#!/usr/bin/env python3
"""Second exact audit of the all-q atomized coset-cochain criterion.

This file is intentionally independent of the candidate verifier.  It compares
the integral cochain equations directly with an exact circular event sweep,
then reconstructs every positive cochain through a generalized CRT and checks
the resulting real-line Helly intersection.  It also checks full owner covers,
strict capacity endpoints, core-speed typing, and the q=6 CRT cover clutter.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from math import gcd
from pathlib import Path


F = Fraction
EXPECTED_SEMANTIC = "bb394b88aa71ebe9a892c370a071227a3e016adb8d550163efc73dba70995ddc"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lcm(a, b):
    return a // gcd(a, b) * b


def mod_one(x):
    return x % 1


def integer_distance(x):
    r = mod_one(x)
    return min(r, 1 - r)


def image_order(q, u):
    return q // gcd(q, u)


def atom_mask(q, atom):
    u, k = atom
    m = image_order(q, u)
    return sum(1 << j for j in range(q) if j % m == k)


def atom_fires(q, atom, t):
    u, k = atom
    return 14 * integer_distance(u * (t + F(k, q))) < 1


def event_endpoints(q, atoms):
    endpoints = set()
    for u, k in atoms:
        for tooth in range(u):
            centre = F(tooth, u) - F(k, q)
            radius = F(1, 14 * u)
            endpoints.add(mod_one(centre - radius))
            endpoints.add(mod_one(centre + radius))
    ordered = tuple(sorted(endpoints))
    require(ordered, ("no endpoints", q, atoms))
    return ordered


def event_midpoints(q, atoms):
    ordered = event_endpoints(q, atoms)
    mids = []
    for i, left in enumerate(ordered):
        right = ordered[(i + 1) % len(ordered)]
        if i + 1 == len(ordered):
            right += 1
        mids.append(mod_one((left + right) / 2))
    return tuple(mids)


def direct_common_phase(q, atoms):
    for t in event_midpoints(q, atoms):
        if all(atom_fires(q, atom, t) for atom in atoms):
            return t
    return None


def allowed_gaps(q, left, right):
    u, k = left
    v, ell = right
    modulus = q * gcd(u, v)
    residue = ((ell - k) * u * v) % modulus
    bound = (q * (u + v) - 1) // 14
    return tuple(
        p for p in range(-bound, bound + 1)
        if (p - residue) % modulus == 0
    )


def oriented(edge_values, i, j):
    if i < j:
        return edge_values[(i, j)]
    return -edge_values[(j, i)]


def cochain_witness(q, atoms):
    edges = tuple(combinations(range(len(atoms)), 2))
    choices = tuple(allowed_gaps(q, atoms[i], atoms[j]) for i, j in edges)
    if any(not values for values in choices):
        return None
    for values in product(*choices):
        candidate = dict(zip(edges, values))
        closes = True
        for i, j, k in combinations(range(len(atoms)), 3):
            ui = atoms[i][0]
            uj = atoms[j][0]
            uk = atoms[k][0]
            circulation = (
                uk * oriented(candidate, i, j)
                + ui * oriented(candidate, j, k)
                + uj * oriented(candidate, k, i)
            )
            if circulation != 0:
                closes = False
                break
        if closes:
            return candidate
    return None


def all_cochains(q, atoms):
    edges = tuple(combinations(range(len(atoms)), 2))
    choices = tuple(allowed_gaps(q, atoms[i], atoms[j]) for i, j in edges)
    if any(not values for values in choices):
        return ()
    witnesses = []
    for values in product(*choices):
        candidate = dict(zip(edges, values))
        if all(
            atoms[k][0] * oriented(candidate, i, j)
            + atoms[i][0] * oriented(candidate, j, k)
            + atoms[j][0] * oriented(candidate, k, i) == 0
            for i, j, k in combinations(range(len(atoms)), 3)
        ):
            witnesses.append(candidate)
    return tuple(witnesses)


def star_compressed_witness(q, atoms):
    """Use only base-star variables; derive every non-star edge exactly."""
    if len(atoms) <= 1:
        return {}
    u0 = atoms[0][0]
    star_choices = tuple(allowed_gaps(q, atoms[0], atoms[i]) for i in range(1, len(atoms)))
    if any(not values for values in star_choices):
        return None
    for values in product(*star_choices):
        candidate = {(0, i): values[i - 1] for i in range(1, len(atoms))}
        valid = True
        for i, j in combinations(range(1, len(atoms)), 2):
            numerator = atoms[i][0] * candidate[(0, j)] - atoms[j][0] * candidate[(0, i)]
            if numerator % u0:
                valid = False
                break
            derived = numerator // u0
            if derived not in allowed_gaps(q, atoms[i], atoms[j]):
                valid = False
                break
            candidate[(i, j)] = derived
        if valid:
            return candidate
    return None


def merge_congruences(a, m, b, n):
    """Return the least solution and modulus, using generalized CRT."""
    g = gcd(m, n)
    require((b - a) % g == 0, ("incompatible CRT", a, m, b, n))
    m1 = m // g
    n1 = n // g
    if n1 == 1:
        step = 0
    else:
        step = ((b - a) // g * pow(m1, -1, n1)) % n1
    modulus = m * n1
    return (a + m * step) % modulus, modulus


def reconstruct_by_crt(q, atoms, witness):
    """Independently turn a positive cochain into coherent tooth centres."""
    if len(atoms) == 1:
        potentials = (F(0),)
    else:
        u0 = atoms[0][0]
        potentials = (F(0),) + tuple(
            -F(oriented(witness, 0, i), q * u0 * atoms[i][0])
            for i in range(1, len(atoms))
        )

    for i, j in combinations(range(len(atoms)), 2):
        expected = F(oriented(witness, i, j), q * atoms[i][0] * atoms[j][0])
        require(potentials[i] - potentials[j] == expected,
                ("potential sign", q, atoms, i, j, witness))

    residues = tuple(-F(k, q) - z for (u, k), z in zip(atoms, potentials))
    scale = 1
    for (u, k), residue in zip(atoms, residues):
        scale = lcm(scale, u)
        scale = lcm(scale, residue.denominator)

    integral_residues = []
    moduli = []
    for (u, k), residue in zip(atoms, residues):
        scaled = scale * residue
        require(scaled.denominator == 1, ("nonintegral scaled residue", scaled))
        integral_residues.append(scaled.numerator)
        moduli.append(scale // u)

    x = integral_residues[0] % moduli[0]
    modulus = moduli[0]
    for residue, next_modulus in zip(integral_residues[1:], moduli[1:]):
        x, modulus = merge_congruences(x, modulus, residue, next_modulus)
    shift = F(x, scale)
    centres = tuple(z + shift for z in potentials)

    for centre, (u, k) in zip(centres, atoms):
        tooth = u * (centre + F(k, q))
        require(tooth.denominator == 1, ("not a tooth centre", q, atoms, tooth))

    left = max(centre - F(1, 14 * atom[0]) for centre, atom in zip(centres, atoms))
    right = min(centre + F(1, 14 * atom[0]) for centre, atom in zip(centres, atoms))
    require(left < right, ("Helly failed", q, atoms, witness, left, right))
    t = mod_one((left + right) / 2)
    require(all(atom_fires(q, atom, t) for atom in atoms),
            ("reduced phase failed", q, atoms, t))
    common_gcd = 0
    for u, k in atoms:
        common_gcd = gcd(common_gcd, u)
    return t, left, right, common_gcd


def floor_fraction(x):
    return x.numerator // x.denominator


def open_interval_mod_one_contains(t, left, right):
    lower = left - t
    upper = right - t
    candidate = floor_fraction(lower) + 1
    return candidate < upper


def audit_constructive_locus():
    """Compare the full interval union, not merely nonemptiness."""
    families = samples = cochains_seen = repeated_families = 0
    for q in range(2, 9):
        atoms = atom_universe(q, range(1, 5))
        for rank in range(1, 4):
            for family in combinations(atoms, rank):
                witnesses = all_cochains(q, family)
                geometries = []
                for witness in witnesses:
                    t0, left, right, common_gcd = reconstruct_by_crt(q, family, witness)
                    geometries.append((left, right, common_gcd))
                    cochains_seen += 1
                test_points = event_endpoints(q, family) + event_midpoints(q, family)
                for t in test_points:
                    direct = all(atom_fires(q, atom, t) for atom in family)
                    reconstructed = any(
                        open_interval_mod_one_contains(
                            t,
                            left + F(n, common_gcd),
                            right + F(n, common_gcd),
                        )
                        for left, right, common_gcd in geometries
                        for n in range(common_gcd)
                    )
                    require(direct == reconstructed,
                            ("constructive locus mismatch", q, family, t, direct, geometries))
                    samples += 1
                if len({u for u, k in family}) < rank:
                    repeated_families += 1
                families += 1
    return families, samples, cochains_seen, repeated_families


def atom_universe(q, speeds):
    return tuple(
        (u, k)
        for u in speeds
        for k in range(image_order(q, u))
    )


def audit_atom_tuples():
    total = positive = repeated = core_typed = 0
    rank_profile = []
    for q in range(2, 10):
        atoms = atom_universe(q, range(1, 7))
        for rank in range(1, 4):
            local = local_positive = local_repeated = 0
            for family in combinations(atoms, rank):
                direct = direct_common_phase(q, family)
                witness = cochain_witness(q, family)
                star = star_compressed_witness(q, family)
                require((direct is not None) == (witness is not None),
                        ("cochain iff mismatch", q, family, direct, witness))
                require((witness is not None) == (star is not None),
                        ("star compression mismatch", q, family, witness, star))
                if witness is not None:
                    reconstruct_by_crt(q, family, witness)
                    positive += 1
                    local_positive += 1
                if len({u for u, k in family}) < rank:
                    repeated += 1
                    local_repeated += 1
                if any(u % q == 0 for u, k in family):
                    core_typed += 1
                total += 1
                local += 1
            rank_profile.append((q, rank, local, local_positive, local_repeated))

    # A smaller exhaustive rank-four universe supplies independent triangle
    # interaction beyond a single triangle.
    for q in range(2, 9):
        atoms = atom_universe(q, range(1, 4))
        local = local_positive = local_repeated = 0
        for family in combinations(atoms, 4):
            direct = direct_common_phase(q, family)
            witness = cochain_witness(q, family)
            star = star_compressed_witness(q, family)
            require((direct is not None) == (witness is not None),
                    ("rank-four mismatch", q, family, direct, witness))
            require((witness is not None) == (star is not None),
                    ("rank-four star mismatch", q, family, witness, star))
            if witness is not None:
                reconstruct_by_crt(q, family, witness)
                positive += 1
                local_positive += 1
            if len({u for u, k in family}) < 4:
                repeated += 1
                local_repeated += 1
            if any(u % q == 0 for u, k in family):
                core_typed += 1
            total += 1
            local += 1
        rank_profile.append((q, 4, local, local_positive, local_repeated))

    # Explicitly audit q|u beyond the bounded small-speed universe.
    extra = 0
    for q in range(2, 13):
        atoms = atom_universe(q, (q, 2 * q, q + 1))
        for rank in range(1, 4):
            for family in combinations(atoms, rank):
                direct = direct_common_phase(q, family)
                witness = cochain_witness(q, family)
                star = star_compressed_witness(q, family)
                require((direct is not None) == (witness is not None),
                        ("core typing mismatch", q, family, direct, witness))
                require((witness is not None) == (star is not None),
                        ("core star mismatch", q, family, witness, star))
                if witness is not None:
                    reconstruct_by_crt(q, family, witness)
                extra += 1
    return total, positive, repeated, core_typed, extra, tuple(rank_profile)


def audit_capacity_endpoints():
    profile = []
    checks = 0
    for m in range(1, 43):
        endpoints = set()
        for r in range(m):
            centre = -F(r, m)
            endpoints.add(mod_one(centre - F(1, 14)))
            endpoints.add(mod_one(centre + F(1, 14)))
        ordered = sorted(endpoints)
        samples = list(ordered)
        for i, left in enumerate(ordered):
            right = ordered[(i + 1) % len(ordered)]
            if i + 1 == len(ordered):
                right += 1
            samples.append(mod_one((left + right) / 2))
        counts = {
            sum(14 * integer_distance(theta + F(r, m)) < 1 for r in range(m))
            for theta in samples
        }
        ceiling = (m + 6) // 7
        require(counts == {ceiling - 1, ceiling},
                ("capacity endpoint mismatch", m, counts, ceiling))
        profile.append((m, tuple(sorted(counts))))
        checks += len(samples)

    # Verify atom pullback and q|u multiplicity for actual multipliers.
    pullbacks = 0
    for q in range(2, 21):
        for u in range(1, 2 * q + 1):
            m = image_order(q, u)
            for numerator in range(0, 28):
                t = F(numerator, 28 * u)
                firing = tuple(k for k in range(m) if atom_fires(q, (u, k), t))
                mask = 0
                for k in firing:
                    mask |= atom_mask(q, (u, k))
                direct_mask = sum(
                    1 << sheet
                    for sheet in range(q)
                    if 14 * integer_distance(u * (t + F(sheet, q))) < 1
                )
                require(mask == direct_mask, ("pullback mismatch", q, u, t))
                require(mask.bit_count() == gcd(q, u) * len(firing),
                        ("atom multiplicity mismatch", q, u, t))
                pullbacks += 1
    return checks, pullbacks, tuple(profile)


def owner_cover_direct(q, owners):
    atoms = atom_universe(q, owners)
    full = (1 << q) - 1
    for t in event_midpoints(q, atoms):
        mask = 0
        for atom in atoms:
            if atom_fires(q, atom, t):
                mask |= atom_mask(q, atom)
        if mask == full:
            return t
    return None


def owner_atom_options(q, u):
    m = image_order(q, u)
    cap = (m + 6) // 7
    return tuple(
        tuple((u, k) for k in labels)
        for size in range(cap + 1)
        for labels in combinations(range(m), size)
    )


def owner_cover_cochain(q, owners, cache):
    full = (1 << q) - 1
    for grouped in product(*(owner_atom_options(q, u) for u in owners)):
        family = tuple(atom for group in grouped for atom in group)
        if not family:
            continue
        mask = 0
        for atom in family:
            mask |= atom_mask(q, atom)
        if mask != full:
            continue
        if family not in cache:
            cache[family] = cochain_witness(q, family)
        if cache[family] is not None:
            return family
    return None


def audit_owner_covers():
    checks = positives = 0
    profile = []
    for q in range(2, 9):
        cache = {}
        base = (1, 2, 3, 4)
        for size in range(1, len(base) + 1):
            for owners in combinations(base, size):
                direct = owner_cover_direct(q, owners)
                family = owner_cover_cochain(q, owners, cache)
                require((direct is not None) == (family is not None),
                        ("owner cover mismatch", q, owners, direct, family))
                checks += 1
                positives += direct is not None
                profile.append((q, owners, direct is not None))
    return checks, positives, tuple(profile)


def q6_minimal_clutter():
    q = 6
    shapes = tuple(
        [("cell", k, 1 << k) for k in range(6)]
        + [("column", k, sum(1 << j for j in range(6) if j % 3 == k)) for k in range(3)]
        + [("row", k, sum(1 << j for j in range(6) if j % 2 == k)) for k in range(2)]
    )
    full = (1 << q) - 1
    minimal = []
    for size in range(1, len(shapes) + 1):
        for indices in combinations(range(len(shapes)), size):
            mask = 0
            for i in indices:
                mask |= shapes[i][2]
            if mask != full:
                continue
            irredundant = True
            for removed in indices:
                reduced = 0
                for j in indices:
                    if j != removed:
                        reduced |= shapes[j][2]
                if reduced == full:
                    irredundant = False
                    break
            if irredundant:
                minimal.append(indices)
    weights = [sum(shapes[i][2].bit_count() for i in indices) for indices in minimal]
    require(len(minimal) == 23, ("q6 clutter count", len(minimal)))
    require(weights.count(6) == 11, ("q6 partition count", weights))
    require(weights.count(7) == 6 and weights.count(8) == 6,
            ("q6 overlap strata", weights))
    return len(minimal), tuple(sorted((weight, weights.count(weight)) for weight in set(weights)))


def q6_correlated_controls():
    rows = []
    for owners in ((2, 8, 14), (2, 8, 10)):
        pairwise = closing = 0
        for labels in product(range(3), repeat=3):
            family = tuple(zip(owners, labels))
            mask = 0
            for atom in family:
                mask |= atom_mask(6, atom)
            if mask != 63:
                continue
            if all(allowed_gaps(6, family[i], family[j]) for i, j in combinations(range(3), 2)):
                pairwise += 1
            if cochain_witness(6, family) is not None:
                closing += 1
        direct = owner_cover_direct(6, owners) is not None
        require(direct == (closing > 0), ("q6 correlated mismatch", owners, direct, closing))
        rows.append((owners, pairwise, closing, direct))
    require(rows[0][1] > 0 and rows[0][2] == 0, ("hostile lost", rows[0]))
    require(rows[1][2] > 0, ("positive lost", rows[1]))
    return tuple(rows)


def source_hygiene():
    source = Path(__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    float_nodes = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(tree)
    )
    require(assert_nodes == 0 and float_nodes == 0, ("source hygiene", assert_nodes, float_nodes))
    return assert_nodes, float_nodes


def main():
    hygiene = source_hygiene()
    capacity = audit_capacity_endpoints()
    tuples = audit_atom_tuples()
    owner_covers = audit_owner_covers()
    constructive_locus = audit_constructive_locus()
    clutter = q6_minimal_clutter()
    q6_controls = q6_correlated_controls()
    semantic_fields = (
        capacity,
        tuples,
        owner_covers,
        constructive_locus,
        clutter,
        q6_controls,
        ("pair_congruence", "triangle_coboundary", "generalized_CRT",
         "coherent_real_lifts", "strict_open_Helly", "base_star_compression"),
        ("body_pointwise", "B_subset_A_C"),
        ("body_aligned_grid", "B_minus_A_C_subset_Gamma_D"),
    )
    semantic = sha256(repr(semantic_fields).encode()).hexdigest()
    if EXPECTED_SEMANTIC != "TO_BE_FROZEN":
        require(semantic == EXPECTED_SEMANTIC, ("semantic drift", semantic))

    source_bytes = Path(__file__).read_bytes().replace(b"\r\n", b"\n")
    print("ALL-Q ATOMIZED COSET-COCHAIN SECOND INDEPENDENT AUDIT")
    print("source_sha256_lf=" + sha256(source_bytes).hexdigest())
    print("hygiene_assert_nodes=%d;float_literals=%d" % hygiene)
    print("capacity_phase_checks=%d;pullback_checks=%d;profile=%r" % capacity)
    print("atom_tuple_checks=%d;positive=%d;repeated_owner=%d;core_typed=%d;extra_core_checks=%d" % tuples[:5])
    print("atom_tuple_rank_profile=%r" % (tuples[5],))
    print("owner_cover_checks=%d;positive=%d;profile=%r" % owner_covers)
    print("constructive_locus_families=%d;samples=%d;cochains=%d;repeated_owner_families=%d;coordinate=y=q*t" % constructive_locus)
    print("q6_minimal_clutter=%d;weight_profile=%r" % clutter)
    print("q6_correlated_controls=%r" % (q6_controls,))
    print("body_sidecars=pointwise:B_subset_A_C;aligned_grid:B_minus_A_C_subset_Gamma_D")
    print("semantic_sha256=" + semantic)
    print("verdict=PASS")


if __name__ == "__main__":
    main()

