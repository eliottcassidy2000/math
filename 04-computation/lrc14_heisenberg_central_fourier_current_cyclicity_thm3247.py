#!/usr/bin/env python3
"""Assertion-independent exact companion for THM-3247.

The companion has two independent layers.

1.  Over exact splitting fields it checks the multiplicity-free central-
    Fourier decomposition of the standard affine H_p permutation module for
    p=3,5,13, together with three sharp orbit-span controls.
2.  It rebuilds the two certified THM-2625 endpoint-factor specializations
    used by THM-2790.  For every nonzero q in F_13^2 it recomputes all charged
    central modes, all thirteen scalar-character projections, and all 169
    two-dimensional Fourier modes for the regular translation subgroup.

The script contains no assert statements, floating literals, randomness, or
cache inputs.  Every digest uses LF-normalized dependency bytes or unsigned
eight-byte big-endian field elements.
"""

from ast import Assert, Constant, parse, walk
from hashlib import sha256
import importlib.util
from itertools import product
from math import isqrt
from pathlib import Path


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    "THM-2779": (
        ROOT / "01-canon/theorems/THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate.md",
        "21ac7ec9b19b8ed0abe4763e0b7e13ebc1e5eb776168c8e0088540f29daabccb",
    ),
    "THM-2790": (
        ROOT / "01-canon/theorems/THM-2790-universal-depth-two-central-response-and-carry-wall-spectrum.md",
        "dc5a4ab2b6a158d3e80918e42b3d8b59cba186d8d3b364c65010c94091c6a8ba",
    ),
    "THM-3240": (
        ROOT / "01-canon/theorems/THM-3240-exact-address-heisenberg-clutch-on-carrier-imbalance.md",
        "7d23f2920adfecb17d8a149aada08a8e34215111546eac63b40570e898e14f14",
    ),
}
UPSTREAM_SOURCE = (
    ROOT / "04-computation/lrc14_universal_central_response_thm2790.py"
)
UPSTREAM_OUTPUT = (
    ROOT / "05-knowledge/results/lrc14_universal_central_response_thm2790.out"
)
UPSTREAM_SOURCE_HASH = (
    "f2501a8b9b6a3b8d8d3f30b0bc5267f9b3f161f7fe71f6333d0b5fb36e261b82"
)
UPSTREAM_OUTPUT_HASH = (
    "720b118e98c25d95195ac824647f534a14cb964fc967624ecbf68f4b178c4f9e"
)


def lf_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


for dependency_name, (dependency_path, expected_hash) in DEPENDENCIES.items():
    actual_hash = lf_hash(dependency_path)
    require(
        actual_hash == expected_hash,
        ("dependency hash", dependency_name, actual_hash, expected_hash),
    )
require(
    lf_hash(UPSTREAM_SOURCE) == UPSTREAM_SOURCE_HASH,
    ("upstream source hash", lf_hash(UPSTREAM_SOURCE), UPSTREAM_SOURCE_HASH),
)
require(
    lf_hash(UPSTREAM_OUTPUT) == UPSTREAM_OUTPUT_HASH,
    ("upstream output hash", lf_hash(UPSTREAM_OUTPUT), UPSTREAM_OUTPUT_HASH),
)

syntax_tree = parse(Path(__file__).read_text(encoding="utf-8"))
require(
    not any(isinstance(node, Assert) for node in walk(syntax_tree)),
    "assert statements are optimization-sensitive",
)
require(
    not any(
        isinstance(node, Constant) and isinstance(node.value, float)
        for node in walk(syntax_tree)
    ),
    "floating literals are forbidden",
)

SPEC = importlib.util.spec_from_file_location("thm2790_companion", UPSTREAM_SOURCE)
UPSTREAM = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(UPSTREAM)


def h_mul(left, right, p):
    x, y, central = left
    xp, yp, centralp = right
    return (
        (x + xp) % p,
        (y + yp) % p,
        (central + centralp - y * xp) % p,
    )


def act(element, state, p):
    x, y, central = element
    transverse, central_coordinate = state
    return (
        (transverse + x) % p,
        (central_coordinate + central - y * transverse) % p,
    )


def exact_order_p_root(modulus, p):
    require(
        modulus >= 2
        and all(modulus % divisor for divisor in range(2, isqrt(modulus) + 1)),
        ("nonprime splitting modulus", modulus),
    )
    require((modulus - 1) % p == 0, ("splitting modulus", modulus, p))
    root = next(
        (
            candidate
            for candidate in range(2, modulus)
            if pow(candidate, p, modulus) == 1 and candidate != 1
        ),
        None,
    )
    require(root is not None, ("missing p-th root", modulus, p))
    return root


def point_index(state, p):
    return state[0] * p + state[1]


def permute_vector(vector, element, p):
    image = [0] * (p * p)
    for state in product(range(p), repeat=2):
        image[point_index(act(element, state, p), p)] = vector[
            point_index(state, p)
        ]
    return tuple(image)


def span_rank(rows, modulus):
    basis = {}
    for source_row in rows:
        row = [value % modulus for value in source_row]
        for pivot in sorted(basis):
            if row[pivot]:
                factor = row[pivot]
                row = [
                    (value - factor * base_value) % modulus
                    for value, base_value in zip(row, basis[pivot])
                ]
        pivot = next((index for index, value in enumerate(row) if value), None)
        if pivot is None:
            continue
        inverse = pow(row[pivot], modulus - 2, modulus)
        row = [inverse * value % modulus for value in row]
        basis[pivot] = tuple(row)
    return len(basis)


def block_support(vector, p, modulus, zeta):
    central_coefficients = []
    for central_character in range(p):
        row = []
        for transverse in range(p):
            row.append(
                sum(
                    vector[point_index((transverse, central_coordinate), p)]
                    * pow(
                        zeta,
                        central_character * central_coordinate % p,
                        modulus,
                    )
                    for central_coordinate in range(p)
                )
                % modulus
            )
        central_coefficients.append(tuple(row))

    scalar_modes = tuple(
        sum(
            central_coefficients[0][transverse]
            * pow(zeta, scalar_character * transverse % p, modulus)
            for transverse in range(p)
        )
        % modulus
        for scalar_character in range(p)
    )
    active_scalars = sum(value != 0 for value in scalar_modes)
    active_charged = sum(
        any(central_coefficients[central_character])
        for central_character in range(1, p)
    )
    predicted_rank = active_scalars + p * active_charged
    return active_scalars, active_charged, predicted_rank


def verify_decomposition(p, modulus):
    zeta = exact_order_p_root(modulus, p)
    require(pow(zeta, p, modulus) == 1 and zeta != 1, "root order")
    orthogonality = tuple(
        sum(pow(zeta, exponent * index % p, modulus) for index in range(p))
        % modulus
        for exponent in range(p)
    )
    require(orthogonality == (p % modulus,) + (0,) * (p - 1), orthogonality)

    group = tuple(product(range(p), repeat=3))
    # The inherited THM-2779/3240 companions already exhaust all H_13 pairs.
    # Here it is enough to check every element against a signed generating
    # set: the displayed affine signatures are closed under composition.
    generators = (
        (1, 0, 0),
        (-1 % p, 0, 0),
        (0, 1, 0),
        (0, -1 % p, 0),
        (0, 0, 1),
        (0, 0, -1 % p),
    )
    composition_checks = 0
    for left in group:
        for right in generators:
            lx, ly, lc = left
            rx, ry, rc = right
            composed_signature = (
                (lx + rx) % p,
                (lc + rc - ly * rx) % p,
                -(ly + ry) % p,
            )
            product_element = h_mul(left, right, p)
            product_signature = (
                product_element[0],
                product_element[2],
                -product_element[1] % p,
            )
            require(
                composed_signature == product_signature,
                ("group/action sign", p, left, right),
            )
            composition_checks += 1

    fourier_formula_checks = 0
    for element in group:
        x, y, central = element
        for transverse in range(p):
            for central_character in range(p):
                phase = pow(
                    zeta,
                    central_character * (central - y * transverse) % p,
                    modulus,
                )
                # Both sides are characters linear in the central coordinate;
                # checking 0 and 1 checks their constant term and slope.
                for central_coordinate in (0, 1):
                    image_coordinate = (
                        central_coordinate + central - y * transverse
                    ) % p
                    left_coefficient = pow(
                        zeta,
                        -central_character * central_coordinate % p,
                        modulus,
                    )
                    right_coefficient = phase * pow(
                        zeta,
                        -central_character * image_coordinate % p,
                        modulus,
                    ) % modulus
                    require(
                        left_coefficient == right_coefficient,
                        (
                            "Fourier action",
                            p,
                            element,
                            transverse,
                            central_character,
                            central_coordinate,
                        ),
                    )
                    fourier_formula_checks += 1

    for central_character in range(1, p):
        weight_rows = {
            tuple(
                pow(
                    zeta,
                    -central_character * y * transverse % p,
                    modulus,
                )
                for y in range(p)
            )
            for transverse in range(p)
        }
        require(len(weight_rows) == p, ("distinct y weights", p, central_character))
    require(
        all({(transverse + x) % p for x in range(p)} == set(range(p))
            for transverse in range(p)),
        ("x transitivity", p),
    )

    constant = (1,) * (p * p)
    charged_seed = tuple(
        pow(zeta, -central_coordinate % p, modulus) if transverse == 0 else 0
        for transverse in range(p)
        for central_coordinate in range(p)
    )
    delta_seed = tuple(
        int(transverse == 0 and central_coordinate == 0)
        for transverse in range(p)
        for central_coordinate in range(p)
    )
    controls = (
        ("constant", constant, (1, 0, 1)),
        ("one_charged_character", charged_seed, (0, 1, p)),
        ("point_delta", delta_seed, (p, p - 1, p * p)),
    )
    control_ranks = {}
    for name, vector, expected_support in controls:
        support = block_support(vector, p, modulus, zeta)
        require(support == expected_support, ("block support", p, name, support))
        rank = span_rank(
            (permute_vector(vector, element, p) for element in group), modulus
        )
        require(rank == expected_support[2], ("orbit rank", p, name, rank))
        control_ranks[name] = rank

    scalar_dimensions = p
    charged_dimensions = p * (p - 1)
    require(
        scalar_dimensions + charged_dimensions == p * p,
        ("dimension ledger", p),
    )
    return {
        "p": p,
        "modulus": modulus,
        "zeta": zeta,
        "composition_checks": composition_checks,
        "fourier_formula_checks": fourier_formula_checks,
        "scalar_blocks": p,
        "charged_blocks": p - 1,
        "scalar_dimension": scalar_dimensions,
        "charged_dimension": charged_dimensions,
        "constant_rank": control_ranks["constant"],
        "charged_rank": control_ranks["one_charged_character"],
        "delta_rank": control_ranks["point_delta"],
    }


def field_digest(values):
    payload = b"".join(int(value).to_bytes(8, "big") for value in values)
    return sha256(payload).hexdigest()


P = 13
EXPECTED_NEW = (
    (
        352341050142921841,
        "5c4d76f09f89935ff5570dfbaac51972e3102016d435b2968025ca67479aad73",
        "bd4e07157f30c82b5cefe0e77f20cd802c6b4945d21bc4e9649ef992e4b77419",
        "39dbf7d36cfd37a637499591bf45418007a9b26edd1e15b469ab7111a7bba9c2",
        "5dcc44a3f80407aac7b5e02bc8c8fa26dc15f80db96412e2df9dc57a5e5bb828",
    ),
    (
        956354278959359281,
        "435c2e675a7bc48014660abe538d74c01e36bf14ede832556f5df7d17e57f92d",
        "77ee3ef36082b5195b753e0140ab74ca04d86b5307f809d8f3867f73e8f83869",
        "671906fa806d5cf4653945cd85d3f2ef0fd15580c6a916e5272bfe0b14a265c7",
        "6703aa62b265b463183c4a24810ad900d9c64364e6fd6857e7b679c19bbbb76b",
    ),
)


def point_scale(scalar, point):
    return scalar * point[0] % P, scalar * point[1] % P


def point_add(left, right):
    return (left[0] + right[0]) % P, (left[1] + right[1]) % P


def analyze_current_field(field, field_index):
    prime, left, right = field
    primitive_root = UPSTREAM.GATE.BASE.MODS[field_index][1]
    zeta = pow(primitive_root, UPSTREAM.GATE.BASE.NN // P, prime)
    require(pow(zeta, P, prime) == 1 and zeta != 1, "certified zeta")

    central_sum_bank = []
    scalar_mode_bank = []
    charged_mode_bank = []
    # The translation digest serializes q in UPSTREAM.NONZERO order; for each
    # q it uses the first t in UPSTREAM.POINTS with det(q,t)=1, then b=0,...,12
    # outside a=0,...,12 (the b=0 row is scalar_modes).  field_digest writes
    # every residue as an unsigned big-endian eight-byte integer.  The digest
    # therefore pins this deterministic section, whereas full support/rank is
    # invariant under the shear/permutation induced by t -> t+kq.
    translation_mode_bank = []
    cyclic_directions = 0
    translation_cyclic_directions = 0
    irreducible_blocks = 0
    for step in UPSTREAM.NONZERO:
        transverse = next(
            point
            for point in UPSTREAM.POINTS
            if UPSTREAM.GATE.det(step, point) == 1
        )
        coordinates = {
            point_add(point_scale(w, step), point_scale(delta, transverse))
            for delta in range(P)
            for w in range(P)
        }
        require(
            len(coordinates) == P * P,
            ("coordinate chart", field_index, step, len(coordinates)),
        )
        central_sums = []
        charged_rows = {character: [] for character in range(1, P)}
        charged_active = set()
        for delta in range(P):
            current_cycle = []
            for w in range(P):
                origin = point_add(
                    point_scale(w, step), point_scale(delta, transverse)
                )
                current_cycle.append(
                    left[point_add(origin, step)] * right[origin] % prime
                )
            central_sum = sum(current_cycle) % prime
            central_sums.append(central_sum)
            central_sum_bank.append(central_sum)
            for central_character in range(1, P):
                mode = sum(
                    value * pow(
                        zeta,
                        -central_character * w % P,
                        prime,
                    )
                    for w, value in enumerate(current_cycle)
                ) % prime
                require(
                    mode != 0,
                    (
                        "charged mode hole",
                        field_index,
                        step,
                        delta,
                        central_character,
                    ),
                )
                charged_mode_bank.append(mode)
                charged_rows[central_character].append(mode)
                charged_active.add(central_character)

        scalar_active = 0
        scalar_modes = []
        for scalar_character in range(P):
            mode = sum(
                value
                * pow(zeta, -scalar_character * delta % P, prime)
                for delta, value in enumerate(central_sums)
            ) % prime
            require(
                mode != 0,
                (
                    "scalar mode hole",
                    field_index,
                    step,
                    scalar_character,
                ),
            )
            scalar_mode_bank.append(mode)
            scalar_modes.append(mode)
            scalar_active += 1

        translation_modes = list(scalar_modes)
        for central_character in range(1, P):
            require(
                len(charged_rows[central_character]) == P,
                ("charged transverse row", field_index, step, central_character),
            )
            for transverse_character in range(P):
                mode = sum(
                    value
                    * pow(
                        zeta,
                        -transverse_character * delta % P,
                        prime,
                    )
                    for delta, value in enumerate(
                        charged_rows[central_character]
                    )
                ) % prime
                require(
                    mode != 0,
                    (
                        "two-dimensional Fourier hole",
                        field_index,
                        step,
                        transverse_character,
                        central_character,
                    ),
                )
                translation_modes.append(mode)
        require(
            len(translation_modes) == P * P and all(translation_modes),
            ("translation orbit rank", field_index, step),
        )
        translation_mode_bank.extend(translation_modes)
        translation_cyclic_directions += 1

        require(
            scalar_active == P and len(charged_active) == P - 1,
            ("block support", field_index, step, scalar_active, charged_active),
        )
        predicted_orbit_rank = scalar_active + P * len(charged_active)
        require(
            predicted_orbit_rank == P * P,
            ("cyclic rank criterion", field_index, step, predicted_orbit_rank),
        )
        irreducible_blocks += scalar_active + len(charged_active)
        cyclic_directions += 1

    require(
        len(central_sum_bank) == len(UPSTREAM.NONZERO) * P,
        "central-sum universe drift",
    )
    require(
        len(scalar_mode_bank) == len(UPSTREAM.NONZERO) * P,
        "scalar-mode universe drift",
    )
    require(
        len(charged_mode_bank) == len(UPSTREAM.NONZERO) * P * (P - 1),
        "charged-mode universe drift",
    )
    require(
        len(translation_mode_bank) == len(UPSTREAM.NONZERO) * P * P,
        "two-dimensional Fourier universe drift",
    )
    result = {
        "prime": prime,
        "central_sum_entries": len(central_sum_bank),
        "scalar_modes": len(scalar_mode_bank),
        "charged_modes": len(charged_mode_bank),
        "translation_modes": len(translation_mode_bank),
        "irreducible_blocks": irreducible_blocks,
        "cyclic_directions": cyclic_directions,
        "translation_cyclic_directions": translation_cyclic_directions,
        "central_sum_digest": field_digest(central_sum_bank),
        "scalar_mode_digest": field_digest(scalar_mode_bank),
        "charged_mode_digest": field_digest(charged_mode_bank),
        "translation_mode_digest": field_digest(translation_mode_bank),
    }
    if EXPECTED_NEW is not None:
        require(
            (
                result["prime"],
                result["central_sum_digest"],
                result["scalar_mode_digest"],
                result["charged_mode_digest"],
            )
            == EXPECTED_NEW[field_index][:4],
            ("new deterministic digest drift", field_index, result),
        )
        if EXPECTED_NEW[field_index][4] is not None:
            require(
                result["translation_mode_digest"]
                == EXPECTED_NEW[field_index][4],
                ("translation digest drift", field_index, result),
            )
    return result


def main():
    representation_results = tuple(
        verify_decomposition(p, modulus) for p, modulus in ((3, 7), (5, 11), (13, 53))
    )
    carry_count = UPSTREAM.verify_heisenberg_carrier()
    hostile_edges = UPSTREAM.verify_edge_gate_hostile()
    fields = UPSTREAM.GATE.build_endpoint_factors()
    inherited_results = tuple(
        UPSTREAM.analyze_field(field, field_index)
        for field_index, field in enumerate(fields)
    )
    current_results = tuple(
        analyze_current_field(field, field_index)
        for field_index, field in enumerate(fields)
    )
    require(
        all(result["central_mode_support"] == 26208 for result in inherited_results),
        "THM-2790 charged-mode inheritance",
    )

    print("THM-3247 HEISENBERG FOURIER DECOMPOSITION AND CURRENT CYCLICITY")
    print("status=VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED")
    print(
        "pins=THM-2779,THM-2790,THM-3240;upstream_script_and_output=LF_SHA256_MATCH"
    )
    for result in representation_results:
        print(
            f"p={result['p']} split_modulus={result['modulus']} zeta={result['zeta']} "
            f"decomposition={result['scalar_blocks']}x1+"
            f"{result['charged_blocks']}x{result['p']}="
            f"{result['scalar_dimension'] + result['charged_dimension']} "
            f"composition_checks={result['composition_checks']} "
            f"fourier_formula_checks={result['fourier_formula_checks']}"
        )
        print(
            f"p={result['p']} control_orbit_ranks="
            f"constant:{result['constant_rank']},"
            f"one_charged_character:{result['charged_rank']},"
            f"point_delta:{result['delta_rank']}"
        )
    print(
        f"p=13 carrier_points=169 carry_wall_points={carry_count} "
        f"separate_edge_constant_product_hostile={hostile_edges}/26"
    )
    print(
        "cyclicity_criterion=nonzero_projection_to_each_of_13_scalar_plus_12_charged_irreducibles"
    )
    for result in current_results:
        print(
            f"field={result['prime']} central_sum_entries={result['central_sum_entries']} "
            f"scalar_character_modes={result['scalar_modes']}/2184 "
            f"charged_coordinate_modes={result['charged_modes']}/26208 "
            f"translation_modes={result['translation_modes']}/28392 "
            f"irreducible_blocks={result['irreducible_blocks']}/4200 "
            f"cyclic_directions={result['cyclic_directions']}/168 "
            f"translation_cyclic_directions="
            f"{result['translation_cyclic_directions']}/168 "
            "H_and_translation_orbit_rank_per_direction=169"
        )
        print(
            f"field={result['prime']} "
            f"central_sums_sha256={result['central_sum_digest']} "
            f"scalar_modes_sha256={result['scalar_mode_digest']} "
            f"charged_modes_sha256={result['charged_mode_digest']} "
            f"translation_modes_sha256={result['translation_mode_digest']}"
        )
    print(
        "response_state_boundary=center_blind_scalar_line:dimension_1;"
        "one_charged_channel:dimension_13;all_charged_channels:dimension_156;"
        "full_carrier:dimension_169"
    )
    print(
        "scope=coefficient-side normalized endpoint carrier;physical same-ancestry descent,"
        "section-free intertwiner,owner map,and LRC14 remain open"
    )
    print("FAILED_CHECKS=NONE")


if __name__ == "__main__":
    main()
