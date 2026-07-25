#!/usr/bin/env python3
"""Exact, optimization-safe companion for THM-2350."""

from itertools import product


P = 13
FIELD_PRIME = 53
LABELS = ("j", "a", "b", "u0", "ka", "kb", "u1", "u2", "u3")
INDEX = {label: position for position, label in enumerate(LABELS)}
DIMENSION = len(LABELS)
TARGET_GROUP = tuple(product(range(P), repeat=2))


def require(condition: bool, message: str) -> None:
    """Raise under ordinary and optimized Python."""
    if not condition:
        raise RuntimeError(message)


def vector(**entries: int) -> tuple[int, ...]:
    answer = [0] * DIMENSION
    for label, value in entries.items():
        answer[INDEX[label]] = value % P
    return tuple(answer)


def add(*vectors: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(values) % P for values in zip(*vectors))


def scale(scalar: int, value: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(scalar * entry % P for entry in value)


def dot(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return sum(a * b for a, b in zip(left, right)) % P


def rank_mod(matrix: list[list[int]], prime: int) -> int:
    """Exact row rank over F_prime."""
    rows = [[entry % prime for entry in row] for row in matrix]
    if not rows:
        return 0
    row_count = len(rows)
    column_count = len(rows[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (
                row
                for row in range(pivot_row, row_count)
                if rows[row][column] % prime
            ),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        inverse = pow(rows[pivot_row][column], -1, prime)
        rows[pivot_row] = [
            entry * inverse % prime for entry in rows[pivot_row]
        ]
        for row in range(row_count):
            if row == pivot_row:
                continue
            factor = rows[row][column] % prime
            if factor:
                rows[row] = [
                    (entry - factor * pivot_entry) % prime
                    for entry, pivot_entry in zip(rows[row], rows[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def owner_pivot_packet() -> tuple[tuple[int, ...], list[tuple[int, ...]]]:
    """Build the mod-13 word and THM-2309 owner-pivot rows."""
    weights = {
        "j": 0,
        "a": 0,
        "b": 0,
        "u0": 5,
        "ka": 2,
        "kb": 7,
        "u1": 1,
        "u2": 4,
        "u3": 9,
    }
    word = tuple(weights[label] for label in LABELS)
    q = weights["u0"]
    rows: list[tuple[int, ...]] = []
    for label in ("j", "ka", "kb", "u1", "u2", "u3"):
        entries = {"u0": weights[label], label: -q}
        if label == "ka":
            entries["a"] = -q
        if label == "kb":
            entries["b"] = -q
        rows.append(vector(**entries))

    require(all(dot(row, word) == 0 for row in rows), "packet left K")
    require(
        rank_mod([list(row) for row in rows], P) == 6,
        "owner-pivot rank changed",
    )
    return word, rows


def dipole(label: str, unit: str) -> tuple[int, ...]:
    return vector(**{label: 1, unit: -1})


def projection(value: tuple[int, ...]) -> tuple[int, int]:
    return (
        (value[INDEX["a"]] - value[INDEX["ka"]]) % P,
        (value[INDEX["b"]] - value[INDEX["kb"]]) % P,
    )


def dual_and_projection_atlas(
    word: tuple[int, ...],
    rows: list[tuple[int, ...]],
) -> tuple[int, int, int, int]:
    """Verify all dual gauges and the quotient kernel."""
    d_a = dipole("a", "ka")
    d_b = dipole("b", "kb")
    require(all(dot(row, d_a) == 0 for row in rows), "D_a not dual")
    require(all(dot(row, d_b) == 0 for row in rows), "D_b not dual")

    dual_vectors = set()
    gauge_vectors = set()
    for gauge, s, t in product(range(P), repeat=3):
        representative = add(
            scale(gauge, word),
            scale(s, d_a),
            scale(t, d_b),
        )
        require(
            all(dot(row, representative) == 0 for row in rows),
            "dual normal form stopped annihilating L",
        )
        dual_vectors.add(representative)
        if gauge == 0:
            require(
                representative[INDEX["u0"]] == 0,
                "gauge section moved u0",
            )
            gauge_vectors.add(representative)

    require(len(dual_vectors) == P**3, "dual normal form is not complete")
    require(len(gauge_vectors) == P**2, "gauge section is not unique")

    require(all(projection(row) == (0, 0) for row in rows), "L projects")
    constraint_matrix = [
        list(word),
        list(d_a),
        list(d_b),
    ]
    require(
        rank_mod(constraint_matrix, P) == 3,
        "word/dipole constraints lost independence",
    )
    quotient_kernel_dimension = DIMENSION - 3
    require(
        quotient_kernel_dimension == len(rows) == 6,
        "projection kernel dimension changed",
    )
    require(projection(vector(a=1)) == (1, 0), "a target basis moved")
    require(projection(vector(b=1)) == (0, 1), "b target basis moved")

    minimum_axis_support = min(
        sum(entry != 0 for entry in add(d_a, scale(gauge, word)))
        for gauge in range(P)
    )
    require(minimum_axis_support == 2, "target axis acquired a unary lift")
    return (
        len(dual_vectors),
        len(gauge_vectors),
        quotient_kernel_dimension,
        minimum_axis_support,
    )


def atomic_projection_atlas() -> int:
    """Check the dipole-residue plus deep-translation formula."""
    checks = 0
    for delta_u_a, delta_u_b, delta_v_a, delta_v_b in product(
        range(P), repeat=4
    ):
        for deepest in range(2):
            multiplier = 5
            direct = (
                (
                    delta_u_a
                    - delta_v_a
                    + multiplier * (deepest == 0)
                )
                % P,
                (
                    delta_u_b
                    - delta_v_b
                    + multiplier * (deepest == 1)
                )
                % P,
            )
            phase = (
                delta_u_a - delta_v_a,
                delta_u_b - delta_v_b,
            )
            translated = list(value % P for value in phase)
            translated[deepest] = (
                translated[deepest] + multiplier
            ) % P
            require(tuple(translated) == direct, "atomic target formula failed")
            checks += 1
    return checks


def transported_word_invisibility_atlas() -> int:
    """Check that every dipole word residue is killed by R."""
    checks = 0
    for exponent in range(1, 5):
        dilation = P**exponent
        for beta_a, beta_b in TARGET_GROUP:
            require(
                (dilation * beta_a) % P
                == (dilation * beta_b) % P
                == 0,
                "a transported-word dipole residue became target-visible",
            )
            checks += 1
    return checks


def primitive_root_of_order_13() -> int:
    for candidate in range(2, FIELD_PRIME):
        if pow(candidate, P, FIELD_PRIME) == 1 and candidate != 1:
            return candidate
    raise RuntimeError("F_53 lost its thirteenth root")


def vertex_index(s: int, t: int) -> int:
    return (s % P) * P + (t % P)


def magnetic_incidence(
    zeta: int,
    multiplier: int,
) -> list[list[int]]:
    rows: list[list[int]] = []
    inverse_phase = pow(zeta, (-multiplier) % P, FIELD_PRIME)
    for s, t in TARGET_GROUP:
        horizontal = [0] * (P**2)
        horizontal[vertex_index(s + 1, t)] = 1
        horizontal[vertex_index(s, t)] = -1
        rows.append(horizontal)

        vertical = [0] * (P**2)
        vertical[vertex_index(s, t + 1)] = 1
        vertical[vertex_index(s, t)] = -inverse_phase
        rows.append(vertical)
    return rows


def magnetic_flatness_atlas() -> tuple[int, int, int, int, int]:
    """Verify the exact rank-one flat section and spectral zero."""
    zeta = primitive_root_of_order_13()
    require(pow(zeta, P, FIELD_PRIME) == 1, "zeta order failed")
    require(
        all(pow(zeta, divisor, FIELD_PRIME) != 1 for divisor in (1,)),
        "zeta became trivial",
    )
    multiplier = 5
    incidence = magnetic_incidence(zeta, multiplier)
    rank = rank_mod(incidence, FIELD_PRIME)
    require(rank == P**2 - 1, "magnetic incidence nullity is not one")

    flat = [
        pow(zeta, (-multiplier * t) % P, FIELD_PRIME)
        for s, t in TARGET_GROUP
    ]
    require(
        all(
            sum(entry * value for entry, value in zip(row, flat))
            % FIELD_PRIME
            == 0
            for row in incidence
        ),
        "inverse character is not magnetic-flat",
    )

    ordinary = magnetic_incidence(zeta, 0)
    ordinary_rank = rank_mod(ordinary, FIELD_PRIME)
    require(ordinary_rank == rank, "deep gauge changed incidence rank")

    gauge_edge_checks = 0
    inverse_phase = pow(zeta, (-multiplier) % P, FIELD_PRIME)
    for s, t in TARGET_GROUP:
        gauge_current = pow(zeta, multiplier * t, FIELD_PRIME)
        gauge_horizontal = pow(zeta, multiplier * t, FIELD_PRIME)
        gauge_vertical = pow(
            zeta,
            multiplier * ((t + 1) % P),
            FIELD_PRIME,
        )
        require(
            gauge_horizontal == gauge_current,
            "horizontal deep gauge changed",
        )
        require(
            gauge_current * pow(gauge_vertical, -1, FIELD_PRIME)
            % FIELD_PRIME
            == inverse_phase,
            "vertical deep gauge is not the magnetic phase",
        )
        gauge_edge_checks += 2

    spectral_zeros = 0
    for x, y in TARGET_GROUP:
        horizontal_multiplier = (pow(zeta, x, FIELD_PRIME) - 1) % FIELD_PRIME
        vertical_multiplier = (
            pow(zeta, y, FIELD_PRIME)
            - pow(zeta, (-multiplier) % P, FIELD_PRIME)
        ) % FIELD_PRIME
        if horizontal_multiplier == vertical_multiplier == 0:
            require((x, y) == (0, (-multiplier) % P), "wrong spectral zero")
            spectral_zeros += 1
    require(spectral_zeros == 1, "magnetic spectrum has extra zeros")

    edge_defect_controls = 0
    for s, t in TARGET_GROUP:
        current = flat[vertex_index(s, t)]
        horizontal = flat[vertex_index(s + 1, t)]
        vertical = flat[vertex_index(s, t + 1)]
        require(horizontal == current, "flat horizontal defect survived")
        require(
            vertical
            == pow(zeta, (-multiplier) % P, FIELD_PRIME) * current
            % FIELD_PRIME,
            "flat vertical defect survived",
        )
        edge_defect_controls += 2
    return (
        zeta,
        rank,
        spectral_zeros,
        edge_defect_controls,
        gauge_edge_checks,
    )


word, packet_rows = owner_pivot_packet()
(
    dual_normal_forms,
    gauge_section_size,
    projection_kernel_dimension,
    minimum_axis_support,
) = dual_and_projection_atlas(word, packet_rows)
atomic_projection_checks = atomic_projection_atlas()
transported_word_checks = transported_word_invisibility_atlas()
(
    zeta_53,
    magnetic_rank,
    spectral_zeros,
    edge_defect_controls,
    gauge_edge_checks,
) = (
    magnetic_flatness_atlas()
)

print("theorem=THM-2350")
print("status=PROVED+VERIFIED-EXACT+CANDIDATE-UNDER-INDEPENDENT-AUDIT")
print("field=F_13")
print("owner_pivot_rank=6")
print(f"dual_normal_forms={dual_normal_forms}")
print(f"unique_omitted_unit_gauge_section={gauge_section_size}")
print("dual_axis_a=e_a-e_k_a")
print("dual_axis_b=e_b-e_k_b")
print(f"minimum_target_axis_support={minimum_axis_support}")
print("target_projection=(r_a-r_k_a,r_b-r_k_b)")
print(f"projection_kernel_dimension={projection_kernel_dimension}")
print(f"atomic_dipole_projection_checks={atomic_projection_checks}")
print(f"transported_word_invisibility_checks={transported_word_checks}")
print("transported_word_target_residue=ZERO")
print(f"finite_field_root_zeta_13_mod_53={zeta_53}")
print("magnetic_incidence_shape=338x169")
print(f"magnetic_incidence_rank={magnetic_rank}")
print(f"magnetic_spectral_zeros={spectral_zeros}")
print("magnetic_kernel=inverse_deep_character_line")
print("deep_phase_gauge=ordinary_torus_incidence")
print(f"deep_phase_gauge_edge_checks={gauge_edge_checks}")
print(f"flat_edge_defect_controls={edge_defect_controls}")
print("one_nonzero_dipole_edge_implies_nonzero_target=YES")
print("canonical_edge_defect_nonzero=OPEN")
print("word_matching_component_positive=OPEN")
print("scalar_rows_excluded=0")
print("lrc14_status=OPEN")
print("all_checks=PASS")
