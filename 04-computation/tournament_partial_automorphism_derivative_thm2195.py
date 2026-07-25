#!/usr/bin/env python3
"""Exact indexing referee for THM-2195's partial-automorphism derivative.

The theorem is proved by partitioning unordered vertex pairs.  This companion
exhausts every labelled tournament of orders 2 through 5, every nonempty
induced vertex set, and every automorphism of that induced subtournament.  It
checks the quotient-pair formula against literal expanded block pairs, then
checks the Hamming, symmetric-difference, orbit-transition, extension-kernel,
and right-invariant-pseudometric descriptions independently.
"""

from itertools import combinations, permutations


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def tournament(q: int, mask: int) -> tuple[tuple[bool, ...], ...]:
    out = [[False] * q for _ in range(q)]
    for bit, (i, j) in enumerate(combinations(range(q), 2)):
        winner, loser = (i, j) if (mask >> bit) & 1 else (j, i)
        out[winner][loser] = True
    return tuple(tuple(row) for row in out)


def automorphisms(
    out: tuple[tuple[bool, ...], ...], vertices: tuple[int, ...]
) -> tuple[tuple[int, ...], ...]:
    answer = []
    for image in permutations(vertices):
        rho = dict(zip(vertices, image))
        if all(
            out[i][j] == out[rho[i]][rho[j]]
            for i in vertices
            for j in vertices
            if i != j
        ):
            answer.append(image)
    return tuple(answer)


def map_from_image(
    vertices: tuple[int, ...], image: tuple[int, ...]
) -> dict[int, int]:
    return dict(zip(vertices, image))


def compose(
    vertices: tuple[int, ...],
    left: tuple[int, ...],
    right: tuple[int, ...],
) -> tuple[int, ...]:
    left_map = map_from_image(vertices, left)
    right_map = map_from_image(vertices, right)
    return tuple(left_map[right_map[c]] for c in vertices)


def weighted_signature_distance(
    out: tuple[tuple[bool, ...], ...],
    vertices: tuple[int, ...],
    left: tuple[int, ...],
    right: tuple[int, ...],
    weights: tuple[int, ...],
    scale: int,
) -> int:
    left_map = map_from_image(vertices, left)
    right_map = map_from_image(vertices, right)
    exterior = tuple(t for t in range(len(out)) if t not in vertices)
    return scale * sum(
        weights[t]
        for t in exterior
        for c in vertices
        if out[left_map[c]][t] != out[right_map[c]][t]
    )


def cycles(
    vertices: tuple[int, ...], rho: dict[int, int]
) -> tuple[tuple[int, ...], ...]:
    unseen = set(vertices)
    answer = []
    while unseen:
        start = min(unseen)
        orbit = []
        c = start
        while c in unseen:
            unseen.remove(c)
            orbit.append(c)
            c = rho[c]
        require(c == start, "orbit did not close")
        answer.append(tuple(orbit))
    return tuple(answer)


def check_case(
    out: tuple[tuple[bool, ...], ...],
    vertices: tuple[int, ...],
    image: tuple[int, ...],
) -> None:
    q = len(out)
    rho = map_from_image(vertices, image)
    sigma = {i: rho[i] if i in rho else i for i in range(q)}
    scale = 3
    weights = tuple(t + 2 for t in range(q))
    sizes = tuple(scale if i in vertices else weights[i] for i in range(q))

    quotient_sum = sum(
        sizes[i] * sizes[j]
        for i, j in combinations(range(q), 2)
        if out[i][j] != out[sigma[i]][sigma[j]]
    )
    literal_expansion = sum(
        1
        for i, j in combinations(range(q), 2)
        for _x in range(sizes[i])
        for _y in range(sizes[j])
        if out[i][j] != out[sigma[i]][sigma[j]]
    )
    require(quotient_sum == literal_expansion, "expanded pair mismatch")

    exterior = tuple(t for t in range(q) if t not in vertices)
    hamming_sum = scale * sum(
        weights[t]
        for t in exterior
        for c in vertices
        if out[c][t] != out[rho[c]][t]
    )
    require(quotient_sum == hamming_sum, "external Hamming mismatch")

    set_sum = 0
    transition_sum = 0
    orbit_constant = True
    for t in exterior:
        support = {c for c in vertices if out[c][t]}
        inverse_support = {c for c in vertices if rho[c] in support}
        set_row = len(support.symmetric_difference(inverse_support))
        intersection_row = 2 * (len(support) - len(support & inverse_support))
        require(set_row == intersection_row, "symmetric-difference mismatch")
        set_sum += scale * weights[t] * set_row

        row_transitions = 0
        for orbit in cycles(vertices, rho):
            bits = tuple(out[c][t] for c in orbit)
            row_transitions += sum(
                bits[k] != bits[(k + 1) % len(bits)]
                for k in range(len(bits))
            )
            orbit_constant &= len(set(bits)) == 1
        transition_sum += scale * weights[t] * row_transitions

        if (
            len(vertices) == 3
            and set(image) == set(vertices)
            and all(rho[rho[rho[c]]] == c for c in vertices)
            and all(rho[c] != c for c in vertices)
        ):
            require(set_row in (0, 2), "triangle transition collapse failed")

    require(hamming_sum == set_sum, "set formula mismatch")
    require(hamming_sum == transition_sum, "orbit transition mismatch")
    require(hamming_sum % (2 * scale) == 0, "parity mismatch")

    extends = all(
        out[i][j] == out[sigma[i]][sigma[j]]
        for i, j in combinations(range(q), 2)
    )
    require((hamming_sum == 0) == extends, "extension kernel mismatch")
    require((hamming_sum == 0) == orbit_constant, "orbit kernel mismatch")


def check_pseudometric(
    out: tuple[tuple[bool, ...], ...],
    vertices: tuple[int, ...],
    autos: tuple[tuple[int, ...], ...],
) -> None:
    q = len(out)
    weights = tuple(t + 2 for t in range(q))
    scale = 3
    for alpha in autos:
        require(
            weighted_signature_distance(
                out, vertices, alpha, alpha, weights, scale
            )
            == 0,
            "pseudometric diagonal mismatch",
        )
        for beta in autos:
            dab = weighted_signature_distance(
                out, vertices, alpha, beta, weights, scale
            )
            dba = weighted_signature_distance(
                out, vertices, beta, alpha, weights, scale
            )
            require(dab == dba, "pseudometric symmetry mismatch")
            for gamma in autos:
                dac = weighted_signature_distance(
                    out, vertices, alpha, gamma, weights, scale
                )
                dcb = weighted_signature_distance(
                    out, vertices, gamma, beta, weights, scale
                )
                require(dab <= dac + dcb, "triangle inequality mismatch")

                alpha_gamma = compose(vertices, alpha, gamma)
                beta_gamma = compose(vertices, beta, gamma)
                translated = weighted_signature_distance(
                    out,
                    vertices,
                    alpha_gamma,
                    beta_gamma,
                    weights,
                    scale,
                )
                require(dab == translated, "right invariance mismatch")


def main() -> None:
    labelled_tournaments = 0
    partial_automorphisms = 0
    for q in range(2, 6):
        edge_count = q * (q - 1) // 2
        for mask in range(1 << edge_count):
            out = tournament(q, mask)
            labelled_tournaments += 1
            for subset_mask in range(1, 1 << q):
                vertices = tuple(
                    i for i in range(q) if (subset_mask >> i) & 1
                )
                autos = automorphisms(out, vertices)
                partial_automorphisms += len(autos)
                for image in autos:
                    check_case(out, vertices, image)
                check_pseudometric(out, vertices, autos)

    require(labelled_tournaments == 1098, "tournament census drift")
    require(partial_automorphisms == 41026, "case census drift")
    print("labelled_tournaments_q2_to_q5=1098")
    print("nonempty_partial_automorphism_cases=41026")
    print("expanded_block_pair_formula=PASS")
    print("hamming_symmetric_difference_transition_forms=PASS")
    print("extension_kernel_and_parity=PASS")
    print("right_invariant_pseudometric=PASS")
    print("status=THM-2195_EXACT_INDEXING_REFEREE_ONLY")


if __name__ == "__main__":
    main()
