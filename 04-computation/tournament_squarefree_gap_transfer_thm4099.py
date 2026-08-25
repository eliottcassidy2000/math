#!/usr/bin/env python3
"""Exact primary referee for THM-4099.

The theorem proves a squarefree gap polynomial for arbitrary finite inserted
sets.  This companion specializes to two inserted vertices and exhausts every
labeled tournament through final order six.  It checks the closed local gap
coefficients against direct block enumeration, contracts the four-state
transfer matrix, and compares all four coefficients with an independent
subset dynamic program for Hamiltonian-path counts.

The same universe checks the sequential derivative/orphan formula, the full
first-deletion fibre decomposition, the degree/defect filtration, and the
minimal mixed-coefficient hostile.  Only the Python standard library is used.
Every executable gate goes through ``require`` so optimization cannot disable
the audit.
"""

from hashlib import sha256
from itertools import permutations
import json


def require(condition, label):
    if not condition:
        raise RuntimeError("CHECK FAILED: " + label)


PERMUTATIONS = {
    n: tuple(permutations(range(n)))
    for n in range(1, 7)
}


def pair_bit(n, i, j):
    require(0 <= i < j < n, "pair_bit domain")
    return i * (2 * n - i - 1) // 2 + (j - i - 1)


def edge(mask, n, i, j):
    """Return whether the labeled tournament encoded by mask has i -> j."""
    require(0 <= i < n and 0 <= j < n and i != j, "edge domain")
    if i > j:
        return not edge(mask, n, j, i)
    return bool((mask >> pair_bit(n, i, j)) & 1)


def mask_from_arc_set(n, arcs):
    arcs = set(arcs)
    mask = 0
    for i in range(n):
        for j in range(i + 1, n):
            forward = (i, j) in arcs
            backward = (j, i) in arcs
            require(forward != backward, "exactly one orientation per pair")
            if forward:
                mask |= 1 << pair_bit(n, i, j)
    return mask


def is_directed_path(mask, n, word):
    return all(edge(mask, n, word[i], word[i + 1])
               for i in range(len(word) - 1))


def hamilton_paths_prefix(mask, n, prefix_size):
    require(1 <= prefix_size <= n, "Hamilton prefix size")
    return tuple(
        word
        for word in PERMUTATIONS[prefix_size]
        if is_directed_path(mask, n, word)
    )


def subset_hamilton_counts(mask, n):
    """Independent DP: counts Hamilton paths on every induced vertex subset."""
    dp = [[0] * n for _ in range(1 << n)]
    for vertex in range(n):
        dp[1 << vertex][vertex] = 1
    for subset in range(1, 1 << n):
        if subset & (subset - 1) == 0:
            continue
        for last in range(n):
            if not (subset >> last) & 1:
                continue
            previous = subset ^ (1 << last)
            total = 0
            for penultimate in range(n):
                if ((previous >> penultimate) & 1
                        and edge(mask, n, penultimate, last)):
                    total += dp[previous][penultimate]
            dp[subset][last] = total
    counts = [0] * (1 << n)
    counts[0] = 1
    for subset in range(1, 1 << n):
        counts[subset] = sum(dp[subset])
    return tuple(counts)


def multiply_squarefree(left, right):
    """Multiply in Z[U,V]/(U^2,V^2), basis (1,U,V,UV)."""
    a, b, c, d = left
    e, f, g, h = right
    return (
        a * e,
        a * f + b * e,
        a * g + c * e,
        a * h + b * g + c * f + d * e,
    )


def matrix_step(gap, state):
    """Apply the THM-4099 four-state gap matrix to a coefficient column."""
    a, b, c, d = gap
    s0, su, sv, suv = state
    return (
        a * s0,
        b * s0 + a * su,
        c * s0 + a * sv,
        d * s0 + c * su + b * sv + a * suv,
    )


def block_is_legal(mask, n, left, block, right):
    word = (() if left is None else (left,)) + block
    word += () if right is None else (right,)
    return is_directed_path(mask, n, word)


def local_coefficients_bruteforce(mask, n, left, right, u, v):
    """Enumerate the four local inserted subsets independently."""
    return (
        int(block_is_legal(mask, n, left, (), right)),
        int(block_is_legal(mask, n, left, (u,), right)),
        int(block_is_legal(mask, n, left, (v,), right)),
        int(block_is_legal(mask, n, left, (u, v), right))
        + int(block_is_legal(mask, n, left, (v, u), right)),
    )


def local_coefficients_closed(mask, n, left, right, u, v):
    """Closed coefficients a+bU+cV+dUV from the theorem."""
    if left is None:
        require(right is not None, "nonempty base left boundary")
        return (
            1,
            int(edge(mask, n, u, right)),
            int(edge(mask, n, v, right)),
            int(edge(mask, n, u, v) and edge(mask, n, v, right))
            + int(edge(mask, n, v, u) and edge(mask, n, u, right)),
        )
    if right is None:
        return (
            1,
            int(edge(mask, n, left, u)),
            int(edge(mask, n, left, v)),
            int(edge(mask, n, left, u) and edge(mask, n, u, v))
            + int(edge(mask, n, left, v) and edge(mask, n, v, u)),
        )
    return (
        int(edge(mask, n, left, right)),
        int(edge(mask, n, left, u) and edge(mask, n, u, right)),
        int(edge(mask, n, left, v) and edge(mask, n, v, right)),
        int(edge(mask, n, left, u)
            and edge(mask, n, u, v)
            and edge(mask, n, v, right))
        + int(edge(mask, n, left, v)
              and edge(mask, n, v, u)
              and edge(mask, n, u, right)),
    )


def gaps(word):
    require(word, "nonempty base word")
    return (
        ((None, word[0]),)
        + tuple((word[i], word[i + 1]) for i in range(len(word) - 1))
        + ((word[-1], None),)
    )


def word_gap_polynomial(mask, n, word, local_table):
    state = (1, 0, 0, 0)
    for left, right in gaps(word):
        closed = local_table[(left, right)]
        via_matrix = matrix_step(closed, state)
        via_ring = multiply_squarefree(state, closed)
        require(via_matrix == via_ring, "four-state matrix orientation")
        state = via_matrix
    defect = sum(
        not edge(mask, n, word[i], word[i + 1])
        for i in range(len(word) - 1)
    )
    degrees = (0, 1, 1, 2)
    require(
        all(coefficient == 0 or degree >= defect
            for coefficient, degree in zip(state, degrees)),
        "degree is at least the number of bad old adjacencies",
    )
    return state, defect


def transfer_polynomial(mask, n, base_size, layer_examples):
    require(n == base_size + 2, "two inserted vertices")
    u, v = base_size, base_size + 1
    local_table = {}
    gap_keys = (
        tuple((None, x) for x in range(base_size))
        + tuple((x, None) for x in range(base_size))
        + tuple(
            (x, y)
            for x in range(base_size)
            for y in range(base_size)
            if x != y
        )
    )
    for left, right in gap_keys:
        closed = local_coefficients_closed(mask, n, left, right, u, v)
        brute = local_coefficients_bruteforce(mask, n, left, right, u, v)
        require(closed == brute, "closed/local-block coefficient agreement")
        local_table[(left, right)] = closed
    total = [0, 0, 0, 0]
    local_checks = len(local_table)
    defect_histogram = {}
    exact_degree_mass = {}
    degrees = (0, 1, 1, 2)
    for word in PERMUTATIONS[base_size]:
        polynomial, defect = word_gap_polynomial(mask, n, word, local_table)
        defect_histogram[defect] = defect_histogram.get(defect, 0) + 1
        exact_mass = sum(
            coefficient
            for coefficient, degree in zip(polynomial, degrees)
            if degree == defect
        )
        exact_degree_mass[defect] = exact_degree_mass.get(defect, 0) + exact_mass
        if defect in (0, 1, 2) and exact_mass and defect not in layer_examples:
            layer_examples[defect] = {
                "final_n": n,
                "mask": mask,
                "base_word": list(word),
                "defect": defect,
                "polynomial": list(polynomial),
            }
        for index, coefficient in enumerate(polynomial):
            total[index] += coefficient
    return (
        tuple(total),
        local_checks,
        tuple(sorted(defect_histogram.items())),
        tuple(sorted(exact_degree_mass.items())),
    )


def signature_tau(mask, n, word, v):
    signature = tuple(int(edge(mask, n, v, x)) for x in word)
    tau = sum(signature[i] == 1 and signature[i + 1] == 0
              for i in range(len(signature) - 1))
    return signature, tau


def insert_at(word, vertex, position):
    return word[:position] + (vertex,) + word[position:]


def first_deletion_fibres(mask, n, u):
    """Full THM-4094 decomposition for A=T+u, with A a prefix of W."""
    require(u == n - 2, "first inserted vertex label")
    old_paths = hamilton_paths_prefix(mask, n, u)
    new_paths = hamilton_paths_prefix(mask, n, u + 1)
    fibres = {path: [] for path in old_paths}
    orphans = []
    for path in new_paths:
        reduced = tuple(x for x in path if x != u)
        if reduced in fibres:
            fibres[reduced].append(path)
        else:
            orphans.append(path)
    require(all(fibres[path] for path in old_paths), "first deletion left-total")
    require(
        sum(len(fibres[path]) for path in old_paths) + len(orphans)
        == len(new_paths),
        "full first deletion partition",
    )
    return old_paths, new_paths, fibres, tuple(orphans)


def sequential_components(mask, n, expected_h):
    """Check H(A+v)=L(1)+L'(1)+F and its exact witness partition."""
    require(n >= 3, "successive insertion base is nonempty")
    u, v = n - 2, n - 1
    _, paths_a, fibres, first_orphans = first_deletion_fibres(mask, n, u)

    histogram = {}
    nonorphan_targets = set()
    legal_position_total = 0
    for path in paths_a:
        _, tau = signature_tau(mask, n, path, v)
        histogram[tau] = histogram.get(tau, 0) + 1
        legal_positions = []
        for position in range(len(path) + 1):
            target = insert_at(path, v, position)
            if is_directed_path(mask, n, target):
                legal_positions.append(position)
                nonorphan_targets.add(target)
        require(len(legal_positions) == 1 + tau, "binary transition balance")
        legal_position_total += len(legal_positions)
    require(len(nonorphan_targets) == legal_position_total,
            "legal insertions have unique deletion source")

    fibre_histogram = {}
    for path_list in tuple(fibres.values()) + (first_orphans,):
        for path in path_list:
            _, tau = signature_tau(mask, n, path, v)
            fibre_histogram[tau] = fibre_histogram.get(tau, 0) + 1
    require(histogram == fibre_histogram,
            "full first-deletion fibres distribute L exactly")

    repaired_targets = set()
    repaired_words = 0
    for word in PERMUTATIONS[n - 1]:
        bad = tuple(
            i
            for i in range(n - 2)
            if not edge(mask, n, word[i], word[i + 1])
        )
        if len(bad) != 1:
            continue
        cut = bad[0]
        if edge(mask, n, word[cut], v) and edge(mask, n, v, word[cut + 1]):
            target = insert_at(word, v, cut + 1)
            require(is_directed_path(mask, n, target), "repaired cut is Hamiltonian")
            repaired_targets.add(target)
            repaired_words += 1
    require(len(repaired_targets) == repaired_words,
            "failed-cut repair has unique deletion word")
    require(nonorphan_targets.isdisjoint(repaired_targets),
            "nonorphan and repaired-orphan witnesses are disjoint")

    l_at_one = sum(histogram.values())
    l_derivative_at_one = sum(tau * count for tau, count in histogram.items())
    formula = l_at_one + l_derivative_at_one + repaired_words
    require(formula == expected_h, "sequential derivative formula")
    require(len(nonorphan_targets | repaired_targets) == expected_h,
            "sequential witnesses exhaust independent DP count")
    return {
        "H_A": len(paths_a),
        "L_histogram": tuple(sorted(histogram.items())),
        "L_prime_1": l_derivative_at_one,
        "F": repaired_words,
        "H_final": expected_h,
        "first_orphans": len(first_orphans),
        "legal_targets": len(nonorphan_targets),
    }


def first_insertion_profile(mask, n):
    require(n == 3, "minimal hostile order")
    u = 1
    old, new, fibres, orphans = first_deletion_fibres(mask, n, u)
    return (
        len(old),
        tuple(len(fibres[path]) for path in old),
        len(orphans),
        len(new),
    )


def main():
    ledger = {
        "theorem": "THM-4099",
        "ring_basis": ["1", "U", "V", "UV"],
        "independent_counter": "endpoint subset DP",
    }

    print("SQUAREFREE GAP TRANSFER EXHAUSTIVE AUDIT")
    squarefree_rows = []
    sequential_rows = []
    layer_examples = {}
    for base_size in range(1, 5):
        n = base_size + 2
        tournament_count = 1 << (n * (n - 1) // 2)
        coefficient_checks = 0
        local_checks = 0
        defect_word_totals = {}
        exact_degree_mass = {}
        sequential_totals = {
            "intermediate_paths": 0,
            "legal_targets": 0,
            "failed_repairs": 0,
            "first_orphans": 0,
        }
        for mask in range(tournament_count):
            dp_counts = subset_hamilton_counts(mask, n)
            polynomial, checks, defect_histogram, degree_mass = transfer_polynomial(
                mask, n, base_size, layer_examples
            )
            local_checks += checks
            for defect, count in defect_histogram:
                defect_word_totals[defect] = defect_word_totals.get(defect, 0) + count
            for defect, mass in degree_mass:
                exact_degree_mass[defect] = exact_degree_mass.get(defect, 0) + mass

            base = (1 << base_size) - 1
            u_bit = 1 << base_size
            v_bit = 1 << (base_size + 1)
            expected = (
                dp_counts[base],
                dp_counts[base | u_bit],
                dp_counts[base | v_bit],
                dp_counts[base | u_bit | v_bit],
            )
            require(polynomial == expected,
                    "all squarefree coefficients equal independent DP")
            coefficient_checks += 4

            components = sequential_components(mask, n, expected[3])
            sequential_totals["intermediate_paths"] += components["H_A"]
            sequential_totals["legal_targets"] += components["legal_targets"]
            sequential_totals["failed_repairs"] += components["F"]
            sequential_totals["first_orphans"] += components["first_orphans"]

        row = {
            "base_n": base_size,
            "final_n": n,
            "tournaments": tournament_count,
            "coefficient_checks": coefficient_checks,
            "local_checks": local_checks,
            "defect_words": tuple(sorted(defect_word_totals.items())),
            "exact_degree_mass": tuple(sorted(exact_degree_mass.items())),
            "failures": 0,
        }
        squarefree_rows.append(row)
        sequential_row = {
            "final_n": n,
            "tournaments": tournament_count,
            **sequential_totals,
            "failures": 0,
        }
        sequential_rows.append(sequential_row)
        print(
            "base_n=", base_size,
            "final_n=", n,
            "tournaments=", tournament_count,
            "coefficient_checks=", coefficient_checks,
            "local_checks=", local_checks,
            "defect_words=", row["defect_words"],
            "failures=", 0,
        )

    require(
        [row["tournaments"] for row in squarefree_rows] == [8, 64, 1024, 32768],
        "frozen exhaustive tournament universes",
    )
    require(set(layer_examples) == {0, 1, 2}, "all live defect layers witnessed")
    ledger["squarefree_rows"] = squarefree_rows
    ledger["layer_examples"] = layer_examples

    print("\nSEQUENTIAL DERIVATIVE / FULL-FIBRE AUDIT")
    for row in sequential_rows:
        print(
            "final_n=", row["final_n"],
            "tournaments=", row["tournaments"],
            "intermediate_paths=", row["intermediate_paths"],
            "legal_targets=", row["legal_targets"],
            "failed_repairs=", row["failed_repairs"],
            "first_orphans=", row["first_orphans"],
            "failures=", row["failures"],
        )
    ledger["sequential_rows"] = sequential_rows

    print("\nLIVE DEFECT-LAYER CONTROLS")
    for defect in (0, 1, 2):
        print("defect=", defect, "example=", layer_examples[defect])

    transitive = mask_from_arc_set(3, {(0, 1), (0, 2), (2, 1)})
    cyclic = mask_from_arc_set(3, {(1, 0), (0, 2), (2, 1)})
    hostile_rows = []
    print("\nMINIMAL MIXED-COEFFICIENT HOSTILE")
    for name, mask in (("transitive", transitive), ("C3", cyclic)):
        counts = subset_hamilton_counts(mask, 3)
        polynomial, _, _, _ = transfer_polynomial(mask, 3, 1, {})
        components = sequential_components(mask, 3, counts[7])
        profile = first_insertion_profile(mask, 3)
        row = {
            "name": name,
            "mask": mask,
            "polynomial": polynomial,
            "first_profile": profile,
            "sequential": components,
        }
        hostile_rows.append(row)
        print(
            "name=", name,
            "mask=", mask,
            "Z=", polynomial,
            "first_profile=", profile,
            "L_histogram=", components["L_histogram"],
            "F=", components["F"],
            "H_final=", components["H_final"],
        )
    require(hostile_rows[0]["first_profile"] == hostile_rows[1]["first_profile"],
            "hostile first-insertion count profiles agree")
    require(hostile_rows[0]["polynomial"][:3]
            == hostile_rows[1]["polynomial"][:3] == (1, 1, 1),
            "hostile proper-face coefficients agree")
    require(hostile_rows[0]["polynomial"][3] == 1,
            "transitive mixed coefficient")
    require(hostile_rows[1]["polynomial"][3] == 3,
            "cyclic mixed coefficient")
    require(edge(transitive, 3, 0, 2) and edge(cyclic, 3, 0, 2),
            "same future arc 0 to v")
    require(edge(transitive, 3, 2, 1) and edge(cyclic, 3, 2, 1),
            "same future arc v to u")
    order_two_h = sorted({
        subset_hamilton_counts(mask, 2)[3]
        for mask in range(2)
    })
    require(order_two_h == [1], "no smaller mixed hostile")
    print("order_two_H_support=", order_two_h, "so_final_order_3_is_minimal=True")
    ledger["minimal_hostile"] = hostile_rows
    ledger["order_two_H_support"] = order_two_h

    semantic_payload = json.dumps(
        ledger,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    semantic_digest = sha256(semantic_payload).hexdigest()
    print("\nsemantic_sha256=", semantic_digest)
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()
