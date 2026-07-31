#!/usr/bin/env python3
"""Exact rooted-forest leading-monomial audit for THM-2765.

For a parent-before-child order, the injection Vandermonde contributes
X_j^(j-1), while every earlier edge contributes both mirror factors to the
new edge and hence X_j^2.  The resulting coefficient-grid certificate gives
linear-range labels with distinct absolute edge differences.  All checks use
explicit exceptions; there are no truth-bearing Python assertions.
"""

from itertools import product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def parent_arrays(vertex_count, root_count=1):
    """All ordered rooted forests with the first root_count vertices roots.

    Every later vertex j chooses one parent with smaller index.  For the tree
    census root_count=1, this is the full recursive-tree parent-array atlas.
    """

    if vertex_count < root_count:
        return
    if vertex_count == root_count:
        yield tuple([-1] * root_count)
        return
    choices = [range(j) for j in range(root_count, vertex_count)]
    for tail in product(*choices):
        yield tuple([-1] * root_count) + tail


def target_exponents(vertex_count, root_count=1):
    exponents = [j for j in range(root_count)]
    for edge_rank, child in enumerate(range(root_count, vertex_count), start=1):
        exponents.append(child + 2 * (edge_rank - 1))
    return tuple(exponents)


def add_term(terms, exponent, coefficient):
    if coefficient == 0:
        return
    terms[exponent] = terms.get(exponent, 0) + coefficient
    if terms[exponent] == 0:
        del terms[exponent]


def unit_exponent(vertex_count, vertex, power=1):
    exponent = [0] * vertex_count
    exponent[vertex] = power
    return tuple(exponent)


def vandermonde_factor(vertex_count, first, second):
    return {
        unit_exponent(vertex_count, second): 1,
        unit_exponent(vertex_count, first): -1,
    }


def squared_edge_difference_factor(vertex_count, parents, first, second):
    """Expand (Y_second)^2-(Y_first)^2 exactly."""

    terms = {}
    for child, sign in ((second, 1), (first, -1)):
        parent = parents[child]
        add_term(terms, unit_exponent(vertex_count, child, 2), sign)
        add_term(terms, unit_exponent(vertex_count, parent, 2), sign)
        mixed = [0] * vertex_count
        mixed[child] += 1
        mixed[parent] += 1
        add_term(terms, tuple(mixed), -2 * sign)
    return terms


def coefficient_of_target(vertex_count, parents, root_count=1):
    """Truncated exact expansion of the claimed top monomial coefficient."""

    target = target_exponents(vertex_count, root_count)
    state = {tuple([0] * vertex_count): 1}
    factors = []
    for first in range(vertex_count):
        for second in range(first + 1, vertex_count):
            factors.append(vandermonde_factor(vertex_count, first, second))
    children = list(range(root_count, vertex_count))
    for first_index in range(len(children)):
        for second_index in range(first_index + 1, len(children)):
            factors.append(squared_edge_difference_factor(
                vertex_count,
                parents,
                children[first_index],
                children[second_index],
            ))

    for factor in factors:
        next_state = {}
        for old_exponent, old_coefficient in state.items():
            for increment, factor_coefficient in factor.items():
                new_exponent = tuple(
                    old_exponent[index] + increment[index]
                    for index in range(vertex_count)
                )
                if any(new_exponent[index] > target[index]
                       for index in range(vertex_count)):
                    continue
                add_term(
                    next_state,
                    new_exponent,
                    old_coefficient * factor_coefficient,
                )
        state = next_state
    return state.get(target, 0), len(state)


def find_labeling(parents, root_count=1):
    """Find one point in the exact coefficient grids by backtracking."""

    vertex_count = len(parents)
    exponents = target_exponents(vertex_count, root_count)
    labels = [None] * vertex_count
    used_labels = set()
    used_differences = set()
    nodes = 0

    def search(vertex):
        nonlocal nodes
        if vertex == vertex_count:
            return True
        for label in range(exponents[vertex] + 1):
            nodes += 1
            if label in used_labels:
                continue
            difference = None
            if vertex >= root_count:
                difference = abs(label - labels[parents[vertex]])
                if difference == 0 or difference in used_differences:
                    continue
            labels[vertex] = label
            used_labels.add(label)
            if difference is not None:
                used_differences.add(difference)
            if search(vertex + 1):
                return True
            if difference is not None:
                used_differences.remove(difference)
            used_labels.remove(label)
            labels[vertex] = None
        return False

    require(search(0), "a rooted coefficient grid had no valid labeling")
    return tuple(labels), nodes


def edge_differences(parents, labels, root_count=1):
    return tuple(
        abs(labels[child] - labels[parents[child]])
        for child in range(root_count, len(parents))
    )


def main():
    coefficient_arrays = 0
    max_truncated_states = 0
    for vertex_count in range(2, 7):
        expected = tuple([0] + [3 * j - 2 for j in range(1, vertex_count)])
        require(target_exponents(vertex_count) == expected,
                "the tree exponent ladder changed")
        expected_degree = (
            vertex_count * (vertex_count - 1) // 2
            + (vertex_count - 1) * (vertex_count - 2)
        )
        require(sum(expected) == expected_degree,
                "the leading monomial stopped having full total degree")
        for parents in parent_arrays(vertex_count):
            coefficient, state_count = coefficient_of_target(
                vertex_count, parents
            )
            coefficient_arrays += 1
            max_truncated_states = max(max_truncated_states, state_count)
            require(coefficient == 1,
                    "the rooted leading coefficient stopped being one")

    recursive_trees = 0
    search_nodes = 0
    largest_search = 0
    largest_label_seen = 0
    for vertex_count in range(2, 9):
        for parents in parent_arrays(vertex_count):
            labels, nodes = find_labeling(parents)
            differences = edge_differences(parents, labels)
            require(len(set(labels)) == vertex_count,
                    "a certificate labeling repeated a vertex value")
            require(len(set(differences)) == vertex_count - 1
                    and min(differences) > 0,
                    "a certificate labeling repeated an absolute edge value")
            require(max(labels) <= 3 * vertex_count - 5,
                    "a certificate labeling escaped the linear range")
            recursive_trees += 1
            search_nodes += nodes
            largest_search = max(largest_search, nodes)
            largest_label_seen = max(largest_label_seen, max(labels))

    # Forest exponent control, including isolated roots.  The predicted maximum
    # is n+2m-3 for m>0 and n-1 for the edgeless case.
    forest_controls = 0
    for vertex_count in range(2, 9):
        for root_count in range(1, vertex_count + 1):
            edge_count = vertex_count - root_count
            exponents = target_exponents(vertex_count, root_count)
            expected_max = (
                vertex_count - 1 if edge_count == 0
                else vertex_count + 2 * edge_count - 3
            )
            require(max(exponents) == expected_max,
                    "the rooted-forest range formula changed")
            forest_controls += 1

    # Both parts of the product are necessary.
    mirror_hostile_parents = (-1, 0, 0)
    mirror_hostile_labels = (1, 0, 2)
    oriented = tuple(
        mirror_hostile_labels[child]
        - mirror_hostile_labels[mirror_hostile_parents[child]]
        for child in range(1, 3)
    )
    require(oriented == (-1, 1)
            and len(set(oriented)) == 2
            and len(set(abs(value) for value in oriented)) == 1,
            "the oriented-versus-absolute hostile changed")

    injectivity_hostile_parents = (-1, 0, 1, 2)
    injectivity_hostile_labels = (0, 1, 3, 0)
    injectivity_differences = edge_differences(
        injectivity_hostile_parents,
        injectivity_hostile_labels,
    )
    require(len(set(injectivity_hostile_labels)) < 4
            and injectivity_differences == (1, 2, 3),
            "the edge-distinct-but-noninjective hostile changed")

    print("ROOTED FOREST NULLSTELLENSATZ RANGE AUDIT")
    print("tree_top_monomial=product_j>=2 X_j^(3j-5) coefficient=1")
    print(f"exact_coefficient_parent_arrays_n2_to_n6={coefficient_arrays}")
    print(f"max_truncated_coefficient_states={max_truncated_states}")
    print(f"recursive_tree_labeling_controls_n2_to_n8={recursive_trees}")
    print(
        f"backtracking_nodes_total={search_nodes} "
        f"max_per_tree={largest_search} max_label_seen={largest_label_seen}"
    )
    print(f"forest_exponent_controls={forest_controls}")
    print("tree_range=0..3n-5 forest_range=0..n+2m-3_for_m>0")
    print("mirror_hostile=P3_labels(1,0,2)_oriented=(-1,1)_absolute_collision")
    print("injectivity_hostile=P4_labels(0,1,3,0)_edge_differences=(1,2,3)")
    print("SCOPE: linear-range distinct-edge labeling; not graceful-tree closure")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
