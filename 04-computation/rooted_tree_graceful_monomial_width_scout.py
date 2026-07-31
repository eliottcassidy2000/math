#!/usr/bin/env python3
"""Exact bounded-monomial width scout for the rooted graceful polynomial.

This is a finite scout, not a theorem dependency.  It asks for the smallest R
such that the vertex Vandermonde times the squared edge-difference
discriminant has a nonzero full-degree monomial with every exponent at most R.
The coefficient-grid theorem would then use a common (R+1)-point label set.
"""

from collections import Counter

from rooted_forest_nullstellensatz_range_thm2765 import (
    add_term,
    parent_arrays,
    squared_edge_difference_factor,
    vandermonde_factor,
)


def factors(vertex_count, parents):
    out = []
    for first in range(vertex_count):
        for second in range(first + 1, vertex_count):
            out.append(vandermonde_factor(vertex_count, first, second))
    for first in range(1, vertex_count):
        for second in range(first + 1, vertex_count):
            out.append(squared_edge_difference_factor(
                vertex_count, parents, first, second
            ))
    return tuple(out)


def bounded_state(vertex_count, parents, radius):
    state = {tuple([0] * vertex_count): 1}
    peak = 1
    for factor in factors(vertex_count, parents):
        next_state = {}
        for old_exponent, old_coefficient in state.items():
            for increment, factor_coefficient in factor.items():
                new_exponent = tuple(
                    old_exponent[index] + increment[index]
                    for index in range(vertex_count)
                )
                if max(new_exponent) > radius:
                    continue
                add_term(
                    next_state,
                    new_exponent,
                    old_coefficient * factor_coefficient,
                )
        state = next_state
        peak = max(peak, len(state))
        if not state:
            break
    return state, peak


def unrooted_tree_code(parents):
    vertex_count = len(parents)
    adjacency = [set() for _ in range(vertex_count)]
    for child in range(1, vertex_count):
        parent = parents[child]
        adjacency[child].add(parent)
        adjacency[parent].add(child)

    remaining = set(range(vertex_count))
    while len(remaining) > 2:
        leaves = {
            vertex for vertex in remaining
            if len(adjacency[vertex] & remaining) <= 1
        }
        remaining -= leaves

    def rooted_code(vertex, parent):
        children = sorted(
            rooted_code(neighbor, vertex)
            for neighbor in adjacency[vertex]
            if neighbor != parent
        )
        return "(" + "".join(children) + ")"

    centers = sorted(remaining)
    if len(centers) == 1:
        return rooted_code(centers[0], -1)
    halves = sorted((
        rooted_code(centers[0], centers[1]),
        rooted_code(centers[1], centers[0]),
    ))
    return "[" + "".join(halves) + "]"


def main():
    print("ROOTED TREE GRACEFUL MONOMIAL WIDTH SCOUT")
    grand_peak = 0
    total_trees = 0
    for vertex_count in range(2, 7):
        total_degree = (
            vertex_count * (vertex_count - 1) // 2
            + (vertex_count - 1) * (vertex_count - 2)
        )
        lower = (total_degree + vertex_count - 1) // vertex_count
        representatives = {}
        recursive_count = 0
        for parents in parent_arrays(vertex_count):
            recursive_count += 1
            representatives.setdefault(unrooted_tree_code(parents), parents)
        histogram = Counter()
        details = []
        for parents in representatives.values():
            total_trees += 1
            first_nonzero = None
            for radius in range(lower, 3 * vertex_count - 4):
                state, peak = bounded_state(vertex_count, parents, radius)
                grand_peak = max(grand_peak, peak)
                if state:
                    first_nonzero = radius
                    break
            if first_nonzero is None:
                raise RuntimeError("the certified 3n-5 monomial disappeared")
            histogram[first_nonzero] += 1
            degrees = [0] * vertex_count
            for child in range(1, vertex_count):
                degrees[child] += 1
                degrees[parents[child]] += 1
            signature = "".join(str(value) for value in sorted(degrees, reverse=True))
            details.append((signature, first_nonzero))
        encoded = ",".join(
            f"{radius}:{histogram[radius]}" for radius in sorted(histogram)
        )
        print(
            f"n={vertex_count} degree={total_degree} lower={lower} "
            f"recursive={recursive_count} shapes={len(representatives)} "
            f"width_hist={encoded}"
        )
        print("shape_widths=" + ",".join(
            f"{signature}:{radius}" for signature, radius in sorted(details)
        ))
    print(f"unrooted_shape_tests={total_trees} peak_states={grand_peak}")
    print("SCOPE: finite coefficient scout only; no improved uniform bound proved")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
