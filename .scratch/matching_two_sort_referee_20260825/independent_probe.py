"""Independent semantic and sort-flow probe for the THM-4090 candidate.

This is deliberately independent of the bit-mask implementation in
arxiv_260813306_audit_20260825/two_sort_boundary_probe.py: denotations are
computed as Python frozensets, and the rule check is stated as a reachability
test on the exact two-sort signature.
"""

from itertools import chain, combinations, product


def powerset(carrier: frozenset[int]) -> tuple[frozenset[int], ...]:
    points = tuple(carrier)
    return tuple(
        frozenset(choice)
        for size in range(len(points) + 1)
        for choice in combinations(points, size)
    )


def gamma_denotation(
    b: frozenset[int],
    a: frozenset[int],
    f: dict[int, frozenset[int]],
) -> frozenset[int]:
    """Denotation of forall x:b forall y:b. f(x and y), of sort a."""

    values = []
    for x, y in product(b, repeat=2):
        x_value = frozenset({x})
        y_value = frozenset({y})
        meet = x_value & y_value
        # Pointwise extension of a unary symbol; union(empty family) = empty.
        values.append(frozenset(chain.from_iterable(f[z] for z in meet)))
    return frozenset.intersection(*values)


def phi_denotation(b: frozenset[int]) -> frozenset[int]:
    """Denotation of forall x:b. x, of sort b."""

    return frozenset.intersection(*(frozenset({x}) for x in b))


def reflexive_transitive_closure(
    sorts: frozenset[str], edges: frozenset[tuple[str, str]]
) -> frozenset[tuple[str, str]]:
    closure = set(edges) | {(sort, sort) for sort in sorts}
    changed = True
    while changed:
        changed = False
        for source, middle in tuple(closure):
            for middle_again, target in tuple(closure):
                if middle == middle_again and (source, target) not in closure:
                    closure.add((source, target))
                    changed = True
    return frozenset(closure)


def main() -> None:
    total_interpretations = 0
    gamma_models = 0
    entailment_failures = []
    gamma_counts: dict[tuple[int, int], int] = {}

    for b_size in range(1, 5):
        b = frozenset(range(b_size))
        for a_size in range(1, 4):
            a = frozenset(range(a_size))
            subsets_a = powerset(a)
            count = 0
            for images in product(subsets_a, repeat=b_size):
                total_interpretations += 1
                f = dict(zip(sorted(b), images, strict=True))
                gamma_total = gamma_denotation(b, a, f) == a
                phi_total = phi_denotation(b) == b
                if gamma_total:
                    gamma_models += 1
                    count += 1
                    if not phi_total:
                        entailment_failures.append((b_size, a_size, images))

                # Closed-form semantic identities used in the proof.
                assert gamma_total == (b_size == 1 and images[0] == a)
                assert phi_total == (b_size == 1)
            gamma_counts[(b_size, a_size)] = count

    # f:b->a is the only positive-arity symbol.  Primitive framing therefore
    # contributes b->a; all other premise-bearing primitive rules preserve sort.
    sorts = frozenset({"b", "a"})
    symbol_edges = frozenset({("b", "a")})
    feeds = reflexive_transitive_closure(sorts, symbol_edges)
    assert ("a", "b") not in feeds
    assert feeds == frozenset({("a", "a"), ("b", "b"), ("b", "a")})

    # Positive and hostile controls.
    assert gamma_denotation(frozenset({0}), frozenset({0, 1}), {0: frozenset({0, 1})}) == frozenset({0, 1})
    assert phi_denotation(frozenset({0})) == frozenset({0})
    assert phi_denotation(frozenset({0, 1})) == frozenset()

    assert total_interpretations == 5050
    assert gamma_models == 3
    assert not entailment_failures
    assert all(count == (1 if b_size == 1 else 0) for (b_size, _), count in gamma_counts.items())

    print("independent two-sort referee probe: PASS")
    print(f"interpretations checked: {total_interpretations}")
    print(f"Gamma models: {gamma_models}")
    print("sort-flow closure:", sorted(feeds))
    print("a does not feed b; Gamma-dependent sort-b lines are impossible")


if __name__ == "__main__":
    main()
