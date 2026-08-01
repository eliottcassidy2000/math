#!/usr/bin/env python3
"""Exact verifier for THM-3037's cusp-braid S4 lift dichotomy."""

from collections import Counter
from itertools import permutations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def compose(p, q):
    """Return the permutation p after q."""
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p):
    ans = [None] * len(p)
    for i, j in enumerate(p):
        ans[j] = i
    return tuple(ans)


def power(p, exponent):
    ans = tuple(range(len(p)))
    for _ in range(exponent):
        ans = compose(p, ans)
    return ans


def order(p):
    identity = tuple(range(len(p)))
    ans = identity
    for exponent in range(1, 1 + 4 * len(p)):
        ans = compose(p, ans)
        if ans == identity:
            return exponent
    raise RuntimeError((p, "order bound"))


def cycle_type(p):
    seen = set()
    sizes = []
    for i in range(len(p)):
        if i in seen:
            continue
        j = i
        size = 0
        while j not in seen:
            seen.add(j)
            size += 1
            j = p[j]
        sizes.append(size)
    return tuple(sorted(sizes, reverse=True))


def fixed_points(p):
    return frozenset(i for i, j in enumerate(p) if i == j)


def generated_group(generators):
    identity = tuple(range(len(generators[0])))
    moves = tuple(generators) + tuple(inverse(g) for g in generators)
    group = {identity}
    frontier = [identity]
    while frontier:
        old = frontier.pop()
        for move in moves:
            new = compose(move, old)
            if new not in group:
                group.add(new)
                frontier.append(new)
    return frozenset(group)


def conjugate(v, p):
    return compose(compose(v, p), inverse(v))


S4 = tuple(permutations(range(4)))
S3 = tuple(permutations(range(3)))
ID4 = tuple(range(4))
ID3 = tuple(range(3))

MATCHINGS = (
    frozenset((frozenset((0, 1)), frozenset((2, 3)))),
    frozenset((frozenset((0, 2)), frozenset((1, 3)))),
    frozenset((frozenset((0, 3)), frozenset((1, 2)))),
)
MATCHING_INDEX = {matching: i for i, matching in enumerate(MATCHINGS)}


def matching_image(p):
    image = []
    for matching in MATCHINGS:
        moved = frozenset(
            frozenset((p[min(edge)], p[max(edge)])) for edge in matching
        )
        image.append(MATCHING_INDEX[moved])
    return tuple(image)


# Verify the quotient itself, rather than assuming the standard S4/V4 model.
require({matching_image(p) for p in S4} == set(S3), "matching quotient image")
for p in S4:
    for q in S4:
        require(
            matching_image(compose(p, q))
            == compose(matching_image(p), matching_image(q)),
            (p, q, "quotient homomorphism"),
        )
KERNEL = tuple(p for p in S4 if matching_image(p) == ID3)
require(len(KERNEL) == 4, "V4 kernel size")
require(Counter(cycle_type(p) for p in KERNEL) == Counter({(1, 1, 1, 1): 1, (2, 2): 3}), "V4 kernel types")

TRANSPOSITIONS = tuple(p for p in S3 if cycle_type(p) == (2, 1))
require(len(TRANSPOSITIONS) == 3, "S3 transpositions")


def lift_records(tau_1, tau_2):
    require(tau_1 != tau_2, "distinct meridians")
    require(
        compose(compose(tau_1, tau_2), tau_1)
        == compose(compose(tau_2, tau_1), tau_2),
        "downstairs braid",
    )
    pre_1 = tuple(p for p in S4 if matching_image(p) == tau_1)
    pre_2 = tuple(p for p in S4 if matching_image(p) == tau_2)
    require(len(pre_1) == len(pre_2) == 4, "quotient fibre size")
    records = []
    for x in pre_1:
        for y in pre_2:
            if compose(compose(x, y), x) != compose(compose(y, x), y):
                continue
            modular_a = compose(compose(x, y), x)
            modular_b = compose(x, y)
            parabolic = compose(modular_a, modular_b)
            group = generated_group((x, y))
            common = fixed_points(x) & fixed_points(y)
            require(power(modular_a, 2) == ID4, (x, y, "A squared"))
            require(power(modular_b, 3) == ID4, (x, y, "B cubed"))
            require(power(compose(x, y), 3) == ID4, (x, y, "braid centre"))
            require(parabolic == inverse(x), (x, y, "AB=X inverse"))
            require(order(parabolic) == order(x), (x, y, "cusp width"))
            records.append(
                (
                    cycle_type(x),
                    cycle_type(y),
                    len(group),
                    len(common),
                    order(parabolic),
                    4 - len(cycle_type(x)),
                    3 - len(cycle_type(tau_1)),
                    x,
                    y,
                    tuple(sorted(common)),
                )
            )
    require(len(records) == 8, "eight braided lifts")
    summary = Counter(record[:7] for record in records)
    expected = Counter(
        {
            ((2, 1, 1), (2, 1, 1), 6, 1, 2, 1, 1): 4,
            ((4,), (4,), 24, 0, 4, 3, 1): 4,
        }
    )
    require(summary == expected, (summary, expected))
    require(all(record[0] == record[1] for record in records), "no mixed branch")

    # The branch bit is a literal V4-valued square defect.  Split lifts have
    # zero square; full lifts have two distinct nonzero squares which generate
    # the whole matching kernel.
    for record in records:
        x, y = record[7], record[8]
        x_square, y_square = power(x, 2), power(y, 2)
        if record[0] == (2, 1, 1):
            require(x_square == y_square == ID4, (x, y, "split square defect"))
        else:
            require(x_square in KERNEL and y_square in KERNEL, (x, y, "V4 squares"))
            require(x_square != ID4 and y_square != ID4, (x, y, "nonzero squares"))
            require(x_square != y_square, (x, y, "distinct squares"))
            require(generated_group((x_square, y_square)) == frozenset(KERNEL), (x, y, "squares generate V4"))

    # Simultaneous conjugation by the kernel is free and transitive on each branch.
    pairs = frozenset((record[7], record[8]) for record in records)
    unseen = set(pairs)
    orbits = []
    while unseen:
        pair = next(iter(unseen))
        orbit = frozenset(
            (conjugate(v, pair[0]), conjugate(v, pair[1])) for v in KERNEL
        )
        require(orbit <= pairs, (pair, "kernel orbit leaves braid set"))
        orbits.append(orbit)
        unseen -= orbit
    require(sorted(map(len, orbits)) == [4, 4], "two free V4 orbits")
    require(
        {cycle_type(next(iter(orbit))[0]) for orbit in orbits}
        == {(2, 1, 1), (4,)},
        "one orbit per branch",
    )
    return tuple(records), summary, tuple(orbits)


all_summaries = Counter()
ordered_pairs = 0
braided_pairs = 0
for tau_1 in TRANSPOSITIONS:
    for tau_2 in TRANSPOSITIONS:
        if tau_1 == tau_2:
            continue
        records, summary, _orbits = lift_records(tau_1, tau_2)
        ordered_pairs += 1
        braided_pairs += len(records)
        all_summaries.update(summary)

require(ordered_pairs == 6, "ordered downstairs pairs")
require(braided_pairs == 48, "all braided lifts")
require(
    all_summaries
    == Counter(
        {
            ((2, 1, 1), (2, 1, 1), 6, 1, 2, 1, 1): 24,
            ((4,), (4,), 24, 0, 4, 3, 1): 24,
        }
    ),
    "universal branch census",
)

# Canonical representatives and the global-owner hostile.
X_SPLIT = (1, 0, 2, 3)  # (1 2)
Y_SPLIT = (0, 2, 1, 3)  # (2 3)
Z_REMOTE = (3, 1, 2, 0)  # (1 4)
X_FULL = (1, 3, 0, 2)  # (1 2 4 3)
Y_FULL = (2, 3, 1, 0)  # (1 3 2 4)
require(
    compose(compose(X_SPLIT, Y_SPLIT), X_SPLIT)
    == compose(compose(Y_SPLIT, X_SPLIT), Y_SPLIT),
    "split representative braid",
)
require(
    compose(compose(X_FULL, Y_FULL), X_FULL)
    == compose(compose(Y_FULL, X_FULL), Y_FULL),
    "full representative braid",
)
require(len(generated_group((X_SPLIT, Y_SPLIT))) == 6, "split S3")
require(fixed_points(X_SPLIT) & fixed_points(Y_SPLIT) == {3}, "split owner")
require(len(generated_group((X_FULL, Y_FULL))) == 24, "full S4")
require(not (fixed_points(X_FULL) & fixed_points(Y_FULL)), "full no owner")
require(len(generated_group((X_SPLIT, Y_SPLIT, Z_REMOTE))) == 24, "remote monodromy hostile")

print("THM3037 cusp-braid S4 lift dichotomy exact certificate")
print(f"s4_order={len(S4)};s3_order={len(S3)};kernel_order={len(KERNEL)}")
print(f"downstairs_ordered_distinct_transposition_pairs={ordered_pairs}")
print("candidate_lift_pairs_per_downstairs_pair=16")
print("braided_lift_pairs_per_downstairs_pair=8")
print(f"all_braided_lift_pairs={braided_pairs}")
for key in sorted(all_summaries, key=repr):
    print(f"universal_class={key};count={all_summaries[key]}")
print("simultaneous_kernel_orbits_per_downstairs_pair=2;sizes=4,4")
print("split_branch=transposition,transposition;group=S3;common_fixed_sheet=1;cusp_width=2;disc_exponents=1,1")
print("full_branch=4cycle,4cycle;group=S4;common_fixed_sheet=0;cusp_width=4;disc_exponents=3,1")
print("mixed_braided_branch_count=0")
print("square_defect=0,0_or_distinct_nonzero_V4_generators")
print("braid_centre=(XY)^3=1;modular_relations=A^2=B^3=1;parabolic_AB=X^-1")
print("global_hostile=<X_split,Y_split,(1 4)>=S4;local_common_sheet_not_global")
print("all_exact_controls=PASS")
