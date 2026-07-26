#!/usr/bin/env python3
"""Exact controls for THM-2422's summand-parity and twin-center ancestry.

Universe
--------
* all integers and primes through ``CENTER_LIMIT + 1``;
* every A014574 center ``c <= CENTER_LIMIT``;
* the normalized nonexceptional necklace
      K = {k >= 1 : 6*k-1 and 6*k+1 are prime};
* summand targets through ``SUMMAND_LIMIT``; and
* the local residue carriers modulo the displayed primes.

The sieve, ancestry search, summand witnesses, and residue enumeration use
only exact integer/Boolean operations.  Every executable check calls
``require`` explicitly, so ``python3`` and ``python3 -O`` execute the same
verification.
"""

from collections import Counter
from math import isqrt


CENTER_LIMIT = 100_000_000
SUMMAND_LIMIT = 10_000
LOCAL_PRIMES = (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def prime_sieve(limit: int) -> bytearray:
    prime = bytearray(b"\x01") * (limit + 1)
    prime[0:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if prime[p]:
            start = p * p
            prime[start : limit + 1 : p] = b"\x00" * (
                (limit - start) // p + 1
            )
    return prime


prime = prime_sieve(CENTER_LIMIT + 1)
centers = [
    c
    for c in range(3, CENTER_LIMIT + 1)
    if prime[c - 1] and prime[c + 1]
]
require(
    centers[:21]
    == [
        4,
        6,
        12,
        18,
        30,
        42,
        60,
        72,
        102,
        108,
        138,
        150,
        180,
        192,
        198,
        228,
        240,
        270,
        282,
        312,
        348,
    ],
    "A014574 startup prefix drifted",
)
require(
    all(c == 4 or c % 6 == 0 for c in centers),
    "a nonexceptional twin center was not divisible by six",
)

center_set = set(centers)
immediate_successes = []
immediate_failures = []
for index in range(1, len(centers)):
    gap = centers[index] - centers[index - 1]
    record = (index + 1, centers[index], centers[index - 1], gap)
    if gap in center_set:
        immediate_successes.append(record)
    else:
        immediate_failures.append(record)

require(
    immediate_failures[0] == (2, 6, 4, 2),
    "the exceptional 4-to-6 startup transition changed",
)
require(
    immediate_failures[1] == (21, 348, 312, 36),
    "the first post-startup failure of immediate ancestry changed",
)
require(
    all(
        centers[index] - centers[index - 1] in center_set
        for index in range(2, 20)
    ),
    "immediate ancestry failed before center 348",
)
require(36 not in center_set, "the 348=312+36 hostile disappeared")

# Remove the exceptional center 4 before scaling.
necklace = [c // 6 for c in centers if c % 6 == 0]
necklace_set = set(necklace)
require(
    necklace[:20]
    == [
        1,
        2,
        3,
        5,
        7,
        10,
        12,
        17,
        18,
        23,
        25,
        30,
        32,
        33,
        38,
        40,
        45,
        47,
        52,
        58,
    ],
    "normalized twin-center prefix drifted",
)
require(
    all(k == 1 or k % 5 in {0, 2, 3} for k in necklace),
    "the modulo-five twin-center filter failed",
)
mod5_parent_channels = {
    residue: sum(
        1
        for left in {0, 2, 3}
        for right in {0, 2, 3}
        if (left + right) % 5 == residue
    )
    for residue in range(5)
}
require(
    mod5_parent_channels == {0: 3, 1: 1, 2: 2, 3: 2, 4: 1},
    "the modulo-five parent-channel profile changed",
)
known_hyp1994_holes = (16, 67, 86, 131, 151, 186, 191, 211, 226, 541, 701)
require(
    sum(hole % 5 == 1 for hole in known_hyp1994_holes) == 10,
    "the modulo-five profile of the recorded HYP-1994 holes changed",
)

# Canonical parent: the least a in K with a < k-a and k-a in K.
# Positivity makes both parents earlier than the target automatically.
canonical_parents = {}
largest_search = (0, 0, 0)
for k in necklace:
    parent = None
    tested = 0
    for a in necklace:
        if 2 * a >= k:
            break
        tested += 1
        if k - a in necklace_set:
            parent = (a, k - a)
            break
    if parent is not None:
        canonical_parents[k] = parent
    if tested > largest_search[0]:
        largest_search = (tested, k, 0 if parent is None else parent[0])

require(
    set(necklace) - set(canonical_parents) == {1, 2},
    "a checked normalized center k>=3 lacked distinct K+K parents",
)
require(
    canonical_parents[3] == (1, 2),
    "the first distinct parent pair changed",
)
require(
    canonical_parents[58] == (18, 40),
    "the canonical repair 58=18+40 changed",
)
require(
    58 - 52 == 6 and 6 not in necklace_set,
    "the normalized immediate-ancestry hostile changed",
)

# Independent traversal: scan possible second parents from the latest
# necklace vertex backward.  This has the opposite order from the primary
# least-first-parent scan and must recover the same canonical pair.
reverse_canonical_parents = {}
reverse_repair_depths = {}
for target_index, k in enumerate(necklace):
    for parent_index in range(target_index - 1, -1, -1):
        b = necklace[parent_index]
        a = k - b
        if a >= b:
            break
        if a in necklace_set:
            reverse_canonical_parents[k] = (a, b)
            reverse_repair_depths[k] = target_index - 1 - parent_index
            break

require(
    reverse_canonical_parents == canonical_parents,
    "forward and reverse canonical-parent scans diverged",
)

# The least first parent is equivalently the latest possible second parent.
# Its distance behind the immediate predecessor is the canonical repair depth.
necklace_index = {value: index for index, value in enumerate(necklace)}
repair_depths = {
    k: necklace_index[k] - 1 - necklace_index[b]
    for k, (_, b) in canonical_parents.items()
}
repair_depth_histogram = Counter(repair_depths.values())
max_repair_depth = max(repair_depths.values())
max_repair_targets = [
    k for k, depth in repair_depths.items() if depth == max_repair_depth
]
require(repair_depths[58] == 3, "the 58 repair depth changed")
require(
    reverse_repair_depths == repair_depths,
    "forward and reverse repair depths diverged",
)
require(
    repair_depth_histogram[0] + 1 == len(immediate_successes),
    "depth-zero ancestry and the immediate-gap count diverged beyond k=2=1+1",
)

# The distinct-summand closure from seeds {2,3}.
missing = {1, 4, 6}
for z in range(1, SUMMAND_LIMIT + 1):
    if z in {2, 3}:
        continue
    if z == 5:
        witness = (2, 3)
    elif z == 7:
        witness = (2, 5)
    elif z == 8:
        witness = (3, 5)
    elif z >= 9:
        witness = (2, z - 2)
    else:
        witness = None

    if z in missing:
        require(witness is None, f"missing target {z} received a witness")
        for a in range(1, z):
            b = z - a
            if a < b:
                require(
                    a in missing or b in missing,
                    f"missing target {z} has live parents {(a, b)}",
                )
    else:
        require(witness is not None, f"live target {z} lacked a witness")
        a, b = witness
        require(
            1 <= a < b < z
            and a + b == z
            and a not in missing
            and b not in missing,
            f"invalid closure witness for {z}: {witness}",
        )

# Exact period-two lattice count for the initial-segment filtration.
summand_total = 0
for z in range(1, SUMMAND_LIMIT + 1):
    direct_count = sum(1 for a in range(1, z) if a < z - a)
    formula_count = (z - 1) // 2
    parity_formula = (2 * z - 3 - (-1) ** z) // 4
    require(
        direct_count == formula_count == parity_formula,
        f"summand fibre count failed at {z}",
    )
    summand_total += direct_count
    require(
        summand_total == ((z - 1) * (z - 1)) // 4,
        f"cumulative summand count failed at {z}",
    )

# Once z>=13, deleting the missing dependency module {1,4,6} deletes
# exactly three labelled parent pairs and no others.
for z in range(13, SUMMAND_LIMIT + 1):
    restricted_count = sum(
        1
        for a in range(1, z)
        if a < z - a and a not in missing and z - a not in missing
    )
    require(
        restricted_count == (z - 1) // 2 - 3,
        f"restricted summand fibre failed at {z}",
    )

# Weak/strict commutative fibres expose the swap-fixed diagonal.  Addition
# sees the even targets; multiplication sees the perfect squares.  Deleting
# the multiplicative unit then makes primes exactly the zero proper-divisor
# witnesses, but the diagonal must be retained to distinguish 9 from a prime.
for z in range(2, SUMMAND_LIMIT + 1):
    root = isqrt(z)
    square = int(root * root == z)
    tau = sum(
        1 if divisor * divisor == z else 2
        for divisor in range(1, root + 1)
        if z % divisor == 0
    )
    weak_additive = sum(1 for a in range(1, z) if a <= z - a)
    strict_additive = sum(1 for a in range(1, z) if a < z - a)
    weak_multiplicative = sum(
        1 for a in range(1, root + 1) if z % a == 0
    )
    strict_multiplicative = weak_multiplicative - square
    proper_weak_multiplicative = weak_multiplicative - 1
    proper_strict_multiplicative = strict_multiplicative - 1

    require(
        weak_additive == z // 2
        and strict_additive == z // 2 - int(z % 2 == 0),
        f"weak/strict additive fibres failed at {z}",
    )
    require(
        2 * weak_multiplicative == tau + square
        and 2 * strict_multiplicative == tau - square,
        f"weak/strict multiplicative fibres failed at {z}",
    )
    require(
        2 * proper_weak_multiplicative == tau - 2 + square
        and 2 * proper_strict_multiplicative == tau - 2 - square,
        f"proper multiplicative fibres failed at {z}",
    )
    require(
        (tau - 2 == 0) == bool(prime[z]),
        f"prime zero-defect test failed at {z}",
    )
    require(
        (proper_weak_multiplicative == 0) == bool(prime[z]),
        f"weak proper-factor atom test failed at {z}",
    )
    require(
        (proper_strict_multiplicative == 0)
        == (bool(prime[z]) or bool(square and prime[root])),
        f"strict proper-factor prime-square ambiguity failed at {z}",
    )

for square_hostile in (4, 9):
    require(
        not prime[square_hostile]
        and sum(
            1
            for a in range(2, isqrt(square_hostile) + 1)
            if a < square_hostile // a and square_hostile % a == 0
        )
        == 0,
        f"the square-diagonal hostile at {square_hostile} changed",
    )

# The Fibonacci chain is an embedded ancestry path, not the full closure.
fibonacci_spine = [2, 3]
while fibonacci_spine[-1] + fibonacci_spine[-2] <= SUMMAND_LIMIT:
    fibonacci_spine.append(fibonacci_spine[-1] + fibonacci_spine[-2])
for index in range(2, len(fibonacci_spine)):
    require(
        fibonacci_spine[index]
        == fibonacci_spine[index - 1] + fibonacci_spine[index - 2],
        "Fibonacci spine recurrence failed",
    )
    require(
        fibonacci_spine[index] not in missing,
        "Fibonacci spine entered the missing dependency module",
    )

parents_of_13 = [
    (a, 13 - a)
    for a in range(1, 13)
    if a < 13 - a and a not in missing and 13 - a not in missing
]
require(
    parents_of_13 == [(2, 11), (3, 10), (5, 8)],
    "the first displayed nonunique Fibonacci ancestry changed",
)

# Synchronous generation has one further startup tooth and then becomes an
# exact dyadic interval frontier.  Here S_(t+1) uses only pairs already in S_t.
synchronous_stages = [{2, 3}]
for _ in range(10):
    old_stage = synchronous_stages[-1]
    values = sorted(old_stage)
    next_stage = set(old_stage)
    for left_index, left in enumerate(values):
        for right in values[left_index + 1 :]:
            next_stage.add(left + right)
    synchronous_stages.append(next_stage)

require(synchronous_stages[0] == {2, 3}, "synchronous stage S0 drifted")
require(synchronous_stages[1] == {2, 3, 5}, "synchronous stage S1 drifted")
require(
    synchronous_stages[2] == {2, 3, 5, 7, 8},
    "synchronous stage S2 drifted",
)
require(
    synchronous_stages[3]
    == {2, 3, 5, 7, 8, 9, 10, 11, 12, 13, 15},
    "synchronous stage S3 drifted",
)
require(
    synchronous_stages[4] == {2, 3, 5} | set(range(7, 29)),
    "the exceptional S3-to-S4 interval fill drifted",
)
for stage_index in range(4, len(synchronous_stages)):
    frontier = 27 * (2 ** (stage_index - 4)) + 1
    require(
        synchronous_stages[stage_index]
        == {2, 3, 5} | set(range(7, frontier + 1)),
        f"dyadic synchronous frontier failed at stage {stage_index}",
    )

# Parent triples a,b,a+b in K are prime sextuples
#   6a+-1, 6b+-1, 6(a+b)+-1.
# After the invertible scaling x=6a, y=6b modulo p>=5, the forbidden
# lines are x=+-1, y=+-1, x+y=+-1.  Their complement has p^2-6p+12 points.
local_counts = []
for p in LOCAL_PRIMES:
    allowed = 0
    for x in range(p):
        for y in range(p):
            if (
                x not in {1, p - 1}
                and y not in {1, p - 1}
                and (x + y) % p not in {1, p - 1}
            ):
                allowed += 1
    require(
        allowed == p * p - 6 * p + 12,
        f"local prime-sextuple count failed modulo {p}",
    )
    target_profile = {}
    for target in range(p):
        if target in {1, p - 1}:
            continue
        target_count = sum(
            1
            for x in range(p)
            if x not in {1, p - 1}
            and (target - x) % p not in {1, p - 1}
        )
        if target == 0:
            expected = p - 2
        elif target in {2, p - 2}:
            expected = p - 3
        else:
            expected = p - 4
        require(
            target_count == expected,
            f"targetwise sextuple count failed modulo {p} at {target}",
        )
        target_profile[target] = target_count
    require(
        sum(target_profile.values()) == allowed,
        f"targetwise and aggregate sextuple counts diverged modulo {p}",
    )
    local_counts.append((p, allowed))

parent_preview = ",".join(
    f"{k}={a}+{b}"
    for k, (a, b) in list(canonical_parents.items())[:20]
)
local_preview = ",".join(f"{p}:{count}" for p, count in local_counts)
repair_preview = ",".join(
    f"{depth}:{repair_depth_histogram[depth]}" for depth in range(11)
)

print("THM-2422 twin-center additive ancestry exact controls")
print(f"center_limit={CENTER_LIMIT}")
print(f"A014574_centers={len(centers)} last_center={centers[-1]}")
print(f"normalized_K={len(necklace)} last_k={necklace[-1]}")
print(
    "immediate_ancestry="
    f"successes:{len(immediate_successes)} "
    f"failures:{len(immediate_failures)}"
)
print("first_raw_failure=6=4+2 with 2_not_A014574")
print("first_post_startup_failure=348=312+36 with 36_not_A014574")
print("normalized_hostile=58=52+6 with 6_not_K")
print("mod5_K_filter=K_minus_1_subset_0_2_3")
print("mod5_parent_channels=0:3,1:1,2:2,3:2,4:1")
print("HYP1994_recorded_holes_mod5_1=10_of_11")
print(
    "distinct_K_plus_K="
    f"all_checked_k_ge_3:{len(necklace) - 2} "
    "exceptions:1,2"
)
print("first_hostile_repair=58=18+40; center_form:348=108+240")
print(f"canonical_parent_preview={parent_preview}")
print(
    "largest_canonical_search="
    f"candidates:{largest_search[0]} "
    f"target_k:{largest_search[1]} "
    f"least_parent:{largest_search[2]}"
)
print(
    "repair_depth_histogram_0_to_10="
    f"{repair_preview};gt10:"
    f"{sum(count for depth, count in repair_depth_histogram.items() if depth > 10)}"
)
print(
    "max_repair_depth="
    f"{max_repair_depth} targets:"
    + ",".join(str(target) for target in max_repair_targets)
)
print(
    "reverse_parent_scan="
    f"matched:{len(reverse_canonical_parents)} max_depth:{max_repair_depth}"
)
print(
    "summand_closure_from_2_3="
    f"[1,{SUMMAND_LIMIT}]_minus_1_4_6"
)
print(
    "summand_fibre="
    "r(z)=floor((z-1)/2)=(2z-3-(-1)^z)/4"
)
print("summand_cumulative=R(N)=floor((N-1)^2/4)")
print("restricted_summand_fibre_for_z_ge_13=r(z)-3")
print(
    "swap_fixed_diagonals="
    "addition:even_targets multiplication:perfect_squares"
)
print("multiplicative_atoms=prime_iff_tau(z)-2=0")
print("strict_factor_empty=prime_or_prime_square")
print("diagonal_hostiles=4_and_9_have_no_strict_proper_factor_pair")
print(
    "fibonacci_spine="
    + ",".join(str(value) for value in fibonacci_spine)
)
print("fibonacci_target_13_parents=2+11,3+10,5+8")
print(
    "synchronous_stage_max="
    + ",".join(
        f"{stage_index}:{max(stage)}"
        for stage_index, stage in enumerate(synchronous_stages)
    )
)
print(
    "synchronous_stage_size="
    + ",".join(
        f"{stage_index}:{len(stage)}"
        for stage_index, stage in enumerate(synchronous_stages)
    )
)
print("dyadic_frontier=M_t=27*2^(t-4)+1 for 4<=t<=10")
print(f"local_sextuple_counts={local_preview}")
print("all_checks=PASS")
