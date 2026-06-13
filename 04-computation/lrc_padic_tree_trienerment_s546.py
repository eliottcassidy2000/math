#!/usr/bin/env python3
"""p-adic tree-cover trienerments for LRC sieve channels.

This is the concrete cover-core companion to the recent p-adic tree threads.
For a speed set V with |V|=n-1, every integer q with 2 <= q <= n is a
product p-adic zero-branch node: the residue ball 0 mod q in the product of
the p-adic trees for primes dividing q.  Its mass is

    z_q(V) = #{v in V : q divides v}.

If z_q=0 and q<n, then t=1/q is an open THM-369 sieve witness because every
runner is at distance at least 1/q > 1/n from the observer.  If z_n=0, the same
probe is a compactified wall witness at t=1/n.

Tournament Analysis declaration:
  * vertices: product p-adic zero-branch obligations q=2,...,n;
  * pairwise observable: zero-branch cover mass z_q plus divisibility depth;
  * switch/gauge: lower cover mass wins; equal-mass comparable nodes orient
    toward the deeper branch; equal-mass incomparable nodes are ties;
  * tie Hamiltonian path: increasing q;
  * fingerprints: tie count, ternary B1, score histogram, strict 3-cycles,
    SCCs, Hamiltonian-path count, and cover-core/min-set-cover statistics.

Assumption challenge:
  The vertices are not runners, arcs, or Gabor cells.  They are proof
  obligations in the p-adic product tree.  The quotient preserves the THM-369
  sieve predicate and the compactified q=n wall predicate.  It destroys exact
  event order, circular gap shape, runner identity except for carrier roles,
  and non-sieve LRC witnesses.  The challenged assumption is that the p-adic
  tree only records denominator/channel rank; here it also records which
  individual speeds are forced to carry the zero-branch cover.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from math import gcd
from statistics import mean
from typing import Iterable


CONFIGS = [(6, 8), (8, 10), (10, 12), (12, 14), (14, 16), (16, 18), (18, 20)]


@dataclass(frozen=True)
class TreeFingerprint:
    node_count: int
    edge_count: int
    tie_count: int
    ternary_b1: int
    score_histogram: tuple[tuple[int, int], ...]
    strict_3cycles: int
    scc_count: int
    largest_scc: int
    hamiltonian_paths: int


@dataclass(frozen=True)
class CoverRecord:
    n: int
    speeds: tuple[int, ...]
    masses: tuple[int, ...]
    open_empty: tuple[int, ...]
    compact_empty: tuple[int, ...]
    open_survivor: bool
    compact_survivor: bool
    singleton_nodes: tuple[int, ...]
    singleton_carriers: tuple[tuple[int, int], ...]
    forced_carrier_count: int
    max_singleton_load: int
    min_open_cover_size: int | None
    prime_power_empty: tuple[int, ...]
    fingerprint: TreeFingerprint


def gcd_many(values: Iterable[int]) -> int:
    result = 0
    for value in values:
        result = gcd(result, value)
    return result


def prime_factorization(value: int) -> tuple[tuple[int, int], ...]:
    factors: list[tuple[int, int]] = []
    divisor = 2
    n = value
    while divisor * divisor <= n:
        if n % divisor == 0:
            exp = 0
            while n % divisor == 0:
                n //= divisor
                exp += 1
            factors.append((divisor, exp))
        divisor += 1 if divisor == 2 else 2
    if n > 1:
        factors.append((n, 1))
    return tuple(factors)


def is_prime_power(value: int) -> bool:
    return len(prime_factorization(value)) == 1


def omega(value: int) -> int:
    return len(prime_factorization(value))


def nodes(n: int, compact: bool = True) -> tuple[int, ...]:
    end = n if compact else n - 1
    return tuple(range(2, end + 1))


def primitive_speed_sets(n: int, max_speed: int) -> Iterable[tuple[int, ...]]:
    for speeds in combinations(range(1, max_speed + 1), n - 1):
        if gcd_many(speeds) == 1:
            yield speeds


def zero_branch_masses(n: int, speeds: tuple[int, ...], compact: bool = True) -> tuple[int, ...]:
    return tuple(sum(1 for speed in speeds if speed % q == 0) for q in nodes(n, compact))


def compare_nodes(q_left: int, mass_left: int, q_right: int, mass_right: int) -> int:
    """Return 1 if left beats right, -1 if right beats left, 0 for tie."""
    if mass_left != mass_right:
        return 1 if mass_left < mass_right else -1
    if q_left != q_right:
        if q_right % q_left == 0:
            return -1
        if q_left % q_right == 0:
            return 1
    return 0


def tree_fingerprint(n: int, masses: tuple[int, ...]) -> TreeFingerprint:
    qs = nodes(n, compact=True)
    vertex_count = len(qs)
    edge_count = vertex_count * (vertex_count - 1) // 2
    adjacency = [0] * vertex_count
    strict_orientation: dict[tuple[int, int], int] = {}
    tie_count = 0

    for i, j in combinations(range(vertex_count), 2):
        cmp_value = compare_nodes(qs[i], masses[i], qs[j], masses[j])
        if cmp_value > 0:
            adjacency[i] |= 1 << j
            strict_orientation[(i, j)] = 1
        elif cmp_value < 0:
            adjacency[j] |= 1 << i
            strict_orientation[(i, j)] = -1
        else:
            tie_count += 1
            adjacency[i] |= 1 << j
            adjacency[j] |= 1 << i
            strict_orientation[(i, j)] = 0

    outdegrees = [adjacency[i].bit_count() for i in range(vertex_count)]
    scc_sizes = strongly_connected_component_sizes(adjacency, vertex_count)
    return TreeFingerprint(
        node_count=vertex_count,
        edge_count=edge_count,
        tie_count=tie_count,
        ternary_b1=2 * edge_count - 3 * tie_count,
        score_histogram=tuple(sorted(Counter(outdegrees).items())),
        strict_3cycles=strict_3cycle_count(vertex_count, strict_orientation),
        scc_count=len(scc_sizes),
        largest_scc=max(scc_sizes) if scc_sizes else 0,
        hamiltonian_paths=hamiltonian_path_count(adjacency, vertex_count),
    )


def strict_3cycle_count(vertex_count: int, strict_orientation: dict[tuple[int, int], int]) -> int:
    cycles = 0
    for a, b, c in combinations(range(vertex_count), 3):
        ab = strict_orientation[(a, b)]
        ac = strict_orientation[(a, c)]
        bc = strict_orientation[(b, c)]
        if ab == 0 or ac == 0 or bc == 0:
            continue
        if (ab, ac, bc) in ((1, -1, 1), (-1, 1, -1)):
            cycles += 1
    return cycles


def strongly_connected_component_sizes(adjacency: list[int], vertex_count: int) -> list[int]:
    reverse = [0] * vertex_count
    for source, mask in enumerate(adjacency):
        remaining = mask
        while remaining:
            low_bit = remaining & -remaining
            target = low_bit.bit_length() - 1
            reverse[target] |= 1 << source
            remaining ^= low_bit

    seen = 0
    order: list[int] = []

    def dfs_forward(start: int) -> None:
        nonlocal seen
        stack = [(start, False)]
        while stack:
            node, expanded = stack.pop()
            if expanded:
                order.append(node)
                continue
            if seen & (1 << node):
                continue
            seen |= 1 << node
            stack.append((node, True))
            remaining = adjacency[node] & ~seen
            while remaining:
                low_bit = remaining & -remaining
                stack.append((low_bit.bit_length() - 1, False))
                remaining ^= low_bit

    for node in range(vertex_count):
        if not (seen & (1 << node)):
            dfs_forward(node)

    assigned = 0
    sizes: list[int] = []
    for start in reversed(order):
        if assigned & (1 << start):
            continue
        size = 0
        stack = [start]
        assigned |= 1 << start
        while stack:
            node = stack.pop()
            size += 1
            remaining = reverse[node] & ~assigned
            while remaining:
                low_bit = remaining & -remaining
                target = low_bit.bit_length() - 1
                assigned |= low_bit
                stack.append(target)
                remaining ^= low_bit
        sizes.append(size)
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(adjacency: list[int], vertex_count: int) -> int:
    full_mask = (1 << vertex_count) - 1
    dp = [[0] * vertex_count for _ in range(1 << vertex_count)]
    for vertex in range(vertex_count):
        dp[1 << vertex][vertex] = 1
    for mask in range(1 << vertex_count):
        for end in range(vertex_count):
            count = dp[mask][end]
            if not count:
                continue
            available = adjacency[end] & ~mask
            while available:
                low_bit = available & -available
                nxt = low_bit.bit_length() - 1
                dp[mask | low_bit][nxt] += count
                available ^= low_bit
    return sum(dp[full_mask])


def min_set_cover_size(open_qs: tuple[int, ...], speeds: tuple[int, ...]) -> int | None:
    full = (1 << len(open_qs)) - 1
    masks = []
    for speed in speeds:
        mask = 0
        for idx, q in enumerate(open_qs):
            if speed % q == 0:
                mask |= 1 << idx
        masks.append(mask)
    if 0 not in (full,):
        covered = 0
        for mask in masks:
            covered |= mask
        if covered != full:
            return None
    large = len(speeds) + 1
    dp = [large] * (1 << len(open_qs))
    dp[0] = 0
    for mask in masks:
        if mask == 0:
            continue
        for state in range(full, -1, -1):
            if dp[state] < large:
                nxt = state | mask
                if dp[state] + 1 < dp[nxt]:
                    dp[nxt] = dp[state] + 1
    return dp[full] if dp[full] < large else None


def cover_record(
    n: int,
    speeds: tuple[int, ...],
    fp_cache: dict[tuple[int, tuple[int, ...]], TreeFingerprint],
) -> CoverRecord:
    compact_qs = nodes(n, compact=True)
    open_qs = nodes(n, compact=False)
    masses = zero_branch_masses(n, speeds, compact=True)
    open_masses = masses[: len(open_qs)]
    open_empty = tuple(q for q, mass in zip(open_qs, open_masses) if mass == 0)
    compact_empty = tuple(q for q, mass in zip(compact_qs, masses) if mass == 0)
    singleton_nodes = tuple(q for q, mass in zip(compact_qs, masses) if mass == 1)

    carrier_counter: Counter[int] = Counter()
    singleton_carriers = []
    for q in singleton_nodes:
        carriers = [speed for speed in speeds if speed % q == 0]
        if len(carriers) == 1:
            carrier_counter[carriers[0]] += 1
            singleton_carriers.append((q, carriers[0]))

    fp_key = (n, masses)
    if fp_key not in fp_cache:
        fp_cache[fp_key] = tree_fingerprint(n, masses)
    prime_power_empty = tuple(q for q in compact_empty if is_prime_power(q))

    return CoverRecord(
        n=n,
        speeds=speeds,
        masses=masses,
        open_empty=open_empty,
        compact_empty=compact_empty,
        open_survivor=not open_empty,
        compact_survivor=not compact_empty,
        singleton_nodes=singleton_nodes,
        singleton_carriers=tuple(singleton_carriers),
        forced_carrier_count=len(carrier_counter),
        max_singleton_load=max(carrier_counter.values()) if carrier_counter else 0,
        min_open_cover_size=min_set_cover_size(open_qs, speeds),
        prime_power_empty=prime_power_empty,
        fingerprint=fp_cache[fp_key],
    )


def scan_config(n: int, max_speed: int) -> str:
    fp_cache: dict[tuple[int, tuple[int, ...]], TreeFingerprint] = {}
    records = [cover_record(n, speeds, fp_cache) for speeds in primitive_speed_sets(n, max_speed)]
    open_survivors = [record for record in records if record.open_survivor]
    compact_survivors = [record for record in records if record.compact_survivor]

    lines = [
        f"=== n={n}, max_speed={max_speed} ===",
        f"primitive speed sets: {len(records)}",
        f"channel rank omega(n/2): {omega(n // 2)}",
        f"open sieve survivors q=2..{n-1}: {len(open_survivors)}",
        f"compact q=2..{n} survivors: {len(compact_survivors)}",
        f"best open witness q distribution: {best_witness_distribution(records)}",
        f"prime-power empty distribution: {histogram(len(r.prime_power_empty) for r in records)}",
    ]
    lines.extend(record_summary("all records", records))
    lines.extend(record_summary("open survivors", open_survivors))
    lines.extend(fingerprint_summary(records, open_survivors))
    lines.extend(named_examples(n, max_speed, fp_cache))
    return "\n".join(lines)


def best_witness_distribution(records: list[CoverRecord]) -> dict[str, int]:
    counter: Counter[str] = Counter()
    for record in records:
        if record.open_empty:
            counter[str(min(record.open_empty))] += 1
        elif record.compact_empty:
            counter["wall"] += 1
        else:
            counter["none"] += 1
    return dict(sorted(counter.items(), key=lambda item: (item[0] == "none", item[0] == "wall", item[0])))


def histogram(values: Iterable[int | None]) -> dict[str, int]:
    counter: Counter[str] = Counter("None" if value is None else str(value) for value in values)
    return dict(sorted(counter.items(), key=lambda item: (item[0] == "None", int(item[0]) if item[0] != "None" else -1)))


def record_summary(label: str, records: list[CoverRecord]) -> list[str]:
    if not records:
        return [f"{label}: none"]
    return [
        f"{label}:",
        f"  open empty count avg/range: {avg_range([len(r.open_empty) for r in records])}",
        f"  compact empty count avg/range: {avg_range([len(r.compact_empty) for r in records])}",
        f"  singleton nodes avg/range: {avg_range([len(r.singleton_nodes) for r in records])}",
        f"  forced carrier count avg/range: {avg_range([r.forced_carrier_count for r in records])}",
        f"  max singleton load avg/range: {avg_range([r.max_singleton_load for r in records])}",
        f"  min open cover size distribution: {histogram(r.min_open_cover_size for r in records)}",
        f"  tree tie count avg/range: {avg_range([r.fingerprint.tie_count for r in records])}",
        f"  tree B1 avg/range: {avg_range([r.fingerprint.ternary_b1 for r in records])}",
        f"  strict 3-cycles avg/range: {avg_range([r.fingerprint.strict_3cycles for r in records])}",
        f"  HP count avg/range: {avg_range([r.fingerprint.hamiltonian_paths for r in records])}",
    ]


def avg_range(values: list[int]) -> str:
    return f"{mean(values):.3f} / {min(values)}..{max(values)}"


def fingerprint_key(record: CoverRecord) -> tuple[object, ...]:
    fp = record.fingerprint
    return (
        len(record.open_empty),
        len(record.compact_empty),
        len(record.singleton_nodes),
        record.forced_carrier_count,
        record.max_singleton_load,
        record.min_open_cover_size,
        fp.tie_count,
        fp.ternary_b1,
        fp.score_histogram,
        fp.scc_count,
        fp.hamiltonian_paths,
    )


def fingerprint_summary(records: list[CoverRecord], open_survivors: list[CoverRecord]) -> list[str]:
    classes = defaultdict(list)
    for record in records:
        classes[fingerprint_key(record)].append(record)
    survivor_classes = {fingerprint_key(record) for record in open_survivors}
    mixed = 0
    for key, members in classes.items():
        outcomes = {member.open_survivor for member in members}
        if len(outcomes) > 1:
            mixed += 1
    return [
        "tree-trienerment fingerprint classes:",
        f"  total classes: {len(classes)}",
        f"  open-survivor classes: {len(survivor_classes)}",
        f"  mixed survivor/non-survivor classes: {mixed}",
    ]


def named_examples(
    n: int,
    max_speed: int,
    fp_cache: dict[tuple[int, tuple[int, ...]], TreeFingerprint],
) -> list[str]:
    candidates = {
        "AP": tuple(range(1, n)),
        "shifted": tuple(range(2, n + 1)),
    }
    if n % 2 == 0:
        apex = n // 2
        candidates["apex_replaced"] = tuple(sorted(set(range(1, n)) - {apex} | {n}))

    lines = ["named examples:"]
    for name, speeds in candidates.items():
        if len(speeds) != n - 1 or max(speeds) > max_speed or gcd_many(speeds) != 1:
            lines.append(f"  {name}: not in scan box")
            continue
        record = cover_record(n, speeds, fp_cache)
        lines.append(
            f"  {name}: open_survivor={record.open_survivor}, compact_survivor={record.compact_survivor}, "
            f"open_empty={record.open_empty}, compact_empty={record.compact_empty}, "
            f"singletons={record.singleton_nodes}, forced_carriers={record.singleton_carriers}, "
            f"min_cover={record.min_open_cover_size}, ties={record.fingerprint.tie_count}, "
            f"HP={record.fingerprint.hamiltonian_paths}"
        )
    return lines


def main() -> None:
    print("LRC p-adic tree-cover trienerment scan, S546")
    print("Nodes q=2..n are product p-adic zero branches; z_q counts speeds divisible by q.")
    print("Empty q<n gives open THM-369 witness t=1/q; empty q=n is compactified wall.")
    print()
    for config in CONFIGS:
        print(scan_config(*config))
        print()


if __name__ == "__main__":
    main()
