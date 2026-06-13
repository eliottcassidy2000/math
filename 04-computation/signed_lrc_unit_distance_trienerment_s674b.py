#!/usr/bin/env python3
"""S674b: signed LRC gauges and unit-distance trienerment edges.

This session tests two reframes from the prompt.

LRC signed speeds:
  The observer predicate is invariant under independent sign changes
  v_i -> +/- v_i, because ||v_i t|| = ||-v_i t||.  Runner-runner
  tournaments are not invariant: alternating signs turn many pair clocks
  from differences into sums.

Unit distance:
  Distances are treated as a three-state pair observable:
      S = short  (< 1), U = unit (= 1), L = long (> 1).
  S and L are oriented in opposite directions relative to a base order.
  U is a trienerment/tie layer, resolved either positive or negative.

Tournament Analysis / assumption challenge:
  For LRC, vertices considered included runners, signed runners, pair clocks,
  residues, fixed circle sections, wall crossings, Fourier modes, and proof
  obligations.  This script chooses pair-clock and lens vertices because the
  observer predicate is sign-invariant while the pair tournament is not.
  For unit distance, vertices considered included points, distance states,
  unit graph edges, nonunit edges, Hamiltonian path obligations, deletion
  owners, and construction moves.  This script chooses distance states and
  mapping lenses because the preserved predicate is whether the equality
  layer, not just a binary flip, carries the Hamiltonian path.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd


def bitcount(x: int) -> int:
    return x.bit_count()


def circular_distance(x: Fraction) -> Fraction:
    y = x % 1
    return min(y, 1 - y)


def tournament_from_positions(positions: tuple[Fraction, ...]) -> tuple[int, ...]:
    """Orient by clockwise half-circle; exact collisions use the base order."""

    n = len(positions)
    adj = [0] * n
    half = Fraction(1, 2)
    for i, j in combinations(range(n), 2):
        delta = (positions[j] - positions[i]) % 1
        if delta == 0 or delta == half or delta < half:
            adj[i] |= 1 << j
        else:
            adj[j] |= 1 << i
    return tuple(adj)


def directed_3cycles(adj: tuple[int, ...]) -> int:
    total = 0
    for a, b, c in combinations(range(len(adj)), 3):
        local = []
        for v in (a, b, c):
            s = 0
            for w in (a, b, c):
                if v != w and ((adj[v] >> w) & 1):
                    s += 1
            local.append(s)
        if sorted(local) == [1, 1, 1]:
            total += 1
    return total


def scc_sizes(adj: tuple[int, ...]) -> tuple[int, ...]:
    n = len(adj)
    rev = [0] * n
    for i, row in enumerate(adj):
        x = row
        while x:
            bit = x & -x
            j = bit.bit_length() - 1
            x -= bit
            rev[j] |= 1 << i

    seen = [False] * n
    order: list[int] = []

    def dfs1(v: int) -> None:
        seen[v] = True
        x = adj[v]
        while x:
            bit = x & -x
            w = bit.bit_length() - 1
            x -= bit
            if not seen[w]:
                dfs1(w)
        order.append(v)

    for v in range(n):
        if not seen[v]:
            dfs1(v)

    seen = [False] * n
    sizes = []

    def dfs2(v: int) -> int:
        seen[v] = True
        total = 1
        x = rev[v]
        while x:
            bit = x & -x
            w = bit.bit_length() - 1
            x -= bit
            if not seen[w]:
                total += dfs2(w)
        return total

    for v in reversed(order):
        if not seen[v]:
            sizes.append(dfs2(v))
    return tuple(sorted(sizes, reverse=True))


def ham_path_count(adj: tuple[int, ...]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            avail = adj[last] & ~mask
            while avail:
                bit = avail & -avail
                nxt = bit.bit_length() - 1
                avail -= bit
                dp[mask | bit][nxt] += count
    return sum(dp[-1])


def redei_insert_path(adj: tuple[int, ...]) -> tuple[int, ...]:
    path: list[int] = []
    for v in range(len(adj)):
        pos = 0
        while pos < len(path) and not ((adj[v] >> path[pos]) & 1):
            pos += 1
        path.insert(pos, v)
    return tuple(path)


def tournament_fingerprint(adj: tuple[int, ...], hp_limit: int = 9) -> dict[str, object]:
    outdegrees = tuple(bitcount(row) for row in adj)
    out = {
        "n": len(adj),
        "score_hist": tuple(sorted(Counter(outdegrees).items())),
        "c3": directed_3cycles(adj),
        "scc": scc_sizes(adj),
    }
    if len(adj) <= hp_limit:
        out["H"] = ham_path_count(adj)
    return out


def fmt_counter(counter: Counter, limit: int | None = None) -> str:
    items = sorted(counter.items())
    if limit is not None and len(items) > limit:
        head = items[:limit]
        return "{" + ", ".join(f"{k}:{v}" for k, v in head) + ", ...}"
    return "{" + ", ".join(f"{k}:{v}" for k, v in items) + "}"


def fp_key(adj: tuple[int, ...]) -> tuple[object, ...]:
    fp = tournament_fingerprint(adj, hp_limit=0)
    return (fp["score_hist"], fp["c3"], fp["scc"])


def sign_patterns(m: int) -> dict[str, tuple[int, ...]]:
    return {
        "all_plus": tuple(1 for _ in range(m)),
        "all_minus": tuple(-1 for _ in range(m)),
        "alternating_plus": tuple(1 if i % 2 == 0 else -1 for i in range(m)),
        "alternating_minus": tuple(-1 if i % 2 == 0 else 1 for i in range(m)),
    }


def signed_positions(speeds: tuple[int, ...], signs: tuple[int, ...], t: Fraction) -> tuple[Fraction, ...]:
    return tuple((signs[i] * speeds[i] * t) % 1 for i in range(len(speeds)))


def clearances(speeds: tuple[int, ...], signs: tuple[int, ...], times: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return tuple(min(circular_distance(signs[i] * speeds[i] * t) for i in range(len(speeds))) for t in times)


def pair_frequency_summary(
    speeds: tuple[int, ...], signs: tuple[int, ...], modulus: int
) -> dict[str, object]:
    abs_counts: Counter[int] = Counter()
    mod_counts: Counter[int] = Counter()
    same = 0
    opposite = 0
    same_abs: Counter[int] = Counter()
    opposite_abs: Counter[int] = Counter()
    for i, j in combinations(range(len(speeds)), 2):
        raw = signs[j] * speeds[j] - signs[i] * speeds[i]
        freq = abs(raw)
        abs_counts[freq] += 1
        mod_counts[raw % modulus] += 1
        if signs[i] == signs[j]:
            same += 1
            same_abs[freq] += 1
        else:
            opposite += 1
            opposite_abs[freq] += 1
    return {
        "same_pairs": same,
        "opposite_pairs": opposite,
        "abs_counts": abs_counts,
        "mod_counts": mod_counts,
        "same_abs": same_abs,
        "opposite_abs": opposite_abs,
        "gcd_mod_counts": Counter(gcd(k, modulus) for k in mod_counts for _ in range(mod_counts[k])),
    }


def lrc_signed_probe() -> list[str]:
    total_n = 14
    modulus = 2 * total_n - 1
    times = tuple(Fraction(k, 2 * modulus * total_n) for k in range(2 * modulus * total_n))
    families = {
        "AP13": tuple(range(1, 14)),
        "Vstar": (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24),
        "2AP13": tuple(2 * x for x in range(1, 14)),
    }
    lines = []
    lines.append("SIGNED LRC GAUGE PROBE")
    lines.append(f"total n={total_n}, moving runners=13, C=2n-1={modulus}, sample_times={len(times)}")
    lines.append("")

    for family_name, speeds in families.items():
        patterns = sign_patterns(len(speeds))
        base_clear = clearances(speeds, patterns["all_plus"], times)
        lines.append(f"{family_name}: speeds={speeds}")
        for pname, signs in patterns.items():
            seq = clearances(speeds, signs, times)
            invariant = seq == base_clear
            unique_fp: set[tuple[object, ...]] = set()
            unique_adj: set[tuple[int, ...]] = set()
            c3_values: Counter[int] = Counter()
            scc_values: Counter[tuple[int, ...]] = Counter()
            for t in times:
                adj = tournament_from_positions(signed_positions(speeds, signs, t))
                unique_adj.add(adj)
                key = fp_key(adj)
                unique_fp.add(key)
                c3_values[key[1]] += 1
                scc_values[key[2]] += 1
            freq = pair_frequency_summary(speeds, signs, modulus)
            lines.append(
                "  "
                + pname
                + f": observer_invariant={invariant}, unique_tournaments={len(unique_adj)}, "
                + f"unique_fingerprints={len(unique_fp)}, c3_values={fmt_counter(c3_values, 8)}, "
                + f"scc_shapes={fmt_counter(scc_values, 5)}"
            )
            lines.append(
                "    pair clocks: "
                + f"same={freq['same_pairs']}, opposite={freq['opposite_pairs']}, "
                + f"abs={fmt_counter(freq['abs_counts'], 16)}, mod{modulus}={fmt_counter(freq['mod_counts'], 16)}"
            )
            if pname.startswith("alternating"):
                lines.append(
                    "    alternating split: same-sign differences="
                    + fmt_counter(freq["same_abs"], 12)
                    + "; opposite-sign sums="
                    + fmt_counter(freq["opposite_abs"], 16)
                )

        mismatches = 0
        c3_delta = Counter()
        alt = patterns["alternating_plus"]
        plus = patterns["all_plus"]
        for t in times:
            a = tournament_from_positions(signed_positions(speeds, plus, t))
            b = tournament_from_positions(signed_positions(speeds, alt, t))
            ka = fp_key(a)
            kb = fp_key(b)
            if ka != kb:
                mismatches += 1
            c3_delta[kb[1] - ka[1]] += 1
        lines.append(
            f"  all_plus vs alternating_plus: fingerprint_mismatches={mismatches}/{len(times)}, "
            + f"c3_delta={fmt_counter(c3_delta, 16)}"
        )
        lines.append("")

    lines.append("LRC reading: signed direction is a gauge symmetry for the observer,")
    lines.append("but alternating signs convert many runner-runner clocks from differences")
    lines.append("into pair sums.  For AP13 the opposite-sign clocks are exactly the odd")
    lines.append("interior sums 3,5,...,25 inside C=27; the missing boundary shells 1 and")
    lines.append("27 are the address tax that the observer projection erases.")
    lines.append("")
    return lines


@dataclass(frozen=True)
class MetricPointSet:
    name: str
    points: tuple[tuple[Fraction, Fraction], ...]
    metric: str


def sqdist_cart(a: tuple[Fraction, Fraction], b: tuple[Fraction, Fraction]) -> Fraction:
    dx = b[0] - a[0]
    dy = b[1] - a[1]
    return dx * dx + dy * dy


def sqdist_axial(a: tuple[Fraction, Fraction], b: tuple[Fraction, Fraction]) -> Fraction:
    dq = b[0] - a[0]
    dr = b[1] - a[1]
    return dq * dq + dq * dr + dr * dr


def hex_patch(radius: int) -> tuple[tuple[Fraction, Fraction], ...]:
    pts = []
    for q in range(-radius, radius + 1):
        for r in range(-radius, radius + 1):
            if max(abs(q), abs(r), abs(q + r)) <= radius:
                pts.append((Fraction(q), Fraction(r)))
    return tuple(sorted(pts, key=lambda p: (sqdist_axial((0, 0), p), p[0], p[1])))


def distance_states(pset: MetricPointSet) -> dict[tuple[int, int], str]:
    sqdist = sqdist_cart if pset.metric == "cart" else sqdist_axial
    states = {}
    for i, j in combinations(range(len(pset.points)), 2):
        d2 = sqdist(pset.points[i], pset.points[j])
        if d2 < 1:
            state = "S"
        elif d2 == 1:
            state = "U"
        else:
            state = "L"
        states[(i, j)] = state
    return states


def unit_graph_from_states(n: int, states: dict[tuple[int, int], str]) -> tuple[int, ...]:
    adj = [0] * n
    for (i, j), state in states.items():
        if state == "U":
            adj[i] |= 1 << j
            adj[j] |= 1 << i
    return tuple(adj)


def find_undirected_hamiltonian_path(adj: tuple[int, ...], max_nodes: int = 1_000_000) -> tuple[int, ...] | None:
    n = len(adj)
    degrees = [bitcount(row) for row in adj]
    if any(d == 0 for d in degrees) and n > 1:
        return None
    calls = 0

    def ordered_neighbors(v: int, used: int) -> list[int]:
        out = []
        x = adj[v] & ~used
        while x:
            bit = x & -x
            w = bit.bit_length() - 1
            x -= bit
            onward = bitcount(adj[w] & ~used)
            out.append((onward, degrees[w], w))
        return [w for _, __, w in sorted(out)]

    starts = sorted(range(n), key=lambda v: (degrees[v], v))
    for start in starts:
        path = [start]
        used = 1 << start

        def dfs(v: int, used_mask: int) -> bool:
            nonlocal calls
            calls += 1
            if calls > max_nodes:
                return False
            if len(path) == n:
                return True
            for w in ordered_neighbors(v, used_mask):
                path.append(w)
                if dfs(w, used_mask | (1 << w)):
                    return True
                path.pop()
            return False

        if dfs(start, used):
            return tuple(path)
    return None


def reorder_states(states: dict[tuple[int, int], str], order: tuple[int, ...]) -> dict[tuple[int, int], str]:
    pos = {old: new for new, old in enumerate(order)}
    out = {}
    for (i, j), state in states.items():
        a = pos[i]
        b = pos[j]
        if a < b:
            out[(a, b)] = state
        else:
            out[(b, a)] = state
    return out


def trienerment_tournament(n: int, states: dict[tuple[int, int], str], unit_positive: bool) -> tuple[int, ...]:
    adj = [0] * n
    for i, j in combinations(range(n), 2):
        state = states[(i, j)]
        if state == "S":
            adj[i] |= 1 << j
        elif state == "L":
            adj[j] |= 1 << i
        else:
            if unit_positive:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
    return tuple(adj)


def path_state_counts(path: tuple[int, ...], states: dict[tuple[int, int], str]) -> Counter[str]:
    counts: Counter[str] = Counter()
    for a, b in zip(path, path[1:]):
        key = (a, b) if a < b else (b, a)
        counts[states[key]] += 1
    return counts


def unit_distance_probe() -> list[str]:
    square_center = MetricPointSet(
        "square_center_short_unit_long",
        (
            (Fraction(0), Fraction(0)),
            (Fraction(1), Fraction(0)),
            (Fraction(1), Fraction(1)),
            (Fraction(0), Fraction(1)),
            (Fraction(1, 2), Fraction(1, 2)),
        ),
        "cart",
    )
    hex7 = MetricPointSet("eisenstein_hex7", hex_patch(1), "axial")
    hex19 = MetricPointSet("eisenstein_hex19", hex_patch(2), "axial")
    hex21_points = tuple(list(hex_patch(2)) + [(Fraction(3), Fraction(0)), (Fraction(0), Fraction(3))])
    hex21 = MetricPointSet("eisenstein_hex19_plus_two_fringe_n21", hex21_points, "axial")
    point_sets = (square_center, hex7, hex19, hex21)

    lines = []
    lines.append("UNIT-DISTANCE SHORT/UNIT/LONG TRIENERMENT PROBE")
    lines.append("Rule: S(<1) follows base order, L(>1) reverses base order, U(=1) is a tie layer.")
    lines.append("")

    for pset in point_sets:
        raw_states = distance_states(pset)
        raw_unit_graph = unit_graph_from_states(len(pset.points), raw_states)
        unit_spine_raw = find_undirected_hamiltonian_path(raw_unit_graph)
        if unit_spine_raw is None:
            order = tuple(range(len(pset.points)))
            ordered_states = raw_states
            order_note = "no_unit_spine_found"
        else:
            order = unit_spine_raw
            ordered_states = reorder_states(raw_states, order)
            order_note = "base_order_is_unit_spine"
        state_counts = Counter(ordered_states.values())
        base_counts = path_state_counts(tuple(range(len(pset.points))), ordered_states)
        lines.append(
            f"{pset.name}: n={len(pset.points)}, order={order_note}, "
            + f"pair_states={fmt_counter(state_counts)}, base_path_states={fmt_counter(base_counts)}"
        )
        if unit_spine_raw is not None:
            lines.append(f"  unit spine in original labels: {unit_spine_raw}")

        for unit_positive in (True, False):
            adj = trienerment_tournament(len(pset.points), ordered_states, unit_positive)
            fp = tournament_fingerprint(adj)
            redei = redei_insert_path(adj)
            redei_counts = path_state_counts(redei, ordered_states)
            label = "unit_positive" if unit_positive else "unit_negative"
            extra = f", H={fp['H']}" if "H" in fp else ""
            lines.append(
                f"  {label}: score_hist={fp['score_hist']}, c3={fp['c3']}, "
                + f"scc={fp['scc']}{extra}, redei_path_states={fmt_counter(redei_counts)}"
            )
        lines.append("")

    lines.append("Unit-distance reading: the equality layer behaves like a wall, not like")
    lines.append("ordinary nonunit mass.  On Eisenstein patches, the U layer can carry a")
    lines.append("unit spine all the way to the n=21 fringe toy.  In the square-center toy,")
    lines.append("short edges isolate the center from the U graph, so any mandatory path")
    lines.append("must pay a short/long address.  This is the small-size impairment to use")
    lines.append("before asking whether non-lattice optimal rows keep a unit spine.")
    lines.append("")
    return lines


def route_tournament(lines: list[str]) -> None:
    routes = [
        ("signed_observer_gauge", (5, 5, 4, 5, 5)),
        ("alternating_pair_sum_clock", (4, 5, 5, 5, 4)),
        ("short_unit_long_trienerment", (5, 5, 5, 4, 5)),
        ("unit_spine_order", (4, 4, 5, 4, 5)),
        ("owner_address_derivative", (4, 5, 4, 5, 5)),
        ("raw_phase_or_edge_count", (1, 1, 2, 1, 4)),
    ]
    n = len(routes)
    adj = [0] * n
    for i, j in combinations(range(n), 2):
        vi = routes[i][1]
        vj = routes[j][1]
        wins_i = sum(1 for a, b in zip(vi, vj) if a > b)
        wins_j = sum(1 for a, b in zip(vi, vj) if b > a)
        if wins_i >= wins_j:
            adj[i] |= 1 << j
        else:
            adj[j] |= 1 << i
    fp = tournament_fingerprint(tuple(adj))
    order = sorted(range(n), key=lambda i: bitcount(adj[i]), reverse=True)
    lines.append("CROSS-LENS TOURNAMENT ANALYSIS")
    lines.append("Observable vector: predicate preservation, side-channel exposure, exact evidence, transfer to LRC14/UD21, overclaim safety.")
    lines.append(
        "fingerprint: "
        + f"score_hist={fp['score_hist']}, c3={fp['c3']}, scc={fp['scc']}, H={fp['H']}"
    )
    lines.append("outscore order:")
    for idx in order:
        lines.append(f"  {routes[idx][0]} -> {bitcount(adj[idx])}")
    lines.append("")


def main() -> None:
    lines: list[str] = []
    lines.extend(lrc_signed_probe())
    lines.extend(unit_distance_probe())
    route_tournament(lines)
    lines.append("SESSION THESIS")
    lines.append("1. Negative speeds do not change LRC loneliness; they change the pair-clock")
    lines.append("   address.  Alternating signs are a controlled way to inject the pair-sum")
    lines.append("   sieve into a runner-runner tournament while leaving the observer shadow")
    lines.append("   identical.")
    lines.append("2. Unit distance should keep S/U/L as a trienerment object.  Collapsing U")
    lines.append("   into either side erases exactly the wall that decides whether a unit")
    lines.append("   Hamiltonian spine is geometric or just a tiling-order artifact.")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
