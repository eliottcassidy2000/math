#!/usr/bin/env python3
"""Exact replay for the unbounded full-active-packet holonomy family.

For H >= 2 put F = 49H+1 and w_j = F-7j, 0 <= j <= 7.  On

    I_H = [x_(6H), x_(7H-1)],   x_m = (m+1/2)/F,

the full global wall word is

    x_m, y_(7,m), y_(6,m), ..., y_(1,m), x_(m+1),

for m = 6H,...,7H-2, where y_(j,m) is the w_j-wall of index
m+1-j.  Every wall is covered, so I_H contains H-1 consecutive active
fastest periods, each carrying the same full seven-owner visitor packet.

This verifier uses only integer arithmetic and fractions.Fraction.  It:

* independently enumerates every wall in I_H and compares the complete list
  with the displayed word for every H=2,...,200;
* checks all wall rainbows and every intervening chamber exactly;
* checks the wall-count, active-period, extent, and ratio formulas;
* checks that H>=5 refutes L < 1/g + 2/f; and
* certifies the elementary marked-cut and circular R_7 tournament
  fingerprints, including edge flips and Hamiltonian-path count.
"""

from fractions import Fraction as Q
from itertools import combinations, permutations


P = 7
H_MIN = 2
H_MAX = 200


def inv7(w):
    return pow(w % P, -1, P)


def floor_q(x):
    return x.numerator // x.denominator


def ceil_q(x):
    return -((-x.numerator) // x.denominator)


def nearest_integer_or_none(x):
    """Nearest integer to x, or None at an exact half-integer."""
    q, r = divmod(x.numerator, x.denominator)
    twice = 2 * r
    if twice == x.denominator:
        return None
    return q + (twice > x.denominator)


def token(w, x):
    n = nearest_integer_or_none(w * x)
    if n is None:
        return None
    return (-inv7(w) * n) % P


def family(H):
    F = 49 * H + 1
    return F, tuple(F - 7 * j for j in range(8))


def f_wall(F, m):
    return Q(2 * m + 1, 2 * F)


def visitor_wall(F, j, m):
    w = F - 7 * j
    n = m + 1 - j
    return Q(2 * n + 1, 2 * w)


def expected_events(H):
    """Closed-form full event word on I_H.

    Event tuples are (coordinate, owner index j, speed, local wall index).
    """
    F, W = family(H)
    events = [(f_wall(F, 6 * H), 0, F, 6 * H)]
    for m in range(6 * H, 7 * H - 1):
        for j in range(7, 0, -1):
            n = m + 1 - j
            events.append((visitor_wall(F, j, m), j, W[j], n))
        events.append((f_wall(F, m + 1), 0, F, m + 1))
    return events


def enumerated_events(H):
    """Independently enumerate every owner wall in the closed interval I_H."""
    F, W = family(H)
    lo = f_wall(F, 6 * H)
    hi = f_wall(F, 7 * H - 1)
    events = []
    for j, w in enumerate(W):
        # (n+1/2)/w lies in [lo,hi] iff n lies in this exact interval.
        n_lo = max(0, ceil_q(w * lo - Q(1, 2)))
        n_hi = min(w - 1, floor_q(w * hi - Q(1, 2)))
        for n in range(n_lo, n_hi + 1):
            x = Q(2 * n + 1, 2 * w)
            assert lo <= x <= hi
            events.append((x, j, w, n))
    events.sort()
    return events


def wall_is_covered(W, event):
    x, owner, _, _ = event
    assert token(W[owner], x) is None
    other_tokens = [token(w, x) for j, w in enumerate(W) if j != owner]
    return None not in other_tokens and sorted(other_tokens) == list(range(P))


def chamber_is_covered(W, left, right):
    x = (left + right) / 2
    tokens = [token(w, x) for w in W]
    return None not in tokens and set(tokens) == set(range(P))


def audit_family(H):
    F, W = family(H)
    g = W[1]
    expected = expected_events(H)
    actual = enumerated_events(H)

    assert W == tuple(sorted(W, reverse=True))
    assert len(set(W)) == 8
    assert min(W) > 0
    assert all(w % P == 1 for w in W)

    # This simultaneously certifies that there are no missing events or ties.
    assert actual == expected
    assert all(actual[i][0] < actual[i + 1][0] for i in range(len(actual) - 1))
    assert len(actual) == 8 * H - 7

    assert all(wall_is_covered(W, event) for event in actual)
    assert all(
        chamber_is_covered(W, actual[i][0], actual[i + 1][0])
        for i in range(len(actual) - 1)
    )

    f_positions = [i for i, event in enumerate(actual) if event[1] == 0]
    assert len(f_positions) == H
    assert [actual[i][3] for i in f_positions] == list(range(6 * H, 7 * H))

    packets = []
    for left, right in zip(f_positions, f_positions[1:]):
        packet = tuple(actual[i][1] for i in range(left + 1, right))
        packets.append(packet)
    assert packets == [(7, 6, 5, 4, 3, 2, 1)] * (H - 1)
    assert len(packets) == H - 1
    assert all(sum(inv7(W[j]) for j in packet) % P == 0 for packet in packets)

    extent = actual[-1][0] - actual[0][0]
    assert extent == Q(H - 1, F)
    assert (F + g - 1) // g == 2  # ceil(f/g)=2
    assert F - g == 7

    proposed_bound = Q(1, g) + Q(2, F)
    if H >= 5:
        assert extent > proposed_bound

    return {
        "events": len(actual),
        "chambers": len(actual) - 1,
        "active_periods": len(packets),
        "extent": extent,
        "bound": proposed_bound,
    }


def tournament_edge(adj, a, b):
    assert a != b
    return adj[a][b]


def score_sequence(adj):
    return tuple(sum(row) for row in adj)


def directed_triangle_count(adj):
    count = 0
    for triple in combinations(range(len(adj)), 3):
        scores = [sum(adj[a][b] for b in triple if b != a) for a in triple]
        count += sorted(scores) == [1, 1, 1]
    return count


def reachable(adj, source, target):
    seen = {source}
    stack = [source]
    while stack:
        a = stack.pop()
        for b, edge in enumerate(adj[a]):
            if edge and b not in seen:
                seen.add(b)
                stack.append(b)
    return target in seen


def scc_sizes(adj):
    remaining = set(range(len(adj)))
    sizes = []
    while remaining:
        a = min(remaining)
        component = {b for b in remaining if reachable(adj, a, b) and reachable(adj, b, a)}
        sizes.append(len(component))
        remaining -= component
    return tuple(sorted(sizes, reverse=True))


def hamiltonian_paths(adj):
    return tuple(
        path
        for path in permutations(range(len(adj)))
        if all(tournament_edge(adj, path[i], path[i + 1]) for i in range(len(path) - 1))
    )


def marked_cut_tournament(phase):
    # Visitor owner a=0,...,6 corresponds to j=a+1 and has token a-phase.
    tokens = tuple((a - phase) % P for a in range(P))
    return tuple(
        tuple(a != b and tokens[a] < tokens[b] for b in range(P))
        for a in range(P)
    )


def circular_tournament(phase):
    tokens = tuple((a - phase) % P for a in range(P))
    return tuple(
        tuple(a != b and (tokens[b] - tokens[a]) % P in (1, 2, 3) for b in range(P))
        for a in range(P)
    )


def edge_flip_count(left, right):
    return sum(
        left[a][b] != right[a][b]
        for a in range(len(left))
        for b in range(a + 1, len(left))
    )


def audit_tournaments():
    marked = tuple(marked_cut_tournament(phase) for phase in range(P))
    circular = tuple(circular_tournament(phase) for phase in range(P))

    assert all(tuple(sorted(score_sequence(T))) == tuple(range(P)) for T in marked)
    assert all(directed_triangle_count(T) == 0 for T in marked)
    assert all(scc_sizes(T) == (1, 1, 1, 1, 1, 1, 1) for T in marked)
    marked_hps = tuple(hamiltonian_paths(T) for T in marked)
    assert all(len(paths) == 1 for paths in marked_hps)

    # Global sheet translation changes the marked cut but not the circular gauge.
    assert all(edge_flip_count(marked[p], marked[(p + 1) % P]) == 6 for p in range(P))
    assert all(T == circular[0] for T in circular)
    assert all(edge_flip_count(marked[p], circular[p]) == 6 for p in range(P))

    R7 = circular[0]
    R7_paths = hamiltonian_paths(R7)
    assert score_sequence(R7) == (3, 3, 3, 3, 3, 3, 3)
    assert directed_triangle_count(R7) == 14
    assert scc_sizes(R7) == (7,)
    assert len(R7_paths) == 175
    assert min(R7_paths) == (0, 1, 2, 3, 4, 5, 6)

    # The chronological packet 6,5,...,0 is a Hamiltonian path in the converse
    # time-compatible circular gauge.
    chronological = tuple(range(6, -1, -1))
    assert all(R7[chronological[i + 1]][chronological[i]] for i in range(6))

    return {
        "marked_scores": tuple(sorted(score_sequence(marked[0]))),
        "marked_triangles": directed_triangle_count(marked[0]),
        "marked_sccs": scc_sizes(marked[0]),
        "marked_hamiltonian_paths": len(marked_hps[0]),
        "translation_edge_flips": edge_flip_count(marked[0], marked[1]),
        "gauge_edge_flips": edge_flip_count(marked[0], R7),
        "circular_scores": score_sequence(R7),
        "circular_triangles": directed_triangle_count(R7),
        "circular_sccs": scc_sizes(R7),
        "circular_hamiltonian_paths": len(R7_paths),
        "circular_tie_path": min(R7_paths),
        "chronological_converse_path": chronological,
    }


def main():
    rows = []
    for H in range(H_MIN, H_MAX + 1):
        rows.append((H, audit_family(H)))
    tournament = audit_tournaments()

    total_events = sum(row["events"] for _, row in rows)
    total_chambers = sum(row["chambers"] for _, row in rows)
    total_active = sum(row["active_periods"] for _, row in rows)
    refutations = sum(row["extent"] > row["bound"] for H, row in rows if H >= 5)

    print("FULL ACTIVE-PACKET HOLONOMY -- EXACT REPLAY")
    print("family: F=49H+1; w_j=F-7j (0<=j<=7); all w_j=1 (mod 7)")
    print("interval: [x_(6H),x_(7H-1)], x_m=(m+1/2)/F")
    print("packet word per active period: (7,6,5,4,3,2,1)")
    print()
    print(f"audited H range: {H_MIN}..{H_MAX} ({len(rows)} families)")
    print(f"complete event-list comparisons: {len(rows)} passed")
    print(f"covered walls checked: {total_events}")
    print(f"covered chambers checked: {total_chambers}")
    print(f"active periods checked: {total_active}")
    print("per-family formulas: walls=8H-7; active=H-1; extent=(H-1)/(49H+1)")
    print("ratio/gap formulas: ceil(f/g)=2; f-g=7")
    print(f"extent refutations for H>=5: {refutations}/{H_MAX - 4}")
    print()
    for H in (2, 5, 10, 50, 200):
        F, W = family(H)
        row = dict(rows)[H]
        print(
            f"sample H={H}: W={tuple(reversed(W))}; walls={row['events']}; "
            f"active={row['active_periods']}; extent={row['extent']}; "
            f"proposed_bound={row['bound']}; refutes={row['extent'] > row['bound']}"
        )
    print()
    print("TOURNAMENT FINGERPRINTS")
    print(f"marked-cut scores: {tournament['marked_scores']}")
    print(f"marked-cut directed triangles: {tournament['marked_triangles']}")
    print(f"marked-cut SCC sizes: {tournament['marked_sccs']}")
    print(f"marked-cut Hamiltonian paths: {tournament['marked_hamiltonian_paths']}")
    print(f"marked-cut edge flips per global translation: {tournament['translation_edge_flips']}")
    print(f"marked-cut/circular gauge edge flips: {tournament['gauge_edge_flips']}")
    print(f"circular R7 scores: {tournament['circular_scores']}")
    print(f"circular R7 directed triangles: {tournament['circular_triangles']}")
    print(f"circular R7 SCC sizes: {tournament['circular_sccs']}")
    print(f"circular R7 Hamiltonian paths: {tournament['circular_hamiltonian_paths']}")
    print(f"circular R7 lexicographic tie path: {tournament['circular_tie_path']}")
    print(f"chronological packet path in converse gauge: {tournament['chronological_converse_path']}")
    print()
    print("SUMMARY: ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
