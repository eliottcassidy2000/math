#!/usr/bin/env python3
"""Classify pair-good decoys in the LRC14 binding-pair carrier.

HYP-3019 counted pair-good decoy times: THM-524 candidate times where at least
one tied pair is at distance >= 1/14, but another runner makes the true minimum
< 1/14.  This script classifies their exact generators.

For a pair source (a,b), lane D=a+b or |a-b|, and candidate t=k/D, a blocker c
kills the pair-good time exactly when

    min((c*k mod D), D-(c*k mod D)) / D < 1/14.

After reducing t=p/q, this is the residue-tooth inequality

    14 * min((c*p mod q), q-(c*p mod q)) < q.

The point is to replace raw decoy counts by a small set of modular tooth
generators: source pair lane, blocker role, zero/near tooth, and active blocker
deck.  Tournament Analysis uses generator/proof carriers as vertices, not
runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from math import gcd


N = 14
THRESHOLD = F(1, N)
HALF = F(1, 2)


def fmt(x: F) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def dist(v: int, t: F) -> F:
    r = (v * t) % 1
    return min(r, 1 - r)


def candidate_times(row: tuple[int, ...]) -> tuple[F, ...]:
    times: set[F] = {HALF}
    for a, b in combinations(row, 2):
        for d in {a + b, abs(a - b)}:
            if d <= 0:
                continue
            for k in range(1, d):
                t = F(k, d)
                if t <= HALF:
                    times.add(t)
    return tuple(sorted(times))


def primitive(row: tuple[int, ...]) -> bool:
    g = 0
    for v in row:
        g = gcd(g, v)
    return g == 1


def residue_depth(v: int, t: F) -> tuple[int, int, int]:
    """Return (depth numerator, reduced denominator, signed residue)."""
    p, q = t.numerator, t.denominator
    r = (v * p) % q
    return min(r, q - r), q, r


def role(row: tuple[int, ...], v: int) -> str:
    if v <= 13:
        return "core"
    if v % N == 0:
        return "tail_multiple14"
    if v % N == 13:
        return "tail_minus1"
    return "tail_other"


def lane_modes(a: int, b: int, t: F) -> tuple[str, ...]:
    modes: list[str] = []
    if ((a + b) * t).denominator == 1:
        modes.append("sum")
    d = abs(a - b)
    if d and (d * t).denominator == 1:
        modes.append("diff")
    if t == HALF and not modes:
        modes.append("halfturn")
    return tuple(modes or ("same_distance",))


def shell_label(a: int, b: int) -> str:
    s = (a + b) % N
    d = abs(a - b)
    flags: list[str] = []
    if s == 0:
        flags.append("sum0")
    if a + b == N:
        flags.append("literal14")
    if gcd(a + b, N) > 1:
        flags.append(f"sum_g{gcd(a + b, N)}")
    if d and gcd(d, N) > 1:
        flags.append(f"diff_g{gcd(d, N)}")
    return "+".join(flags) if flags else f"sum{s}"


@dataclass(frozen=True)
class DecoyTime:
    t: F
    actual_gap: F
    good_pairs: tuple[tuple[int, int, tuple[str, ...], str], ...]
    blockers: tuple[int, ...]
    blocker_depths: tuple[tuple[int, int, int, int], ...]

    @property
    def zero_blocked(self) -> bool:
        return any(depth == 0 for _v, depth, _q, _r in self.blocker_depths)

    @property
    def blocker_roles(self) -> tuple[str, ...]:
        return tuple(sorted({role((), v) for v in self.blockers}))


@dataclass(frozen=True)
class RowReport:
    name: str
    family: str
    row: tuple[int, ...]
    candidate_count: int
    decoys: tuple[DecoyTime, ...]
    witnesses: int
    boundary: int


def decoy_times(row: tuple[int, ...]) -> tuple[DecoyTime, ...]:
    out: list[DecoyTime] = []
    for t in candidate_times(row):
        distances = {v: dist(v, t) for v in row}
        actual = min(distances.values())
        good_pairs: list[tuple[int, int, tuple[str, ...], str]] = []
        for a, b in combinations(row, 2):
            if distances[a] == distances[b] and distances[a] >= THRESHOLD:
                good_pairs.append((a, b, lane_modes(a, b, t), shell_label(a, b)))
        if not good_pairs or actual >= THRESHOLD:
            continue
        blockers = tuple(v for v in row if distances[v] == actual)
        depths = tuple((v, *residue_depth(v, t)) for v in blockers)
        out.append(
            DecoyTime(
                t=t,
                actual_gap=actual,
                good_pairs=tuple(good_pairs),
                blockers=blockers,
                blocker_depths=depths,
            )
        )
    return tuple(out)


def named_rows() -> tuple[tuple[str, tuple[int, ...], str], ...]:
    return (
        ("AP13", tuple(range(1, 14)), "AP/GW boundary"),
        ("GW 12->24", tuple(list(range(1, 12)) + [13, 24]), "AP/GW boundary"),
        ("petal 10->20", tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 20]), "unit petal"),
        ("petal 13->26", tuple(list(range(1, 13)) + [26]), "unit petal"),
        ("K33 12->36", tuple(list(range(1, 12)) + [13, 36]), "K33/state-lift"),
        ("K33 family drop12,13 add26,36", tuple(list(range(1, 12)) + [26, 36]), "K33/state-lift"),
        ("covering 12->84", tuple(list(range(1, 12)) + [13, 84]), "covering"),
        ("single far tail 12->200", tuple(list(range(1, 13)) + [200]), "covering"),
        ("drop6 core add180", tuple([1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 180]), "drop-core"),
        ("two-far doublet 20,21", tuple(list(range(1, 12)) + [20, 21]), "genuine-wide doublet"),
    )


def row_report(name: str, row: tuple[int, ...], family: str) -> RowReport:
    assert len(row) == 13 and primitive(row)
    times = candidate_times(row)
    decoys = decoy_times(row)
    witnesses = 0
    boundary = 0
    for t in times:
        g = min(dist(v, t) for v in row)
        if g >= THRESHOLD:
            witnesses += 1
        if g == THRESHOLD:
            boundary += 1
    return RowReport(name, family, row, len(times), decoys, witnesses, boundary)


def role_signature(row: tuple[int, ...], d: DecoyTime) -> str:
    roles = sorted(role(row, v) for v in d.blockers)
    return "+".join(roles)


def depth_bucket(d: DecoyTime) -> str:
    if d.zero_blocked:
        return "zero_tooth"
    min_depth, q = min((depth, q) for _v, depth, q, _r in d.blocker_depths)
    clearance = q - N * min_depth
    if clearance == 1:
        return "near_wall_slack1"
    if clearance <= 3:
        return "near_wall_slack2_3"
    return "deep_tooth"


def source_lane_bucket(d: DecoyTime) -> str:
    modes = sorted({mode for _a, _b, modes, _shell in d.good_pairs for mode in modes})
    if d.t == HALF:
        return "halfturn"
    return "+".join(modes)


def generator_keys(row: tuple[int, ...], d: DecoyTime) -> set[tuple[str, str, str, str]]:
    keys: set[tuple[str, str, str, str]] = set()
    blocker_sig = role_signature(row, d)
    bucket = depth_bucket(d)
    for a, b, modes, shell in d.good_pairs:
        for mode in modes:
            keys.add((mode, shell, blocker_sig, bucket))
    return keys


def compressed_cover(row: tuple[int, ...], decoys: tuple[DecoyTime, ...]) -> list[tuple[tuple[str, str, str, str], int]]:
    uncovered = set(range(len(decoys)))
    key_to_ids: dict[tuple[str, str, str, str], set[int]] = defaultdict(set)
    for idx, d in enumerate(decoys):
        for key in generator_keys(row, d):
            key_to_ids[key].add(idx)
    chosen: list[tuple[tuple[str, str, str, str], int]] = []
    while uncovered:
        key, ids = max(key_to_ids.items(), key=lambda item: len(item[1] & uncovered))
        gain = len(ids & uncovered)
        if gain == 0:
            break
        chosen.append((key, gain))
        uncovered -= ids
    return chosen


@dataclass(frozen=True)
class Carrier:
    name: str
    score: tuple[int, ...]


def tournament_fingerprint() -> dict[str, object]:
    carriers = (
        Carrier("blocker_residue_tooth_inequality", (6, 6, 6, 6, 5, 6)),
        Carrier("source_pair_lane_plus_blocker_deck", (6, 6, 5, 5, 6, 5)),
        Carrier("normal_fan_active_owner_support", (6, 5, 5, 6, 5, 5)),
        Carrier("family_compressed_generator_cover", (5, 6, 4, 5, 5, 6)),
        Carrier("zero_tooth_divisibility_subcase", (4, 5, 6, 3, 4, 6)),
        Carrier("near_wall_clearance_margin", (5, 5, 4, 5, 4, 5)),
        Carrier("raw_decoy_count", (1, 1, 1, 1, 1, 2)),
        Carrier("raw_pair_gap", (1, 0, 2, 1, 0, 1)),
    )
    n = len(carriers)
    out = [0] * n
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        si, sj = carriers[i].score, carriers[j].score
        wi = sum(a > b for a, b in zip(si, sj))
        wj = sum(b > a for a, b in zip(si, sj))
        if wi > wj or (wi == wj and i < j):
            a, b = i, j
        else:
            a, b = j, i
        adj[a][b] = True
        out[a] += 1
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or (adj[i][k] and adj[k][j] and adj[j][i]):
            c3 += 1
    # Hamiltonian path DP.
    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            ways = dp.get((mask, last), 0)
            if not ways:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + ways
    full = (1 << n) - 1
    hp = sum(dp.get((full, i), 0) for i in range(n))
    order = sorted(range(n), key=lambda i: carriers[i].score, reverse=True)
    return {
        "score_hist": sorted(Counter(out).items()),
        "directed_3cycles": c3,
        "hamiltonian_path_count": hp,
        "tie_path": " > ".join(carriers[i].name for i in order),
    }


def print_assumption_challenge() -> None:
    print("Assumption challenge")
    print("  Candidate vertex sets considered:")
    print("    runners, gaps, fixed sections, section boundaries, wall crossings,")
    print("    residues, cover arcs, Fourier modes, matroid circuits, proof obligations,")
    print("    active normal-fan owners, and decoy generator triples.")
    print("  Chosen vertex set:")
    print("    generator triples (source pair lane, blocker residue tooth, blocker deck).")
    print("  Preserved LRC predicate:")
    print("    pair-good candidate is harmless exactly when every blocker residue tooth")
    print("    is outside threshold; decoy means at least one tooth inequality is strict.")
    print("  Destroyed by raw counts:")
    print("    which blocker tooth generated the failure and whether many times share one mechanism.")
    print()


def main() -> None:
    print("LRC14 pair-good decoy classifier")
    print("=" * 72)
    print_assumption_challenge()
    print("[0] Exact generator lemma")
    print("  For t=p/q in lowest terms, runner c blocks iff")
    print("    14 * min(c*p mod q, q-(c*p mod q)) < q.")
    print("  A pair-good decoy is exactly a source pair lane plus at least one such")
    print("  blocker tooth.  Counts are therefore shadows of modular tooth families.")
    print()

    reports = tuple(row_report(name, row, family) for name, row, family in named_rows())
    print("[1] Row-level decoy compression")
    for rep in reports:
        decoys = rep.decoys
        zero = sum(1 for d in decoys if d.zero_blocked)
        role_hist = Counter(role_signature(rep.row, d) for d in decoys)
        depth_hist = Counter(depth_bucket(d) for d in decoys)
        lane_hist = Counter(source_lane_bucket(d) for d in decoys)
        cover = compressed_cover(rep.row, decoys)[:6]
        print(f"  {rep.name}")
        print(
            f"    candidates={rep.candidate_count}; witnesses={rep.witnesses}; "
            f"boundary={rep.boundary}; decoys={len(decoys)}; zero_blocked={zero}"
        )
        print(f"    blocker_roles={dict(role_hist)}")
        print(f"    depth_buckets={dict(depth_hist)}")
        print(f"    source_lanes={dict(lane_hist)}")
        print("    greedy generator cover:")
        for key, gain in cover:
            mode, shell, blocker_sig, bucket = key
            print(f"      +{gain:4d}  mode={mode:<4s} shell={shell:<18s} blockers={blocker_sig:<20s} tooth={bucket}")
    print()

    print("[2] Global generator census over named rows")
    global_keys = Counter()
    role_keys = Counter()
    tooth_keys = Counter()
    for rep in reports:
        for d in rep.decoys:
            role_keys[role_signature(rep.row, d)] += 1
            tooth_keys[depth_bucket(d)] += 1
            for key in generator_keys(rep.row, d):
                global_keys[key] += 1
    print(f"  total_decoy_times={sum(len(r.decoys) for r in reports)}")
    print(f"  distinct_generator_keys={len(global_keys)}")
    print(f"  blocker_role_hist={dict(role_keys)}")
    print(f"  tooth_bucket_hist={dict(tooth_keys)}")
    print("  top generator keys:")
    for (mode, shell, blocker_sig, bucket), count in global_keys.most_common(12):
        print(f"    {count:5d}  mode={mode:<4s} shell={shell:<18s} blockers={blocker_sig:<20s} tooth={bucket}")
    print()

    print("[3] Why the count is not the invariant")
    print("  The largest rows have thousands of decoys because many source pairs cross")
    print("  inside the same small set of blocker tooth regimes.  The proof object is")
    print("  the blocker residue inequality, not the number of candidate times hitting it.")
    print("  Packet sidecar fields suggested:")
    print("    decoy_source_lane, decoy_source_shell, decoy_blocker_role,")
    print("    decoy_blocker_depth_bucket, decoy_zero_tooth_flag,")
    print("    decoy_generator_cover_size, decoy_top_generator_key.")
    print()

    print("[4] Tournament Analysis")
    fp = tournament_fingerprint()
    print("  vertices are generator/proof carriers, not runners.")
    print("  pairwise observable:")
    print("    predicate exactness, modular checkability, denominator retention,")
    print("    normal-fan compatibility, family compression, and anti-count guard.")
    print("  switch:")
    print("    majority comparison of these coordinates; ties use the listed path.")
    print(f"  score_hist={fp['score_hist']}")
    print(f"  directed_3cycles={fp['directed_3cycles']}")
    print(f"  hamiltonian_path_count={fp['hamiltonian_path_count']}")
    print(f"  tie_hamiltonian_path={fp['tie_path']}")


if __name__ == "__main__":
    main()
