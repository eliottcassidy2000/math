#!/usr/bin/env python3
"""Exact tooth-label winding atlas for the twelve-speed LRC frontier.

This is a computational scout, not a uniform proof of the n=12 sporadic
branch.  It studies completions ``A=P union {w}`` at the threshold L=1/13.
Every geometric quantity below is a ``fractions.Fraction``.

The main bank is a transparent exhaustive max-peel slice: every primitive
12-subset A of {1,...,20} for which

    w=max(A),  M(A\\{w})>1/12,  and  M(A)<=1/10.

The inequalities are selected exactly by maximal masks on the cells cut out
by all threshold endpoints.  Eleven eligible deletions of {1,...,12} are
added as the tight/all-covered controls.  For every row the script computes
all open components of E_{1/13}(P), their endpoint owners, exact tooth slack

    sigma(J,w) = 1/13 - ||w c_J|| - w h_J,

and the cyclic word of nearest w-tooth labels.

Tournament Analysis uses safe components, rather than runners, as vertices.
For components i,j the pairwise observable is the signed circular displacement
of the within-tooth phases theta_i={w c_i}.  Clockwise half-circle orientation
is the switch/gauge.  Chronological component order 0->1->... is the fixed tie
Hamiltonian path for coincident or antipodal phases.  We report score
histograms, directed 3-cycles, SCCs, edge flips from the tie path, and exact
Hamiltonian-path counts when the component count is at most 12.

The SHA-256 digest commits to every exact row and component, so the concise
output still audits the full bank.  Replay with

    python3 04-computation/lrc13_n12_tooth_winding_atlas_codex_S3.py
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import comb, gcd


L = F(1, 13)
CORE_FLOOR = F(1, 12)
NEAR_CEILING = F(1, 10)
N = 20
HP_COMPONENT_CAP = 12


def ftext(x: F) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def circle_distance(x: F) -> F:
    r = x % 1
    return min(r, 1 - r)


def set_mask(values: tuple[int, ...]) -> int:
    return sum(1 << (v - 1) for v in values)


def maximal_strict_safe_masks(n: int, level: F) -> tuple[int, ...]:
    """Maximal speed masks simultaneously strict-safe somewhere at ``level``.

    All predicates ||v t||>level are constant on the open cells between their
    rational endpoints.  Sampling the exact midpoint of every global cell is
    therefore lossless.  A set S has M(S)>level iff its mask is contained in
    one of the returned maximal masks.
    """

    endpoints = {F(0), F(1)}
    for v in range(1, n + 1):
        for k in range(v):
            endpoints.add((F(k, v) - level / v) % 1)
            endpoints.add((F(k, v) + level / v) % 1)
    points = sorted(endpoints)
    raw: set[int] = set()
    for lo, hi in zip(points, points[1:]):
        if lo == hi:
            continue
        t = (lo + hi) / 2
        raw.add(
            sum(
                1 << (v - 1)
                for v in range(1, n + 1)
                if circle_distance(v * t) > level
            )
        )

    maximal: list[int] = []
    for mask in sorted(raw, key=lambda m: (-m.bit_count(), m)):
        if not any(mask & ~larger == 0 for larger in maximal):
            maximal.append(mask)
    return tuple(maximal)


def contained_in_some(mask: int, containers: tuple[int, ...]) -> bool:
    return any(mask & ~container == 0 for container in containers)


def intersect_interval_unions(
    left: list[tuple[F, F]], right: list[tuple[F, F]]
) -> list[tuple[F, F]]:
    out: list[tuple[F, F]] = []
    i = j = 0
    while i < len(left) and j < len(right):
        a, b = left[i]
        c, d = right[j]
        lo, hi = max(a, c), min(b, d)
        if hi > lo:
            out.append((lo, hi))
        if b <= d:
            i += 1
        else:
            j += 1
    return out


def strict_safe_components(speeds: tuple[int, ...], level: F = L) -> tuple[tuple[F, F], ...]:
    """The exact open components of {t in R/Z: min_v ||v t||>level}."""

    current = [(F(0), F(1))]
    for v in speeds:
        runner_safe = [
            ((k + level) / v, (k + 1 - level) / v) for k in range(v)
        ]
        current = intersect_interval_unions(current, runner_safe)
        if not current:
            break
    return tuple(current)


def exact_m(speeds: tuple[int, ...]) -> F:
    """Exact max_t min_v ||v t|| from all pair crossings and self-cusps."""

    denominators = {2 * v for v in speeds}
    for i, v in enumerate(speeds):
        for u in speeds[i + 1 :]:
            denominators.add(u + v)
            denominators.add(u - v)

    best_num, best_den = 0, 1
    for q in denominators:
        for p in range(1, q):
            clearance_num = min(
                min((v * p) % q, q - ((v * p) % q)) for v in speeds
            )
            if clearance_num * best_den > best_num * q:
                best_num, best_den = clearance_num, q
    return F(best_num, best_den)


def nearest_integer(x: F) -> int | None:
    floor = x.numerator // x.denominator
    remainder = x - floor
    if 2 * remainder < 1:
        return floor
    if 2 * remainder > 1:
        return floor + 1
    return None


@dataclass(frozen=True)
class Component:
    lo: F
    hi: F
    center: F
    half_width: F
    left_owners: tuple[int, ...]
    right_owners: tuple[int, ...]
    tooth_label: int | None
    tooth_phase: F
    sigma: F


@dataclass(frozen=True)
class TournamentFingerprint:
    score_histogram: tuple[tuple[int, int], ...]
    directed_triangles: int
    scc_sizes: tuple[int, ...]
    edge_flips: int
    phase_ties: int
    hamiltonian_paths: int | None


@dataclass(frozen=True)
class Row:
    family: str
    A: tuple[int, ...]
    P: tuple[int, ...]
    w: int
    core_m: F
    full_m: F
    components: tuple[Component, ...]
    labels: tuple[int, ...] | None
    increments: tuple[int, ...] | None
    winding: int | None
    tournament: TournamentFingerprint

    @property
    def all_covered(self) -> bool:
        return all(component.sigma >= 0 for component in self.components)

    @property
    def endpoint_pure(self) -> bool:
        return all(
            len(component.left_owners) == len(component.right_owners) == 1
            for component in self.components
        )

    @property
    def minimum_sigma(self) -> F:
        return min(component.sigma for component in self.components)

    @property
    def zero_slacks(self) -> int:
        return sum(component.sigma == 0 for component in self.components)

    @property
    def negative_slacks(self) -> int:
        return sum(component.sigma < 0 for component in self.components)


def endpoint_owners(speeds: tuple[int, ...], point: F) -> tuple[int, ...]:
    return tuple(v for v in speeds if circle_distance(v * point) == L)


def strongly_connected_component_sizes(out_bits: tuple[int, ...]) -> tuple[int, ...]:
    n = len(out_bits)
    seen = [False] * n
    finish: list[int] = []

    def dfs1(v: int) -> None:
        seen[v] = True
        bits = out_bits[v]
        while bits:
            bit = bits & -bits
            u = bit.bit_length() - 1
            bits ^= bit
            if not seen[u]:
                dfs1(u)
        finish.append(v)

    for v in range(n):
        if not seen[v]:
            dfs1(v)

    reverse = [0] * n
    for v, bits in enumerate(out_bits):
        while bits:
            bit = bits & -bits
            u = bit.bit_length() - 1
            bits ^= bit
            reverse[u] |= 1 << v

    seen = [False] * n
    sizes: list[int] = []

    def dfs2(v: int) -> int:
        seen[v] = True
        size = 1
        bits = reverse[v]
        while bits:
            bit = bits & -bits
            u = bit.bit_length() - 1
            bits ^= bit
            if not seen[u]:
                size += dfs2(u)
        return size

    for v in reversed(finish):
        if not seen[v]:
            sizes.append(dfs2(v))
    return tuple(sorted(sizes, reverse=True))


def hamiltonian_path_count(out_bits: tuple[int, ...]) -> int:
    """Exact directed Hamiltonian-path count by subset DP."""

    n = len(out_bits)
    states = 1 << n
    dp = [[0] * n for _ in range(states)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, states):
        if mask & (mask - 1) == 0:
            continue
        last_bits = mask
        while last_bits:
            last_bit = last_bits & -last_bits
            last = last_bit.bit_length() - 1
            last_bits ^= last_bit
            previous = mask ^ last_bit
            predecessor_bits = previous
            total = 0
            while predecessor_bits:
                predecessor_bit = predecessor_bits & -predecessor_bits
                predecessor = predecessor_bit.bit_length() - 1
                predecessor_bits ^= predecessor_bit
                if out_bits[predecessor] & last_bit:
                    total += dp[previous][predecessor]
            dp[mask][last] = total
    return sum(dp[-1])


def component_tournament(components: tuple[Component, ...]) -> TournamentFingerprint:
    n = len(components)
    out = [0] * n
    ties = 0
    flips = 0
    for i in range(n):
        for j in range(i + 1, n):
            delta = (components[j].tooth_phase - components[i].tooth_phase) % 1
            if delta == 0 or delta == F(1, 2):
                # Fixed tie Hamiltonian path: chronological i -> j.
                winner, loser = i, j
                ties += 1
            elif delta < F(1, 2):
                winner, loser = i, j
            else:
                winner, loser = j, i
                flips += 1
            out[winner] |= 1 << loser

    out_tuple = tuple(out)
    scores = tuple(bits.bit_count() for bits in out_tuple)
    score_histogram = tuple(sorted(Counter(scores).items()))
    directed_triangles = comb(n, 3) - sum(comb(score, 2) for score in scores)
    scc_sizes = strongly_connected_component_sizes(out_tuple)
    hp = hamiltonian_path_count(out_tuple) if n <= HP_COMPONENT_CAP else None
    return TournamentFingerprint(
        score_histogram=score_histogram,
        directed_triangles=directed_triangles,
        scc_sizes=scc_sizes,
        edge_flips=flips,
        phase_ties=ties,
        hamiltonian_paths=hp,
    )


def audit_row(family: str, A: tuple[int, ...], w: int) -> Row:
    P = tuple(v for v in A if v != w)
    if len(P) != 11:
        raise AssertionError("completion speed must occur exactly once")

    components: list[Component] = []
    for lo, hi in strict_safe_components(P):
        center = (lo + hi) / 2
        half_width = (hi - lo) / 2
        unwrapped_label = nearest_integer(w * center)
        tooth_label = None if unwrapped_label is None else unwrapped_label % w
        phase = (w * center) % 1
        sigma = L - circle_distance(w * center) - w * half_width
        components.append(
            Component(
                lo=lo,
                hi=hi,
                center=center,
                half_width=half_width,
                left_owners=endpoint_owners(P, lo),
                right_owners=endpoint_owners(P, hi),
                tooth_label=tooth_label,
                tooth_phase=phase,
                sigma=sigma,
            )
        )

    component_tuple = tuple(components)
    if any(component.tooth_label is None for component in component_tuple):
        labels = increments = None
        winding = None
    else:
        labels = tuple(int(component.tooth_label) for component in component_tuple)
        increments = tuple(
            (labels[(i + 1) % len(labels)] - labels[i]) % w
            for i in range(len(labels))
        )
        increment_sum = sum(increments)
        if increment_sum % w:
            raise AssertionError("cyclic label increments have nonintegral winding")
        winding = increment_sum // w

    row = Row(
        family=family,
        A=A,
        P=P,
        w=w,
        core_m=exact_m(P),
        full_m=exact_m(A),
        components=component_tuple,
        labels=labels,
        increments=increments,
        winding=winding,
        tournament=component_tournament(component_tuple),
    )
    if not row.components:
        raise AssertionError("a super-1/12 core must have nonempty E_1/13")
    # THM-765(A), checked independently against exact M on every row.
    if row.all_covered != (row.full_m <= L):
        raise AssertionError("component containment disagrees with exact full maximum")
    if row.all_covered and any(component.tooth_label is None for component in row.components):
        raise AssertionError("an all-covered component must have a unique tooth label")
    return row


def canonical_row_lines(row: Row) -> list[str]:
    t = row.tournament
    lines = [
        "|".join(
            [
                row.family,
                ",".join(map(str, row.A)),
                str(row.w),
                ftext(row.core_m),
                ftext(row.full_m),
                str(row.labels),
                str(row.increments),
                str(row.winding),
                str(t.score_histogram),
                str(t.directed_triangles),
                str(t.scc_sizes),
                str(t.edge_flips),
                str(t.phase_ties),
                str(t.hamiltonian_paths),
            ]
        )
    ]
    for component in row.components:
        lines.append(
            "|".join(
                [
                    ftext(component.lo),
                    ftext(component.hi),
                    ftext(component.center),
                    ftext(component.half_width),
                    ",".join(map(str, component.left_owners)),
                    ",".join(map(str, component.right_owners)),
                    str(component.tooth_label),
                    ftext(component.tooth_phase),
                    ftext(component.sigma),
                ]
            )
        )
    return lines


def compact_row(row: Row) -> str:
    t = row.tournament
    return (
        f"A={row.A} w={row.w} M(P)={ftext(row.core_m)} M(A)={ftext(row.full_m)} "
        f"k={len(row.components)} covered={row.all_covered} "
        f"zero/negative={row.zero_slacks}/{row.negative_slacks} "
        f"min_sigma={ftext(row.minimum_sigma)} "
        f"wind={row.winding} owner_pure={row.endpoint_pure} "
        f"T(score={t.score_histogram},C3={t.directed_triangles},SCC={t.scc_sizes},"
        f"flips={t.edge_flips},ties={t.phase_ties},H={t.hamiltonian_paths})"
    )


def print_detailed_row(title: str, row: Row) -> None:
    print(f"\n{title}")
    print("  " + compact_row(row))
    print(f"  labels={row.labels}")
    print(f"  increments={row.increments}; winding_checksum={row.winding}")
    for i, component in enumerate(row.components):
        print(
            f"  J{i:02d}=({ftext(component.lo)},{ftext(component.hi)}) "
            f"c={ftext(component.center)} h={ftext(component.half_width)} "
            f"owners={component.left_owners}->{component.right_owners} "
            f"label={component.tooth_label} phase={ftext(component.tooth_phase)} "
            f"sigma={ftext(component.sigma)}"
        )


def basic_unlabelled_signature(row: Row) -> tuple[object, ...]:
    """A deliberately coarse winding/owner/tournament summary for collision tests."""

    owner_sizes = Counter()
    for component in row.components:
        owner_sizes[len(component.left_owners)] += 1
        owner_sizes[len(component.right_owners)] += 1
    increment_hist = None if row.increments is None else tuple(sorted(Counter(row.increments).items()))
    t = row.tournament
    return (
        len(row.components),
        row.winding,
        increment_hist,
        tuple(sorted(owner_sizes.items())),
        t.score_histogram,
        t.scc_sizes,
        t.edge_flips,
        t.phase_ties,
    )


def main() -> None:
    masks_12 = maximal_strict_safe_masks(N, CORE_FLOOR)
    masks_10 = maximal_strict_safe_masks(N, NEAR_CEILING)

    primitive_sets = 0
    super_max_cores = 0
    branch_specs: list[tuple[tuple[int, ...], int]] = []
    for A in combinations(range(1, N + 1), 12):
        if gcd(*A) != 1:
            continue
        primitive_sets += 1
        w = A[-1]
        P = A[:-1]
        p_mask = set_mask(P)
        a_mask = p_mask | (1 << (w - 1))
        if not contained_in_some(p_mask, masks_12):
            continue
        super_max_cores += 1
        if contained_in_some(a_mask, masks_10):
            continue
        branch_specs.append((A, w))

    branch = tuple(audit_row("near-max-peel", A, w) for A, w in branch_specs)

    ap = tuple(range(1, 13))
    controls = tuple(
        audit_row("AP-control", ap, w)
        for w in ap
        if exact_m(tuple(v for v in ap if v != w)) > CORE_FLOOR
    )
    rows = controls + branch

    # Independent selector assertions after the exact maxima have been computed.
    if not all(row.core_m > CORE_FLOOR and row.full_m <= NEAR_CEILING for row in branch):
        raise AssertionError("mask-selected near branch failed exact-M replay")
    if not all(row.full_m == L and row.all_covered for row in controls):
        raise AssertionError("AP controls must be exactly tight and all-covered")

    digest = sha256()
    for row in rows:
        for line in canonical_row_lines(row):
            digest.update(line.encode("ascii"))
            digest.update(b"\n")

    print("LRC(13), n=12: exact tooth-label winding atlas (codex-S3)")
    print("STATUS: bounded exact evidence only; this is NOT a uniform sporadic-branch proof.")
    print(f"threshold L={ftext(L)}; N={N}; near ceiling={ftext(NEAR_CEILING)}")
    print(
        "bank: exhaustive primitive A subset [1,20], |A|=12, w=max(A), "
        "M(A\\{w})>1/12, M(A)<=1/10"
    )
    print(
        f"selector: primitive_sets={primitive_sets}; super_max_cores={super_max_cores}; "
        f"near_rows={len(branch)}; maximal_masks@1/12={len(masks_12)}; "
        f"maximal_masks@1/10={len(masks_10)}"
    )
    print(f"controls: eligible AP deletions={len(controls)} (w=1,...,11; w=12 has M(P)=1/12)")
    print(f"full_exact_payload_sha256={digest.hexdigest()}")

    print("\nExact global checks")
    print(
        f"  component/M equivalence mismatches=0/{len(rows)}; "
        f"all-covered labels nonunique=0/{sum(row.all_covered for row in rows)}"
    )
    print(
        f"  all-covered: controls={sum(row.all_covered for row in controls)}/{len(controls)}, "
        f"near branch={sum(row.all_covered for row in branch)}/{len(branch)}"
    )
    print(
        "  M(A) histogram (near branch): "
        + str(tuple((ftext(value), count) for value, count in sorted(Counter(row.full_m for row in branch).items())))
    )
    print(
        "  component-count histogram: "
        + str(tuple(sorted(Counter(len(row.components) for row in branch).items())))
    )
    print(
        "  completion-w histogram: "
        + str(tuple(sorted(Counter(row.w for row in branch).items())))
    )
    print(
        "  negative-component histogram: "
        + str(tuple(sorted(Counter(row.negative_slacks for row in branch).items())))
    )

    winding_hist = Counter(row.winding for row in branch)
    endpoint_pure_count = sum(row.endpoint_pure for row in branch)
    print("\nTooth-label winding scout")
    print(f"  winding histogram (near branch)={tuple(sorted(winding_hist.items(), key=lambda x: str(x[0])))}")
    print(f"  endpoint-pure rows={endpoint_pure_count}/{len(branch)}")
    print(
        "  raw proposed separator W1: endpoint-pure + one cyclic tooth winding => all-covered."
    )
    w1_liars = [row for row in branch if row.endpoint_pure and row.winding == 1]
    if w1_liars:
        print("  W1 REFUTED; first lexicographic liar: " + compact_row(w1_liars[0]))
    else:
        print("  W1 has no liar in this bank (evidence only).")

    print(
        "  tournament refinement T1: W1 + transitive phase tournament (C3=0,H=1) => all-covered."
    )
    t1_liars = [
        row
        for row in w1_liars
        if row.tournament.directed_triangles == 0
        and row.tournament.hamiltonian_paths == 1
    ]
    if t1_liars:
        print("  T1 REFUTED; first lexicographic liar: " + compact_row(t1_liars[0]))
    else:
        print("  T1 has no liar in this bank (evidence only).")

    control_signatures: dict[tuple[object, ...], Row] = {}
    for row in controls:
        control_signatures.setdefault(basic_unlabelled_signature(row), row)
    signature_collision: tuple[Row, Row] | None = None
    for row in branch:
        control = control_signatures.get(basic_unlabelled_signature(row))
        if control is not None:
            signature_collision = (control, row)
            break
    if signature_collision is None:
        print(
            "  coarse combined winding/owner/tournament signature separates the AP controls "
            "from this bank (bounded evidence only)."
        )
    else:
        control, liar = signature_collision
        print(
            f"  coarse combined signature COLLISION: AP deletion w={control.w} and "
            f"escaping A={liar.A}."
        )

    tournament_profiles = Counter(
        (row.tournament.score_histogram, row.tournament.scc_sizes) for row in branch
    )
    hp_rows = [row for row in branch if row.tournament.hamiltonian_paths is not None]
    print("\nTournament Analysis over safe components")
    print(
        "  observable=within-tooth midpoint phase displacement; gauge=clockwise half-circle; "
        "tie path=chronological component order."
    )
    print(f"  distinct (score histogram,SCC) profiles={len(tournament_profiles)}")
    print(
        f"  directed-3-cycle range={min(row.tournament.directed_triangles for row in branch)}.."
        f"{max(row.tournament.directed_triangles for row in branch)}; "
        f"edge-flip range={min(row.tournament.edge_flips for row in branch)}.."
        f"{max(row.tournament.edge_flips for row in branch)}"
    )
    print(
        "  SCC-size histogram="
        + str(tuple(sorted(Counter(row.tournament.scc_sizes for row in branch).items())))
    )
    print(
        f"  exact H(T) evaluated for k<={HP_COMPONENT_CAP}: {len(hp_rows)}/{len(branch)} rows; "
        f"range={min(row.tournament.hamiltonian_paths for row in hp_rows)}.."
        f"{max(row.tournament.hamiltonian_paths for row in hp_rows)}; "
        f"distinct={len({row.tournament.hamiltonian_paths for row in hp_rows})}"
    )

    print("\nAP tight/all-covered controls")
    for row in controls:
        print("  " + compact_row(row))

    closest_by_m = min(branch, key=lambda row: (row.full_m, row.A))
    closest_by_tooth = min(branch, key=lambda row: (-row.minimum_sigma, row.A))
    deepest_tooth = min(branch, key=lambda row: (row.minimum_sigma, row.A))
    max_components = max(branch, key=lambda row: (len(row.components), tuple(-v for v in row.A)))
    print("\nExtremal branch summaries")
    print("  smallest M(A): " + compact_row(closest_by_m))
    print("  smallest worst-tooth deficit: " + compact_row(closest_by_tooth))
    print("  largest worst-tooth deficit: " + compact_row(deepest_tooth))
    print("  most components (lex first tie): " + compact_row(max_components))

    print_detailed_row("Detailed tight control (richest AP deletion, w=11)", controls[-1])
    print_detailed_row("Detailed first W1 liar", w1_liars[0])
    if closest_by_m.A != w1_liars[0].A:
        print_detailed_row("Detailed closest escaping row by M(A)", closest_by_m)

    print("\nInterpretation")
    print("  * Exact sigma, not the raw cyclic label winding, detects tooth escape.")
    print("  * The first liar already has one winding and pure endpoint ownership; metric width/phase is indispensable.")
    print("  * Component tournaments retain within-tooth phase order but still quotient away sigma's width term.")
    print("  * Challenge to the vertex assumption: components are more faithful than runners, yet a proof object must")
    print("    carry (component interval, endpoint owners, tooth phase, width, deletion deck), not only a tournament.")


if __name__ == "__main__":
    main()
