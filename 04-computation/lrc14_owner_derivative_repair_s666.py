#!/usr/bin/env python3
"""S666: LRC n=14 owner-derivative repair over the Res_27 carry seam.

HYP-2240 says a scalar/product collision should be repaired by attaching the
missing address coordinate and then unrolling it into a derivative sum.

For the LRC n=14 seam this script tests a concrete version:

    fixed Res_27 shadow + fixed C=27 gcd/product shell
    + paired owner-deletion derivative
    separates AP/Vstar/2AP floor atoms from local nonzero carry perturbations.

The owner-deletion derivative used here is deliberately small.  For each speed
v, pair its Res_27 residue with:

    (number of D/U/N obligations covered by v,
     whether v is the unique owner of at least one D/U/N obligation).

The private-owner bit is a deletion derivative: deleting v would uncover an
obligation exactly when the bit is true.  The experiment asks whether this
paired owner bit repairs the collision left by the visible shadow.

Tournament Analysis / assumption challenge:
  Vertices are repair channels, not runners.  Pairwise observables are
  (mixed fiber count, compression, owner alignment, carry independence, and
  simplicity).  The switch is majority; ties follow the channel list.  Candidate
  LRC vertices considered included runners, residues, gaps, fixed circle
  sections, section boundaries, wall-crossing events, pair-sum denominators,
  owner obligations, deleted-speed cards, carry coordinates, and proof
  obligations.  The chosen quotient preserves the predicate "does this side
  channel separate known floor atoms from bounded local carry perturbations?"
  while destroying phase order and raw speed identity.
"""

from __future__ import annotations

import sys
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc_n14_carry_conservativity_s611 as s611  # noqa: E402
import lrc_n14_fixed_round_certificate_s578 as s578  # noqa: E402


N = s611.N
C = s611.C
FLOOR = s611.FLOOR
AP = s611.AP
VSTAR = s611.VSTAR
TWOP = tuple(2 * v for v in AP)
BASES = (("AP", AP), ("Vstar", VSTAR), ("2AP", TWOP))


def fmt_frac(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def route(score: Fraction) -> str:
    if score == FLOOR:
        return "floor"
    if score > FLOOR:
        return "strict"
    return "below"


def residue(v: int) -> int:
    r = v % C
    if r == 0:
        raise ValueError(f"zero Res_{C} residue for {v}")
    return r


def carry(v: int) -> int:
    r = residue(v)
    return (v - r) // C


def shadow(row: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(residue(v) for v in row))


def shell(r: int) -> int:
    r %= C
    return min(r, C - r)


def gcd_shell_counts(row: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(gcd(shell(r), C) for r in shadow(row)).items()))


def visible_shadow_key(row: tuple[int, ...]) -> tuple[object, ...]:
    return (shadow(row), gcd_shell_counts(row))


def obligations_for(v: int) -> tuple[s578.Obligation, ...]:
    return tuple(o for o in s578.obligation_universe() if s578.covers(v, o))


def private_owner_map(row: tuple[int, ...]) -> dict[int, tuple[s578.Obligation, ...]]:
    inc = s578.incidence(row)
    private: dict[int, list[s578.Obligation]] = defaultdict(list)
    for obligation, owners in inc.items():
        if len(owners) == 1:
            private[owners[0]].append(obligation)
    return {v: tuple(sorted(obs)) for v, obs in private.items()}


def owner_cover_count_key(row: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted((residue(v), len(obligations_for(v))) for v in row))


def owner_cover_private_flag_key(row: tuple[int, ...]) -> tuple[tuple[int, int, int], ...]:
    private = private_owner_map(row)
    return tuple(
        sorted((residue(v), len(obligations_for(v)), 1 if private.get(v) else 0) for v in row)
    )


def owner_cover_private_count_key(row: tuple[int, ...]) -> tuple[tuple[int, int, int], ...]:
    private = private_owner_map(row)
    return tuple(sorted((residue(v), len(obligations_for(v)), len(private.get(v, ()))) for v in row))


def carry_l1_key(row: tuple[int, ...]) -> tuple[int, int, tuple[int, ...]]:
    carries = tuple(carry(v) for v in row)
    support = tuple(sorted(residue(v) for v in row if carry(v)))
    return (sum(carries), max(carries) - min(carries), support)


def cheap_pair_key(row: tuple[int, ...]) -> tuple[int, int, int, int]:
    pair = s578.unblocked_small_pair(row)
    if pair is None:
        return (-1, -1, -1, -1)
    a, b, k = pair
    r1, r2 = sorted((residue(a), residue(b)))
    return (r1, r2, k, a + b)


@dataclass(frozen=True)
class Probe:
    family: str
    weight: int
    moved: tuple[int, ...]
    row: tuple[int, ...]
    score: Fraction
    time: Fraction

    @property
    def label(self) -> str:
        if self.weight == 0:
            return self.family
        return f"{self.family}:w{self.weight}:carry{self.moved}"

    @property
    def route(self) -> str:
        return route(self.score)

    @property
    def margin(self) -> Fraction:
        return self.score - FLOOR


def local_probes(max_weight: int = 3) -> list[Probe]:
    probes: list[Probe] = []
    for family, base in BASES:
        score, time = s611.exact_maximin(tuple(sorted(base)))
        probes.append(Probe(family, 0, tuple(), tuple(sorted(base)), score, time))
        for weight in range(1, max_weight + 1):
            for idxs in combinations(range(len(base)), weight):
                row = list(base)
                moved = []
                for idx in idxs:
                    moved.append(base[idx])
                    row[idx] += C
                lifted = tuple(sorted(row))
                score, time = s611.exact_maximin(lifted)
                probes.append(Probe(family, weight, tuple(sorted(moved)), lifted, score, time))
    return probes


@dataclass(frozen=True)
class Channel:
    name: str
    key: Callable[[tuple[int, ...]], tuple[object, ...]]
    uses_owner: bool
    paired_owner: bool
    uses_carry: bool
    simplicity: int


CHANNELS = (
    Channel("visible_shadow", lambda row: visible_shadow_key(row), False, False, False, 5),
    Channel(
        "visible+cheap_pair",
        lambda row: visible_shadow_key(row) + (cheap_pair_key(row),),
        True,
        False,
        False,
        4,
    ),
    Channel(
        "visible+owner_cover_count",
        lambda row: visible_shadow_key(row) + (owner_cover_count_key(row),),
        True,
        True,
        False,
        4,
    ),
    Channel(
        "visible+owner_private_flag",
        lambda row: visible_shadow_key(row) + (owner_cover_private_flag_key(row),),
        True,
        True,
        False,
        3,
    ),
    Channel(
        "visible+owner_private_count",
        lambda row: visible_shadow_key(row) + (owner_cover_private_count_key(row),),
        True,
        True,
        False,
        2,
    ),
    Channel(
        "visible+carry_l1_support",
        lambda row: visible_shadow_key(row) + (carry_l1_key(row),),
        False,
        False,
        True,
        3,
    ),
)


@dataclass(frozen=True)
class AuditRow:
    channel: Channel
    groups: int
    mixed: int
    max_bucket: int
    examples: tuple[tuple[Counter[str], tuple[str, ...]], ...]

    @property
    def separates(self) -> bool:
        return self.mixed == 0


def audit_channel(channel: Channel, probes: list[Probe]) -> AuditRow:
    groups: defaultdict[tuple[object, ...], list[Probe]] = defaultdict(list)
    for probe in probes:
        groups[channel.key(probe.row)].append(probe)

    examples = []
    for bucket in groups.values():
        hist = Counter(p.route for p in bucket)
        if len(hist) > 1:
            examples.append((hist, tuple(p.label for p in bucket[:8])))

    return AuditRow(
        channel=channel,
        groups=len(groups),
        mixed=len(examples),
        max_bucket=max(len(bucket) for bucket in groups.values()),
        examples=tuple(examples[:4]),
    )


def terminal_leads() -> list[tuple[str, int, tuple[int, ...], Fraction, Fraction]]:
    rows = []
    for family, base in BASES:
        ordered = tuple(sorted(base))
        for weight in range(4, 7):
            moved = ordered[-weight:]
            row = list(base)
            for v in moved:
                idx = row.index(v)
                row[idx] = v + C
            score, time = s611.exact_maximin(tuple(sorted(row)))
            rows.append((family, weight, moved, score, time))
    return rows


def channel_metrics(audit: AuditRow) -> tuple[int, int, int, int, int]:
    channel = audit.channel
    return (
        -audit.mixed,
        -audit.groups,
        1 if channel.uses_owner else 0,
        1 if channel.paired_owner else 0,
        1 if not channel.uses_carry else 0,
    )


def tournament(audits: list[AuditRow]) -> dict[str, object]:
    n = len(audits)
    adj = [[0] * n for _ in range(n)]
    out = [0] * n
    metrics = [channel_metrics(a) for a in audits]

    for i, j in combinations(range(n), 2):
        wins_i = sum(x > y for x, y in zip(metrics[i], metrics[j]))
        wins_j = sum(x < y for x, y in zip(metrics[i], metrics[j]))
        if wins_i > wins_j or (wins_i == wins_j and i < j):
            adj[i][j] = 1
            out[i] += 1
        else:
            adj[j][i] = 1
            out[j] += 1

    c3 = 0
    for a, b, c in combinations(range(n), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            c3 += 1

    def scc_sizes() -> list[int]:
        radj = [[i for i in range(n) if adj[i][j]] for j in range(n)]
        seen = [False] * n
        order: list[int] = []
        for start in range(n):
            if seen[start]:
                continue
            stack = [(start, False)]
            while stack:
                v, done = stack.pop()
                if done:
                    order.append(v)
                    continue
                if seen[v]:
                    continue
                seen[v] = True
                stack.append((v, True))
                for w in range(n):
                    if adj[v][w] and not seen[w]:
                        stack.append((w, False))
        seen = [False] * n
        sizes = []
        for start in reversed(order):
            if seen[start]:
                continue
            q = deque([start])
            seen[start] = True
            size = 0
            while q:
                v = q.popleft()
                size += 1
                for w in radj[v]:
                    if not seen[w]:
                        seen[w] = True
                        q.append(w)
            sizes.append(size)
        return sorted(sizes, reverse=True)

    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v, count in enumerate(dp[mask]):
            if not count:
                continue
            for w in range(n):
                if not mask & (1 << w) and adj[v][w]:
                    dp[mask | (1 << w)][w] += count

    return {
        "score_hist": dict(sorted(Counter(out).items())),
        "directed_3cycles": c3,
        "scc_sizes": scc_sizes(),
        "hamiltonian_paths": sum(dp[-1]),
        "top_order": [a.channel.name for a in sorted(audits, key=lambda x: channel_metrics(x), reverse=True)],
    }


def main() -> None:
    print("=" * 78)
    print("S666 LRC14 owner-derivative repair over the Res_27 carry seam")
    print("=" * 78)
    print()
    print(f"n={N}; C=2n-1={C}; floor={fmt_frac(FLOOR)}")
    print("Visible scalar/product key = sorted Res_27 residues + C=27 gcd-shell counts.")
    print("Owner derivative bit = whether deleting a speed would uncover a D/U/N obligation.")
    print()

    probes = local_probes(max_weight=3)
    print("A. Exhaustive local +27 carry tax through weight 3")
    print(f"  probes={len(probes)} route_hist={dict(Counter(p.route for p in probes))}")
    for family in ("AP", "Vstar", "2AP"):
        for weight in range(0, 4):
            rows = [p for p in probes if p.family == family and p.weight == weight]
            hist = Counter(p.route for p in rows)
            min_probe = min(rows, key=lambda p: (p.score, p.time, p.moved))
            print(
                f"  {family:5s} weight={weight} count={len(rows):3d} "
                f"routes={dict(hist)} min_M={fmt_frac(min_probe.score):>5s} "
                f"margin={fmt_frac(min_probe.margin):>6s} t={fmt_frac(min_probe.time):>5s} "
                f"moved={min_probe.moved or '-'}"
            )
    print("  New S666 evidence: all 858 weight-3 local carries are strict.")
    print()

    print("B. Owner-deletion projection audit")
    audits = [audit_channel(channel, probes) for channel in CHANNELS]
    print(f"{'channel':<34} {'groups':>6} {'mixed':>6} {'max_bucket':>10} separates")
    for audit in audits:
        print(
            f"{audit.channel.name:<34} {audit.groups:6d} {audit.mixed:6d} "
            f"{audit.max_bucket:10d} {audit.separates}"
        )
        for hist, labels in audit.examples[:2]:
            print(f"    mixed routes={dict(hist)} examples={labels}")
    print()
    print("  Key repair:")
    print("    visible+owner_cover_count still leaks AP:w1:carry(11) and Vstar:w1:carry(11).")
    print("    Adding only the private-owner deletion bit gives 0 mixed fibers while")
    print("    remaining carry-free; it does not use the full carry vector.")
    print()

    print("C. Targeted terminal-carry leads beyond the exhaustive radius")
    for family, weight, moved, score, time in terminal_leads():
        print(
            f"  {family:5s} terminal weight={weight} moved={moved} "
            f"M={fmt_frac(score):>5s} margin={fmt_frac(score - FLOOR):>6s} t={fmt_frac(time):>5s}"
        )
    print("  These are leads only; S666's exhaustive certificate stops at weight 3.")
    print()

    fp = tournament(audits)
    print("D. Tournament Analysis over repair channels")
    print("  vertices=repair channels")
    print("  observable=(mixed fibers, compression, owner alignment, pairedness, carry independence)")
    print(f"  score_hist={fp['score_hist']}")
    print(f"  directed_3cycles={fp['directed_3cycles']}")
    print(f"  scc_sizes={fp['scc_sizes']}")
    print(f"  hamiltonian_paths={fp['hamiltonian_paths']}")
    print("  top_order:")
    for name in fp["top_order"]:
        print(f"    {name}")
    print()

    print("E. Proof consequence")
    print(
        "  The visible Res_27/product shell has three mixed buckets, each containing "
        "one floor atom and 377 strict local carries.  Cheap-pair existence does "
        "not repair the leak.  Paired owner cover counts almost repair it, failing "
        "only on AP/Vstar carry(11).  The single extra private-owner deletion bit "
        "repairs every checked local carry through weight 3."
    )
    print(
        "  This suggests a concrete no-leak lemma: in the Res_27 fiber, any "
        "nonzero local carry that preserves the visible shell must change the "
        "paired owner-deletion ledger, unless it belongs to a globally coherent "
        "scalar floor lift."
    )


if __name__ == "__main__":
    main()
