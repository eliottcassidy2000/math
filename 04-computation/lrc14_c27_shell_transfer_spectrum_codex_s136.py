#!/usr/bin/env python3
"""S136: C=27 shell-transfer spectrum for the LRC14 low-gap frontier.

This follows S133--S134 and the new binary-relational mandate.  The vertex set
is not the runners.  It is the thirteen antipodal summand shells

    P_a = {a, 27-a}, 1 <= a <= 13.

The quotient preserves exactly the data used by the S571 second-gap bridge:
which shells are unit-visible, which nonunit gcd stratum they occupy, and where
a row has a hole/double transfer relative to the AP lower transversal.

The main computational question is deliberately modest:

    in the S130 AP/GW/single-replacement atlas, what are the marked C=27
    transfers on the low-gap frontier M <= 2/27?

The script also builds a shell-priority tournament for each selected row.  Its
edges are proof-priority edges between shell pairs, not clock geometry: holes
and surplus shells outrank ordinary occupied shells; unit-visible shells outrank
nonunit shells; ties use the lower shell label as a Hamiltonian path.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
from pathlib import Path
import argparse
import sys


REPO = Path(__file__).resolve().parents[1]
N = 14
C = 2 * N - 1


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s134 = load_module(
    "s134_bigraded_for_s135",
    REPO / "04-computation" / "lrc14_bigraded_relation_signature_codex_s134.py",
)
s127 = s134.s127


@dataclass(frozen=True)
class Shell:
    a: int
    gcd27: int
    lower_count: int
    upper_count: int

    @property
    def total(self) -> int:
        return self.lower_count + self.upper_count

    @property
    def side_word(self) -> str:
        if self.total == 0:
            return "hole"
        if self.total > 2:
            return "multi"
        if self.total == 2:
            return "double"
        if self.lower_count:
            return "lower"
        return "upper"

    @property
    def detail(self) -> str:
        return f"{self.a}:{self.gcd27}:{self.lower_count}{self.upper_count}"


@dataclass(frozen=True)
class RowItem:
    label: str
    speeds: tuple[int, ...]
    M: F
    branch: str
    q_threshold: int
    shells: tuple[Shell, ...]
    zero_count: int
    shell_mask: int
    shell_edge_flips_from_ap: int

    @property
    def holes(self) -> tuple[Shell, ...]:
        return tuple(s for s in self.shells if s.total == 0)

    @property
    def doubles(self) -> tuple[Shell, ...]:
        return tuple(s for s in self.shells if s.total >= 2)


def shell_profile(speeds: tuple[int, ...]) -> tuple[tuple[Shell, ...], int]:
    residues = Counter(v % C for v in speeds)
    shells = []
    for a in range(1, (C + 1) // 2):
        shells.append(
            Shell(
                a=a,
                gcd27=gcd(a, C),
                lower_count=residues[a],
                upper_count=residues[C - a],
            )
        )
    return tuple(shells), residues[0]


def transfer_label(item: RowItem, detailed: bool = True) -> str:
    if not item.holes and not item.doubles and item.zero_count == 0:
        return "perfect-transversal"
    h = ",".join(f"{s.a}:g{s.gcd27}" for s in item.holes) or "-"
    if detailed:
        d = ",".join(s.detail for s in item.doubles) or "-"
    else:
        d = ",".join(f"{s.a}:g{s.gcd27}" for s in item.doubles) or "-"
    z = str(item.zero_count) if item.zero_count else "-"
    return f"H[{h}] D[{d}] Z[{z}]"


def coarse_transfer_label(item: RowItem) -> str:
    if not item.holes and not item.doubles and item.zero_count == 0:
        return "perfect"
    h = tuple(s.gcd27 for s in item.holes)
    d = tuple(s.gcd27 for s in item.doubles)
    z = item.zero_count
    return f"Hg{h or '-'}->Dg{d or '-'};z={z}"


def shell_priority(shell: Shell) -> tuple[int, int, int, int]:
    """Role score used for the shell-priority tournament.

    This is a proof-carrier order, not a geometric order.  Holes are most urgent
    because a missed unit shell creates a second-gap witness.  Doubles are the
    compensating transfer packet.  Unit shells outrank nonunit shells because
    multiplication by units can move them to the observer.
    """
    status_rank = {
        "hole": 5,
        "multi": 4,
        "double": 4,
        "upper": 2,
        "lower": 1,
    }[shell.side_word]
    visibility_rank = {1: 3, 3: 2, 9: 1}[shell.gcd27]
    side_rank = shell.upper_count - shell.lower_count
    return (status_rank, visibility_rank, side_rank, -shell.a)


def tournament_from_shells(shells: tuple[Shell, ...]) -> int:
    mask = 0
    bit = 0
    scores = [shell_priority(s) for s in shells]
    for i, j in combinations(range(len(shells)), 2):
        if scores[i] >= scores[j]:
            mask |= 1 << bit
        bit += 1
    return mask


def edge_flips(mask_a: int, mask_b: int, n: int) -> int:
    return sum(
        1
        for i, j in combinations(range(n), 2)
        if s127.edge(mask_a, n, i, j) != s127.edge(mask_b, n, i, j)
    )


def build_items(max_replacement: int) -> list[RowItem]:
    rows = s134.s130.candidate_rows(max_replacement)
    ap_shells, _ = shell_profile(tuple(range(1, 14)))
    ap_mask = tournament_from_shells(ap_shells)
    items: list[RowItem] = []
    for row in rows:
        sig = s134.signature(row)
        shells, zero_count = shell_profile(row.speeds)
        mask = tournament_from_shells(shells)
        items.append(
            RowItem(
                label=row.label,
                speeds=row.speeds,
                M=sig.M,
                branch=sig.farey_branch,
                q_threshold=sig.q_threshold,
                shells=shells,
                zero_count=zero_count,
                shell_mask=mask,
                shell_edge_flips_from_ap=edge_flips(ap_mask, mask, len(shells)),
            )
        )
    return items


def selected_labels() -> tuple[str, ...]:
    return (
        "AP",
        "GW 12->24",
        "near-miss 12->36",
        "swap 10->20",
        "swap 13->26",
        "residue-liar 12->26",
        "swap 12->60",
    )


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge / quotient declaration")
    print("  considered vertices:")
    print("    runners, residues, C=27 shell pairs, shell transfers, Farey branches,")
    print("    Kpq incidence owners, visible folds, hidden balanced relations, and")
    print("    proof obligations.")
    print("  chosen vertices:")
    print("    antipodal summand shells P_a={a,27-a}.")
    print("  quotient preserves:")
    print("    unit visibility, gcd stratum, AP-lower transversal holes/doubles, and")
    print("    the marked transfer that S571 says multiplication can or cannot see.")
    print("  quotient destroys:")
    print("    exact speed size, off-shell denominators, and full time geometry; those")
    print("    are reattached through exact M(S), q-threshold, and Farey branch.")
    print("  challenged assumption:")
    print("    a runner-vertex tournament is the right first quotient.  For the C=27")
    print("    branch, shell-pair vertices preserve the missing-pair predicate that")
    print("    runner tournaments blur.")


def print_selected(items: list[RowItem]) -> None:
    by_label = {item.label: item for item in items}
    print()
    print("[1] Named C=27 shell-transfer signatures")
    print(
        f"  {'row':20s} {'M':>7s} {'branch':23s} {'qthr':>4s} "
        f"{'flips':>5s} {'transfer':45s} shell fp"
    )
    for label in selected_labels():
        item = by_label.get(label)
        if item is None:
            continue
        fp = s127.tournament_fingerprint(item.shell_mask, len(item.shells))
        shell_fp = f"c3={fp['c3']},hp={fp['hp']}"
        print(
            f"  {item.label:20s} {str(item.M):>7s} {item.branch:23s} "
            f"{item.q_threshold:4d} {item.shell_edge_flips_from_ap:5d} "
            f"{transfer_label(item):45s} {shell_fp}"
        )
    print()
    print("  Reading:")
    print("    AP is the perfect lower transversal.  GW has the nonunit transfer")
    print("    H[12:g3] -> D[3:g3].  The near-miss keeps the same hole but moves")
    print("    surplus to the unique gcd-9 shell, H[12:g3] -> D[9:g9], matching")
    print("    the first K33/Farey child.  The two 2/27 petals are the unit-hole")
    print("    second-gap cases.")


def print_frontier(items: list[RowItem], cutoff: F, name: str) -> None:
    front = sorted((item for item in items if item.M <= cutoff), key=lambda x: (x.M, x.label))
    print()
    print(f"[2] Low-gap frontier {name}: M <= {cutoff}")
    print(f"  rows={len(front)}")
    print(f"  branch_counts={dict(sorted(Counter(item.branch for item in front).items()))}")
    for item in front:
        print(
            f"    {item.label:18s} M={str(item.M):>5s} branch={item.branch:23s} "
            f"qthr={item.q_threshold:2d} flips={item.shell_edge_flips_from_ap:2d} "
            f"{transfer_label(item)}"
        )


def print_bank_summary(items: list[RowItem]) -> None:
    print()
    print("[3] Whole-bank shell-transfer summary")
    print(f"  audited rows={len(items)}")
    print(f"  branch_counts={dict(sorted(Counter(item.branch for item in items).items()))}")
    low = [item for item in items if item.M <= F(2, 27)]
    print(f"  rows with M<=2/27={len(low)}")
    print("  low-frontier coarse transfer counts:")
    for label, count in sorted(Counter(coarse_transfer_label(item) for item in low).items()):
        examples = [item.label for item in low if coarse_transfer_label(item) == label]
        print(f"    {count:2d} {label:24s} examples={examples}")
    print()
    print("  transfer labels that mix tight/K33/C27 with loose rows:")
    by_label: dict[str, list[RowItem]] = defaultdict(list)
    for item in items:
        by_label[transfer_label(item)].append(item)
    mixed = []
    for label, group in by_label.items():
        branches = {item.branch for item in group}
        if len(branches) > 1 and any(item.M <= F(2, 27) for item in group):
            mixed.append((label, group))
    for label, group in sorted(mixed, key=lambda x: (min(i.M for i in x[1]), x[0]))[:8]:
        branches = Counter(item.branch for item in group)
        examples = ", ".join(f"{item.label}:{item.M}" for item in sorted(group, key=lambda x: (x.M, x.label))[:5])
        print(f"    {label:45s} branches={dict(branches)} examples={examples}")
    print()
    print("  Guardrail:")
    print("    marked shell transfer is necessary relation data but not sufficient")
    print("    by itself.  Exact M/Farey branch must remain attached, because some")
    print("    transfer labels recur in safely loose rows.")


def print_channel_tournament() -> None:
    print()
    print("[4] Tournament Analysis on proof quotients")
    channels = [
        ("exact M/Farey branch", (5, 5, 4, 5, 4)),
        ("marked C27 shell transfer", (4, 5, 5, 4, 5)),
        ("unit-visible shell holes", (3, 5, 5, 3, 4)),
        ("nonunit 3-adic depth", (3, 4, 4, 3, 5)),
        ("Kpq/K33 incidence owners", (3, 4, 3, 5, 4)),
        ("shell-priority tournament", (2, 4, 4, 3, 4)),
        ("raw visible fold count", (1, 2, 3, 2, 2)),
        ("raw runner residues", (0, 1, 1, 1, 1)),
    ]
    mask = 0
    bit = 0
    wins = {i: set() for i in range(len(channels))}
    for i, j in combinations(range(len(channels)), 2):
        if channels[i][1] >= channels[j][1]:
            winner, loser = i, j
            mask |= 1 << bit
        else:
            winner, loser = j, i
        wins[winner].add(loser)
        bit += 1
    fp = s127.tournament_fingerprint(mask, len(channels))
    order = sorted(range(len(channels)), key=lambda i: channels[i][1], reverse=True)
    print("  vertices: proof quotients/channels, not runners.")
    print("  pair observable:")
    print("    theorem scale, shell-predicate retention, visibility/sign retention,")
    print("    HYP-2908 handoff strength, scalar-decoy resistance.")
    print("  switch/gauge: lexicographically larger role score wins.")
    print(f"  fingerprint: score_hist={fp['score_hist']} c3={fp['c3']} scc={fp['scc']} hp={fp['hp']}")
    print("  Hamiltonian path:")
    print("    " + " > ".join(channels[i][0] for i in order))


def print_proof_readout() -> None:
    print()
    print("[5] Proof readout")
    print("  S136 turns the C=27 branch into a marked transfer problem.")
    print("  Bounded frontier evidence through the S130 single-replacement atlas says:")
    print("    M<=3/41 gives exactly AP, GW, and the 12->36 K33 near-miss;")
    print("    M<=2/27 adds exactly the two unit-hole petals 10->20 and 13->26.")
    print("  This suggests a sharper lemma target:")
    print("    after standard reductions, any low-gap non-AP/GW atom must either")
    print("    have a unit-visible C=27 hole, hence belongs to the second-gap/petal")
    print("    branch, or move a gcd-3 nonunit hole into the gcd-9 shell, hence")
    print("    carries the first K33/Farey-child state-lift packet.")
    print("  The warning is equally important: shell transfer alone is not the proof.")
    print("  The exact Farey branch remains attached to prevent safely loose rows from")
    print("  masquerading as frontier rows.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-replacement", type=int, default=140)
    args = parser.parse_args()

    print("S136 LRC14 C=27 SHELL-TRANSFER SPECTRUM")
    print("=" * 78)
    print(f"max_replacement={args.max_replacement}")
    items = build_items(args.max_replacement)
    print_assumption_challenge()
    print_selected(items)
    print_frontier(items, F(3, 41), "floor plus first K33 child")
    print_frontier(items, F(2, 27), "second-gap shell")
    print_bank_summary(items)
    print_channel_tournament()
    print_proof_readout()


if __name__ == "__main__":
    main()
