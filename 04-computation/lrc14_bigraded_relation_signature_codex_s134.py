#!/usr/bin/env python3
"""S134: bigraded summand/multiplicand relation signatures for LRC14.

This follows S130--S133 and reconnects them to the older S560/S571/S26/S28
summand/multiplicand graph work.

The central rule is:

  addition creates relation shells, multiplication tests their visibility.

For a speed row S we compute:

  * pair-sum support and excess over the AP minimum,
  * observer-visible folds a+b=c in S,
  * balanced pair-sum collisions a+b=c+d split by visible/hidden shell,
  * denominator clearance blockers for the exact optimum M(S)=p/q,
  * C=27 antipodal shell profile for the n=14 second-gap graph,
  * the Farey branch of M(S).

The point is not to prove LRC14.  It is to identify the next proof interface:
S133 separates the p=2 C=27 branch from the p>=3 K33 branch; S134 asks which
typed relation-channel should carry the remaining rigidity inside each branch.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from math import gcd
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[1]
N = 14
C_SECOND = 2 * N - 1
THR = F(1, N)


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s130 = load_module(
    "s130_mutated_farey_for_s134",
    REPO / "04-computation" / "lrc14_mutated_farey_tournament_codex_s130.py",
)
s127 = s130.s127


@dataclass(frozen=True)
class Signature:
    label: str
    family: str
    M: F
    q_threshold: int
    farey_e: int
    farey_branch: str
    sumset: int
    sumset_excess: int
    visible_folds: int
    visible_collision_pairs: int
    hidden_collision_pairs: int
    exact_denominator_blockers: int
    c27_lower_only: int
    c27_upper_only: int
    c27_both: int
    c27_empty: int
    c27_nonunit_events: int
    c27_zero_hits: int
    product_edges: int
    k_rank: str

    @property
    def relation_visibility(self) -> tuple[int, int, int]:
        return (self.visible_folds, self.visible_collision_pairs, -self.hidden_collision_pairs)


def exact_M(row: tuple[int, ...]) -> F:
    M, _pts = s130.s124.M_exact(row)
    return M


def farey_branch(M: F) -> str:
    p, q = M.numerator, M.denominator
    e = N * p - q
    if e == 0:
        return "tight-floor"
    if e == 1 and p == 1:
        return "q-parent-star"
    if e == 1 and p == 2 and q == C_SECOND:
        return "C27-petal-two-block"
    if e == 1 and p >= 3:
        return "K33-unit-excess"
    if e > 1:
        return "nonunit-excess"
    return "below-or-other"


def k_rank(p: int, q: int) -> str:
    m = min(p, q)
    if m >= 3:
        return "K33-wall"
    if m == 2:
        return "two-block"
    return "star"


def pair_sum_counter(speeds: tuple[int, ...]) -> Counter[int]:
    ctr: Counter[int] = Counter()
    for a, b in combinations(speeds, 2):
        ctr[a + b] += 1
    return ctr


def additive_signature(speeds: tuple[int, ...]) -> tuple[int, int, int, int, int]:
    ctr = pair_sum_counter(speeds)
    speed_set = set(speeds)
    sumset = len(ctr)
    ap_min = 2 * len(speeds) - 3
    visible_folds = sum(v for c, v in ctr.items() if c in speed_set)
    visible_collisions = 0
    hidden_collisions = 0
    for c, v in ctr.items():
        if v < 2:
            continue
        pairs = v * (v - 1) // 2
        if c in speed_set:
            visible_collisions += pairs
        else:
            hidden_collisions += pairs
    return sumset, sumset - ap_min, visible_folds, visible_collisions, hidden_collisions


def c27_shell_profile(speeds: tuple[int, ...]) -> tuple[int, int, int, int, int, int]:
    residues = {s % C_SECOND for s in speeds}
    lower_only = upper_only = both = empty = nonunit_events = 0
    zero_hits = 1 if 0 in residues else 0
    for a in range(1, (C_SECOND + 1) // 2):
        b = C_SECOND - a
        has_a = a in residues
        has_b = b in residues
        if has_a and has_b:
            both += 1
        elif has_a:
            lower_only += 1
        elif has_b:
            upper_only += 1
        else:
            empty += 1
        if gcd(a, C_SECOND) > 1 and (has_a or has_b):
            nonunit_events += int(has_a) + int(has_b)
    return lower_only, upper_only, both, empty, nonunit_events, zero_hits


def signature(row: s130.Row) -> Signature:
    speeds = row.speeds
    M = exact_M(speeds)
    p, q = M.numerator, M.denominator
    e = N * p - q
    sumset, excess, folds, vis_col, hid_col = additive_signature(speeds)
    c27 = c27_shell_profile(speeds)
    return Signature(
        label=row.label,
        family=row.family,
        M=M,
        q_threshold=s130.s124.q_threshold(speeds),
        farey_e=e,
        farey_branch=farey_branch(M),
        sumset=sumset,
        sumset_excess=excess,
        visible_folds=folds,
        visible_collision_pairs=vis_col,
        hidden_collision_pairs=hid_col,
        exact_denominator_blockers=sum(1 for s in speeds if s % q == 0),
        c27_lower_only=c27[0],
        c27_upper_only=c27[1],
        c27_both=c27[2],
        c27_empty=c27[3],
        c27_nonunit_events=c27[4],
        c27_zero_hits=c27[5],
        product_edges=p * q,
        k_rank=k_rank(p, q),
    )


def selected_rows() -> list[s130.Row]:
    ap = tuple(range(1, 14))
    rows = [
        s130.Row("AP", ap, "known tight"),
        s130.Row("GW 12->24", tuple(list(range(1, 12)) + [13, 24]), "known tight"),
        s130.Row("residue-liar 12->26", tuple(list(range(1, 12)) + [13, 26]), "q-threshold loose"),
        s130.Row("petal 10->20", tuple(sorted((set(ap) - {10}) | {20})), "C27 p=2 petal"),
        s130.Row("petal 13->26", tuple(sorted((set(ap) - {13}) | {26})), "C27 p=2 petal"),
        s130.Row("near-miss 12->36", tuple(list(range(1, 12)) + [13, 36]), "K33 near miss"),
        s130.Row("shiftAP +10", tuple(range(10, 23)), "translation decoy"),
        s130.Row("committed 84", tuple(list(range(1, 12)) + [13, 84]), "q-covering tail"),
    ]
    return rows


def print_selected_signatures() -> None:
    print("[1] Named-row bigraded signatures")
    print("    addition supplies shells; multiplication tests visibility.")
    header = (
        f"  {'row':22s} {'M':>7s} {'br':23s} {'sum+':>5s} {'fold':>5s} "
        f"{'visC':>5s} {'hidC':>5s} {'qblk':>4s} {'C27 L/U/B/E/N/Z':>19s} {'K':>9s}"
    )
    print(header)
    for sig in map(signature, selected_rows()):
        c27_text = (
            f"{sig.c27_lower_only}/{sig.c27_upper_only}/{sig.c27_both}/"
            f"{sig.c27_empty}/{sig.c27_nonunit_events}/{sig.c27_zero_hits}"
        )
        print(
            f"  {sig.label:22s} {str(sig.M):>7s} {sig.farey_branch:23s} "
            f"{sig.sumset_excess:5d} {sig.visible_folds:5d} "
            f"{sig.visible_collision_pairs:5d} {sig.hidden_collision_pairs:5d} "
            f"{sig.exact_denominator_blockers:4d} {c27_text:>19s} {sig.k_rank:>9s}"
        )
    print()
    print("  C27 columns are lower/upper/both/empty/nonunit-events/zero-hits.")
    print("  The shiftAP row is the calibration from HYP-2639: same AP-style")
    print("  sumset and energy, but no observer-visible low folds in the original")
    print("  LRC coordinate frame, so it is very safe despite identical raw additive")
    print("  structure.")


def row_bank(max_replacement: int = 70) -> list[Signature]:
    return [signature(row) for row in s130.candidate_rows(max_replacement)]


def average_int(values: list[int]) -> F:
    return F(sum(values), len(values)) if values else F(0)


def print_bank_summary(sigs: list[Signature]) -> None:
    print()
    print("[2] S130 row-bank branch summary")
    by_branch: dict[str, list[Signature]] = defaultdict(list)
    for sig in sigs:
        by_branch[sig.farey_branch].append(sig)
    print(
        f"  {'branch':23s} {'rows':>4s} {'avg sum+':>9s} {'avg fold':>9s} "
        f"{'avg hidC':>9s} {'avg qblk':>9s} {'K-ranks':>24s}"
    )
    for branch in sorted(by_branch):
        group = by_branch[branch]
        ranks = Counter(sig.k_rank for sig in group)
        print(
            f"  {branch:23s} {len(group):4d} "
            f"{str(average_int([s.sumset_excess for s in group])):>9s} "
            f"{str(average_int([s.visible_folds for s in group])):>9s} "
            f"{str(average_int([s.hidden_collision_pairs for s in group])):>9s} "
            f"{str(average_int([s.exact_denominator_blockers for s in group])):>9s} "
            f"{dict(sorted(ranks.items()))!s:>24s}"
        )
    print()
    print("  Reading:")
    print("    the exact M-branch already splits rows better than any raw additive")
    print("    statistic, but the branch summary says which relation channel should")
    print("    be kept after the split.  The C27 p=2 rows carry a two-block")
    print("    multiplicand tag; K33 unit-excess rows carry three-owner incidence.")


def print_extremes(sigs: list[Signature]) -> None:
    print()
    print("[3] Extremes that explain why scalarizing fails")
    print("  most visible folds:")
    for sig in sorted(sigs, key=lambda s: (-s.visible_folds, s.sumset_excess, s.label))[:6]:
        print(
            f"    {sig.label:18s} M={str(sig.M):>6s} branch={sig.farey_branch:23s} "
            f"fold={sig.visible_folds:3d} sum+={sig.sumset_excess:3d} qblk={sig.exact_denominator_blockers}"
        )
    print("  highest hidden balanced collision payload:")
    for sig in sorted(sigs, key=lambda s: (-s.hidden_collision_pairs, s.visible_folds, s.label))[:6]:
        print(
            f"    {sig.label:18s} M={str(sig.M):>6s} branch={sig.farey_branch:23s} "
            f"hidC={sig.hidden_collision_pairs:3d} fold={sig.visible_folds:3d} sum+={sig.sumset_excess:3d}"
        )
    print("  C27 nonunit shell events among q=14 rows:")
    candidates = [s for s in sigs if s.q_threshold == N and s.c27_nonunit_events]
    for sig in sorted(candidates, key=lambda s: (-s.c27_nonunit_events, s.farey_branch, s.label))[:8]:
        print(
            f"    {sig.label:18s} M={str(sig.M):>6s} branch={sig.farey_branch:23s} "
            f"C27_nonunit={sig.c27_nonunit_events:2d} L/U/B/E="
            f"{sig.c27_lower_only}/{sig.c27_upper_only}/{sig.c27_both}/{sig.c27_empty}"
        )


@dataclass(frozen=True)
class Channel:
    name: str
    theorem_scale: int
    branch_separation: int
    sign_visibility: int
    multiplicand_clearance: int
    old_repo_maturity: int
    scalar_decoy_resistance: int

    def score(self) -> tuple[int, int, int, int, int, int]:
        return (
            self.theorem_scale,
            self.branch_separation,
            self.sign_visibility,
            self.multiplicand_clearance,
            self.old_repo_maturity,
            self.scalar_decoy_resistance,
        )


def channel_tournament() -> None:
    print()
    print("[4] Tournament Analysis on relation channels")
    channels = [
        Channel("q/Farey branch", 5, 5, 2, 2, 4, 5),
        Channel("C27 typed shell", 4, 5, 4, 3, 5, 4),
        Channel("visible folds", 2, 3, 5, 2, 4, 4),
        Channel("hidden balanced shells", 1, 2, 4, 1, 4, 3),
        Channel("multiplicand clearance", 3, 4, 2, 5, 5, 4),
        Channel("Kpq incidence", 2, 4, 1, 4, 3, 3),
        Channel("raw sumset/energy", 0, 1, 0, 0, 4, 0),
        Channel("raw runner vertices", 0, 0, 0, 0, 1, 0),
    ]
    wins: dict[int, set[int]] = {i: set() for i in range(len(channels))}
    bitmask = 0
    bit = 0
    raw_sumset_flips = 0
    for i, j in combinations(range(len(channels)), 2):
        ci, cj = channels[i], channels[j]
        if ci.score() >= cj.score():
            winner, loser = i, j
        else:
            winner, loser = j, i
        wins[winner].add(loser)
        if winner == i:
            bitmask |= 1 << bit
        sumset_pref = i if ci.old_repo_maturity >= cj.old_repo_maturity else j
        if sumset_pref != winner:
            raw_sumset_flips += 1
        bit += 1

    fp = s127.tournament_fingerprint(bitmask, len(channels))
    ham = 0
    first_path: tuple[str, ...] | None = None
    for perm in permutations(range(len(channels))):
        if all(perm[t + 1] in wins[perm[t]] for t in range(len(perm) - 1)):
            ham += 1
            if first_path is None:
                first_path = tuple(channels[i].name for i in perm)
    hist = Counter(len(v) for v in wins.values())
    print("  vertices considered/challenged:")
    print("    runners, exact fractions, C27 shells, folds, balanced collisions,")
    print("    divisor blockers, Kpq incidence, raw energy, and proof obligations.")
    print("  chosen vertices: typed relation channels.")
    print("  pair observable: role score")
    print("    (theorem scale, branch separation, sign visibility, multiplicand")
    print("     clearance, old-repo maturity, scalar-decoy resistance)")
    print("  switch/gauge: lexicographically larger role score wins.")
    print(
        f"  fingerprint: score_hist={dict(sorted(hist.items()))} "
        f"c3={fp['c3']} scc={fp['scc']} hp={ham}"
    )
    print(f"  first Hamiltonian path: {' > '.join(first_path or ())}")
    print(f"  old-repo-maturity-only gauge would flip {raw_sumset_flips} pairs.")


def proof_readout() -> None:
    print()
    print("[5] Proof readout")
    print("  The old summand/multiplicand graph does not replace the S130-S133")
    print("  Farey split.  It refines it.")
    print("  Refined proof interface:")
    print("    (1) keep exact q/Farey branch as theorem scale;")
    print("    (2) inside p=2, use C27 typed shell + multiplicand clearance;")
    print("    (3) inside p>=3, use Kpq/K33 incidence + owner packets;")
    print("    (4) inside relation-rich residuals, separate visible folds from")
    print("        hidden balanced shells before applying additive-energy bounds.")
    print("  Proposed next lemma:")
    print("    every remaining non-AP/GW q=14 atom either has a C27 typed-shell")
    print("    defect handled by petal/lift rigidity, or has a three-owner K33")
    print("    incidence packet with enough sign-visible relation mass to feed the")
    print("    HYP-2908 tournament-state lift.")


def main() -> None:
    print("S134 LRC14 BIGRADED SUMMAND x MULTIPLICAND RELATION SIGNATURES")
    print("=" * 78)
    print_selected_signatures()
    sigs = row_bank()
    print_bank_summary(sigs)
    print_extremes(sigs)
    channel_tournament()
    proof_readout()


if __name__ == "__main__":
    main()
