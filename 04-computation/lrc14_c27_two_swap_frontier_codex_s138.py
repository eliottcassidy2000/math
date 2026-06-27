#!/usr/bin/env python3
"""S138: exact C=27 two-swap frontier audit for LRC14.

S136 made the C=27 shell-transfer quotient explicit on the single-replacement
AP/Goddyn-Wong atlas.  This script pushes that quotient one local layer farther:
all primitive rows obtained from AP={1,...,13} by deleting two AP entries and
adding two values <= B.

The main theorem-shaping question is:

    after the rigorous q-threshold prefilter, do any two-swap rows create new
    low-gap atoms at M <= 2/27?

Vertices for Tournament Analysis are deliberately not runners.  The main
binary relation is between proof obligations/channels, with C=27 shell transfers
and exact Farey branch attached as typed carrier data.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
from pathlib import Path
import argparse
import os
import sys


try:
    sys.stdout.reconfigure(line_buffering=True)
except AttributeError:
    pass


REPO = Path(__file__).resolve().parents[1]
N = 14
C = 27
AP = tuple(range(1, 14))
LOW_341 = F(3, 41)
LOW_227 = F(2, 27)


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s124 = load_module(
    "s124_two_swap_frontier",
    REPO / "04-computation" / "lrc14_ap_gw_common_conditions_codex_s124.py",
)


@dataclass(frozen=True)
class Shell:
    a: int
    gcd27: int
    lower: int
    upper: int

    @property
    def total(self) -> int:
        return self.lower + self.upper

    @property
    def state(self) -> str:
        if self.total == 0:
            return "hole"
        if self.total > 2:
            return "multi"
        if self.total == 2:
            return "double"
        if self.upper:
            return "upper"
        return "lower"

    @property
    def detail(self) -> str:
        return f"{self.a}:g{self.gcd27}:{self.lower}{self.upper}"


@dataclass(frozen=True)
class RowAudit:
    speeds: tuple[int, ...]
    label: str
    M: F
    denoms: tuple[int, ...]
    q_threshold: int
    branch: str
    farey_excess: int
    shells: tuple[Shell, ...]
    zero_count: int
    transfer: str
    coarse_transfer: str
    shell_mask: int
    shell_flips_from_ap: int


def exact_m_worker(S: tuple[int, ...]) -> tuple[tuple[int, ...], F, tuple[int, ...]]:
    M, pts = s124.M_exact(S)
    return S, M, tuple(sorted({t.denominator for t in pts}))


def row_label(S: tuple[int, ...]) -> str:
    s = set(S)
    removed = tuple(sorted(set(AP) - s))
    added = tuple(sorted(s - set(AP)))
    if not removed and not added:
        return "AP"
    if removed == (12,) and added == (24,):
        return "GW 12->24"
    if removed == (12,) and added == (36,):
        return "near-miss 12->36"
    if len(removed) == 1 and len(added) == 1:
        return f"swap {removed[0]}->{added[0]}"
    return f"drop{removed}->add{added}"


def farey_branch(M: F) -> tuple[str, int]:
    p, q = M.numerator, M.denominator
    e = N * p - q
    if e == 0:
        return "tight-floor", e
    if e == 1:
        if p == 1:
            return "right-parent-star", e
        if p == 2 and q == C:
            return "C27-petal-two-block", e
        if p >= 3:
            return "K33-unit-excess", e
    if e > 1:
        return f"nonunit-excess-e{e}", e
    return f"below-floor-e{e}", e


def shell_profile(S: tuple[int, ...]) -> tuple[tuple[Shell, ...], int]:
    residues = Counter(v % C for v in S)
    shells = []
    for a in range(1, 14):
        shells.append(
            Shell(
                a=a,
                gcd27=gcd(a, C),
                lower=residues[a],
                upper=residues[C - a],
            )
        )
    return tuple(shells), residues[0]


def transfer_labels(shells: tuple[Shell, ...], zero_count: int) -> tuple[str, str]:
    holes = tuple(s for s in shells if s.total == 0)
    doubles = tuple(s for s in shells if s.total >= 2)
    if not holes and not doubles and zero_count == 0:
        return "perfect-transversal", "perfect"
    h = ",".join(f"{s.a}:g{s.gcd27}" for s in holes) or "-"
    d = ",".join(s.detail for s in doubles) or "-"
    z = str(zero_count) if zero_count else "-"
    detailed = f"H[{h}] D[{d}] Z[{z}]"
    coarse = f"Hg{tuple(s.gcd27 for s in holes) or '-'}->Dg{tuple(s.gcd27 for s in doubles) or '-'};z={zero_count}"
    return detailed, coarse


def shell_priority(shell: Shell) -> tuple[int, int, int, int]:
    status_rank = {
        "hole": 5,
        "multi": 4,
        "double": 4,
        "upper": 2,
        "lower": 1,
    }[shell.state]
    visibility_rank = {1: 3, 3: 2, 9: 1}[shell.gcd27]
    side_rank = shell.upper - shell.lower
    return (status_rank, visibility_rank, side_rank, -shell.a)


def tournament_mask(scores: list[tuple[int, ...]]) -> int:
    mask = 0
    bit = 0
    for i, j in combinations(range(len(scores)), 2):
        if scores[i] >= scores[j]:
            mask |= 1 << bit
        bit += 1
    return mask


def edge(mask: int, n: int, i: int, j: int) -> bool:
    if i == j:
        raise ValueError("loop")
    if i > j:
        return not edge(mask, n, j, i)
    bit = 0
    for a, b in combinations(range(n), 2):
        if a == i and b == j:
            return bool((mask >> bit) & 1)
        bit += 1
    raise AssertionError("unreachable")


def tournament_fingerprint(mask: int, n: int) -> dict[str, object]:
    outdeg = [sum(1 for j in range(n) if j != i and edge(mask, n, i, j)) for i in range(n)]
    c3 = 0
    for a, b, c in combinations(range(n), 3):
        if edge(mask, n, a, b) and edge(mask, n, b, c) and edge(mask, n, c, a):
            c3 += 1
        if edge(mask, n, a, c) and edge(mask, n, c, b) and edge(mask, n, b, a):
            c3 += 1
    scc = strongly_connected_components(mask, n)
    return {
        "score_hist": dict(sorted(Counter(outdeg).items())),
        "c3": c3,
        "scc": tuple(sorted((len(comp) for comp in scc), reverse=True)),
        "hp": hamiltonian_path_count(mask, n),
    }


def strongly_connected_components(mask: int, n: int) -> list[set[int]]:
    graph = [[j for j in range(n) if i != j and edge(mask, n, i, j)] for i in range(n)]
    rgraph = [[j for j in range(n) if i != j and edge(mask, n, j, i)] for i in range(n)]
    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)

    comps: list[set[int]] = []
    seen.clear()
    for v in reversed(order):
        if v in seen:
            continue
        comp = set()
        stack = [v]
        seen.add(v)
        while stack:
            x = stack.pop()
            comp.add(x)
            for w in rgraph[x]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        comps.append(comp)
    return comps


def hamiltonian_path_count(mask: int, n: int) -> int:
    dp: dict[tuple[int, int], int] = {(1 << i, i): 1 for i in range(n)}
    for size in range(1, n):
        nxt: dict[tuple[int, int], int] = defaultdict(int)
        for (used, last), count in dp.items():
            if used.bit_count() != size:
                nxt[(used, last)] += count
                continue
            for v in range(n):
                if used & (1 << v):
                    continue
                if edge(mask, n, last, v):
                    nxt[(used | (1 << v), v)] += count
            nxt[(used, last)] += count
        dp = nxt
    full = (1 << n) - 1
    return sum(count for (used, _), count in dp.items() if used == full)


def audit_row(S: tuple[int, ...], M: F, denoms: tuple[int, ...], ap_mask: int) -> RowAudit:
    branch, e = farey_branch(M)
    shells, zero_count = shell_profile(S)
    transfer, coarse = transfer_labels(shells, zero_count)
    mask = tournament_mask([shell_priority(s) for s in shells])
    flips = sum(
        1
        for i, j in combinations(range(len(shells)), 2)
        if edge(mask, len(shells), i, j) != edge(ap_mask, len(shells), i, j)
    )
    return RowAudit(
        speeds=S,
        label=row_label(S),
        M=M,
        denoms=denoms,
        q_threshold=s124.q_threshold(S),
        branch=branch,
        farey_excess=e,
        shells=shells,
        zero_count=zero_count,
        transfer=transfer,
        coarse_transfer=coarse,
        shell_mask=mask,
        shell_flips_from_ap=flips,
    )


def compute_exact(rows: list[tuple[int, ...]], workers: int, progress_every: int) -> dict[tuple[int, ...], tuple[F, tuple[int, ...]]]:
    if workers <= 1:
        out = {}
        for idx, S in enumerate(rows, 1):
            _, M, denoms = exact_m_worker(S)
            out[S] = (M, denoms)
            if progress_every and idx % progress_every == 0:
                print(f"  exact progress {idx}/{len(rows)}")
        return out

    out: dict[tuple[int, ...], tuple[F, tuple[int, ...]]] = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(exact_m_worker, S): S for S in rows}
        for idx, fut in enumerate(as_completed(futures), 1):
            S, M, denoms = fut.result()
            out[S] = (M, denoms)
            if progress_every and idx % progress_every == 0:
                print(f"  exact progress {idx}/{len(rows)}")
    return out


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge / quotient declaration")
    print("  considered vertices:")
    print("    runners, AP deletion sites, added tails, residues, C27 shell pairs,")
    print("    transfer labels, exact Farey branches, q-threshold filters, and")
    print("    proof obligations.")
    print("  chosen primary vertices:")
    print("    two-swap row obligations, with C27 antipodal shells attached as")
    print("    a typed quotient; this is not a runner tournament.")
    print("  quotient preserves:")
    print("    q-threshold necessity, exact M/Farey branch, unit/nonunit shell holes,")
    print("    and AP-lower-transversal transfer defects.")
    print("  quotient destroys:")
    print("    full off-shell time geometry and large-speed magnitude; exact M(S)")
    print("    is therefore carried alongside every shell label.")
    print("  challenged assumption:")
    print("    single replacements already expose the whole local frontier.  The")
    print("    two-swap bank tests whether interactions of two AP holes manufacture")
    print("    a new low-gap atom.")


def print_bank_readout(limit: int, rows: list[tuple[int, ...]], q14: list[tuple[int, ...]]) -> None:
    print()
    print("[1] Bank and rigorous q-threshold filter")
    print(f"  two-swap primitive rows with added values <= {limit}: {len(rows)}")
    print(f"  rows with q_threshold >= 14: {len(q14)}")
    print("  reason this is lossless for M<=2/27:")
    print("    if q_threshold(S)=q<=13, then t=1/q gives M(S)>=1/q>=1/13>2/27.")
    print("    Hence every row below the second-gap cutoff must pass the q>=14 filter.")


def print_frontier(audits: list[RowAudit], cutoff: F, title: str) -> None:
    front = sorted((a for a in audits if a.M <= cutoff), key=lambda a: (a.M, a.label, a.speeds))
    print()
    print(f"[2] Low-gap frontier {title}: M <= {cutoff}")
    print(f"  rows={len(front)}")
    print(f"  branch_counts={dict(sorted(Counter(a.branch for a in front).items()))}")
    print(f"  transfer_counts={dict(sorted(Counter(a.coarse_transfer for a in front).items()))}")
    for a in front:
        print(
            f"    {a.label:28s} M={str(a.M):>5s} e={a.farey_excess:>2d} "
            f"qthr={a.q_threshold:>2d} denoms={list(a.denoms)!s:14s} "
            f"flips={a.shell_flips_from_ap:>2d} {a.branch:21s} {a.transfer}"
        )


def print_transfer_context(audits: list[RowAudit]) -> None:
    print()
    print("[3] Transfer-label context across q>=14 two-swap rows")
    print(f"  audited exact rows={len(audits)}")
    print(f"  exact branch counts={dict(sorted(Counter(a.branch for a in audits).items()))}")
    low = [a for a in audits if a.M <= LOW_227]
    print(f"  low rows at M<=2/27={len(low)}")
    by_transfer: dict[str, list[RowAudit]] = defaultdict(list)
    for a in audits:
        by_transfer[a.transfer].append(a)
    mixed = []
    for label, group in by_transfer.items():
        branches = {g.branch for g in group}
        if len(branches) > 1 and any(g.M <= LOW_227 for g in group):
            mixed.append((label, group))
    if not mixed:
        print("  no low transfer label is shared with a different Farey branch in this bank.")
    else:
        print("  low transfer labels that also occur in other Farey branches:")
        for label, group in sorted(mixed, key=lambda x: (min(g.M for g in x[1]), x[0]))[:12]:
            branches = Counter(g.branch for g in group)
            examples = ", ".join(
                f"{g.label}:{g.M}" for g in sorted(group, key=lambda x: (x.M, x.label, x.speeds))[:6]
            )
            print(f"    {label:48s} branches={dict(branches)} examples={examples}")
    print()
    print("  most common coarse transfers in q>=14 bank:")
    for label, count in Counter(a.coarse_transfer for a in audits).most_common(10):
        low_count = sum(1 for a in low if a.coarse_transfer == label)
        print(f"    {count:5d} total, {low_count:2d} low  {label}")


def print_new_atom_diagnosis(audits: list[RowAudit]) -> None:
    print()
    print("[4] Two-swap interaction diagnosis")
    low = sorted((a for a in audits if a.M <= LOW_227), key=lambda a: (a.M, a.label, a.speeds))
    singles = [a for a in low if len(set(AP) - set(a.speeds)) <= 1]
    genuine = [a for a in low if len(set(AP) - set(a.speeds)) >= 2]
    print(f"  low rows inherited from AP/single-swap forms: {len(singles)}")
    print(f"  genuine two-hole low rows: {len(genuine)}")
    if genuine:
        print("  genuine two-hole rows:")
        for a in genuine:
            removed = tuple(sorted(set(AP) - set(a.speeds)))
            added = tuple(sorted(set(a.speeds) - set(AP)))
            print(
                f"    remove={removed} add={added} M={a.M} branch={a.branch} "
                f"transfer={a.transfer}"
            )
    else:
        print("  no genuine two-hole row reaches M<=2/27 at this ceiling.")

    print()
    print("  branch-local interpretation:")
    if genuine:
        print("    the frontier is not only the single-swap S136 list; at least one")
        print("    two-hole interaction enters the second-gap neighborhood and needs")
        print("    its own discharge lemma.")
    else:
        print("    the second-gap neighborhood is stable under this two-swap test;")
        print("    local proof effort can focus on the known AP/GW/C27/K33 atoms.")


def proof_channel_tournament() -> None:
    print()
    print("[5] Tournament Analysis: proof-channel relation")
    channels = [
        ("q>=14 threshold prefilter", (5, 5, 4, 5, 5, 4)),
        ("exact M/Farey branch", (5, 5, 5, 5, 4, 5)),
        ("two-hole interaction flag", (4, 5, 4, 4, 5, 4)),
        ("marked C27 shell transfer", (4, 5, 5, 4, 4, 5)),
        ("unit-hole petal discharge", (3, 5, 5, 3, 4, 5)),
        ("gcd3-to-gcd9 K33 packet", (3, 4, 5, 5, 4, 5)),
        ("AP-tail theorem S124/S35", (3, 4, 3, 4, 4, 3)),
        ("raw runner-residue tournament", (1, 2, 2, 1, 1, 1)),
    ]
    scores = [score for _, score in channels]
    mask = tournament_mask(scores)
    fp = tournament_fingerprint(mask, len(channels))
    order = sorted(range(len(channels)), key=lambda i: scores[i], reverse=True)
    print("  vertices: proof channels and obligations, not runners.")
    print("  pair observable:")
    print("    theorem-scale retention, shell predicate retention, interaction")
    print("    visibility, finite-certifiability, state-lift fit, anti-scalar guard.")
    print("  switch/gauge:")
    print("    lexicographically larger role vector wins; ties use listed Hamiltonian path.")
    print(f"  fingerprint: score_hist={fp['score_hist']} c3={fp['c3']} scc={fp['scc']} hp={fp['hp']}")
    print("  Hamiltonian path:")
    print("    " + " > ".join(channels[i][0] for i in order))


def print_proof_readout(limit: int, audits: list[RowAudit]) -> None:
    print()
    print("[6] Proof-route readout")
    low341 = [a for a in audits if a.M <= LOW_341]
    low227 = [a for a in audits if a.M <= LOW_227]
    genuine = [a for a in low227 if len(set(AP) - set(a.speeds)) >= 2]
    print(f"  Through two-swap ceiling {limit}:")
    print(f"    q>=14 exact rows audited: {len(audits)}")
    print(f"    M<=3/41 rows: {len(low341)}")
    print(f"    M<=2/27 rows: {len(low227)}")
    print(f"    genuine two-hole low rows: {len(genuine)}")
    print("  What this buys:")
    print("    the C27 shell-transfer branch is now stress-tested against interacting")
    print("    AP holes, not just one-petal accelerations.")
    print("  Next lemma target:")
    print("    prove that any genuine two-hole row either has a q<14 direct witness,")
    print("    a unit C27 hole discharged by petal/two-block rigidity, or a nonunit")
    print("    gcd3->gcd9 transfer that state-lifts to the K33/HYP-2908 packet.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=40, help="maximum added value in the AP two-swap bank")
    parser.add_argument(
        "--workers",
        type=int,
        default=0,
        help="parallel exact-M workers; 0 means min(cpu_count, 8)",
    )
    parser.add_argument("--progress-every", type=int, default=1000)
    args = parser.parse_args()

    workers = args.workers
    if workers == 0:
        workers = max(1, min(os.cpu_count() or 1, 8))

    print("S138 LRC14 C=27 TWO-SWAP FRONTIER AUDIT")
    print("=" * 78)
    print(f"limit={args.limit}, workers={workers}")
    print_assumption_challenge()

    rows = s124.bank_two_swaps(args.limit)
    q14 = [S for S in rows if s124.q_threshold(S) >= 14]
    print_bank_readout(args.limit, rows, q14)

    ap_shells, _ = shell_profile(AP)
    ap_mask = tournament_mask([shell_priority(s) for s in ap_shells])
    exact = compute_exact(q14, workers=workers, progress_every=args.progress_every)
    audits = [
        audit_row(S, M, denoms, ap_mask)
        for S, (M, denoms) in exact.items()
    ]
    audits.sort(key=lambda a: (a.M, a.label, a.speeds))

    print_frontier(audits, LOW_341, "floor plus first K33 child")
    print_frontier(audits, LOW_227, "second-gap shell")
    print_transfer_context(audits)
    print_new_atom_diagnosis(audits)
    proof_channel_tournament()
    print_proof_readout(args.limit, audits)


if __name__ == "__main__":
    main()
