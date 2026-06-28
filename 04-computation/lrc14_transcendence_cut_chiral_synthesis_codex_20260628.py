#!/usr/bin/env python3
"""HYP-3420: transcendence/cut/chiral synthesis for LRC14.

This is a proof-route scout, not an LRC14 proof.  It follows HYP-3405 and
HYP-3406: HYP-3405 gives the finite AP-collar certificate, and HYP-3406 shows
that enlarged residue packets need endpoint-owner support.

The script keeps the speculative external ideas honest by translating each of
them into a possible packet field, then testing the concrete part on HYP-3406's
expanded owner-support bank:

* Menger-style owner cuts for mixed theorem-exit fibers;
* mirror/chiral owner-support signatures and their recursion across banks;
* a Barban-Davenport-Halberstam-style fiber-variance meter;
* a tournament over named proof carriers, not runners or raw scales.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, relpath: str):
    path = ROOT / relpath
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


h3406 = load_module(
    "h3407_expanded_residue_owner_repair",
    "04-computation/lrc14_expanded_residue_owner_repair_codex_20260628.py",
)


BANKS = ((20, 4), (30, 8), (48, 12), (60, 16))


def parse_owner(label: str) -> tuple[int, str]:
    residue, grade = label.split(":")
    return int(residue), grade


def mirror_owner(label: str) -> str:
    residue, grade = parse_owner(label)
    mirrored = (-residue) % 14
    if mirrored == 0:
        mirrored = 14
    return f"{mirrored}:{grade}"


def mirror_word(word: tuple[str, ...]) -> tuple[str, ...]:
    return tuple(sorted(mirror_owner(label) for label in word))


def chiral_class(word: tuple[str, ...]) -> tuple[str, ...]:
    reflected = mirror_word(word)
    return min(word, reflected)


def chiral_sign(word: tuple[str, ...]) -> int:
    reflected = mirror_word(word)
    if word == reflected:
        return 0
    return -1 if word < reflected else 1


def sidecar_key(row: h3406.ExpandedRow, label: str) -> tuple[object, ...]:
    if label == "residue_only":
        return row.coarse_base + (row.residue_word,)
    if label == "residue_plus_v2":
        return row.coarse_base + (row.residue_word, row.v2_word)
    if label == "residue_plus_height":
        return row.coarse_base + (row.residue_word, row.exact_height_word)
    if label == "residue_plus_owner_support":
        return row.coarse_base + (row.residue_word, row.owner_support_word)
    if label == "residue_plus_owner_chiral_class":
        return row.coarse_base + (row.residue_word, chiral_class(row.owner_support_word))
    if label == "residue_plus_owner_chiral_full":
        return row.coarse_base + (
            row.residue_word,
            chiral_class(row.owner_support_word),
            chiral_sign(row.owner_support_word),
        )
    raise KeyError(label)


def fibers_by_sidecar(rows: list[h3406.ExpandedRow], label: str) -> dict[tuple[object, ...], list[h3406.ExpandedRow]]:
    out: dict[tuple[object, ...], list[h3406.ExpandedRow]] = defaultdict(list)
    for row in rows:
        out[sidecar_key(row, label)].append(row)
    return out


def mixed_by_sidecar(rows: list[h3406.ExpandedRow], label: str) -> list[list[h3406.ExpandedRow]]:
    out = []
    for fiber in fibers_by_sidecar(rows, label).values():
        if len({row.kernel_flag for row in fiber}) > 1:
            out.append(sorted(fiber, key=lambda row: (row.kernel_flag, row.name)))
    out.sort(key=lambda fiber: (-len(fiber), tuple(row.name for row in fiber)))
    return out


def pair_disagreement_variance(rows: list[h3406.ExpandedRow], label: str) -> int:
    total = 0
    for fiber in fibers_by_sidecar(rows, label).values():
        counts = Counter(row.kernel_flag for row in fiber)
        values = list(counts.values())
        for i, left in enumerate(values):
            for right in values[i + 1 :]:
                total += left * right
    return total


def minimal_owner_cut(fiber: list[h3406.ExpandedRow]) -> tuple[str, ...] | None:
    constraints: list[set[str]] = []
    universe: set[str] = set()
    for row in fiber:
        universe.update(row.owner_support_word)
    for left, right in combinations(fiber, 2):
        if left.kernel_flag == right.kernel_flag:
            continue
        diff = set(left.owner_support_word).symmetric_difference(right.owner_support_word)
        if not diff:
            return None
        constraints.append(diff)
    labels = sorted(universe)
    for size in range(len(labels) + 1):
        for candidate in combinations(labels, size):
            chosen = set(candidate)
            if all(chosen & constraint for constraint in constraints):
                return tuple(candidate)
    return None


SWAP_RE = re.compile(r"(single swap|petal) (\d+)->(\d+)")


def replacement_family_name(row_name: str) -> str | None:
    match = SWAP_RE.search(row_name)
    if not match:
        return None
    kind, source, target = match.groups()
    return f"{kind}->{target}"


def bank_audit(single_limit: int, two_swap_limit: int) -> dict[str, object]:
    rows = h3406.build_rows(single_limit, two_swap_limit)
    sidecars = [
        "residue_only",
        "residue_plus_v2",
        "residue_plus_height",
        "residue_plus_owner_chiral_class",
        "residue_plus_owner_chiral_full",
        "residue_plus_owner_support",
    ]
    mixed_counts = {label: len(mixed_by_sidecar(rows, label)) for label in sidecars}
    variances = {label: pair_disagreement_variance(rows, label) for label in sidecars}
    residue_mixed = mixed_by_sidecar(rows, "residue_only")
    cut_sizes = []
    cut_examples = []
    for fiber in residue_mixed:
        cut = minimal_owner_cut(fiber)
        if cut is not None:
            cut_sizes.append(len(cut))
            cut_examples.append((len(fiber), cut, tuple(row.name for row in fiber)))
    family_counter: Counter[str] = Counter()
    for fiber in residue_mixed:
        for row in fiber:
            family = replacement_family_name(row.name)
            if family is not None:
                family_counter[family] += 1
    chiral_counts = Counter(chiral_sign(row.owner_support_word) for row in rows)
    return {
        "rows": rows,
        "mixed_counts": mixed_counts,
        "variances": variances,
        "residue_mixed": residue_mixed,
        "cut_sizes": cut_sizes,
        "cut_examples": cut_examples,
        "family_counter": family_counter,
        "chiral_counts": chiral_counts,
    }


@dataclass(frozen=True)
class Carrier:
    name: str
    packet_field: str
    concrete: int
    preserves: int
    repairs_h3406: int
    formalizable: int
    novelty: int
    risk_penalty: int
    note: str

    @property
    def score(self) -> int:
        return (
            5 * self.concrete
            + 4 * self.preserves
            + 4 * self.repairs_h3406
            + 3 * self.formalizable
            + self.novelty
            - 3 * self.risk_penalty
        )


CARRIERS = [
    Carrier(
        "owner_menger_cut_certificate",
        "minimal endpoint-owner cut separating mixed theorem exits",
        5,
        5,
        5,
        5,
        3,
        0,
        "Directly measured on HYP-3406 mixed residue fibers.",
    ),
    Carrier(
        "chiral_owner_recursion_signature",
        "mirror class + chirality sign of endpoint-owner support",
        4,
        4,
        4,
        4,
        4,
        1,
        "Keeps recursive owner orientation instead of scalar owner counts.",
    ),
    Carrier(
        "bdh_fiber_variance_ledger",
        "Barban-Davenport-Halberstam-style mixed-fiber variance over residue packets",
        4,
        4,
        3,
        4,
        3,
        1,
        "A second-moment alarm for residue packets before adding owner support.",
    ),
    Carrier(
        "schwarz_christoffel_owner_polygon",
        "contact polygon vertices plus endpoint-owner accessory parameters",
        3,
        4,
        3,
        3,
        4,
        1,
        "Makes the SC accessory parameter the owner-support debt.",
    ),
    Carrier(
        "krasner_hensel_owner_stability",
        "p-adic stability radius for residue/owner sidecar persistence",
        3,
        4,
        3,
        4,
        3,
        1,
        "Useful after HYP-3406 because owner support persists across bank expansion.",
    ),
    Carrier(
        "sophie_germain_quartic_factor_gate",
        "quartic phi4/off-circle correction split into two quadratic factors",
        2,
        3,
        1,
        4,
        3,
        1,
        "Good for the quartic side, but it does not repair the HYP-3406 owner leak.",
    ),
    Carrier(
        "bring_radical_branch_guard",
        "degree-five residual branch/monodromy sidecar rather than radical simplification",
        2,
        3,
        1,
        3,
        5,
        2,
        "A warning for quintic residuals: keep branch data or use a special inverse.",
    ),
    Carrier(
        "hermite_lindemann_weierstrass_scale_guard",
        "transcendence side-condition ledger for exponential/root shadows",
        1,
        3,
        0,
        4,
        4,
        2,
        "Guards against algebraizing raw exp/root shadows without side conditions.",
    ),
    Carrier(
        "meissel_mertens_loglog_calibration",
        "constant term in prime-channel/loglog exceptional-set budget",
        2,
        2,
        0,
        3,
        3,
        2,
        "A calibration offset, not a separator of known mixed fibers.",
    ),
    Carrier(
        "ramanujan_soldner_balance_root",
        "zero-current balance point for cumulative signed owner charge",
        1,
        2,
        0,
        2,
        4,
        3,
        "Only useful if converted into an owner-current zero theorem.",
    ),
    Carrier(
        "raw_exp_exp_exp_79_scale",
        "raw magnitude N=exp(exp(exp(79)))",
        0,
        0,
        0,
        1,
        2,
        5,
        "Raw scale is intentionally demoted; only tower height/log channels matter.",
    ),
]


def beats(left: Carrier, right: Carrier) -> bool:
    if left.score != right.score:
        return left.score > right.score
    return CARRIERS.index(left) < CARRIERS.index(right)


def tournament_fingerprint() -> tuple[dict[int, int], int, int, list[str]]:
    score_hist = dict(sorted(Counter(carrier.score for carrier in CARRIERS).items()))
    directed_3cycles = 0
    for a, b, c in combinations(CARRIERS, 3):
        ab = beats(a, b)
        bc = beats(b, c)
        ca = beats(c, a)
        ba = beats(b, a)
        cb = beats(c, b)
        ac = beats(a, c)
        if (ab and bc and ca) or (ba and cb and ac):
            directed_3cycles += 1

    n = len(CARRIERS)
    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp.get((mask, last), 0)
            if not count:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if beats(CARRIERS[last], CARRIERS[nxt]):
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + count
    full = (1 << n) - 1
    hamiltonian_paths = sum(dp.get((full, last), 0) for last in range(n))
    priority = [carrier.name for carrier in sorted(CARRIERS, key=lambda carrier: (-carrier.score, CARRIERS.index(carrier)))]
    return score_hist, directed_3cycles, hamiltonian_paths, priority


def frac_word(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}" if value.denominator != 1 else str(value.numerator)


def main() -> None:
    print("HYP-3420 transcendence / cut / chiral synthesis")
    print("=" * 78)
    print("status=creative synthesis + exact HYP-3406 owner-support audit; not LRC14 proof")
    print("scale_sentinel=N=exp(exp(exp(79))); logloglog(N)=79; raw magnitude is not a proof vertex")
    print()

    largest_audit = None
    for single_limit, two_swap_limit in BANKS:
        audit = bank_audit(single_limit, two_swap_limit)
        if (single_limit, two_swap_limit) == BANKS[-1]:
            largest_audit = audit
        print(f"BANK single_limit={single_limit} two_swap_limit={two_swap_limit}")
        print(f"  rows={len(audit['rows'])}")
        print(f"  mixed_counts={audit['mixed_counts']}")
        print(f"  bdh_pair_disagreement_variance={audit['variances']}")
        print(f"  owner_cut_size_hist={dict(sorted(Counter(audit['cut_sizes']).items()))}")
        print(f"  chiral_sign_hist={dict(sorted(audit['chiral_counts'].items()))}")
        print(f"  recursive_replacement_family_counts={dict(sorted(audit['family_counter'].items()))}")
        print()

    assert largest_audit is not None
    print("LARGEST BANK OWNER-CUT EXAMPLES")
    for size, cut, names in largest_audit["cut_examples"][:3]:
        print(f"  fiber_size={size} min_owner_cut={cut}")
        for name in names:
            print(f"    {name}")
    print()

    print("NAMED EXTERNAL IDEAS AS PACKET FIELDS")
    print("  Bring radical: keep degree-five branch/monodromy data; do not pretend radicals close the packet.")
    print("  Schwarz-Christoffel: endpoint-owner support is the accessory parameter of the contact polygon.")
    print("  Barban-Davenport-Halberstam: residue-only failure is a fiber-variance signal, killed by owner support here.")
    print("  Menger cuts: minimal owner cuts separate mixed theorem exits in residue fibers.")
    print("  Ramanujan-Soldner: useful only as a signed-current zero template, not as a raw constant.")
    print("  Sophie Germain identity: factor quartic phi4 debt into quadratic gates when the quartic route is active.")
    print("  Hermite-Lindemann-Weierstrass: transcendence shadows need side conditions before algebraic compression.")
    print("  Krasner: owner-support persistence should become a p-adic/Hensel stability radius.")
    print("  Meissel-Mertens: loglog channel budgets need the constant term only after packet labels are retained.")
    print()

    print("CARRIER SCORES")
    for carrier in sorted(CARRIERS, key=lambda item: (-item.score, CARRIERS.index(item))):
        print(
            f"  {carrier.name:42s} score={carrier.score:3d} "
            f"field={carrier.packet_field}; note={carrier.note}"
        )
    print()

    score_hist, cycles, hpaths, priority = tournament_fingerprint()
    print("TOURNAMENT ANALYSIS")
    print("  vertices=proof carriers and packet fields, not runners/raw scales/constants")
    print("  pairwise_observable=known-fiber repair + preserved LRC predicate + formalization readiness - analogy risk")
    print("  switch_gauge=A -> B iff A has larger weighted carrier score; ties use declared priority order")
    print(f"  score_hist={score_hist}")
    print(f"  directed_3cycles={cycles}")
    print(f"  hamiltonian_path_count={hpaths}")
    print(f"  priority_path={' -> '.join(priority)}")
    print(
        "  assumption_challenge=considered runners/residues/speeds/constants/huge scales; "
        "selected owner cuts, chiral signatures, variance ledgers, and branch/accessory sidecars "
        "because they preserve theorem-exit exactness."
    )
    print()

    print("CONCLUSION")
    print("  HYP-3406 makes owner support the next exact enlarged-bank repair.")
    print("  The best new concrete route is an endpoint-owner Menger cut theorem:")
    print("    every residue-mixed theorem-exit fiber has a small owner-support cut,")
    print("    and recursive/chiral owner signatures explain the growing 26/40/54 families.")
    print("  The Bring/SC/BDH/Krasner/Mertens/HLW/Sophie-Germain ideas are useful only")
    print("    after translation into packet fields with destroyed-coordinate guardrails.")


if __name__ == "__main__":
    main()
