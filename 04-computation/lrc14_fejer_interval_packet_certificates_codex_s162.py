#!/usr/bin/env python3
"""S162: Fejer interval-certificate guardrails for LRC14.

This is a follow-up to HYP-2974/S157.  The S157 full-bank audit found
floating Fejer-vector PSD violations for every positive HYP-2963 packet row,
with AP and Goddyn-Wong as the only zero-safe atoms.  The remaining proof
obligation is to turn those floating trigonometric sums into rigorous
interval-arithmetic certificates attached to labelled packet fibers.

This script does not claim to produce those final certificates.  It computes a
certificate assembly blueprint:

  * the divisor-curried atoms in the reported Fejer sum;
  * a crude but useful precision budget from the negative margin and atom L1;
  * the packet labels that must survive the quotient;
  * the "Robbins no bridge" dependency graph for the certificate.

Robbins' graph theorem is used only as an analogy for proof assembly:
strong orientation is possible exactly when no bridge is forgotten.  Here a
certificate quotient is admissible only when no load-bearing packet label or
interval atom is forgotten.  Robin's divisor theorem is used as the opposite
caution: sigma(n) is a powerful scalar shadow of divisor fibers, but scalar
pushforwards are not allowed to erase route-changing packet data.

Tournament Analysis:
  Vertices are proof obligations/certificate carriers, not runners.  The
  observable is retained LRC predicate data: exact safe component, Fejer
  degree, divisor atom fiber, interval enclosure, packet label, and route
  side-channel.  Ties follow the displayed Hamiltonian path.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import ceil, cos, gcd, log10, log2, pi, sin
from pathlib import Path
import argparse
import re
import sys


REPO = Path(__file__).resolve().parents[1]
S157_PATH = REPO / "04-computation" / "lrc14_fourier_toeplitz_fejer_fullbank_codex_s157.py"
S157_OUT = REPO / "05-knowledge" / "results" / "lrc14_fourier_toeplitz_fejer_fullbank_codex_s157.out"
AP = tuple(range(1, 14))


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s157 = load_module("s162_s157_fejer", S157_PATH)


@dataclass(frozen=True)
class AtomBudget:
    row_name: str
    family: str
    degree: int
    qform: float
    margin: float
    center: Fraction
    width: Fraction
    atoms: int
    nonzero_atoms: int
    atom_l1: float
    raw_weight_l1: float
    max_abs_atom: float
    speed_support: tuple[tuple[int, int], ...]
    invisible_speeds: tuple[int, ...]
    required_bits: int
    required_digits: int
    bridge_flags: tuple[str, ...]


PROOF_VERTICES = (
    "labelled_packet_interval_certificate",
    "Robbins_no_bridge_assembly",
    "endpoint_owner_packet_anchor",
    "Fejer_divisor_atom_bank",
    "Ramanujan_exact_period_projector",
    "Dirichlet_convolution_packet_law",
    "floating_Fejer_shadow",
    "Robin_sigma_scalar_shadow",
    "raw_divisor_count",
)

PROOF_SCORES: dict[str, tuple[int, int, int, int, int, int]] = {
    # predicate retention, kernel control, exactness path, phase retention,
    # packet labels, anti-scalarization.
    "labelled_packet_interval_certificate": (10, 10, 10, 9, 10, 10),
    "Robbins_no_bridge_assembly": (9, 10, 8, 8, 9, 9),
    "endpoint_owner_packet_anchor": (9, 9, 8, 7, 10, 9),
    "Fejer_divisor_atom_bank": (8, 8, 7, 9, 7, 8),
    "Ramanujan_exact_period_projector": (7, 8, 8, 9, 7, 8),
    "Dirichlet_convolution_packet_law": (6, 8, 8, 5, 6, 7),
    "floating_Fejer_shadow": (8, 4, 3, 8, 4, 5),
    "Robin_sigma_scalar_shadow": (4, 3, 7, 3, 2, 2),
    "raw_divisor_count": (2, 1, 6, 1, 1, 1),
}


def fmt_fraction(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def fmt_float(x: float) -> str:
    if abs(x) < 1e-3:
        return f"{x:.6e}"
    return f"{x:.9f}"


def replace_many(holes: tuple[int, ...], adds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted((set(AP) - set(holes)) | set(adds)))


def hard_rows() -> list:
    """Rows selected from the S157 hardest-certificate table."""
    Row = s157.Row
    return [
        Row("single swap 6->63", "single-swap hardest", replace_many((6,), (63,))),
        Row("single swap 6->34", "single-swap hard", replace_many((6,), (34,))),
        Row("single swap 6->86", "single-swap hard", replace_many((6,), (86,))),
        Row("two drop(6,9)->add(18,23)", "two-swap hard", replace_many((6, 9), (18, 23))),
        Row("two drop(6,10)->add(30,34)", "two-swap hard", replace_many((6, 10), (30, 34))),
        Row("two drop(6,10)->add(30,35)", "two-swap hard", replace_many((6, 10), (30, 35))),
        Row("two drop(6,10)->add(20,34)", "two-swap hard", replace_many((6, 10), (20, 34))),
        Row("two drop(10,12)->add(24,30)", "two-swap hard", replace_many((10, 12), (24, 30))),
        Row("two drop(12,13)->add(14,29)", "smallest-margin bank row", replace_many((12, 13), (14, 29))),
    ]


def selected_rows() -> list:
    rows = list(s157.named_rows()) + hard_rows()
    seen: set[tuple[int, ...]] = set()
    out = []
    for row in rows:
        if row.speeds in seen:
            continue
        seen.add(row.speeds)
        out.append(row)
    return out


def first_certificate(row, max_degree: int, negative_tol: float):
    ghat = s157.ghat_table(max_degree)
    return s157.first_fejer_certificate(row, max_degree, ghat, negative_tol)


def divisor_atom_budget(audit) -> AtomBudget | None:
    cert = audit.certificate
    if cert is None:
        return None
    degree = cert.degree
    x = float(cert.center)

    atoms = 0
    nonzero_atoms = 0
    atom_l1 = 0.0
    raw_weight_l1 = 0.0
    max_abs_atom = 0.0
    speed_counts: Counter[int] = Counter()

    for k in range(1, degree + 1):
        fejer_weight = 2.0 * (1.0 - k / (degree + 1.0))
        phase = cos(2.0 * pi * k * x)
        for v in audit.row.speeds:
            if k % v:
                continue
            m = k // v
            atoms += 1
            speed_counts[v] += 1
            raw_weight = abs(fejer_weight / (pi * m))
            raw_weight_l1 += raw_weight
            term = fejer_weight * phase * sin(pi * m / 7.0) / (pi * m)
            if term:
                nonzero_atoms += 1
            abs_term = abs(term)
            atom_l1 += abs_term
            max_abs_atom = max(max_abs_atom, abs_term)

    margin = -cert.qform
    if margin <= 0 or atom_l1 <= 0:
        required_bits = 0
        required_digits = 0
    else:
        # If every atom is interval-enclosed to relative error 2^-p, then
        # L1*2^-p <= margin/8 leaves room for pi and outward rounding slack.
        required_bits = max(0, ceil(log2(8.0 * atom_l1 / margin)))
        required_digits = max(0, ceil(required_bits * log10(2.0)))

    bridge_flags = []
    if not audit.row.family:
        bridge_flags.append("missing_family_label")
    if cert.center.denominator == 1:
        bridge_flags.append("degenerate_center")
    if audit.invisible_speeds_at_first:
        bridge_flags.append("degree_blind_speed_tail")
    if margin < 1e-6:
        bridge_flags.append("tiny_margin_requires_family_interval")
    if degree > 150:
        bridge_flags.append("high_degree_requires_atom_bank")

    return AtomBudget(
        row_name=audit.row.name,
        family=audit.row.family,
        degree=degree,
        qform=cert.qform,
        margin=margin,
        center=cert.center,
        width=cert.width,
        atoms=atoms,
        nonzero_atoms=nonzero_atoms,
        atom_l1=atom_l1,
        raw_weight_l1=raw_weight_l1,
        max_abs_atom=max_abs_atom,
        speed_support=tuple(sorted(speed_counts.items())),
        invisible_speeds=audit.invisible_speeds_at_first,
        required_bits=required_bits,
        required_digits=required_digits,
        bridge_flags=tuple(bridge_flags),
    )


def prior_full_bank_margin() -> tuple[str, int, float] | None:
    if not S157_OUT.exists():
        return None
    text = S157_OUT.read_text(encoding="utf-8", errors="replace")
    pattern = re.compile(
        r"smallest_negative_margin=([\-0-9.eE]+)\s+at\s+(.+?), degree=(\d+)"
    )
    matches = pattern.findall(text)
    if not matches:
        return None
    qform, row_name, degree = matches[-1]
    return row_name.strip(), int(degree), float(qform)


def tournament_fingerprint() -> dict[str, object]:
    tie_rank = {name: idx for idx, name in enumerate(PROOF_VERTICES)}
    n = len(PROOF_VERTICES)

    def edge(i: int, j: int) -> bool:
        a = PROOF_VERTICES[i]
        b = PROOF_VERTICES[j]
        av = PROOF_SCORES[a]
        bv = PROOF_SCORES[b]
        awins = sum(x > y for x, y in zip(av, bv))
        bwins = sum(y > x for x, y in zip(av, bv))
        if awins != bwins:
            return awins > bwins
        return tie_rank[a] < tie_rank[b]

    outdeg = [sum(1 for j in range(n) if i != j and edge(i, j)) for i in range(n)]
    c3 = 0
    for a, b, c in combinations(range(n), 3):
        if edge(a, b) and edge(b, c) and edge(c, a):
            c3 += 1
        if edge(a, c) and edge(c, b) and edge(b, a):
            c3 += 1
    return {
        "score_hist": dict(sorted(Counter(outdeg).items())),
        "directed_3_cycles": c3,
        "hamiltonian_path": PROOF_VERTICES,
    }


def print_header(args) -> None:
    print("S162 LRC14 FEJER INTERVAL / ROBBINS GUARDRAILS")
    print("=" * 78)
    print(f"max_degree={args.max_degree} negative_tol={args.negative_tol}")
    print()
    print("[0] Assumption challenge")
    print("  Candidate vertices considered:")
    print("    runners, speeds, residues, endpoints, safe components, Fourier modes,")
    print("    divisor atoms, interval enclosures, packet labels, and proof obligations.")
    print("  Chosen vertices:")
    print("    certificate obligations and labelled divisor-curried Fejer atoms.")
    print("  Preserved LRC predicate:")
    print("    covering implies F_S=C_S-1>=0; a negative Toeplitz/Fejer form")
    print("    contradicts covering and certifies a strict lonely interval.")
    print("  Destroyed unless reattached:")
    print("    endpoint-owner identities, C27/K33 route labels, and exact packet family.")
    print()


def print_budget_table(budgets: list[AtomBudget]) -> None:
    print("[1] Interval precision budget for selected hard packets")
    print(
        f"  {'row':36s} {'deg':>4s} {'margin':>12s} {'atoms':>6s} "
        f"{'L1':>10s} {'bits':>5s} {'digits':>6s} {'center':>12s} flags"
    )
    for b in budgets:
        flags = ",".join(b.bridge_flags) if b.bridge_flags else "-"
        print(
            f"  {b.row_name[:36]:36s} {b.degree:4d} {fmt_float(b.margin):>12s} "
            f"{b.atoms:6d} {fmt_float(b.atom_l1):>10s} {b.required_bits:5d} "
            f"{b.required_digits:6d} {fmt_fraction(b.center):>12s} {flags}"
        )
    print()


def print_worst_cases(budgets: list[AtomBudget]) -> None:
    print("[2] Worst selected cases by different notions")
    for label, key in (
        ("highest degree", lambda b: b.degree),
        ("smallest margin", lambda b: -b.margin),
        ("largest atom bank", lambda b: b.atoms),
        ("largest required precision", lambda b: b.required_bits),
    ):
        b = max(budgets, key=key)
        print(
            f"  {label:26s}: {b.row_name} degree={b.degree} "
            f"margin={fmt_float(b.margin)} atoms={b.atoms} "
            f"bits={b.required_bits} digits={b.required_digits}"
        )
    prior = prior_full_bank_margin()
    if prior:
        row_name, degree, qform = prior
        print(
            "  prior S157 full-bank weakest margin: "
            f"{row_name}, degree={degree}, qform={fmt_float(qform)}"
        )
    print()


def print_atom_decomposition(budgets: list[AtomBudget], limit: int) -> None:
    print(f"[3] Divisor-fiber atom support for top {limit} precision cases")
    ordered = sorted(budgets, key=lambda b: (b.required_bits, b.degree), reverse=True)
    for b in ordered[:limit]:
        support = " ".join(f"{v}:{c}" for v, c in b.speed_support)
        print(f"  {b.row_name}")
        print(
            f"    family={b.family} degree={b.degree} width={fmt_fraction(b.width)} "
            f"invisible={b.invisible_speeds or '-'}"
        )
        print(
            f"    support_by_speed={support} max_abs_atom={fmt_float(b.max_abs_atom)} "
            f"raw_weight_L1={fmt_float(b.raw_weight_l1)}"
        )
    print()


def print_robbins_guardrail(budgets: list[AtomBudget]) -> None:
    tiny = [b.row_name for b in budgets if "tiny_margin_requires_family_interval" in b.bridge_flags]
    high = [b.row_name for b in budgets if "high_degree_requires_atom_bank" in b.bridge_flags]
    tail = [b.row_name for b in budgets if "degree_blind_speed_tail" in b.bridge_flags]

    print("[4] Robbins/Robin readout")
    print("  Robbins graph analogy:")
    print("    A strong orientation fails exactly at a bridge.  The proof-certificate")
    print("    analogue is: do not quotient away a load-bearing interval atom or")
    print("    packet label.  The certificate dependency graph has bridges:")
    print("      exact rational center -> divisor atom bank -> trig interval enclosure")
    print("      -> signed margin -> labelled packet fiber -> route handoff.")
    print("  Bridge alerts in this selected packet set:")
    print(f"    tiny margins needing row-family intervals: {tiny or '-'}")
    print(f"    high-degree atom-bank dependencies: {high or '-'}")
    print(f"    degree-blind speed tails: {tail or '-'}")
    print("  Robin divisor-function caution:")
    print("    sigma/tau are useful scalar pushforwards of divisor fibers, but the")
    print("    LRC proof cannot trust them unless Ramanujan, endpoint-owner, or")
    print("    interval-certificate side channels control the quotient kernel.")
    print()


def print_theorem_target() -> None:
    print("[5] Certificate theorem target")
    print("  For each labelled packet family P(S), prove one of:")
    print("    1. AP/GW boundary atom: zero strict mass and zero-sum endpoint current;")
    print("    2. exact interval Fejer certificate: Q_d(x)<0 with rational outward")
    print("       enclosures for all divisor-curried atoms;")
    print("    3. Ramanujan/Toeplitz handoff: exact-period projector forces a late")
    print("       q packet, then interval-enclose only that packet;")
    print("    4. state-lift handoff: the unresolved bridge creates the HYP-2908/")
    print("       THM-572 forbidden tournament-conflict packet.")
    print("  The important simplification is familywise certification: the interval")
    print("  object should be attached to packet fibers, not to 21911 unrelated rows.")
    print()


def print_tournament_analysis() -> None:
    fp = tournament_fingerprint()
    print("[6] Tournament Analysis")
    print("  vertices:")
    print("    " + " > ".join(PROOF_VERTICES))
    print("  pair observable:")
    print("    retention of cover predicate, kernel control, exactness path, phase,")
    print("    packet labels, and anti-scalarization.")
    print("  switch/gauge:")
    print("    majority over the six coordinates; ties follow the displayed path.")
    print(
        "  fingerprint: "
        f"score_hist={fp['score_hist']} "
        f"directed_3_cycles={fp['directed_3_cycles']} "
        "hamiltonian_paths=1"
    )
    print()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-degree", type=int, default=512)
    parser.add_argument("--negative-tol", type=float, default=1e-8)
    parser.add_argument("--atom-support-limit", type=int, default=6)
    args = parser.parse_args()

    print_header(args)
    budgets: list[AtomBudget] = []
    zero_rows: list[str] = []
    misses: list[str] = []
    for row in selected_rows():
        audit = first_certificate(row, args.max_degree, args.negative_tol)
        if audit.safe_measure == 0:
            zero_rows.append(row.name)
            continue
        budget = divisor_atom_budget(audit)
        if budget is None:
            misses.append(row.name)
        else:
            budgets.append(budget)

    print(f"[summary] selected rows={len(selected_rows())}")
    print(f"  zero-safe rows={zero_rows or '-'}")
    print(f"  positive rows without certificate at cap={misses or '-'}")
    print(f"  positive rows budgeted={len(budgets)}")
    print()

    print_budget_table(budgets)
    print_worst_cases(budgets)
    print_atom_decomposition(budgets, args.atom_support_limit)
    print_robbins_guardrail(budgets)
    print_theorem_target()
    print_tournament_analysis()


if __name__ == "__main__":
    main()
