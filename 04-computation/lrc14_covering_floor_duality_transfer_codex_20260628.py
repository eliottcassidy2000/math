#!/usr/bin/env python3
"""HYP-3424: covering-floor duality transfer scout for LRC14.

This scout is a proof-routing computation, not a numerical LRC solver.  It
looks back at the even/odd, odd/even, add/mult, and chart-change work and asks
which pieces survive the S259 correction, now with HYP-3423's topology-to-
magnitude guardrail, HYP-3422's two-adic interval relocation, and HYP-3421's
off-grid transparency/Rprime branch as floor-facing companions:

    the covering floor is 2-adic; even speeds are the binding obstruction.

The useful duality is therefore not "odd/coprime witnesses prove the floor".
That route was refuted by the 0/400 coprime-transparency test.  The useful
duality is a transfer protocol:

    even fold / 2-adic descent
    -> SPEC decorrelation floor
    -> odd phase-cover debt or finite owner-current sidecar
    -> topology-to-magnitude legality guardrail.

Tournament Analysis uses proof obligations and transfer gates as vertices,
not runners, arcs, residues, or named constants.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from math import comb


OBLIGATIONS = (
    "critical_path_SPEC_floor",
    "two_adic_even_binding",
    "odd_phase_cover_debt",
    "quotient_legality",
    "topology_magnitude_guardrail",
    "finite_owner_sidecar",
    "add_mult_energy_penalty",
    "off_path_7adic_filter",
)


@dataclass(frozen=True)
class Signal:
    name: str
    source: str
    readout: str
    transfer: str
    verdict: str


@dataclass(frozen=True)
class Carrier:
    code: str
    name: str
    keeps: tuple[str, ...]
    destroys: tuple[str, ...]
    obligations: frozenset[str]
    floor_feed: int
    two_adic: int
    odd_resurrection: int
    quotient_legality: int
    finite_sidecar: int
    add_mult: int
    scalar_risk: int
    off_path_risk: int

    @property
    def score(self) -> int:
        return (
            5 * self.floor_feed
            + 4 * self.two_adic
            + 3 * self.odd_resurrection
            + 3 * self.quotient_legality
            + 2 * self.finite_sidecar
            + 2 * self.add_mult
            - 4 * self.scalar_risk
            - 5 * self.off_path_risk
        )


SIGNALS = (
    Signal(
        "coprime_witness_failure",
        "S259/HYP-3418",
        "0/400 covering sets satisfy the naive coprime-to-14 transparency test",
        "odd/coprime data is a phase reservoir, not a standalone witness",
        "reject odd-only witness reduction",
    ),
    Signal(
        "even_binding_obstruction",
        "S259/HYP-3418",
        "binding speeds at the optimum: even=376 versus odd=279; near t=1/2 only 58/300",
        "the floor witness is pulled by even speeds away from the odd-friendly point",
        "promote 2-adic even descent",
    ),
    Signal(
        "two_adic_interval_relocation",
        "HYP-3422",
        "E_safe intersects the union of two odd half-lift filters on all 24 audited covering rows",
        "the even fold should become an interval-overlap lemma rather than a slogan",
        "route the 2-adic branch through finite-ruler/Helly overlap",
    ),
    Signal(
        "topology_magnitude_guardrail",
        "HYP-3423",
        "q-uniform C2/Borsuk-Ulam and C6 topology cannot prove q-specific magnitude floors",
        "topological/Galois data is legal only after a q-specific floor sidecar is retained",
        "block topology-to-magnitude scalarization",
    ),
    Signal(
        "offgrid_rprime_companion",
        "HYP-3421",
        "resonant danger sits on the 14-grid while checked full optima are off-grid; canonical 84m rows feed Rprime",
        "signed SPEC needs the off-grid transparency/Rprime closure, not a raw coprime reduction",
        "route SPEC branch through the all-packet transparency classifier",
    ),
    Signal(
        "even_good_window_wall",
        "S558 even-fold measure reduction",
        "AP and V* have safe slack 0; random samples have min safe slack about 0.10448",
        "even-good windows exist, but odd danger can blanket them at the wall",
        "replace measure counting by positional odd-phase theorem",
    ),
    Signal(
        "frontier_even_owner_label",
        "HYP-3420/HYP-3417/HYP-3419",
        "HYP-3420 owner/chiral classes are exit-exact on scanned banks; HYP-3417 frontier current uses {2:g2,11:g1,13:g1}",
        "finite sidecars already expose one even-cover label plus binding labels",
        "route finite failures to owner-current/Menger sidecar",
    ),
    Signal(
        "add_mult_energy_bridge",
        "HYP-2272/HYP-2128/HYP-2129",
        "8*C(14,2)+1=27^2 and additive pair-sum face is multiplicative energy",
        "additive abundance should become a SPEC or odd-phase debt, not a scalar endpoint",
        "use add/mult only after floor payload is retained",
    ),
    Signal(
        "c3_residue_off_path_filter",
        "HYP-3411/HYP-3413/HYP-3415",
        "C3/Galois/census closes residue/equality structure but not the covering floor",
        "7-adic residue structure is admissible only as a sidecar feeding the floor",
        "demote terminal 7-adic proof attempts",
    ),
)


CARRIERS = (
    Carrier(
        "D00",
        "two_adic_even_floor_descent",
        (
            "even speed self-similarity e=2e'",
            "HYP-3422 E_safe half-lift interval relocation",
            "folded LRC<=13 floor window",
            "binding-speed payload",
        ),
        ("raw odd/coprime witness at t=1/2", "apex-7 census equality story"),
        frozenset(
            {
                "critical_path_SPEC_floor",
                "two_adic_even_binding",
                "odd_phase_cover_debt",
                "off_path_7adic_filter",
            }
        ),
        5,
        5,
        3,
        3,
        1,
        1,
        0,
        0,
    ),
    Carrier(
        "D01",
        "signed_SPEC_decorrelation_floor",
        (
            "|SPEC| < product",
            "14-grid resonance modes",
            "HYP-3421 off-grid/Rprime transparency",
            "Parseval/tail budget",
        ),
        ("finite owner labels unless sidecarred", "raw pointwise witness geometry"),
        frozenset(
            {
                "critical_path_SPEC_floor",
                "two_adic_even_binding",
                "quotient_legality",
                "add_mult_energy_penalty",
            }
        ),
        5,
        4,
        2,
        3,
        0,
        3,
        0,
        0,
    ),
    Carrier(
        "D02",
        "even_good_odd_phase_cover",
        (
            "even-good window G",
            "odd danger arcs restricted to G",
            "pinch/phase positions",
        ),
        ("union-bound total measure", "odd-only lonely point at 1/2"),
        frozenset(
            {
                "two_adic_even_binding",
                "odd_phase_cover_debt",
                "critical_path_SPEC_floor",
            }
        ),
        4,
        4,
        5,
        2,
        0,
        1,
        0,
        0,
    ),
    Carrier(
        "D03",
        "recursive_quotient_sidecar_router",
        (
            "theorem-exit purity",
            "first destroyed coordinate",
            "named terminal debt",
        ),
        ("row order", "untyped scalar equality"),
        frozenset(
            {
                "quotient_legality",
                "topology_magnitude_guardrail",
                "finite_owner_sidecar",
                "odd_phase_cover_debt",
            }
        ),
        3,
        2,
        3,
        5,
        4,
        1,
        0,
        0,
    ),
    Carrier(
        "D04",
        "owner_current_even_cover_sidecar",
        (
            "endpoint-owner current",
            "2:g2 even-cover label",
            "binding labels 11:g1 and 13:g1",
        ),
        ("global SPEC inequality if used terminally",),
        frozenset({"finite_owner_sidecar", "two_adic_even_binding", "quotient_legality"}),
        2,
        4,
        1,
        5,
        5,
        0,
        1,
        0,
    ),
    Carrier(
        "D05",
        "additive_multiplicative_energy_transfer",
        (
            "pair-sum shell",
            "multiplicative energy of roots of unity",
            "odd-square triangular modulus 27",
        ),
        ("2-adic depth", "odd phase position inside G"),
        frozenset({"add_mult_energy_penalty", "critical_path_SPEC_floor"}),
        3,
        1,
        2,
        2,
        0,
        5,
        1,
        0,
    ),
    Carrier(
        "D06",
        "parity_redei_odd_sector_guard",
        (
            "addition certifies multiplication mod 2",
            "self-converse odd worry sector",
            "flip-invariant parity",
        ),
        ("metric floor", "SPEC magnitude", "even binding displacement"),
        frozenset({"odd_phase_cover_debt", "add_mult_energy_penalty"}),
        1,
        1,
        4,
        1,
        0,
        4,
        2,
        0,
    ),
    Carrier(
        "D07",
        "c3_galois_residue_skeleton_filter",
        (
            "unit C3 binding skeleton",
            "Q(sqrt(-7)) transverse sidecar",
            "GW q == 1 mod 3 gate",
        ),
        ("2-adic magnitude floor", "covering optimum displacement"),
        frozenset(
            {
                "off_path_7adic_filter",
                "quotient_legality",
                "topology_magnitude_guardrail",
            }
        ),
        1,
        0,
        1,
        3,
        0,
        1,
        2,
        4,
    ),
    Carrier(
        "D08",
        "raw_7adic_or_scalar_shadow",
        ("apex-7 names", "famous constants", "one-number summaries"),
        ("SPEC floor", "even binding", "odd phase cover", "owner sidecars"),
        frozenset(),
        0,
        0,
        0,
        0,
        0,
        0,
        5,
        5,
    ),
)


def orientation(a: Carrier, b: Carrier) -> tuple[str, str]:
    if (a.score, -int(a.code[1:])) >= (b.score, -int(b.code[1:])):
        return a.code, b.code
    return b.code, a.code


def directed_3cycles(carriers: tuple[Carrier, ...]) -> list[tuple[str, str, str]]:
    wins = {c.code: set() for c in carriers}
    for a, b in combinations(carriers, 2):
        winner, loser = orientation(a, b)
        wins[winner].add(loser)
    cycles = []
    for a, b, c in combinations((x.code for x in carriers), 3):
        if b in wins[a] and c in wins[b] and a in wins[c]:
            cycles.append((a, b, c))
        if c in wins[a] and b in wins[c] and a in wins[b]:
            cycles.append((a, c, b))
    return cycles


def minimal_obligation_covers(carriers: tuple[Carrier, ...]) -> list[tuple[str, ...]]:
    target = set(OBLIGATIONS)
    for size in range(1, len(carriers) + 1):
        winners = []
        for subset in combinations(carriers, size):
            covered = set().union(*(c.obligations for c in subset))
            if target <= covered:
                winners.append(tuple(c.code for c in subset))
        if winners:
            return winners
    return []


def carrier(code: str) -> Carrier:
    return next(c for c in CARRIERS if c.code == code)


def main() -> None:
    print("HYP-3424 COVERING-FLOOR DUALITY TRANSFER SCOUT")
    print("status=SYNTHESIS / proof-routing scout; not an LRC14 proof")
    print("source=HYP-3238/HYP-3234/HYP-2272/HYP-2128/HYP-2129 + S259/HYP-3418 + HYP-3423 + HYP-3422 + HYP-3421")
    print()

    print("## Exact Anchors")
    pairs = comb(14, 2)
    odd_square = 8 * pairs + 1
    print(f"C(14,2)={pairs}")
    print(f"8*C(14,2)+1={odd_square}=27^2")
    print("HYP-3418 binding counts at optima: even=376 odd=279; near_half=58/300")
    print("HYP-3418 coprime-transparency test: 0/400")
    print("S558 even-fold wall: AP/V* safe slack=0; random min safe slack about 0.10448")
    print("HYP-3423 topology guardrail: q-uniform topology cannot close q-specific magnitude floors")
    print("HYP-3422 interval relocation: E_safe half-lift branch overlap holds 24/24 on audited rows")
    print("HYP-3421 off-grid/Rprime branch: resonance transparency feeds the signed SPEC floor")
    print("HYP-3420 owner/chiral sidecar: residue+owner_chiral_class leaves 0 mixed fibers on scanned banks")
    print("HYP-3417 frontier current: {2:g2,11:g1,13:g1}")
    print()

    print("## Duality Transfer Signals")
    for signal in SIGNALS:
        print(f"- {signal.name} [{signal.source}]")
        print(f"  readout: {signal.readout}")
        print(f"  transfer: {signal.transfer}")
        print(f"  verdict: {signal.verdict}")
    print()

    print("## Carrier Scores")
    ordered = tuple(sorted(CARRIERS, key=lambda c: (-c.score, c.code)))
    for c in ordered:
        print(f"{c.code} score={c.score:3d} {c.name}")
        print(f"  covers={tuple(sorted(c.obligations))}")
        print(f"  keeps={c.keeps}")
        print(f"  destroys={c.destroys}")
    print()

    covers = minimal_obligation_covers(CARRIERS)
    print("## Minimal Obligation Covers")
    print(f"obligations={OBLIGATIONS}")
    print(f"minimum_cover_size={len(covers[0]) if covers else 0}")
    for cover in covers[:8]:
        names = tuple(carrier(code).name for code in cover)
        print(f"cover={cover} names={names}")
    print()

    cycles = directed_3cycles(CARRIERS)
    print("## Tournament Analysis")
    print("vertices=proof obligations and transfer gates, not runners/arcs/residues/constants")
    print(
        "pairwise_observable="
        "floor-feed + 2-adic payload + odd resurrection + quotient/topology legality + finite sidecar"
    )
    print("switch=higher weighted proof-facing score; ties by declared code order")
    print(f"vertex_count={len(CARRIERS)}")
    print(f"score_hist={dict(sorted(Counter(c.score for c in CARRIERS).items()))}")
    print(f"directed_3cycles={len(cycles)}")
    path = " -> ".join(c.code for c in ordered)
    print(f"hamiltonian_path={path}")
    print()

    print("## Candidate Lemma Shape")
    print(
        "For every primitive LRC14 covering packet after the q-witness and "
        "LRC<=13 induction split, either:"
    )
    print("  1. the even-speed 2-adic descent supplies a signed-SPEC floor packet with |SPEC| < product;")
    print("  2. HYP-3422 interval relocation supplies an E_safe half-lift overlap;")
    print("  3. the even-good window remains but odd danger arcs emit a positional phase-cover debt;")
    print("  4. the first mixed quotient is discharged by owner-current/Menger sidecar, with even-cover labels named;")
    print("  5. HYP-3421 off-grid transparency routes the signed-SPEC/Rprime branch;")
    print("  6. additive/multiplicative energy concentration converts to a low-SPEC penalty or named odd debt;")
    print("  7. HYP-3423 blocks q-uniform topology/Galois data from closing a q-specific floor;")
    print("  8. the packet is off-path 7-adic/census data and must be demoted unless it feeds one of the above.")
    print()

    print("## Assumption Challenge")
    print(
        "Considered vertices: runners, gaps, fixed circle sections, section boundaries, "
        "wall-crossing events, residues, cover arcs, Fourier modes/SPEC bins, parity "
        "sectors, additive-pair shells, multiplicative-energy atoms, owner currents, "
        "quotient states, and proof obligations."
    )
    print(
        "Chosen vertices: proof obligations and duality-transfer gates.  The quotient "
        "preserves the LRC covering-floor predicate only when it feeds the SPEC floor, "
        "the 2-adic descent, odd phase debt, finite owner sidecars, or an explicit "
        "HYP-3423 legality/off-path filter.  It destroys raw row identity, terminal "
        "7-adic numerology, and scalar shadows."
    )


if __name__ == "__main__":
    main()
