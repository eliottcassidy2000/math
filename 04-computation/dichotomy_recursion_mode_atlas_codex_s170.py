#!/usr/bin/env python3
"""Dichotomy / recursion mode atlas.

This is a synthesis scout, not a proof engine.  It compares recurring binary
and operator dichotomies in the repository as proof carriers for LRC-style
quotient guardrails.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from functools import lru_cache


@dataclass(frozen=True)
class Mode:
    name: str
    dichotomy: str
    recursion_law: str
    preserves: str
    destroys_if_scalarized: str
    anchors: str
    features: tuple[int, ...]


FEATURES = (
    "retained_address",
    "quotient_guardrail",
    "lrc_transfer",
    "recursion_depth_control",
    "side_channel_exposure",
    "proof_maturity",
    "unification_breadth",
)


MODES = [
    Mode(
        "sign_cut",
        "positive/negative",
        "epsilon_i in {+1,-1}; opposite signs turn differences into sums",
        "observer loneliness M is invariant, while pair clocks expose cut edges",
        "pair-sum shell partners and zero clocks",
        "signed-lrc-theory-sign-is-a-cut-and-the-worryset-splits-s699.md; T764",
        (5, 5, 5, 3, 5, 5, 4),
    ),
    Mode(
        "parity_fold",
        "odd/even",
        "odd is unit-propagating; even is the first 2-adic projection defect",
        "2-adic depth, unit/nonunit shell split, clean/wall geometry",
        "CRT block where doubling stops being an automorphism",
        "lrc-parity-ladder-proof-program-s576.md; lrc-orbit-functor-rigidity-s590.md",
        (5, 5, 5, 4, 4, 4, 5),
    ),
    Mode(
        "additive_pair_sum_face",
        "addition/multiplication",
        "a+b=D creates denominator pincers and shell-partner clocks",
        "endpoint owners, pair-sum denominators, additive witness resolution",
        "unit-orbit product structure and local factor labels",
        "THM-401; lrc-the-pair-sum-sieve-modulus-is-2n-1-an-identity-s590.md",
        (5, 5, 5, 3, 5, 5, 5),
    ),
    Mode(
        "multiplicative_unit_orbit",
        "multiplication/addition",
        "unit multiplication transports witnesses inside a residue orbit",
        "unit-shell visibility, CRT localization, orbit functor action",
        "observer-coupled additive folds and endpoint labels",
        "lrc-orbit-functor-rigidity-s590.md; HYP-3001",
        (5, 4, 5, 4, 3, 4, 5),
    ),
    Mode(
        "dyadic_doubling_tree",
        "*2/+2",
        "n -> 2n is vertical descent in the 2-adic tree",
        "valuation depth, shell descent, phase plateau after first gate",
        "horizontal archimedean ordering and plus-two ruler weave",
        "the-line-and-the-trees-why-N-has-plus2-and-times2-the-adelic-treebolic-geometry-s547o.md; lrc-h-loneliness-recursion-s505.md",
        (4, 4, 5, 5, 4, 4, 5),
    ),
    Mode(
        "plus_two_line_motion",
        "+2/*2",
        "n -> n+2 is horocyclic line motion weaving through 2-adic levels",
        "archimedean order, successor shell, additive progression",
        "clean valuation descent and multiplicative orbit closure",
        "the-line-and-the-trees-why-N-has-plus2-and-times2-the-adelic-treebolic-geometry-s547o.md",
        (3, 3, 3, 3, 3, 3, 5),
    ),
    Mode(
        "farey_sum_affine_check",
        "sum/product/power",
        "(p+q)/q = 1+p/q preserves the true Farey order",
        "exact M=p/q as an affine theorem-safe check",
        "branch family after the unit-excess gate",
        "HYP-3001; HYP-3000; HYP-2999",
        (5, 5, 4, 3, 3, 4, 4),
    ),
    Mode(
        "farey_product_scheduler",
        "product/sum",
        "after e=14p-q=1, (p*q)/q=p routes p=1,2,>=3 packet families",
        "unit-excess branch scheduler for q-parent, C27, K33/state-lift",
        "global order; product-value has many row-bank flips if used alone",
        "HYP-3001-lrc14-farey-mutation-certificate-scheduler.md",
        (5, 5, 5, 4, 4, 4, 4),
    ),
    Mode(
        "zeckendorf_path_normal_form",
        "abundance/uniqueness",
        "no adjacent Fibonacci carries give a confluent path normal form",
        "path-rank row fibers and carry-width debt",
        "smoothing entropy and product incidence",
        "HYP-3000; summand-graph-fermat-zeckendorf.md",
        (5, 5, 4, 5, 4, 4, 5),
    ),
    Mode(
        "collatz_affine_halving",
        "+1,/2",
        "3n+1 then divide by 2^k; formal logarithm turns +1 into bounded defect",
        "valuation-drop recursion and iterated-linearization scale",
        "exact integer residue dependence under scalar drift models",
        "collatz-iterated-log-tower.md; lrc-iterated-logs-are-the-inverse-hyperoperation-tower-s597.md",
        (4, 4, 3, 5, 4, 3, 5),
    ),
    Mode(
        "triune_fraction_recursion",
        "sum/product/fraction/recursion",
        "sum -> product -> fraction -> recursion -> next sum",
        "boundary memory: owner, branch, lift, deletion, carry, route",
        "raw scalar and static sum-product shadows",
        "triune-cycle-everywhere-s664.md; pi_sum_product_fraction_s662.py",
        (5, 5, 4, 5, 5, 4, 5),
    ),
    Mode(
        "smoothing_switchboard",
        "global/local smoothing",
        "packet route first, smoothing kernel second",
        "which Fejer/Ramanujan/Kaczynski/large-sieve kernel may forget which label",
        "raw scalar smoothing route without endpoint/packet side channels",
        "lrc14-smoothing-switchboard-codex-s164.md; LTI-151",
        (5, 5, 5, 4, 5, 4, 4),
    ),
]

TIE_PATH = [mode.name for mode in MODES]
INDEX = {mode.name: i for i, mode in enumerate(MODES)}


def compare(a: Mode, b: Mode) -> str:
    wins_a = 0
    wins_b = 0
    for av, bv in zip(a.features, b.features):
        if av > bv:
            wins_a += 1
        elif bv > av:
            wins_b += 1
    if wins_a > wins_b:
        return a.name
    if wins_b > wins_a:
        return b.name
    return a.name if INDEX[a.name] < INDEX[b.name] else b.name


def build_edges() -> dict[tuple[str, str], str]:
    edges = {}
    for a, b in combinations(MODES, 2):
        winner = compare(a, b)
        loser = b.name if winner == a.name else a.name
        edges[(winner, loser)] = winner
    return edges


EDGES = build_edges()


def beats(a: str, b: str) -> bool:
    if a == b:
        return False
    return (a, b) in EDGES


def score(mode: Mode) -> int:
    return sum(1 for other in MODES if beats(mode.name, other.name))


def directed_3cycles() -> list[tuple[str, str, str]]:
    cycles = []
    for a, b, c in combinations([m.name for m in MODES], 3):
        if beats(a, b) and beats(b, c) and beats(c, a):
            cycles.append((a, b, c))
        elif beats(a, c) and beats(c, b) and beats(b, a):
            cycles.append((a, c, b))
    return cycles


def sccs() -> list[list[str]]:
    vertices = [m.name for m in MODES]
    seen = set()
    components = []

    def reachable(start: str) -> set[str]:
        stack = [start]
        out = set()
        while stack:
            v = stack.pop()
            if v in out:
                continue
            out.add(v)
            stack.extend(w for w in vertices if beats(v, w) and w not in out)
        return out

    for v in vertices:
        if v in seen:
            continue
        rv = reachable(v)
        comp = sorted([w for w in vertices if v in reachable(w) and w in rv], key=INDEX.get)
        seen.update(comp)
        components.append(comp)
    return components


@lru_cache(None)
def hp_count(last: str | None, remaining: tuple[str, ...]) -> int:
    if not remaining:
        return 1
    total = 0
    for i, v in enumerate(remaining):
        if last is None or beats(last, v):
            total += hp_count(v, remaining[:i] + remaining[i + 1 :])
    return total


def print_mode_table() -> None:
    print("DICHOTOMY / RECURSION MODE ATLAS S170")
    print("=" * 72)
    print("Status: synthesis scout; not a proof.")
    print("Vertex set: proof carriers / recursion modes, not runners.")
    print("Pairwise observable:", ", ".join(FEATURES))
    print("Gauge: coordinate majority over the observable vector.")
    print("Tie Hamiltonian path:", " > ".join(TIE_PATH))
    print()
    for mode in MODES:
        print(f"[{mode.name}]")
        print(f"  dichotomy: {mode.dichotomy}")
        print(f"  recursion: {mode.recursion_law}")
        print(f"  preserves: {mode.preserves}")
        print(f"  destroys if scalarized: {mode.destroys_if_scalarized}")
        print(f"  anchors: {mode.anchors}")
        print(f"  feature vector: {dict(zip(FEATURES, mode.features))}")
        print()


def print_synthesis() -> None:
    print("SYNTHESIS")
    print("=" * 72)
    claims = [
        ("odd/even", "not mere parity; it decides whether doubling is an automorphism or a projection defect"),
        ("positive/negative", "observer-blind but pair-visible; signs are cuts selecting difference or sum clocks"),
        ("addition/multiplication", "addition creates folds, denominators, and endpoint pincers; multiplication transports along unit orbits and product shells"),
        ("+1,/2", "affine residue noise becomes bounded after the right logarithmic coordinate; the valuation drop is the recursion counter"),
        ("*2,+2", "doubling is vertical 2-adic descent; plus-two is horizontal line motion weaving through ruler-function depths"),
        ("sum/product/fraction/recursion", "the stable grammar is projection repair: aggregate, factor, remember boundary, unroll to the next aggregate"),
    ]
    for left, right in claims:
        print(f"- {left}: {right}.")
    print()
    print("LRC14 readout:")
    print("- Keep exact M=p/q and e=14p-q before applying any Farey mutation.")
    print("- Use product collapse p only after e=1; otherwise it is a false scalar quotient.")
    print("- Treat signs as a side-channel cut, not as an observer witness.")
    print("- Treat parity as the first CRT block where doubling may leak labels.")
    print("- Route every packet by proof currency before picking a smoothing kernel.")
    print()


def print_tournament() -> None:
    print("TOURNAMENT ANALYSIS")
    print("=" * 72)
    scores = {mode.name: score(mode) for mode in MODES}
    hist: dict[int, int] = {}
    for value in scores.values():
        hist[value] = hist.get(value, 0) + 1
    ordered = sorted(MODES, key=lambda m: (-scores[m.name], INDEX[m.name]))
    cycles = directed_3cycles()
    components = sccs()
    print("score_hist:", dict(sorted(hist.items())))
    print("scores:")
    for mode in ordered:
        print(f"  {mode.name}: {scores[mode.name]}")
    print("directed_3cycles:", len(cycles), cycles[:12])
    print("scc_sizes:", [len(c) for c in components])
    print("sccs:", components)
    print("hamiltonian_path_count:", hp_count(None, tuple(m.name for m in MODES)))
    print("one high-retention path:", " > ".join(m.name for m in ordered))
    print()
    print("Most useful next finite check:")
    print("  Add six fields to the HYP-2963 packet bank:")
    print("  parity_block, sign_cut_status, additive_pair_sum_owner,")
    print("  multiplicative_unit_orbit, recursion_boundary_state, smoothing_route.")
    print("  Then test whether every hard non-AP/GW packet is classified by exactly")
    print("  one primary mode plus named side-channel debts.")


def main() -> None:
    print_mode_table()
    print_synthesis()
    print_tournament()


if __name__ == "__main__":
    main()
