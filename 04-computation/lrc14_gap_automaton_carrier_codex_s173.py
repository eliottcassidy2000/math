#!/usr/bin/env python3
"""Gap/automaton carrier scout for LRC14.

This is a synthesis computation, not a proof.  It tests whether the newest
gap-language prompts should become packet fields for the LRC14 proof stack.

Tournament Analysis is over proof carriers and automaton/gap languages, not
runners.  The pairwise observable is the retained proof payload; several
priority gauges are used, and induced n-subtournaments are canonicalized up to
isomorphism to expose which tournament classes this carrier surface can
generate.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
from typing import Iterable


LIMIT = 4096
MAX_P = 384


def digits(n: int, base: int) -> tuple[int, ...]:
    if n == 0:
        return (0,)
    out: list[int] = []
    while n:
        out.append(n % base)
        n //= base
    return tuple(reversed(out))


def digit_string(n: int, base: int) -> str:
    alphabet = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    return "".join(alphabet[d] for d in digits(n, base))


def is_fibbinary(n: int) -> bool:
    bits = digit_string(n, 2)
    return "11" not in bits


def is_moser_de_bruijn(n: int) -> bool:
    return all(d in (0, 1) for d in digits(n, 4))


def seq_fibbinary(limit: int = LIMIT) -> list[int]:
    return [n for n in range(limit + 1) if is_fibbinary(n)]


def seq_moser(limit: int = LIMIT) -> list[int]:
    return [n for n in range(limit + 1) if is_moser_de_bruijn(n)]


def residue_hist(values: Iterable[int], modulus: int = 14) -> dict[int, int]:
    hist = Counter(v % modulus for v in values)
    return {r: hist.get(r, 0) for r in range(modulus)}


def mixed_residue_count(predicate, limit: int = LIMIT, modulus: int = 14) -> int:
    mixed = 0
    for r in range(modulus):
        vals = [predicate(n) for n in range(r, limit + 1, modulus)]
        if any(vals) and not all(vals):
            mixed += 1
    return mixed


def closure_violations(values: set[int], multiplier: int, limit: int = LIMIT) -> int:
    return sum(1 for n in values if multiplier * n <= limit and multiplier * n not in values)


@dataclass(frozen=True)
class Automaton:
    name: str
    alphabet: tuple[str, ...]
    states: tuple[str, ...]
    start: str
    accept: frozenset[str]
    transition: dict[tuple[str, str], str]
    interpretation: str

    def accepts(self, word: str) -> bool:
        state = self.start
        for ch in word:
            state = self.transition[(state, ch)]
        return state in self.accept


def fibbinary_automaton() -> Automaton:
    transition = {
        ("last0", "0"): "last0",
        ("last0", "1"): "last1",
        ("last1", "0"): "last0",
        ("last1", "1"): "dead",
        ("dead", "0"): "dead",
        ("dead", "1"): "dead",
    }
    return Automaton(
        name="fibbinary_no_adjacent_one",
        alphabet=("0", "1"),
        states=("last0", "last1", "dead"),
        start="last0",
        accept=frozenset(("last0", "last1")),
        transition=transition,
        interpretation="binary normal forms with no adjacent 1s; doubling appends 0",
    )


def moser_automaton() -> Automaton:
    transition: dict[tuple[str, str], str] = {}
    for ch in "0123":
        transition[("ok", ch)] = "ok" if ch in "01" else "dead"
        transition[("dead", ch)] = "dead"
    return Automaton(
        name="moser_base4_01_digit_language",
        alphabet=("0", "1", "2", "3"),
        states=("ok", "dead"),
        start="ok",
        accept=frozenset(("ok",)),
        transition=transition,
        interpretation="sums of distinct powers of 4; multiply by 4 appends 0",
    )


@dataclass(frozen=True)
class Carrier:
    name: str
    finite_state: int
    gap_support: int
    doubling_transfer: int
    valuation_budget: int
    lrc_packet_transfer: int
    boundary_guardrail: int
    quotient_safety: int
    proof_maturity: int
    destroyed: str

    def vector(self) -> dict[str, int]:
        return {
            "finite_state": self.finite_state,
            "gap_support": self.gap_support,
            "doubling_transfer": self.doubling_transfer,
            "valuation_budget": self.valuation_budget,
            "lrc_packet_transfer": self.lrc_packet_transfer,
            "boundary_guardrail": self.boundary_guardrail,
            "quotient_safety": self.quotient_safety,
            "proof_maturity": self.proof_maturity,
        }


GAUGES = {
    "proof_priority": (
        "lrc_packet_transfer",
        "quotient_safety",
        "boundary_guardrail",
        "proof_maturity",
        "finite_state",
        "gap_support",
        "doubling_transfer",
        "valuation_budget",
    ),
    "automaton_priority": (
        "finite_state",
        "doubling_transfer",
        "quotient_safety",
        "lrc_packet_transfer",
        "gap_support",
        "boundary_guardrail",
        "valuation_budget",
        "proof_maturity",
    ),
    "gap_boundary_priority": (
        "gap_support",
        "boundary_guardrail",
        "lrc_packet_transfer",
        "quotient_safety",
        "finite_state",
        "doubling_transfer",
        "valuation_budget",
        "proof_maturity",
    ),
    "dyadic_valuation_priority": (
        "doubling_transfer",
        "valuation_budget",
        "finite_state",
        "lrc_packet_transfer",
        "quotient_safety",
        "gap_support",
        "boundary_guardrail",
        "proof_maturity",
    ),
    "exception_budget_priority": (
        "valuation_budget",
        "boundary_guardrail",
        "quotient_safety",
        "lrc_packet_transfer",
        "finite_state",
        "gap_support",
        "doubling_transfer",
        "proof_maturity",
    ),
}


MAJORITY_FIELDS = (
    "finite_state",
    "gap_support",
    "doubling_transfer",
    "valuation_budget",
    "lrc_packet_transfer",
    "boundary_guardrail",
    "quotient_safety",
    "proof_maturity",
)


def carrier_bank() -> list[Carrier]:
    return [
        Carrier(
            "labelled_lrc_gap_automaton_packet",
            finite_state=3,
            gap_support=3,
            doubling_transfer=3,
            valuation_budget=3,
            lrc_packet_transfer=4,
            boundary_guardrail=3,
            quotient_safety=4,
            proof_maturity=2,
            destroyed="nothing if automaton state, gap support, valuation budget, and packet route are retained",
        ),
        Carrier(
            "two_adic_littlewood_hurwitz_doubling",
            finite_state=2,
            gap_support=2,
            doubling_transfer=4,
            valuation_budget=4,
            lrc_packet_transfer=2,
            boundary_guardrail=2,
            quotient_safety=3,
            proof_maturity=3,
            destroyed="endpoint ownership and exact LRC route unless recoupled to packet labels",
        ),
        Carrier(
            "fibbinary_no_adjacent_language",
            finite_state=4,
            gap_support=2,
            doubling_transfer=4,
            valuation_budget=1,
            lrc_packet_transfer=2,
            boundary_guardrail=1,
            quotient_safety=3,
            proof_maturity=2,
            destroyed="magnitude and endpoint geometry if only membership is kept",
        ),
        Carrier(
            "moser_base4_digit_language",
            finite_state=4,
            gap_support=3,
            doubling_transfer=2,
            valuation_budget=1,
            lrc_packet_transfer=2,
            boundary_guardrail=2,
            quotient_safety=3,
            proof_maturity=2,
            destroyed="2-adic transition data; multiply-by-2 is not the native closure",
        ),
        Carrier(
            "ostrowski_hadamard_lacunary_boundary",
            finite_state=1,
            gap_support=4,
            doubling_transfer=2,
            valuation_budget=1,
            lrc_packet_transfer=1,
            boundary_guardrail=4,
            quotient_safety=3,
            proof_maturity=3,
            destroyed="finite packet identity unless the lacunary support is attached to a source family",
        ),
        Carrier(
            "fermat_catalan_exponent_budget",
            finite_state=1,
            gap_support=1,
            doubling_transfer=1,
            valuation_budget=4,
            lrc_packet_transfer=1,
            boundary_guardrail=2,
            quotient_safety=3,
            proof_maturity=2,
            destroyed="which additive or multiplicative fiber generated the exponent budget",
        ),
        Carrier(
            "skolem_mahler_lech_zero_set_guard",
            finite_state=3,
            gap_support=2,
            doubling_transfer=2,
            valuation_budget=2,
            lrc_packet_transfer=1,
            boundary_guardrail=3,
            quotient_safety=3,
            proof_maturity=2,
            destroyed="local endpoint owners and Haar/Fejer certificate data",
        ),
        Carrier(
            "large_stick_potato_visibility_core",
            finite_state=1,
            gap_support=1,
            doubling_transfer=0,
            valuation_budget=0,
            lrc_packet_transfer=2,
            boundary_guardrail=3,
            quotient_safety=2,
            proof_maturity=1,
            destroyed="number-theoretic clocks; keeps only largest visible/inscribed core shape",
        ),
        Carrier(
            "raw_sequence_membership_scalar",
            finite_state=0,
            gap_support=0,
            doubling_transfer=0,
            valuation_budget=0,
            lrc_packet_transfer=0,
            boundary_guardrail=0,
            quotient_safety=0,
            proof_maturity=0,
            destroyed="all packet, endpoint, scale, automaton-state, and proof-route labels",
        ),
    ]


def build_tournament(carriers: list[Carrier], priority: tuple[str, ...]) -> list[list[bool]]:
    n = len(carriers)
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        vi = carriers[i].vector()
        vj = carriers[j].vector()
        winner: int | None = None
        for key in priority:
            if vi[key] > vj[key]:
                winner = i
                break
            if vi[key] < vj[key]:
                winner = j
                break
        if winner is None:
            winner = i
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj


def build_majority_tournament(carriers: list[Carrier]) -> list[list[bool]]:
    """Orient by fieldwise majority over retained proof dimensions."""
    n = len(carriers)
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        vi = carriers[i].vector()
        vj = carriers[j].vector()
        wi = sum(vi[key] > vj[key] for key in MAJORITY_FIELDS)
        wj = sum(vj[key] > vi[key] for key in MAJORITY_FIELDS)
        winner = i if wi >= wj else j
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj


def tournament_for_gauge(carriers: list[Carrier], gauge_name: str) -> list[list[bool]]:
    if gauge_name == "majority_payload":
        return build_majority_tournament(carriers)
    return build_tournament(carriers, GAUGES[gauge_name])


def gauge_names() -> tuple[str, ...]:
    return tuple(GAUGES) + ("majority_payload",)


def directed_3cycles(adj: list[list[bool]]) -> int:
    out = 0
    for i, j, k in combinations(range(len(adj)), 3):
        degs = []
        for v in (i, j, k):
            degs.append(sum(1 for w in (i, j, k) if adj[v][w]))
        if sorted(degs) == [1, 1, 1]:
            out += 1
    return out


def scc_components(adj: list[list[bool]]) -> list[list[int]]:
    n = len(adj)
    index = 0
    stack: list[int] = []
    on_stack = [False] * n
    indices = [-1] * n
    lowlink = [0] * n
    comps: list[list[int]] = []

    def strongconnect(v: int) -> None:
        nonlocal index
        indices[v] = index
        lowlink[v] = index
        index += 1
        stack.append(v)
        on_stack[v] = True
        for w in range(n):
            if not adj[v][w]:
                continue
            if indices[w] == -1:
                strongconnect(w)
                lowlink[v] = min(lowlink[v], lowlink[w])
            elif on_stack[w]:
                lowlink[v] = min(lowlink[v], indices[w])
        if lowlink[v] == indices[v]:
            comp: list[int] = []
            while True:
                w = stack.pop()
                on_stack[w] = False
                comp.append(w)
                if w == v:
                    break
            comps.append(sorted(comp))

    for v in range(n):
        if indices[v] == -1:
            strongconnect(v)
    return sorted(comps, key=lambda comp: (-len(comp), comp))


def scc_sizes(adj: list[list[bool]]) -> list[int]:
    return [len(comp) for comp in scc_components(adj)]


def hamiltonian_path_count(adj: list[list[bool]]) -> int:
    n = len(adj)

    @lru_cache(maxsize=None)
    def count(mask: int, last: int) -> int:
        if mask == (1 << n) - 1:
            return 1
        total = 0
        for nxt in range(n):
            if mask & (1 << nxt):
                continue
            if adj[last][nxt]:
                total += count(mask | (1 << nxt), nxt)
        return total

    return sum(count(1 << start, start) for start in range(n))


def word_for_order(adj: list[list[bool]], order: tuple[int, ...]) -> str:
    bits: list[str] = []
    for a in range(len(order)):
        for b in range(a + 1, len(order)):
            bits.append("1" if adj[order[a]][order[b]] else "0")
    return "".join(bits)


def canonical_word(adj: list[list[bool]], subset: tuple[int, ...]) -> str:
    return min(word_for_order(adj, order) for order in permutations(subset))


def tournament_fingerprint(adj: list[list[bool]], carriers: list[Carrier]) -> dict[str, object]:
    outdeg = [sum(row) for row in adj]
    path_order = sorted(range(len(carriers)), key=lambda i: (-outdeg[i], carriers[i].name))
    components = scc_components(adj)
    return {
        "score_hist": dict(sorted(Counter(outdeg).items())),
        "directed_3cycles": directed_3cycles(adj),
        "scc_sizes": [len(comp) for comp in components],
        "nontrivial_sccs": [
            tuple(carriers[i].name for i in comp) for comp in components if len(comp) > 1
        ],
        "hamiltonian_path_count": hamiltonian_path_count(adj),
        "outdegree_path": " > ".join(carriers[i].name for i in path_order),
    }


def induced_class_census(carriers: list[Carrier], sizes: tuple[int, ...] = (4, 5, 6)) -> dict[int, dict[str, object]]:
    out: dict[int, dict[str, object]] = {}
    for size in sizes:
        class_to_witness: dict[str, tuple[str, tuple[str, ...]]] = {}
        class_to_gauges: defaultdict[str, set[str]] = defaultdict(set)
        for gauge_name in gauge_names():
            adj = tournament_for_gauge(carriers, gauge_name)
            for subset in combinations(range(len(carriers)), size):
                word = canonical_word(adj, subset)
                class_to_gauges[word].add(gauge_name)
                if word not in class_to_witness:
                    class_to_witness[word] = (
                        gauge_name,
                        tuple(carriers[i].name for i in subset),
                    )
        out[size] = {
            "distinct_classes": len(class_to_witness),
            "witnesses": class_to_witness,
            "gauges": class_to_gauges,
        }
    return out


def reciprocal_sum(triple: tuple[int, int, int]) -> Fraction:
    p, q, r = triple
    return Fraction(1, p) + Fraction(1, q) + Fraction(1, r)


def fermat_catalan_budget(max_exp: int = 10) -> dict[str, object]:
    triples = [(p, q, r) for p, q, r in product(range(2, max_exp + 1), repeat=3) if p <= q <= r]
    hyperbolic = [t for t in triples if reciprocal_sum(t) < 1]
    boundary = [t for t in triples if reciprocal_sum(t) == 1]
    elliptic = [t for t in triples if reciprocal_sum(t) > 1]
    first_hyper = sorted(hyperbolic, key=lambda t: (sum(t), t))[:12]
    return {
        "search_box": f"2 <= p <= q <= r <= {max_exp}",
        "hyperbolic_count": len(hyperbolic),
        "boundary_count": len(boundary),
        "elliptic_count": len(elliptic),
        "first_hyperbolic_triples": first_hyper,
    }


def unit_excess_hits(fib: set[int], moser: set[int]) -> list[tuple[int, int, str]]:
    hits: list[tuple[int, int, str]] = []
    for p in range(1, MAX_P + 1):
        q = 14 * p - 1
        tags: list[str] = []
        if q in fib:
            tags.append("fibbinary")
        if q in moser:
            tags.append("moser")
        if tags:
            hits.append((p, q, "+".join(tags)))
    return hits


def format_hist(hist: dict[int, int]) -> str:
    return "{" + ", ".join(f"{k}:{v}" for k, v in hist.items()) + "}"


def main() -> None:
    fib = seq_fibbinary()
    moser = seq_moser()
    fib_set = set(fib)
    moser_set = set(moser)
    both = sorted(fib_set & moser_set)
    autos = [fibbinary_automaton(), moser_automaton()]
    carriers = carrier_bank()

    print("=== LRC14 gap/automaton carrier scout S173 ===")
    print("Status: synthesis computation and quotient guardrail; not a proof")
    print()

    print("External source notes")
    print("  arXiv:2506.04110: 2-adic Littlewood work uses doubling of continued-fraction data as a live obstruction carrier; for this pass it contributes the retained field 'doubling transition state'.")
    print("  Ostrowski-Hadamard gap theorem: lacunary exponent support with ratio bounded above 1 creates natural-boundary behavior; for LRC it is a boundary warning, not a scalar certificate.")
    print("  DOI 10.5555/1109557.1109610: the ACM record resolves to a computational-geometry 'large sticks and potatoes' paper, not Fermat-Catalan; here it contributes only a visibility/inscribed-core analogy.")
    print("  Fermat-Catalan: the useful abstraction is the hyperbolic exponent budget 1/p + 1/q + 1/r < 1 as a finite-exception ledger.")
    print()

    print("Finite automata")
    for auto in autos:
        print(f"  {auto.name}: states={list(auto.states)}, start={auto.start}, accept={sorted(auto.accept)}")
        print(f"    interpretation: {auto.interpretation}")
    print()

    print(f"Finite sequence audit through N={LIMIT}")
    print(f"  fibbinary_count={len(fib)} first20={fib[:20]}")
    print(f"  moser_count={len(moser)} first20={moser[:20]}")
    print(f"  intersection_count={len(both)} first20={both[:20]}")
    print(f"  fibbinary_residues_mod14={format_hist(residue_hist(fib))}")
    print(f"  moser_residues_mod14={format_hist(residue_hist(moser))}")
    print(f"  fibbinary_mixed_residue_classes_mod14={mixed_residue_count(is_fibbinary)}")
    print(f"  moser_mixed_residue_classes_mod14={mixed_residue_count(is_moser_de_bruijn)}")
    print(f"  fibbinary_double_closure_violations={closure_violations(fib_set, 2)}")
    print(f"  moser_times4_closure_violations={closure_violations(moser_set, 4)}")
    print(f"  moser_double_closure_violations={closure_violations(moser_set, 2)}")
    print("  lacunary_support_ratios: powers_of_2=2, powers_of_4=4; full fibbinary and Moser languages are finite-state sparse but not Hadamard-lacunary as full sequences.")
    print()

    print("LRC14 unit-excess chain q=14p-1 against automaton languages")
    hits = unit_excess_hits(fib_set, moser_set)
    print(f"  hits_with_p<={MAX_P}: {len(hits)}")
    for p, q, tags in hits[:18]:
        print(f"    p={p:3d} q={q:5d} tags={tags:<16s} bin={digit_string(q, 2):>14s} base4={digit_string(q, 4):>8s}")
    if len(hits) > 18:
        print(f"    ... {len(hits) - 18} more")
    print("  readout: membership alone mixes every residue class; the packet must retain automaton state plus exact M/qdiv/endpoint labels.")
    print()

    budget = fermat_catalan_budget()
    print("Fermat-Catalan style exponent budget")
    print(f"  box={budget['search_box']}")
    print(f"  counts: hyperbolic={budget['hyperbolic_count']} boundary={budget['boundary_count']} elliptic={budget['elliptic_count']}")
    print(f"  first_hyperbolic_triples={budget['first_hyperbolic_triples']}")
    print("  LRC use: if a residual packet needs simultaneous additive, multiplicative, and valuation exponents in the hyperbolic region, record it as finite-exception debt rather than a free scalar.")
    print()

    print("Proposed packet fields")
    fields = [
        "automaton_language_id",
        "automaton_state_word",
        "gap_support_ratio_label",
        "hadamard_boundary_warning",
        "doubling_transition_state",
        "base4_digit_mask",
        "zeckendorf_carry_state",
        "valuation_exponent_budget",
        "finite_exception_budget",
        "visibility_core_label",
    ]
    for field in fields:
        print(f"  - {field}")
    print()

    print("Tournament Analysis over proof/gap carriers")
    for gauge_name in gauge_names():
        adj = tournament_for_gauge(carriers, gauge_name)
        fp = tournament_fingerprint(adj, carriers)
        print(f"  gauge={gauge_name}")
        if gauge_name == "majority_payload":
            print(f"    rule=fieldwise majority over {MAJORITY_FIELDS}; ties follow declared carrier order")
        else:
            print(f"    priority={GAUGES[gauge_name]}")
        print(f"    score_hist={fp['score_hist']}")
        print(f"    directed_3cycles={fp['directed_3cycles']}")
        print(f"    scc_sizes={fp['scc_sizes']}")
        if fp["nontrivial_sccs"]:
            print(f"    nontrivial_sccs={fp['nontrivial_sccs']}")
        print(f"    hamiltonian_path_count={fp['hamiltonian_path_count']}")
        print(f"    outdegree_path={fp['outdegree_path']}")
    print()

    print("Induced tournament isomorphism-class census")
    census = induced_class_census(carriers)
    for size, data in census.items():
        print(f"  n={size}: distinct_classes_generated={data['distinct_classes']}")
        witnesses = data["witnesses"]
        for idx, (word, (gauge, names)) in enumerate(sorted(witnesses.items())[:6], start=1):
            print(f"    class{idx}: word={word} gauge={gauge} subset={names}")
        if data["distinct_classes"] > 6:
            print(f"    ... {data['distinct_classes'] - 6} more classes")
    print("  readout: this is an achievable-subset census for the chosen carrier surface, not a universal tournament enumeration.  It gives future agents exact class words to compare against when a new proof carrier is added.")
    print()

    print("LRC14 proof readout")
    print("  Necessary condition S173-A: any sequence/gap proof quotient must retain the automaton state or prove that the LRC predicate is constant on each automaton-language fiber.")
    print("  Necessary condition S173-B: doubling-compatible evidence must say whether the native transition is x -> 2x, x -> 4x, or a recentered 2-adic continued-fraction/Hurwitz move.")
    print("  Necessary condition S173-C: lacunary gap claims are boundary-obstruction claims unless attached to a finite source-packet family with exact endpoints.")
    print("  Necessary condition S173-D: Fermat-Catalan style exponent budgets are useful only as finite-exception ledgers over labelled packets; they cannot replace exact M/qdiv/route labels.")
    print("  Candidate theorem target: a primitive non-AP/GW LRC14 residual whose exact packet data are finite-state under all declared carry/doubling languages must either emit a known q/Fejer/Ramanujan/Haar/state-lift certificate or expose a named nonregular/natural-boundary residual sector.")


if __name__ == "__main__":
    main()
