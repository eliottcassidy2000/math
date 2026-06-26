#!/usr/bin/env python3
"""Automatic gap carriers for the LRC14 proof search.

This is a finite, exact scout prompted by the Moser-de Bruijn sequence,
fibbinary numbers, finite automata, Ostrowski-Hadamard gaps, and the recent
2-adic Littlewood / Fermat-Catalan discussion.  It does not prove LRC14.

The useful question is quotient safety: which automatic-language coordinates
survive 2-adic shifts, carry moves, and residue projections before we try to
use them as LRC packet labels?

Tournament Analysis uses proof carriers as vertices, not runners or arcs.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations


LIMIT_BITS = 14
LIMIT = (1 << LIMIT_BITS) - 1


def bits_lsb(n: int, width: int | None = None) -> list[int]:
    if n == 0:
        bits = [0]
    else:
        bits = []
        value = n
        while value:
            bits.append(value & 1)
            value >>= 1
    if width is not None and len(bits) < width:
        bits.extend([0] * (width - len(bits)))
    return bits


def digits_base4_msb(n: int) -> list[int]:
    if n == 0:
        return [0]
    digits = []
    value = n
    while value:
        digits.append(value % 4)
        value //= 4
    return list(reversed(digits))


def is_moser(n: int) -> bool:
    """Moser-de Bruijn set: binary 1s only at even positions from the right."""
    return all(bit == 0 for pos, bit in enumerate(bits_lsb(n)) if pos % 2 == 1)


def is_moser_base4(n: int) -> bool:
    return all(digit in (0, 1) for digit in digits_base4_msb(n))


def is_fibbinary(n: int) -> bool:
    """Fibbinary set: no adjacent 1s in the binary expansion."""
    return (n & (n >> 1)) == 0


@dataclass(frozen=True)
class DFA:
    name: str
    alphabet: tuple[int, ...]
    start: str
    accept: frozenset[str]
    transitions: dict[tuple[str, int], str]
    read_direction: str

    def run(self, word: list[int]) -> tuple[bool, str]:
        state = self.start
        for digit in word:
            state = self.transitions[(state, digit)]
        return state in self.accept, state

    def transition_rows(self) -> list[str]:
        states = sorted({self.start, *self.accept, *(s for s, _ in self.transitions)})
        rows = []
        for state in states:
            cells = [f"{digit}->{self.transitions[(state, digit)]}" for digit in self.alphabet]
            mark = "*" if state in self.accept else " "
            rows.append(f"  {mark}{state}: " + ", ".join(cells))
        return rows


MOSER_BINARY_DFA = DFA(
    name="moser_binary_even_lsb",
    alphabet=(0, 1),
    start="even",
    accept=frozenset({"even", "odd"}),
    read_direction="LSB-first",
    transitions={
        ("even", 0): "odd",
        ("even", 1): "odd",
        ("odd", 0): "even",
        ("odd", 1): "dead",
        ("dead", 0): "dead",
        ("dead", 1): "dead",
    },
)

FIBBINARY_DFA = DFA(
    name="fibbinary_no_adjacent_lsb",
    alphabet=(0, 1),
    start="last0",
    accept=frozenset({"last0", "last1"}),
    read_direction="LSB-first",
    transitions={
        ("last0", 0): "last0",
        ("last0", 1): "last1",
        ("last1", 0): "last0",
        ("last1", 1): "dead",
        ("dead", 0): "dead",
        ("dead", 1): "dead",
    },
)

MOSER_BASE4_DFA = DFA(
    name="moser_base4_digits_0_1",
    alphabet=(0, 1, 2, 3),
    start="ok",
    accept=frozenset({"ok"}),
    read_direction="MSB-first",
    transitions={
        ("ok", 0): "ok",
        ("ok", 1): "ok",
        ("ok", 2): "dead",
        ("ok", 3): "dead",
        ("dead", 0): "dead",
        ("dead", 1): "dead",
        ("dead", 2): "dead",
        ("dead", 3): "dead",
    },
)


def fib_numbers(limit: int) -> list[int]:
    values = [0, 1]
    while len(values) <= limit:
        values.append(values[-1] + values[-2])
    return values


FIB = fib_numbers(40)


def count_by_length() -> list[dict[str, int]]:
    rows = []
    for width in range(1, LIMIT_BITS + 1):
        high = 1 << width
        moser = sum(1 for n in range(high) if is_moser(n))
        fibbinary = sum(1 for n in range(high) if is_fibbinary(n))
        both = sum(1 for n in range(high) if is_moser(n) and is_fibbinary(n))
        rows.append(
            {
                "bits": width,
                "moser": moser,
                "moser_formula": 1 << ((width + 1) // 2),
                "fibbinary": fibbinary,
                "fibbinary_formula": FIB[width + 2],
                "intersection": both,
            }
        )
    return rows


def first_values(predicate, count: int) -> list[int]:
    values = []
    n = 0
    while len(values) < count:
        if predicate(n):
            values.append(n)
        n += 1
    return values


def closure_stats(name: str, predicate, operation_name: str, operation) -> dict[str, object]:
    domain = [n for n in range(LIMIT + 1) if predicate(n) and operation(n) <= LIMIT]
    preserved = [n for n in domain if predicate(operation(n))]
    failures = [n for n in domain if n not in preserved]
    return {
        "language": name,
        "operation": operation_name,
        "domain": len(domain),
        "preserved": len(preserved),
        "failed": len(failures),
        "first_failures": failures[:8],
    }


def residue_profile(predicate, modulus: int = 14) -> Counter[int]:
    return Counter(n % modulus for n in range(LIMIT + 1) if predicate(n))


def prefix_eventual_period_stress(predicate, max_period: int = 64, min_cutoff: int = 512) -> str:
    """Finite stress only: looks for a small eventual period in the prefix."""
    prefix = [1 if predicate(n) else 0 for n in range(LIMIT + 1)]
    for cutoff in range(min_cutoff + 1):
        for period in range(1, max_period + 1):
            ok = True
            for i in range(cutoff + period, len(prefix)):
                if prefix[i] != prefix[i - period]:
                    ok = False
                    break
            if ok:
                return f"period={period}, cutoff={cutoff}"
    return f"none with period<={max_period}, cutoff<={min_cutoff}, prefix<=2^{LIMIT_BITS}-1"


@dataclass(frozen=True)
class Carrier:
    name: str
    role: str
    preserves: tuple[str, ...]
    destroys: tuple[str, ...]
    vector: tuple[int, ...]


FEATURES = (
    "packet_label_retention",
    "2adic_shift_control",
    "carry_boundary_memory",
    "gap_lacunarity",
    "sml_eventual_period_guard",
    "finite_automaton_exactness",
    "lrc_transfer_safety",
)

CARRIERS = (
    Carrier(
        "product_phase_automaton",
        "product DFA keeping both Moser parity phase and fibbinary previous-bit state",
        ("2-adic phase", "carry boundary", "regular-language exactness"),
        ("endpoint geometry and exact safe interval lengths",),
        (5, 5, 5, 4, 3, 5, 5),
    ),
    Carrier(
        "fibbinary_2adic_shift_dfa",
        "no-adjacent-1 language, closed under n -> 2n",
        ("2-adic shift closure", "Zeckendorf-style carry boundary"),
        ("base-4 lacunary atom phase", "endpoint owner geometry"),
        (4, 5, 5, 2, 2, 5, 4),
    ),
    Carrier(
        "fermat_catalan_valuation_gate",
        "p-adic valuation / local-to-global gate before importing recurrence zeros",
        ("valuation depth", "local prime labels", "named arithmetic obstruction"),
        ("automatic language fiber unless attached explicitly",),
        (4, 4, 4, 3, 5, 3, 4),
    ),
    Carrier(
        "sml_eventual_periodic_zero_gate",
        "linear-recurrence zero-set discipline: finite union of arithmetic progressions",
        ("eventual periodic index constraint", "anti-lacunary warning"),
        ("raw automatic value-set structure", "endpoint geometry"),
        (3, 3, 2, 1, 5, 2, 3),
    ),
    Carrier(
        "moser_base4_gap_dfa",
        "base-4 digit language 0/1, stable under n -> 4n but phase-breaking under n -> 2n",
        ("Hadamard-lacunary atom memory", "base-4 support", "Moser subset of fibbinary"),
        ("odd/even 2-adic phase if read as scalar values",),
        (4, 2, 3, 5, 1, 5, 3),
    ),
    Carrier(
        "ostrowski_hadamard_atom_gap",
        "primitive atom gap p_{j+1}/p_j=4 natural-boundary heuristic",
        ("lacunary atom warning",),
        ("finite carry language", "LRC packet labels", "endpoint owners"),
        (2, 1, 1, 5, 1, 2, 2),
    ),
    Carrier(
        "stick_potato_safe_geometry",
        "largest safe component / inner approximation guardrail from polygon geometry",
        ("safe-set geometry", "largest-component objective"),
        ("arithmetic valuation", "automatic carry language"),
        (4, 1, 1, 1, 1, 2, 3),
    ),
    Carrier(
        "raw_sequence_scalar_shadow",
        "membership count or first-values list without retained automaton state",
        ("quick weather report",),
        ("2-adic phase", "carry boundary", "valuation labels", "endpoint geometry"),
        (1, 1, 1, 1, 1, 1, 1),
    ),
)

TIE_PATH = [carrier.name for carrier in CARRIERS]
INDEX = {name: i for i, name in enumerate(TIE_PATH)}


def compare(a: Carrier, b: Carrier) -> str:
    wins_a = 0
    wins_b = 0
    for av, bv in zip(a.vector, b.vector):
        if av > bv:
            wins_a += 1
        elif bv > av:
            wins_b += 1
    if wins_a > wins_b:
        return a.name
    if wins_b > wins_a:
        return b.name
    return a.name if INDEX[a.name] < INDEX[b.name] else b.name


EDGES: set[tuple[str, str]] = set()
for left, right in combinations(CARRIERS, 2):
    win = compare(left, right)
    lose = right.name if win == left.name else left.name
    EDGES.add((win, lose))


def beats(a: str, b: str) -> bool:
    return (a, b) in EDGES


def score(name: str) -> int:
    return sum(1 for other in TIE_PATH if beats(name, other))


def directed_3cycles() -> list[tuple[str, str, str]]:
    cycles = []
    for a, b, c in combinations(TIE_PATH, 3):
        if beats(a, b) and beats(b, c) and beats(c, a):
            cycles.append((a, b, c))
        elif beats(a, c) and beats(c, b) and beats(b, a):
            cycles.append((a, c, b))
    return cycles


def sccs() -> list[list[str]]:
    def reachable(start: str) -> set[str]:
        stack = [start]
        out = set()
        while stack:
            vertex = stack.pop()
            if vertex in out:
                continue
            out.add(vertex)
            stack.extend(other for other in TIE_PATH if beats(vertex, other))
        return out

    seen: set[str] = set()
    components = []
    for vertex in TIE_PATH:
        if vertex in seen:
            continue
        rv = reachable(vertex)
        comp = sorted([other for other in TIE_PATH if vertex in reachable(other) and other in rv], key=INDEX.get)
        seen.update(comp)
        components.append(comp)
    return components


@lru_cache(None)
def hp_count(mask: int, last: int) -> int:
    if mask == (1 << len(TIE_PATH)) - 1:
        return 1
    total = 0
    last_name = TIE_PATH[last]
    for nxt, next_name in enumerate(TIE_PATH):
        if mask & (1 << nxt):
            continue
        if beats(last_name, next_name):
            total += hp_count(mask | (1 << nxt), nxt)
    return total


def hamiltonian_path_count() -> int:
    return sum(hp_count(1 << i, i) for i in range(len(TIE_PATH)))


def print_section(title: str) -> None:
    print(f"\n## {title}")


def main() -> None:
    print("LRC14 Moser/Fibbinary Automatic Gap Carrier")
    print(f"finite audit range: 0..{LIMIT} (2^{LIMIT_BITS}-1)")

    print_section("Finite Automata")
    for dfa in (MOSER_BINARY_DFA, FIBBINARY_DFA, MOSER_BASE4_DFA):
        print(f"{dfa.name} ({dfa.read_direction})")
        for row in dfa.transition_rows():
            print(row)

    print_section("First Values")
    print("Moser-de Bruijn:", first_values(is_moser, 24))
    print("Fibbinary:", first_values(is_fibbinary, 32))
    print("Moser subset of fibbinary through limit:", all(not is_moser(n) or is_fibbinary(n) for n in range(LIMIT + 1)))
    print("Binary DFA agrees with base-4 DFA through limit:", all(is_moser(n) == is_moser_base4(n) for n in range(LIMIT + 1)))

    print_section("Length Counts")
    print("bits  moser/formula  fibbinary/formula  intersection")
    for row in count_by_length():
        print(
            f"{row['bits']:>4}  "
            f"{row['moser']:>5}/{row['moser_formula']:<5}  "
            f"{row['fibbinary']:>8}/{row['fibbinary_formula']:<8}  "
            f"{row['intersection']:>5}"
        )

    print_section("2-adic And Carry Closure")
    for predicate_name, predicate in (("moser", is_moser), ("fibbinary", is_fibbinary)):
        for operation_name, operation in (
            ("times2", lambda n: 2 * n),
            ("times4", lambda n: 4 * n),
            ("plus1", lambda n: n + 1),
        ):
            stats = closure_stats(predicate_name, predicate, operation_name, operation)
            print(
                f"{stats['language']:>9} {stats['operation']:>7}: "
                f"{stats['preserved']:>4}/{stats['domain']:<4} preserved, "
                f"failed={stats['failed']:<4} first_failures={stats['first_failures']}"
            )

    print_section("Residues Mod 14")
    for predicate_name, predicate in (("moser", is_moser), ("fibbinary", is_fibbinary)):
        profile = residue_profile(predicate)
        ordered = [profile[r] for r in range(14)]
        print(f"{predicate_name:>9}: counts={ordered} min={min(ordered)} max={max(ordered)}")

    print_section("Eventual-Period Stress")
    print("This is only a finite stress test against treating raw value sets as SML zero sets.")
    print("moser:", prefix_eventual_period_stress(is_moser))
    print("fibbinary:", prefix_eventual_period_stress(is_fibbinary))

    print_section("Tournament Analysis")
    print("vertex-set challenge:")
    print("  considered runners, gaps, fixed sections, residues, cover arcs, Fourier modes,")
    print("  matroid circuits, and proof obligations; chose proof-language carriers because")
    print("  the LRC predicate being tested is quotient safety under lost automaton state.")
    print("features:", ", ".join(FEATURES))
    print("tie Hamiltonian path:", " > ".join(TIE_PATH))
    scores = Counter(score(name) for name in TIE_PATH)
    print("score_hist:", sorted(scores.items()))
    print("directed_3cycles:", len(directed_3cycles()))
    print("scc_sizes:", [len(comp) for comp in sccs()])
    print("hamiltonian_path_count:", hamiltonian_path_count())
    ordered = sorted(TIE_PATH, key=lambda name: (-score(name), INDEX[name]))
    print("score_order:", " > ".join(ordered))

    print_section("Proof-Order Takeaway")
    print("1. Moser-de Bruijn is an exact automatic subset of fibbinary, but its natural")
    print("   base-4 gap coordinate is stable under n->4n, not n->2n.  LRC14 must retain")
    print("   the even/odd 2-adic phase if this carrier is used.")
    print("2. Fibbinary is the safer 2-adic normal-form language: n->2n preserves it,")
    print("   while +1 exposes the carry boundary.  That matches the Zeckendorf/path")
    print("   normal-form lane from HYP-3000/HYP-3003.")
    print("3. Raw automatic value sets fail the finite SML-period stress here.  If a")
    print("   Fermat-Catalan or linear-recurrence zero-set argument is imported, its")
    print("   ultimately-periodic index coordinate must be explicit rather than inferred")
    print("   from Moser/fibbinary membership.")
    print("4. The ACM stick/potato paper is useful only as a largest-safe-component")
    print("   geometry analogy; it does not supply arithmetic or automaton closure.")


if __name__ == "__main__":
    main()
