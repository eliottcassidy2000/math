#!/usr/bin/env python3
"""LRC14 lacunary finite-automaton carrier scout.

This is a proof-interface note with exact finite checks, not a proof of LRC14.

Prompt sources used as carriers:
  * 2-adic Littlewood / Hurwitz multiplication-by-2 finite transducer:
    https://arxiv.org/abs/2506.04110
  * Ostrowski-Hadamard gap theorem:
    https://en.wikipedia.org/wiki/Ostrowski%E2%80%93Hadamard_gap_theorem
  * "Finding large sticks and potatoes in polygons" as the geometric analogue
    of computing the largest safe component:
    https://doi.org/10.1145/1109557.1109610
  * Moser-de Bruijn sequence: base-4 digits in {0,1}
  * Fibbinary numbers: binary words with no adjacent 1s

Tournament Analysis declaration:
  vertices: proof carriers / finite automata / geometric exact certificates,
    not runners or arcs;
  pairwise observable: retained LRC predicate, exact scale retention,
    finite-state checkability, carry-collision control, quotient-loss warning;
  switch/gauge: majority comparison of the observable vector;
  tie Hamiltonian path: the listed carrier order.

Assumption challenge:
  Alternate vertex sets considered: runners, gaps, residue classes, fixed circle
  sections, section boundaries, wall-crossing events, cover arcs, Fourier modes,
  matroid circuits, automaton states, and proof obligations.  This script uses
  proof obligations because the LRC predicate is "there is a safe time at exact
  threshold 1/14"; raw runner or sequence vertices destroy exact Farey scale and
  endpoint ownership.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations


N = 14
THRESHOLD = Fraction(1, N)


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def dist0(x: Fraction) -> Fraction:
    r = frac_part(x)
    return min(r, 1 - r)


def exact_score(row: tuple[int, ...], t: Fraction) -> Fraction:
    return min(dist0(Fraction(v) * t) for v in row)


def candidate_times(row: tuple[int, ...]) -> set[Fraction]:
    out: set[Fraction] = set()
    for i, a in enumerate(row):
        for m in range(a):
            t = Fraction(2 * m + 1, 2 * a)
            if 0 < t < 1:
                out.add(t)
        for b in row[i + 1 :]:
            for den in (a + b, abs(a - b)):
                if den <= 0:
                    continue
                for m in range(1, den):
                    out.add(Fraction(m, den))
    return out


def exact_maximin(row: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    best = Fraction(-1)
    best_t = Fraction(0)
    for t in candidate_times(tuple(sorted(row))):
        val = exact_score(row, t)
        if val > best or (val == best and t < best_t):
            best = val
            best_t = t
    return best, best_t


def danger_intervals(row: tuple[int, ...], delta: Fraction = THRESHOLD) -> list[tuple[Fraction, Fraction]]:
    intervals: list[tuple[Fraction, Fraction]] = []
    for v in row:
        half = delta / v
        for k in range(v):
            center = Fraction(k, v)
            lo = center - half
            hi = center + half
            if lo < 0:
                intervals.append((lo + 1, Fraction(1)))
                intervals.append((Fraction(0), hi))
            elif hi > 1:
                intervals.append((lo, Fraction(1)))
                intervals.append((Fraction(0), hi - 1))
            else:
                intervals.append((lo, hi))
    return intervals


def safe_components(row: tuple[int, ...], delta: Fraction = THRESHOLD) -> tuple[Fraction, Fraction, int]:
    intervals = sorted(danger_intervals(row, delta))
    merged: list[list[Fraction]] = []
    for lo, hi in intervals:
        if not merged or lo > merged[-1][1]:
            merged.append([lo, hi])
        elif hi > merged[-1][1]:
            merged[-1][1] = hi

    safe: list[tuple[Fraction, Fraction]] = []
    cursor = Fraction(0)
    for lo, hi in merged:
        if lo > cursor:
            safe.append((cursor, lo))
        if hi > cursor:
            cursor = hi
    if cursor < 1:
        safe.append((cursor, Fraction(1)))

    lengths = [hi - lo for lo, hi in safe]
    return sum(lengths, Fraction(0)), (max(lengths) if lengths else Fraction(0)), len(safe)


def is_fibbinary(x: int) -> bool:
    return x > 0 and (x & (x << 1)) == 0


def is_moser_de_bruijn(x: int) -> bool:
    if x <= 0:
        return False
    y = x
    while y:
        if y % 4 not in (0, 1):
            return False
        y //= 4
    return True


def generate(predicate, count: int) -> tuple[int, ...]:
    out = []
    x = 1
    while len(out) < count:
        if predicate(x):
            out.append(x)
        x += 1
    return tuple(out)


def first_residue_transversal(predicate, modulus: int = N) -> tuple[int, ...]:
    found: dict[int, int] = {}
    x = 1
    while len(found) < modulus - 1:
        if predicate(x):
            r = x % modulus
            if r and r not in found:
                found[r] = x
        x += 1
        if x > 2_000_000:
            raise RuntimeError("residue transversal search exceeded limit")
    return tuple(found[r] for r in range(1, modulus))


def fibbinary_count_by_width(max_width: int) -> list[tuple[int, int]]:
    rows = []
    for width in range(1, max_width + 1):
        lo = 1 << (width - 1)
        hi = 1 << width
        rows.append((width, sum(1 for x in range(lo, hi) if is_fibbinary(x))))
    return rows


def moser_count_by_width(max_width: int) -> list[tuple[int, int]]:
    rows = []
    for width in range(1, max_width + 1):
        lo = 1 << (width - 1)
        hi = 1 << width
        rows.append((width, sum(1 for x in range(lo, hi) if is_moser_de_bruijn(x))))
    return rows


def residue_profile(row: tuple[int, ...]) -> str:
    counts = Counter(v % N for v in row)
    missing = [r for r in range(1, N) if counts[r] == 0]
    doubled = [r for r in range(1, N) if counts[r] > 1]
    zero = counts[0]
    return f"missing={missing or '-'} doubled={doubled or '-'} zero_count={zero}"


def fmt_fraction(x: Fraction) -> str:
    return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)


@dataclass(frozen=True)
class RowRecord:
    name: str
    row: tuple[int, ...]
    M: Fraction
    t: Fraction
    safe_mass: Fraction
    longest_safe: Fraction
    components: int

    @property
    def margin(self) -> Fraction:
        return self.M - THRESHOLD


def row_record(name: str, row: tuple[int, ...]) -> RowRecord:
    row = tuple(sorted(row))
    M, t = exact_maximin(row)
    safe_mass, longest, components = safe_components(row)
    return RowRecord(name, row, M, t, safe_mass, longest, components)


@dataclass(frozen=True)
class Carrier:
    name: str
    exact_scale: int
    lrc_predicate: int
    finite_state: int
    carry_collision: int
    quotient_warning: int
    note: str

    @property
    def vector(self) -> tuple[int, int, int, int, int]:
        return (
            self.exact_scale,
            self.lrc_predicate,
            self.finite_state,
            self.carry_collision,
            self.quotient_warning,
        )


CARRIERS = [
    Carrier(
        "farey_exact_M_scheduler",
        5,
        5,
        2,
        4,
        5,
        "keeps M=p/q and e=14p-q before any quotient",
    ),
    Carrier(
        "large_stick_safe_component",
        5,
        5,
        2,
        2,
        4,
        "computes exact safe intervals / largest safe stick",
    ),
    Carrier(
        "two_adic_hurwitz_transducer",
        4,
        4,
        5,
        5,
        4,
        "finite-state dyadic carry lane suggested by 2LC",
    ),
    Carrier(
        "fibbinary_no_adjacent_carry",
        3,
        3,
        5,
        4,
        3,
        "path-width-one carry language; Zeckendorf-compatible",
    ),
    Carrier(
        "moser_de_bruijn_base4_sparse",
        2,
        2,
        5,
        5,
        5,
        "base-4 no-carry language; lacunary generator warning",
    ),
    Carrier(
        "fermat_catalan_power_collision",
        2,
        2,
        1,
        4,
        4,
        "scarcity of primitive power collisions; stress lane only",
    ),
    Carrier(
        "ostrowski_hadamard_gap_guard",
        1,
        1,
        1,
        3,
        5,
        "Hadamard gaps warn that lacunary analytic continuation is brittle",
    ),
    Carrier(
        "raw_sequence_scalar",
        1,
        1,
        2,
        1,
        1,
        "sequence membership alone forgets scale, residues, owners",
    ),
]


def beats(a: Carrier, b: Carrier) -> bool:
    av = a.vector
    bv = b.vector
    greater = sum(x > y for x, y in zip(av, bv))
    less = sum(x < y for x, y in zip(av, bv))
    if greater != less:
        return greater > less
    return CARRIERS.index(a) < CARRIERS.index(b)


def tournament_fingerprint() -> tuple[dict[int, int], int, list[int], int, list[str]]:
    n = len(CARRIERS)
    adj = [[False] * n for _ in range(n)]
    for i, a in enumerate(CARRIERS):
        for j, b in enumerate(CARRIERS):
            if i != j:
                adj[i][j] = beats(a, b)

    scores = [sum(adj[i][j] for j in range(n) if j != i) for i in range(n)]
    score_hist = dict(sorted(Counter(scores).items()))

    directed_3cycles = 0
    for a, b, c in combinations(range(n), 3):
        edges = [adj[a][b], adj[b][c], adj[c][a]]
        redges = [adj[a][c], adj[c][b], adj[b][a]]
        if all(edges) or all(redges):
            directed_3cycles += 1

    # Kosaraju SCCs without external dependencies.
    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w in range(n):
            if adj[v][w] and w not in seen:
                dfs(w)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)

    seen.clear()
    scc_sizes: list[int] = []

    def rdfs(v: int) -> int:
        seen.add(v)
        size = 1
        for w in range(n):
            if adj[w][v] and w not in seen:
                size += rdfs(w)
        return size

    for v in reversed(order):
        if v not in seen:
            scc_sizes.append(rdfs(v))

    # Hamiltonian path count in the carrier tournament.
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    hp = sum(dp[(1 << n) - 1])

    path = [CARRIERS[i].name for i in sorted(range(n), key=lambda i: -scores[i])]
    return score_hist, directed_3cycles, sorted(scc_sizes), hp, path


def print_dfa_tables() -> None:
    print("Finite automata used as carry languages")
    print("=" * 88)
    print("fibbinary DFA, bits read left-to-right after leading 1")
    print("  state Z(last bit 0/start): 0->Z, 1->O")
    print("  state O(last bit 1):       0->Z, 1->DEAD")
    print("  accepted states: Z,O")
    print()
    print("Moser-de Bruijn DFA, bits read least-significant first")
    print("  state E(even bit position): 0->O, 1->O")
    print("  state O(odd bit position):  0->E, 1->DEAD")
    print("  accepted states: E,O after at least one 1")
    print("  note: this language is a strict sparse sublanguage of fibbinary.")


def print_sequence_tables() -> None:
    print()
    print("Sequence growth and residue transversals")
    print("=" * 88)
    fib = generate(is_fibbinary, 13)
    moser = generate(is_moser_de_bruijn, 13)
    fib_trans = first_residue_transversal(is_fibbinary)
    moser_trans = first_residue_transversal(is_moser_de_bruijn)
    print(f"first 13 fibbinary:        {fib}")
    print(f"first 13 Moser-de Bruijn:  {moser}")
    print(f"first fibbinary residues:  {fib_trans}")
    print(f"first Moser residues:      {moser_trans}")
    print()
    print("bit-width counts (exact)")
    print(f"{'width':>5} {'fibbinary':>10} {'Moser':>10}")
    for (w, fc), (_, mc) in zip(fibbinary_count_by_width(12), moser_count_by_width(12)):
        print(f"{w:5d} {fc:10d} {mc:10d}")


def print_lrc_rows() -> None:
    print()
    print("Exact LRC14 row audit")
    print("=" * 88)
    rows = [
        ("AP boundary", tuple(range(1, 14))),
        ("GW boundary 12->24", tuple(list(range(1, 12)) + [13, 24])),
        ("first13 fibbinary", generate(is_fibbinary, 13)),
        ("first13 Moser", generate(is_moser_de_bruijn, 13)),
        ("fibbinary residue transversal", first_residue_transversal(is_fibbinary)),
        ("Moser residue transversal", first_residue_transversal(is_moser_de_bruijn)),
    ]
    print(
        f"{'row':<31} {'max':>5} {'M':>10} {'margin':>10} "
        f"{'best t':>10} {'safe_mass':>12} {'longest':>10} {'comp':>5}"
    )
    records: list[RowRecord] = []
    for name, row in rows:
        rec = row_record(name, row)
        records.append(rec)
        print(
            f"{name:<31} {max(row):5d} {fmt_fraction(rec.M):>10} "
            f"{fmt_fraction(rec.margin):>10} {fmt_fraction(rec.t):>10} "
            f"{fmt_fraction(rec.safe_mass):>12} {fmt_fraction(rec.longest_safe):>10} "
            f"{rec.components:5d}"
        )
        print(f"  residues: {residue_profile(rec.row)}")
        print(f"  row: {rec.row}")
    print()
    print("Readout:")
    print("  * Fibbinary and Moser transversals can restore the mod-14 one-hole support,")
    print("    but their exact M values are strict, so the automata do not create new")
    print("    tight atoms in this scout.")
    print("  * Moser is a useful negative-control: it is extremely finite-state and")
    print("    no-carry, yet too magnitude-lacunary to preserve the LRC14 boundary.")
    print("  * The largest-safe-component column is the circle analogue of the")
    print("    large-stick/potato geometry prompt: it is exact and predicate-preserving.")


def print_carrier_tournament() -> None:
    print()
    print("Tournament Analysis over proof carriers")
    print("=" * 88)
    print(
        f"{'carrier':<34} {'scale':>5} {'LRC':>4} {'DFA':>4} "
        f"{'carry':>5} {'warn':>5} note"
    )
    for c in CARRIERS:
        print(
            f"{c.name:<34} {c.exact_scale:5d} {c.lrc_predicate:4d} "
            f"{c.finite_state:4d} {c.carry_collision:5d} {c.quotient_warning:5d} {c.note}"
        )
    score_hist, cycles, sccs, hp, path = tournament_fingerprint()
    print()
    print(f"score_hist={score_hist}")
    print(f"directed_3cycles={cycles}")
    print(f"scc_sizes={sccs}")
    print(f"hamiltonian_path_count={hp}")
    print("one Hamiltonian path by score:")
    for i, name in enumerate(path, 1):
        print(f"  {i}. {name}")


def print_synthesis() -> None:
    print()
    print("Synthesis for the LRC14 proof stack")
    print("=" * 88)
    print(
        "The Fermat-Catalan analogy should be used as scarcity discipline: "
        "power-like collisions are not a dense certificate source unless the "
        "exact primitive side channels are retained."
    )
    print(
        "The Ostrowski-Hadamard analogy is a warning, not a proof engine: "
        "Hadamard-lacunary supports are naturally boundary-rigid, so analytic "
        "continuation across forgotten packet labels is the wrong expectation."
    )
    print(
        "The actionable finite-automaton import is local: keep a small DFA or "
        "transducer for dyadic/Ostrowski carry states, but attach it to "
        "M=p/q, e=14p-q, residue ownership, and exact safe components."
    )
    print(
        "Next executable target: add automaton_state and carry_language fields "
        "to HYP-2963/HYP-3001 packet records, then test whether zero-open "
        "non-AP/GW residuals emit a nontrivial dyadic carry state or collapse "
        "to the known AP/GW boundary skeleton."
    )


def main() -> None:
    print_dfa_tables()
    print_sequence_tables()
    print_lrc_rows()
    print_carrier_tournament()
    print_synthesis()


if __name__ == "__main__":
    main()
