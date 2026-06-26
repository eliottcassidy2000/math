#!/usr/bin/env python3
"""LRC14 automatic-gap / power-lift ledger scout.

This is a proof-interface script, not a certificate for LRC14.  It tests a
small named packet bank against three sequence languages that the current
proof stack keeps circling:

* Moser-de Bruijn: binary ones only in even positions, equivalently a carry-free
  base-4 square channel with the unique x + 2y representation.
* Fibbinary / Zeckendorf: no adjacent binary ones, a path-independent-set
  channel for endpoint debt.
* Hadamard/Ostrowski gap: lacunary frequency tails whose large ratios should
  not be smoothed away unless their address coordinate is retained.

Tournament Analysis is performed on proof-carrier vertices, not runners.  This
records the challenged assumption required by AGENTS.md: the load-bearing
vertices here are packet languages, gap carriers, and power-lift guards.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations
import math


@dataclass(frozen=True)
class Row:
    name: str
    speeds: tuple[int, ...]
    route_hint: str


ROWS: tuple[Row, ...] = (
    Row("AP13", tuple(range(1, 14)), "AP/GW boundary atom"),
    Row("GW_12_to_24", tuple(list(range(1, 12)) + [13, 24]), "AP/GW boundary atom"),
    Row("K33_12_to_36", tuple(list(range(1, 12)) + [13, 36]), "K33 state-lift / Fejer route"),
    Row("petal_10_to_20", tuple(list(range(1, 10)) + [11, 12, 13, 20]), "C27 petal route"),
    Row("petal_13_to_26", tuple(list(range(1, 13)) + [26]), "C27 petal route"),
    Row("cover_12_to_84", tuple(list(range(1, 12)) + [13, 84]), "covering moment / q-witness route"),
    Row("q27_probe_tail260", (7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 91, 260), "Res_27 carry-owner probe"),
)


def is_fibbinary(n: int) -> bool:
    return n > 0 and (n & (n >> 1)) == 0


def is_moser(n: int) -> bool:
    if n <= 0:
        return False
    bit = 0
    m = n
    while m:
        if (m & 1) and (bit % 2 == 1):
            return False
        bit += 1
        m >>= 1
    return True


def language_letter(n: int) -> str:
    if is_moser(n):
        return "M"
    if is_fibbinary(n):
        return "F"
    return "C"


def zeckendorf_atoms_upto(n: int) -> list[int]:
    atoms = [1, 2]
    while atoms[-1] + atoms[-2] <= n:
        atoms.append(atoms[-1] + atoms[-2])
    return atoms


def zeckendorf_indices(n: int) -> tuple[int, ...]:
    atoms = zeckendorf_atoms_upto(n)
    out: list[int] = []
    rem = n
    for idx in range(len(atoms) - 1, -1, -1):
        a = atoms[idx]
        if a <= rem:
            out.append(idx + 2)  # atoms are F_2=1, F_3=2, ...
            rem -= a
    return tuple(out)


def zeckendorf_gap_word(n: int) -> str:
    idx = zeckendorf_indices(n)
    if len(idx) <= 1:
        return "single"
    gaps = [idx[i] - idx[i + 1] for i in range(len(idx) - 1)]
    return ".".join(str(g) for g in gaps)


def moser_decompose(n: int) -> tuple[int, int]:
    """Return the unique Moser pair x,y with n=x+2y.

    The even binary positions go to x; the odd binary positions are shifted
    down by one and go to y.  Both x and y are in the Moser-de Bruijn sequence.
    """

    x = 0
    y = 0
    bit = 0
    m = n
    while m:
        if m & 1:
            if bit % 2 == 0:
                x |= 1 << bit
            else:
                y |= 1 << (bit - 1)
        bit += 1
        m >>= 1
    return x, y


def integer_nth_root(n: int, e: int) -> int:
    lo, hi = 1, int(round(n ** (1.0 / e))) + 3
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**e <= n:
            lo = mid
        else:
            hi = mid
    return lo


def perfect_power(n: int) -> tuple[int, int] | None:
    if n <= 1:
        return None
    for e in range(int(math.log2(n)) + 1, 1, -1):
        b = integer_nth_root(n, e)
        if b > 1 and b**e == n:
            return b, e
    return None


def hadamard_ratios(row: tuple[int, ...]) -> tuple[Fraction, ...]:
    vals = sorted(row)
    return tuple(Fraction(vals[i + 1], vals[i]) for i in range(len(vals) - 1))


def max_gap_ratio(row: tuple[int, ...]) -> Fraction:
    return max(hadamard_ratios(row))


def tail_gap_ratio(row: tuple[int, ...]) -> Fraction:
    vals = sorted(row)
    return Fraction(vals[-1], vals[-2])


def unit_excess_payloads(pmax: int = 8) -> list[dict[str, object]]:
    out = []
    for p in range(1, pmax + 1):
        q = 14 * p - 1
        payloads = {
            "p": p,
            "q": q,
            "p_plus_q": p + q,
            "p_times_q": p * q,
        }
        powers = {k: perfect_power(v) for k, v in payloads.items() if isinstance(v, int)}
        out.append({"payloads": payloads, "powers": {k: v for k, v in powers.items() if v}})
    return out


def row_summary(row: Row) -> dict[str, object]:
    word = "".join(language_letter(v) for v in row.speeds)
    counts = {letter: word.count(letter) for letter in "MFC"}
    moser_parts = [moser_decompose(v) for v in row.speeds]
    powers = [(v, perfect_power(v)) for v in row.speeds if perfect_power(v)]
    return {
        "name": row.name,
        "route": row.route_hint,
        "word": word,
        "counts": counts,
        "tail_gap": tail_gap_ratio(row.speeds),
        "max_gap": max_gap_ratio(row.speeds),
        "powers": powers,
        "zeck_gaps": [zeckendorf_gap_word(v) for v in row.speeds],
        "moser_x2y": moser_parts,
    }


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, int, int, int, int, int]


CARRIERS: tuple[Carrier, ...] = (
    Carrier("labelled_packet_bank", (5, 5, 5, 5, 5, 5)),
    Carrier("fibbinary_zeckendorf_fsm", (3, 5, 4, 3, 4, 5)),
    Carrier("moser_square_x_plus_2y", (3, 5, 4, 4, 3, 5)),
    Carrier("ostrowski_hadamard_gap", (4, 3, 5, 3, 4, 4)),
    Carrier("fermat_catalan_power_guard", (4, 3, 3, 5, 4, 4)),
    Carrier("hurwitz_doubling_cf_state", (3, 4, 4, 4, 4, 4)),
    Carrier("potato_visibility_guard", (2, 3, 3, 2, 3, 3)),
    Carrier("raw_scalar_speed_word", (1, 1, 1, 1, 1, 1)),
)

SCORE_LABELS = (
    "preserves_lrc_predicate",
    "finite_state_checkable",
    "keeps_packet_labels",
    "guards_power_or_carry_lift",
    "connects_to_existing_routes",
    "resists_scalarization",
)


def carrier_total(c: Carrier) -> int:
    return sum(c.scores)


def orient(a: Carrier, b: Carrier) -> int:
    if carrier_total(a) != carrier_total(b):
        return 1 if carrier_total(a) > carrier_total(b) else -1
    if a.scores != b.scores:
        return 1 if a.scores > b.scores else -1
    return 1 if a.name < b.name else -1


def tournament_edges() -> dict[tuple[str, str], str]:
    edges = {}
    for a, b in combinations(CARRIERS, 2):
        if orient(a, b) > 0:
            edges[(a.name, b.name)] = a.name
        else:
            edges[(a.name, b.name)] = b.name
    return edges


def score_histogram() -> dict[int, int]:
    hist: dict[int, int] = {}
    for c in CARRIERS:
        outdeg = 0
        for d in CARRIERS:
            if c is not d and orient(c, d) > 0:
                outdeg += 1
        hist[outdeg] = hist.get(outdeg, 0) + 1
    return dict(sorted(hist.items()))


def directed_3cycles() -> int:
    count = 0
    for a, b, c in combinations(CARRIERS, 3):
        wins = {
            (a.name, b.name): orient(a, b) > 0,
            (b.name, c.name): orient(b, c) > 0,
            (c.name, a.name): orient(c, a) > 0,
        }
        if all(wins.values()) or not any(wins.values()):
            count += 1
    return count


def is_hamiltonian_path(order: tuple[Carrier, ...]) -> bool:
    return all(orient(order[i], order[i + 1]) > 0 for i in range(len(order) - 1))


def hamiltonian_path_count() -> int:
    return sum(1 for order in permutations(CARRIERS) if is_hamiltonian_path(order))


def hamiltonian_path_example() -> tuple[str, ...]:
    ranked = sorted(CARRIERS, key=lambda c: (carrier_total(c), c.scores, c.name), reverse=True)
    return tuple(c.name for c in ranked)


def print_sources() -> None:
    print("EXTERNAL SOURCE READOUT")
    print("- 2-adic Littlewood paper arXiv:2506.04110: use Hurwitz multiplication-by-2 continued-fraction states as a dyadic packet analogy; reported lower-bound shift C>=15.")
    print("- Ostrowski-Hadamard gap theorem: lacunary exponent support with p_{j+1}/p_j > lambda > 1 should be treated as a boundary/natural-frontier label, not a smooth tail.")
    print("- Moser-de Bruijn sequence: powers-of-4 support / even binary digits / unique x+2y decomposition.")
    print("- Fibbinary numbers: no adjacent binary ones, matching the path independent-set / Zeckendorf endpoint-debt language.")
    print("- Fermat-Catalan input: use power equations as no-lift guards requiring cyclotomic/p-adic labels; do not scalarize perfect-power collisions.")
    print("- ACM DOI 10.5555/1109557.1109610 is the polygon sticks/potatoes paper; here it contributes only an approximation/visibility guardrail.")


def print_rows() -> None:
    print("\nNAMED ROW AUTOMATIC-LANGUAGE LEDGER")
    for row in ROWS:
        s = row_summary(row)
        print(f"{s['name']}: route={s['route']}")
        print(f"  word={s['word']} counts={s['counts']} tail_gap={s['tail_gap']} max_gap={s['max_gap']}")
        print(f"  perfect_power_speeds={s['powers']}")
        print(f"  zeckendorf_gap_words={s['zeck_gaps']}")
        print(f"  moser_x_plus_2y_first_last={s['moser_x2y'][0]} -> {s['moser_x2y'][-1]}")


def print_payloads() -> None:
    print("\nUNIT-EXCESS p/(14p-1) PAYLOAD POWER CHECK")
    for item in unit_excess_payloads(8):
        payloads = item["payloads"]
        powers = item["powers"]
        print(f"p={payloads['p']}: q={payloads['q']} p+q={payloads['p_plus_q']} p*q={payloads['p_times_q']} powers={powers}")


def print_tournament() -> None:
    print("\nTOURNAMENT ANALYSIS OVER PROOF CARRIERS")
    print("vertices_are=proof carriers, not runners")
    print("score_labels=" + ", ".join(SCORE_LABELS))
    for c in sorted(CARRIERS, key=lambda x: carrier_total(x), reverse=True):
        print(f"  {c.name}: total={carrier_total(c)} scores={c.scores}")
    print(f"score_histogram={score_histogram()}")
    print(f"directed_3cycles={directed_3cycles()}")
    print(f"hamiltonian_path_count={hamiltonian_path_count()}")
    print("hamiltonian_path_example=" + " > ".join(hamiltonian_path_example()))


def print_proof_targets() -> None:
    print("\nPROOF-TARGET CONSEQUENCES")
    print("1. Add automatic_language_class, fibbinary_carry_status, moser_even_bit_status, lacunary_gap_ratio, and power_lift_guard to HYP-2963 packets.")
    print("2. A quotient may forget the automaton state only if the LRC predicate is fiber-constant, the state is reconstructible, a dual certificate annihilates it, or the packet is routed to AP/GW, C27, K33, Res_27, or covering residual debt.")
    print("3. Moser/even-bit support is stricter than fibbinary/no-adjacent support: use it as a no-hidden-square-carry subfiber.")
    print("4. Hadamard-lacunary tails should be certified by a boundary/frontier label before Fejer or Haar smoothing treats them as harmless.")
    print("5. Fermat-Catalan style perfect-power payloads are guards against illegal lifts, not standalone LRC certificates.")


def main() -> None:
    print("LRC14 AUTOMATIC GAP / POWER-LIFT LEDGER SCOUT (S178)")
    print_sources()
    print_rows()
    print_payloads()
    print_tournament()
    print_proof_targets()


if __name__ == "__main__":
    main()
