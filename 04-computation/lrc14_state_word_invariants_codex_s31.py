#!/usr/bin/env python3
"""Scout state-word invariants for LRC14 rows.

The computation treats a row E as a weighted cyclic word over the 6 inner
sectors: each atom of the wall arrangement records which sectors are missed by
the multiples of E.  Scalar quantities such as p0, L_y, far-element plateaus,
fold profiles, and signed AP-to-defect transports are then projections of this
state word or of couplings between two state words.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import log2
from typing import Iterable


SectorState = tuple[int, ...]
Atom = tuple[Fraction, Fraction, Fraction, SectorState]
Row = tuple[int, ...]


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def sector(x: Fraction) -> int:
    """Return the 0..6 sector for x in [0,1)."""

    y = frac_part(x)
    return (y.numerator * 7) // y.denominator


def missed_state(row: Row, x: Fraction) -> SectorState:
    seen = {sector(x * r) for r in row}
    return tuple(s for s in range(1, 7) if s not in seen)


def breakpoints(row: Row) -> list[Fraction]:
    bps = {Fraction(0), Fraction(1)}
    for r in row:
        if r == 0:
            continue
        for a in range(0, 7 * r + 1):
            bps.add(Fraction(a, 7 * r))
    return sorted(bps)


def state_atoms(row: Row) -> list[Atom]:
    bps = breakpoints(row)
    atoms: list[Atom] = []
    for lo, hi in zip(bps, bps[1:]):
        if lo == hi:
            continue
        mid = (lo + hi) / 2
        atoms.append((lo, hi, hi - lo, missed_state(row, mid)))
    return atoms


def g_weight(k: int, t: int) -> Fraction:
    if k == 8:
        return Fraction((t - 1) * (t - 2) * (t - 4) * (t - 5), 40)
    if k in (9, 10):
        return Fraction(-(t - 2) * (t - 3) * (t - 6), 36)
    raise ValueError(f"no g-weight formula recorded for k={k}")


def fmtq(q: Fraction, digits: int = 8) -> str:
    return f"{q} ({float(q):.{digits}f})"


def fmtf(x: float, digits: int = 5) -> str:
    return f"{x:.{digits}f}"


def hamming_state(a: SectorState, b: SectorState) -> int:
    return len(set(a).symmetric_difference(b))


def entropy(counter: Counter[SectorState] | Counter[tuple[SectorState, SectorState]]) -> float:
    out = 0.0
    for mass in counter.values():
        p = float(mass)
        if p:
            out -= p * log2(p)
    return out


def visible_fold_profile(row: Row) -> Counter[int]:
    """Nontrivial pair-sum folds by target c: 0<a<b and a+b=c in row."""

    out: Counter[int] = Counter()
    present = set(row)
    for a, b in combinations(row, 2):
        if a == 0:
            continue
        c = a + b
        if c in present:
            out[c] += 1
    return out


@dataclass(frozen=True)
class StateSignature:
    label: str
    row: Row
    ly: Fraction
    p_by_count: tuple[Fraction, ...]
    plateau: Fraction
    state_mass: Counter[SectorState]
    transition_jumps: Counter[int]
    single_sector_mass: tuple[Fraction, ...]
    fold_target_count: int
    fold_reciprocal: Fraction

    @property
    def p0(self) -> Fraction:
        return self.p_by_count[0]

    @property
    def p1(self) -> Fraction:
        return self.p_by_count[1]

    @property
    def state_entropy(self) -> float:
        return entropy(self.state_mass)

    @property
    def effective_states(self) -> float:
        return 2 ** self.state_entropy

    @property
    def single_tv(self) -> Fraction:
        mean = self.p1 / 6
        return sum(abs(v - mean) for v in self.single_sector_mass) / 2

    @property
    def support(self) -> int:
        return len(self.state_mass)

    @property
    def jumps(self) -> int:
        return sum(self.transition_jumps.values())


def state_signature(label: str, row_iter: Iterable[int]) -> StateSignature:
    row = tuple(sorted(row_iter))
    atoms = state_atoms(row)
    by_count = [Fraction(0) for _ in range(7)]
    state_mass: Counter[SectorState] = Counter()
    single = [Fraction(0) for _ in range(6)]
    for _lo, _hi, mass, state in atoms:
        by_count[len(state)] += mass
        state_mass[state] += mass
        if len(state) == 1:
            single[state[0] - 1] += mass
    k = len(row)
    ly = sum(by_count[t] * g_weight(k, t) for t in range(7))
    jumps: Counter[int] = Counter()
    for a, b in zip(atoms, atoms[1:] + atoms[:1]):
        if a[3] != b[3]:
            jumps[hamming_state(a[3], b[3])] += 1
    folds = visible_fold_profile(row)
    fold_reciprocal = sum(Fraction(v, target) for target, v in folds.items())
    return StateSignature(
        label=label,
        row=row,
        ly=ly,
        p_by_count=tuple(by_count),
        plateau=by_count[0] + by_count[1] / 7,
        state_mass=state_mass,
        transition_jumps=jumps,
        single_sector_mass=tuple(single),
        fold_target_count=sum(folds.values()),
        fold_reciprocal=fold_reciprocal,
    )


@dataclass(frozen=True)
class TransportSignature:
    label: str
    row_a: Row
    row_b: Row
    positive: Fraction
    negative: Fraction
    neutral: Fraction
    signed: Fraction
    pair_mass: Counter[tuple[SectorState, SectorState]]
    count_delta_mass: Counter[int]

    @property
    def support(self) -> int:
        return len(self.pair_mass)

    @property
    def pair_entropy(self) -> float:
        return entropy(self.pair_mass)


def transport_signature(label: str, row_a_iter: Iterable[int], row_b_iter: Iterable[int]) -> TransportSignature:
    row_a = tuple(sorted(row_a_iter))
    row_b = tuple(sorted(row_b_iter))
    k = len(row_a)
    assert len(row_a) == len(row_b)
    bps = sorted(set(breakpoints(row_a)).union(breakpoints(row_b)))
    positive = Fraction(0)
    negative = Fraction(0)
    neutral = Fraction(0)
    pair_mass: Counter[tuple[SectorState, SectorState]] = Counter()
    count_delta_mass: Counter[int] = Counter()
    for lo, hi in zip(bps, bps[1:]):
        if lo == hi:
            continue
        mid = (lo + hi) / 2
        state_a = missed_state(row_a, mid)
        state_b = missed_state(row_b, mid)
        mass = hi - lo
        delta = g_weight(k, len(state_b)) - g_weight(k, len(state_a))
        pair_mass[(state_a, state_b)] += mass
        count_delta_mass[len(state_b) - len(state_a)] += mass
        if delta > 0:
            positive += mass * delta
        elif delta < 0:
            negative += mass * (-delta)
        else:
            neutral += mass
    return TransportSignature(
        label=label,
        row_a=row_a,
        row_b=row_b,
        positive=positive,
        negative=negative,
        neutral=neutral,
        signed=positive - negative,
        pair_mass=pair_mass,
        count_delta_mass=count_delta_mass,
    )


def top_states(sig: StateSignature, n: int = 5) -> str:
    items = sig.state_mass.most_common(n)
    return "; ".join(f"{state}:{mass}" for state, mass in items)


def top_pairs(sig: TransportSignature, n: int = 5) -> str:
    items = sig.pair_mass.most_common(n)
    return "; ".join(f"{a}->{b}:{mass}" for (a, b), mass in items)


def histogram_string(counter: Counter[int]) -> str:
    return ", ".join(f"{key}:{counter[key]}" for key in sorted(counter))


def print_state_table(signatures: list[StateSignature]) -> None:
    print("single-row state signatures")
    print(
        "label | k | L_y | p0 | p1 | p0+p1/7 | states | H | eff | "
        "jumps | jump types | p1_TV | folds | fold_rec"
    )
    for sig in signatures:
        print(
            f"{sig.label} | {len(sig.row)} | {fmtq(sig.ly, 7)} | "
            f"{fmtq(sig.p0, 7)} | {fmtq(sig.p1, 7)} | {fmtq(sig.plateau, 7)} | "
            f"{sig.support} | {fmtf(sig.state_entropy)} | {fmtf(sig.effective_states)} | "
            f"{sig.jumps} | {histogram_string(sig.transition_jumps)} | "
            f"{fmtq(sig.single_tv, 7)} | {sig.fold_target_count} | {sig.fold_reciprocal}"
        )
        print(f"  top states: {top_states(sig)}")
        print(f"  single-sector mass: {[str(v) for v in sig.single_sector_mass]}")


def print_transport_table(transports: list[TransportSignature]) -> None:
    print()
    print("common-wall state transports")
    print(
        "label | signed | positive | negative | neutral measure | support pairs | "
        "pair H | count-delta mass"
    )
    for sig in transports:
        print(
            f"{sig.label} | {fmtq(sig.signed, 8)} | {fmtq(sig.positive, 8)} | "
            f"{fmtq(sig.negative, 8)} | {fmtq(sig.neutral, 8)} | "
            f"{sig.support} | {fmtf(sig.pair_entropy)} | "
            f"{histogram_string(sig.count_delta_mass)}"
        )
        print(f"  top pairs: {top_pairs(sig)}")


def main() -> None:
    rows = [
        ("AP8", range(8)),
        ("AP9", range(9)),
        ("AP10", range(10)),
        ("D9_top_plus1", (0, 1, 2, 3, 4, 5, 6, 7, 9)),
        ("D9_top_plus2", (0, 1, 2, 3, 4, 5, 6, 7, 10)),
        ("M9_spectrum_2k", (1, 2, 3, 4, 5, 6, 7, 8, 18)),
        ("D9_gap_8", (0, 1, 2, 3, 4, 5, 6, 8, 9)),
        ("odd9", (1, 3, 5, 7, 9, 11, 13, 15, 17)),
        ("sporadic_1347", (0, 1, 3, 4, 7, 8, 9, 10, 11)),
        ("kps_pocket", (3, 5, 16, 28, 30, 33, 35, 36, 37)),
    ]
    signatures = [state_signature(label, row) for label, row in rows]
    transports = [
        transport_signature("AP9 -> D9_top_plus1", range(9), (0, 1, 2, 3, 4, 5, 6, 7, 9)),
        transport_signature("AP9 -> D9_top_plus2", range(9), (0, 1, 2, 3, 4, 5, 6, 7, 10)),
        transport_signature("AP9 -> M9_spectrum_2k", range(9), (1, 2, 3, 4, 5, 6, 7, 8, 18)),
        transport_signature("AP9 -> D9_gap_8", range(9), (0, 1, 2, 3, 4, 5, 6, 8, 9)),
        transport_signature("AP9 -> odd9", range(9), (1, 3, 5, 7, 9, 11, 13, 15, 17)),
    ]
    print("HYP-2648 / T895 state-word invariant scout")
    print()
    print_state_table(signatures)
    print_transport_table(transports)
    print()
    print("invariant reading")
    print("1. The measured cyclic state word determines p0, L_y, p0+p1/7, single-sector bias, fold")
    print("   reciprocals, and all common-wall transfer shadows by quotienting or coupling.")
    print("2. Neutral measure can dominate a transfer while the signed result is decided by a small")
    print("   addressed residue of state pairs.  This matches HYP-2647's warning that balanced")
    print("   moving-sector marginals do not determine the obstruction.")
    print("3. The far-element plateau quotient only retains (p0,p1).  It is useful because it is")
    print("   monotone under KPS-style x-cell tails, but it discards transition adjacency and")
    print("   sector labels that the AP-to-defect wall proof needs.")
    print("4. Fold targets and coset/reciprocal signatures should be treated as addresses on the")
    print("   same finite state object, not as independent explanations.")
    print()
    print("tournament-analysis quotient")
    print("vertices: measured state word, signed state transport, missed-set histogram,")
    print("p0/p1 plateau, fold-target profile, mod-7 reciprocal phase, raw relation rank, raw runners")
    print("edge rule: A -> B when B is a deterministic quotient or valuation of A")
    print("Hamiltonian path: state word -> transport -> histogram -> plateau -> fold profile ->")
    print("mod-7 phase -> relation rank -> runners")
    print("challenged assumption: tournament vertices need not be runners or arcs; here they are")
    print("finite addresses and quotient obligations that preserve the LRC bad-set measure.")


if __name__ == "__main__":
    main()
