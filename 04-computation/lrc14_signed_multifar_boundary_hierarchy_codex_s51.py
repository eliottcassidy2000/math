#!/usr/bin/env python3
"""HYP-2680 / T919: signed multi-far boundary hierarchy for LRC(14).

This scout extends THM-548/HYP-2679 from the two-far curvature to the signed
three-far and higher Newton differences.  It keeps the direct sector-cover
quantity from HYP-2679:

    p0(E) = measure{x : all six inner sectors 1..6 are hit by E*x}.

The other "p0" convention in the repo, the all-runners-survive danger-arc
measure, is not used here.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import comb, gcd

from lrc14_wide_branch_ridge_codex_s47 import (
    CAP,
    additive_energy,
    missed_distribution,
    p0,
    primitive,
    sumset_excess,
)


Row = tuple[int, ...]


def fmt(q: F) -> str:
    return f"{q} ({float(q):.9f})"


@lru_cache(maxsize=None)
def p0_cached(row: Row) -> F:
    return p0(tuple(sorted(set(row))))


@lru_cache(maxsize=None)
def profile_cached(row: Row) -> tuple[F, ...]:
    return missed_distribution(tuple(sorted(set(row))))


@lru_cache(maxsize=None)
def stirling2(n: int, k: int) -> int:
    if n == k == 0:
        return 1
    if n == 0 or k == 0 or k > n:
        return 0
    if k == 1 or k == n:
        return 1
    return k * stirling2(n - 1, k) + stirling2(n - 1, k - 1)


def phi_coeffs(s: int) -> list[int]:
    """Integer coefficients c_t with Phi_s = 7^-s sum_t c_t p_t."""

    return [((-1) ** (s - t)) * factorial(t) * stirling2(s, t) for t in range(1, s + 1)]


@lru_cache(maxsize=None)
def factorial(n: int) -> int:
    out = 1
    for a in range(2, n + 1):
        out *= a
    return out


def phi_s(core: Row, s: int) -> F:
    prof = profile_cached(core)
    total = F(0)
    for t in range(1, min(6, s) + 1):
        coeff = ((-1) ** (s - t)) * factorial(t) * stirling2(s, t)
        total += coeff * prof[t]
    return total / (7**s)


def c_t(t: int, r: int) -> F:
    return sum(((-1) ** i) * comb(t, i) * F(7 - i, 7) ** r for i in range(t + 1))


def boundary_value_direct(core: Row, r: int) -> F:
    prof = profile_cached(core)
    return sum(prof[t] * c_t(t, r) for t in range(7))


def boundary_value_newton(core: Row, r: int) -> F:
    return profile_cached(core)[0] + sum(F(comb(r, s)) * phi_s(core, s) for s in range(1, r + 1))


def delta_s(core: Row, far: Row) -> F:
    far = tuple(far)
    s = len(far)
    total = F(0)
    for q in range(s + 1):
        sign = (-1) ** (s - q)
        for sub in combinations(far, q):
            total += sign * p0_cached(tuple(sorted(core + sub)))
    return total


def relation_distance(vals: Row, height: int = 4) -> tuple[int, tuple[int, ...], int]:
    """Return min |c.vals|, primitive coefficient vector, and l1 height."""

    best_val: int | None = None
    best_coeff: tuple[int, ...] | None = None
    best_l1: int | None = None
    for coeff in product(range(-height, height + 1), repeat=len(vals)):
        if all(c == 0 for c in coeff):
            continue
        g = 0
        for c in coeff:
            g = gcd(g, abs(c))
        if g != 1:
            continue
        val = abs(sum(c * v for c, v in zip(coeff, vals)))
        l1 = sum(abs(c) for c in coeff)
        key = (val, l1, coeff)
        if best_val is None or key < (best_val, best_l1, best_coeff):
            best_val = val
            best_coeff = coeff
            best_l1 = l1
    assert best_val is not None and best_coeff is not None and best_l1 is not None
    return best_val, best_coeff, best_l1


def pair_relation_distance(vals: Row, height: int = 4) -> tuple[int, tuple[int, int], tuple[int, int]]:
    best = None
    for i, j in combinations(range(len(vals)), 2):
        d, coeff, l1 = relation_distance((vals[i], vals[j]), height)
        item = (d, l1, (i, j), coeff)
        if best is None or item < best:
            best = item
    assert best is not None
    d, _l1, pair, coeff = best
    return d, pair, coeff


@dataclass(frozen=True)
class TripleReport:
    core: Row
    far: Row
    delta: F
    phi: F
    full_p0: F

    @property
    def row(self) -> Row:
        return tuple(sorted(self.core + self.far))

    @property
    def dev(self) -> F:
        return self.delta - self.phi

    @property
    def abs_dev(self) -> F:
        return abs(self.dev)

    @property
    def cap(self) -> F | None:
        return CAP.get(len(self.row))

    @property
    def margin(self) -> F | None:
        cap = self.cap
        return None if cap is None else cap - self.full_p0


def make_triple(core: Row, far: Row) -> TripleReport:
    core = tuple(sorted(core))
    far = tuple(sorted(far))
    return TripleReport(
        core=core,
        far=far,
        delta=delta_s(core, far),
        phi=phi_s(core, len(far)),
        full_p0=p0_cached(tuple(sorted(core + far))),
    )


def scan_core_triples(core: Row, far_values: range, keep: int = 8) -> dict[str, list[TripleReport]]:
    tops: dict[str, list[TripleReport]] = {
        "abs_dev": [],
        "positive_dev": [],
        "negative_dev": [],
        "p0": [],
    }

    def push(name: str, rep: TripleReport, key) -> None:
        tops[name].append(rep)
        tops[name].sort(key=key, reverse=True)
        del tops[name][keep:]

    for far in combinations(far_values, 3):
        rep = make_triple(core, far)
        push("abs_dev", rep, lambda r: (r.abs_dev, -relation_distance(r.far)[0], r.full_p0, r.far))
        if rep.dev > 0:
            push("positive_dev", rep, lambda r: (r.dev, r.full_p0, r.far))
        if rep.dev < 0:
            push("negative_dev", rep, lambda r: (abs(r.dev), r.full_p0, r.far))
        push("p0", rep, lambda r: (r.full_p0, r.abs_dev, r.far))
    return tops


def print_coeff_table(max_s: int = 8) -> None:
    print("SIGNED NEWTON / STIRLING COEFFICIENTS")
    print("Phi_s(B) = 7^-s * sum_t coeff[s,t] * p_t(B)")
    for s in range(1, max_s + 1):
        coeff = [((-1) ** (s - t)) * factorial(t) * stirling2(s, t) for t in range(1, min(6, s) + 1)]
        print(f"  s={s}: denom=7^{s}, coeff_t1..t{min(6,s)}={coeff}")
    print()


def verify_boundary_identities() -> None:
    print("BOUNDARY VALUE IDENTITY CHECKS")
    cores = [
        (0, 4, 6, 8, 10, 12, 14),
        (0, 1, 2, 3, 4, 5, 6, 7),
        (0, 1, 2, 4, 8, 9, 11),
    ]
    for core in cores:
        print(f"core={core}")
        print(f"  profile={profile_cached(core)}")
        for r in range(1, 7):
            direct = boundary_value_direct(core, r)
            newton = boundary_value_newton(core, r)
            print(f"  r={r}: P_r direct={fmt(direct)} ; Newton sum={fmt(newton)} ; equal={direct == newton}")
    print()


def print_triple(label: str, rep: TripleReport) -> None:
    trip_d, trip_coeff, trip_l1 = relation_distance(rep.far, 4)
    pair_d, pair, pair_coeff = pair_relation_distance(rep.far, 4)
    cap_text = "n/a" if rep.cap is None else fmt(rep.cap)
    margin_text = "n/a" if rep.margin is None else fmt(rep.margin)
    print(f"{label}: core={rep.core}, far={rep.far}, row={rep.row}")
    print(f"  Delta_3={fmt(rep.delta)}")
    print(f"  Phi_3={fmt(rep.phi)}")
    print(f"  deviation={fmt(rep.dev)} abs={fmt(rep.abs_dev)}")
    print(f"  p0(full)={fmt(rep.full_p0)} cap={cap_text} margin={margin_text}")
    print(f"  triple_resdist={trip_d} coeff={trip_coeff} l1={trip_l1}; pair_resdist={pair_d} pair={pair} coeff={pair_coeff}")
    print(f"  row_exc={sumset_excess(rep.row)} energy={additive_energy(rep.row)} primitive={primitive(rep.row)}")


def named_triple_reports() -> None:
    print("NAMED THREE-FAR REPORTS")
    named = [
        ("dilated-core consecutive triple", (0, 4, 6, 8, 10, 12, 14), (15, 16, 17)),
        ("dilated-core AP triple", (0, 4, 6, 8, 10, 12, 14), (15, 18, 21)),
        ("consec8 consecutive triple", (0, 1, 2, 3, 4, 5, 6, 7), (15, 16, 17)),
        ("consec8 separated triple", (0, 1, 2, 3, 4, 5, 6, 7), (17, 23, 31)),
        ("third-pocket first triple", (0, 3, 5), (16, 28, 30)),
        ("third-pocket active triple", (0, 3, 5, 16, 28), (30, 33, 35)),
    ]
    for label, core, far in named:
        print_triple(label, make_triple(core, far))
    print()


def selected_core_scans() -> list[TripleReport]:
    print("SELECTED-CORE FAR-TRIPLE SCANS (far 15..22)")
    selected = [
        ("dilated_core", (0, 4, 6, 8, 10, 12, 14)),
        ("consec8", (0, 1, 2, 3, 4, 5, 6, 7)),
        ("mixed_sparse", (0, 1, 2, 4, 8, 9, 11)),
    ]
    global_top: list[TripleReport] = []
    far_values = range(15, 23)
    for label, core in selected:
        tops = scan_core_triples(core, far_values, keep=5)
        print(f"core label={label}, core={core}, triples={comb(len(tuple(far_values)),3)}")
        for name in ("abs_dev", "positive_dev", "negative_dev", "p0"):
            print(f"  top {name}:")
            for rep in tops[name][:3]:
                trip_d, trip_coeff, _ = relation_distance(rep.far, 4)
                print(
                    f"    far={rep.far} dev={rep.dev} abs={rep.abs_dev} "
                    f"p0={rep.full_p0} margin={rep.margin if rep.margin is not None else 'n/a'} res={trip_d}:{trip_coeff}"
                )
        global_top.extend(tops["abs_dev"])
    print()
    return sorted(global_top, key=lambda r: (r.abs_dev, r.full_p0, r.far), reverse=True)[:8]


def small_exhaustive_bank() -> list[TripleReport]:
    print("SMALL EXACT BANK: all 7-core rows in [0,14], far triple (15,16,17)")
    far_values = range(15, 18)
    top_abs: list[TripleReport] = []
    top_p0: list[TripleReport] = []
    stats = Counter()
    raw = primitive_count = 0

    def push(store: list[TripleReport], rep: TripleReport, key, keep: int = 10) -> None:
        store.append(rep)
        store.sort(key=key, reverse=True)
        del store[keep:]

    for rest in combinations(range(1, 15), 6):
        core = (0,) + rest
        for far in combinations(far_values, 3):
            raw += 1
            row = tuple(sorted(core + far))
            if not primitive(row):
                continue
            primitive_count += 1
            rep = make_triple(core, far)
            trip_d, _coeff, _l1 = relation_distance(rep.far, 4)
            pair_d, _pair, _pc = pair_relation_distance(rep.far, 4)
            stats["triple_exact_relation"] += int(trip_d == 0)
            stats["pair_exact_relation"] += int(pair_d == 0)
            stats["dev_positive"] += int(rep.dev > 0)
            stats["dev_negative"] += int(rep.dev < 0)
            stats["dev_zero"] += int(rep.dev == 0)
            push(top_abs, rep, lambda r: (r.abs_dev, -relation_distance(r.far)[0], r.full_p0, r.row))
            push(top_p0, rep, lambda r: (r.full_p0, r.abs_dev, r.row))

    print(f"  raw={raw} primitive={primitive_count}")
    print(f"  stats={dict(stats)}")
    print("  top abs deviations:")
    for rep in top_abs[:8]:
        trip_d, trip_coeff, _ = relation_distance(rep.far, 4)
        print(f"    row={rep.row} Delta3={rep.delta} Phi3={rep.phi} dev={rep.dev} abs={rep.abs_dev} res={trip_d}:{trip_coeff} p0={rep.full_p0} margin={rep.margin if rep.margin is not None else 'n/a'}")
    print("  top direct p0:")
    for rep in top_p0[:8]:
        trip_d, trip_coeff, _ = relation_distance(rep.far, 4)
        print(f"    row={rep.row} p0={rep.full_p0} margin={rep.margin if rep.margin is not None else 'n/a'} dev={rep.dev} abs={rep.abs_dev} res={trip_d}:{trip_coeff}")
    print()
    return top_abs


def tournament_analysis(reports: list[TripleReport]) -> None:
    print("TOURNAMENT ANALYSIS ON THREE-FAR PROOF OBLIGATIONS")
    # Vertices are top deviation rows/proof obligations, not runners.
    reps = reports[:8]
    n = len(reps)
    wins = [0] * n
    cycles = 0

    def beats(a: TripleReport, b: TripleReport) -> bool:
        ka = (a.abs_dev, -relation_distance(a.far)[0], a.full_p0, a.row)
        kb = (b.abs_dev, -relation_distance(b.far)[0], b.full_p0, b.row)
        return ka > kb

    for i, j in combinations(range(n), 2):
        if beats(reps[i], reps[j]):
            wins[i] += 1
        else:
            wins[j] += 1

    for i, j, k in combinations(range(n), 3):
        eij = beats(reps[i], reps[j])
        ejk = beats(reps[j], reps[k])
        eki = beats(reps[k], reps[i])
        if (eij and ejk and eki) or ((not eij) and (not ejk) and (not eki)):
            cycles += 1

    order = sorted(range(n), key=lambda i: wins[i], reverse=True)
    print("  vertices: top three-far deviation rows")
    print("  observable: larger |Delta_3-Phi_3| wins; tie by smaller relation distance, then p0")
    print(f"  score_hist={dict(Counter(wins))}")
    print(f"  directed_3cycles={cycles}")
    print("  Hamiltonian path:")
    for rank, idx in enumerate(order, 1):
        rep = reps[idx]
        d, coeff, _ = relation_distance(rep.far, 4)
        print(f"    {rank}. far={rep.far} core={rep.core} abs_dev={rep.abs_dev} res={d}:{coeff} p0={rep.full_p0}")
    print()


def main() -> None:
    print("HYP-2680 / T919 -- signed multi-far boundary hierarchy\n")
    print_coeff_table(8)
    verify_boundary_identities()
    named_triple_reports()
    selected_top = selected_core_scans()
    bank_top = small_exhaustive_bank()
    tournament_analysis(sorted(selected_top + bank_top, key=lambda r: (r.abs_dev, r.full_p0, r.row), reverse=True))
    print("SYNTHESIS")
    print("  Phi_3=(p1-6p2+6p3)/343 is the correct signed three-far boundary target.")
    print("  Exact deviations concentrate on low relation-distance triples; AP triples like (15,18,21)")
    print("  expose three-body relations that pairwise screening alone can miss.")
    print("  The next proof step is a signed Abel bound for Delta_3-Phi_3 stratified by")
    print("  the full relation lattice, then induction to higher Newton orders.")


if __name__ == "__main__":
    main()
