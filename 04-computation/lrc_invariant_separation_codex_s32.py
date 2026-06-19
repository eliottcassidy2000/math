#!/usr/bin/env python3
"""Invariant separation scout for Lonely Runner structure.

This complements HYP-2648's state-word invariant.  It runs two small exact
labs:

1. Max-min lab for M(S)=max_t min_v ||v t||.
   We compare coarse speed invariants against the optimal-time clearance word
   (denominator, folded residues, active runners, source crossings).

2. LRC14 sector lab for offset rows E.
   We compare additive summaries/fold profiles against the measured missed
   sector state word used by HYP-2648.

Readout: the structure is determined only after retaining an address on the
wall/crossing arrangement.  Scalar shadows are useful routing labels, not
determinants.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from math import gcd, log2


def set_gcd(S: tuple[int, ...]) -> int:
    g = 0
    for v in S:
        g = gcd(g, v)
    return g


def fold_distance(v: int, t: F) -> int:
    num, den = t.numerator, t.denominator
    r = (v * num) % den
    return min(r, den - r)


def candidate_times(S: tuple[int, ...]) -> set[F]:
    cands: set[F] = set()
    for i, v in enumerate(S):
        for a in range(v):
            cands.add(F(2 * a + 1, 2 * v))
        for w in S[i + 1 :]:
            for den in (v + w, abs(v - w)):
                if den == 0:
                    continue
                for num in range(1, den):
                    cands.add(F(num, den))
    return {t for t in cands if 0 < t < 1}


def exact_M(S: tuple[int, ...]) -> tuple[F, tuple[F, ...]]:
    best = F(0)
    times: list[F] = []
    for t in candidate_times(S):
        den = t.denominator
        val = F(min(fold_distance(v, t) for v in S), den)
        if val > best:
            best = val
            times = [t]
        elif val == best:
            times.append(t)
    return best, tuple(sorted(times))


def crossing_sources(S: tuple[int, ...], t: F) -> tuple[str, ...]:
    src: set[str] = set()
    for i, v in enumerate(S):
        if (t * 2 * v).denominator == 1 and int(t * 2 * v) % 2 == 1:
            src.add(f"peak:{v}")
        for w in S[i + 1 :]:
            if ((v + w) * t).denominator == 1:
                src.add(f"sum:{v}+{w}")
            if v != w and (abs(v - w) * t).denominator == 1:
                src.add(f"diff:{v}-{w}")
    return tuple(sorted(src))


@dataclass(frozen=True)
class ClearanceRecord:
    S: tuple[int, ...]
    M: F
    t: F
    folds: tuple[int, ...]
    active: tuple[int, ...]
    sources: tuple[str, ...]

    @property
    def q(self) -> int:
        return self.t.denominator

    @property
    def j(self) -> int:
        return min(self.folds)

    @property
    def clearance_word(self) -> tuple[int, tuple[int, ...]]:
        return (self.q, self.folds)

    @property
    def active_word(self) -> tuple[int, int, int, tuple[int, ...]]:
        return (self.q, self.j, len(self.active), tuple(sorted(self.active)))


def clearance_record(S: tuple[int, ...]) -> ClearanceRecord:
    M, times = exact_M(S)
    # Canonicalize by smallest denominator, then smallest numerator.
    t = min(times, key=lambda x: (x.denominator, x.numerator))
    folds_by_runner = tuple(fold_distance(v, t) for v in S)
    j = min(folds_by_runner)
    active = tuple(v for v, d in zip(S, folds_by_runner) if d == j)
    return ClearanceRecord(
        S=S,
        M=M,
        t=t,
        folds=tuple(sorted(folds_by_runner)),
        active=active,
        sources=crossing_sources(S, t),
    )


def speed_features(S: tuple[int, ...]) -> dict[str, object]:
    sums = {a + b for a in S for b in S if a <= b}
    folds = Counter()
    present = set(S)
    for a, b in combinations(S, 2):
        c = a + b
        if c in present:
            folds[c] += 1
    diffs = tuple(b - a for a, b in zip(S, S[1:]))
    k = len(S)
    return {
        "sumset_excess": len(sums) - (2 * k - 1),
        "fold_count": sum(folds.values()),
        "fold_profile": tuple(sorted(folds.items())),
        "gap_pattern": diffs,
    }


def mixed_fiber_report(records: list[ClearanceRecord], feature_name: str, feature_fn) -> tuple[int, int, str]:
    fibers: dict[object, list[ClearanceRecord]] = defaultdict(list)
    for rec in records:
        fibers[feature_fn(rec)].append(rec)
    mixed = [rows for rows in fibers.values() if len({r.M for r in rows}) > 1]
    if not mixed:
        return 0, max((len(v) for v in fibers.values()), default=0), "no mixed fibers"
    rows = max(mixed, key=len)
    by_m: dict[F, ClearanceRecord] = {}
    for rec in rows:
        by_m.setdefault(rec.M, rec)
        if len(by_m) >= 2:
            break
    ex = " ; ".join(f"S={r.S}, M={r.M}, t={r.t}, active={r.active}" for r in by_m.values())
    return len(mixed), len(rows), ex


def run_maxmin_lab() -> None:
    print("=" * 88)
    print("A. Exact max-min clearance invariant lab")
    print("=" * 88)
    bank: list[tuple[int, int]] = [(5, 13), (6, 11)]
    records: list[ClearanceRecord] = []
    for k, B in bank:
        for S in combinations(range(1, B + 1), k):
            if set_gcd(S) == 1:
                records.append(clearance_record(tuple(S)))
    print(f"records: {len(records)} primitive speed sets from banks {bank}")
    print("feature | mixed M fibers | largest mixed fiber | example")
    feature_table = [
        ("sumset_excess", lambda r: speed_features(r.S)["sumset_excess"]),
        ("fold_count", lambda r: speed_features(r.S)["fold_count"]),
        ("fold_profile", lambda r: speed_features(r.S)["fold_profile"]),
        ("gap_pattern", lambda r: speed_features(r.S)["gap_pattern"]),
        ("optimal_denominator_q", lambda r: r.q),
        ("q_and_active_count", lambda r: (r.q, len(r.active))),
        ("q_j_active_values", lambda r: r.active_word),
        ("clearance_word_q_folds", lambda r: r.clearance_word),
    ]
    for name, fn in feature_table:
        mixed_count, largest, example = mixed_fiber_report(records, name, fn)
        print(f"{name:24s} | {mixed_count:14d} | {largest:19d} | {example}")

    print()
    print("sample clearance words near the AP floor")
    close = sorted(records, key=lambda r: (float(r.M - F(1, len(r.S) + 1)), max(r.S)))[:10]
    for rec in close:
        print(
            f"S={rec.S}, M={rec.M}, t={rec.t}, q={rec.q}, j={rec.j}, "
            f"folds={rec.folds}, active={rec.active}, sources={rec.sources}"
        )


SectorState = tuple[int, ...]
StateWord = tuple[tuple[F, SectorState], ...]


def frac_part(x: F) -> F:
    return x - (x.numerator // x.denominator)


def sector(x: F) -> int:
    y = frac_part(x)
    return (y.numerator * 7) // y.denominator


def missed_state(row: tuple[int, ...], x: F) -> SectorState:
    seen = {sector(e * x) for e in row}
    return tuple(s for s in range(1, 7) if s not in seen)


def breakpoints(row: tuple[int, ...]) -> list[F]:
    bps = {F(0), F(1)}
    for e in row:
        if e == 0:
            continue
        for a in range(0, 7 * e + 1):
            bps.add(F(a, 7 * e))
    return sorted(bps)


def state_word(row: tuple[int, ...]) -> StateWord:
    bps = breakpoints(row)
    word: list[tuple[F, SectorState]] = []
    for lo, hi in zip(bps, bps[1:]):
        if lo == hi:
            continue
        mid = (lo + hi) / 2
        word.append((hi - lo, missed_state(row, mid)))
    return tuple(word)


def state_mass(word: StateWord) -> tuple[tuple[SectorState, F], ...]:
    c: Counter[SectorState] = Counter()
    for mass, state in word:
        c[state] += mass
    return tuple(sorted(c.items()))


def missed_histogram(word: StateWord) -> tuple[F, ...]:
    out = [F(0) for _ in range(7)]
    for mass, state in word:
        out[len(state)] += mass
    return tuple(out)


def g_weight(k: int, t: int) -> F:
    if k == 8:
        return F((t - 1) * (t - 2) * (t - 4) * (t - 5), 40)
    if k in (9, 10):
        return F(-(t - 2) * (t - 3) * (t - 6), 36)
    raise ValueError(k)


def L_y(row: tuple[int, ...]) -> F:
    hist = missed_histogram(state_word(row))
    k = len(row)
    return sum(hist[t] * g_weight(k, t) for t in range(7))


def state_entropy(word: StateWord) -> float:
    c = Counter(dict(state_mass(word)))
    ent = 0.0
    for mass in c.values():
        p = float(mass)
        if p:
            ent -= p * log2(p)
    return ent


def transition_signature(word: StateWord) -> tuple[tuple[int, int], ...]:
    def hamming(a: SectorState, b: SectorState) -> int:
        return len(set(a).symmetric_difference(b))

    c: Counter[int] = Counter()
    states = [state for _mass, state in word]
    for a, b in zip(states, states[1:] + states[:1]):
        if a != b:
            c[hamming(a, b)] += 1
    return tuple(sorted(c.items()))


def sector_features(row: tuple[int, ...]) -> dict[str, object]:
    features = speed_features(tuple(x for x in row if x != 0))
    word = state_word(row)
    hist = missed_histogram(word)
    return {
        "sumset_excess": features["sumset_excess"],
        "fold_count": features["fold_count"],
        "fold_profile": features["fold_profile"],
        "missed_histogram": hist,
        "state_mass": state_mass(word),
        "transition_signature": transition_signature(word),
        "full_state_word": word,
    }


def sector_mixed_report(rows: list[tuple[int, ...]], feature_name: str, feature_fn) -> tuple[int, int, str]:
    fibers: dict[object, list[tuple[int, ...]]] = defaultdict(list)
    for row in rows:
        fibers[feature_fn(row)].append(row)
    mixed = [rs for rs in fibers.values() if len({L_y(r) for r in rs}) > 1]
    if not mixed:
        return 0, max((len(v) for v in fibers.values()), default=0), "no mixed L_y fibers"
    rs = max(mixed, key=len)
    by_ly: dict[F, tuple[int, ...]] = {}
    for row in rs:
        by_ly.setdefault(L_y(row), row)
        if len(by_ly) >= 2:
            break
    ex = " ; ".join(f"E={row}, L_y={ly}" for ly, row in by_ly.items())
    return len(mixed), len(rs), ex


def run_sector_lab() -> None:
    print()
    print("=" * 88)
    print("B. LRC14 sector state-word separation lab")
    print("=" * 88)
    rows = [tuple((0,) + c) for c in combinations(range(1, 14), 8) if set_gcd(tuple(c)) == 1]
    print(f"rows: {len(rows)} primitive k=9 offset rows E={{0}}+8-subsets of [1,13]")
    print("feature | mixed L_y fibers | largest mixed fiber | example")
    feature_table = [
        ("sumset_excess", lambda r: sector_features(r)["sumset_excess"]),
        ("fold_count", lambda r: sector_features(r)["fold_count"]),
        ("fold_profile", lambda r: sector_features(r)["fold_profile"]),
        ("transition_signature", lambda r: sector_features(r)["transition_signature"]),
        ("state_mass", lambda r: sector_features(r)["state_mass"]),
        ("missed_histogram", lambda r: sector_features(r)["missed_histogram"]),
        ("full_state_word", lambda r: sector_features(r)["full_state_word"]),
    ]
    for name, fn in feature_table:
        mixed_count, largest, example = sector_mixed_report(rows, name, fn)
        print(f"{name:24s} | {mixed_count:15d} | {largest:19d} | {example}")

    top = sorted(rows, key=lambda r: L_y(r), reverse=True)[:8]
    print()
    print("top k=9 rows by L_y with state metrics")
    print("row | L_y | p_hist | states | entropy | transitions | fold_count")
    for row in top:
        word = state_word(row)
        print(
            f"{row} | {L_y(row)} | {missed_histogram(word)} | {len(state_mass(word))} | "
            f"{state_entropy(word):.5f} | {transition_signature(word)} | "
            f"{sector_features(row)['fold_count']}"
        )


def run_synthesis() -> None:
    print()
    print("=" * 88)
    print("C. Invariant reading")
    print("=" * 88)
    print("1. Additive summaries are routing labels: sumset excess, fold count,")
    print("   fold profile, and gap pattern all have mixed fibers for exact M or L_y.")
    print("2. For max-min M, the local determinant is the optimal-time clearance word:")
    print("   denominator q, folded residues, active runners, and crossing source.")
    print("   It is a finite code on the THM-524 envelope, not a raw speed invariant.")
    print("3. For LRC14 sector measures, the missed-count histogram determines L_y,")
    print("   but HYP-2648's measured state word determines the richer structure:")
    print("   transition complexity, sector bias, and signed wall transport.")
    print("4. The common object is an addressed wall/crossing sheaf: atoms plus labels.")
    print("   Scalar invariants become trustworthy only after the address is retained.")
    print()
    print("Tournament Analysis")
    print("vertices: addressed_wall_sheaf, optimal_clearance_word, measured_state_word,")
    print("missed_histogram, additive_freiman_labels, fold_profile, raw_speeds")
    print("edge rule: A -> B when B is a quotient of A that can mix target fibers")
    print("Hamiltonian path:")
    print("  addressed_wall_sheaf > optimal_clearance_word > measured_state_word >")
    print("  missed_histogram > additive_freiman_labels > fold_profile > raw_speeds")
    print("challenged assumption: runners and arcs are not the determinant; the determinant")
    print("is the finite address carried by the wall arrangement before scalarization.")


def main() -> None:
    print("HYP-2650 / T897 LRC invariant separation scout")
    run_maxmin_lab()
    run_sector_lab()
    run_synthesis()


if __name__ == "__main__":
    main()
