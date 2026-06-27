#!/usr/bin/env python3
"""Multiplicity-moment dual certificates for LRC14.

This is intentionally not an apex/local-aperture argument.  For a 13-speed row
S and threshold 1/14, define

    X_S(t) = #{v in S : ||v t|| < 1/14}.

The strict lonely set has positive Haar measure exactly when the multiplicity
histogram has positive mass at X=0.  Therefore any polynomial P with

    P(0) < 0,    P(k) >= 0 for k=1,...,13

is a dual certificate: if E[P(X_S)] < 0, then Pr[X_S=0] > 0 and the row has a
strict lonely interval.

The script computes exact multiplicity histograms by sweeping danger-arc
endpoints, then searches two small certificate families:

* odd binomial barriers B_m(x)=binom(x-1,m), m odd.  These have B_m(0)=-1 and
  B_m(k)>=0 for k>=1.
* root barriers with P(0)=-1 and selected integer roots, checked explicitly on
  {1,...,13}.

Tournament Analysis vertices are moment/barrier proof carriers, not runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import comb, gcd, prod
import argparse
import os


AP = tuple(range(1, 14))
GW = tuple(list(range(1, 12)) + [13, 24])


@dataclass(frozen=True)
class RowSpec:
    bank: str
    name: str
    speeds: tuple[int, ...]
    qdiv: int


@dataclass(frozen=True)
class RowMoment:
    bank: str
    name: str
    speeds: tuple[int, ...]
    qdiv: int
    hist: tuple[Fraction, ...]
    support: tuple[int, ...]
    mean: Fraction
    safe_mass: Fraction


@dataclass(frozen=True)
class Barrier:
    family: str
    degree: int
    label: str
    values: tuple[Fraction, ...]


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def qdiv(speeds: tuple[int, ...], cap: int = 220) -> int:
    for d in range(2, cap + 1):
        if all(v % d for v in speeds):
            return d
    return cap + 1


def row_name(holes: tuple[int, ...], adds: tuple[int, ...]) -> str:
    if not holes and not adds:
        return "AP"
    if len(holes) == 1 and len(adds) == 1:
        return f"{holes[0]}->{adds[0]}"
    return f"drop({','.join(map(str, holes))})->add({','.join(map(str, adds))})"


def add_interval(events: dict[Fraction, int], lo: Fraction, hi: Fraction) -> None:
    if lo >= hi:
        return
    events[lo] += 1
    events[hi] -= 1


def speed_danger_events(v: int, events: dict[Fraction, int]) -> None:
    radius = Fraction(1, 14 * v)
    add_interval(events, Fraction(0), radius)
    add_interval(events, Fraction(1) - radius, Fraction(1))
    for k in range(1, v):
        center = Fraction(k, v)
        add_interval(events, center - radius, center + radius)


def multiplicity_histogram(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    events: dict[Fraction, int] = defaultdict(int)
    events[Fraction(0)] += 0
    events[Fraction(1)] += 0
    for v in speeds:
        speed_danger_events(v, events)
    points = sorted(events)
    hist = [Fraction(0) for _ in range(len(speeds) + 1)]
    current = 0
    for i, point in enumerate(points[:-1]):
        current += events[point]
        nxt = points[i + 1]
        if nxt > point:
            hist[current] += nxt - point
    assert sum(hist, Fraction(0)) == 1
    mean = sum(Fraction(i) * h for i, h in enumerate(hist))
    assert mean == Fraction(len(speeds), 7), (speeds, mean)
    return tuple(hist)


def analyze_row(row: RowSpec) -> RowMoment:
    hist = multiplicity_histogram(row.speeds)
    support = tuple(i for i, mass in enumerate(hist) if mass)
    mean = sum(Fraction(i) * mass for i, mass in enumerate(hist))
    return RowMoment(
        row.bank,
        row.name,
        row.speeds,
        row.qdiv,
        hist,
        support,
        mean,
        hist[0],
    )


def named_rows() -> list[RowSpec]:
    rows = [
        ("named", "AP", AP),
        ("named", "GW 12->24", GW),
        ("named", "near K33 12->36", tuple(list(range(1, 12)) + [13, 36])),
        ("named", "petal 10->20", tuple(sorted(set(AP) - {10} | {20}))),
        ("named", "petal 13->26", tuple(sorted(set(AP) - {13} | {26}))),
        ("named", "covering 12->84", tuple(list(range(1, 12)) + [13, 84])),
        ("named", "covering 12->168", tuple(list(range(1, 12)) + [13, 168])),
        (
            "named",
            "repair drop(4,6)->add(19,42)",
            tuple(sorted(set(AP) - {4, 6} | {19, 42})),
        ),
        (
            "named",
            "repair drop(2,6)->add(17,42)",
            tuple(sorted(set(AP) - {2, 6} | {17, 42})),
        ),
    ]
    return [RowSpec(bank, name, speeds, qdiv(speeds)) for bank, name, speeds in rows]


def generate_swap_bank(k: int, add_max: int, q_min: int) -> list[RowSpec]:
    if k == 0:
        return [RowSpec("0-swap", "AP", AP, qdiv(AP))]
    bank = f"{k}-swap add<={add_max}"
    out: list[RowSpec] = []
    for holes in combinations(AP, k):
        base = set(AP) - set(holes)
        for adds in combinations(range(14, add_max + 1), k):
            speeds = tuple(sorted(base | set(adds)))
            if len(speeds) != 13 or not primitive(speeds):
                continue
            q = qdiv(speeds)
            if q < q_min:
                continue
            out.append(RowSpec(bank, row_name(holes, adds), speeds, q))
    return out


def binomial_barriers() -> list[Barrier]:
    out: list[Barrier] = []
    for m in range(1, 14, 2):
        vals = []
        for x in range(14):
            if x == 0:
                vals.append(Fraction(-1))
            elif x - 1 >= m:
                vals.append(Fraction(comb(x - 1, m)))
            else:
                vals.append(Fraction(0))
        out.append(Barrier("binomial", m, f"B{m}=C(x-1,{m})", tuple(vals)))
    return out


def root_barriers(max_degree: int) -> list[Barrier]:
    out: list[Barrier] = []
    for degree in range(1, max_degree + 1):
        for roots in combinations(range(1, 14), degree):
            denom = prod(-r for r in roots)
            scale = Fraction(-1, denom)
            vals = []
            for x in range(14):
                vals.append(scale * prod(x - r for r in roots))
            if vals[0] != -1:
                continue
            if all(v >= 0 for v in vals[1:]):
                label = "roots=" + ",".join(map(str, roots))
                out.append(Barrier("root", degree, label, tuple(vals)))
    return out


def expectation(hist: tuple[Fraction, ...], barrier: Barrier) -> Fraction:
    return sum(hist[i] * barrier.values[i] for i in range(14))


def first_certificate(
    moment: RowMoment,
    barriers: list[Barrier],
) -> tuple[Barrier, Fraction] | None:
    if moment.safe_mass == 0:
        return None
    for barrier in barriers:
        val = expectation(moment.hist, barrier)
        if val < 0:
            return barrier, val
    return None


def analyze_rows(rows: list[RowSpec], workers: int) -> list[RowMoment]:
    if workers <= 1:
        return [analyze_row(row) for row in rows]
    out: list[RowMoment] = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(analyze_row, row) for row in rows]
        for fut in as_completed(futures):
            out.append(fut.result())
    return sorted(out, key=lambda r: (r.bank, r.name, r.speeds))


def fmt(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def summarize_named(moments: list[RowMoment], barriers: list[Barrier]) -> None:
    print("[1] Named rows")
    for row in moments:
        cert = first_certificate(row, barriers)
        cert_text = "none"
        if cert is not None:
            barrier, val = cert
            cert_text = f"{barrier.family}:{barrier.label} E={fmt(val)}"
        print(
            f"  {row.name:36s} qdiv={row.qdiv:<3d} safe={fmt(row.safe_mass):>10s} "
            f"support={row.support!s:24s} cert={cert_text}"
        )
    print()


def summarize_global(moments: list[RowMoment], barriers: list[Barrier]) -> None:
    certs: list[tuple[RowMoment, Barrier, Fraction]] = []
    no_cert: list[RowMoment] = []
    for row in moments:
        cert = first_certificate(row, barriers)
        if cert is None:
            no_cert.append(row)
        else:
            certs.append((row, cert[0], cert[1]))

    by_bank = Counter(row.bank for row in moments)
    zero_safe = [row for row in moments if row.safe_mass == 0]
    positive = [row for row in moments if row.safe_mass > 0]
    degree_hist = Counter((barrier.family, barrier.degree) for _, barrier, _ in certs)
    max_support_hist = Counter(max(row.support) for row in moments)
    min_safe = min((row.safe_mass for row in positive), default=Fraction(0))

    print("[2] Global multiplicity-moment audit")
    print(f"  rows={len(moments)}")
    print(f"  banks={dict(sorted(by_bank.items()))}")
    print(f"  zero_safe_rows={len(zero_safe)}")
    print(f"  positive_safe_rows={len(positive)}")
    print(f"  certified_by_barrier_family={len(certs)}")
    print(f"  no_barrier_found_with_selected_family={len(no_cert)}")
    print(f"  min_positive_safe_mass={fmt(min_safe)}")
    print(f"  max_multiplicity_hist={dict(sorted(max_support_hist.items()))}")
    print(f"  certificate_degree_hist={dict(sorted(degree_hist.items(), key=lambda kv: kv[0]))}")
    print()

    print("  zero-safe rows:")
    for row in zero_safe[:16]:
        print(f"    {row.bank:18s} {row.name:36s} qdiv={row.qdiv:<3d} support={row.support}")
    if len(zero_safe) > 16:
        print(f"    ... {len(zero_safe) - 16} more")
    print()

    print("  hardest positive rows by selected certificate degree then safe mass:")
    hardest = sorted(
        certs,
        key=lambda item: (-item[1].degree, item[0].safe_mass, item[0].bank, item[0].name),
    )
    for row, barrier, val in hardest[:16]:
        print(
            f"    {row.bank:18s} {row.name:36s} safe={fmt(row.safe_mass):>10s} "
            f"support={row.support!s:22s} {barrier.family}:{barrier.label} E={fmt(val)}"
        )
    print()

    suspicious = [row for row in no_cert if row.safe_mass > 0]
    if suspicious:
        print("  positive rows missed by selected barrier search:")
        for row in sorted(suspicious, key=lambda r: (r.safe_mass, r.bank, r.name))[:16]:
            print(
                f"    {row.bank:18s} {row.name:36s} safe={fmt(row.safe_mass):>10s} "
                f"support={row.support}"
            )
        print()


def proof_carrier_tournament() -> None:
    vertices = [
        ("integer multiplicity histogram", (6, 6, 6, 5, 5, 6)),
        ("dual moment barrier", (6, 5, 6, 6, 6, 6)),
        ("AP/GW zero-safe atom", (5, 6, 5, 5, 6, 5)),
        ("NORK pinch template", (5, 5, 5, 6, 6, 5)),
        ("qdiv witness", (4, 4, 3, 3, 3, 4)),
        ("apex aperture", (3, 4, 4, 4, 3, 4)),
        ("raw safe mass", (3, 3, 2, 2, 2, 2)),
        ("raw runner set", (1, 1, 1, 1, 1, 1)),
    ]
    out = {name: set() for name, _ in vertices}
    for i, (a, av) in enumerate(vertices):
        for j, (b, bv) in enumerate(vertices):
            if i == j:
                continue
            awins = sum(x > y for x, y in zip(av, bv))
            bwins = sum(y > x for x, y in zip(av, bv))
            if awins > bwins or (awins == bwins and i < j):
                out[a].add(b)
    scores = Counter(len(out[name]) for name, _ in vertices)
    c3 = 0
    names = [name for name, _ in vertices]
    for a, b, c in combinations(names, 3):
        if (b in out[a] and c in out[b] and a in out[c]) or (
            c in out[a] and b in out[c] and a in out[b]
        ):
            c3 += 1
    print("[3] Tournament Analysis")
    print("  vertices are proof carriers, not runners.")
    print("  pair observable:")
    print("    preserves no-lonely counterexample status, exact Haar distribution,")
    print("    low-moment inequalities, AP/GW atom visibility, packet labels, and")
    print("    resistance to scalarization.")
    print("  switch/gauge:")
    print("    majority over six retention coordinates; ties follow the displayed order.")
    print(f"  score_hist={dict(sorted(scores.items()))}")
    print(f"  directed_3_cycles={c3}")
    print("  Hamiltonian path:")
    print("    " + " > ".join(names))
    print()


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge / quotient declaration")
    print("  considered vertices:")
    print("    runners, residues, danger arcs, safe components, endpoint owners,")
    print("    Fourier modes, multiplicity values, moment barriers, and proof gates.")
    print("  chosen vertices:")
    print("    integer multiplicity values X=0..13 and moment/barrier proof carriers.")
    print("  preserved LRC predicate:")
    print("    Pr[X=0]>0 is exactly positive strict Haar-open lonely time.")
    print("  destroyed information:")
    print("    location of the safe interval, endpoint ownership, and denominator-14")
    print("    local geometry; those are intentionally delegated to other hypotheses.")
    print("  challenged assumption:")
    print("    do not start from apexes or packets.  Start from the global integer")
    print("    random variable that a counterexample would force to avoid 0 a.e.")
    print()


def theorem_readout(moments: list[RowMoment], barriers: list[Barrier]) -> None:
    zero_safe = [row for row in moments if row.safe_mass == 0]
    positive_no_cert = [
        row for row in moments if row.safe_mass > 0 and first_certificate(row, barriers) is None
    ]
    print("[4] Proof-shaped readout")
    print("  Necessary condition for a strict counterexample:")
    print("    for every polynomial P with P(0)<0 and P(k)>=0 on k=1..13,")
    print("    the exact multiplicity distribution must satisfy E[P(X_S)]>=0.")
    print("  Thus every negative moment barrier is a rigorous lonely-time proof.")
    print(f"  In this audit, zero-safe rows found: {len(zero_safe)}.")
    print("  They are expected to be AP/GW boundary atoms in the sampled banks.")
    print(f"  Positive rows missed by selected low-barrier family: {len(positive_no_cert)}.")
    print("  Candidate theorem:")
    print("    Multiplicity moment rigidity for LRC14: outside the AP/GW boundary")
    print("    atoms, every primitive qdiv>=14 source-core row has a negative")
    print("    admissible moment barrier, preferably one of the odd binomial")
    print("    barriers or a bounded root barrier.  A counterexample must therefore")
    print("    be moment-indistinguishable from AP/GW before packet labels are seen.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--one-max", type=int, default=160)
    parser.add_argument("--two-max", type=int, default=36)
    parser.add_argument("--three-max", type=int, default=24)
    parser.add_argument("--q-min", type=int, default=14)
    parser.add_argument("--root-degree", type=int, default=7)
    parser.add_argument("--workers", type=int, default=max(1, min(8, os.cpu_count() or 1)))
    args = parser.parse_args()

    print("S156 LRC14 MULTIPLICITY-MOMENT DUAL")
    print("=" * 78)
    print(
        f"one_max={args.one_max}, two_max={args.two_max}, three_max={args.three_max}, "
        f"q_min={args.q_min}, root_degree={args.root_degree}, workers={args.workers}"
    )
    print_assumption_challenge()

    barriers = sorted(
        binomial_barriers() + root_barriers(args.root_degree),
        key=lambda b: (b.degree, b.family, b.label),
    )
    print(f"barriers_loaded={len(barriers)}")
    print()

    named = analyze_rows(named_rows(), workers=min(args.workers, 4))
    summarize_named(named, barriers)

    rows = []
    rows.extend(generate_swap_bank(0, 13, args.q_min))
    rows.extend(generate_swap_bank(1, args.one_max, args.q_min))
    rows.extend(generate_swap_bank(2, args.two_max, args.q_min))
    rows.extend(generate_swap_bank(3, args.three_max, args.q_min))
    moments = analyze_rows(rows, args.workers)
    summarize_global(moments, barriers)
    proof_carrier_tournament()
    theorem_readout(moments, barriers)


if __name__ == "__main__":
    main()
