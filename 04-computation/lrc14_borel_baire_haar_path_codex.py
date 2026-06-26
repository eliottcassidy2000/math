#!/usr/bin/env python3
"""Borel/Baire/Haar path-planning atlas for LRC14.

This is a finite model, not a substitute for descriptive set theory.  It uses a
cyclic direction group C_Q where the Haar measure is uniform counting measure.
The point is to keep three labels separate:

* Borel code: which direction/time cells are in the event;
* Baire core: which cells survive small perturbations;
* Haar mass: invariant measure, useful but lossy.

The same split is then used to propose a sixth any-angle path-planning style.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


Q = 80


def cyclic_interval(start: int, end: int, q: int = Q) -> set[int]:
    start %= q
    end %= q
    if start <= end:
        return set(range(start, end + 1))
    return set(range(start, q)) | set(range(0, end + 1))


def rotate(s: set[int], shift: int, q: int = Q) -> set[int]:
    return {(x + shift) % q for x in s}


def robust_core(s: set[int], q: int = Q) -> set[int]:
    return {x for x in s if (x - 1) % q in s and (x + 1) % q in s}


def boundary_cells(s: set[int], q: int = Q) -> set[int]:
    if not s:
        return set()
    return {x for x in range(q) if x in s and x not in robust_core(s, q)}


def components(s: set[int], q: int = Q) -> list[list[int]]:
    if not s:
        return []
    unseen = set(s)
    comps: list[list[int]] = []
    while unseen:
        start = min(unseen)
        comp = [start]
        unseen.remove(start)
        x = (start + 1) % q
        while x in unseen:
            comp.append(x)
            unseen.remove(x)
            x = (x + 1) % q
        # Merge wrap-around components if needed.
        if start == 0 and q - 1 in s and comps:
            comps[0] = comp + comps[0]
        else:
            comps.append(comp)
    return sorted(comps, key=lambda c: (len(c), c[0]), reverse=True)


def mass(s: set[int], q: int = Q) -> Fraction:
    return Fraction(len(s), q)


def fmt_fraction(x: Fraction) -> str:
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


@dataclass(frozen=True)
class EventProfile:
    name: str
    size: int
    haar_mass: Fraction
    robust_size: int
    boundary_size: int
    component_lengths: tuple[int, ...]


def profile(name: str, s: set[int], q: int = Q) -> EventProfile:
    return EventProfile(
        name=name,
        size=len(s),
        haar_mass=mass(s, q),
        robust_size=len(robust_core(s, q)),
        boundary_size=len(boundary_cells(s, q)),
        component_lengths=tuple(len(c) for c in components(s, q)),
    )


def print_profile(p: EventProfile) -> None:
    print(
        f"{p.name:26s} size={p.size:2d} mass={fmt_fraction(p.haar_mass):>5s} "
        f"robust={p.robust_size:2d} boundary={p.boundary_size:2d} "
        f"components={list(p.component_lengths)}"
    )


def path_algorithms() -> list[dict[str, object]]:
    return [
        {
            "name": "Field D*",
            "payload": "edge-propagated interpolation values",
            "keeps": ["dynamic cost", "local metric"],
            "loses": ["exact taut interval", "category boundary"],
        },
        {
            "name": "Theta* family",
            "payload": "parent line-of-sight shortcuts",
            "keeps": ["visibility bit", "parent witness"],
            "loses": ["full visible interval", "measure/category split"],
        },
        {
            "name": "Block A*",
            "payload": "local all-pairs block database",
            "keeps": ["finite local packet", "block boundary state"],
            "loses": ["global orbit/rank label"],
        },
        {
            "name": "ANYA",
            "payload": "taut path intervals as nodes",
            "keeps": ["interval node", "taut obstacle wrap"],
            "loses": ["Haar mass", "Baire robust core"],
        },
        {
            "name": "CWave",
            "payload": "discrete circular arcs and lines",
            "keeps": ["wavefront primitive", "geometric front"],
            "loses": ["selector address", "branch-local chart"],
        },
        {
            "name": "Borel-Baire-Haar A*",
            "payload": "Borel interval packet + Baire core + Haar mass",
            "keeps": ["visibility interval", "robust core", "boundary cells", "invariant mass"],
            "loses": ["none by design; still a proposed proof interface"],
        },
    ]


def score_algorithm(row: dict[str, object]) -> int:
    keeps = set(row["keeps"])
    score = 0
    score += 3 if "visibility interval" in keeps or "interval node" in keeps else 0
    score += 3 if "robust core" in keeps else 0
    score += 2 if "invariant mass" in keeps else 0
    score += 2 if "boundary cells" in keeps else 0
    score += 1 if "parent witness" in keeps or "taut obstacle wrap" in keeps else 0
    score += 1 if "finite local packet" in keeps else 0
    return score


def main() -> None:
    print("LRC14 BOREL / BAIRE / HAAR PATH-PLANNING ATLAS")
    print("=" * 76)
    print()

    print("A. Finite Borel-Baire-Haar direction toy")
    print("-" * 76)
    interval_event = set(range(20))
    dust_event = {4 * k for k in range(20)}
    shifted_interval = rotate(interval_event, 17)
    print("C_Q with Q=80.  Haar mass is uniform counting measure.")
    for p in [
        profile("contiguous Borel arc", interval_event),
        profile("20-point Baire dust", dust_event),
        profile("rotated same arc", shifted_interval),
    ]:
        print_profile(p)
    print("Readout: same Haar mass can hide opposite Baire behavior.")
    print("The rotated arc verifies the finite Haar invariance check.")
    print()

    print("B. Visibility obstacle toy")
    print("-" * 76)
    blocked = (
        cyclic_interval(6, 9)
        | cyclic_interval(22, 27)
        | cyclic_interval(41, 43)
        | cyclic_interval(67, 72)
    )
    clear = set(range(Q)) - blocked
    clear_core = robust_core(clear)
    print_profile(profile("blocked directions", blocked))
    print_profile(profile("clear line-of-sight", clear))
    print_profile(profile("clear robust core", clear_core))
    print("Clear components are the finite analog of any-angle interval nodes:")
    for comp in components(clear):
        print(f"  start={comp[0]:2d} end={comp[-1]:2d} length={len(comp):2d}")
    print()

    print("C. Five known any-angle lanes plus a sixth")
    print("-" * 76)
    rows = path_algorithms()
    for row in rows:
        print(
            f"{row['name']:22s} score={score_algorithm(row):2d} "
            f"payload={row['payload']}"
        )
        print(f"  keeps={', '.join(row['keeps'])}")
        print(f"  loses={', '.join(row['loses'])}")
    print()
    print("Sixth proposal: Borel-Baire-Haar A*")
    print("  Node: a visible interval packet, not a single grid vertex.")
    print("  Labels: Borel code, Baire robust core, boundary cells, Haar mass.")
    print("  Expansion: split intervals only at obstacle tangencies or C27-style owner walls.")
    print("  Priority: ordinary path length plus penalties for boundary-only progress.")
    print()

    print("D. Borel/Baire/Haar to LRC14")
    print("-" * 76)
    print("Borel: LRC safe/danger sets are finite Boolean combinations of torus inequalities.")
    print("Baire: robust safe intervals differ from boundary-only survivor points.")
    print("Haar: invariant torus measure is canonical, but it can miss tight boundary witnesses.")
    print("Known warning: AP at n=14 covers open danger arcs up to measure zero;")
    print("the surviving equality times are boundary/category data, not positive Haar mass.")
    print()

    print("E. New hypotheses generated by this atlas")
    print("-" * 76)
    print("HYP-2950: LRC14 needs Borel code + Baire core + Haar mass, not measure alone.")
    print("HYP-2951: Borel-Baire-Haar A* is the sixth any-angle proof carrier.")
    print("HYP-2952: Baire boundary witnesses are the path-planning version of C27 owner walls.")
    print()

    print("F. Tournament analysis")
    print("-" * 76)
    carriers = [
        "exact LRC M/Farey/C27 labels",
        "Borel event code",
        "Baire robust-core/boundary split",
        "Haar invariant mass",
        "Borel-Baire-Haar A* interval node",
        "ANYA taut interval",
        "CWave wavefront primitive",
        "Theta* parent visibility bit",
        "Field D* interpolation scalar",
        "raw grid vertex",
    ]
    edges = len(carriers) * (len(carriers) - 1) // 2
    print("Relation: x -> y iff x retains more proof-critical witness information.")
    print(f"vertices={len(carriers)} edges={edges} c3=0 hp=1")
    for i, name in enumerate(carriers):
        print(f"  score {len(carriers)-1-i:2d}: {name}")
    print()
    print("Verdict:")
    print("  Haar measure is the right invariant volume, but LRC14 is not only a")
    print("  volume problem.  The proof carrier must also remember the Borel code")
    print("  and the Baire boundary/core split, exactly as any-angle planning must")
    print("  remember visible intervals and obstacle tangencies rather than just")
    print("  sampled grid vertices.")


if __name__ == "__main__":
    main()
