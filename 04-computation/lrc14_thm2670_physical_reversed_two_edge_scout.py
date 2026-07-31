#!/usr/bin/env python3
"""Exact scout for THM-2670 equation (35).

For one reversed clock triple ``a <- b <- c`` this script materializes every
positive sharp-graph atom

    E^(s0)_(b,a;j,h)  and  E^(s1)_(c,b;h,k)

before rail/source/state labels are forgotten, and tests the physical overlap

    E^(s0)_(b,a;j,h) intersect D^(-1)E^(s1)_(c,b;h,k),
    D(x)={13x}.

The high-clock words are not expanded into 13^6 copies.  With z={13^6 x},
their product is the single finite word Q0(z)Q1(13z).  All remaining carrier
pieces are put on the denominator 13*T, where pulling the second edge back by
D requires only thirteen copies.  Arithmetic is exact integer interval
arithmetic throughout.

This is deliberately a scout rather than canon: it currently audits a chosen
clock triple and sector pair.  It prints the complete 12-by-12 source-pair
census of formal and physical state triples and retains a labelled atom-pair
witness for every positive source/state key.
"""

from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
import argparse
import sys


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import lrc14_cross_time_target_future_diagonal_thm2616 as cross
import lrc14_fallback_rail_digit_diagonal_pullback_thm2592 as old
import lrc14_guard_cospan_successor_private_clock_collapse as cospan


P = 13
Q7 = 7
INV2 = 7
PREDECESSOR_SPEED = P**5
SUCCESSOR_SPEED = P**6
T = old.T
TD = P * T
SECTORS = ("safe", "danger", "guard_free")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


@dataclass(frozen=True)
class Atom:
    sector: str
    rail_index: int
    source: int
    rail_clock: int
    shallow_clock: int
    rail_digit: int
    j: int
    h: int
    epsilon: int
    kappa: int
    old_numerator: int
    pieces: tuple


def prefix_intervals(prefix):
    starts, lens, _ = prefix
    return tuple((left, left + length) for left, length in zip(starts, lens))


def make_prefix(intervals):
    starts = [left for left, _ in intervals]
    lens = [right - left for left, right in intervals]
    pref = [0]
    for length in lens:
        pref.append(pref[-1] + length)
    return starts, lens, pref


def phi_at(x, starts, lens, pref):
    from bisect import bisect_right

    i = bisect_right(starts, x) - 1
    if i < 0:
        return 0
    return pref[i] + min(x - starts[i], lens[i])


def delayed_weighted_numerator(pieces, prefix, denominator=TD):
    """Numerator over SUCCESSOR_SPEED*denominator of int w(x)1_Q(Rx)."""
    starts, lens, pref = prefix
    length_q = pref[-1]
    weighted_len = 0
    acc_r = 0
    acc_phi = 0
    rred = SUCCESSOR_SPEED % denominator
    for left, right, weight in pieces:
        rleft = left * rred % denominator
        rright = right * rred % denominator
        weighted_len += weight * (right - left)
        acc_r += weight * (rright - rleft)
        acc_phi += weight * (
            phi_at(rright, starts, lens, pref)
            - phi_at(rleft, starts, lens, pref)
        )
    floor_numerator = SUCCESSOR_SPEED * weighted_len - acc_r
    require(
        floor_numerator % denominator == 0,
        "generalized weighted floor count is not integral",
    )
    result = length_q * (floor_numerator // denominator) + acc_phi
    require(result >= 0, "negative generalized delayed overlap")
    return result


def lift_profile(pieces):
    return tuple((P * a, P * b, w) for a, b, w in pieces)


def pullback_profile(pieces):
    """Represent x -> profile({13x}) on the denominator 13*T."""
    return tuple(
        (a + sheet * T, b + sheet * T, w)
        for sheet in range(P)
        for a, b, w in pieces
    )


def pullback_arrival_profile(pieces):
    arrival_left = 6 * T // P
    arrival_right = 7 * T // P
    require(
        all(arrival_left <= left < right <= arrival_right
            for left, right, _ in pieces),
        "weighted base event escaped arrival digit 6",
    )
    return tuple((6 * T + a, 6 * T + b, w) for a, b, w in pieces)


def lift_intervals(intervals):
    return tuple((P * a, P * b) for a, b in intervals)


def pullback_intervals(intervals):
    """Represent z -> 1_Q({13z}) on the denominator 13*T."""
    return tuple(
        (a + sheet * T, b + sheet * T)
        for sheet in range(P)
        for a, b in intervals
    )


def pullback_arrival_intervals(intervals):
    """Exact D-pullback on the only sheet meeting the retained arrival cell."""
    arrival_left = 6 * T // P
    arrival_right = 7 * T // P
    require(
        all(arrival_left <= left < right <= arrival_right
            for left, right in intervals),
        "base event escaped the inherited arrival digit 6",
    )
    # The first lifted event lies in [6T,7T) on the 13T grid.  A pulled-back
    # second arrival event can meet it only in sheet 6.
    return tuple((6 * T + left, 6 * T + right) for left, right in intervals)


def intersect_intervals(left, right):
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        a = max(left[i][0], right[j][0])
        b = min(left[i][1], right[j][1])
        if a < b:
            if out and out[-1][1] == a:
                out[-1] = (out[-1][0], b)
            else:
                out.append((a, b))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(out)


def merge_intervals(intervals):
    ordered = sorted(intervals)
    if not ordered:
        return ()
    out = [list(ordered[0])]
    for left, right in ordered[1:]:
        if left <= out[-1][1]:
            out[-1][1] = max(out[-1][1], right)
        else:
            out.append([left, right])
    return tuple((left, right) for left, right in out)


def intersect_profiles(left, right):
    """Multiply two sorted nonnegative step profiles."""
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        a = max(left[i][0], right[j][0])
        b = min(left[i][1], right[j][1])
        if a < b:
            weight = left[i][2] * right[j][2]
            if weight:
                if out and out[-1][1] == a and out[-1][2] == weight:
                    out[-1] = (out[-1][0], b, weight)
                else:
                    out.append((a, b, weight))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(out)


def build_sector_prefixes(module, inherited_safe):
    prefixes, _, _, _ = cospan.build_guard_cospan(module, inherited_safe)
    return prefixes


def build_edge_atoms(
    module,
    rails,
    present,
    starts,
    sector_prefixes,
    sector,
    rail_clock,
    shallow_clock,
):
    require(sector in SECTORS, "unknown delayed sector")
    d = (shallow_clock - rail_clock) % Q7
    require(d != 0, "edge clocks must be distinct")
    atoms = []
    candidates = 0
    for rail_index, (source, ell, rail_digit, pieces) in enumerate(rails):
        if ell != rail_clock:
            continue
        for h in range(12):
            r = (-h - 1) % P
            require(r != 0, "sharp graph entered the absent r=0 sheet")
            overlap = old.intersect_weighted_union(
                pieces,
                present[shallow_clock, (-h) % P],
                starts[shallow_clock, (-h) % P],
            )
            for epsilon, left, right in (
                (0, 14 * r, 14 * r + 13),
                (1, 14 * r - 13, 14 * r),
            ):
                half_tooth = old.intersect_weighted_comb(
                    overlap, module.C3, 182, left, right
                )
                for kappa in (0, 1):
                    candidates += 1
                    j = INV2 * (r - epsilon - kappa) % P
                    require(
                        h == (-2 * j - kappa - epsilon - 1) % P,
                        "predecessor/successor state equation failed",
                    )
                    predecessor = old.intersect_weighted_comb(
                        half_tooth,
                        PREDECESSOR_SPEED,
                        P,
                        j,
                        j + 1,
                    )
                    carry = old.intersect_weighted_comb(
                        predecessor,
                        SUCCESSOR_SPEED,
                        2,
                        kappa,
                        kappa + 1,
                    )
                    if not carry:
                        continue
                    numerator = old.delayed_weighted_numerator(
                        carry, sector_prefixes[sector][shallow_clock][h]
                    )
                    if not numerator:
                        continue
                    atoms.append(
                        Atom(
                            sector,
                            rail_index,
                            source,
                            rail_clock,
                            shallow_clock,
                            rail_digit,
                            j,
                            h,
                            epsilon,
                            kappa,
                            numerator,
                            tuple(carry),
                        )
                    )
    return tuple(atoms), candidates


def combined_word_prefix(prefix0, prefix1):
    q0 = lift_intervals(prefix_intervals(prefix0))
    q1 = pullback_intervals(prefix_intervals(prefix1))
    return make_prefix(intersect_intervals(q0, q1))


def single_edge_normalization_controls(atoms, prefixes, which, limit=12):
    """Check both denominator lifts reproduce Haar invariance exactly."""
    checked = 0
    whole_prefix = make_prefix(((0, TD),))
    for atom in atoms:
        if checked >= limit:
            break
        prefix = prefixes[atom.sector][atom.shallow_clock][atom.h]
        if which == "first":
            carrier = lift_profile(atom.pieces)
            word = make_prefix(lift_intervals(prefix_intervals(prefix)))
        else:
            carrier = pullback_profile(atom.pieces)
            word = make_prefix(pullback_intervals(prefix_intervals(prefix)))
        lifted = delayed_weighted_numerator(carrier, word)
        require(lifted == P * atom.old_numerator,
                f"{which}-edge denominator/Haar normalization failed")
        # A harmless whole-word control exercises the same generalized sweep.
        require(delayed_weighted_numerator(carrier, whole_prefix) >= lifted,
                f"{which}-edge whole-word domination failed")
        checked += 1
    require(checked == min(limit, len(atoms)), "normalization control shortfall")
    return checked


def source_matrix(counter):
    return tuple(
        tuple(counter.get((s0, s1), 0) for s1 in range(1, P))
        for s0 in range(1, P)
    )


def build_global_atom_bank(
    module, rails, present, starts, prefixes, sector
):
    bank = defaultdict(list)
    candidates = 0
    for rail in range(Q7):
        for shallow in range(Q7):
            if rail == shallow:
                continue
            atoms, count = build_edge_atoms(
                module, rails, present, starts, prefixes,
                sector, rail, shallow,
            )
            candidates += count
            for atom in atoms:
                key = (
                    atom.shallow_clock,
                    atom.rail_clock,
                    atom.source,
                    atom.j,
                    atom.h,
                )
                bank[key].append(atom)
    return {key: tuple(value) for key, value in bank.items()}, candidates


def global_physical_census(
    sector0, sector1, verbose=True, audit_atom_pairs=False
):
    """Exhaust all distinct labelled-state chains for one sector ordering."""
    (
        module,
        inherited_safe,
        _,
        _,
        rails,
        present,
        starts,
    ) = cross.build_carrier_data()
    prefixes = build_sector_prefixes(module, inherited_safe)
    bank0, candidates0 = build_global_atom_bank(
        module, rails, present, starts, prefixes, sector0
    )
    bank1 = bank0
    candidates1 = candidates0
    if sector1 != sector0:
        bank1, candidates1 = build_global_atom_bank(
            module, rails, present, starts, prefixes, sector1
        )

    support0 = {
        key: merge_intervals(
            (left, right)
            for atom in atoms
            for left, right, weight in atom.pieces
            if weight
        )
        for key, atoms in bank0.items()
    }
    support1 = support0 if bank1 is bank0 else {
        key: merge_intervals(
            (left, right)
            for atom in atoms
            for left, right, weight in atom.pieces
            if weight
        )
        for key, atoms in bank1.items()
    }
    lifted = {key: lift_intervals(value) for key, value in support0.items()}
    pulled = {
        key: pullback_arrival_intervals(value)
        for key, value in support1.items()
    }

    second_by_middle = defaultdict(list)
    for key1 in bank1:
        shallow, rail, source, predecessor, output = key1
        second_by_middle[shallow, predecessor].append(key1)

    totals = Counter()
    by_clock = Counter()
    by_clock_base = Counter()
    by_clock_physical = Counter()
    by_source_physical = Counter()
    by_clock_source_physical = Counter()
    physical_state_keys = set()
    word_cache = {}
    first_physical = None
    for key0 in bank0:
        a, b, source0, j, h = key0
        for key1 in second_by_middle.get((b, h), ()):
            b1, c, source1, h1, k = key1
            require(b1 == b and h1 == h, "global middle label mismatch")
            totals["formal"] += 1
            by_clock[a, b, c] += 1
            base = intersect_intervals(lifted[key0], pulled[key1])
            if not base:
                totals["base_zero"] += 1
                continue
            totals["base_positive"] += 1
            by_clock_base[a, b, c] += 1
            positive_atom_profiles = None
            if audit_atom_pairs:
                positive_atom_profiles = []
                for atom0 in bank0[key0]:
                    for atom1 in bank1[key1]:
                        atom_base = intersect_profiles(
                            lift_profile(atom0.pieces),
                            pullback_arrival_profile(atom1.pieces),
                        )
                        if atom_base:
                            positive_atom_profiles.append(atom_base)
                require(positive_atom_profiles,
                        "positive union event has no positive atom pair")
                totals["base_atom_pairs"] += len(positive_atom_profiles)
            word_key = (a, h, b, k)
            prefix = word_cache.get(word_key)
            if prefix is None:
                prefix = combined_word_prefix(
                    prefixes[sector0][a][h], prefixes[sector1][b][k]
                )
                word_cache[word_key] = prefix
            if prefix[2][-1] == 0:
                totals["word_empty"] += 1
                continue
            if audit_atom_pairs:
                totals["word_compatible_atom_pairs"] += len(
                    positive_atom_profiles
                )
            value = delayed_weighted_numerator(
                tuple((left, right, 1) for left, right in base), prefix
            )
            if not value:
                totals["skew_zero"] += 1
                continue
            totals["physical"] += 1
            physical_state_keys.add(
                (a, b, c, source0, source1, j, h, k)
            )
            by_clock_physical[a, b, c] += 1
            by_source_physical[source0, source1] += 1
            by_clock_source_physical[(a, b, c), (source0, source1)] += 1
            if audit_atom_pairs:
                atom_pairs = sum(
                    bool(delayed_weighted_numerator(atom_base, prefix))
                    for atom_base in positive_atom_profiles
                )
                require(atom_pairs > 0,
                        "positive union event has no positive atom pair")
                totals["physical_atom_pairs"] += atom_pairs
                totals["atom_pair_skew_zero"] += (
                    len(positive_atom_profiles) - atom_pairs
                )
            if first_physical is None:
                first_physical = (
                    (a, b, c, source0, source1, j, h, k), value
                )

    require(
        totals["formal"]
        == totals["base_zero"] + totals["word_empty"]
        + totals["skew_zero"] + totals["physical"],
        "global chain partition failed",
    )
    result = {
        "sector_pair": (sector0, sector1),
        "atom_counts": (sum(map(len, bank0.values())),
                        sum(map(len, bank1.values()))),
        "candidate_atom_counts": (candidates0, candidates1),
        "totals": dict(totals),
        "formal_by_clock": dict(by_clock),
        "base_by_clock": dict(by_clock_base),
        "physical_by_clock": dict(by_clock_physical),
        "physical_by_source": dict(by_source_physical),
        "physical_by_clock_source": dict(by_clock_source_physical),
        "physical_state_keys": physical_state_keys,
        "first_physical": first_physical,
    }
    if verbose:
        print("THM-2670 global reversed two-edge census -- exact")
        print(f"sector ordering={sector0}/{sector1}")
        print(
            f"atoms={result['atom_counts'][0]}/{result['atom_counts'][1]} "
            f"candidate atoms={candidates0}/{candidates1}"
        )
        print(f"totals={dict(sorted(totals.items()))}")
        print(f"formal by clock={tuple(sorted(by_clock.items()))}")
        print(f"base-positive by clock={tuple(sorted(by_clock_base.items()))}")
        print(f"physical by clock={tuple(sorted(by_clock_physical.items()))}")
        print(
            "physical source-pair support="
            f"{len(by_source_physical)}/144; counts="
            f"{tuple(sorted(by_source_physical.items()))}"
        )
        print(f"first physical={first_physical}")
    return result


def print_matrix(label, matrix):
    print(label)
    for source, row in enumerate(matrix, 1):
        print(f"  s0={source:2d}: " + " ".join(f"{value:3d}" for value in row))


def run_probe(a, b, c, sector0, sector1, diagnose_relaxed=False):
    require(len({a, b}) == 2 and len({b, c}) == 2,
            "each consecutive clock pair must be distinct")
    (
        module,
        inherited_safe,
        _,
        _,
        rails,
        present,
        starts,
    ) = cross.build_carrier_data()
    prefixes = build_sector_prefixes(module, inherited_safe)

    first, first_candidates = build_edge_atoms(
        module, rails, present, starts, prefixes,
        sector0, b, a,
    )
    second, second_candidates = build_edge_atoms(
        module, rails, present, starts, prefixes,
        sector1, c, b,
    )
    require(first_candidates > 0 and second_candidates > 0,
            "edge atom universes unexpectedly empty")
    controls0 = single_edge_normalization_controls(first, prefixes, "first")
    controls1 = single_edge_normalization_controls(second, prefixes, "second")

    second_by_input = defaultdict(list)
    pulled_second = {}
    for atom in second:
        second_by_input[atom.j].append(atom)
        pulled_second[atom] = pullback_profile(atom.pieces)
    lifted_first = {atom: lift_profile(atom.pieces) for atom in first}

    word_cache = {}
    carrier_cache = {}
    formal_keys = set()
    physical_keys = set()
    witnesses = {}
    formal_occurrences = 0
    carrier_occurrences = 0
    nonempty_word_occurrences = 0
    delayed_tests = 0
    physical_occurrences = 0
    positive_numerators = []
    formal_by_source = Counter()
    physical_by_source = Counter()

    for atom0 in first:
        for atom1 in second_by_input[atom0.h]:
            formal_occurrences += 1
            key = (atom0.source, atom1.source, atom0.j, atom0.h, atom1.h)
            formal_keys.add(key)
            formal_by_source[atom0.source, atom1.source] += 1

            carrier_key = (atom0, atom1)
            carrier = carrier_cache.get(carrier_key)
            if carrier is None:
                carrier = intersect_profiles(
                    lifted_first[atom0], pulled_second[atom1]
                )
                carrier_cache[carrier_key] = carrier
            if not carrier:
                continue
            carrier_occurrences += 1

            word_key = (atom0.h, atom1.h)
            prefix = word_cache.get(word_key)
            if prefix is None:
                prefix = combined_word_prefix(
                    prefixes[sector0][a][atom0.h],
                    prefixes[sector1][b][atom1.h],
                )
                word_cache[word_key] = prefix
            if prefix[2][-1] == 0:
                continue
            nonempty_word_occurrences += 1
            delayed_tests += 1
            numerator = delayed_weighted_numerator(carrier, prefix)
            if not numerator:
                continue

            physical_occurrences += 1
            positive_numerators.append(numerator)
            if key not in physical_keys:
                physical_keys.add(key)
                physical_by_source[atom0.source, atom1.source] += 1
                witnesses[key] = (
                    atom0.rail_index,
                    atom0.rail_digit,
                    atom0.epsilon,
                    atom0.kappa,
                    atom1.rail_index,
                    atom1.rail_digit,
                    atom1.epsilon,
                    atom1.kappa,
                    numerator,
                )

    print("THM-2670 physical reversed two-edge scout -- exact")
    print(f"clocks a<-b<-c = {a}<-{b}<-{c}")
    print(f"sectors first/second = {sector0}/{sector1}")
    print(f"T={T}  TD=13*T={TD}  R=13^6={SUCCESSOR_SPEED}")
    print(
        f"edge atoms first={len(first)}/{first_candidates} "
        f"second={len(second)}/{second_candidates}"
    )
    print(
        "carrier piece lengths first(min/max/sum)="
        f"{min(map(lambda x: len(x.pieces), first))}/"
        f"{max(map(lambda x: len(x.pieces), first))}/"
        f"{sum(len(x.pieces) for x in first)} "
        "second(min/max/sum)="
        f"{min(map(lambda x: len(x.pieces), second))}/"
        f"{max(map(lambda x: len(x.pieces), second))}/"
        f"{sum(len(x.pieces) for x in second)}"
    )
    print(f"normalization controls first={controls0} second={controls1}")
    print(
        f"formal occurrence pairs={formal_occurrences} "
        f"formal source/state keys={len(formal_keys)}"
    )
    print(
        f"physical positive occurrence pairs={physical_occurrences} "
        f"physical source/state keys={len(physical_keys)}"
    )
    print(
        f"nonempty carrier pairs={carrier_occurrences} "
        f"nonempty paired-word pairs={nonempty_word_occurrences} "
        f"delayed tests={delayed_tests}"
    )
    print(f"lost formal source/state keys={len(formal_keys - physical_keys)}")
    if positive_numerators:
        print(
            f"positive numerator min/max={min(positive_numerators)}/"
            f"{max(positive_numerators)} over denominator R*13*T"
        )
    print_matrix("formal occurrence-pair counts by (s0,s1):", source_matrix(formal_by_source))
    print_matrix("physical state-triple counts by (s0,s1):", source_matrix(physical_by_source))
    print("first ten labelled physical witnesses:")
    for key in sorted(witnesses)[:10]:
        print(f"  key(s0,s1,j,h,k)={key} witness={witnesses[key]}")

    if diagnose_relaxed:
        relaxed = Counter()
        relaxed_pairs = 0
        example = {}
        for atom0 in first:
            for atom1 in second:
                if intersect_profiles(lifted_first[atom0], pulled_second[atom1]):
                    delta = (atom1.j - atom0.h) % P
                    relaxed[delta] += 1
                    relaxed_pairs += 1
                    example.setdefault(delta, (atom0, atom1))
        print(f"relaxed full-carrier pairs (state equality dropped)={relaxed_pairs}")
        print(f"relaxed predecessor mismatch histogram={dict(sorted(relaxed.items()))}")
        for delta in sorted(example)[:3]:
            atom0, atom1 = example[delta]
            print(
                "  relaxed example delta="
                f"{delta}: first=(rail={atom0.rail_index},s={atom0.source},"
                f"j={atom0.j},h={atom0.h},eps={atom0.epsilon},kap={atom0.kappa}) "
                f"second=(rail={atom1.rail_index},s={atom1.source},"
                f"j={atom1.j},h={atom1.h},eps={atom1.epsilon},kap={atom1.kappa})"
            )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--a", type=int, default=2)
    parser.add_argument("--b", type=int, default=1)
    parser.add_argument("--c", type=int, default=0)
    parser.add_argument("--sector0", choices=SECTORS, default="safe")
    parser.add_argument("--sector1", choices=SECTORS, default="safe")
    parser.add_argument("--diagnose-relaxed", action="store_true")
    parser.add_argument("--global-census", action="store_true")
    args = parser.parse_args()
    require(all(0 <= clock < Q7 for clock in (args.a, args.b, args.c)),
            "clock labels must lie in F_7")
    if args.global_census:
        global_physical_census(args.sector0, args.sector1)
    else:
        run_probe(
            args.a, args.b, args.c, args.sector0, args.sector1,
            diagnose_relaxed=args.diagnose_relaxed,
        )


if __name__ == "__main__":
    main()
