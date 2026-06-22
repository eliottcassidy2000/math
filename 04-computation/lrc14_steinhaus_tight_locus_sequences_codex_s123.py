#!/usr/bin/env python3
"""S123: exact Steinhaus/gap-sequence probes for the LRC14 tight locus.

This script targets the remaining three-gap rigidity/census obligation:

    primitive 13-speed tight rows M(S)=1/14 should be only AP and Goddyn-Wong.

It does not claim that theorem.  It records exact sequence invariants that make
the obstruction sharper:

* the denominator-14 cyclic gap word, with collisions retained as zero gaps;
* the residue defect sequence against the perfect AP clock;
* the apex-lock filter: all denominator-14 unit witnesses are blocked down to
  exactly 1/14, so no unit point immediately proves M(S)>1/14;
* the off-apex escape value, computed by the exact rational lower-envelope
  critical set.

The important separation is:

    apex-locked  !=  tight.

Many perturbations have AP/GW-like residue defects at denominator 14, but only
AP and GW survive the global exact-M check in the bounded banks below.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import sys


sys.stdout.reconfigure(line_buffering=True)


N = 14
THR = F(1, N)
UNITS14 = (1, 3, 5, 9, 11, 13)
AP = tuple(range(1, 14))
GW = tuple(list(range(1, 12)) + [13, 24])


def norm1(x: F) -> F:
    r = x % 1
    return min(r, 1 - r)


def primitive(S: tuple[int, ...]) -> bool:
    g = 0
    for s in S:
        g = gcd(g, s)
    return g == 1


def candidate_taus(S: tuple[int, ...]) -> set[F]:
    """All lower-envelope vertices in (0,1/2].

    The max of min_s ||s t|| occurs at a tent peak or a crossing of two tent
    sides, so these rational candidates are exact.  The function is even and
    1-periodic, so the half-interval is enough for M(S).
    """
    S = tuple(sorted(set(S)))
    cands: set[F] = set()
    for i, a in enumerate(S):
        k = 0
        while True:
            t = F(2 * k + 1, 2 * a)
            if t > F(1, 2):
                break
            cands.add(t)
            k += 1
        for b in S[i + 1 :]:
            for d in (a + b, b - a):
                if d <= 0:
                    continue
                k = 1
                while True:
                    t = F(k, d)
                    if t > F(1, 2):
                        break
                    cands.add(t)
                    k += 1
    cands.add(F(1, 2))
    return cands


_M_CACHE: dict[tuple[int, ...], tuple[F, tuple[F, ...]]] = {}


def M_exact(S: tuple[int, ...]) -> tuple[F, tuple[F, ...]]:
    S = tuple(sorted(set(S)))
    if S in _M_CACHE:
        return _M_CACHE[S]
    best = F(0)
    pts: list[F] = []
    for t in candidate_taus(S):
        if not (0 < t < 1):
            continue
        val = min(norm1(s * t) for s in S)
        if val > best:
            best = val
            pts = [t]
        elif val == best:
            pts.append(t)
    out = (best, tuple(sorted(pts)))
    _M_CACHE[S] = out
    return out


def binders(S: tuple[int, ...], t: F, value: F) -> tuple[int, ...]:
    return tuple(s for s in S if norm1(s * t) == value)


def residue_counts_with_observer(S: tuple[int, ...], a: int) -> tuple[int, ...]:
    counts = [0] * N
    counts[0] = 1  # observer
    for s in S:
        counts[(a * s) % N] += 1
    return tuple(counts)


def cyclic_gap_word_from_counts(counts: tuple[int, ...]) -> tuple[int, ...]:
    positions: list[int] = []
    for r, c in enumerate(counts):
        positions.extend([r] * c)
    positions.sort()
    gaps: list[int] = []
    for i, r in enumerate(positions):
        nxt = positions[(i + 1) % len(positions)]
        if i == len(positions) - 1:
            nxt += N
        gaps.append(nxt - r)
    return tuple(gaps)


def gap_partition(gaps: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(gaps).items()))


def residue_defect(counts: tuple[int, ...]) -> tuple[int, ...]:
    # AP clock has exactly one occupant at every residue: observer at 0, one
    # runner at each nonzero residue.
    return tuple(c - 1 for c in counts)


def steinhaus_energy(gaps: tuple[int, ...]) -> int:
    # Zero on AP's uniform 14-clock.  GW gives one zero gap and one two-gap,
    # hence energy 2.
    return sum((g - 1) * (g - 1) for g in gaps)


def apex_unit_profile(S: tuple[int, ...], a: int) -> dict[str, object]:
    counts = residue_counts_with_observer(S, a)
    gaps = cyclic_gap_word_from_counts(counts)
    defect = residue_defect(counts)
    speed_residues = tuple((a * s) % N for s in S)
    residue_distances = [min(r, N - r) for r in speed_residues]
    return {
        "a": a,
        "min_residue_distance": min(residue_distances),
        "counts": counts,
        "gap_word": gaps,
        "gap_partition": gap_partition(gaps),
        "distinct_gap_lengths": len(set(gaps)),
        "energy": steinhaus_energy(gaps),
        "missing": tuple(i for i, d in enumerate(defect) if d < 0),
        "extra": tuple((i, d) for i, d in enumerate(defect) if d > 0),
        "l1_defect": sum(abs(d) for d in defect),
    }


def apex_profiles(S: tuple[int, ...]) -> tuple[dict[str, object], ...]:
    return tuple(apex_unit_profile(S, a) for a in UNITS14)


def best_apex_min_distance(S: tuple[int, ...]) -> int:
    return max(p["min_residue_distance"] for p in apex_profiles(S))  # type: ignore[arg-type]


def is_apex_locked(S: tuple[int, ...]) -> bool:
    # If a unit has min distance >=2 residues, then t=a/14 proves strict slack.
    # Tight rows must have best apex distance exactly 1.  Rows with a multiple
    # of 14 have distance 0 and are not apex-boundary candidates.
    return best_apex_min_distance(S) == 1


def canonical_profile_key(S: tuple[int, ...]) -> tuple[tuple[tuple[int, int], int, int], ...]:
    """A unit-invariant sequence key for denominator-14 gap behavior."""
    parts = []
    for p in apex_profiles(S):
        parts.append(
            (
                p["gap_partition"],  # type: ignore[arg-type]
                p["energy"],  # type: ignore[arg-type]
                p["min_residue_distance"],  # type: ignore[arg-type]
            )
        )
    return tuple(sorted(parts))


def row_summary(name: str, S: tuple[int, ...]) -> None:
    M, pts = M_exact(S)
    print(f"{name}: S={S}")
    print(f"  M={M} ({float(M):.9f}), argmax_count={len(pts)}, denominators={sorted({t.denominator for t in pts})}")
    print(f"  argmaxes={pts[:12]}{' ...' if len(pts) > 12 else ''}")
    print(f"  apex_best_residue_distance={best_apex_min_distance(S)}, apex_locked={is_apex_locked(S)}")
    print("  denominator-14 unit gap sequences:")
    for p in apex_profiles(S):
        a = p["a"]
        t = F(a, N)
        b = binders(S, t, min(norm1(s * t) for s in S))
        print(
            f"    a={a:2d}: minResid={p['min_residue_distance']}, "
            f"gapPart={p['gap_partition']}, E={p['energy']}, "
            f"missing={p['missing']}, extra={p['extra']}, binders={b}"
        )
    print()


def scan_single_swaps(limit: int = 300) -> dict[str, object]:
    rows = []
    profile_counts: Counter[tuple[tuple[tuple[int, int], int, int], ...]] = Counter()
    tight = []
    below = []
    apex_locked_loose = []
    apex_slack_direct = []
    apex_blocked_exact = []
    for rem in AP:
        kept = [x for x in AP if x != rem]
        for v in range(1, limit + 1):
            S = tuple(sorted(set(kept + [v])))
            if len(S) != 13 or not primitive(S):
                continue
            best_dist = best_apex_min_distance(S)
            locked = best_dist == 1
            key = canonical_profile_key(S)
            profile_counts[key] += 1
            if best_dist > 1:
                # Some denominator-14 unit point is already strictly lonely.
                # Exact M is unnecessary for the tight/below classification.
                apex_slack_direct.append((rem, v, S, None, (), key, best_dist))
                rows.append((rem, v, S, None, (), key, best_dist))
                continue
            M, pts = M_exact(S)
            rec = (rem, v, S, M, pts, key, best_dist)
            rows.append(rec)
            if M == THR:
                tight.append(rec)
            elif M < THR:
                below.append(rec)
            elif locked:
                apex_locked_loose.append(rec)
            else:
                apex_blocked_exact.append(rec)
    return {
        "rows": rows,
        "profile_counts": profile_counts,
        "tight": tight,
        "below": below,
        "apex_locked_loose": apex_locked_loose,
        "apex_slack_direct": apex_slack_direct,
        "apex_blocked_exact": apex_blocked_exact,
    }


def scan_two_swaps(limit: int = 36) -> dict[str, object]:
    """A moderate exact AP two-swap bank.

    Replace two AP elements by two values in [1,limit].  This is not a proof
    search; it is a stress test for whether the sequence filters identify a
    wider tight family.
    """
    tight = []
    below = []
    apex_locked = []
    count = 0
    by_profile: Counter[tuple[tuple[tuple[int, int], int, int], ...]] = Counter()
    direct_slack = 0
    blocked_exact = 0
    values = range(1, limit + 1)
    for rems in combinations(AP, 2):
        kept = set(AP) - set(rems)
        for adds in combinations(values, 2):
            if kept & set(adds):
                continue
            S = tuple(sorted(kept | set(adds)))
            if len(S) != 13 or not primitive(S):
                continue
            count += 1
            key = canonical_profile_key(S)
            by_profile[key] += 1
            best_dist = best_apex_min_distance(S)
            if best_dist > 1:
                direct_slack += 1
                continue
            M, pts = M_exact(S)
            blocked_exact += 1
            rec = (rems, adds, S, M, pts, key, best_dist)
            if M == THR:
                tight.append(rec)
            elif M < THR:
                below.append(rec)
            if is_apex_locked(S):
                apex_locked.append(rec)
    return {
        "count": count,
        "tight": tight,
        "below": below,
        "apex_locked": apex_locked,
        "by_profile": by_profile,
        "direct_slack": direct_slack,
        "blocked_exact": blocked_exact,
    }


def format_rec(rec: tuple[object, ...]) -> str:
    # Supports both single and two-swap records.
    S = rec[2]  # type: ignore[assignment]
    M = rec[3]  # type: ignore[assignment]
    pts = rec[4]  # type: ignore[assignment]
    return f"S={S}, M={M}, argmax={pts[:4]}{' ...' if len(pts) > 4 else ''}"


def q_covering_necessary(S: tuple[int, ...]) -> bool:
    return all(any(s % q == 0 for s in S) for q in range(2, 15))


def scan_small_covering_box(top: int = 19) -> tuple[int, F, tuple[int, ...], tuple[F, ...], int, int]:
    count = 0
    tight = 0
    below = 0
    best = F(1)
    best_S: tuple[int, ...] = ()
    best_pts: tuple[F, ...] = ()
    for S in combinations(range(1, top + 1), 13):
        if not q_covering_necessary(S):
            continue
        count += 1
        M, pts = M_exact(S)
        if M < best:
            best, best_S, best_pts = M, S, pts
        if M == THR:
            tight += 1
        if M < THR:
            below += 1
    return count, best, best_S, best_pts, tight, below


def main() -> None:
    print("S123 LRC14 STEINHAUS TIGHT-LOCUS SEQUENCE PROBE")
    print("=" * 78)
    print("All values are exact Fractions.  No line below proves the global census.")
    print()

    row_summary("AP tight clock", AP)
    row_summary("Goddyn-Wong tight clock", GW)
    row_summary("near sporadic escape 24->36", tuple(list(range(1, 12)) + [13, 36]))
    row_summary("covering AP-core plus 84", tuple(list(range(1, 12)) + [13, 84]))

    print("[A] Exact single-swap bank around AP")
    print("  bound: replacement v<=80 (sequence-focused; broader v<=300 census exists in S51)")
    single = scan_single_swaps(limit=80)
    tight = single["tight"]
    below = single["below"]
    locked_loose = single["apex_locked_loose"]
    direct = single["apex_slack_direct"]
    blocked_exact = single["apex_blocked_exact"]
    unique_tight = {rec[2]: rec for rec in tight}
    unique_non_ap_tight = [rec for S, rec in unique_tight.items() if S != AP]
    print(f"  rows checked={len(single['rows'])}")
    print(f"  profile keys={len(single['profile_counts'])}")
    print(f"  tight records={len(tight)}; unique tight sets={len(unique_tight)}")
    print(f"  below-threshold rows={len(below)}")
    print(f"  apex-locked but globally loose rows={len(locked_loose)}")
    print(f"  denom-14 direct-slack rows={len(direct)}")
    print(f"  apex-blocked nonlocked rows exact-checked={len(blocked_exact)}")
    print("  non-AP tight rows:")
    if not unique_non_ap_tight:
        print("    none")
    for rec in unique_non_ap_tight:
        print(f"    remove {rec[0]}, add {rec[1]} -> {format_rec(rec)}")
    print("  first apex-locked loose escapes:")
    for rec in locked_loose[:8]:
        print(f"    remove {rec[0]}, add {rec[1]} -> {format_rec(rec)}")
    print()

    print("[B] Exact two-swap AP-neighborhood stress bank")
    print("  bound: replacement values <=18")
    two = scan_two_swaps(limit=18)
    print(f"  rows checked={two['count']}")
    print(f"  profile keys={len(two['by_profile'])}")
    print(f"  denom-14 direct-slack rows={two['direct_slack']}")
    print(f"  apex-blocked rows exact-checked={two['blocked_exact']}")
    print(f"  apex-locked rows={len(two['apex_locked'])}")
    unique_two_tight = {rec[2]: rec for rec in two["tight"]}
    unique_two_non_ap = [rec for S, rec in unique_two_tight.items() if S != AP]
    print(f"  tight records={len(two['tight'])}; unique tight sets={len(unique_two_tight)}")
    print(f"  below-threshold rows={len(two['below'])}")
    print("  non-AP tight rows found:")
    if not unique_two_non_ap:
        print("    none")
    else:
        for rec in unique_two_non_ap[:12]:
            print(f"    {format_rec(rec)}")
    print("  first apex-locked loose two-swap escapes:")
    for rec in two["apex_locked"][:8]:
        if rec not in two["tight"]:
            print(f"    remove {rec[0]}, add {rec[1]} -> {format_rec(rec)}")
    print()

    print("[C] Bounded covering checksum")
    count, best, best_S, best_pts, tight_count, below_count = scan_small_covering_box(top=19)
    print(f"  q-covering rows in [1,19]: {count}")
    print(f"  min M={best} ({float(best):.9f}), row={best_S}, argmax={best_pts[:8]}")
    print(f"  tight rows={tight_count}, below-threshold rows={below_count}")
    print()

    print("[D] Sequence-level readout")
    print("  AP gap word: every unit has partition ((1,14)), energy 0.")
    print("  GW gap words: denominator-14 units keep one collision/missing packet in")
    print("  a small Steinhaus-energy class; global exact-M, not residues alone,")
    print("  isolates 12->24 from many residue-similar perturbations.")
    print("  Proposed sharp target:")
    print("    (1) apex-lock is a finite residue/gap-sequence filter;")
    print("    (2) every non-AP/GW apex-locked reduced row has an off-apex escape")
    print("        M(S)>1/14;")
    print("    (3) scale-separated rows need the existing analytic comb/Weyl input.")
    print("  This is a clearer form of the three-gap rigidity census, not a closure.")


if __name__ == "__main__":
    main()
