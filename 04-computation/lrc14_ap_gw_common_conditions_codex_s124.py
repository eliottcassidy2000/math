#!/usr/bin/env python3
"""S124: AP/Goddyn-Wong common-condition atlas for LRC14.

Goal: derive exact necessary-condition candidates by asking what the two known
tight rows share:

    AP = {1,...,13}
    GW = {1,...,11,13,24}

The first filters below are rigorous necessary conditions for a primitive
tight row, using existing repo theorems:
  * punctured q-cover: multiples for every q=2..13, but no multiple of 14;
  * apex unit cover: every unit residue mod 14 is represented.

The later filters are stronger fingerprints suggested by AP/GW:
  * exact unit shell and exact half-step residue 7;
  * cofinite nonzero residue support and maximal antipodal zsum;
  * AP complement binders have actual sums 14, not merely 0 mod 14;
  * residue defect is either zero (AP) or one even-shell dipole (GW-like);
  * one-petal acceleration: the dipole is made by replacing m with 2m.

These later filters are not claimed as theorems.  They are proposed proof
targets / obstructions: any future tight row must either satisfy them or teach
us exactly where AP/GW intuition fails.
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
UNITS = (1, 3, 5, 9, 11, 13)
EVENS = (2, 4, 6, 8, 10, 12)
ODD_SKELETON = (1, 3, 5, 7, 9, 11, 13)
AP = tuple(range(1, 14))
GW = tuple(list(range(1, 12)) + [13, 24])
RESIDUE_LIAR_26 = tuple(list(range(1, 12)) + [13, 26])
NEAR_MISS_36 = tuple(list(range(1, 12)) + [13, 36])


def norm1(x: F) -> F:
    r = x % 1
    return min(r, 1 - r)


def primitive(S: tuple[int, ...]) -> bool:
    g = 0
    for s in S:
        g = gcd(g, s)
    return g == 1


def residue_counts(S: tuple[int, ...]) -> tuple[int, ...]:
    c = [0] * N
    for s in S:
        c[s % N] += 1
    return tuple(c)


def q_punctured_cover(S: tuple[int, ...]) -> bool:
    return all(any(s % q == 0 for s in S) for q in range(2, 14)) and all(s % 14 for s in S)


def q_threshold(S: tuple[int, ...]) -> int:
    """q(S) = first divisor denominator with no multiple in S."""
    q = 2
    while any(s % q == 0 for s in S):
        q += 1
    return q


def missing_divisors(S: tuple[int, ...], top: int = 14) -> tuple[int, ...]:
    return tuple(q for q in range(2, top + 1) if not any(s % q == 0 for s in S))


def same_residue_multiset_as_ap(S: tuple[int, ...]) -> bool:
    return residue_counts(S) == residue_counts(AP)


def unit_cover(S: tuple[int, ...]) -> bool:
    c = residue_counts(S)
    return all(c[u] >= 1 for u in UNITS)


def exact_unit_shell(S: tuple[int, ...]) -> bool:
    c = residue_counts(S)
    return c[0] == 0 and all(c[u] == 1 for u in UNITS)


def exact_odd_skeleton(S: tuple[int, ...]) -> bool:
    c = residue_counts(S)
    return c[0] == 0 and all(c[r] == 1 for r in ODD_SKELETON)


def gcd_shell_profile(S: tuple[int, ...]) -> bool:
    c = residue_counts(S)
    return (
        c[0] == 0
        and sum(c[u] for u in UNITS) == 6
        and sum(c[e] for e in EVENS) == 6
        and c[7] == 1
    )


def cofinite_residue_set(S: tuple[int, ...]) -> bool:
    """AP/GW kind: nonzero residue support misses at most one residue."""
    c = residue_counts(S)
    return c[0] == 0 and sum(1 for r in range(1, N) if c[r]) >= 12


def zsum_max(S: tuple[int, ...]) -> bool:
    """Maximal antipodal pair count modulo 14.

    AP has the six literal antipodal pairs.  GW loses (2,12) but the doubled
    residue 10 adds (4,24), so the total stays six.
    """
    return sum(1 for i, a in enumerate(S) for b in S[i + 1 :] if (a + b) % N == 0) == 6


def complement_sum14(S: tuple[int, ...]) -> bool:
    """The unit-shell binders are the literal AP complement pairs.

    THM-568/HYP-2909 naturally gives sums divisible by 14.  AP and GW share the
    stronger literal fact that the unique speed in residue r and the unique speed
    in residue 14-r sum to exactly 14 for r=1,3,5.
    """
    by_res: dict[int, list[int]] = {r: [] for r in range(N)}
    for s in S:
        by_res[s % N].append(s)
    for r in (1, 3, 5):
        a = by_res[r]
        b = by_res[N - r]
        if len(a) != 1 or len(b) != 1:
            return False
        if a[0] + b[0] != 14:
            return False
    return True


def even_dipole_data(S: tuple[int, ...]) -> tuple[str, tuple[int, ...], tuple[int, ...]]:
    """Return AP/GW-style defect type against the AP residue clock."""
    c = residue_counts(S)
    defect = [c[r] - (1 if r else 0) for r in range(N)]
    support = tuple(r for r, d in enumerate(defect) if d)
    missing = tuple(r for r, d in enumerate(defect) if d < 0)
    extra = tuple(r for r, d in enumerate(defect) if d > 0)
    if not support:
        return ("zero", missing, extra)
    if set(support).issubset(set(EVENS)) and len(missing) == 1 and len(extra) == 1 and sum(abs(defect[r]) for r in support) == 2:
        return ("single_even_dipole", missing, extra)
    return ("other", missing, extra)


def single_even_dipole_or_zero(S: tuple[int, ...]) -> bool:
    return even_dipole_data(S)[0] in {"zero", "single_even_dipole"}


def one_petal_or_ap(S: tuple[int, ...]) -> bool:
    """AP or the GW-style acceleration m -> 2m.

    This is deliberately strong and specific.  For GW, residue 12 is missing,
    residue 10 is doubled, and the out-of-window representative is 24=2*12.
    The near-miss 36=3*12 fails here.
    """
    kind, missing, extra = even_dipole_data(S)
    if kind == "zero":
        return True
    if kind != "single_even_dipole":
        return False
    m = missing[0]
    e = extra[0]
    outside = sorted(s for s in S if s > 13 and s % 14 == e)
    return len(outside) == 1 and outside[0] == 2 * m


def minimal_reps_except_one_petal(S: tuple[int, ...]) -> bool:
    """All residue representatives are minimal, except possibly one m -> 2m petal."""
    kind, missing, extra = even_dipole_data(S)
    if kind == "zero":
        return set(S) == set(AP)
    if kind != "single_even_dipole":
        return False
    m = missing[0]
    e = extra[0]
    expected = (set(AP) - {m}) | {2 * m}
    return set(S) == expected and (2 * m) % 14 == e


def top_petal_or_ap(S: tuple[int, ...]) -> bool:
    """The known LRC14 top petal: AP or 12 -> 24."""
    return set(S) == set(AP) or set(S) == set(GW)


FILTERS = (
    ("punctured_q_cover", q_punctured_cover),
    ("unit_cover", unit_cover),
    ("exact_unit_shell", exact_unit_shell),
    ("exact_odd_skeleton", exact_odd_skeleton),
    ("gcd_shell_profile", gcd_shell_profile),
    ("cofinite_residue_set", cofinite_residue_set),
    ("zsum_max", zsum_max),
    ("complement_sum14", complement_sum14),
    ("even_dipole_or_zero", single_even_dipole_or_zero),
    ("one_petal_or_ap", one_petal_or_ap),
    ("minimal_one_petal", minimal_reps_except_one_petal),
    ("top_petal_or_ap", top_petal_or_ap),
)


def first_failed_filter(S: tuple[int, ...]) -> str:
    for name, pred in FILTERS:
        if not pred(S):
            return name
    return "PASS"


def gw_window(v: int) -> tuple[int, int]:
    """Goddyn-Wong acceleration window for AP {1,...,13} at site v."""
    return (N - v, 2 * N - 1 - 2 * v)


def gw_window_coprime_blockers(v: int) -> tuple[int, ...]:
    lo, hi = gw_window(v)
    return tuple(x for x in range(lo, hi + 1) if gcd(x, v) == 1)


def gw_gate_passes(v: int) -> bool:
    return not gw_window_coprime_blockers(v)


def candidate_taus(S: tuple[int, ...]) -> set[F]:
    S = tuple(sorted(set(S)))
    cands: set[F] = {F(1, 2)}
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
    return cands


_M_CACHE: dict[tuple[int, ...], tuple[F, tuple[F, ...]]] = {}


def M_exact(S: tuple[int, ...]) -> tuple[F, tuple[F, ...]]:
    S = tuple(sorted(set(S)))
    if S in _M_CACHE:
        return _M_CACHE[S]
    best = F(0)
    pts: list[F] = []
    for t in candidate_taus(S):
        val = min(norm1(s * t) for s in S)
        if val > best:
            best = val
            pts = [t]
        elif val == best:
            pts.append(t)
    out = (best, tuple(sorted(pts)))
    _M_CACHE[S] = out
    return out


def row_signature(S: tuple[int, ...]) -> str:
    M, pts = M_exact(S)
    kind, missing, extra = even_dipole_data(S)
    q = q_threshold(S)
    return (
        f"S={S}, q={q}, M={M}, denoms={sorted({t.denominator for t in pts})}, "
        f"defect={kind}, missing={missing}, extra={extra}"
    )


def bank_single_swaps(limit: int = 300) -> list[tuple[int, ...]]:
    rows = set()
    for rem in AP:
        kept = set(AP) - {rem}
        for v in range(1, limit + 1):
            S = tuple(sorted(kept | {v}))
            if len(S) == 13 and primitive(S):
                rows.add(S)
    return sorted(rows)


def bank_two_swaps(limit: int = 40) -> list[tuple[int, ...]]:
    rows = set()
    vals = range(1, limit + 1)
    for rems in combinations(AP, 2):
        kept = set(AP) - set(rems)
        for adds in combinations(vals, 2):
            S = tuple(sorted(kept | set(adds)))
            if len(S) == 13 and primitive(S):
                rows.add(S)
    return sorted(rows)


def bank_box(top: int = 19) -> list[tuple[int, ...]]:
    return [S for S in combinations(range(1, top + 1), 13) if primitive(S)]


def filter_counts(rows: list[tuple[int, ...]]) -> tuple[list[set[tuple[int, ...]]], list[tuple[int, ...]]]:
    survivors: list[set[tuple[int, ...]]] = []
    current = set(rows)
    for _, pred in FILTERS:
        current = {S for S in current if pred(S)}
        survivors.append(set(current))
    return survivors, sorted(current)


def report_bank(name: str, rows: list[tuple[int, ...]], exact_final: bool = True) -> None:
    print(f"[{name}]")
    print(f"  total primitive rows={len(rows)}")
    survivors, final_rows = filter_counts(rows)
    prev = len(rows)
    for (fname, _), surv in zip(FILTERS, survivors):
        print(f"  after {fname:<20}: {len(surv):5d}  (drop {prev - len(surv):5d})")
        prev = len(surv)
    print(f"  final survivors={len(final_rows)}")
    for S in final_rows[:12]:
        print(f"    {row_signature(S) if exact_final else S}")
    if len(final_rows) > 12:
        print(f"    ... {len(final_rows)-12} more")
    if exact_final:
        tight = []
        loose = []
        below = []
        for S in final_rows:
            M, _ = M_exact(S)
            if M == THR:
                tight.append(S)
            elif M < THR:
                below.append(S)
            else:
                loose.append(S)
        print(f"  final exact M: tight={len(tight)}, loose={len(loose)}, below={len(below)}")
    print()


def report_q_threshold_rows() -> None:
    print("[q-threshold / residue-liar audit]")
    print("  q(S)=min d>=2 with no multiple of d in S; the q-witness gives M(S)>=1/q(S).")
    print("  rows with the same mod-14 residues can split immediately by divisibility.")
    rows = (
        ("AP", AP),
        ("GW 12->24", GW),
        ("residue_liar 12->26", RESIDUE_LIAR_26),
        ("Farey near_miss 12->36", NEAR_MISS_36),
    )
    for label, S in rows:
        M, pts = M_exact(S)
        q = q_threshold(S)
        det = 1 * 41 - 3 * 14
        farey = " yes" if M == F(3, 41) else " no"
        print(
            f"  {label:<24} q={q:<2} q_witness=1/{q:<2} M={str(M):<5} "
            f"same_residues_as_AP={same_residue_multiset_as_ap(S)} "
            f"missing<=14={missing_divisors(S)} first_fail={first_failed_filter(S)} "
            f"denoms={sorted({t.denominator for t in pts})} farey_3_41={farey}"
        )
        if M == F(3, 41):
            print(f"    Farey check: det[[1,3],[14,41]]={det}.")
    print()


def report_petal_ledger() -> None:
    print("[Minimal AP-doubling petal ledger]")
    print("  Rows are AP with one valid replacement v -> 2v.  The GW gate is the")
    print("  interval [14-v, 27-2v]; it passes exactly when no integer in the")
    print("  interval is coprime to v.")
    print(
        f"  {'v':>2} {'2v':>3} {'q':>2} {'gate':>5} {'window':>9} "
        f"{'coprime blockers':>19} {'M':>7} {'denoms':>10} {'first failed filter'}"
    )
    for v in AP:
        new = 2 * v
        S = tuple(sorted((set(AP) - {v}) | {new}))
        if len(S) != 13 or new in AP:
            continue
        M, pts = M_exact(S)
        lo, hi = gw_window(v)
        blockers = gw_window_coprime_blockers(v)
        print(
            f"  {v:>2} {new:>3} {q_threshold(S):>2} {str(gw_gate_passes(v)):>5} "
            f"{str((lo, hi)):>9} {str(blockers):>19} {str(M):>7} "
            f"{str(sorted({t.denominator for t in pts})):>10} {first_failed_filter(S)}"
        )
    print("  Only v=12 passes the Jacobsthal/Goddyn-Wong gate and has M=1/14.")
    print()


def survivors_after(rows: list[tuple[int, ...]], filter_name: str) -> list[tuple[int, ...]]:
    current = set(rows)
    for name, pred in FILTERS:
        current = {S for S in current if pred(S)}
        if name == filter_name:
            return sorted(current)
    raise ValueError(f"unknown filter {filter_name!r}")


def report_terminal_census_shrink() -> None:
    print("[Terminal finite-census shrink]")
    print("  The broad banks test the filters, but the final exact-M census is much")
    print("  smaller: after the minimal-one-petal premise, both local banks have the")
    print("  same four-row terminal core.")

    single_wide = bank_single_swaps(300)
    two_wide = bank_two_swaps(40)
    terminal_single = survivors_after(single_wide, "minimal_one_petal")
    terminal_two = survivors_after(two_wide, "minimal_one_petal")
    terminal = sorted(set(terminal_single) | set(terminal_two))
    print(f"  single replacements v<=300 terminal rows: {len(terminal_single)}")
    print(f"  two replacements values<=40 terminal rows: {len(terminal_two)}")
    print(f"  union terminal core: {len(terminal)} rows, max element {max(max(S) for S in terminal)}")
    for S in terminal:
        M, pts = M_exact(S)
        verdict = "TIGHT" if M == THR else "loose"
        print(
            f"    {row_signature(S)} first_fail={first_failed_filter(S)} "
            f"argdenoms={sorted({t.denominator for t in pts})} {verdict}"
        )

    print("  Stabilization by replacement ceiling:")
    print("    limit | single terminal | two-replacement terminal | terminal max")
    for limit in (16, 20, 24, 30, 40):
        s_rows = survivors_after(bank_single_swaps(limit), "minimal_one_petal")
        t_rows = survivors_after(bank_two_swaps(limit), "minimal_one_petal")
        union = sorted(set(s_rows) | set(t_rows))
        max_seen = max((max(S) for S in union), default=None)
        print(f"    {limit:5d} | {len(s_rows):15d} | {len(t_rows):24d} | {max_seen}")
    print(
        "  Therefore, inside the AP/GW minimal-petal proof target, the exact terminal"
    )
    print(
        "  census is bounded by max speed 24; the larger banks only stress-test"
    )
    print("  earlier necessary filters and near-miss escape modes.")
    print()


def condition_tournament(rows: list[tuple[int, ...]]) -> None:
    """Tournament Analysis over conditions.

    Vertices are condition filters.  For filters A,B, orient A -> B when A is
    stronger than B on the sampled bank (A's pass set is a proper subset of B's
    pass set).  Ties use the listed order as Hamiltonian path.
    """
    pass_sets = []
    for _, pred in FILTERS:
        pass_sets.append({S for S in rows if pred(S)})
    n = len(FILTERS)
    adj = [[False] * n for _ in range(n)]
    flips = 0
    incomparable = 0
    for i in range(n):
        for j in range(i + 1, n):
            Ai, Aj = pass_sets[i], pass_sets[j]
            if Ai < Aj:
                adj[i][j] = True
            elif Aj < Ai:
                adj[j][i] = True
                flips += 1
            elif Ai == Aj:
                adj[i][j] = True
            else:
                incomparable += 1
                # Orient by smaller pass-set first, then tie order.
                if (len(Ai), i) <= (len(Aj), j):
                    adj[i][j] = True
                else:
                    adj[j][i] = True
                    flips += 1
    scores = [sum(row) for row in adj]
    cyc3 = 0
    for a, b, c in combinations(range(n), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            cyc3 += 1
        if adj[a][c] and adj[c][b] and adj[b][a]:
            cyc3 += 1
    print("[Tournament Analysis: condition-strength carrier]")
    print("  vertices:", ", ".join(name for name, _ in FILTERS))
    print("  observable: pass-set inclusion over single-swap bank v<=300")
    print("  switch: stronger condition -> weaker condition; ties use listed path")
    print(f"  score_hist={dict(sorted(Counter(scores).items()))}")
    print(f"  directed_3cycles={cyc3}, incomparable_pairs={incomparable}, opposite_to_list_flips={flips}")
    print()


def main() -> None:
    print("S124 AP/GW COMMON NECESSARY-CONDITION ATLAS")
    print("=" * 78)
    print("Known tight rows:")
    for label, S in (
        ("AP", AP),
        ("GW", GW),
        ("residue_liar_26", RESIDUE_LIAR_26),
        ("near_miss_36", NEAR_MISS_36),
    ):
        flags = ", ".join(name for name, pred in FILTERS if pred(S))
        print(f"  {label:<12} {row_signature(S)}")
        print(f"    passes: {flags}")
    print()

    report_q_threshold_rows()
    report_petal_ledger()
    report_terminal_census_shrink()

    single = bank_single_swaps()
    report_bank("AP single replacements v<=300", single)
    report_bank("AP two replacements values<=40", bank_two_swaps())
    report_bank("bounded primitive 13-subsets [1,19]", bank_box())
    condition_tournament(single)

    print("[Derived necessary-condition stack]")
    print("  Theorem-level: q=2..13 covered and q=14 omitted; otherwise a q-grid witness")
    print("  gives M>1/14, or the denominator-14 apex is blocked.")
    print("  Theorem-level: unit residues must be represented; otherwise some a/14 is")
    print("  strictly lonely.")
    print("  AP/GW fingerprint: cofinite nonzero residue support, maximal antipodal")
    print("  zsum, exact odd skeleton, AP complement sums 14, and all residue defect")
    print("  lives as zero or one even-shell dipole.")
    print("  Strong speculative fingerprint: residue representatives are minimal except")
    print("  for one acceleration petal m -> 2m; the only tight petal seen is the top")
    print("  one, 12 -> 24.  Other petals survive coarse filters but escape off apex.")
    print("  In these banks the full top-petal stack leaves only AP/GW (or AP alone")
    print("  when 24 is outside the bank), with no below-threshold row.")


if __name__ == "__main__":
    main()
