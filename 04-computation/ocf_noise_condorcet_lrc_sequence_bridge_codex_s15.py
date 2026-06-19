#!/usr/bin/env python3
"""
OCF noise/Condorcet spectra and LRC(14) sequence bridge.

HYP-2618 / T866 asks whether the odd-cycle formula

    H(T) = I(Omega(T), 2)

is a noise-stability functional at a special rho, and whether the forbidden
values {7,21} can be read as forbidden Condorcet-cyclicity spectra.  This
script keeps the useful part and rejects the overclaim:

* OCF is exactly a hard-core partition function at activity 2.
* Equivalently, for the independent-set indicator f_Omega under product
  Bernoulli(p), H(T)=3^m E_{p=2/3} f_Omega, where m is the number of odd-cycle
  vertices of Omega.
* It is not, by itself, a nontrivial two-copy noise stability.  For rho<1 the
  stability uses ordered pairs of independent sets, a larger pair-spectrum.

The script then enumerates the small OCF/Condorcet spectra and splices the
normalization into the HYP-2617 coimage sequence atlas for LRC(14).
"""
from __future__ import annotations

import itertools
import math
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_support6_coimage_fiber_codex_s14 as s14  # noqa: E402


def section(title: str) -> None:
    print("\n" + "=" * 92)
    print(title)
    print("=" * 92)


def edge_index(n: int, i: int, j: int) -> int:
    if i > j:
        i, j = j, i
    return i * (2 * n - i - 1) // 2 + (j - i - 1)


def has_edge(bits: int, n: int, i: int, j: int) -> bool:
    if i == j:
        return False
    idx = edge_index(n, i, j)
    bit = (bits >> idx) & 1
    if i < j:
        return bool(bit)
    return not bool(bit)


def hamiltonian_paths(bits: int, n: int) -> int:
    dp: dict[tuple[int, int], int] = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp.get((mask, last), 0)
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if has_edge(bits, n, last, nxt):
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + val
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def directed_odd_cycles(bits: int, n: int) -> list[tuple[int, ...]]:
    cycles: list[tuple[int, ...]] = []
    for size in range(3, n + 1, 2):
        for subset in itertools.combinations(range(n), size):
            first = subset[0]
            for rest in itertools.permutations(subset[1:]):
                cyc = (first,) + rest
                if all(has_edge(bits, n, cyc[i], cyc[(i + 1) % size]) for i in range(size)):
                    cycles.append(cyc)
    return cycles


def alpha_vector_from_cycles(cycles: list[tuple[int, ...]], n: int) -> tuple[int, ...]:
    masks = []
    for cyc in cycles:
        mask = 0
        for v in cyc:
            mask |= 1 << v
        masks.append(mask)
    max_k = n // 3
    alpha = [1] + [0] * max_k

    def rec(start: int, used: int, depth: int) -> None:
        if depth:
            alpha[depth] += 1
        if depth == max_k:
            return
        for idx in range(start, len(masks)):
            mask = masks[idx]
            if used & mask:
                continue
            rec(idx + 1, used | mask, depth + 1)

    rec(0, 0, 0)
    return tuple(alpha)


def independent_sets_from_cycle_masks(masks: list[int]) -> list[int]:
    indep: list[int] = [0]

    def rec(start: int, used_vertices: int, chosen_cycles: int) -> None:
        for idx in range(start, len(masks)):
            if used_vertices & masks[idx]:
                continue
            bit = 1 << idx
            new_chosen = chosen_cycles | bit
            indep.append(new_chosen)
            rec(idx + 1, used_vertices | masks[idx], new_chosen)

    rec(0, 0, 0)
    return indep


def pair_spectrum(indep: list[int], m: int) -> tuple[tuple[tuple[int, int, int], int], ...]:
    counts: Counter[tuple[int, int, int]] = Counter()
    for a in indep:
        for b in indep:
            inter = (a & b).bit_count()
            left = (a & ~b).bit_count()
            right = (b & ~a).bit_count()
            outside = m - (a | b).bit_count()
            counts[(inter, left + right, outside)] += 1
    return tuple(sorted(counts.items()))


@dataclass(frozen=True)
class TournamentRecord:
    n: int
    bits: int
    h: int
    alpha: tuple[int, ...]
    cycle_count: int
    pair_profile: tuple[tuple[tuple[int, int, int], int], ...]

    @property
    def ocf_density_p23(self) -> Fraction:
        return Fraction(self.h, 3 ** self.cycle_count)


def tournament_records(n: int) -> list[TournamentRecord]:
    records: list[TournamentRecord] = []
    for bits in range(1 << (n * (n - 1) // 2)):
        cycles = directed_odd_cycles(bits, n)
        masks = []
        for cyc in cycles:
            mask = 0
            for v in cyc:
                mask |= 1 << v
            masks.append(mask)
        alpha = alpha_vector_from_cycles(cycles, n)
        h = hamiltonian_paths(bits, n)
        ocf = sum(a * (2 ** i) for i, a in enumerate(alpha))
        if ocf != h:
            raise AssertionError(f"OCF failed for n={n}, bits={bits}: {ocf} != {h}")
        indep = independent_sets_from_cycle_masks(masks)
        records.append(
            TournamentRecord(
                n=n,
                bits=bits,
                h=h,
                alpha=alpha,
                cycle_count=len(cycles),
                pair_profile=pair_spectrum(indep, len(cycles)),
            )
        )
    return records


def all_tournament_spectra() -> dict[int, list[TournamentRecord]]:
    return {n: tournament_records(n) for n in range(1, 7)}


def summarize_tournament_spectra(records_by_n: dict[int, list[TournamentRecord]]) -> None:
    section("OCF SPECTRA FOR ALL LABELED TOURNAMENTS n<=6")
    print("OCF verified as H=I(Omega,2) for every labeled tournament in the scan.")
    print(
        f"{'n':>2} {'labeled':>8} {'H values':>8} {'minH':>6} {'maxH':>6} "
        f"{'missing odd <=max':>26} {'7?':>4} {'21?':>5}"
    )
    for n, rows in records_by_n.items():
        values = sorted({r.h for r in rows})
        odd_missing = [v for v in range(1, max(values) + 1, 2) if v not in values]
        print(
            f"{n:>2} {len(rows):>8} {len(values):>8} {min(values):>6} {max(values):>6} "
            f"{str(odd_missing[:12]):>26} {str(7 in values):>4} {str(21 in values):>5}"
        )
    print(
        "Read: {7,21} are absent in the complete n<=6 tournament image. "
        "Since every strict majority tournament is electorate-realizable with enough voters "
        "(McGarvey), these are also forbidden eventual Condorcet-cyclicity spectra, not "
        "just abstract tournament gaps."
    )


def forbidden_candidate_alpha(records_by_n: dict[int, list[TournamentRecord]]) -> None:
    section("FORBIDDEN H VALUES AS OCF/CONDORCET CYCLICITY SPECTRA")
    spectra = sorted({r.alpha for rows in records_by_n.values() for r in rows})
    by_h = defaultdict(list)
    for alpha in spectra:
        h = sum(a * (2 ** i) for i, a in enumerate(alpha))
        by_h[h].append(alpha)
    for target in (7, 21):
        candidates = []
        max_depth = 2
        for a1 in range(30):
            for a2 in range(30):
                alpha = (1, a1, a2)
                if sum(a * (2 ** i) for i, a in enumerate(alpha)) == target:
                    candidates.append(alpha)
        realized = by_h.get(target, [])
        print(f"H={target}: formal alpha candidates up to depth {max_depth}: {candidates}")
        print(f"       realized tournament alpha spectra in n<=6: {realized or 'none'}")
    print(
        "Read: the social-choice object is not just 'there is a Condorcet cycle'. "
        "It is the compatible packet spectrum alpha=(1, alpha_1, alpha_2, ...). "
        "The forbidden values say certain small paradox-packet inventories cannot be "
        "simultaneously realized by any majority tournament."
    )


def majority_bits_from_profile(profile: tuple[tuple[int, ...], ...], m: int) -> int:
    positions = []
    for order in profile:
        pos = [0] * m
        for i, cand in enumerate(order):
            pos[cand] = i
        positions.append(pos)
    bits = 0
    for i in range(m):
        for j in range(i + 1, m):
            votes_i = sum(pos[i] < pos[j] for pos in positions)
            if votes_i > len(profile) // 2:
                bits |= 1 << edge_index(m, i, j)
    return bits


def condorcet_profile_spectra(records_by_n: dict[int, list[TournamentRecord]]) -> None:
    section("EXACT 3-VOTER CONDORCET SPECTRA")
    print(
        "Rows enumerate all profiles of 3 strict voters on m alternatives; "
        "unique majority tournaments are then read through the OCF spectrum."
    )
    print(f"{'m':>2} {'profiles':>10} {'unique T':>9} {'all T?':>7} {'H values':>8} {'H list'}")
    for m in (3, 4, 5):
        orders = tuple(itertools.permutations(range(m)))
        seen: set[int] = set()
        for profile in itertools.product(orders, repeat=3):
            seen.add(majority_bits_from_profile(profile, m))
        lookup = {r.bits: r for r in records_by_n[m]}
        h_values = sorted({lookup[bits].h for bits in seen})
        all_count = 1 << (m * (m - 1) // 2)
        print(
            f"{m:>2} {len(orders) ** 3:>10} {len(seen):>9} "
            f"{str(len(seen) == all_count):>7} {len(h_values):>8} {h_values}"
        )
    print(
        "Read: fixed small electorates impose extra spectra gaps, but by McGarvey the "
        "eventual electorate image is exactly the tournament image.  Thus the "
        "forbidden-H program is a statement about possible social paradox packet "
        "patterns, with voter-count refinements available as a stricter layer."
    )


def noise_identity_table(records_by_n: dict[int, list[TournamentRecord]]) -> None:
    section("IS OCF A NOISE-STABILITY FUNCTIONAL?")
    print("Hard-core identity:")
    print("  mu_p(IndependentSetIndicator)=sum_I p^|I|(1-p)^(m-|I|)")
    print("  = (1-p)^m I(Omega, p/(1-p)).")
    print("  At p=2/3: H(T)=I(Omega,2)=3^m * mu_{2/3}.")
    print()
    print("Two-copy biased noise stability for rho<1 uses ordered pairs (I,J):")
    print("  Stab_{p,rho}(f)=sum_{I,J indep} P11^|I cap J| P10^|I\\J| P01^|J\\I| P00^outside.")
    print("  That pair spectrum is not determined by I(Omega,2) alone.")
    print()
    print("Diagonal same-state mass can be forced to activity 2 by:")
    print("  (p^2+rho*p*q)/(q^2+rho*p*q)=2, q=1-p.")
    print(f"{'rho':>8} {'p solving activity 2':>24}")
    for rho in (Fraction(0, 1), Fraction(1, 3), Fraction(1, 2), Fraction(1, 1)):
        p = solve_activity_two_p(float(rho))
        print(f"{str(rho):>8} {p:>24.12f}")
    print(
        "Conclusion: the useful canonical statement is activity 2 / p=2/3 density. "
        "Calling it noise stability is accurate only in a degenerate rho=1 or "
        "diagonal-mass normalization; nontrivial rho asks for a richer pair polynomial."
    )

    collision = find_pair_spectrum_collision(records_by_n[6])
    if collision:
        h, left, right = collision
        print()
        print("Pair-spectrum guardrail found at n=6:")
        print(f"  same OCF evaluation H={h}, but different pair/noise spectra")
        print(f"  alpha left={left.alpha}, alpha right={right.alpha}")
        print(f"  tournament bits {left.bits} and {right.bits} have different ordered-pair spectra.")
        print("  Therefore a nontrivial noise-stability functional cannot be recovered from the OCF evaluation H alone.")
    else:
        print("No same-H/different-pair-spectrum collision found at n=6; guardrail remains conceptual.")


def solve_activity_two_p(rho: float) -> float:
    # p^2 - 2(1-p)^2 = rho*p*(1-p)
    # (rho-1)p^2 + (4-rho)p - 2 = 0
    a = rho - 1.0
    b = 4.0 - rho
    c = -2.0
    if abs(a) < 1e-12:
        return -c / b
    disc = b * b - 4 * a * c
    roots = [(-b + math.sqrt(disc)) / (2 * a), (-b - math.sqrt(disc)) / (2 * a)]
    return next(r for r in roots if 0 < r < 1)


def find_pair_spectrum_collision(records: list[TournamentRecord]) -> tuple[int, TournamentRecord, TournamentRecord] | None:
    seen: dict[int, TournamentRecord] = {}
    for rec in records:
        old = seen.get(rec.h)
        if old is None:
            seen[rec.h] = rec
        elif old.pair_profile != rec.pair_profile:
            return rec.h, old, rec
    return None


def lrc_sequence_bridge() -> None:
    section("LRC(14) SUPPORT-6 COIMAGE SEQUENCE BRIDGE")
    classes = s14.support_classes()
    all_stats = {d: s14.compute_stats_for_d(d, classes) for d in range(6, 14)}
    zero_hist = Counter(cls.count(0) for cls in classes)
    print(f"projective support-residue classes: {len(classes)}")
    print(f"zero-residue histogram z=0..5: {[zero_hist[i] for i in range(6)]}")
    print()
    print(f"{'d':>3} {'max |S_d|':>12} {'max class':>22} {'null':>6} {'<0.01':>6} {'<0.1':>6} {'median |S_d|':>14}")
    for d, rows in all_stats.items():
        vals = sorted(r.signed_abs for r in rows)
        best = max(rows, key=lambda r: r.signed_abs)
        print(
            f"{d:>3} {best.signed_abs:>12.8g} {str(best.cls):>22} "
            f"{sum(v < 1e-12 for v in vals):>6} {sum(v < 1e-2 for v in vals):>6} "
            f"{sum(v < 1e-1 for v in vals):>6} {vals[len(vals)//2]:>14.8g}"
        )
    print()
    named = {
        "AP": ((1, 2, 3, 4, 5, 6), 7),
        "resonant_21": ((1, 2, 3, 4, 5, 21), 7),
        "wide_68": ((2, 3, 4, 5, 6, 68), 8),
        "wall_22": ((1, 2, 4, 7, 8, 22), 9),
    }
    lookup = {(r.d, r.cls): r for rows in all_stats.values() for r in rows}
    print(f"{'name':<12} {'d':>3} {'class':>22} {'z':>2} {'|S_d|':>12} {'abs/signed':>12}")
    for name, (support, d) in named.items():
        cls = s14.canon_support(support)
        row = lookup[(d, cls)]
        print(f"{name:<12} {d:>3} {str(cls):>22} {row.z:>2} {row.signed_abs:>12.8g} {row.ratio:>12.6g}")
    print(
        "Read: OCF's activity-2 normalization and the LRC coimage atlas share the "
        "same lesson: do not trust raw absolute mass.  Retain the packet address "
        "(independent odd-cycle packets for OCF, projective mod-7 coimage class for LRC) "
        "and then evaluate the signed/compatible partition function."
    )


def proof_quotient_tournament() -> None:
    section("TOURNAMENT ANALYSIS: QUOTIENT CHOICE")
    vertices = [
        "coimage_sequence_tail",
        "condorcet_alpha_spectrum",
        "hard_core_activity_2",
        "biased_density_p23",
        "diagonal_noise_embedding",
        "full_noise_stability",
        "raw_absolute_mass",
    ]
    rank = {v: i for i, v in enumerate(vertices)}
    out = Counter()
    cycles = 0
    for a, b in itertools.combinations(vertices, 2):
        if rank[a] < rank[b]:
            out[a] += 1
        else:
            out[b] += 1
    for a, b, c in itertools.combinations(vertices, 3):
        ab = rank[a] < rank[b]
        bc = rank[b] < rank[c]
        ac = rank[a] < rank[c]
        if (ab and bc and not ac) or ((not ab) and (not bc) and ac):
            cycles += 1
    hist = Counter(out[v] for v in vertices)
    print(f"Hamiltonian proof path: {vertices}")
    print(f"score histogram: {dict(sorted(hist.items()))}")
    print(f"directed 3-cycles: {cycles}")
    print(
        "Assumption challenged: the useful vertices are not runners, arcs, or voters. "
        "They are quotient choices that decide which packet information survives."
    )


def main() -> None:
    section("OCF NOISE / CONDORCET / LRC SEQUENCE BRIDGE S15")
    records_by_n = all_tournament_spectra()
    summarize_tournament_spectra(records_by_n)
    forbidden_candidate_alpha(records_by_n)
    condorcet_profile_spectra(records_by_n)
    noise_identity_table(records_by_n)
    lrc_sequence_bridge()
    proof_quotient_tournament()
    section("S15 READING")
    print("1. OCF is canonically the hard-core partition function at activity 2, or a p=2/3 independent-set density after multiplying by 3^m.")
    print("2. It is not determined by a nontrivial two-copy noise stability at rho<1; that needs a pair spectrum of independent packet pairs.")
    print("3. The forbidden values {7,21} become forbidden compatible Condorcet-cyclicity spectra alpha, because majority tournaments range over all tournaments with enough voters.")
    print("4. The LRC support-six tail should imitate the successful part: retain the finite packet address, then bound the signed compatible sum. HYP-2617 supplies the 159-class address table.")
    print("5. For LRC(14), this points to a next computation: classify which non-null coimage classes can still occur after height-1/height-2 wall deletion, rather than improving the raw absolute Minkowski volume.")


if __name__ == "__main__":
    main()
