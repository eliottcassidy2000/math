#!/usr/bin/env python3
"""S663: apply the sum/product/fraction trinity to LRC and sibling problems.

HYP-2238 says a value should often be stored as a triune carrier:

    (sum packets, product factors, fraction boundary state).

For LRC n=14 this becomes:

    sum face      = active wall / pair-sum packets,
    product face  = C=27 gcd-shell and local obstruction data,
    fraction face = carry vector k in v = r + 27 k, plus owner boundary state.

The finite experiment below tests the slogan on the known floor atoms AP,
Vstar, and nonprimitive 2AP.  Local +27 carry perturbations keep the same
Res_27 shadow and product/gcd shell, but become strict.  The fraction/carry
face is the missing side channel that separates them.
"""

from __future__ import annotations

import sys
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lrc_n14_carry_conservativity_s611 import (  # noqa: E402
    AP,
    C,
    FLOOR,
    N,
    VSTAR,
    exact_maximin,
    exact_score,
    fmt_frac,
)


TWOP = tuple(2 * v for v in AP)


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def dist(x: Fraction) -> Fraction:
    r = frac_part(x)
    return min(r, 1 - r)


def residue_shadow(row: tuple[int, ...]) -> tuple[int, ...]:
    residues = []
    for v in row:
        r = v % C
        if r == 0:
            raise ValueError(f"zero Res_{C} residue in {row}")
        residues.append(r)
    return tuple(sorted(residues))


def carry_pairs(row: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    pairs = []
    for v in row:
        r = v % C
        if r == 0:
            raise ValueError(f"zero Res_{C} residue in {row}")
        pairs.append((r, (v - r) // C))
    return tuple(sorted(pairs))


def shell(r: int) -> int:
    r %= C
    if r == 0:
        return 0
    return min(r, C - r)


def gcd_shell_counts(shadow: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(gcd(shell(r), C) for r in shadow).items()))


def gcd_shell_mass(shadow: tuple[int, ...]) -> int:
    return sum(g * count for g, count in gcd_shell_counts(shadow))


def candidate_times(row: tuple[int, ...]) -> set[Fraction]:
    candidates: set[Fraction] = set()
    for i, a in enumerate(row):
        for m in range(a):
            t = Fraction(2 * m + 1, 2 * a)
            if 0 < t < 1:
                candidates.add(t)
        for b in row[i + 1 :]:
            for den in (a + b, abs(a - b)):
                if den <= 0:
                    continue
                for m in range(1, den):
                    candidates.add(Fraction(m, den))
    return candidates


def maximin_times(row: tuple[int, ...]) -> tuple[Fraction, tuple[Fraction, ...]]:
    best, _ = exact_maximin(tuple(sorted(row)))
    best_times = []
    for t in candidate_times(tuple(sorted(row))):
        if exact_score(tuple(sorted(row)), t) == best:
            best_times.append(t)
    return best, tuple(sorted(best_times))


def active_packets(row: tuple[int, ...], limit: int = 8) -> tuple[tuple[Fraction, tuple[int, ...], tuple[int, ...]], ...]:
    score, times = maximin_times(row)
    packets = []
    for t in times[:limit]:
        active = tuple(v for v in row if dist(Fraction(v) * t) == score)
        pair_sums = tuple(sorted({a + b for a, b in combinations(active, 2)}))
        packets.append((t, active, pair_sums))
    return tuple(packets)


def active_packet_key(row: tuple[int, ...]) -> tuple[tuple[int, tuple[int, ...]], ...]:
    """Small key: denominator plus active pair sums at best times."""

    packets = active_packets(row, limit=12)
    return tuple((t.denominator, sums) for t, _, sums in packets)


def continuant(partials: tuple[int, ...]) -> int:
    """Continuant K(a_1,...,a_m) with K()=1 and K(a_1)=a_1."""

    if not partials:
        return 1
    km2 = 1
    km1 = partials[0]
    for a in partials[1:]:
        km2, km1 = km1, a * km1 + km2
    return km1


def fraction_state(row: tuple[int, ...]) -> dict[str, object]:
    pairs = carry_pairs(row)
    carries = tuple(k for _, k in pairs)
    partials = tuple(k + 1 for k in carries)
    mod14_word = tuple((r - k) % N for r, k in pairs)
    parity_word = tuple((r + k) & 1 for r, k in pairs)
    apex = tuple(r for r, k in pairs if (r - k) % N == 0)
    pair_apex = []
    for (r1, k1), (r2, k2) in combinations(pairs, 2):
        if (r1 + r2 - k1 - k2) % N == 0:
            pair_apex.append((r1, r2))
    return {
        "carry_pairs": pairs,
        "carry_support": tuple(r for r, k in pairs if k),
        "carry_l1": sum(carries),
        "carry_span": max(carries) - min(carries) if carries else 0,
        "continuant": continuant(partials),
        "mod14_word": mod14_word,
        "parity_word": parity_word,
        "apex": apex,
        "pair_apex_count": len(pair_apex),
    }


def route(score: Fraction) -> str:
    if score == FLOOR:
        return "floor"
    if score > FLOOR:
        return "strict"
    return "below"


@dataclass(frozen=True)
class TriuneRecord:
    name: str
    row: tuple[int, ...]
    score: Fraction
    time: Fraction
    shadow_score: Fraction
    shadow_time: Fraction
    shadow: tuple[int, ...]
    active_key: tuple[tuple[int, tuple[int, ...]], ...]
    product_key: tuple[tuple[int, int], ...]
    product_mass: int
    fraction: tuple[tuple[str, object], ...]

    @property
    def route(self) -> str:
        return route(self.score)

    @property
    def shadow_route(self) -> str:
        return route(self.shadow_score)

    @property
    def sum_product_key(self) -> tuple[object, ...]:
        return (self.shadow, self.active_key, self.product_key, self.product_mass)

    @property
    def full_triune_key(self) -> tuple[object, ...]:
        return self.sum_product_key + (self.fraction,)


def triune_record(name: str, row: tuple[int, ...]) -> TriuneRecord:
    row = tuple(sorted(row))
    shadow = residue_shadow(row)
    score, time = exact_maximin(row)
    shadow_score, shadow_time = exact_maximin(shadow)
    frac_state = tuple(sorted(fraction_state(row).items()))
    return TriuneRecord(
        name=name,
        row=row,
        score=score,
        time=time,
        shadow_score=shadow_score,
        shadow_time=shadow_time,
        shadow=shadow,
        active_key=active_packet_key(shadow),
        product_key=gcd_shell_counts(shadow),
        product_mass=gcd_shell_mass(shadow),
        fraction=frac_state,
    )


def perturbations(base_name: str, base: tuple[int, ...], max_weight: int = 2) -> list[TriuneRecord]:
    records = []
    for weight in range(1, max_weight + 1):
        for idxs in combinations(range(len(base)), weight):
            row = list(base)
            moved = []
            for idx in idxs:
                moved.append(base[idx])
                row[idx] += C
            label = f"{base_name}+27({','.join(map(str, moved))})"
            records.append(triune_record(label, tuple(row)))
    return records


def mixed_groups(records: list[TriuneRecord], key_name: str) -> list[tuple[object, Counter[str], list[str]]]:
    groups: defaultdict[object, list[TriuneRecord]] = defaultdict(list)
    for rec in records:
        key = rec.sum_product_key if key_name == "sum_product" else rec.full_triune_key
        groups[key].append(rec)
    mixed = []
    for key, rows in groups.items():
        hist = Counter(rec.route for rec in rows)
        if len(hist) > 1:
            mixed.append((key, hist, [rec.name for rec in rows[:8]]))
    return mixed


@dataclass(frozen=True)
class Application:
    name: str
    sum_face: str
    product_face: str
    fraction_face: str
    next_move: str
    fit: int
    finite_test: int
    proof_leverage: int


APPLICATIONS = [
    Application(
        "LRC14 carry-owner theorem",
        "odd-wall and pair-sum packets",
        "C=27 gcd shells, Pillai mass, local obstructions",
        "carry word k in v=r+27k, owner route, continuant",
        "triune perturbation atlas, then no-leak lemma",
        5,
        5,
        5,
    ),
    Application(
        "OCF / H(T)",
        "odd-cycle packet coefficients",
        "strong-component product and local forbidden H packets",
        "deletion/substitution boundary continuant for Hamiltonian paths",
        "encode macro-word DP as a continuant",
        5,
        4,
        5,
    ),
    Application(
        "Tournament decks",
        "deleted-card loss sums",
        "deck product/scissors component packets",
        "paired card-owner derivative",
        "replace unpaired deck multiset by card+boundary state",
        4,
        5,
        4,
    ),
    Application(
        "Unit distance",
        "unit-edge spine plus tile/bulk packets",
        "direction support and norm/unit-shell factors",
        "point-deletion frontier/ear owner state",
        "triune record for 21-core extension candidates",
        4,
        4,
        5,
    ),
    Application(
        "Finite-field Kakeya/Falconer",
        "distance sums and pinned-distance packets",
        "line-direction local products and concurrency factors",
        "owner of pins / line-choice recursion",
        "separate full-distance-support twins by owner recursion",
        4,
        4,
        4,
    ),
    Application(
        "Goldbach/Lemoine",
        "prime-pair sums E=p+q, O=p+2q",
        "singular-series local obstruction products",
        "ordered pair reconstruction q=O-E, p=2E-O",
        "treat pair reconstruction as fraction/owner face",
        4,
        3,
        3,
    ),
    Application(
        "CH / forcing",
        "cardinal/equinumerosity sums",
        "model-local consistency factors",
        "generic extension boundary state",
        "use forcing state as fraction face in scalar twins",
        3,
        2,
        4,
    ),
    Application(
        "pi/e trace-norm",
        "trace S=e+pi and power sums",
        "norm P=e*pi and discriminant products",
        "branch sheet of T^2-S*T+P",
        "triune trace/norm/branch ledger",
        4,
        3,
        4,
    ),
]


def tournament_from_applications(apps: list[Application]) -> dict[str, object]:
    n = len(apps)
    metrics = [(a.fit, a.finite_test, a.proof_leverage) for a in apps]
    adj = [[0] * n for _ in range(n)]
    out = [0] * n
    for i, j in combinations(range(n), 2):
        wins_i = sum(x > y for x, y in zip(metrics[i], metrics[j]))
        wins_j = sum(x < y for x, y in zip(metrics[i], metrics[j]))
        if wins_i > wins_j or (wins_i == wins_j and apps[i].name < apps[j].name):
            adj[i][j] = 1
            out[i] += 1
        else:
            adj[j][i] = 1
            out[j] += 1

    c3 = 0
    for a, b, c in combinations(range(n), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            c3 += 1

    def sccs() -> list[list[str]]:
        radj = [[i for i in range(n) if adj[i][j]] for j in range(n)]
        seen = [False] * n
        order: list[int] = []
        for start in range(n):
            if seen[start]:
                continue
            stack = [(start, False)]
            while stack:
                v, done = stack.pop()
                if done:
                    order.append(v)
                    continue
                if seen[v]:
                    continue
                seen[v] = True
                stack.append((v, True))
                for w in range(n):
                    if adj[v][w] and not seen[w]:
                        stack.append((w, False))
        seen = [False] * n
        comps: list[list[str]] = []
        for start in reversed(order):
            if seen[start]:
                continue
            comp = []
            q = deque([start])
            seen[start] = True
            while q:
                v = q.popleft()
                comp.append(apps[v].name)
                for w in radj[v]:
                    if not seen[w]:
                        seen[w] = True
                        q.append(w)
            comps.append(comp)
        return comps

    # Hamiltonian path count by DP.
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            val = dp[mask][v]
            if not val:
                continue
            for w in range(n):
                if not mask & (1 << w) and adj[v][w]:
                    dp[mask | (1 << w)][w] += val
    return {
        "outscores": {apps[i].name: out[i] for i in range(n)},
        "score_hist": dict(sorted(Counter(out).items())),
        "directed_3cycles": c3,
        "scc_sizes": sorted([len(c) for c in sccs()], reverse=True),
        "hamiltonian_paths": sum(dp[-1]),
        "top_order": [a.name for a in sorted(apps, key=lambda a: (-a.fit, -a.finite_test, -a.proof_leverage, a.name))],
    }


def main() -> None:
    print("=" * 78)
    print("S663 triune carrier applications: LRC first, then sibling problems")
    print("=" * 78)
    print()
    print("Triune translation")
    print("  sum face      = additive wall / pair-sum packets")
    print("  product face  = C=27 gcd-shell / local obstruction product")
    print("  fraction face = carry-owner boundary state k in v = r + 27k")
    print()

    base_rows = [("AP", AP), ("Vstar", VSTAR), ("2AP", TWOP)]
    base_records = [triune_record(name, row) for name, row in base_rows]
    perturb_records: list[TriuneRecord] = []
    for name, row in base_rows:
        perturb_records.extend(perturbations(name, row, max_weight=2))
    all_records = base_records + perturb_records

    print("A. Base floor triune records")
    for rec in base_records:
        frac_state = dict(rec.fraction)
        print(f"  {rec.name}: M={fmt_frac(rec.score)} t={fmt_frac(rec.time)} route={rec.route}")
        print(f"    shadow={rec.shadow}")
        print(f"    sum active key={rec.active_key[:4]}{' ...' if len(rec.active_key) > 4 else ''}")
        print(f"    product gcd shells={rec.product_key} mass={rec.product_mass}")
        print(
            "    fraction carry_support="
            f"{frac_state['carry_support']} l1={frac_state['carry_l1']} "
            f"span={frac_state['carry_span']} continuant={frac_state['continuant']} "
            f"pair_apex_count={frac_state['pair_apex_count']}"
        )
    print()

    print("B. Local +27 carry perturbation atlas")
    by_family: defaultdict[tuple[str, int], list[TriuneRecord]] = defaultdict(list)
    for rec in perturb_records:
        base = rec.name.split("+", 1)[0]
        weight = rec.name.count(",") + 1
        by_family[(base, weight)].append(rec)
    for key, rows in sorted(by_family.items()):
        hist = Counter(r.route for r in rows)
        min_row = min(rows, key=lambda r: (r.score, r.name))
        print(
            f"  {key[0]:5s} weight={key[1]} count={len(rows):3d} routes={dict(hist)} "
            f"min_M={fmt_frac(min_row.score)} margin={fmt_frac(min_row.score - FLOOR)} "
            f"via {min_row.name}"
        )
    print()

    print("C. Projection-collision audit")
    mixed_sp = mixed_groups(all_records, "sum_product")
    mixed_full = mixed_groups(all_records, "full_triune")
    print(
        "  using only shadow additive packets + product/gcd shell: "
        f"{len(mixed_sp)} mixed floor/strict groups"
    )
    for _, hist, names in mixed_sp[:5]:
        print(f"    mixed routes={dict(hist)} examples={names}")
    print(f"  using full triune key including carry-continuant state: {len(mixed_full)} mixed groups")
    print("  Interpretation: +27 perturbations preserve the Res_27 sum/product shadow,")
    print("  but the carry/fraction face splits floor atoms from strict carries.")
    print()

    print("D. Cross-problem triune application table")
    print(f"{'application':<34} | {'sum face':<34} | {'product face':<38} | fraction face")
    for app in APPLICATIONS:
        print(f"{app.name:<34} | {app.sum_face:<34} | {app.product_face:<38} | {app.fraction_face}")
    print()

    fp = tournament_from_applications(APPLICATIONS)
    print("E. Tournament Analysis over application routes")
    print("  vertices=application routes")
    print("  observable=triune fit, finite testability, proof leverage")
    print(f"  outscores={fp['outscores']}")
    print(f"  score_hist={fp['score_hist']}")
    print(f"  directed_3cycles={fp['directed_3cycles']}")
    print(f"  scc_sizes={fp['scc_sizes']}")
    print(f"  hamiltonian_paths={fp['hamiltonian_paths']}")
    print("  top_order:")
    for name in fp["top_order"]:
        print(f"    {name}")
    print()

    print("F. Working hypotheses")
    print("  1. LRC14: fixed odd-wall packets + fixed C=27 product shell are not enough;")
    print("     the carry/owner continuant is the fraction face that restores sufficiency.")
    print("  2. OCF/tournaments: search for a continuant encoding of deletion/substitution")
    print("     boundary state beside the independence-polynomial packet carrier.")
    print("  3. Unit distance: pair edge-count sums and direction/norm products with")
    print("     point-deletion frontier owners before trusting a candidate core.")
    print("  4. Number theory: Goldbach/Lemoine pair reconstruction and pi/e branch sheets")
    print("     are fraction-face analogues, not just algebraic afterthoughts.")


if __name__ == "__main__":
    main()
