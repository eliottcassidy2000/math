#!/usr/bin/env python3
"""
S678: rank easier proof targets near LRC n=14 and unit-distance n=21.

This is a comparative proof-target atlas, not a proof of any new LRC or
unit-distance case.  The question is whether there are neighboring or later
values whose proof carriers look easier than the headline frontiers.

The answer depends on the predicate:

* LRC full theorem difficulty still grows with runner count, but the local
  THM-401 shell at C=2n-1 can be cleaner after n=14.
* Unit-distance lower-bound/spine certificates can be easier than exact
  maximum proofs.  THM-408 makes n=13 and n=14 one-slab Moser rows, while
  n=21 is the first minus-family two-slab row.

Tournament Analysis / assumption challenge:
  Vertices are proof obligations rather than runners or point sets.  Alternate
  vertex sets considered here include runners, residues mod 2n-1, gcd-shells,
  fixed clocks, Moser points, Moser slabs, unit directions, ears, centered
  Eisenstein shells, and exact-upper-bound obstructions.  The chosen quotient
  preserves "which proof obligation should be easier" and destroys raw
  embedding and continuous-time phase data.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from itertools import combinations
from math import gcd


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    p = 3
    while p * p <= n:
        if n % p == 0:
            return False
        p += 2
    return True


def factor_dict(n: int) -> dict[int, int]:
    out: dict[int, int] = {}
    x = n
    p = 2
    while p * p <= x:
        while x % p == 0:
            out[p] = out.get(p, 0) + 1
            x //= p
        p += 1 if p == 2 else 2
    if x > 1:
        out[x] = out.get(x, 0) + 1
    return out


def fmt_factor(n: int) -> str:
    fac = factor_dict(n)
    if not fac:
        return "1"
    return " * ".join(str(p) if e == 1 else f"{p}^{e}" for p, e in fac.items())


def folded_orbits(C: int) -> tuple[tuple[int, ...], ...]:
    """Orbits of nonzero residues under x -> 2x and x -> -x mod C."""
    unseen = set(range(1, C))
    orbits: list[tuple[int, ...]] = []
    while unseen:
        root = min(unseen)
        orbit: set[int] = set()
        queue = deque([root])
        while queue:
            x = queue.popleft() % C
            if x == 0 or x in orbit:
                continue
            orbit.add(x)
            for y in ((2 * x) % C, (-x) % C):
                if y and y not in orbit:
                    queue.append(y)
        orbits.append(tuple(sorted(orbit)))
        unseen -= orbit
    return tuple(orbits)


def centered_hex_radius(n: int) -> int | None:
    # Centered triangular/Eisenstein hex shells: 1 + 3r(r+1).
    r = 0
    while 1 + 3 * r * (r + 1) <= n:
        if 1 + 3 * r * (r + 1) == n:
            return r
        r += 1
    return None


@dataclass(frozen=True)
class LrcTarget:
    n: int
    C: int
    C_factor: str
    C_omega: int
    C_max_exp: int
    orbit_count: int
    gcd_strata: int
    n_odd_prime: bool
    C_prime: bool
    three_shell: bool
    doubling_sporadic: bool
    proven_window: bool
    local_score: int
    full_score: int
    reason: str


def lrc_target(n: int) -> LrcTarget:
    C = 2 * n - 1
    fac = factor_dict(C)
    orbits = folded_orbits(C)
    gcd_strata = len({gcd(orbit[0], C) for orbit in orbits})
    n_odd_prime = n % 2 == 1 and is_prime(n)
    C_prime = is_prime(C)
    three_shell = C % 3 == 0
    doubling_sporadic = n % 2 == 0 and n % 3 == 2
    proven_window = n <= 13

    # Local carrier score: lower is easier in the repo's current residue/carry
    # machinery.  It intentionally does not decide the whole LRC theorem.
    local_score = 0
    local_score += 3 if not n_odd_prime else -2
    local_score += 1 if n % 2 == 0 else 0
    local_score += 3 if not C_prime else -2
    local_score += 2 * max(0, len(fac) - 1)
    local_score += max(0, max(fac.values()) - 1) if fac else 0
    local_score += max(0, gcd_strata - 1)
    local_score += 3 if three_shell else 0
    local_score += 5 if doubling_sporadic else 0
    local_score = max(local_score, -4)

    # Full-proof score adds size and frontier status.  Smaller proved rows are
    # honestly easier, while later clean carriers are "easier subproblems" only.
    full_score = local_score + max(0, n - 13)
    if proven_window:
        full_score -= 8
    if n == 14:
        full_score += 4

    reasons: list[str] = []
    if proven_window:
        reasons.append("already in proved window")
    if n_odd_prime:
        reasons.append("odd-prime n keeps the paper-style tight shortcut alive")
    if C_prime:
        reasons.append(f"C={C} prime gives no gcd-strata side channel")
    elif C == 27:
        reasons.append("C=27=3^3 is the active Res_27 wall")
    elif three_shell:
        reasons.append("C has a 3-shell")
    if doubling_sporadic:
        reasons.append("Vstar doubling-sporadic wall is present")
    if not reasons:
        reasons.append("mixed but not carrying the n=14 triadic wall")

    return LrcTarget(
        n=n,
        C=C,
        C_factor=fmt_factor(C),
        C_omega=len(fac),
        C_max_exp=max(fac.values()) if fac else 0,
        orbit_count=len(orbits),
        gcd_strata=gcd_strata,
        n_odd_prime=n_odd_prime,
        C_prime=C_prime,
        three_shell=three_shell,
        doubling_sporadic=doubling_sporadic,
        proven_window=proven_window,
        local_score=local_score,
        full_score=full_score,
        reason="; ".join(reasons),
    )


@dataclass(frozen=True)
class UnitDistanceTarget:
    n: int
    family: str
    m: int | None
    expected_edges: int | None
    spine_edges: int | None
    bulk_edges: int | None
    centered_hex_radius: int | None
    nearest_moser: str
    lower_spine_score: int
    exact_upper_score: int
    reason: str


def moser_family(n: int) -> tuple[str, int, int] | None:
    if n >= 6 and (n - 6) % 8 == 0:
        m = (n - 6) // 8
        edges = 8 if m == 0 else 27 * m + 6
        return ("P_m^+", m, edges)
    if n >= 5 and (n - 5) % 8 == 0:
        m = (n - 5) // 8
        edges = 6 if m == 0 else 27 * m + 3
        return ("P_m^-", m, edges)
    return None


def nearest_moser_label(n: int) -> tuple[str, int]:
    candidates: list[tuple[int, str]] = []
    for m in range(0, 5):
        candidates.append((abs(n - (8 * m + 6)), f"P_{m}^+"))
        candidates.append((abs(n - (8 * m + 5)), f"P_{m}^-"))
    dist, label = min(candidates, key=lambda item: (item[0], item[1]))
    return label, dist


def unit_distance_target(n: int) -> UnitDistanceTarget:
    fam = moser_family(n)
    hex_r = centered_hex_radius(n)
    nearest, ear_debt = nearest_moser_label(n)

    if fam:
        family, m, edges = fam
        spine = n - 1
        bulk = edges - spine
    else:
        family, m, edges, spine, bulk = ("none", None, None, None, None)

    lower_spine_score = 0
    exact_upper_score = 0
    reasons: list[str] = []

    if fam:
        lower_spine_score += 2 * (m or 0)
        exact_upper_score += 3 * (m or 0)
        reasons.append(f"explicit THM-408 {family} row")
        if (m or 0) <= 1:
            reasons.append("only one full slab or less")
        else:
            lower_spine_score += 2
            exact_upper_score += 6
            reasons.append("two-or-more slabs expose ear/core side channels")
    else:
        lower_spine_score += 5 + ear_debt
        exact_upper_score += 6 + 2 * ear_debt
        reasons.append(f"{ear_debt} vertices from nearest Moser row {nearest}")

    if hex_r is not None:
        lower_spine_score -= 3
        exact_upper_score -= 1
        reasons.append(f"centered Eisenstein shell radius {hex_r}")

    if n <= 14:
        lower_spine_score -= 3
        exact_upper_score -= 5
        reasons.append("small exact-witness regime already heavily audited")

    if n in (21, 22):
        exact_upper_score += 5
        reasons.append("frontier two-slab exact-upper obstruction")

    lower_spine_score += max(0, n - 21) // 4
    if n > 21:
        # Exact planar upper bounds should pay a serious size/embedding tax.
        # A larger row may have a pretty lower-bound spine, but that does not
        # make its extremal upper-bound problem easier than n=21.
        exact_upper_score += 10 + (n - 21)

    return UnitDistanceTarget(
        n=n,
        family=family,
        m=m,
        expected_edges=edges,
        spine_edges=spine,
        bulk_edges=bulk,
        centered_hex_radius=hex_r,
        nearest_moser=nearest,
        lower_spine_score=lower_spine_score,
        exact_upper_score=exact_upper_score,
        reason="; ".join(reasons),
    )


@dataclass(frozen=True)
class Route:
    name: str
    preserves_predicate: int
    side_channel_retention: int
    computational_leverage: int
    theorem_value: int
    overclaim_risk: int


ROUTES = [
    Route("split_local_carrier_from_full_theorem", 5, 5, 4, 5, 0),
    Route("LRC_clean_C_later_values", 5, 4, 5, 4, 1),
    Route("UD_one_slab_spine_values", 5, 4, 5, 4, 1),
    Route("UD_centered_Eisenstein_shells", 4, 5, 3, 3, 2),
    Route("raw_smaller_n_is_easier", 3, 2, 5, 3, 1),
    Route("literal_14_21_numerology", 1, 1, 5, 1, 5),
]


def route_score(route: Route) -> int:
    return (
        3 * route.preserves_predicate
        + 2 * route.side_channel_retention
        + 2 * route.theorem_value
        + route.computational_leverage
        - 3 * route.overclaim_risk
    )


def beats(a: Route, b: Route) -> bool:
    sa = route_score(a)
    sb = route_score(b)
    if sa != sb:
        return sa > sb
    return a.name < b.name


def directed_triangles(routes: list[Route]) -> int:
    total = 0
    for a, b, c in combinations(range(len(routes)), 3):
        wins = []
        for x, y in ((a, b), (a, c), (b, c)):
            wins.append((x, y) if beats(routes[x], routes[y]) else (y, x))
        if sorted(Counter(x for x, _ in wins).values()) == [1, 1, 1]:
            total += 1
    return total


def hamiltonian_paths(routes: list[Route]) -> int:
    n = len(routes)
    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp.get((mask, last), 0)
            if not count:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if beats(routes[last], routes[nxt]):
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def yesno(flag: bool) -> str:
    return "yes" if flag else "no"


def print_lrc() -> None:
    targets = [lrc_target(n) for n in range(3, 25)]
    base = next(t for t in targets if t.n == 14)
    print("LRC target atlas")
    print("----------------")
    print(f"baseline n=14: C={base.C}={base.C_factor}, local_score={base.local_score}, full_score={base.full_score}")
    print()
    print("n | C | C factor | orbits | gcd strata | odd-prime n | C prime | 3-shell | Vstar wall | local | full | note")
    print("-- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | --")
    for t in targets:
        marker = "*" if t.n == 14 else ""
        print(
            f"{t.n}{marker} | {t.C} | {t.C_factor} | {t.orbit_count} | {t.gcd_strata} | "
            f"{yesno(t.n_odd_prime)} | {yesno(t.C_prime)} | {yesno(t.three_shell)} | "
            f"{yesno(t.doubling_sporadic)} | {t.local_score} | {t.full_score} | {t.reason}"
        )
    print()
    print("Locally cleaner later LRC carriers than n=14")
    print("n | C | local_score | why it is interesting | warning")
    print("-- | -- | -- | -- | --")
    later = [t for t in targets if t.n > 14 and t.local_score < base.local_score]
    for t in sorted(later, key=lambda x: (x.local_score, x.n))[:10]:
        warning = "larger full LRC theorem, not an immediate frontier shortcut"
        print(f"{t.n} | {t.C}={t.C_factor} | {t.local_score} | {t.reason} | {warning}")
    print()


def print_unit_distance() -> None:
    targets = [unit_distance_target(n) for n in range(3, 31)]
    base = next(t for t in targets if t.n == 21)
    print("Unit-distance target atlas")
    print("--------------------------")
    print(
        "baseline n=21: "
        f"family={base.family}, m={base.m}, edges={base.expected_edges}, "
        f"lower_spine_score={base.lower_spine_score}, exact_upper_score={base.exact_upper_score}"
    )
    print()
    print("n | family | m | edges | spine | bulk | hex shell | lower-spine | exact-upper | note")
    print("-- | -- | -- | -- | -- | -- | -- | -- | -- | --")
    for t in targets:
        marker = "*" if t.n == 21 else ""
        print(
            f"{t.n}{marker} | {t.family} | {t.m if t.m is not None else '-'} | "
            f"{t.expected_edges if t.expected_edges is not None else '-'} | "
            f"{t.spine_edges if t.spine_edges is not None else '-'} | "
            f"{t.bulk_edges if t.bulk_edges is not None else '-'} | "
            f"{t.centered_hex_radius if t.centered_hex_radius is not None else '-'} | "
            f"{t.lower_spine_score} | {t.exact_upper_score} | {t.reason}"
        )
    print()
    print("Easier unit-distance values than n=21 by predicate")
    print("predicate | values | reading")
    print("-- | -- | --")
    easier_spine = [t.n for t in targets if t.n != 21 and t.lower_spine_score < base.lower_spine_score]
    easier_exact = [t.n for t in targets if t.n != 21 and t.exact_upper_score < base.exact_upper_score]
    print(
        "explicit/lower-bound unit spine | "
        f"{easier_spine} | n=13 and n=14 are the clean one-slab Moser rows; n=19 is a symmetric Eisenstein shell"
    )
    print(
        "exact-upper style proof | "
        f"{easier_exact} | smaller audited rows are easier; n=22 is not easier here because the 60/61 side channel remains"
    )
    print()


def print_tournament_analysis() -> None:
    routes = sorted(ROUTES, key=lambda r: (-route_score(r), r.name))
    scores = {route.name: route_score(route) for route in routes}
    outdegrees = {
        r.name: sum(1 for other in routes if other != r and beats(r, other))
        for r in routes
    }
    print("Tournament Analysis over target-selection routes")
    print("------------------------------------------------")
    print(f"vertices={len(routes)}")
    print(f"score_hist={dict(sorted(Counter(outdegrees.values()).items()))}")
    print(f"directed_3cycles={directed_triangles(routes)}")
    print(f"hamiltonian_paths={hamiltonian_paths(routes)}")
    print("tie Hamiltonian path:")
    for r in routes:
        print(f"  outdegree={outdegrees[r.name]} route_score={scores[r.name]} {r.name}")
    print()


def main() -> None:
    print("S678 easier frontier target atlas")
    print("=================================")
    print()
    print_lrc()
    print_unit_distance()
    print("Synthesis")
    print("---------")
    print("LRC: n<=13 are the genuinely easier/proved values.  Among later values,")
    print("n=15,16,18,19,22,24 have cleaner local C=2n-1 carriers than n=14 in")
    print("this scorecard; n=19 is the prettiest because n and C=37 are both prime.")
    print("The caution is that this is local-carrier ease, not a proof that LRC(19)")
    print("is globally easier than finishing the immediate n=14 frontier.")
    print()
    print("Unit distance: for lower-bound or unit-spine certificates, n=13 and n=14")
    print("are easier than n=21 because they are the one-slab THM-408 rows.  The")
    print("centered Eisenstein shell n=19 is a promising symmetric side target.")
    print("For exact upper-bound work, n=22 is not easier than n=21 despite having")
    print("an explicit P_2^+ spine, because the repo's n=22 obstacle is precisely")
    print("the extra endpoint/ear side channel.")
    print()
    print_tournament_analysis()
    print("Assumption challenge")
    print("--------------------")
    print("This session does not use runners or points as tournament vertices.")
    print("It uses proof obligations: clean C-shells, odd-prime shortcut rows,")
    print("Moser slab count, Eisenstein symmetry, and exact-upper side channels.")
    print("The quotient preserves which theorem route becomes simpler and destroys")
    print("the raw continuous LRC orbit and the concrete planar embedding.")


if __name__ == "__main__":
    main()
