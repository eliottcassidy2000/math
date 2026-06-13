#!/usr/bin/env python3
"""Build an add/multiply gauge stack for LRC Tournament Analysis.

codex-2026-06-01 S513

This script extends the S506/S509 Lonely Runner gauge work in the direction
requested by the user: odd/even prime representation fibers, the x+2 / x*2
operation grid, product-sum collision equations, A000568 odd-cycle survival,
and LRC endpoint debt are treated as competing arc criteria.

Tournament Analysis declarations:

Denominator tournaments
    vertices: selected LRC denominators N.
    pairwise observables:
        additive prime fiber size
            even N: unordered Goldbach pairs N=p+q
            odd N: Levy/Lemoine pairs N=p+2q
        multiplicative branch labels: N=2^h * odd_core, divisor count,
            largest proper divisor, and phi(N) endpoint debt.
        product-sum collision data: two-factor solutions of
            (a-1)(b-1)=N-1 and small product-sum seed minima.
        A000568 shadow: number/fraction of odd integer partitions of N.
    switches/gauges:
        scalar gauges orient larger observable to smaller, with increasing N
        as the Hamiltonian tie path.
        majority gauges orient by a pairwise majority vote over named
        observables, again using increasing N for exact ties.
    tie Hamiltonian path: increasing denominator order.

Route tournament
    vertices: candidate LRC arc criteria.
    observable: hand-scored vector of virtues (LRC signal, dynamic data,
        operation-grid relevance, A000568 relevance, proof potential) and
        projection risk.
    switch/gauge: majority vote; lower projection risk wins its coordinate.
    tie Hamiltonian path: listed order of criteria.

Stored output:
    05-knowledge/results/lrc_add_mult_gauge_stack_s513.out
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd


DENOMINATORS = (10, 12, 14, 15, 16, 18, 20, 21, 22, 24, 26, 27, 28, 30)
MAX_GRID_ODD = 19
MAX_GRID_HEIGHT = 4


def fmt_frac(x: Fraction) -> str:
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def fmt_dec(x: Fraction | float, places: int = 3) -> str:
    return f"{float(x):.{places}f}"


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0:
        return False
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def primes_upto(n: int) -> list[int]:
    return [k for k in range(2, n + 1) if is_prime(k)]


def v2(n: int) -> int:
    h = 0
    while n % 2 == 0:
        h += 1
        n //= 2
    return h


def odd_core(n: int) -> int:
    return n >> v2(n)


def divisors(n: int) -> list[int]:
    out: list[int] = []
    for d in range(1, int(n**0.5) + 1):
        if n % d == 0:
            out.append(d)
            if d * d != n:
                out.append(n // d)
    return sorted(out)


def largest_proper_divisor(n: int) -> int:
    if n <= 1:
        return 1
    return max(d for d in divisors(n) if d < n)


def euler_phi(n: int) -> int:
    out = n
    m = n
    p = 2
    while p * p <= m:
        if m % p == 0:
            out -= out // p
            while m % p == 0:
                m //= p
        p += 1 if p == 2 else 2
    if m > 1:
        out -= out // m
    return out


def partition_count(n: int, odd_only: bool = False) -> int:
    parts = range(1, n + 1, 2) if odd_only else range(1, n + 1)
    dp = [0] * (n + 1)
    dp[0] = 1
    for part in parts:
        for total in range(part, n + 1):
            dp[total] += dp[total - part]
    return dp[n]


def goldbach_pairs(n: int) -> list[tuple[int, int]]:
    if n % 2:
        return []
    out = []
    for p in primes_upto(n // 2):
        q = n - p
        if p <= q and is_prime(q):
            out.append((p, q))
    return out


def lemoine_pairs(n: int) -> list[tuple[int, int]]:
    """Prime pairs (p,q) with odd n = p + 2q."""
    if n % 2 == 0:
        return []
    out = []
    for q in primes_upto(n // 2):
        p = n - 2 * q
        if is_prime(p):
            out.append((p, q))
    return out


def two_factor_solutions_for_arity(k: int) -> list[tuple[int, int, int]]:
    """Solutions (a,b,product) to (a-1)(b-1)=k-1, a<=b."""
    out = []
    target = k - 1
    for d in divisors(target):
        e = target // d
        a, b = d + 1, e + 1
        if a <= b:
            out.append((a, b, a * b))
    return out


def enumerate_product_sum_minima(max_k: int, max_product: int = 10000) -> dict[int, tuple[int, tuple[int, ...]]]:
    """Small product-sum minima via nondecreasing nonunit seed search."""
    best: dict[int, tuple[int, tuple[int, ...]]] = {}

    def visit(start: int, product: int, total: int, seed: tuple[int, ...]) -> None:
        if seed:
            defect = product - total
            if defect >= 0:
                arity = len(seed) + defect
                if 2 <= arity <= max_k:
                    current = best.get(arity)
                    if current is None or (product, seed) < current:
                        best[arity] = (product, seed)
        for factor in range(start, max_product + 1):
            new_product = product * factor
            if new_product > max_product:
                break
            new_total = total + factor
            # Once the arity defect is already much too large, larger factors
            # in this nondecreasing branch will not rescue it.
            if seed and new_product - new_total + len(seed) + 1 > max_k + 12:
                break
            visit(factor, new_product, new_total, seed + (factor,))

    visit(2, 1, 0, ())
    return best


PRODUCT_SUM_MINIMA = enumerate_product_sum_minima(max(DENOMINATORS))


@dataclass(frozen=True)
class DenominatorRecord:
    n: int
    parity_kind: str
    add_count: int
    add_examples: tuple[tuple[int, int], ...]
    dyadic_height: int
    odd_core: int
    divisor_count: int
    lpd: int
    phi: int
    twofactor_count: int
    twofactor_best_product: int | None
    ps_min_product: int | None
    ps_seed: tuple[int, ...]
    partitions: int
    odd_partitions: int

    @property
    def odd_partition_fraction(self) -> Fraction:
        return Fraction(self.odd_partitions, self.partitions)

    @property
    def endpoint_debt_fraction(self) -> Fraction:
        return Fraction(self.phi, self.n)


def denominator_record(n: int) -> DenominatorRecord:
    if n % 2 == 0:
        examples = tuple(goldbach_pairs(n))
        parity_kind = "Goldbach"
    else:
        examples = tuple(lemoine_pairs(n))
        parity_kind = "Lemoine"
    twofactor = two_factor_solutions_for_arity(n)
    ps_min_product, ps_seed = PRODUCT_SUM_MINIMA.get(n, (None, ()))
    return DenominatorRecord(
        n=n,
        parity_kind=parity_kind,
        add_count=len(examples),
        add_examples=examples[:3],
        dyadic_height=v2(n),
        odd_core=odd_core(n),
        divisor_count=len(divisors(n)),
        lpd=largest_proper_divisor(n),
        phi=euler_phi(n),
        twofactor_count=len(twofactor),
        twofactor_best_product=min((p for _, _, p in twofactor), default=None),
        ps_min_product=ps_min_product,
        ps_seed=ps_seed,
        partitions=partition_count(n, odd_only=False),
        odd_partitions=partition_count(n, odd_only=True),
    )


RECORDS = tuple(denominator_record(n) for n in DENOMINATORS)


def compare(a, b) -> int:
    return (a > b) - (a < b)


def orient_by_key(records: tuple[DenominatorRecord, ...], key) -> list[list[int]]:
    m = len(records)
    adj = [[0] * m for _ in range(m)]
    for i, j in combinations(range(m), 2):
        c = compare(key(records[i]), key(records[j]))
        if c == 0:
            winner = i
        else:
            winner = i if c > 0 else j
        loser = j if winner == i else i
        adj[winner][loser] = 1
    return adj


def orient_by_majority(records: tuple[DenominatorRecord, ...], dimensions) -> list[list[int]]:
    m = len(records)
    adj = [[0] * m for _ in range(m)]
    for i, j in combinations(range(m), 2):
        votes_i = 0
        votes_j = 0
        for _name, key, high_wins in dimensions:
            c = compare(key(records[i]), key(records[j]))
            if c == 0:
                continue
            if high_wins:
                votes_i += int(c > 0)
                votes_j += int(c < 0)
            else:
                votes_i += int(c < 0)
                votes_j += int(c > 0)
        if votes_i == votes_j:
            winner = i
        else:
            winner = i if votes_i > votes_j else j
        loser = j if winner == i else i
        adj[winner][loser] = 1
    return adj


def count_hamiltonian_paths(adj: list[list[int]]) -> int:
    n = len(adj)
    dp: dict[tuple[int, int], int] = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            ways = dp.get((mask, v), 0)
            if not ways:
                continue
            for w in range(n):
                if not (mask & (1 << w)) and adj[v][w]:
                    dp[(mask | (1 << w), w)] += ways
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def scores(adj: list[list[int]]) -> tuple[int, ...]:
    return tuple(sum(row) for row in adj)


def score_histogram(adj: list[list[int]]) -> str:
    hist = Counter(scores(adj))
    return ",".join(f"{score}:{hist[score]}" for score in sorted(hist))


def directed_triangles(adj: list[list[int]]) -> int:
    out = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            out += 1
    return out


def scc_sizes(adj: list[list[int]]) -> tuple[int, ...]:
    n = len(adj)
    graph = [[j for j in range(n) if adj[i][j]] for i in range(n)]
    rev = [[i for i in range(n) if adj[i][j]] for j in range(n)]
    seen = [False] * n
    order: list[int] = []

    def dfs(v: int) -> None:
        seen[v] = True
        for w in graph[v]:
            if not seen[w]:
                dfs(w)
        order.append(v)

    for v in range(n):
        if not seen[v]:
            dfs(v)

    seen = [False] * n
    sizes: list[int] = []
    for start in reversed(order):
        if seen[start]:
            continue
        q = deque([start])
        seen[start] = True
        size = 0
        while q:
            v = q.pop()
            size += 1
            for w in rev[v]:
                if not seen[w]:
                    seen[w] = True
                    q.append(w)
        sizes.append(size)
    return tuple(sorted(sizes, reverse=True))


def top_vertices(adj: list[list[int]], limit: int = 5) -> str:
    score_pairs = sorted(
        ((scores(adj)[i], RECORDS[i].n) for i in range(len(adj))),
        reverse=True,
    )
    return " ".join(f"{n}:{score}" for score, n in score_pairs[:limit])


def edge_flip_rate(a: list[list[int]], b: list[list[int]]) -> Fraction:
    flips = 0
    total = 0
    for i, j in combinations(range(len(a)), 2):
        total += 1
        flips += int(a[i][j] != b[i][j])
    return Fraction(flips, total)


def scalar_grid_stack_score(record: DenominatorRecord) -> int:
    max_add = max(r.add_count for r in RECORDS)
    first_even_bonus = 4 if record.dyadic_height == 1 else 0
    pure_dyadic_bonus = 3 if record.odd_core == 1 else 0
    ps_bonus = 0 if record.ps_min_product is None else max(0, 10 - record.ps_min_product // 10)
    return (
        5 * (max_add - record.add_count)
        + record.phi
        + 3 * record.dyadic_height
        + 2 * record.twofactor_count
        + first_even_bonus
        + pure_dyadic_bonus
        + ps_bonus
    )


GAUGES = {
    "additive_prime_fiber": orient_by_key(RECORDS, lambda r: (r.add_count, r.phi, -r.n)),
    "additive_scarcity": orient_by_key(RECORDS, lambda r: (-r.add_count, r.phi, -r.n)),
    "multiplicative_branch_depth": orient_by_key(
        RECORDS, lambda r: (r.dyadic_height, r.divisor_count, r.lpd, -r.odd_core)
    ),
    "product_sum_collision": orient_by_key(
        RECORDS,
        lambda r: (
            r.twofactor_count,
            -(r.twofactor_best_product or 10**9),
            -(r.ps_min_product or 10**9),
        ),
    ),
    "lrc_endpoint_debt": orient_by_key(RECORDS, lambda r: (r.phi, r.endpoint_debt_fraction, -r.n)),
    "a000568_odd_survival": orient_by_key(
        RECORDS, lambda r: (r.odd_partition_fraction, r.odd_partitions, -r.n)
    ),
    "scalar_grid_stack": orient_by_key(RECORDS, lambda r: (scalar_grid_stack_score(r), -r.n)),
}

GAUGES["add_mult_majority"] = orient_by_majority(
    RECORDS,
    (
        ("additive abundance", lambda r: r.add_count, True),
        ("divisor branching", lambda r: r.divisor_count, True),
        ("dyadic height", lambda r: r.dyadic_height, True),
        ("product-sum critical pairs", lambda r: r.twofactor_count, True),
        ("endpoint debt", lambda r: r.phi, True),
        ("A000568 odd survival", lambda r: r.odd_partition_fraction, True),
    ),
)

GAUGES["loneliness_pressure_majority"] = orient_by_majority(
    RECORDS,
    (
        ("additive scarcity", lambda r: r.add_count, False),
        ("endpoint debt", lambda r: r.phi, True),
        ("first doubling height", lambda r: r.dyadic_height, True),
        ("small odd core", lambda r: r.odd_core, False),
        ("product-sum critical pairs", lambda r: r.twofactor_count, True),
        ("small product-sum witness", lambda r: r.ps_min_product or 10**9, False),
    ),
)


@dataclass(frozen=True)
class Route:
    name: str
    lrc_signal: int
    dynamic: int
    operation_grid: int
    a000568: int
    proof_potential: int
    projection_risk: int
    note: str


ROUTES = (
    Route("phase_half_H", 7, 8, 2, 5, 6, 5, "global semicircle/entropy meter"),
    Route("theta_close", 8, 8, 2, 3, 6, 4, "local 1/N crowding switch"),
    Route("open_arc_density", 8, 7, 3, 6, 7, 4, "best S506c bridge gauge"),
    Route("danger_deficit", 9, 9, 2, 3, 8, 3, "dynamic shortfall from threshold"),
    Route("edge_danger_deficit", 9, 8, 6, 4, 8, 3, "pair-cell close-burden meter"),
    Route("pressure_relief_SCC", 8, 8, 3, 3, 9, 3, "blocker-debt dependency gauge"),
    Route("endpoint_owner_debt", 10, 8, 4, 4, 10, 2, "keeps labelled LRC boundary fiber"),
    Route("odd_core_row", 3, 1, 9, 8, 4, 7, "horizontal x+2 branch label"),
    Route("dyadic_height_branch", 5, 2, 9, 7, 6, 6, "vertical x*2 branch label"),
    Route("product_sum_defect", 5, 2, 10, 5, 8, 5, "equation tying product to sum"),
    Route("goldbach_lemoine_fiber", 4, 1, 9, 4, 7, 6, "additive prime abundance fiber"),
    Route("A000568_marked_quotient", 7, 5, 5, 10, 9, 2, "rooted quotient/sheaf route"),
    Route("add_mult_gauge_stack", 9, 6, 10, 8, 9, 2, "hybrid metric proposed here"),
)


def route_tournament(routes: tuple[Route, ...]) -> list[list[int]]:
    dims = (
        ("lrc_signal", lambda r: r.lrc_signal, True),
        ("dynamic", lambda r: r.dynamic, True),
        ("operation_grid", lambda r: r.operation_grid, True),
        ("a000568", lambda r: r.a000568, True),
        ("proof_potential", lambda r: r.proof_potential, True),
        ("projection_risk", lambda r: r.projection_risk, False),
    )
    n = len(routes)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        votes_i = 0
        votes_j = 0
        for _name, key, high_wins in dims:
            c = compare(key(routes[i]), key(routes[j]))
            if c == 0:
                continue
            if high_wins:
                votes_i += int(c > 0)
                votes_j += int(c < 0)
            else:
                votes_i += int(c < 0)
                votes_j += int(c > 0)
        if votes_i == votes_j:
            winner = i
        else:
            winner = i if votes_i > votes_j else j
        loser = j if winner == i else i
        adj[winner][loser] = 1
    return adj


def route_top_vertices(adj: list[list[int]], limit: int = 8) -> str:
    sc = scores(adj)
    pairs = sorted(((sc[i], ROUTES[i].name) for i in range(len(ROUTES))), reverse=True)
    return " | ".join(f"{name}:{score}" for score, name in pairs[:limit])


def print_method() -> None:
    print("LRC ADD/MULT GAUGE STACK - S513")
    print("=" * 76)
    print("Denominator vertices:", " ".join(map(str, DENOMINATORS)))
    print("Tie Hamiltonian path: increasing denominator order.")
    print("Even additive fiber: Goldbach pairs N=p+q.")
    print("Odd additive fiber: Levy/Lemoine pairs N=p+2q.")
    print("Multiplicative grid: N = 2^h * odd_core; horizontal odd_core -> odd_core+2, vertical h -> h+1.")
    print("Product-sum bridge: P-S=D and D ones repair the equation; two-factor layer is (a-1)(b-1)=N-1.")
    print("A000568 bridge: unlabelled tournament Burnside terms survive on odd cycle partitions; here we track odd partition survival as a coarse shadow.")
    print()


def print_operation_grid() -> None:
    print("x+2 / x*2 OPERATION GRID")
    print("=" * 76)
    header = "odd_core | " + " | ".join(f"h={h}" for h in range(MAX_GRID_HEIGHT + 1))
    print(header)
    print("-" * len(header))
    for odd in range(1, MAX_GRID_ODD + 1, 2):
        cells = []
        for h in range(MAX_GRID_HEIGHT + 1):
            n = odd * (2**h)
            if n % 2 == 0:
                count = len(goldbach_pairs(n))
                tag = f"{n}:G{count}"
            else:
                count = len(lemoine_pairs(n))
                tag = f"{n}:L{count}"
            cells.append(tag)
        print(f"{odd:>8} | " + " | ".join(cells))
    print()
    print("Reading: horizontal motion keeps h fixed and adds 2 to the odd core; vertical motion doubles.")
    print("Goldbach/Lemoine counts are additive fibers sitting on these multiplicative addresses.")
    print()


def print_records() -> None:
    print("DENOMINATOR RECORDS")
    print("=" * 76)
    header = (
        "N  core h kind      add ex        phi div lpd tf ps_min seed        oddpart"
    )
    print(header)
    print("-" * len(header))
    for r in RECORDS:
        examples = ",".join(f"{a}+{b}" if r.parity_kind == "Goldbach" else f"{a}+2*{b}" for a, b in r.add_examples)
        examples = examples if examples else "-"
        ps = "-" if r.ps_min_product is None else str(r.ps_min_product)
        seed = "-" if not r.ps_seed else "*".join(map(str, r.ps_seed))
        print(
            f"{r.n:>2} {r.odd_core:>5} {r.dyadic_height:>1} {r.parity_kind:<9} "
            f"{r.add_count:>3} {examples:<13} {r.phi:>3} {r.divisor_count:>3} "
            f"{r.lpd:>3} {r.twofactor_count:>2} {ps:>6} {seed:<10} "
            f"{r.odd_partitions:>3}/{r.partitions}"
        )
    print()


def print_gauge_fingerprints() -> None:
    print("DENOMINATOR TOURNAMENT FINGERPRINTS")
    print("=" * 76)
    header = "gauge                         H        c3  SCCs        score_hist             top scores"
    print(header)
    print("-" * len(header))
    for name, adj in GAUGES.items():
        h = count_hamiltonian_paths(adj)
        c3 = directed_triangles(adj)
        scc = scc_sizes(adj)
        print(
            f"{name:<29} {h:>8} {c3:>5} {str(scc):<11} "
            f"{score_histogram(adj):<22} {top_vertices(adj)}"
        )
    print()
    print("Scalar gauges mostly collapse to transitive orders. Majority gauges keep the operation conflict visible.")
    print()


def print_edge_flips() -> None:
    names = list(GAUGES)
    print("EDGE-FLIP RATES BETWEEN DENOMINATOR GAUGES")
    print("=" * 76)
    print(" " * 29 + " ".join(f"{name[:6]:>6}" for name in names))
    for i, name_a in enumerate(names):
        row = [f"{name_a:<29}"]
        for j, name_b in enumerate(names):
            if i == j:
                row.append("    --")
            else:
                row.append(f"{fmt_dec(edge_flip_rate(GAUGES[name_a], GAUGES[name_b]), 2):>6}")
        print(" ".join(row))
    print()


def print_route_tournament() -> None:
    adj = route_tournament(ROUTES)
    print("CANDIDATE ARC-CRITERIA ROUTE TOURNAMENT")
    print("=" * 76)
    print(f"vertices={len(ROUTES)} H={count_hamiltonian_paths(adj)} c3={directed_triangles(adj)} SCCs={scc_sizes(adj)}")
    print(f"score_hist={score_histogram(adj)}")
    print("top route scores:", route_top_vertices(adj))
    print()
    print("route                         LRC dyn op  A568 proof risk  note")
    print("-" * 76)
    sc = scores(adj)
    for i, r in sorted(enumerate(ROUTES), key=lambda item: (sc[item[0]], item[1].name), reverse=True):
        print(
            f"{r.name:<29} {r.lrc_signal:>3} {r.dynamic:>3} {r.operation_grid:>3} "
            f"{r.a000568:>5} {r.proof_potential:>5} {r.projection_risk:>4}  {r.note}"
        )
    print()


def print_synthesis() -> None:
    print("SYNTHESIS")
    print("=" * 76)
    print("1. Addition is the flattening operation: as a shadow it becomes the total order x<z, but its fibers remember p+q and p+2q abundance.")
    print("2. Multiplication is the branching operation: divisors, dyadic height, and odd core build the sparse vertical skeleton under that order.")
    print("3. Product-sum equations are collision equations between those two operations: multiplication creates surplus D=P-S, and addition repairs it by D ones.")
    print("4. A000568 enters because tournament quotienting forgets labels by odd-cycle fixed-point survival; the LRC analogue must keep marked endpoint fibers over that quotient.")
    print("5. A useful loneliness metric should therefore be a stack, not one number: dynamic LRC threshold gauges, pair-cell danger, endpoint debt, odd/dyadic branch labels, and additive/product-sum fiber entropy.")
    print()
    print("Proposed practical arc criteria:")
    print("  runner level: phase_half_H, theta_close, open_arc_density, danger_deficit, pressure_relief_SCC.")
    print("  pair-cell level: edge_danger_deficit plus odd_core_row, dyadic_height_branch, product_sum_defect labels.")
    print("  denominator level: Goldbach/Lemoine scarcity, phi(N) endpoint debt, A000568 odd survival, and add_mult_gauge_stack majority.")
    print()
    print("Conjectural proof shape: a counterexample row would need to be scarce in additive fibers, high in endpoint debt, trapped in the marked A000568 fiber, and cyclic/nonpeeling in pressure gauges at the same time.")


def main() -> None:
    print_method()
    print_operation_grid()
    print_records()
    print_gauge_fingerprints()
    print_edge_flips()
    print_route_tournament()
    print_synthesis()


if __name__ == "__main__":
    main()
