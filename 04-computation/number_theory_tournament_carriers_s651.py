#!/usr/bin/env python3
"""S651: number-theory tournament carriers.

This script answers a narrow version of "how do tournaments fit into number
theory?" without pretending to solve RH, BSD, abc, Goldbach, or twin primes.

The working rule is:

    number-theory tournament = pairwise quotient of a retained arithmetic carrier.

Three concrete tests are included:

1. Prime-race tournaments: vertices are reduced residue classes mod q, and
   a -> b when pi(X;q,a) currently exceeds pi(X;q,b).
2. Local obstruction ledgers: vertices are local residue channels for prime
   tuples / Goldbach forms, priced by survivor factors.
3. Proof-route tournament: vertices are proof methods/carriers for hard
   number-theory problems, ranked by side-channel retention rather than raw
   scalar numerology.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from math import gcd, prod


MAX_X = 100_000
RACE_QS = [3, 4, 5, 7, 8, 11, 13, 16]
RACE_XS = [1_000, 10_000, 100_000]
LOCAL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]


def primes_upto(n: int) -> list[int]:
    sieve = [True] * (n + 1)
    if n >= 0:
        sieve[0] = False
    if n >= 1:
        sieve[1] = False
    for p in range(2, int(n**0.5) + 1):
        if sieve[p]:
            start = p * p
            sieve[start : n + 1 : p] = [False] * (((n - start) // p) + 1)
    return [i for i, ok in enumerate(sieve) if ok]


PRIMES = primes_upto(MAX_X)


def reduced_residues(q: int) -> list[int]:
    return [a for a in range(q) if gcd(a, q) == 1]


def race_counts(q: int, x: int) -> dict[int, int]:
    residues = reduced_residues(q)
    counts = {a: 0 for a in residues}
    for p in PRIMES:
        if p > x:
            break
        r = p % q
        if r in counts:
            counts[r] += 1
    return counts


def edge_key(i: int, j: int) -> tuple[int, int]:
    return (i, j) if i < j else (j, i)


def race_tournament(labels: list[int], counts: dict[int, int]) -> dict[tuple[int, int], int]:
    """Return edges on index pairs.  Value 1 means lower index beats higher."""

    edges: dict[tuple[int, int], int] = {}
    for i, j in combinations(range(len(labels)), 2):
        a = labels[i]
        b = labels[j]
        ca = counts[a]
        cb = counts[b]
        if ca > cb:
            edges[(i, j)] = 1
        elif cb > ca:
            edges[(i, j)] = 0
        else:
            # Tie Hamiltonian path: lower residue wins.  Ties are reported.
            edges[(i, j)] = 1 if a < b else 0
    return edges


def dynamic_race_tournament(
    labels: list[int], q: int, x: int
) -> tuple[dict[tuple[int, int], int], dict[tuple[int, int], int]]:
    """Pairwise-majority race over prime checkpoints up to x.

    The final-count race is always transitive because one scalar count orders
    every residue.  This dynamic race asks a genuinely pairwise question:
    for each pair of residues, which one led at more prime checkpoints?
    """

    counts = {a: 0 for a in labels}
    margins = {pair: 0 for pair in combinations(range(len(labels)), 2)}
    for p in PRIMES:
        if p > x:
            break
        r = p % q
        if r in counts:
            counts[r] += 1
        for i, j in combinations(range(len(labels)), 2):
            a = labels[i]
            b = labels[j]
            if counts[a] > counts[b]:
                margins[(i, j)] += 1
            elif counts[b] > counts[a]:
                margins[(i, j)] -= 1
            else:
                margins[(i, j)] += 1 if a < b else -1

    final_counts = counts
    edges: dict[tuple[int, int], int] = {}
    for i, j in combinations(range(len(labels)), 2):
        margin = margins[(i, j)]
        if margin > 0:
            edges[(i, j)] = 1
        elif margin < 0:
            edges[(i, j)] = 0
        else:
            # Tie path: final count, then lower residue.
            a = labels[i]
            b = labels[j]
            if final_counts[a] > final_counts[b]:
                edges[(i, j)] = 1
            elif final_counts[b] > final_counts[a]:
                edges[(i, j)] = 0
            else:
                edges[(i, j)] = 1 if a < b else 0
    return edges, margins


def beats(edges: dict[tuple[int, int], int], i: int, j: int) -> bool:
    if i < j:
        return edges[(i, j)] == 1
    return edges[(j, i)] == 0


def score_hist(edges: dict[tuple[int, int], int], n: int) -> dict[int, int]:
    scores = Counter()
    for i in range(n):
        scores[i] = 0
    for i, j in combinations(range(n), 2):
        winner = i if beats(edges, i, j) else j
        scores[winner] += 1
    return dict(sorted(Counter(scores.values()).items()))


def directed_3cycles(edges: dict[tuple[int, int], int], n: int) -> int:
    total = 0
    for a, b, c in combinations(range(n), 3):
        out = Counter()
        for i, j in [(a, b), (a, c), (b, c)]:
            out[i if beats(edges, i, j) else j] += 1
        if sorted(out.values()) == [1, 1, 1]:
            total += 1
    return total


def hamiltonian_count(edges: dict[tuple[int, int], int], n: int) -> int:
    @lru_cache(maxsize=None)
    def dp(mask: int, last: int) -> int:
        if mask == (1 << last):
            return 1
        prev_mask = mask ^ (1 << last)
        total = 0
        for prev in range(n):
            if ((prev_mask >> prev) & 1) and beats(edges, prev, last):
                total += dp(prev_mask, prev)
        return total

    full = (1 << n) - 1
    return sum(dp(full, last) for last in range(n))


def scc_sizes(edges: dict[tuple[int, int], int], n: int) -> list[int]:
    graph = [[] for _ in range(n)]
    rev = [[] for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if beats(edges, i, j):
            graph[i].append(j)
            rev[j].append(i)
        else:
            graph[j].append(i)
            rev[i].append(j)

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
    sizes = []

    def rdfs(v: int) -> int:
        seen[v] = True
        size = 1
        for w in rev[v]:
            if not seen[w]:
                size += rdfs(w)
        return size

    for v in reversed(order):
        if not seen[v]:
            sizes.append(rdfs(v))
    return sorted(sizes, reverse=True)


def edge_flips(a: dict[tuple[int, int], int], b: dict[tuple[int, int], int]) -> int:
    return sum(1 for k in a if a[k] != b[k])


def ranking_from_counts(counts: dict[int, int]) -> list[int]:
    return sorted(counts, key=lambda r: (-counts[r], r))


def inverse_mod(a: int, p: int) -> int:
    return pow(a % p, -1, p)


def linear_survivor_row(forms: tuple[tuple[int, int], ...], p: int) -> tuple[int, int, tuple[int, ...], float]:
    """Forms are a*n+b.  Return p, survivors, forbidden residues, density."""

    forbidden = set()
    for a, b in forms:
        aa = a % p
        bb = b % p
        if aa == 0:
            if bb == 0:
                return p, 0, tuple(range(p)), 0.0
            continue
        forbidden.add((-bb * inverse_mod(aa, p)) % p)
    survivors = p - len(forbidden)
    return p, survivors, tuple(sorted(forbidden)), survivors / p


def local_ledger(forms: tuple[tuple[int, int], ...], primes: list[int]) -> list[dict[str, object]]:
    k = len(forms)
    partial = 1.0
    rows = []
    for p in primes:
        _p, survivors, forbidden, density = linear_survivor_row(forms, p)
        if p == 0:
            normalized = 0.0
        else:
            denom = (1.0 - 1.0 / p) ** k
            normalized = (survivors / p) / denom if denom else 0.0
        partial *= normalized
        rows.append(
            {
                "p": p,
                "survivors": survivors,
                "forbidden": forbidden,
                "density": density,
                "normalized_factor": normalized,
                "partial_product": partial,
            }
        )
    return rows


@dataclass(frozen=True)
class LocalProblem:
    name: str
    forms: tuple[tuple[int, int], ...]
    preserved_predicate: str
    destroyed_data: str


LOCAL_PROBLEMS = [
    LocalProblem(
        "twin_primes_gap_2",
        ((1, 0), (1, 2)),
        "n and n+2 survive every local prime channel",
        "global primality, parity barrier, distribution level",
    ),
    LocalProblem(
        "prime_triplet_0_2_6",
        ((1, 0), (1, 2), (1, 6)),
        "n,n+2,n+6 is locally admissible",
        "global primality and higher-order correlations",
    ),
    LocalProblem(
        "prime_quadruplet_0_2_6_8",
        ((1, 0), (1, 2), (1, 6), (1, 8)),
        "prime quadruplet pattern is locally admissible",
        "global primality and tuple-count asymptotics",
    ),
    LocalProblem(
        "sophie_germain",
        ((1, 0), (2, 1)),
        "n and 2n+1 survive every local prime channel",
        "global primality and safe-prime density",
    ),
    LocalProblem(
        "goldbach_N_210",
        ((1, 0), (-1, 210)),
        "a and 210-a survive every local prime channel",
        "major/minor arcs and representation count",
    ),
    LocalProblem(
        "goldbach_N_2110",
        ((1, 0), (-1, 2110)),
        "a and 2110-a survive every local prime channel",
        "major/minor arcs and representation count",
    ),
]


@dataclass(frozen=True)
class ProblemCarrier:
    name: str
    raw_scalar: str
    hidden_witness: str
    side_channels: tuple[str, ...]
    candidate_vertices: tuple[str, ...]
    tags: tuple[str, ...]


PROBLEM_CARRIERS = [
    ProblemCarrier(
        "Riemann_Hypothesis",
        "all nontrivial zeta zeros on Re(s)=1/2",
        "cancellation in primes via explicit formula",
        ("zero spectrum", "prime races", "Mellin/Fourier modes", "Euler product"),
        ("residue classes", "zeros", "Fourier modes", "proof obligations"),
        ("zero-spectrum", "prime-race", "global-cancellation", "Euler-product"),
    ),
    ProblemCarrier(
        "Generalized_RH_Chebotarev",
        "zeros / Frobenius equidistribution",
        "character and conjugacy-class cancellation",
        ("Dirichlet characters", "Frobenius classes", "Artin L-functions", "prime races"),
        ("characters", "conjugacy classes", "residue classes", "zeros"),
        ("zero-spectrum", "prime-race", "Frobenius", "local-global"),
    ),
    ProblemCarrier(
        "Twin_Primes",
        "gap 2 occurs infinitely often",
        "admissible pair survives local obstructions and parity barrier",
        ("local residues", "singular series", "sieve weights", "level of distribution"),
        ("local primes", "residue classes", "gaps", "sieve weights"),
        ("local-obstruction", "sieve", "parity-barrier", "prime-tuples"),
    ),
    ProblemCarrier(
        "Goldbach",
        "every large even N has p+q=N",
        "additive convolution of primes stays positive",
        ("local residues", "major arcs", "minor arcs", "representation counts"),
        ("residue pairs", "prime-pair edges", "circle arcs", "proof obligations"),
        ("local-obstruction", "additive-cover", "circle-method", "prime-tuples"),
    ),
    ProblemCarrier(
        "BSD",
        "analytic rank equals algebraic rank",
        "L-function vanishing order matches Mordell-Weil/Selmer data",
        ("Frobenius traces", "local root numbers", "Selmer groups", "regulators", "Tamagawa factors"),
        ("primes of good reduction", "Frobenius traces", "Selmer classes", "local conditions"),
        ("Frobenius", "local-global", "zero-spectrum", "height-regulator"),
    ),
    ProblemCarrier(
        "abc",
        "height controlled by radical up to epsilon",
        "valuation-side channel prevents too much additive/multiplicative alignment",
        ("prime valuations", "radical", "height", "conductor/discriminant", "Szpiro channel"),
        ("prime divisors", "valuation strata", "height packets", "elliptic curve conductors"),
        ("height-valuation", "local-global", "discriminant", "side-channel"),
    ),
    ProblemCarrier(
        "Euler_Heegner_Rabinowitsch",
        "finite lucky prime horizon",
        "class-number-one/no-early-collapse side channel",
        ("discriminant", "class group", "norm factorization", "boundary square"),
        ("input slots", "boundary positions", "class group states", "proof obligations"),
        ("finite-window", "discriminant", "local-global", "side-channel"),
    ),
    ProblemCarrier(
        "Collatz",
        "all orbits reach 1",
        "2-adic/3-adic carry defect remains controlled",
        ("parity word", "rapidity defect", "2-adic state", "3-adic inverse tree"),
        ("parity states", "defect packets", "residue classes", "proof obligations"),
        ("arithmetic-dynamics", "height-valuation", "side-channel", "finite-window"),
    ),
]


@dataclass(frozen=True)
class Route:
    name: str
    preserves_predicate: int
    computable_now: int
    cross_problem: int
    proof_strength: int
    risk: int


ROUTES = [
    Route("side_channel_jackknife", 5, 5, 5, 4, 1),
    Route("local_obstruction_product_ledger", 5, 5, 4, 4, 1),
    Route("prime_race_tournament", 4, 5, 4, 3, 2),
    Route("Frobenius_trace_tournament", 5, 3, 5, 5, 2),
    Route("zero_mode_spectral_packet", 5, 3, 5, 5, 3),
    Route("height_valuation_ledger", 4, 4, 4, 4, 2),
    Route("finite_exception_window", 4, 4, 3, 4, 1),
    Route("two_shadow_trace_norm_carrier", 4, 5, 4, 4, 1),
    Route("raw_scalar_numerology", 1, 5, 2, 1, 5),
]

CRITERIA = ["preserves_predicate", "computable_now", "cross_problem", "proof_strength"]


def route_tournament() -> tuple[dict[tuple[int, int], int], Counter[str]]:
    edges: dict[tuple[int, int], int] = {}
    wins = Counter({r.name: 0 for r in ROUTES})
    for i, j in combinations(range(len(ROUTES)), 2):
        a = ROUTES[i]
        b = ROUTES[j]
        av = 0
        bv = 0
        for criterion in CRITERIA:
            if getattr(a, criterion) > getattr(b, criterion):
                av += 1
            elif getattr(b, criterion) > getattr(a, criterion):
                bv += 1
        if a.risk < b.risk:
            av += 1
        elif b.risk < a.risk:
            bv += 1
        winner = a if av >= bv else b
        edges[(i, j)] = 1 if winner is a else 0
        wins[winner.name] += 1
    return edges, wins


def problem_similarity_rows() -> list[tuple[float, str, str, list[str]]]:
    rows = []
    for a, b in combinations(PROBLEM_CARRIERS, 2):
        A = set(a.tags)
        B = set(b.tags)
        shared = sorted(A & B)
        if not shared:
            continue
        score = len(shared) / len(A | B)
        rows.append((score, a.name, b.name, shared))
    return sorted(rows, reverse=True)


def local_problem_tournament() -> tuple[dict[tuple[int, int], int], Counter[str], dict[str, float]]:
    """Orient local problems by how much local mass survives after primes <=31.

    This is a toy diagnostic, not a theorem: a higher product means the local
    obstruction ledger is less restrictive, so it wins the "local abundance"
    gauge.
    """

    products: dict[str, float] = {}
    for problem in LOCAL_PROBLEMS:
        rows = local_ledger(problem.forms, LOCAL_PRIMES)
        products[problem.name] = float(rows[-1]["partial_product"])
    edges: dict[tuple[int, int], int] = {}
    wins = Counter({p.name: 0 for p in LOCAL_PROBLEMS})
    for i, j in combinations(range(len(LOCAL_PROBLEMS)), 2):
        a = LOCAL_PROBLEMS[i]
        b = LOCAL_PROBLEMS[j]
        av = products[a.name]
        bv = products[b.name]
        if av > bv or (av == bv and a.name < b.name):
            edges[(i, j)] = 1
            wins[a.name] += 1
        else:
            edges[(i, j)] = 0
            wins[b.name] += 1
    return edges, wins, products


def print_prime_races(lines: list[str]) -> None:
    lines.append("Prime-Race Tournaments")
    lines.append("----------------------")
    lines.append("Observable: pi(X;q,a)-pi(X;q,b) for reduced residue classes a,b mod q.")
    lines.append("Switch/gauge: orient a -> b when pi(X;q,a) > pi(X;q,b); ties go to lower residue.")
    lines.append("Tie Hamiltonian path: residue classes sorted by (-count, residue).")
    lines.append("Dynamic variant: orient a -> b when a led b at more prime checkpoints up to X.")
    lines.append("")
    for q in RACE_QS:
        labels = reduced_residues(q)
        previous_edges = None
        lines.append(f"q={q}, vertices={labels}")
        for x in RACE_XS:
            counts = race_counts(q, x)
            edges = race_tournament(labels, counts)
            n = len(labels)
            ties = sum(1 for i, j in combinations(labels, 2) if counts[i] == counts[j])
            flips = "-" if previous_edges is None else str(edge_flips(previous_edges, edges))
            previous_edges = edges
            ranking = ranking_from_counts(counts)
            lines.append(
                "  X={:<6} ranking={} counts={} ties={} flips={} score_hist={} c3={} scc={} H={}".format(
                    x,
                    ranking,
                    {r: counts[r] for r in ranking},
                    ties,
                    flips,
                    score_hist(edges, n),
                    directed_3cycles(edges, n),
                    scc_sizes(edges, n),
                    hamiltonian_count(edges, n),
                )
            )
        dynamic_edges, margins = dynamic_race_tournament(labels, q, RACE_XS[-1])
        final_edges = race_tournament(labels, race_counts(q, RACE_XS[-1]))
        n = len(labels)
        strongest = sorted(
            ((abs(v), labels[i], labels[j], v) for (i, j), v in margins.items()),
            reverse=True,
        )[:6]
        lines.append(
            "  dynamic-majority X={} flips_vs_final={} score_hist={} c3={} scc={} H={}".format(
                RACE_XS[-1],
                edge_flips(dynamic_edges, final_edges),
                score_hist(dynamic_edges, n),
                directed_3cycles(dynamic_edges, n),
                scc_sizes(dynamic_edges, n),
                hamiltonian_count(dynamic_edges, n),
            )
        )
        lines.append(
            "  strongest dynamic pair margins: "
            + ", ".join(f"{a}/{b}:{margin}" for _abs, a, b, margin in strongest)
        )
        lines.append("")


def print_local_ledgers(lines: list[str]) -> None:
    lines.append("Local Obstruction Ledgers")
    lines.append("-------------------------")
    lines.append("Observable: survivor-factor comparison across local prime channels.")
    lines.append("Switch/gauge for the local-problem tournament: larger partial normalized product wins.")
    lines.append("Tie Hamiltonian path: descending partial product, then name.")
    lines.append("")
    for problem in LOCAL_PROBLEMS:
        rows = local_ledger(problem.forms, LOCAL_PRIMES)
        lines.append(f"{problem.name}")
        lines.append(f"  preserved predicate: {problem.preserved_predicate}")
        lines.append(f"  destroyed data: {problem.destroyed_data}")
        for row in rows[:7]:
            forbidden = ",".join(str(x) for x in row["forbidden"])
            lines.append(
                "  p={:<2} survivors={:<2} forbidden={:<8} factor={:.6f} partial={:.6f}".format(
                    row["p"],
                    row["survivors"],
                    forbidden,
                    row["normalized_factor"],
                    row["partial_product"],
                )
            )
        lines.append(f"  partial product through p={LOCAL_PRIMES[-1]}: {rows[-1]['partial_product']:.6f}")
        lines.append("")
    edges, wins, products = local_problem_tournament()
    n = len(LOCAL_PROBLEMS)
    order = sorted(range(n), key=lambda i: (wins[LOCAL_PROBLEMS[i].name], products[LOCAL_PROBLEMS[i].name]), reverse=True)
    lines.append("Local-problem Tournament Analysis")
    lines.append(f"vertices={n}")
    lines.append(f"score_hist={score_hist(edges, n)}")
    lines.append(f"directed_3cycles={directed_3cycles(edges, n)}")
    lines.append(f"scc_sizes={scc_sizes(edges, n)}")
    lines.append(f"Hamiltonian_paths={hamiltonian_count(edges, n)}")
    lines.append("ranking:")
    for rank, idx in enumerate(order, start=1):
        name = LOCAL_PROBLEMS[idx].name
        lines.append(f"  {rank}. {name} score={wins[name]} partial={products[name]:.6f}")
    lines.append("")


def print_problem_atlas(lines: list[str]) -> None:
    lines.append("Hard-Problem Carrier Atlas")
    lines.append("--------------------------")
    for problem in PROBLEM_CARRIERS:
        lines.append(f"- {problem.name}")
        lines.append(f"  raw scalar: {problem.raw_scalar}")
        lines.append(f"  hidden witness: {problem.hidden_witness}")
        lines.append(f"  retained side channels: {', '.join(problem.side_channels)}")
        lines.append(f"  plausible tournament vertices: {', '.join(problem.candidate_vertices)}")
    lines.append("")
    lines.append("Problem tag similarities")
    for score, a, b, shared in problem_similarity_rows()[:14]:
        lines.append(f"  {score:.3f} {a} <-> {b}: {shared}")
    lines.append("")


def print_route_tournament(lines: list[str]) -> None:
    lines.append("Proof-Route Tournament")
    lines.append("----------------------")
    lines.append("Vertices: proof routes/carriers, not integers or conjectures.")
    lines.append("Observable: majority over predicate retention, computability, cross-problem transfer, proof strength, and lower risk.")
    lines.append("Switch/gauge: A -> B when A wins the majority of criteria.")
    lines.append("Tie Hamiltonian path: route score descending.")
    edges, wins = route_tournament()
    n = len(ROUTES)
    order = sorted(range(n), key=lambda i: wins[ROUTES[i].name], reverse=True)
    lines.append(f"vertices={n}")
    lines.append(f"score_hist={score_hist(edges, n)}")
    lines.append(f"directed_3cycles={directed_3cycles(edges, n)}")
    lines.append(f"scc_sizes={scc_sizes(edges, n)}")
    lines.append(f"Hamiltonian_paths={hamiltonian_count(edges, n)}")
    lines.append("ranking:")
    for rank, idx in enumerate(order, start=1):
        route = ROUTES[idx]
        lines.append(
            f"  {rank}. {route.name} score={wins[route.name]} "
            f"(preserve={route.preserves_predicate}, compute={route.computable_now}, "
            f"transfer={route.cross_problem}, proof={route.proof_strength}, risk={route.risk})"
        )
    lines.append("")


def main() -> str:
    lines: list[str] = []
    lines.append("S651 Number-Theory Tournament Carriers")
    lines.append("======================================")
    lines.append("")
    lines.append("Thesis")
    lines.append("------")
    lines.append("A tournament fits number theory when a pairwise arithmetic observable is")
    lines.append("collapsed to a binary orientation while the side channel that made the")
    lines.append("observable meaningful is kept nearby.")
    lines.append("")
    lines.append("Assumption Challenge")
    lines.append("--------------------")
    lines.append("I did not assume tournament vertices must be primes, runners, or arcs.")
    lines.append("Alternate vertex sets considered: residue classes, gaps, local primes,")
    lines.append("Dirichlet characters, zero modes, Frobenius classes, Selmer local")
    lines.append("conditions, divisor/valuation strata, height packets, circle-method arcs,")
    lines.append("and proof obligations.")
    lines.append("")
    lines.append("Preserved predicates: pairwise lead in prime races; local admissibility for")
    lines.append("linear prime forms; and route quality under retained side channels.")
    lines.append("Destroyed data: exact prime locations, zero ordinates, analytic estimates,")
    lines.append("global primality, Selmer arithmetic, height inequalities, and embeddings.")
    lines.append("")
    print_prime_races(lines)
    print_local_ledgers(lines)
    print_problem_atlas(lines)
    print_route_tournament(lines)
    lines.append("Reading")
    lines.append("-------")
    lines.append("Prime races are the cleanest literal tournament: residue classes compete,")
    lines.append("and edge flips are zero-driven wall crossings.  This is a finite shadow of")
    lines.append("GRH/explicit-formula data, not a replacement for it.")
    lines.append("")
    lines.append("Twin primes, Goldbach, Sophie Germain primes, and prime tuples become local")
    lines.append("obstruction tournaments only after the local product ledger is retained.")
    lines.append("The product says which quotient survives local tests; it does not decide")
    lines.append("the global conjecture because parity barriers and minor-arc/global")
    lines.append("correlations are destroyed by the quotient.")
    lines.append("")
    lines.append("For RH/BSD/abc, the tournament move is moral progress rather than proof:")
    lines.append("choose vertices as zeros, Frobenius traces, local conditions, valuations,")
    lines.append("or proof obligations, then ask which quotient preserves the hard predicate.")
    lines.append("The route tournament says the strongest transferable first move remains")
    lines.append("side-channel jackknife/local ledgers, with raw scalar numerology last.")
    return "\n".join(lines)


if __name__ == "__main__":
    print(main())
