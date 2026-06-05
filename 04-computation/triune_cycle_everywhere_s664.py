#!/usr/bin/env python3
"""S664: search the repo for the sum -> product -> fraction -> sum cycle.

S662 found the representation-face cycle for pi.  S663 showed that, for LRC
n=14, the fraction/carry face repairs a real floor/strict projection collision
left by the additive/product shadow.  This session asks where the same loop
appears elsewhere.

Interpretation used here:

    sum      = additive aggregation, packets, traces, counts, moments;
    product  = local factors, gcd/norm shells, independent components;
    fraction = boundary state, owner/lift/branch, deletion or recursion data.

The cycle is not just "all three words appear".  The useful pattern is:

    sum -> product    additive packets expose/localize factor data;
    product -> fraction
                      local factors are not enough without an owner/lift;
    fraction -> sum   recursive boundary state re-expands into packet sums.

Tournament Analysis:
  * Face tournament vertices are the four carrier faces.  Pairwise observable
    is domain-weighted dependency support.  The switch is weighted majority.
    The tie Hamiltonian path is sum -> product -> fraction -> raw_scalar.
  * Domain tournament vertices are problem lanes.  Pairwise observable is
    (cycle_strength, exact_collision_evidence, finite_lab, transfer_leverage).
    The tie path is the curated domain order below.

Assumption challenge:
  Candidate vertices considered: faces, domains, files, runners, residues,
  wall pairs, gcd shells, carry words, deck cards, deleted vertices, points,
  prime pairs, models/generics, and proof obligations.  The script uses faces
  and domains because the preserved predicate is "which carrier is needed to
  repair scalar/projection leakage."  It intentionally destroys many local
  labels during the repo scan; curated entries restore the known owner state.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCAN_DIRS = ("00-navigation", "04-computation", "05-knowledge", "07-reflections")
TEXT_SUFFIXES = {".md", ".py", ".out", ".txt", ".tsv"}
MAX_SCAN_BYTES = 1_500_000
SELF_SCAN_EXCLUDE = {"05-knowledge/results/triune_cycle_everywhere_s664.out"}
EXCLUDED_PARTS = {
    ".git",
    ".lake",
    "__pycache__",
    ".mypy_cache",
    ".pytest_cache",
    "node_modules",
    "site-packages",
}

FACES = ("sum", "product", "fraction", "raw_scalar")
CYCLE_EDGES = (("sum", "product"), ("product", "fraction"), ("fraction", "sum"))

FACE_TERMS = {
    "sum": (
        "sum",
        "additive",
        "packet",
        "pair-sum",
        "pair sum",
        "wall",
        "trace",
        "power sum",
        "moment",
        "sigma(",
        "goldbach",
        "lemoine",
        "distance support",
        "burnside",
        "count",
    ),
    "product": (
        "product",
        "factor",
        "gcd",
        "shell",
        "norm",
        "local",
        "singular series",
        "euler product",
        "component",
        "concurrency",
        "prime power",
        "scissors",
        "independence polynomial",
    ),
    "fraction": (
        "fraction",
        "continued",
        "continuant",
        "carry",
        "owner",
        "boundary",
        "deletion",
        "branch",
        "lift",
        "fiber",
        "generic",
        "recursion",
        "quotient",
        "marked",
    ),
}


@dataclass(frozen=True)
class Domain:
    name: str
    query_terms: tuple[str, ...]
    sum_face: str
    product_face: str
    fraction_face: str
    cycle_edges: dict[tuple[str, str], int]
    exact_collision: int
    finite_lab: int
    transfer_leverage: int
    evidence: str
    projection_repair: str
    next_probe: str

    @property
    def cycle_strength(self) -> int:
        return sum(self.cycle_edges.get(edge, 0) for edge in CYCLE_EDGES)

    @property
    def complete_cycle(self) -> bool:
        return all(self.cycle_edges.get(edge, 0) > 0 for edge in CYCLE_EDGES)


DOMAINS = [
    Domain(
        "LRC14 carry-owner",
        ("LRC", "n=14", "C=27", "carry", "odd-wall", "Res_27"),
        "odd-wall packets and active pair sums",
        "C=27 gcd/Pillai shell and local obstruction mass",
        "carry word k in v=r+27k plus owner/continuant state",
        {("sum", "product"): 5, ("product", "fraction"): 5, ("fraction", "sum"): 5},
        5,
        5,
        5,
        "S663: 3 mixed sum/product groups become 0 mixed groups with carry state.",
        "Res_27 additive/product shadow -> carry-continuant full key.",
        "Attach HYP-2165 owner routes to the carry key and test deeper carries.",
    ),
    Domain(
        "Pi representations",
        ("Wallis", "Brouncker", "cotangent", "Leibniz", "Basel", "Machin"),
        "Leibniz, Machin, Basel, cotangent packet sums",
        "Wallis/sine factors, local zeros, Euler-product attitude",
        "Brouncker/cotangent continuants and convergent owners",
        {("sum", "product"): 5, ("product", "fraction"): 4, ("fraction", "sum"): 5},
        4,
        5,
        4,
        "S662 face tournament gives sum -> product -> fraction -> sum exactly.",
        "Raw decimal scalar -> representation trinity.",
        "Use cotangent-pole residues as a finite analog of LRC wall ownership.",
    ),
    Domain(
        "OCF / H(T)",
        ("OCF", "H(T)", "independence polynomial", "Hamiltonian paths", "Omega"),
        "odd-cycle packet coefficients and forward-polynomial sums",
        "I(Omega,2), strong-component multiplication, local forbidden H packets",
        "deletion/contraction/substitution boundary state for Hamiltonian paths",
        {("sum", "product"): 5, ("product", "fraction"): 4, ("fraction", "sum"): 4},
        4,
        4,
        5,
        "H(T)=I(Omega(T),2) packages odd-cycle sums as a partition function.",
        "Component/product OCF shadow -> deletion/substitution continuant.",
        "Build a macro-word continuant toy and compare with Hamiltonian-path DP.",
    ),
    Domain(
        "Tournament decks",
        ("tournament deck", "full deck", "deleted score", "paired", "HYP-2236"),
        "deleted-card loss sums, H-loss and c3-loss ledgers",
        "deck product/scissors packets and component cards",
        "paired card-owner derivative: (card, deleted score)",
        {("sum", "product"): 4, ("product", "fraction"): 5, ("fraction", "sum"): 4},
        5,
        5,
        4,
        "S660: full decks collide through n=6; paired deleted score resolves them.",
        "Unpaired deck/product shadow -> paired card-boundary state.",
        "Push canonical augmentation toward n=7 and test paired derivative minimality.",
    ),
    Domain(
        "Unit distance",
        ("unit distance", "Moser", "Eisenstein", "spine", "direction", "frontier"),
        "unit-edge spine, edge-packet counts, tile flips",
        "direction support, Eisenstein norm shells, unit-vector factors",
        "point-deletion frontier/ear owner and recursive attachment state",
        {("sum", "product"): 4, ("product", "fraction"): 4, ("fraction", "sum"): 4},
        3,
        4,
        5,
        "S622/S626/S628: direction/norm carriers need frontier ownership.",
        "Direction product shadow -> point deletion/ear owner.",
        "For n=21/22 cores, store edge sums + direction products + deletion owners.",
    ),
    Domain(
        "Finite-field Kakeya/Falconer",
        ("Kakeya", "Falconer", "finite-field", "pinned", "concurrency"),
        "distance and pinned-distance support packets",
        "line-direction products, concurrency factors, quadratic distance shells",
        "line-choice recursion and owner of pins",
        {("sum", "product"): 4, ("product", "fraction"): 4, ("fraction", "sum"): 3},
        3,
        5,
        4,
        "S659: full distance support is scalar; pinned/concurrency labels vary.",
        "Direction/distance product shadow -> pin owner state.",
        "Enumerate full-distance twins and ask which owner labels split them.",
    ),
    Domain(
        "Goldbach / Lemoine pair plane",
        ("Goldbach", "Lemoine", "prime pair", "singular series", "twin primes"),
        "E=p+q and O=p+2q additive pair fibers",
        "Hardy-Littlewood singular-series local obstruction products",
        "ordered-pair reconstruction q=O-E, p=2E-O",
        {("sum", "product"): 4, ("product", "fraction"): 3, ("fraction", "sum"): 4},
        3,
        4,
        4,
        "Pair-plane threads: same pair has additive columns and local sieves.",
        "Sum/local sieve shadow -> ordered reconstruction branch.",
        "Track when the same prime pair serves Goldbach and Lemoine columns.",
    ),
    Domain(
        "Pi/e trace-norm",
        ("pi/e", "trace-norm", "T^2-S*T+P", "Schanuel", "branch sheet"),
        "trace S=e+pi and Newton power sums",
        "norm P=e*pi and discriminant/local product shadows",
        "branch sheet of T^2-S*T+P, discriminant sign, root owner",
        {("sum", "product"): 4, ("product", "fraction"): 5, ("fraction", "sum"): 4},
        4,
        4,
        4,
        "S635/S636: any two of trace/norm/branch reconstruct the pair.",
        "Trace/norm shadow -> branch sheet.",
        "Run finite-field Vieta labs as toy models for branch-owner leakage.",
    ),
    Domain(
        "Perfect / aliquot carrier",
        ("perfect number", "aliquot", "divisor", "sigma", "fixed point"),
        "sigma(n) divisor-sum packet",
        "Euler product over prime powers in sigma(n)",
        "aliquot iteration boundary/orbit owner",
        {("sum", "product"): 5, ("product", "fraction"): 4, ("fraction", "sum"): 5},
        3,
        5,
        3,
        "Perfect numbers are fixed points of the divisor-sum/aliquot carrier.",
        "Divisor product formula -> aliquot orbit state.",
        "Use aliquot preperiod/cycle owners as fraction faces for sigma shadows.",
    ),
    Domain(
        "A000568 / marked LRC quotient",
        ("A000568", "marked", "rooted", "observer", "Burnside"),
        "Burnside orbit sums and chamber-count packets",
        "cycle-type products and automorphism/local stabilizer factors",
        "rooted/marked observer fiber and endpoint boundary state",
        {("sum", "product"): 4, ("product", "fraction"): 3, ("fraction", "sum"): 4},
        3,
        4,
        4,
        "Raw A000568 classes are too coarse; marked observer fibers repair them.",
        "Unmarked orbit/product count -> rooted/marked boundary fiber.",
        "Treat every quotient count as needing its marked fraction face.",
    ),
    Domain(
        "CH / forcing",
        ("Continuum Hypothesis", "forcing", "generic", "constructible", "cardinal-shadow"),
        "cardinal/equinumerosity sums and continuum-size shadows",
        "model-local consistency factors and product of independent choices",
        "generic extension boundary state / forcing name",
        {("sum", "product"): 3, ("product", "fraction"): 3, ("fraction", "sum"): 3},
        2,
        2,
        4,
        "S656: scalar cardinal shadow needs model/generic side channel.",
        "Cardinal/product consistency shadow -> generic boundary state.",
        "Use forcing twins as the abstract version of carry-lift leakage.",
    ),
    Domain(
        "Cauldron / Schur games",
        ("cauldron", "Schur", "finite-sums", "block-turn", "schedule"),
        "active forbidden-sum hyperedges",
        "color/resource products and residue-channel survival factors",
        "turn-schedule boundary word / minimax state",
        {("sum", "product"): 4, ("product", "fraction"): 3, ("fraction", "sum"): 3},
        3,
        5,
        3,
        "S618-S621: schedule word acts as a boundary carrier for additive rules.",
        "Sum-free/product resource shadow -> schedule/minimax owner state.",
        "Compare parity vs two-block games with the trinity ledger explicitly.",
    ),
    Domain(
        "Heegner / prime polynomials",
        ("Heegner", "prime generating", "unique factorization", "class number"),
        "polynomial value sums/traces such as n^2+n+A",
        "unique-factorization/product obstruction in quadratic orders",
        "reduced-form/continued-fraction class-cycle boundary",
        {("sum", "product"): 3, ("product", "fraction"): 5, ("fraction", "sum"): 3},
        2,
        3,
        4,
        "Heegner numbers say product factorization controls polynomial primes.",
        "UFD product law -> reduced-form/class-cycle owner.",
        "Build a small discriminant/form-cycle atlas for prime-generating runs.",
    ),
]


def lower_terms(terms: tuple[str, ...]) -> tuple[str, ...]:
    return tuple(term.lower() for term in terms)


def iter_text_files() -> list[tuple[Path, str]]:
    rows: list[tuple[Path, str]] = []
    for dirname in SCAN_DIRS:
        base = ROOT / dirname
        if not base.exists():
            continue
        for path in base.rglob("*"):
            if any(part in EXCLUDED_PARTS for part in path.parts):
                continue
            rel = path.relative_to(ROOT)
            if str(rel) in SELF_SCAN_EXCLUDE:
                continue
            if not path.is_file() or path.suffix not in TEXT_SUFFIXES:
                continue
            try:
                if path.stat().st_size > MAX_SCAN_BYTES:
                    continue
                text = path.read_text(encoding="utf-8", errors="ignore")
            except OSError:
                continue
            rows.append((path, text))
    return rows


def term_count(text_lower: str, terms: tuple[str, ...]) -> int:
    return sum(text_lower.count(term) for term in terms)


def repo_scan(texts: list[tuple[Path, str]]) -> dict[str, object]:
    face_terms = {face: lower_terms(terms) for face, terms in FACE_TERMS.items()}
    top_files = []
    face_file_counts = Counter()
    all_three = 0
    for path, text in texts:
        low = text.lower()
        counts = {face: term_count(low, terms) for face, terms in face_terms.items()}
        present = [face for face, count in counts.items() if count]
        for face in present:
            face_file_counts[face] += 1
        if all(counts[face] for face in ("sum", "product", "fraction")):
            all_three += 1
            total = counts["sum"] + counts["product"] + counts["fraction"]
            rel = path.relative_to(ROOT)
            top_files.append((total, counts, str(rel)))
    top_files.sort(key=lambda row: (-row[0], row[2]))
    return {
        "files_scanned": len(texts),
        "face_file_counts": dict(sorted(face_file_counts.items())),
        "files_with_all_three": all_three,
        "top_files": top_files[:14],
    }


def domain_hits(texts: list[tuple[Path, str]], domain: Domain) -> dict[str, object]:
    qterms = lower_terms(domain.query_terms)
    face_terms = {face: lower_terms(terms) for face, terms in FACE_TERMS.items()}
    file_hits = 0
    face_counts = Counter()
    sample_files: list[tuple[int, str]] = []
    for path, text in texts:
        low = text.lower()
        qcount = term_count(low, qterms)
        if not qcount:
            continue
        file_hits += 1
        rel = str(path.relative_to(ROOT))
        sample_files.append((qcount, rel))
        for face, terms in face_terms.items():
            face_counts[face] += term_count(low, terms)
    sample_files.sort(key=lambda row: (-row[0], row[1]))
    return {
        "files": file_hits,
        "face_counts": dict(sorted(face_counts.items())),
        "sample_files": [name for _, name in sample_files[:4]],
    }


def tournament_fingerprint(
    vertices: tuple[str, ...],
    weights: dict[tuple[str, str], int],
    tie_order: tuple[str, ...],
) -> dict[str, object]:
    n = len(vertices)
    idx = {v: i for i, v in enumerate(vertices)}
    tie_pos = {v: i for i, v in enumerate(tie_order)}
    adj = [[0] * n for _ in range(n)]
    out = [0] * n
    edges = []
    for a, b in combinations(vertices, 2):
        ab = weights.get((a, b), 0)
        ba = weights.get((b, a), 0)
        if ab > ba:
            winner, loser = a, b
        elif ba > ab:
            winner, loser = b, a
        else:
            winner, loser = (a, b) if tie_pos[a] < tie_pos[b] else (b, a)
        adj[idx[winner]][idx[loser]] = 1
        out[idx[winner]] += 1
        edges.append((winner, loser, weights.get((winner, loser), 0), weights.get((loser, winner), 0)))

    c3 = 0
    cycle_triples = []
    for a, b, c in combinations(range(n), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            c3 += 1
            cycle_triples.append((vertices[a], vertices[b], vertices[c]))

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
                comp.append(vertices[v])
                for w in radj[v]:
                    if not seen[w]:
                        seen[w] = True
                        q.append(w)
            comps.append(sorted(comp))
        return comps

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
        "outscores": {vertices[i]: out[i] for i in range(n)},
        "score_hist": dict(sorted(Counter(out).items())),
        "directed_3cycles": c3,
        "cycle_triples": cycle_triples,
        "sccs": sccs(),
        "hamiltonian_paths": sum(dp[-1]),
        "edges": edges,
    }


def face_tournament(domains: list[Domain]) -> tuple[dict[tuple[str, str], int], dict[str, object]]:
    weights: dict[tuple[str, str], int] = defaultdict(int)
    for domain in domains:
        for edge, strength in domain.cycle_edges.items():
            weights[edge] += strength
    for face in ("sum", "product", "fraction"):
        weights[(face, "raw_scalar")] += 3 * len(domains)
    return dict(weights), tournament_fingerprint(FACES, weights, ("sum", "product", "fraction", "raw_scalar"))


def domain_tournament(domains: list[Domain]) -> dict[str, object]:
    weights: dict[tuple[str, str], int] = defaultdict(int)
    metrics = {
        d.name: (d.cycle_strength, d.exact_collision, d.finite_lab, d.transfer_leverage)
        for d in domains
    }
    order = tuple(d.name for d in domains)
    for a, b in combinations(domains, 2):
        wins_a = sum(x > y for x, y in zip(metrics[a.name], metrics[b.name]))
        wins_b = sum(x < y for x, y in zip(metrics[a.name], metrics[b.name]))
        if wins_a > wins_b:
            weights[(a.name, b.name)] = 1
        elif wins_b > wins_a:
            weights[(b.name, a.name)] = 1
        else:
            winner, loser = (a.name, b.name) if order.index(a.name) < order.index(b.name) else (b.name, a.name)
            weights[(winner, loser)] = 1
    fp = tournament_fingerprint(tuple(d.name for d in domains), weights, order)
    fp["metrics"] = metrics
    fp["top_order"] = [
        name for name, _score in sorted(
            fp["outscores"].items(), key=lambda item: (-item[1], order.index(item[0]))
        )
    ]
    return fp


def primes_upto(n: int) -> list[int]:
    sieve = [True] * (n + 1)
    sieve[0:2] = [False, False]
    for p in range(2, int(n**0.5) + 1):
        if sieve[p]:
            for q in range(p * p, n + 1, p):
                sieve[q] = False
    return [p for p in range(n + 1) if sieve[p]]


def goldbach_lemoine_lab(limit: int = 160) -> dict[str, object]:
    primes = primes_upto(limit)
    prime_set = set(primes)
    pairs = [(p, q) for p in primes for q in primes if p <= q and p + q <= limit]
    shared = []
    diagonal = []
    for p, q in pairs:
        even = p + q
        odd = p + 2 * q
        if odd <= limit + 80 and odd % 2 == 1:
            rec_p = 2 * even - odd
            rec_q = odd - even
            if rec_p == p and rec_q == q:
                shared.append((p, q, even, odd))
                if p == q:
                    diagonal.append((p, even, odd))
    even_counts = Counter(p + q for p, q in pairs)
    lemoine_counts = Counter(p + 2 * q for p in primes for q in primes if p + 2 * q <= limit + 80)
    return {
        "prime_count": len(primes),
        "goldbach_columns": len(even_counts),
        "max_goldbach_fiber": max(even_counts.values()),
        "lemoine_columns": len(lemoine_counts),
        "max_lemoine_fiber": max(lemoine_counts.values()),
        "shared_pair_count": len(shared),
        "diagonal_examples": diagonal[:8],
        "shared_examples": shared[:8],
    }


def vieta_branch_lab(p: int = 17) -> dict[str, object]:
    pairs = [(x, y) for x in range(p) for y in range(p)]
    sum_buckets: defaultdict[int, list[tuple[int, int]]] = defaultdict(list)
    product_buckets: defaultdict[int, list[tuple[int, int]]] = defaultdict(list)
    sp_buckets: defaultdict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    branch_buckets: defaultdict[tuple[int, int, int], list[tuple[int, int]]] = defaultdict(list)
    for x, y in pairs:
        s = (x + y) % p
        prod = (x * y) % p
        branch = (x - y) % p
        sum_buckets[s].append((x, y))
        product_buckets[prod].append((x, y))
        sp_buckets[(s, prod)].append((x, y))
        branch_buckets[(s, prod, branch)].append((x, y))
    return {
        "field": p,
        "sum_only_max_fiber": max(len(v) for v in sum_buckets.values()),
        "product_only_max_fiber": max(len(v) for v in product_buckets.values()),
        "sum_product_max_fiber": max(len(v) for v in sp_buckets.values()),
        "sum_product_branch_max_fiber": max(len(v) for v in branch_buckets.values()),
        "sum_product_collision_example": next(
            (key, vals) for key, vals in sp_buckets.items() if len(vals) == 2
        ),
    }


def factor(n: int) -> dict[int, int]:
    out: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def sigma_product(n: int) -> int:
    prod = 1
    for p, a in factor(n).items():
        prod *= (p ** (a + 1) - 1) // (p - 1)
    return prod


def sigma_sum(n: int) -> int:
    return sum(d for d in range(1, n + 1) if n % d == 0)


def aliquot_lab(limit: int = 300) -> dict[str, object]:
    perfect = []
    product_formula_mismatches = []
    two_cycles = []
    aliquot = {n: sigma_product(n) - n for n in range(1, limit + 1)}
    for n in range(1, limit + 1):
        if sigma_sum(n) != sigma_product(n):
            product_formula_mismatches.append(n)
        if aliquot[n] == n:
            perfect.append(n)
    for n in range(2, limit + 1):
        m = aliquot[n]
        if n < m <= limit and aliquot.get(m) == n:
            two_cycles.append((n, m))
    return {
        "limit": limit,
        "sigma_product_mismatches": product_formula_mismatches[:5],
        "perfect_fixed_points": perfect,
        "aliquot_two_cycles": two_cycles[:8],
    }


def print_wrapped(prefix: str, text: str, width: int = 96) -> None:
    words = text.split()
    line = prefix
    for word in words:
        if len(line) + len(word) + 1 > width:
            print(line.rstrip())
            line = " " * len(prefix)
        line += word + " "
    print(line.rstrip())


def main() -> None:
    texts = iter_text_files()
    scan = repo_scan(texts)
    domain_scans = {domain.name: domain_hits(texts, domain) for domain in DOMAINS}
    face_weights, face_fp = face_tournament(DOMAINS)
    domain_fp = domain_tournament(DOMAINS)

    print("=" * 78)
    print("S664 triune-cycle-everywhere atlas")
    print("=" * 78)
    print()
    print("Cycle criterion")
    print("  sum      = additive aggregation: packets, traces, walls, counts, moments")
    print("  product  = local factorization: gcd/norm shells, components, sieves")
    print("  fraction = boundary memory: carry, owner, branch, lift, deletion, recursion")
    print("  useful cycle = sum -> product -> fraction -> sum")
    print()

    print("A. Repo scan")
    print(f"  files_scanned={scan['files_scanned']}")
    print(f"  face_file_counts={scan['face_file_counts']}")
    print(f"  files_with_all_three_faces={scan['files_with_all_three']}")
    print("  top all-three files:")
    for total, counts, rel in scan["top_files"]:
        print(f"    score={total:4d} counts={counts} file={rel}")
    print()

    print("B. Domain atlas")
    header = (
        f"{'domain':<34} {'repo_files':>10} {'cycle':>6} "
        f"{'exact':>5} {'lab':>3} {'xfer':>4} edges"
    )
    print(header)
    for domain in DOMAINS:
        ds = domain_scans[domain.name]
        edges = ",".join(
            f"{a[0]}>{b[0]}:{domain.cycle_edges.get((a, b), 0)}"
            for a, b in CYCLE_EDGES
        )
        print(
            f"{domain.name:<34} {ds['files']:10d} {domain.cycle_strength:6d} "
            f"{domain.exact_collision:5d} {domain.finite_lab:3d} "
            f"{domain.transfer_leverage:4d} {edges}"
        )
    print()

    print("C. Strongest reframings")
    for domain in sorted(DOMAINS, key=lambda d: (-d.cycle_strength, -d.exact_collision, d.name))[:8]:
        print(f"  {domain.name}")
        print_wrapped("    sum:      ", domain.sum_face)
        print_wrapped("    product:  ", domain.product_face)
        print_wrapped("    fraction: ", domain.fraction_face)
        print_wrapped("    repair:   ", domain.projection_repair)
        print_wrapped("    next:     ", domain.next_probe)
    print()

    print("D. Face-dependency Tournament Analysis")
    print("  vertices=sum, product, fraction, raw_scalar")
    print("  observable=domain-weighted dependency support")
    print("  switch=weighted majority; tie path=sum -> product -> fraction -> raw_scalar")
    print(f"  weights={dict(sorted(face_weights.items()))}")
    print(f"  outscores={face_fp['outscores']}")
    print(f"  score_hist={face_fp['score_hist']}")
    print(f"  directed_3cycles={face_fp['directed_3cycles']} {face_fp['cycle_triples']}")
    print(f"  sccs={face_fp['sccs']}")
    print(f"  hamiltonian_paths={face_fp['hamiltonian_paths']}")
    print("  edges:")
    for winner, loser, w_for, w_against in face_fp["edges"]:
        print(f"    {winner} -> {loser}  support={w_for}:{w_against}")
    print()

    print("E. Domain-route Tournament Analysis")
    print("  vertices=problem lanes")
    print("  observable=(cycle_strength, exact_collision_evidence, finite_lab, transfer_leverage)")
    print(f"  score_hist={domain_fp['score_hist']}")
    print(f"  directed_3cycles={domain_fp['directed_3cycles']}")
    print(f"  sccs={domain_fp['sccs']}")
    print(f"  hamiltonian_paths={domain_fp['hamiltonian_paths']}")
    print("  top_order:")
    for name in domain_fp["top_order"]:
        print(f"    {name}: metrics={domain_fp['metrics'][name]}")
    print()

    print("F. Small side labs")
    gl = goldbach_lemoine_lab()
    print("  Goldbach/Lemoine pair plane")
    print(f"    prime_count={gl['prime_count']}")
    print(f"    goldbach_columns={gl['goldbach_columns']} max_fiber={gl['max_goldbach_fiber']}")
    print(f"    lemoine_columns={gl['lemoine_columns']} max_fiber={gl['max_lemoine_fiber']}")
    print(f"    shared_pair_count={gl['shared_pair_count']}")
    print(f"    diagonal_examples(p,2p,3p)={gl['diagonal_examples']}")
    print(f"    shared_examples(p,q,E,O)={gl['shared_examples']}")
    vf = vieta_branch_lab()
    print("  Vieta branch lab over F_17")
    print(f"    sum_only_max_fiber={vf['sum_only_max_fiber']}")
    print(f"    product_only_max_fiber={vf['product_only_max_fiber']}")
    print(f"    sum_product_max_fiber={vf['sum_product_max_fiber']}")
    print(f"    sum_product_branch_max_fiber={vf['sum_product_branch_max_fiber']}")
    print(f"    collision_example((S,P),pairs)={vf['sum_product_collision_example']}")
    al = aliquot_lab()
    print("  Perfect/aliquot carrier")
    print(f"    sigma_product_mismatches={al['sigma_product_mismatches']}")
    print(f"    perfect_fixed_points={al['perfect_fixed_points']}")
    print(f"    aliquot_two_cycles={al['aliquot_two_cycles']}")
    print()

    print("G. Hypothesis upgrades")
    print("  1. The triune cycle is a projection-repair grammar, not a numerology.")
    print("  2. The strongest exact instances are LRC14 and tournament decks: both")
    print("     have visible collisions that disappear only after retaining owner state.")
    print("  3. OCF/H(T) is the algebraic hinge: additive odd-cycle packets are already")
    print("     a product-valued partition function, but deletion/substitution boundaries")
    print("     are the missing fraction face.")
    print("  4. Unit distance should be attacked with the same impairment test:")
    print("     edge-spine sums + direction/norm products + point-deletion owners.")
    print("  5. In number theory, Vieta/Goldbach-Lemoine/perfect-number carriers all say")
    print("     the branch/orbit owner is the fraction face that turns local products")
    print("     back into additive statements.")


if __name__ == "__main__":
    main()
