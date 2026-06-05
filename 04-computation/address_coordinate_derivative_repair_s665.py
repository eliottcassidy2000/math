#!/usr/bin/env python3
"""S665: scalar/product collision -> address coordinate -> derivative sum.

S664 found the broad triune loop:

    sum -> product -> fraction -> sum.

This script sharpens that into the operational search pattern requested by the
user:

    find a scalar/product collision,
    attach the missing fraction/address coordinate,
    unroll it back into a derivative sum.

An address coordinate is the fraction face made explicit: carry word, deleted
card owner, branch sheet, Weyl m-function, continued-fraction coefficient,
activity mark, generic/model state, or point-deletion owner.  A derivative sum
is the result of differentiating the quotient in the right category: log
derivative, deletion derivative, spectral derivative, moment/J-fraction
expansion, owner-loss deck, or recursive parameter sum.

Tournament Analysis:
  Vertices are candidate repair lanes.  Pairwise observable is the vector
  (collision_evidence, address_clarity, derivative_sum_clarity,
   finite_actionability, repo_leverage, source_confidence).  The switch is
  majority; ties follow the curated source order.

Assumption challenge:
  Candidate vertices considered included formulas, domains, files, invariants,
  vertices, edges, deleted cards, residues, runners, primes, roots, eigenvalues,
  moments, bases, models, and proof obligations.  This script uses repair lanes
  because the preserved predicate is "which missing coordinate repairs a
  quotient collision and yields a derivative ledger?"
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Lane:
    name: str
    family: str
    source_kind: str
    source: str
    scalar_product_collision: str
    address_coordinate: str
    derivative_sum: str
    next_probe: str
    collision_evidence: int
    address_clarity: int
    derivative_clarity: int
    finite_actionability: int
    repo_leverage: int
    source_confidence: int

    @property
    def score_vector(self) -> tuple[int, int, int, int, int, int]:
        return (
            self.collision_evidence,
            self.address_clarity,
            self.derivative_clarity,
            self.finite_actionability,
            self.repo_leverage,
            self.source_confidence,
        )

    @property
    def total(self) -> int:
        return sum(self.score_vector)

    @property
    def exactish(self) -> bool:
        return self.collision_evidence >= 4 and self.address_clarity >= 4


LANES = [
    Lane(
        "LRC14 Res27 carry",
        "repo/LRC",
        "repo exact",
        "S663/HYP-2239",
        "same Res27 additive/product shadow mixes floor and strict rows",
        "carry word k in v=r+27k plus owner/continuant state",
        "carry-deletion derivative deck: M, wall pairs, gcd mass, apex and owner losses",
        "attach HYP-2165 owner routes and test deeper carry perturbations",
        5,
        5,
        5,
        5,
        5,
        5,
    ),
    Lane(
        "Tournament deck paired owner",
        "repo/tournaments",
        "repo exact",
        "S660/HYP-2236",
        "full decks and global scalar repairs collide through n=6",
        "paired card address (card isomorphism type, deleted vertex outdegree)",
        "deleted-score, H-loss, and c3-loss sums paired to each card",
        "extend paired derivative minimality to n=7 with canonical augmentation",
        5,
        5,
        5,
        5,
        5,
        5,
    ),
    Lane(
        "Graph reconstruction Kocay matrix",
        "web/graphs",
        "external theorem framework",
        "arXiv:1301.4121; arXiv:2601.00620",
        "many graphs may share coarse deck-derived product/scalar shadows",
        "covering-number or vertex-pair parameter address rows",
        "Kocay-style linear constraints and recursive pair-parameter equations",
        "port paired derivative deck tests from tournaments to ordinary graphs",
        4,
        5,
        5,
        4,
        4,
        5,
    ),
    Lane(
        "Riemann/Weil explicit formula",
        "web/number theory",
        "classical external",
        "explicit formulae for L-functions",
        "Euler product alone hides zero/gamma/infinity address data",
        "zeros, gamma factor, and test-function/Fourier address",
        "logarithmic derivative turns products into prime-power and zero sums",
        "use as the model for product-ledger -> derivative-sum transformations",
        4,
        5,
        5,
        3,
        4,
        4,
    ),
    Lane(
        "Determinant adjugate derivative",
        "web/linear algebra",
        "classical external",
        "Jacobi determinant formula",
        "same determinant/product of eigenvalues can hide directional response",
        "adjugate/inverse row-column address for the perturbation",
        "Jacobi trace/cofactor formula gives the derivative sum",
        "use adjugate-style owner cofactors for tournament path determinant toys",
        4,
        5,
        5,
        5,
        3,
        5,
    ),
    Lane(
        "Ihara/Selberg graph zeta",
        "web/graphs",
        "external theorem family",
        "Bass-Ihara determinant; Selberg trace formula",
        "determinant/product zeta can collide at spectral level",
        "primitive closed geodesic or scattering/geodesic-length address",
        "log derivative/trace formula expands product into closed-orbit sums",
        "compute Ihara/Bartholdi zeta for small tournament-derived graphs",
        4,
        4,
        5,
        4,
        4,
        4,
    ),
    Lane(
        "Sparse phase retrieval Prony",
        "web/signal",
        "external inverse problem",
        "sparse phase retrieval via Prony",
        "Fourier magnitude/product spectrum forgets phase and pair assignment",
        "phase/support-difference owner from Prony reconstruction",
        "inverse Fourier/Prony step reconstructs spike coefficient sums",
        "model LRC Res27 wall packets as a sparse phase-retrieval problem",
        4,
        5,
        5,
        4,
        3,
        5,
    ),
    Lane(
        "Turnpike/beltway homometry",
        "web/inverse geometry",
        "external inverse problem",
        "turnpike and beltway reconstruction",
        "unassigned pairwise distances collide between point sets",
        "which distance belongs to which pair/end owner",
        "placed points re-expand into coordinate and interval sums",
        "use homometric point-set collisions as unit-distance owner stress tests",
        4,
        4,
        4,
        5,
        4,
        5,
    ),
    Lane(
        "Stieltjes moment J-fraction",
        "web/moments",
        "external theorem family",
        "SIAM structured matrices; arXiv:1907.05204",
        "moment/Hankel scalar data can be indeterminate or too coarse",
        "continued-fraction coefficients and Hankel determinant address",
        "J-fraction/orthogonal-polynomial recurrence unrolls into moment sums",
        "build a finite moment-collision toy beside OCF independence packets",
        4,
        5,
        5,
        4,
        3,
        4,
    ),
    Lane(
        "Inverse spectral Weyl m-function",
        "web/spectral",
        "external theorem",
        "Barry Simon inverse spectral theory",
        "spectrum/product data alone needs norming or spectral measure address",
        "Weyl m-function or spectral measure",
        "m-function expansion and A-equation unroll potential into derivative/integral sums",
        "use m-function as archetype for LRC owner address, not just analogy",
        4,
        5,
        5,
        3,
        3,
        5,
    ),
    Lane(
        "Tutte deletion-contraction activities",
        "web/combinatorics",
        "external theorem/tool",
        "Tutte polynomial references",
        "many graph invariants are scalar specializations of one product-like polynomial",
        "edge order plus internal/external activity address",
        "deletion-contraction derivative sums recover spanning tree/activity ledgers",
        "turn tournament H deletion-contraction into an activity-address toy",
        3,
        4,
        5,
        5,
        4,
        4,
    ),
    Lane(
        "Matroid Hodge characteristic polynomial",
        "web/matroids",
        "external theorem",
        "Adiprasito-Huh-Katz",
        "characteristic-polynomial coefficients are scalar shadows of matroid geometry",
        "Chow ring degree/address and Lefschetz operators",
        "Hodge-Riemann inequalities unroll coefficients into intersection/derivative constraints",
        "look for a Chow-ring analogue of OCF side-channel positivity",
        3,
        4,
        4,
        2,
        3,
        5,
    ),
    Lane(
        "Pi/e Vieta branch sheet",
        "repo/number theory",
        "repo exact toy",
        "S635/S636/HYP-2211/HYP-2212",
        "trace/norm pair leaves root order or branch sheet unresolved",
        "discriminant D=e-pi or finite-field branch x-y",
        "Newton power sums and branch derivative split the two sheets",
        "reuse finite-field branch lab as a generic collision detector",
        4,
        5,
        4,
        5,
        4,
        4,
    ),
    Lane(
        "OCF noncommutative address",
        "repo/tournaments",
        "repo/web theorem lead",
        "T121/T131/T1667; Mitrovic",
        "commutative OCF H=I(Omega,2) forgets vertex/set-partition addresses",
        "noncommuting variables and set partitions",
        "deletion-contraction in the richer W_X unrolls to OCF packet sums",
        "test transfer-matrix symmetry by projecting NC identities downward",
        4,
        5,
        5,
        3,
        5,
        4,
    ),
    Lane(
        "Unit distance deletion owner",
        "repo/unit distance",
        "repo finite evidence",
        "S622/S626/S628/HYP-2203",
        "edge count and direction/norm product shadows do not encode construction ownership",
        "point-deletion frontier/ear owner",
        "frontier-gain loss sums under deleting or adding a point",
        "run n=21/22 owner-preserving impairment lab",
        3,
        4,
        4,
        4,
        5,
        4,
    ),
    Lane(
        "Finite-field Kakeya pinned owner",
        "repo/finite field",
        "repo exact finite",
        "S659/HYP-2235",
        "full distance support saturates while union size and pinned minima vary",
        "pin owner, line-choice recursion, concurrency label",
        "pinned-distance and concurrency derivative sums",
        "classify full-distance-support twins by owner labels",
        4,
        4,
        4,
        5,
        4,
        4,
    ),
    Lane(
        "Union-closed frequency pressure",
        "repo/set systems",
        "repo exact finite",
        "S658/HYP-2234",
        "frequency vectors collide at m=4",
        "set pressure address sum_B |A union B|",
        "pressure-loss derivative over added/deleted sets",
        "classify frequency collision buckets by pressure derivative profile",
        4,
        4,
        4,
        5,
        3,
        4,
    ),
    Lane(
        "Goldbach/Lemoine pair plane",
        "repo/number theory",
        "repo exact identity",
        "S642/S643/HYP-2218/HYP-2219",
        "even E=p+q and odd O=p+2q are separate scalar columns",
        "ordered pair address q=O-E, p=2E-O",
        "parent-rich odd/even companion sums and swap-reflection derivative",
        "rank shared prime-pair columns by local sieve plus reconstruction address",
        4,
        5,
        4,
        5,
        4,
        4,
    ),
    Lane(
        "Perfect/aliquot orbit owner",
        "repo/number theory",
        "repo finite exact",
        "S644/S645/HYP-2220",
        "sigma product/abundancy scalar proximity does not decide fixedness",
        "aliquot orbit owner: fixed, amicable, sociable, preperiod",
        "aliquot derivative s(n)-n and orbit-loss sums",
        "attach orbit owner to triangular pair-count and divisor-product rows",
        3,
        4,
        4,
        5,
        4,
        4,
    ),
    Lane(
        "A000568 marked observer fiber",
        "repo/tournaments LRC",
        "repo exact finite",
        "S509/S512/S586/HYP-1977",
        "unmarked tournament classes and scalar counts mix LRC status",
        "rooted/marked observer fiber and endpoint threshold labels",
        "source/deletion/gap-threshold derivative sums",
        "use address-coordinate grammar for every quotient count",
        4,
        5,
        4,
        4,
        5,
        4,
    ),
    Lane(
        "CH forcing generic address",
        "repo/set theory",
        "standard theorem transfer",
        "S656/HYP-2232",
        "continuum cardinal scalar is independent of ZFC model",
        "forcing generic or inner-model address",
        "absoluteness/forcing-extension derivative ledger",
        "use as negative control for scalar shadows needing model labels",
        3,
        5,
        3,
        2,
        3,
        5,
    ),
    Lane(
        "Cauldron schedule word",
        "repo/games",
        "repo exact finite",
        "S618-S621",
        "same additive rule changes outcome under turn-block product schedule",
        "schedule word / player-owner state",
        "minimax value derivative under one schedule-block flip",
        "build schedule-address quotient for two-block variants",
        3,
        4,
        4,
        5,
        3,
        4,
    ),
    Lane(
        "Heegner reduced-form cycle",
        "repo/web number theory",
        "theorem transfer",
        "S651; Heegner/class-number-one thread",
        "prime-generating polynomial scalar run hides UFD/product address",
        "reduced binary quadratic form cycle / class-group address",
        "class-cycle neighbor derivative sums for polynomial values",
        "make a discriminant/form-cycle atlas for Heegner-like prime runs",
        2,
        4,
        3,
        3,
        3,
        4,
    ),
]


REPO_ARTIFACTS = {
    "S663/HYP-2239": [
        "04-computation/triune_carrier_applications_s663.py",
        "05-knowledge/results/triune_carrier_applications_s663.out",
    ],
    "S660/HYP-2236": [
        "04-computation/tournament_deck_derivative_s660.py",
        "05-knowledge/results/tournament_deck_derivative_s660.out",
    ],
    "S635/S636/HYP-2211/HYP-2212": [
        "04-computation/rational_shadow_carrier_s635.py",
        "04-computation/quadratic_pi_e_carrier_s636.py",
    ],
    "S622/S626/S628/HYP-2203": [
        "04-computation/unit_distance_impairment_lab_s622.py",
        "04-computation/unit_distance_unit_spine_tournament_s626.py",
        "04-computation/unit_distance_spine_ladder_s628.py",
    ],
    "S659/HYP-2235": [
        "04-computation/finite_field_kakeya_falconer_s659.py",
        "05-knowledge/results/finite_field_kakeya_falconer_s659.out",
    ],
    "S658/HYP-2234": [
        "04-computation/frontier_small_carriers_s658.py",
        "05-knowledge/results/frontier_small_carriers_s658.out",
    ],
    "S642/S643/HYP-2218/HYP-2219": [
        "04-computation/goldbach_lemoine_pair_projection_s642.py",
        "04-computation/goldbach_lemoine_pair_bridge_s643.py",
    ],
    "S644/S645/HYP-2220": [
        "04-computation/vieta_perfect_aliquot_carriers_s644.py",
        "04-computation/perfect_numbers_aliquot_carrier_s645.py",
    ],
    "S509/S512/S586/HYP-1977": [
        "04-computation/lrc_a000568_iso_analogy_s509.py",
        "04-computation/lrc_iso_class_constraint_s512.py",
        "04-computation/tournament_perspective_observer_coupling_s586.py",
    ],
    "S656/HYP-2232": [
        "04-computation/continuum_hypothesis_carrier_s656.py",
        "05-knowledge/results/continuum_hypothesis_carrier_s656.out",
    ],
    "S618-S621": [
        "04-computation/cauldron_game_s618.py",
        "04-computation/cauldron_two_block_adversarial_s620.py",
        "04-computation/cauldron_block_turn_minimax_s621.py",
    ],
}


WEB_SOURCES = [
    (
        "DLMF Euler products and logarithmic derivatives",
        "https://dlmf.nist.gov/27.4",
        "Euler products unroll into prime-power/Mangoldt-style derivative sums",
    ),
    (
        "explicit formulae for L-functions",
        "https://en.wikipedia.org/wiki/Explicit_formulae_for_L-functions",
        "Euler product/log-derivative prime-power sums plus zero sums",
    ),
    (
        "Jacobi determinant formula",
        "https://people.eecs.berkeley.edu/~wkahan/MathH110/jacobi.pdf",
        "determinant products differentiate through adjugate/cofactor addresses",
    ),
    (
        "graph reconstruction algebraic formulation",
        "https://arxiv.org/abs/1301.4121",
        "deck collisions, Kocay linear constraints, covering-number matrices",
    ),
    (
        "graph reconstruction pair parameters",
        "https://arxiv.org/abs/2601.00620",
        "new reconstructible vertex-pair address parameters",
    ),
    (
        "sparse phase retrieval with Prony",
        "https://www.frontiersin.org/articles/10.3389/fams.2017.00005/full",
        "Fourier intensity/autocorrelation data needs phase and support addresses",
    ),
    (
        "turnpike and beltway reconstruction",
        "https://arxiv.org/abs/1804.02465",
        "unassigned distance products need pair/end-owner addresses",
    ),
    (
        "inverse spectral Weyl m-function",
        "https://arxiv.org/abs/math/9906118",
        "m-function/spectral measure determines potential and has A-equation",
    ),
    (
        "structured matrices and continued fractions",
        "https://epubs.siam.org/doi/10.1137/090781127",
        "Hankel matrices, moment problems, J-fractions, root localization",
    ),
    (
        "continued fractions and Hankel determinants",
        "https://arxiv.org/abs/1907.05204",
        "continued fractions generate Hankel determinant formulae and Lax pairs",
    ),
    (
        "Tutte polynomial deletion-contraction",
        "https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.polynomials.tutte_polynomial.html",
        "deletion-contraction and activity-address expansions",
    ),
    (
        "Hodge theory for combinatorial geometries",
        "https://arxiv.org/abs/1511.02888",
        "matroid Chow-ring address proves log-concavity of characteristic coefficients",
    ),
    (
        "Selberg trace formula",
        "https://mathworld.wolfram.com/SelbergTraceFormula.html",
        "spectrum/geodesic trace formula as closed-orbit derivative sum",
    ),
]


def artifact_status() -> dict[str, tuple[int, int]]:
    rows = {}
    for label, rels in REPO_ARTIFACTS.items():
        exists = sum(1 for rel in rels if (ROOT / rel).exists())
        rows[label] = (exists, len(rels))
    return rows


def tournament(lanes: list[Lane]) -> dict[str, object]:
    n = len(lanes)
    adj = [[0] * n for _ in range(n)]
    out = [0] * n
    edges = []
    for i, j in combinations(range(n), 2):
        vi = lanes[i].score_vector
        vj = lanes[j].score_vector
        wi = sum(a > b for a, b in zip(vi, vj))
        wj = sum(a < b for a, b in zip(vi, vj))
        if wi > wj or (wi == wj and i < j):
            winner, loser = i, j
        else:
            winner, loser = j, i
        adj[winner][loser] = 1
        out[winner] += 1
        edges.append((lanes[winner].name, lanes[loser].name))

    c3 = 0
    for a, b, c in combinations(range(n), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            c3 += 1

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
    sccs: list[list[str]] = []
    for start in reversed(order):
        if seen[start]:
            continue
        comp = []
        q = deque([start])
        seen[start] = True
        while q:
            v = q.popleft()
            comp.append(lanes[v].name)
            for w in radj[v]:
                if not seen[w]:
                    seen[w] = True
                    q.append(w)
        sccs.append(comp)

    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if not (mask & (1 << nxt)) and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val

    return {
        "outscores": {lanes[i].name: out[i] for i in range(n)},
        "score_hist": dict(sorted(Counter(out).items())),
        "directed_3cycles": c3,
        "scc_sizes": sorted([len(c) for c in sccs], reverse=True),
        "hamiltonian_paths": sum(dp[-1]),
        "top_order": [lanes[i].name for i in sorted(range(n), key=lambda k: (-out[k], -lanes[k].total, k))],
        "edges_sample": edges[:20],
    }


def primes_upto(limit: int) -> list[int]:
    sieve = [True] * (limit + 1)
    sieve[0:2] = [False, False]
    for p in range(2, int(limit**0.5) + 1):
        if sieve[p]:
            for q in range(p * p, limit + 1, p):
                sieve[q] = False
    return [i for i in range(limit + 1) if sieve[i]]


def finite_vieta_lab(p: int = 17) -> dict[str, object]:
    sum_product: defaultdict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    full: defaultdict[tuple[int, int, int], list[tuple[int, int]]] = defaultdict(list)
    power_sums = {}
    for x in range(p):
        for y in range(p):
            s = (x + y) % p
            prod = (x * y) % p
            branch = (x - y) % p
            sum_product[(s, prod)].append((x, y))
            full[(s, prod, branch)].append((x, y))
            power_sums[(s, prod, branch)] = tuple(((x**k + y**k) % p) for k in range(1, 7))
    collision_key, collision_vals = next((k, v) for k, v in sum_product.items() if len(v) == 2)
    branches = [((x - y) % p, power_sums[(collision_key[0], collision_key[1], (x - y) % p)]) for x, y in collision_vals]
    return {
        "field": p,
        "sum_product_max_fiber": max(len(v) for v in sum_product.values()),
        "full_address_max_fiber": max(len(v) for v in full.values()),
        "collision_example": (collision_key, collision_vals),
        "branch_power_sums": branches,
    }


def log_derivative_lab(limit: int = 13, terms: int = 20) -> dict[str, object]:
    primes = primes_upto(limit)
    von_mangoldt: dict[int, int] = {}
    for p in primes:
        pk = p
        while pk <= terms:
            von_mangoldt[pk] = p
            pk *= p
    # For product prod_p (1 - x^p)^-1, x d/dx log(product)
    # has coefficient p on every multiple of p.  This is the toy version of
    # "product -> derivative sum".
    coeffs = []
    for n in range(1, terms + 1):
        coeffs.append((n, sum(p for p in primes if n % p == 0), von_mangoldt.get(n, 0)))
    return {
        "primes": primes,
        "coefficients_n_sum_p_dividing_n_and_prime_power_marker": coeffs,
    }


def determinant_adjugate_lab() -> dict[str, object]:
    rows = []
    for a, b in [(2, 6), (3, 4)]:
        rows.append(
            {
                "diag": (a, b),
                "determinant": a * b,
                "derivative_in_E11_direction": b,
                "adjugate_diagonal": (b, a),
            }
        )
    return {
        "collision": "diag(2,6) and diag(3,4) both have determinant 12",
        "rows": rows,
        "lesson": "the adjugate row/column address splits the determinant collision",
    }


def turnpike_homometric_lab() -> dict[str, object]:
    a = (0, 1, 2, 6, 8, 11)
    b = (0, 1, 6, 7, 9, 11)

    def distances(points: tuple[int, ...]) -> tuple[int, ...]:
        return tuple(
            sorted(points[j] - points[i] for i, j in combinations(range(len(points)), 2))
        )

    def adjacent_gaps(points: tuple[int, ...]) -> tuple[int, ...]:
        return tuple(points[i + 1] - points[i] for i in range(len(points) - 1))

    return {
        "sets": (a, b),
        "same_distance_multiset": distances(a) == distances(b),
        "distance_multiset": distances(a),
        "coordinate_sums": (sum(a), sum(b)),
        "adjacent_gaps": (adjacent_gaps(a), adjacent_gaps(b)),
        "lesson": "distance product shadows need pair/end-owner addresses",
    }


def goldbach_lemoine_lab(limit: int = 100) -> dict[str, object]:
    primes = primes_upto(limit)
    shared = []
    for p in primes:
        for q in primes:
            if p > q:
                continue
            even = p + q
            odd = p + 2 * q
            if even <= limit and odd <= 2 * limit:
                shared.append((p, q, even, odd, odd - even, 2 * even - odd))
    diagonals = [(p, even, odd) for p, q, even, odd, _, _ in shared if p == q]
    return {
        "shared_pair_rows": len(shared),
        "first_shared_rows_p_q_E_O_qrec_prec": shared[:10],
        "diagonal_p_2p_3p": diagonals[:10],
    }


def family_summary(lanes: list[Lane]) -> dict[str, object]:
    by_family: defaultdict[str, list[Lane]] = defaultdict(list)
    for lane in lanes:
        by_family[lane.family].append(lane)
    return {
        family: {
            "count": len(rows),
            "best": max(rows, key=lambda lane: (lane.total, lane.collision_evidence)).name,
            "total_score": sum(lane.total for lane in rows),
        }
        for family, rows in sorted(by_family.items())
    }


def main() -> None:
    print("=" * 78)
    print("S665 address-coordinate derivative repair atlas")
    print("=" * 78)
    print()
    print("Grammar")
    print("  1. find a scalar/product collision")
    print("  2. attach the missing fraction/address coordinate")
    print("  3. unroll it back into a derivative sum")
    print()

    print("A. Source anchors")
    for label, url, role in WEB_SOURCES:
        print(f"  {label}: {role}")
        print(f"    {url}")
    print()

    print("B. Repo artifact availability")
    for label, (exists, total) in sorted(artifact_status().items()):
        print(f"  {label:<34} {exists}/{total} artifacts present")
    print()

    print("C. Candidate repair lanes")
    print(f"{'lane':<38} {'family':<22} {'total':>5} vector")
    for lane in sorted(LANES, key=lambda x: (-x.total, x.name)):
        print(f"{lane.name:<38} {lane.family:<22} {lane.total:5d} {lane.score_vector}")
    print()

    print("D. Exact-ish lanes")
    for lane in [lane for lane in LANES if lane.exactish]:
        print(f"  {lane.name}")
        print(f"    collision:  {lane.scalar_product_collision}")
        print(f"    address:    {lane.address_coordinate}")
        print(f"    derivative: {lane.derivative_sum}")
        print(f"    next:       {lane.next_probe}")
    print()

    fp = tournament(LANES)
    print("E. Tournament Analysis over repair lanes")
    print("  vertices=repair lanes")
    print("  observable=(collision, address, derivative, finite, repo, source)")
    print(f"  score_hist={fp['score_hist']}")
    print(f"  directed_3cycles={fp['directed_3cycles']}")
    print(f"  scc_sizes={fp['scc_sizes']}")
    print(f"  hamiltonian_paths={fp['hamiltonian_paths']}")
    print("  top_order:")
    for name in fp["top_order"][:12]:
        print(f"    {name}")
    print()

    print("F. Family summary")
    for family, row in family_summary(LANES).items():
        print(f"  {family:<24} count={row['count']} total={row['total_score']} best={row['best']}")
    print()

    print("G. Finite side labs")
    print("  Vieta branch over F_17")
    vl = finite_vieta_lab()
    print(f"    sum_product_max_fiber={vl['sum_product_max_fiber']}")
    print(f"    full_address_max_fiber={vl['full_address_max_fiber']}")
    print(f"    collision_example={vl['collision_example']}")
    print(f"    branch_power_sums={vl['branch_power_sums']}")
    print("  Toy log derivative of product prod_p (1 - x^p)^-1")
    ld = log_derivative_lab()
    print(f"    primes={ld['primes']}")
    print(f"    coeffs(n, sum_p|n p, prime_power_marker)={ld['coefficients_n_sum_p_dividing_n_and_prime_power_marker']}")
    print("  Determinant/adjugate collision")
    da = determinant_adjugate_lab()
    print(f"    collision={da['collision']}")
    print(f"    rows={da['rows']}")
    print(f"    lesson={da['lesson']}")
    print("  Turnpike homometric collision")
    th = turnpike_homometric_lab()
    print(f"    sets={th['sets']}")
    print(f"    same_distance_multiset={th['same_distance_multiset']}")
    print(f"    distance_multiset={th['distance_multiset']}")
    print(f"    coordinate_sums={th['coordinate_sums']}")
    print(f"    adjacent_gaps={th['adjacent_gaps']}")
    print("  Goldbach/Lemoine reconstruction")
    gl = goldbach_lemoine_lab()
    print(f"    shared_pair_rows={gl['shared_pair_rows']}")
    print(f"    first_rows={gl['first_shared_rows_p_q_E_O_qrec_prec']}")
    print(f"    diagonals={gl['diagonal_p_2p_3p']}")
    print()

    print("H. Synthesis")
    print("  Address coordinates are not decorations.  They are the derivative variable.")
    print("  In exact repo cases, unpaired/global scalar repairs fail; paired/addressed")
    print("  derivatives split the collision.  In external theorem families, the same")
    print("  mechanism appears as log derivatives, adjugates, Weyl m-functions,")
    print("  J-fractions, Prony reconstructions, turnpike owner maps, Kocay matrices,")
    print("  activity expansions, trace formulas, and Hodge operators.")
    print("  Best next builds: LRC owner derivative, tournament n=7 paired deck, OCF")
    print("  noncommutative address projection, and unit-distance deletion-owner lab.")


if __name__ == "__main__":
    main()
