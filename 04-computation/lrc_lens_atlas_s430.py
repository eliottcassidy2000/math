#!/usr/bin/env python3
"""
lrc_lens_atlas_s430.py

codex-2026-05-31 S430

Build an atlas of mathematical lenses for the Lonely Runner Conjecture and
compare each lens to the tournament machinery in this repository.

The goal is not to prove a new case.  It is to answer a structural question:
how many genuinely different-looking problems are isomorphic to the same LRC
core, and which repo/tournament lens does each one activate?
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass


@dataclass(frozen=True)
class Lens:
    key: str
    name: str
    status: str
    atom: str
    certificate: str
    tournament_shadow: str
    sources: tuple[str, ...]


LENSES: tuple[Lens, ...] = (
    Lens(
        "runner_frame",
        "moving runners / stationary runner",
        "exact formulation",
        "relative speed v_i",
        "exists t with all ||v_i t|| >= 1/(k+1)",
        "choose a distinguished vertex and measure all comparisons relative to it",
        ("Wills/Cusick/Goddyn", "THM-357"),
    ),
    Lens(
        "diophantine_norm",
        "simultaneous Diophantine approximation",
        "exact formulation",
        "integer vector v",
        "max_t min_i ||v_i t|| >= 1/(k+1)",
        "score/height functional on a labelled complete decision vector",
        ("Cambridge finite-checking intro", "THM-369"),
    ),
    Lens(
        "forbidden_intervals",
        "forbidden interval cover on R/Z",
        "exact formulation",
        "open interval ||v t|| < 1/n",
        "union fails to be the whole circle",
        "good-cut coverage: an uncovered cut is a witness",
        ("THM-357", "HYP-1855"),
    ),
    Lens(
        "endpoint_hypergraph",
        "endpoint-protection hypergraph",
        "exact finite certificate",
        "forbidden endpoint (n*m+eps)/(n*v)",
        "full-measure cover and every endpoint protected",
        "labelled cut-protection cycles in tournament SCC/good-cut theory",
        ("THM-357", "THM-365", "S386"),
    ),
    Lens(
        "coarse_sieve",
        "small-denominator divisibility sieve",
        "proved necessary certificate",
        "reduced fraction a/q with q <= n",
        "counterexample must pay every q by a divisible speed",
        "cut rows that force at least one crossing backward arc",
        ("THM-366", "THM-369"),
    ),
    Lens(
        "torus_line",
        "one-dimensional subtorus of (R/Z)^k",
        "exact formulation",
        "line t*v modulo Z^k",
        "line hits the central cube",
        "orbit of a tournament under a one-parameter flip/order flow",
        ("Giri-Kravitz spectrum", "HYP-1893"),
    ),
    Lens(
        "view_obstruction",
        "view obstruction by centered cubes",
        "exact geometric formulation",
        "ray through a periodic obstacle field",
        "a line of sight misses all obstacles",
        "visibility through the tournament tiling hypercube",
        ("Cusick", "Henze-Malikiosis"),
    ),
    Lens(
        "distance_graph",
        "distance graph / circular chromatic number",
        "equivalent graph-coloring formulation",
        "integer-distance edge set D",
        "circular coloring threshold exceeds the LRC bound",
        "chromatic/fugacity axis of tournament conflict graphs",
        ("distance graph literature", "chromatic tournament thread"),
    ),
    Lens(
        "flows_matroids",
        "nowhere-zero flows and regular matroids",
        "equivalent implication family",
        "flow value set",
        "reduce arbitrary k values to [1..k]",
        "cycle-space constraints behind tournament OCF/flow shadows",
        ("Bienia-Goddyn-Gvozdjak-Sebo-Tarsi",),
    ),
    Lens(
        "zonotope_radius",
        "lattice zonotope covering radius",
        "exact convex-geometric restatement",
        "LR zonotope with n generators in dimension n-1",
        "covering/intersection condition at center",
        "zonotope/polytope shadow of the tournament staircase",
        ("Henze-Malikiosis", "Malikiosis-Santos-Schymura"),
    ),
    Lens(
        "finite_checking",
        "finite checking / minimal counterexample bound",
        "proof technology",
        "bounded primitive integer velocity vector",
        "check all speeds below a computable threshold",
        "finite labelled tournament enumeration after quotienting symmetries",
        ("Tao", "linearly-exponential checking", "Rosenfeld"),
    ),
    Lens(
        "product_prime_sieve",
        "prime-product counterexample sieve",
        "proof technology",
        "minimal counterexample product",
        "forced prime divisibility outruns the finite-check bound",
        "Lucas-active parity/sieve layers in tournament bucket transport",
        ("Rosenfeld 8", "Trakulthongchai 9/10", "Sungkawichai-Trakulthongchai 11-13"),
    ),
    Lens(
        "integer_program",
        "integer program / set cover",
        "repo proof technology",
        "speed column vs denominator/endpoint row",
        "no n-1 primitive columns cover all rows",
        "tournament good-cut incidence matrix and dual row weights",
        ("HYP-1890", "S420"),
    ),
    Lens(
        "bruhat_tits",
        "p-adic / Bruhat-Tits endpoint frontier",
        "repo proof technology",
        "endpoint denominator depth vector",
        "all-repaired branch has positive frontier mass",
        "valuation skeleton of natural-number multiplication shadow",
        ("HYP-1868", "HYP-1880", "S399-S410"),
    ),
    Lens(
        "adelic_product",
        "adelic gap-debt product",
        "repo proof technology",
        "real gap times p-adic endpoint debt",
        "avoid the simultaneous zero corner",
        "AFM-style conservation: frustration cannot vanish in all projections",
        ("HYP-1866", "HYP-1892"),
    ),
    Lens(
        "spectrum",
        "Lonely Runner spectra",
        "structural exact lens",
        "value D(T) for proper subtori",
        "spectrum stays below the LRC bound",
        "H-spectrum / root-spectrum style invariant for tournaments",
        ("Giri-Kravitz spectrum",),
    ),
    Lens(
        "mixed_thresholds",
        "mixed thresholds / Fourier bands",
        "variant lens",
        "threshold vector d_i",
        "nonempty simultaneous safe band for unequal radii",
        "weighted tournament arcs / nonuniform cut thresholds",
        ("Jensen mixed thresholds",),
    ),
    Lens(
        "shifted_lrc",
        "shifted LRC / coloopless-cosimple zonotopes",
        "variant lens with counterexamples",
        "runners with shifted starts or zonotope translations",
        "covering-radius statement can fail outside the original class",
        "warning: tournament complement/tiling complement analogies can split",
        ("shifted LRC", "coloopless zonotopes"),
    ),
    Lens(
        "natural_modes",
        "natural-number operation modes",
        "repo organizing lens",
        "n = 2^r * odd root",
        "row x*2 and top-row x+2 moves preserve obstruction mass",
        "tournament n->2n blowup and n->n+2 Mode-B recursion",
        ("HYP-1881", "HYP-1891", "THM-362"),
    ),
    Lens(
        "product_sum",
        "product-sum / Egyptian reciprocal ledger",
        "repo exploratory lens",
        "additive slack equals multiplicative factor packing",
        "primitive breaker creates gap or private endpoint leaf",
        "additive shadow is transitive tournament; product shadow is divisibility DAG",
        ("S365", "HYP-1865", "S397"),
    ),
)


FORMULATION_EDGE_PAIRS: tuple[tuple[str, str], ...] = (
    ("runner_frame", "diophantine_norm"),
    ("diophantine_norm", "forbidden_intervals"),
    ("forbidden_intervals", "endpoint_hypergraph"),
    ("diophantine_norm", "torus_line"),
    ("torus_line", "view_obstruction"),
    ("torus_line", "zonotope_radius"),
    ("diophantine_norm", "distance_graph"),
    ("diophantine_norm", "flows_matroids"),
    ("diophantine_norm", "coarse_sieve"),
    ("torus_line", "spectrum"),
)

PROOF_EDGE_PAIRS: tuple[tuple[str, str], ...] = (
    ("finite_checking", "product_prime_sieve"),
    ("product_prime_sieve", "coarse_sieve"),
    ("coarse_sieve", "integer_program"),
    ("endpoint_hypergraph", "integer_program"),
    ("endpoint_hypergraph", "bruhat_tits"),
    ("bruhat_tits", "adelic_product"),
    ("natural_modes", "bruhat_tits"),
    ("natural_modes", "product_sum"),
    ("integer_program", "natural_modes"),
    ("product_sum", "distance_graph"),
)


def component_sizes(edges: tuple[tuple[str, str], ...]) -> list[list[str]]:
    graph: dict[str, set[str]] = defaultdict(set)
    for a, b in edges:
        graph[a].add(b)
        graph[b].add(a)

    seen: set[str] = set()
    components: list[list[str]] = []
    for lens in LENSES:
        if lens.key in seen or lens.key not in graph:
            continue
        queue = deque([lens.key])
        seen.add(lens.key)
        comp = []
        while queue:
            node = queue.popleft()
            comp.append(node)
            for nxt in sorted(graph[node]):
                if nxt not in seen:
                    seen.add(nxt)
                    queue.append(nxt)
        components.append(sorted(comp))
    return sorted(components, key=lambda comp: (-len(comp), comp[0]))


def print_header(title: str) -> None:
    print(title)
    print("=" * 92)


def print_lens_inventory() -> None:
    print_header("1. LRC LENS INVENTORY")
    counts = Counter(lens.status for lens in LENSES)
    print(f"Total lenses recorded: {len(LENSES)}")
    for status, count in sorted(counts.items()):
        print(f"  {status:<34} {count}")
    print()
    for i, lens in enumerate(LENSES, 1):
        print(f"{i:>2}. {lens.name}")
        print(f"    status:      {lens.status}")
        print(f"    atom:        {lens.atom}")
        print(f"    certificate: {lens.certificate}")
        print(f"    tournament:  {lens.tournament_shadow}")
        print(f"    sources:     {', '.join(lens.sources)}")
        print()


def print_equivalence_graph() -> None:
    print_header("2. FORMULATION / TRANSPORT GRAPH")
    print("Edges here mean: exact restatement, standard transport, or necessary")
    print("constraint family living on the same LRC core.  Status lines above")
    print("separate strict equivalences from proof constraints.")
    for a, b in FORMULATION_EDGE_PAIRS:
        print(f"  {a:22} <-> {b}")
    print()
    for comp in component_sizes(FORMULATION_EDGE_PAIRS):
        print(f"component size {len(comp):>2}: {', '.join(comp)}")
    print()

    print_header("3. PROOF-TECHNOLOGY GRAPH")
    print("Edges here mean: same obstruction measured with different proof machinery.")
    for a, b in PROOF_EDGE_PAIRS:
        print(f"  {a:22} <-> {b}")
    print()
    for comp in component_sizes(PROOF_EDGE_PAIRS):
        print(f"component size {len(comp):>2}: {', '.join(comp)}")
    print()


def print_tournament_dictionary() -> None:
    print_header("4. COMMON INCIDENCE CORE WITH TOURNAMENTS")
    rows = [
        (
            "LRC forbidden endpoint",
            "row to be protected",
            "Hamiltonian-path cut / good cut",
        ),
        (
            "LRC speed interval",
            "column/protector",
            "backward arc crossing a cut",
        ),
        (
            "uncovered endpoint",
            "lonely witness",
            "bad cut / non-SC certificate",
        ),
        (
            "all endpoints protected",
            "counterexample certificate",
            "every cut protected inside SCC",
        ),
        (
            "unit endpoint a/n",
            "coarse denominator row",
            "boundary cut forced by a single crossing",
        ),
        (
            "p-adic endpoint depth",
            "row export after repair",
            "projection defect after deleting/adding vertices",
        ),
        (
            "finite-check speed bound",
            "finite labelled search modulo symmetries",
            "finite tournament enumeration modulo S_n/complement",
        ),
        (
            "distance/circular coloring",
            "avoid forbidden difference set",
            "color conflict graph Omega(T) at fugacity/chromatic axes",
        ),
    ]
    print(f"{'LRC object':<28} {'incidence role':<34} tournament analogue")
    print("-" * 92)
    for left, mid, right in rows:
        print(f"{left:<28} {mid:<34} {right}")
    print()


def print_new_questions() -> None:
    print_header("5. QUESTIONS GENERATED BY THE ATLAS")
    questions = [
        "Can the finite-checking prime-product sieve be rewritten as a dual row-weight certificate in the S420 IP matrix?",
        "Does the Giri-Kravitz spectrum recursion S(n) -> S(n-1) have a direct tournament analogue as deletion of a vertex or quotient by a cut?",
        "Can view-obstruction cubes be replaced by a tournament tiling polytope whose facets are good cuts and odd-cycle conflict rows?",
        "Is the row/column natural-number split the common skeleton behind zonotope paving, p-adic endpoint depth, and product-sum factor packing?",
        "Can a counterexample be defined without speeds first: as a leafless labelled protection hypergraph satisfying all coarse sieve rows and then asking realizability?",
    ]
    for i, question in enumerate(questions, 1):
        print(f"{i}. {question}")


def main() -> None:
    print("LRC lens atlas (codex-2026-05-31 S430)")
    print()
    print_lens_inventory()
    print_equivalence_graph()
    print_tournament_dictionary()
    print_new_questions()


if __name__ == "__main__":
    main()
