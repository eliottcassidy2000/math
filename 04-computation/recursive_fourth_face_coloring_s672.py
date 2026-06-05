"""S672: the recursive fourth face behind sum/product/fraction carriers.

The repo's triune carrier program used:

    sum      = additive packets / moments / walls
    product  = local factors / shells / obstruction ledgers
    fraction = boundary states / owners / branches

This scout splits a fourth face out of the fraction face:

    recursion = transition law under iteration or outer extension.

The point is not that recursion is another static formula.  It is the proof
engine that says how a retained boundary state evolves.  Zhou-Markov's
irrationality proofs use integral recurrences this way, and Paris-Harrington
colorings show the extreme version: finite coloring witnesses whose growth
outstrips PA-provably total functions.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, permutations


def banner(title: str) -> None:
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def poly_add(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        out[i] = (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(out)


def poly_scale(a: tuple[int, ...], c: int) -> tuple[int, ...]:
    return tuple(c * x for x in a)


def poly_shift_t(a: tuple[int, ...]) -> tuple[int, ...]:
    return (0,) + a


def poly_fmt(a: tuple[int, ...]) -> str:
    terms: list[str] = []
    for power, coeff in enumerate(a):
        if coeff == 0:
            continue
        sign = "-" if coeff < 0 else "+"
        mag = abs(coeff)
        if power == 0:
            body = str(mag)
        elif power == 1:
            body = "t" if mag == 1 else f"{mag}*t"
        else:
            body = f"t^{power}" if mag == 1 else f"{mag}*t^{power}"
        terms.append((sign, body))
    if not terms:
        return "0"
    first_sign, first_body = terms[0]
    s = ("" if first_sign == "+" else "-") + first_body
    for sign, body in terms[1:]:
        s += f" {sign} {body}"
    return s


def zhou_markov_pi_polys(count: int = 8) -> list[tuple[int, ...]]:
    """Polynomials P_n(t) for I_n=P_n(pi^2) in the pi recurrence.

    In the paper's pi proof, I_0=2, I_1=4 and
        I_n = (4n-2)I_{n-1} - pi^2 I_{n-2}.
    We store t=pi^2.
    """

    polys = [(2,), (4,)]
    for n in range(2, count):
        prev = poly_scale(polys[-1], 4 * n - 2)
        back = poly_scale(poly_shift_t(polys[-2]), -1)
        polys.append(poly_add(prev, back))
    return polys


def has_large_homogeneous(coloring: dict[tuple[int, int], int], n: int) -> bool:
    """Partial check for the baby Paris-Harrington edge-coloring case.

    Vertices are 1..n.  A homogeneous set H is "large" when |H| >= min(H).
    Only fully assigned candidate subsets are tested.
    """

    vertices = range(1, n + 1)
    for size in range(3, n + 1):
        for subset in combinations(vertices, size):
            if size < min(subset):
                continue
            values: list[int] = []
            complete = True
            for edge in combinations(subset, 2):
                if edge not in coloring:
                    complete = False
                    break
                values.append(coloring[edge])
            if complete and len(set(values)) == 1:
                return True
    return False


def count_large_avoiding_edge_colorings(n: int) -> tuple[int, int]:
    """Count 2-colorings of K_n avoiding large homogeneous sets.

    This is only the x=2 baby case of the user's PH function:
    color 2-tuples by 2 colors and seek a large homogeneous set of size 3.
    It is tiny, but it exposes the recursive extension-state viewpoint.
    """

    edges = list(combinations(range(1, n + 1), 2))
    coloring: dict[tuple[int, int], int] = {}
    count = 0
    nodes = 0

    def rec(i: int) -> None:
        nonlocal count, nodes
        nodes += 1
        if i == len(edges):
            count += 1
            return
        edge = edges[i]
        for color in (0, 1):
            coloring[edge] = color
            if not has_large_homogeneous(coloring, n):
                rec(i + 1)
            del coloring[edge]

    rec(0)
    return count, nodes


@dataclass(frozen=True)
class Face:
    name: str
    static_exactness: int
    local_factor_leverage: int
    boundary_memory: int
    outer_extension_survival: int
    growth_witness: int
    algorithmic_use: int
    coloring_transfer: int

    @property
    def vector(self) -> tuple[int, ...]:
        return (
            self.static_exactness,
            self.local_factor_leverage,
            self.boundary_memory,
            self.outer_extension_survival,
            self.growth_witness,
            self.algorithmic_use,
            self.coloring_transfer,
        )


CRITERIA = (
    "static_exactness",
    "local_factor_leverage",
    "boundary_memory",
    "outer_extension_survival",
    "growth_witness",
    "algorithmic_use",
    "coloring_transfer",
)

FACES = [
    Face("sum", 5, 2, 1, 1, 1, 3, 4),
    Face("product", 3, 5, 2, 2, 2, 4, 4),
    Face("fraction", 2, 3, 5, 3, 3, 4, 3),
    Face("recursion", 4, 4, 4, 5, 5, 5, 5),
    Face("raw_scalar", 1, 1, 0, 0, 1, 1, 1),
]


def face_tournament(faces: list[Face]) -> dict[str, set[str]]:
    priority = {name: i for i, name in enumerate(["recursion", "fraction", "product", "sum", "raw_scalar"])}
    wins: dict[str, set[str]] = {f.name: set() for f in faces}
    by_name = {f.name: f for f in faces}
    for a, b in combinations(by_name, 2):
        va, vb = by_name[a].vector, by_name[b].vector
        score_a = sum(x > y for x, y in zip(va, vb))
        score_b = sum(y > x for x, y in zip(va, vb))
        if score_a > score_b:
            wins[a].add(b)
        elif score_b > score_a:
            wins[b].add(a)
        elif priority[a] < priority[b]:
            wins[a].add(b)
        else:
            wins[b].add(a)
    return wins


def count_directed_triangles(wins: dict[str, set[str]]) -> int:
    total = 0
    for a, b, c in combinations(wins, 3):
        ab = b in wins[a]
        bc = c in wins[b]
        ca = a in wins[c]
        ba = a in wins[b]
        cb = b in wins[c]
        ac = c in wins[a]
        if (ab and bc and ca) or (ba and cb and ac):
            total += 1
    return total


def hamiltonian_paths(wins: dict[str, set[str]]) -> int:
    names = list(wins)
    total = 0
    for path in permutations(names):
        if all(path[i + 1] in wins[path[i]] for i in range(len(path) - 1)):
            total += 1
    return total


DOMAIN_ATLAS = [
    {
        "domain": "Zhou-Markov trig irrationality",
        "sum": "small positive/oscillatory integral packets",
        "product": "rational denominator clearing b^n or b^(2n+1)",
        "fraction": "Lambert-style continued fractions become optional boundary memory",
        "recursion": "I_n=(4n-2)I_{n-1}-r^2 I_{n-2}; four-state cos recurrence plus descent",
        "test": "nonzero integer forced into (-1,1)",
    },
    {
        "domain": "Paris-Harrington colorings",
        "sum": "finite coloring/counting of x-tuples",
        "product": "Ramsey parameter product: arity, colors, size, horizon",
        "fraction": "relative largeness |H| >= min(H) is a boundary/observer state",
        "recursion": "minimal witness sigma(n,k) eventually dominates PA-provably total functions",
        "test": "finite coloring theorem with proof-strength growth",
    },
    {
        "domain": "LRC14 carry-owner seam",
        "sum": "odd-wall and pair-sum danger packets",
        "product": "C=27 gcd shells and local residue factors",
        "fraction": "carry-owner continuant v=r+27k",
        "recursion": "coherent +27 outer extensions must preserve no-leak or pay strict tax",
        "test": "floor vs strict fiber purity",
    },
    {
        "domain": "A000568 endpoint enumeration",
        "sum": "class counts and streamed H/c3/SCC payloads",
        "product": "Aut(parent)-orbit incident cubes",
        "fraction": "deleted-card with L/M/U owner side",
        "recursion": "n -> n+1 endpoint extension purity conjecture",
        "test": "Phi injective through n=8",
    },
    {
        "domain": "Friedman/outer-extension toy",
        "sum": "finite tree/color histograms",
        "product": "embedding/downset probes as local factors",
        "fraction": "outer insertion address",
        "recursion": "one-leaf extension preserving homeomorphic color-chain predicates",
        "test": "coarse quotients leak; embedding profiles pure",
    },
    {
        "domain": "unit-distance carrier",
        "sum": "unit-edge spine and edge-gain counts",
        "product": "Eisenstein direction/norm support",
        "fraction": "point-deletion frontier owner",
        "recursion": "Moser slab/ear extension ledger",
        "test": "n=21/22 edge optimality obstruction search",
    },
    {
        "domain": "cauldron/Schur coloring",
        "sum": "A+B=C additive witnesses",
        "product": "color-class residue/sieve patterns",
        "fraction": "active cauldron/removal state",
        "recursion": "online placement game and first-boil/last-boil dynamics",
        "test": "finite-sums and removal variants split ordinary Schur shadows",
    },
]


def main() -> None:
    banner("S672 recursive fourth face scout")
    print("Thesis: fraction is a boundary state; recursion is the transition law")
    print("that proves the state survives, shrinks, or explodes under extension.")

    banner("A. Zhou-Markov recurrence face")
    print("arXiv:0911.1933 uses recurrences of integrals as proof engines.")
    print("For the pi proof, with t=pi^2:")
    for i, poly in enumerate(zhou_markov_pi_polys(9)):
        print(f"  P_{i}(t) = {poly_fmt(poly)}")
    print("The contradiction shape is: denominator clearing makes an integer;")
    print("the integral bound makes it tend to 0; recurrence/descent keeps a")
    print("nonzero subsequence alive.  That is a recursive representation.")

    banner("B. Baby Paris-Harrington coloring recursion")
    print("2-color edge case: avoid homogeneous H with |H|>=3 and |H|>=min(H).")
    print("This baby case collapses at N=6, but the state-machine viewpoint is")
    print("the same one that becomes proof-theoretically huge in full PH.")
    print(f"{'N':>3} {'avoiding_colorings':>20} {'search_nodes':>14}")
    for n in range(3, 7):
        count, nodes = count_large_avoiding_edge_colorings(n)
        print(f"{n:>3} {count:>20} {nodes:>14}")

    banner("C. Four-face carrier atlas")
    for row in DOMAIN_ATLAS:
        print(f"- {row['domain']}")
        print(f"    sum:       {row['sum']}")
        print(f"    product:   {row['product']}")
        print(f"    fraction:  {row['fraction']}")
        print(f"    recursion: {row['recursion']}")
        print(f"    test:      {row['test']}")

    banner("D. Face tournament")
    wins = face_tournament(FACES)
    score_hist = Counter(len(v) for v in wins.values())
    print("vertices=representation faces")
    print("observable=(static exactness, local factors, boundary memory,")
    print("            outer-extension survival, growth witness, algorithmic use,")
    print("            coloring transfer)")
    print(f"criteria={CRITERIA}")
    for face in FACES:
        print(f"  {face.name:>10}: {face.vector}")
    print(f"score_hist={dict(sorted(score_hist.items()))}")
    print(f"directed_3cycles={count_directed_triangles(wins)}")
    print(f"hamiltonian_paths={hamiltonian_paths(wins)}")
    print("top_order:")
    for name in sorted(wins, key=lambda x: (-len(wins[x]), x)):
        print(f"  {name:>10} score={len(wins[name])} beats={sorted(wins[name])}")

    banner("E. New grammar")
    print("Old useful loop:")
    print("  sum -> product -> fraction -> sum")
    print("Sharper recursive loop:")
    print("  sum -> product -> fraction -> recursion -> sum")
    print()
    print("Read it as:")
    print("  sum: aggregate what is visible")
    print("  product: localize obstruction factors")
    print("  fraction: attach owner/boundary/branch state")
    print("  recursion: prove the state is stable under all allowed extensions")
    print("  recursion -> sum: unroll the transition into the next additive ledger")

    banner("F. Concrete next tests")
    tests = [
        "A000568 n=9: test whether half-filter Phi remains pure under endpoint recursion.",
        "LRC14: separate carry-owner continuant state from the recursion law of coherent +27 lifts.",
        "Unit distance n=21/22: keep edge/direction product fixed and vary point-deletion recursion owners.",
        "Cauldron/Schur: record active removal state as fraction, then enumerate game recursion.",
        "OCF/H(T): turn deletion/substitution continuants into a recursive representation of H payloads.",
        "PH micro-lab: build stronger small hypergraph-coloring backtracking with relative-largeness states.",
    ]
    for i, test in enumerate(tests, 1):
        print(f"{i}. {test}")


if __name__ == "__main__":
    main()
