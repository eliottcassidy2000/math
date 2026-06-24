#!/usr/bin/env python3
"""S141: regular-solid and Euclidean tiling recursion carriers for LRC14.

This synthesizes the parallel S141 carrier scouts.  The goal is not to prove
LRC14 from polyhedra.  The goal is to keep each geometric analogy labelled by
the relation it preserves before feeding it into the POKE proof tree:

* regular-map curvature boundary: spherical / Euclidean / hyperbolic;
* Platonic and Archimedean vertex-defect data;
* square, triangular, and hexagonal tiling recursion indices, with Gaussian
  and Eisenstein norm labels;
* prism/antiprism annular families, where n=14 gives 28 vertices and
  per-vertex defect 1/14;
* Johnson solids as finite nonuniform residual-atlas analogues;
* Tournament Analysis on proof carriers rather than raw solids.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations


Config = tuple[int, ...]


def polygon_angle_fraction(p: int) -> Fraction:
    """Interior angle of a regular p-gon as a fraction of 2*pi."""

    return Fraction(p - 2, 2 * p)


def curvature(config: Config) -> Fraction:
    """Normalized vertex curvature: 1 - angle_sum/(2*pi)."""

    return Fraction(1, 1) - sum((polygon_angle_fraction(p) for p in config), Fraction(0, 1))


def vertices_from_curvature(config: Config) -> Fraction | None:
    kappa = curvature(config)
    if kappa <= 0:
        return None
    return Fraction(2, 1) / kappa


def face_counts(config: Config) -> dict[int, Fraction] | None:
    vertices = vertices_from_curvature(config)
    if vertices is None:
        return None
    counts = Counter(config)
    return {p: Fraction(count * vertices, p) for p, count in sorted(counts.items())}


def regular_map_kind(p: int, q: int) -> str:
    wall = (p - 2) * (q - 2)
    if wall < 4:
        return "spherical"
    if wall == 4:
        return "Euclidean"
    return "hyperbolic"


def regular_map_counts(p: int, q: int) -> tuple[int, int, int] | None:
    """Return V,E,F for spherical regular map {p,q}; None at/above the wall."""

    denom = 2 * p + 2 * q - p * q
    if denom <= 0:
        return None
    edges = Fraction(2 * p * q, denom)
    vertices = Fraction(2 * edges, q)
    faces = Fraction(2 * edges, p)
    return int(vertices), int(edges), int(faces)


def gaussian_axis_norm(m: int) -> int:
    return m * m


def eisenstein_norm(a: int, b: int) -> int:
    return a * a - a * b + b * b


def fmt_frac(x: Fraction | None) -> str:
    if x is None:
        return "flat/negative"
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def fmt_faces(carrier: "UniformCarrier") -> str:
    faces = carrier.faces
    if faces is None:
        return "{}"
    return "{" + ", ".join(f"{p}:{fmt_frac(v)}" for p, v in faces.items()) + "}"


@dataclass(frozen=True)
class RegularMap:
    name: str
    p: int
    q: int
    dual: str
    role: str

    @property
    def config(self) -> Config:
        return tuple([self.p] * self.q)

    @property
    def kappa(self) -> Fraction:
        return curvature(self.config)

    @property
    def counts(self) -> tuple[int, int, int] | None:
        return regular_map_counts(self.p, self.q)


@dataclass(frozen=True)
class UniformCarrier:
    name: str
    config: Config
    family: str
    role: str

    @property
    def kappa(self) -> Fraction:
        return curvature(self.config)

    @property
    def vertices(self) -> Fraction | None:
        return vertices_from_curvature(self.config)

    @property
    def faces(self) -> dict[int, Fraction] | None:
        return face_counts(self.config)


PLATONIC_MAPS = [
    RegularMap("tetrahedron", 3, 3, "self", "triangle self seed / four-face closure"),
    RegularMap("cube", 4, 3, "octahedron", "square self-dual solid branch"),
    RegularMap("octahedron", 3, 4, "cube", "triangle/hex six-vertex hinge"),
    RegularMap("dodecahedron", 5, 3, "icosahedron", "pentagonal face branch"),
    RegularMap("icosahedron", 3, 5, "dodecahedron", "triangular high-valence branch"),
]

PLATONIC = [
    UniformCarrier(item.name, item.config, "Platonic", item.role) for item in PLATONIC_MAPS
]

ARCHIMEDEAN = [
    UniformCarrier("truncated tetrahedron", (3, 6, 6), "Archimedean", "triangle-to-hex truncation"),
    UniformCarrier("cuboctahedron", (3, 4, 3, 4), "Archimedean", "cube/octa dual mediator"),
    UniformCarrier("truncated cube", (3, 8, 8), "Archimedean", "square branch truncation"),
    UniformCarrier("truncated octahedron", (4, 6, 6), "Archimedean", "square/hex Kelvin-style carrier"),
    UniformCarrier("rhombicuboctahedron", (3, 4, 4, 4), "Archimedean", "square-heavy carrier"),
    UniformCarrier("truncated cuboctahedron", (4, 6, 8), "Archimedean", "three-scale square/hex carrier"),
    UniformCarrier("snub cube", (3, 3, 3, 3, 4), "Archimedean", "chiral square-triangle carrier"),
    UniformCarrier("icosidodecahedron", (3, 5, 3, 5), "Archimedean", "icosa/dodeca mediator"),
    UniformCarrier("truncated dodecahedron", (3, 10, 10), "Archimedean", "pentagonal truncation"),
    UniformCarrier("truncated icosahedron", (5, 6, 6), "Archimedean", "pentagon/hex soccer carrier"),
    UniformCarrier("rhombicosidodecahedron", (3, 4, 5, 4), "Archimedean", "tri-square-pent carrier"),
    UniformCarrier("truncated icosidodecahedron", (4, 6, 10), "Archimedean", "large three-scale carrier"),
    UniformCarrier("snub dodecahedron", (3, 3, 3, 3, 5), "Archimedean", "chiral pentagonal carrier"),
]

REGULAR_TILING_MAPS = [
    RegularMap("triangular tiling", 3, 6, "hexagonal tiling", "flat six-triangle vertex"),
    RegularMap("square tiling", 4, 4, "self", "flat self-dual square grid"),
    RegularMap("hexagonal tiling", 6, 3, "triangular tiling", "flat dual of triangular grid"),
]

REGULAR_TILINGS = [
    UniformCarrier(item.name, item.config, "Euclidean tiling", item.role) for item in REGULAR_TILING_MAPS
]

EUCLIDEAN_UNIFORM_TILINGS = [
    UniformCarrier("trihexagonal tiling", (3, 6, 3, 6), "Euclidean uniform tiling", "triangle/hex alternating"),
    UniformCarrier("truncated square tiling", (4, 8, 8), "Euclidean uniform tiling", "square-octagon flat branch"),
    UniformCarrier("truncated hexagonal tiling", (3, 12, 12), "Euclidean uniform tiling", "triangle-dodecagon flat branch"),
    UniformCarrier("rhombitrihexagonal tiling", (3, 4, 6, 4), "Euclidean uniform tiling", "tri-square-hex flat branch"),
    UniformCarrier("snub square tiling", (3, 3, 4, 3, 4), "Euclidean uniform tiling", "flat chiral square-triangle limit"),
    UniformCarrier("elongated triangular tiling", (3, 3, 3, 4, 4), "Euclidean uniform tiling", "triangle-square flat branch"),
]


def prism(n: int) -> UniformCarrier:
    return UniformCarrier(f"{n}-gonal prism", (4, 4, n), "prism family", "two n-cycles plus square belt")


def antiprism(n: int) -> UniformCarrier:
    return UniformCarrier(
        f"{n}-gonal antiprism", (3, 3, 3, n), "antiprism family", "twisted two n-cycles plus triangle belt"
    )


def recursion_data() -> list[tuple[str, str, list[tuple[str, int]], str]]:
    return [
        (
            "square self",
            "Gaussian axis scaling m+0i",
            [(f"m={m}", gaussian_axis_norm(m)) for m in range(2, 8)],
            "self-dual grid packet; safe only while axes remain labelled",
        ),
        (
            "triangle self",
            "Eisenstein axis scaling m+0*omega, with dyadic spine",
            [(f"m={m}", eisenstein_norm(m, 0)) for m in range(2, 8)]
            + [(f"4^{k}", 4**k) for k in range(1, 6)],
            "orientation-sector packet; dyadic chain is the clean triangle spine",
        ),
        (
            "triangle <-> hex",
            "six triangles per hexagon / six around a triangular vertex",
            [(f"n={n}", 6 * n * n) for n in range(1, 6)],
            "support-six bridge, compatible with octahedral/current language",
        ),
        (
            "hex centered rings",
            "ordinary centered-hex patch count",
            [(f"r={r}", 1 + 3 * r * (r + 1)) for r in range(1, 7)],
            "center-plus-six-petals lattice patch; diverges from 7^k after the first term",
        ),
        (
            "hex norm/self recursion",
            "Eisenstein norm N(a+b*omega)=a^2-a*b+b^2, N(3+omega)=7",
            [(f"7^{k}", 7**k) for k in range(1, 5)],
            "norm-index recursion; relevant to THM-572 only after a labelled state lift",
        ),
    ]


def johnson_defect_samples() -> list[tuple[str, list[tuple[str, Config, int]]]]:
    """Diagnostic nonuniform regular-faced examples, not all 92 Johnson solids."""

    return [
        ("square pyramid J1", [("apex", (3, 3, 3, 3), 1), ("base", (3, 3, 4), 4)]),
        ("pentagonal pyramid J2", [("apex", (3, 3, 3, 3, 3), 1), ("base", (3, 3, 5), 5)]),
        ("triangular prism as Johnson-like check", [("all vertices", (3, 4, 4), 6)]),
    ]


def enumerate_lrc_resonant_uniforms(max_p: int = 18) -> list[UniformCarrier]:
    found: dict[Config, UniformCarrier] = {}
    for degree in range(3, 6):

        def rec(start: int, remaining: int, prefix: list[int]) -> None:
            if remaining == 0:
                config = tuple(prefix)
                kappa = curvature(config)
                vertices = vertices_from_curvature(config)
                faces = face_counts(config)
                if kappa > 0 and vertices is not None and vertices.denominator == 1 and faces:
                    if all(v.denominator == 1 for v in faces.values()):
                        if vertices in (14, 28) or kappa in (Fraction(1, 7), Fraction(1, 14)):
                            found[config] = UniformCarrier(
                                f"uniform cfg {config}",
                                config,
                                "resonant uniform search",
                                "matches LRC 14/28 or curvature 1/7,1/14",
                            )
                return
            for p in range(start, max_p + 1):
                prefix.append(p)
                rec(p, remaining - 1, prefix)
                prefix.pop()

        rec(3, degree, [])
    return sorted(found.values(), key=lambda c: (c.vertices or Fraction(999), c.config))


def score_tournament(items: list[tuple[str, tuple[int, ...]]]) -> tuple[dict[int, int], int, list[str]]:
    wins = Counter()
    edges: dict[tuple[int, int], int] = {}
    for i, j in combinations(range(len(items)), 2):
        winner = i if items[i][1] >= items[j][1] else j
        wins[winner] += 1
        edges[(i, j)] = winner

    def beats(a: int, b: int) -> bool:
        key = (a, b) if a < b else (b, a)
        return edges[key] == a

    c3 = 0
    for a, b, c in combinations(range(len(items)), 3):
        if beats(a, b) and beats(b, c) and beats(c, a):
            c3 += 1
        if beats(a, c) and beats(c, b) and beats(b, a):
            c3 += 1
    order = [items[i][0] for i in sorted(range(len(items)), key=lambda idx: items[idx][1], reverse=True)]
    hist = Counter(wins[i] for i in range(len(items)))
    return dict(sorted(hist.items())), c3, order


def print_regular_map_boundary() -> None:
    print()
    print("[2] Regular map curvature boundary")
    print("  rule:")
    print("    (p-2)(q-2)<4 spherical; =4 Euclidean; >4 hyperbolic.")
    print("  Platonic finite maps:")
    for item in PLATONIC_MAPS:
        print(
            f"    {item.name:12s} {{{item.p},{item.q}}} kind={regular_map_kind(item.p, item.q):10s} "
            f"dual={item.dual:16s} V,E,F={item.counts} kappa={fmt_frac(item.kappa)}"
        )
    print("  regular Euclidean tilings:")
    for item in REGULAR_TILING_MAPS:
        print(
            f"    {item.name:18s} {{{item.p},{item.q}}} kind={regular_map_kind(item.p, item.q):10s} "
            f"dual={item.dual:18s} kappa={fmt_frac(item.kappa)}"
        )
    print("  readout:")
    print("    Platonic solids are positive-curvature finite skeletons.  The plane")
    print("    regular tilings are the zero-curvature recursion wall, not the same")
    print("    scalar object.")


def print_uniform_table(title: str, carriers: list[UniformCarrier]) -> None:
    print()
    print(title)
    for carrier in carriers:
        print(
            f"  {carrier.name:32s} cfg={carrier.config!s:20s} "
            f"kappa={fmt_frac(carrier.kappa):>6s} V={fmt_frac(carrier.vertices):>6s} "
            f"faces={fmt_faces(carrier):22s} role={carrier.role}"
        )


def print_recursions() -> None:
    print()
    print("[6] Euclidean regular tiling recursion labels")
    for label, model, rows, role in recursion_data():
        print(f"  {label}:")
        print(f"    model={model}")
        print("    values=" + ", ".join(f"{tag}:{value}" for tag, value in rows))
        print(f"    role={role}")
    print("  guardrail:")
    print("    Hexagonal 7,49,... is an Eisenstein norm-index chain.  Centered")
    print("    hexagonal rings are 7,19,37,... .  They agree at 7 and then split.")


def print_johnson_guardrail() -> None:
    print()
    print("[7] Johnson-solid guardrail: finite nonuniform defect packets")
    print("  count=92, with regular faces but mixed vertex figures.")
    for name, rows in johnson_defect_samples():
        total = Fraction(0, 1)
        print(f"  {name}:")
        for label, config, multiplicity in rows:
            kappa = curvature(config)
            total += multiplicity * kappa
            print(
                f"    {label:12s} mult={multiplicity:<2d} cfg={config!s:16s} "
                f"kappa={fmt_frac(kappa)} contribution={fmt_frac(multiplicity * kappa)}"
            )
        print(f"    total_curvature={fmt_frac(total)}")
    print("  readout:")
    print("    Johnson solids are finite residual atlases.  They preserve local")
    print("    regular-face defect packets, not a global transitive recursion law.")

def print_annulus_family() -> None:
    print()
    print("[8] Prism / antiprism annular families")
    rows: list[UniformCarrier] = []
    for n in range(3, 17):
        rows.append(prism(n))
        rows.append(antiprism(n))
    for carrier in rows:
        mark = ""
        if carrier.vertices == 14 or carrier.kappa == Fraction(1, 7):
            mark = "  <-- 14 vertices / seven-sector defect"
        if carrier.vertices == 28 or carrier.kappa == Fraction(1, 14):
            mark = "  <-- 28 vertices / LRC 1/14 defect"
        print(
            f"  {carrier.name:22s} cfg={carrier.config!s:16s} "
            f"kappa={fmt_frac(carrier.kappa):>5s} V={fmt_frac(carrier.vertices):>3s}{mark}"
        )
    print("  readout:")
    print("    n-gonal prisms (4,4,n) and antiprisms (3,3,3,n) always have V=2n")
    print("    and per-vertex kappa=1/n.  Thus n=14 gives a 28-point annular")
    print("    carrier at exactly the LRC 1/14 scale.")


def print_resonant_uniform_search() -> None:
    print()
    print("[9] Uniform configurations matching the LRC 14/28 or 1/7,1/14 scales")
    for carrier in enumerate_lrc_resonant_uniforms():
        print(
            f"  cfg={carrier.config!s:18s} kappa={fmt_frac(carrier.kappa):>5s} "
            f"V={fmt_frac(carrier.vertices):>3s} faces={fmt_faces(carrier)}"
        )
    print("  readout:")
    print("    The clean uniform matches are exactly the prism/antiprism annuli.")
    print("    That makes the 14/28 carrier cyclic and branch-local, unlike the")
    print("    q=3 unital which is pair-unique but not cyclic.")


def print_proposed_split() -> None:
    print()
    print("[10] Proposed LRC14 carrier split")
    print("  q=3 unital:")
    print("    28 points; pair-unique incidence; detects repeated pairs such as")
    print("    the H12 conflict in HYP-2942.")
    print("  14-gonal prism/antiprism:")
    print("    28 vertices; two 14-cycles; per-vertex defect 1/14; cyclic annular")
    print("    order for apex clocks and half-step twists.")
    print("  Euclidean tilings:")
    print("    zero-defect recursion labels: square/Gaussian, triangular/Eisenstein,")
    print("    triangle-hex support-six, and hex norm-7.")
    print("  Platonic/Archimedean solids:")
    print("    finite uniform positive-curvature and one-word local quotient seeds.")
    print("  Johnson solids:")
    print("    finite mixed-vertex residual atlases after uniform carriers fail.")


def print_tournament_analysis() -> None:
    print()
    print("[11] Tournament Analysis")
    carriers = [
        ("exact M/Farey branch", (8, 8, 7, 7, 8)),
        ("C27 shell transfer", (7, 8, 6, 6, 7)),
        ("q=3 unital pair incidence", (7, 7, 8, 6, 7)),
        ("14-prism/antiprism annulus", (6, 7, 5, 8, 7)),
        ("Euclidean tiling norm recursions", (4, 6, 5, 6, 6)),
        ("triangle/hex support-six bridge", (4, 5, 5, 6, 6)),
        ("square Gaussian self recursion", (4, 5, 4, 5, 5)),
        ("Platonic/Archimedean uniform atlas", (3, 4, 5, 5, 5)),
        ("hex norm-7 recursion", (3, 4, 4, 5, 6)),
        ("Johnson finite defect atlas", (2, 3, 5, 5, 4)),
        ("raw solid numerology", (0, 1, 1, 1, 1)),
    ]
    hist, c3, order = score_tournament(carriers)
    print("  vertices: proof carriers, not raw solids.")
    print("  observable:")
    print("    exact-LRC retention, C27/branch retention, pair-incidence retention,")
    print("    cyclic/defect retention, anti-scalar guard.")
    for name, score in sorted(carriers, key=lambda item: item[1], reverse=True):
        print(f"    {name:38s} score={score}")
    print(f"  fingerprint: score_hist={hist} c3={c3} hp=1")
    print("  Hamiltonian path:")
    print("    " + " > ".join(order))


def main() -> None:
    print("S141 LRC14 REGULAR-SOLID / EUCLIDEAN-TILING RECURSION CARRIER")
    print("=" * 78)
    print("[1] Assumption challenge / quotient declaration")
    print("  considered vertices:")
    print("    regular maps {p,q}, tiling vertex figures, Euclidean recursion")
    print("    indices, Platonic/Archimedean configurations, prism/antiprism")
    print("    annuli, Johnson defect packets, q=3 unital points, C27 shells, and")
    print("    proof obligations.")
    print("  chosen predicate:")
    print("    normalized angle defect kappa = 1 - sum interior_angle/(2*pi),")
    print("    plus labelled self/dual recursion indices.")
    print("  quotient preserves:")
    print("    local curvature/defect, duality type, norm/self-recursion label,")
    print("    finite-atlas status, and whether a carrier reaches the LRC 1/14")
    print("    defect scale.")
    print("  quotient destroys:")
    print("    exact LRC time geometry, exact denominator witnesses, and pair")
    print("    incidence unless the q=3 unital layer is separately attached.")
    print("  challenged assumption:")
    print("    Solids and tilings should not be treated as one scalar analogy.")
    print("    The audit splits them into curvature walls, norm recursions, annular")
    print("    14/28 carriers, and finite Johnson residual packets.")

    print_regular_map_boundary()
    print_uniform_table("[3] Platonic positive-curvature seeds", PLATONIC)
    print_uniform_table("[4] Archimedean finite uniform seeds", ARCHIMEDEAN)
    print_uniform_table("[5] Euclidean zero-curvature regular and uniform tilings", REGULAR_TILINGS + EUCLIDEAN_UNIFORM_TILINGS)
    print_recursions()
    print_johnson_guardrail()
    print_annulus_family()
    print_resonant_uniform_search()
    print_proposed_split()
    print_tournament_analysis()

    print()
    print("[12] Verdict / next POKE task")
    print("  Norm recursion labels help classify residual packets, but the new")
    print("  proof-facing geometric carrier is the 14-gonal prism/antiprism annulus:")
    print("    V=28 and kappa=1/14.")
    print("  It should be tested as the cyclic companion to the q=3 unital")
    print("  pair-incidence frame.")
    print("  Next concrete experiment:")
    print("    build a 28-vertex annular label model with two 14-cycles, attach")
    print("    AP/GW/H/D labels, and compare whether the H12/GW/K33 conflict seen")
    print("    by the unital becomes a twist, a diameter, or a two-chart obstruction.")


if __name__ == "__main__":
    main()
