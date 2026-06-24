#!/usr/bin/env python3
"""S141: regular-solid and Euclidean tiling recursion carriers for LRC14.

The prompt asks for Platonic, Archimedean, Johnson-solid, and Euclidean tiling
recursion signals.  This script keeps the geometric data labelled:

* vertex curvature/defect of regular polygon vertex configurations;
* self/dual recursion counts for square, triangular, and hexagonal tilings;
* prism/antiprism annular families, where n=14 gives 28 vertices and
  per-vertex curvature 1/14;
* a finite Johnson-defect guardrail showing why nonuniform regular-faced solids
  belong to finite-atlas work rather than a single global recursion law;
* Tournament Analysis on proof carriers.

The output is not an LRC14 proof.  It is a structural carrier audit.
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


def fmt_frac(x: Fraction | None) -> str:
    if x is None:
        return "flat/negative"
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


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


PLATONIC = [
    UniformCarrier("tetrahedron", (3, 3, 3), "Platonic", "triangle self seed / four-face closure"),
    UniformCarrier("cube", (4, 4, 4), "Platonic", "square self-dual solid branch"),
    UniformCarrier("octahedron", (3, 3, 3, 3), "Platonic", "triangle/hex six-vertex hinge"),
    UniformCarrier("dodecahedron", (5, 5, 5), "Platonic", "pentagonal face branch"),
    UniformCarrier("icosahedron", (3, 3, 3, 3, 3), "Platonic", "triangular high-valence branch"),
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

REGULAR_TILINGS = [
    UniformCarrier("triangular tiling", (3, 3, 3, 3, 3, 3), "Euclidean tiling", "flat six-triangle vertex"),
    UniformCarrier("square tiling", (4, 4, 4, 4), "Euclidean tiling", "flat self-dual square grid"),
    UniformCarrier("hexagonal tiling", (6, 6, 6), "Euclidean tiling", "flat dual of triangular grid"),
]

EUCLIDEAN_UNIFORM_TILINGS = [
    UniformCarrier("trihexagonal tiling", (3, 6, 3, 6), "Euclidean uniform tiling", "triangle/hex alternating"),
    UniformCarrier("truncated square tiling", (4, 8, 8), "Euclidean uniform tiling", "square-octagon flat branch"),
    UniformCarrier("truncated hexagonal tiling", (3, 12, 12), "Euclidean uniform tiling", "triangle-dodecagon flat branch"),
    UniformCarrier("rhombitrihexagonal tiling", (3, 4, 6, 4), "Euclidean uniform tiling", "tri-square-hex flat branch"),
    UniformCarrier("snub square tiling", (3, 3, 3, 3, 4), "Euclidean uniform tiling", "flat chiral square-triangle limit"),
    UniformCarrier("elongated triangular tiling", (3, 3, 3, 4, 4), "Euclidean uniform tiling", "triangle-square flat branch"),
]


def prism(n: int) -> UniformCarrier:
    return UniformCarrier(f"{n}-gonal prism", (4, 4, n), "prism family", "two n-cycles plus square belt")


def antiprism(n: int) -> UniformCarrier:
    return UniformCarrier(
        f"{n}-gonal antiprism", (3, 3, 3, n), "antiprism family", "twisted two n-cycles plus triangle belt"
    )


def self_recursion_data() -> dict[str, list[tuple[str, int]]]:
    out: dict[str, list[tuple[str, int]]] = {
        "square self n^2": [(f"n={n}", n * n) for n in range(2, 8)],
        "triangle self 4^k": [(f"k={k}", 4**k) for k in range(1, 6)],
        "triangle->hex exchange 6*n^2": [(f"n={n}", 6 * n * n) for n in range(1, 6)],
        "hex centered patch 1+3r(r+1)": [(f"r={r}", 1 + 3 * r * (r + 1)) for r in range(1, 7)],
        "hex self hexaflake 7^k": [(f"k={k}", 7**k) for k in range(1, 5)],
    }
    return out


def johnson_defect_samples() -> list[tuple[str, list[tuple[str, Config, int]]]]:
    """Small nonuniform regular-faced examples with exact vertex defects.

    These are not an enumeration of all 92 Johnson solids.  They are diagnostic
    examples showing that Johnson objects are finite defect atlases, not single
    vertex-transitive recursion laws.
    """

    return [
        ("square pyramid J1", [("apex", (3, 3, 3, 3), 1), ("base", (3, 3, 4), 4)]),
        ("pentagonal pyramid J2", [("apex", (3, 3, 3, 3, 3), 1), ("base", (3, 3, 5), 5)]),
        ("triangular prism as Johnson-like check", [("all vertices", (3, 4, 4), 6)]),
    ]


def print_uniform_table(title: str, carriers: list[UniformCarrier]) -> None:
    print()
    print(title)
    for carrier in carriers:
        faces = carrier.faces
        face_text = "{}" if faces is None else "{" + ", ".join(f"{p}:{fmt_frac(v)}" for p, v in faces.items()) + "}"
        print(
            f"  {carrier.name:32s} cfg={carrier.config!s:20s} "
            f"kappa={fmt_frac(carrier.kappa):>6s} V={fmt_frac(carrier.vertices):>6s} "
            f"faces={face_text:22s} role={carrier.role}"
        )


def face_text(carrier: UniformCarrier) -> str:
    faces = carrier.faces
    if faces is None:
        return "{}"
    return "{" + ", ".join(f"{p}:{fmt_frac(v)}" for p, v in faces.items()) + "}"


def enumerate_lrc_resonant_uniforms(max_p: int = 16) -> list[UniformCarrier]:
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


def print_johnson_guardrail() -> None:
    print()
    print("[6] Johnson-solid guardrail: finite nonuniform defect packets")
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
    print("    Johnson solids should be treated like finite residual atlases.  They")
    print("    preserve regular-face local defects, but not a single transitive")
    print("    vertex configuration or global recursion law.")


def main() -> None:
    print("S141 LRC14 REGULAR-SOLID / EUCLIDEAN-TILING RECURSION CARRIER")
    print("=" * 78)
    print("[1] Assumption challenge / quotient declaration")
    print("  considered vertices:")
    print("    regular tiling cells, tiling vertex figures, Platonic/Archimedean")
    print("    vertex configurations, prism/antiprism annuli, Johnson defect packets,")
    print("    q=3 unital points, C27 shells, and LRC proof obligations.")
    print("  chosen predicate:")
    print("    normalized angle defect kappa = 1 - sum interior_angle/(2*pi),")
    print("    plus declared self/dual recursion counts.")
    print("  quotient preserves:")
    print("    local curvature/defect, duality type, self-recursion count, finite")
    print("    atlas status, and whether a carrier has the LRC 1/14 defect scale.")
    print("  quotient destroys:")
    print("    exact LRC time geometry, exact denominator witnesses, and pair")
    print("    incidence unless the q=3 unital layer is separately attached.")
    print("  challenged assumption:")
    print("    Platonic/Archimedean/Johnson solids form one global atlas.  The audit")
    print("    splits them into uniform curvature laws, annular 14/28 carriers, and")
    print("    finite nonuniform Johnson defect packets.")

    print_uniform_table("[2] Platonic positive-curvature seeds", PLATONIC)
    print_uniform_table("[3] Archimedean finite uniform seeds", ARCHIMEDEAN)
    print_uniform_table("[4] Euclidean zero-curvature regular tilings", REGULAR_TILINGS)

    print()
    print("[5] Prompt recursion counts for the three Euclidean regular tilings")
    for label, rows in self_recursion_data().items():
        print(f"  {label:32s} " + ", ".join(f"{tag}:{value}" for tag, value in rows))
    print("  readout:")
    print("    square recursion is n^2/self-dual; triangle self-recursion is dyadic")
    print("    4^k; triangle<->hex exchange is 6*n^2; hex self-recursion has two")
    print("    useful readings: centered lattice patches 1+3r(r+1), and hexaflake")
    print("    self-similar packets 7^k.")

    print_johnson_guardrail()

    print()
    print("[7] Prism / antiprism annular families")
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
    print("    carrier at exactly the LRC 1/14 scale, complementary to the q=3")
    print("    unital's 28-point pair-incidence carrier.")

    print()
    print("[8] Uniform configurations matching the LRC 14/28 or 1/7,1/14 scales")
    for carrier in enumerate_lrc_resonant_uniforms():
        print(
            f"  cfg={carrier.config!s:18s} kappa={fmt_frac(carrier.kappa):>5s} "
            f"V={fmt_frac(carrier.vertices):>3s} faces={face_text(carrier)}"
        )
    print("  readout:")
    print("    The clean uniform matches are exactly the prism/antiprism annuli.")
    print("    That makes the 14/28 carrier cyclic and branch-local, unlike the")
    print("    q=3 unital which is pair-unique but not cyclic.")

    print()
    print("[9] Proposed LRC14 carrier split")
    print("  q=3 unital:")
    print("    preserves pair incidence; 28 points; detects repeated pairs such as")
    print("    the H12 conflict in HYP-2942.")
    print("  14-gonal prism/antiprism:")
    print("    preserves cyclic annular order; 28 vertices; per-vertex defect 1/14;")
    print("    candidate carrier for apex clocks, two 14-cycles, and half-step twists.")
    print("  Euclidean tilings:")
    print("    preserve zero-defect recursion rules; useful as flat limits and")
    print("    self/dual recursion mnemonics, not as finite LRC certificates.")
    print("  Johnson solids:")
    print("    preserve finite nonuniform defect packets; useful as finite residual")
    print("    atlases after uniform/prism carriers fail.")

    print()
    print("[10] Tournament Analysis")
    carriers = [
        ("exact M/Farey branch", (7, 7, 6, 6, 7)),
        ("C27 shell transfer", (6, 7, 5, 5, 6)),
        ("q=3 unital pair incidence", (6, 6, 7, 5, 6)),
        ("14-prism/antiprism annulus", (5, 6, 4, 7, 6)),
        ("Euclidean tiling zero-defect laws", (3, 5, 3, 6, 5)),
        ("triangle/hex dual exchange 6", (3, 4, 3, 5, 5)),
        ("square self n^2 recursion", (3, 4, 2, 4, 4)),
        ("hex self 7^k recursion", (2, 3, 2, 4, 5)),
        ("Platonic/Archimedean uniform atlas", (2, 3, 4, 4, 4)),
        ("Johnson finite defect atlas", (1, 2, 4, 5, 4)),
        ("raw solid numerology", (0, 1, 1, 1, 1)),
    ]
    hist, c3, order = score_tournament(carriers)
    print("  vertices: proof carriers, not raw solids.")
    print("  observable:")
    print("    exact-LRC retention, C27/branch retention, pair-incidence retention,")
    print("    cyclic/defect retention, anti-scalar guard.")
    for name, score in sorted(carriers, key=lambda item: item[1], reverse=True):
        print(f"    {name:36s} score={score}")
    print(f"  fingerprint: score_hist={hist} c3={c3} hp=1")
    print("  Hamiltonian path:")
    print("    " + " > ".join(order))

    print()
    print("[11] Verdict / next POKE task")
    print("  The new high-leverage carrier is the 14-gonal prism/antiprism annulus:")
    print("    V=28 and kappa=1/14, matching both the q=3 unital point count and")
    print("    the LRC14 threshold.  It should be tested as the cyclic companion to")
    print("    the q=3 unital pair-incidence frame.")
    print("  Next concrete experiment:")
    print("    build a 28-vertex annular label model with two 14-cycles, attach")
    print("    AP/GW/H/D labels, and compare whether the H12/GW/K33 conflict seen")
    print("    by the unital becomes a twist, a diameter, or a two-chart obstruction.")


if __name__ == "__main__":
    main()
