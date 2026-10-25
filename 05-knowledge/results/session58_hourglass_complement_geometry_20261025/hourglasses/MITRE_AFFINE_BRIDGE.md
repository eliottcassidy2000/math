# The magic hourglass is a mitre, and two rooted mitres recover affine magic patterns

Status: CITED classical identification; PROVED elementary counting and arithmetic boundary; FINITE-EXACT independent enumeration. No new Steiner-system existence theorem or integer-square solution is claimed.

## Primary-source identification

A. E. Brouwer, *Steiner triple systems without forbidden subconfigurations*, ZW104/77 (1977), §3, printed pp6–7, defines the five-block seven-point configuration

    {0,1,1′}, {0,2,2′}, {0,3,3′}, {1,2,3}, {1′,2′,3′}

as a mitre and identifies it with AG(2,3) after removing two points and all seven lines incident with either point. The direct renaming 0=E, (1,2,3)=(A,B,C), (1′,2′,3′)=(I,H,G) gives exactly the horizontal magic-hourglass line hypergraph. This nomenclature and deletion interpretation are classical, not discoveries of this session. Sources: https://ir.cwi.nl/pub/6839 and https://ir.cwi.nl/pub/6839/6839D.pdf . PDF text/version inspected; a requested screenshot failed with cache-miss, and no diagram-dependent claim is based on it. No full article is copied into the artifact.

The arithmetic predicate remains extra: the abstract mitre has seven vertices and five blocks regardless of whether any values assigned to those vertices are squares or have equal sums.

## Minimal Steiner host and explicit nine-point completion

A Steiner triple system means every unordered pair of distinct points belongs to exactly one three-point block. The hourglass is a maximal partial such system on its seven points: its six uncovered pairs form a six-cycle and cannot contain another block (THEOREM.md §7). Hence no seven-point Steiner triple system can contain it. An eight-point Steiner triple system cannot exist because each point would require (8−1)/2 blocks. Fewer than seven points cannot hold an injective copy. Thus nine is a lower bound for a Steiner host.

For the matching upper bound, take F3² and all triples of distinct points whose vector sum is0. Given two points u≠v, the third point −u−v is different from both and is uniquely determined, so this is a Steiner triple system with9points and12blocks. In coordinates (column,row), place the nine cells A,...,I at the usual3×3 positions with rows/columns indexed −1,0,1. The five horizontal hourglass blocks are affine lines. Deleting D,F removes precisely the other7lines, so the remaining7points induce the hourglass. Thus9 is the minimum Steiner host order. This is an incidence extension, not a way of manufacturing square labels for D,F.

## Exact census and rooted overlap

There are C(9,2)=36 mitres in AG(2,3), one for each deleted pair. Each pair lies on a unique line; the third point c of that line is the mitre's unique degree3 point, its center. Equivalently, for each of9centers, choose one of the4lines through c and remove the other2points on that line. Four choices per center give36.

Completeness: any mitre has7points, leaving one complementary pair. The plane has exactly5lines avoiding those2points (12−(4+4−1)=5), so a five-line mitre on that support must consist of all those lines. This proves the count without reliance on finite enumeration.

Choose TWO distinct mitres with the same center c. Their missing pairs lie in different directions and are disjoint. They therefore share5points and cover9. The common blocks are the two other lines through c; hence their5+5−2=8-block union is the ordinary magic-square incidence pattern. There are9*C(4,2)=54 unordered same-center pairs, giving54 different eight-line patterns.

Distinctness of the54 patterns follows directly: the center is the unique point of degree4 in the union; the four points of degree2 are the two missing pairs; their lines through the center determine the two omitted directions. Thus the pair is reconstructed from its union.

Each pattern has exactly8 labelled grid layouts, related by D4 symmetries. Independent derivation: every injective vector-valued magic array over F3² has the normal form from THEOREM.md with center c and linearly independent offsets x,z. There are9choices of center and(9−1)(9−3)=48ordered independent offset pairs, giving432layouts. Every layout gives one same-center mitre pair; the8ordinary grid symmetries act freely and account for54patterns. This recovers Session57's432layouts/54D4-orbits by a new explicit mitre-pair interpretation, not by equating two unexplained matching counts.

## Which four affine lines are omitted, and why characteristic3 matters

The full affine plane has12lines; a standard3×3 magic pattern has8. The four additional wraparound diagonals, in cell labels, are

    AFH, BDI, BFG, CDH.

For the general magic parametrization with center E and corner offsets x,z, their sums are respectively

    3E+3x, 3E−3x, 3E−3z, 3E+3z.

Consequently all12affine triples have the original total3E iff3x=3z=0 in the additive target group. Over a torsion-free group this forces x=z=0 and all entries equal. In characteristic3 these extra constraints vanish identically.

This is the exact arithmetic boundary of the finite-geometry connection. An affine-plane host for a seven-cell incidence pattern does NOT imply that a nonconstant rational or integer magic square satisfies its four additional lines. Over characteristic3, however, the entire affine-plane incidence is compatible with the additive magic equations. Square membership and injectivity must still be verified in the chosen extension field; the abstract9points of F3² are not automatically9distinct scalar squares in any field.

## Nonzero-center finite hostile and reproducibility

The first anchor hostile over F13 has center0 and is unsuitable for a chart normalized at nonzero center. A chart-safe control occurs over F37 with E=1,x=8,z=27:

    9   3  28
    20  1  19
    11 36  30.

The seven hourglass values are distinct quadratic residues, all five retained line sums are3mod37, and the two forced omitted values20,19 are nonsquares. This center is already normalized to1. An exhaustive exact search of every prime≤37 confirms that37 is the smallest prime with a nonzero-square center, seven distinct square hourglass entries, and BOTH omitted entries nonsquare. This is explicitly not a minimum over extension fields or a global integer statement.

`python3 check_geometry.py` verifies all36deleted-pair mitres, the54same-center pairs, the432normal-form layouts by independent enumeration, their equality of pattern sets and eight-element fibers, the four omitted-line defects for343integer parameter triples, and the stated prime-field hostile minimality. Normal and optimized Python return byte-identical output. It imports only the local exact definitions from check.py and uses explicit exceptions rather than assertions. geometry-certificate.json records the checks and all small-prime outcomes.

The deliverable adds a genuine7→9 incidence bridge with an exact characteristic3 lifting boundary. It does not transfer Keller inverse-fiber structure or resolve square membership over the integers.
