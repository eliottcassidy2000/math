# Independent audit: complete color holdouts and additive edge transport

Status: **PASS**, with the finite-carrier and field/group distinctions below.
Source manifest SHA256:
6b476fd57a203f90fc1aefcd51acd043fd549dc689630d15a513711a05ab048f.
All 17 frozen source files match their recorded sizes and digests. No source
file was changed by this audit.

## 1. Full holdout completeness was independently checked

This audit wrote and ran a second C++ enumerator, count.cpp. It constructs
adjacency by iterating square values and complementary summands, then explores
every unused neighbor from every start. It counts ALL oriented Hamiltonian
paths and divides by two, rather than importing the source enumerator or its
completed-path orientation filter. No pruning is used.

The independent counts are:

- Q33: 1338 oriented, 669 unoriented paths;
- Q34: 8082 oriented, 4041 unoriented paths;
- Q35: 34350 oriented, 17175 unoriented paths;
- Q36: 25920 oriented, 12960 unoriented paths.

Each provided path was separately checked for every vertex exactly once,
square adjacency, smaller-first endpoint orientation, and uniqueness.
Because the valid unique provided sets have exactly the independently counted
cardinalities, they are complete. Thus this audit does establish independent
full-census verification for Q33--Q36, not merely imported-path validation.
The Q32 training completeness remains inherited; its path objects were checked.

## 2. Fixed grammar and every score recomputed

The six vertex features, seven edge features, F2/F3 coefficient sets and
global color-renaming quotient match the earlier frozen square-color packet.
The convention is the ordinary Zeckendorf expansion with F2=1. Deduplication
is on the fixed order-80 vertex/edge domain, not a claim that formulas define
the same coloring at all larger integers.

An independent Python scorer regenerates all 1584 distinct candidates,
selects exactly the 1549 Q32 successes, and computes every exact distinct-label
count on Q32--Q36. It agrees with every ledger entry. All 286 failed candidates
have their actual recorded collision pairs checked, and the 16 separately
expanded representative pairs are also validated.

The surviving counts are exactly:

- vertex2: 4 of 56;
- edge2: 15 of 60;
- vertex3: 281 of 359;
- edge3: 963 of 1074.

Total: 1263 of 1549. The first-failure distributions match the source. The
same endpoint pair and the same ordered color word are required; no endpoint
or path-position information is silently removed. Exact string equality,
not equality of hash values alone, determines a collision.

The two-color Zeckendorf-index-sum survivor is a finite positive finding.
Neither its survival nor the large cohort of other survivors proves any
all-size decoder theorem. Pair-only tests at orders36 and80 have genuinely
been superseded at33--36, not at80 or beyond36.

## 3. All-size additive transport proof and exact lost information

For a fixed coloring c with values in an abelian group A, define the edge
color e(uv)=c(u)+c(v). Given c(v_0), the recurrence

    c(v_(i+1))=e(v_i v_(i+1))-c(v_i)

recovers the entire vertex-color word. The forward map is immediate. When
actual starting vertices are retained and c is fixed, this proves equality
of the two collision partitions on any supplied oriented path family.
It does not prove the paths themselves are recoverable unless their vertex
words already separate them. Arbitrary edge-colorings need not have vertex
lifts.

The field kernel statement is also correct. Along edges a kernel vector
alternates signs. On a connected bipartite component its root value is free;
on a non-bipartite component an odd cycle forces twice that value to vanish.
Thus over F2 every component contributes one dimension, while over F3 only
bipartite components contribute one. Isolated vertices are counted as
bipartite components. The source's independent matrix-rank checker was
replayed normally and under -O, byte-matching both source receipts.

### Audited extension: the common abelian-group universe

The same proof works over ANY abelian group A, without a field hypothesis.
For a connected bipartite component the kernel is isomorphic to A, with
alternating signs on its two parts. For a connected non-bipartite component
the kernel is isomorphic to A[2]={a in A:2a=0}; its assignments are constant,
since a=-a for a in A[2]. Different components contribute independently.

For finite A and a graph with b bipartite components among k total components,

    |kernel| = |A|^b |A[2]|^(k-b).

For A=Z/qZ this is q^b gcd(q,2)^(k-b). This places the two-/three-color
discrepancy in one precise universe: it is surviving 2-torsion, not an
analogy based solely on palette cardinality. The extension is proved by the
same sign propagation, with additional exhaustive tests over Z/4Z and Z/6Z
on every simple graph through four vertices.

## 4. Audit controls and reproduction

Independent check.py validates all path objects and complete color scores,
then checks incidence-kernel counts on all 1099 labeled simple graphs through
five vertices over F2 and F3, and reconstructs every color word of length at
most eight from its first color and additive edge word. Normal and optimized
runs each pass 88002 explicit checks and byte-match. The cyclic 4-/6-color
torsion extension has 150 further byte-matching checks.

To reproduce after publication, first run the source color-holdouts
verify.py --fresh to regenerate its intentionally omitted path files. Then:

    g++ -O3 -std=c++17 count.cpp -o /tmp/s42-count
    /tmp/s42-count > independent-counts.json
    python3 check.py
    python3 -O check.py
    python3 check_pins_witnesses.py
    python3 check_torsion_extension.py
    python3 -O check_torsion_extension.py

The audit expects its sibling source folder named color-holdouts (published)
or session42-color-holdouts (workspace). The count-only program does not need
the generated source path files. No giant path inventory is duplicated in
this audit packet.
