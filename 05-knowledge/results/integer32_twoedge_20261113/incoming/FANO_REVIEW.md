# Incoming Fano lift: independent literal replay and weight audit

Status: **AUDIT PASS / FINITE-EXACT** for the declared128-twist family,
canonical partition, majority/port structure, and orbit classification.
The metagraph degree convention and regular-source obstruction pass the
source/proof review below. Large Hamilton-path counts are **INHERITED,
NOT INDEPENDENTLY REPLAYED**. No all-tournament classification follows.

## Source pin and independent path

Reviewed the four cached files in `/tmp/integer32-incoming/` against their
exact bytes at commit `0eab57c4327b6056480f1943bd2629a32dc6e969`, under
`05-knowledge/results/session321_selective_returns_and_fano_wiring_20261113/`:
FANO_METAGRAPH.md, BSD_WIRING_STENCIL.md, BSD_PARITY.md, WEDGE_WIRING.md.
[SOURCE_PINS.json](SOURCE_PINS.json) records their Git blobs and SHA256 hashes.
The review did not require a worktree rebase or import incoming producer code.

[check_fano_twists.py](check_fano_twists.py) constructs adjacency directly from

    (a,b)(c,d)=(a+2^b c,b+d),
    S0={1,2,3}, S1={0,1,2,4}, S2={1,2,4},
    x -> y iff x^(-1)y belongs to S.

Vertices have index3a+b, with a modulo7 and b modulo3. It reverses exactly
the three internal arcs of triangle F_a when bit a of the mask is1. It uses
no Hamilton-path DP, imported canonicalizer, or incoming implementation.
Normal and optimized runs are byte-identical and pass **72,243 explicit
checks**. The finite universe is all128 masks, all pairs in each tournament,
all5040 quotient permutations, and all21 ordinary plus21 converse candidate
lifts for every mask.

## What the literal replay establishes

Every twist is a21-vertex tournament with outdegree10 at every vertex.
Forward and reverse reachability from one vertex prove strong connectivity.
Every internal triangle remains a directed3-cycle. For every unordered pair
in every twist, common-in and common-out counts agree; common-out count7
holds **exactly** for pairs in the same F_a. Cross pairs have maximum count5.
Thus the triangle partition is intrinsic under both isomorphism and converse.

Between F_a and F_c with c−a in D={1,2,4}, exactly eight arcs point from a
to c and one points back. The unique dissent joins port b to the same port b,
where c−a=2^(b+2) modulo7. All three outgoing quotient arcs select distinct
ports. Hence the majority quotient is the fixed Paley7 tournament, and its
dissent incidences distinguish every vertex within each triangle.

The checker also constructs F8 independently as F2[X]/(X³+X+1). Its seven
translated exponent triples a+D are exactly seven zero-sum field triples,
with a unique such triple through every pair. This verifies the claimed
Fano interpretation without identifying the exponent space F7 with the
additive field F8.

### Why these are full isomorphism classes, not just action orbits

Any isomorphism preserves the intrinsic triangles and induces a Paley quotient
automorphism. Exhausting all5040 quotient permutations gives precisely the21
maps a→t+2^k a. A quotient map has at most one lift: its action on each
outgoing dissent arc forces each of the three ports. The candidate lift

    (a,b) -> (t+2^k a,b+k)

is checked on all210 vertex pairs for every mask and every t,k. It carries a
mask to its affine image. Thus the mask stabilizers are the **full** tournament
automorphism groups in this family.

Similarly every quotient anti-map is a→t−2^k a, with forced port map

    (a,b) -> (t−2^k a,b+k).

Its target mask is the complement of the image mask. Literal adjacency
checks verify all of these anti-isomorphisms. Because the point universe
has odd size7, mask sizes k and7−k never agree. Consequently **none** of the
128 twists is self-converse. The argument is about the odd size of the
seven-point universe; it does not assume that each individual mask has
odd cardinality.

The ordinary and merged representatives are:

    ordinary: 0,1,3,7,11,13,15,23,29,31,63,127  (12 classes)
    merged:   0,1,3,7,11,13                    (6 classes).

The ordinary action has cycle types1^7 once,7 six times, and1·3·3 fourteen
times. Its direct fixed-mask sum agrees with Burnside's12. At mask size3,
representatives7,11,13 are respectively a transitive triple, a Fano-line
directed triangle, and a nonline directed triangle, with orbit sizes21,7,7.
These exhaust all35 triples.

## Metagraph incidence and the factor of two

The incoming weight statement agrees with **THM-4731**, §3,
`01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md`.
For a non-self-converse class T:

    ordinary marked-path orbit count = H(T)/|Aut(T)|,
    fundamental black degree         = H(T)/|Aut(T)|,
    full merged black degree          = 2 H(T)/|Aut(T)|.

The automorphism action on Hamilton words is free: fixing a full word fixes
every vertex. The two ordinary converse classes each supply H/|Aut| marked
incidences. Their sigma twins become a pair of black lines with the same
merged endpoints; quotienting the twin pair gives the fundamental black
convention. No blue presentation can have this chiral source. The source
therefore uses H/21 correctly for untwisted U; **2H/21** is the corresponding
full merged black degree. Neither degree counts distinct adjacent classes.

The regular-source obstruction is independent of the large H values. If T
has order n=2r+1 and outdegree r, its path-complement partner has outdegrees
r+1 at the path start, r−1 at its end, and r elsewhere. This follows by
reversing all nonpath arcs while retaining the directed path. At n=21 these
are11,9,10. Consequently no partner is regular, even after converse merging.
All six Fano merged vertices are therefore loop-free and mutually nonadjacent
in the actual complement-line metagraph. This conclusion needs no Hamilton
count or assumption that different degree totals determine graph wiring.

The arithmetic-only checks of the inherited degree table give:

| Mask | Checked automorphism order | Inherited H, not recomputed | Fundamental black degree | Full merged black degree |
|---:|---:|---:|---:|---:|
|0|21|123376391895435|5875066280735|11750132561470|
|1|3|123381234928431|41127078309477|82254156618954|
|3|1|123384468287435|123384468287435|246768936574870|
|7|1|123386086601285|123386086601285|246772173202570|
|11|3|123385527546237|41128509182079|82257018364158|
|13|3|123386644838787|41128881612929|82257763225858|

Here |Aut| is the full order justified by the canonical-quotient/port proof.
The checker verifies divisibility and sums the ten incoming first-arc values
to H(U)/21. Those arithmetic equalities do **not** independently establish
H(U), any first-arc count, or distinct Hamilton spectra. The six structural
classes above are established independently of these inherited numbers.

## The two Fano families are different objects

This incoming family has **128 masks on one fixed seven-point Fano/Paley
quotient**, modifying internal orientations in a21-vertex port lift. The
current bridge lane's other family consists of **30 labelled Fano line
systems/simplex codes on seven labels**. The latter comes from relabelling
the incidence structure and its pairwise-blind, triple-faithful character
statistics. These have different states, actions, and quotient operations.
The shared named Fano incidence supports comparisons; it does not identify
128 with30 or transfer their automorphism groups, orbit counts, or H values.

## Companion algebra: bounded proof review

The other incoming texts supply compatible constraints, not a second proof
of the orbit classification.

- **BSD_PARITY.md:** the displayed selective-attachment block matrix retains
  the actual cut c and its complement r. Over F4, the TT_m block minus lambda
  is invertible because its diagonal lies in F2 while lambda does not. Its
  inverse on the all-ones vector is (eps−lambda)^(-1), so the Schur complement
  depends on m only through parity. Squarefreeness of X²+X+1 and Frobenius
  conjugacy justify the doubled eigenspace dimension. The algebra passes;
  the listed large-family nullities and arithmetic realizations were not
  replayed in this checker.
- **BSD_WIRING_STENCIL.md:** with row-zero A and binary undirected path
  Laplacian L, B=A+J+(n mod2)I+L gives
  Phi(B)+Phi(A)=AL+LA+L²+L. The converse identity cancels J, and JL=LJ=0.
  Symmetry uses A^T=A+J+I; the diagonal vanishes because diag(AL+LA)=diag(L)
  and diag(L²)=0. The interior distance-at-least-three stencil follows;
  endpoint and short-diagonal exceptions must retain the full formula.
  This is a necessary aligned-line constraint, not a nullity-preservation
  law or a complete wiring criterion.
- **WEDGE_WIRING.md:** the signed edge disagreement integrates from wedge
  differences up to complement on connected L(K_n). Testing either support
  for a coherent spanning directed path supplies the exact aligned line
  criterion. Its uniqueness bound n≥5 is correct; at n=4 complementary
  three-edge paths may both pass. The ordinary alignment count divides by
  both source and target automorphism orders, with the actual polarity fixed.
  A converse-merged count must aggregate distinct ordinary fibre elements,
  never double a self-converse singleton. This agrees with the source-weight
  distinction above and with the current reconstruction lane.

These are source/proof audits of the displayed finite algebra. No claim
about elliptic-curve ranks, BSD, Bott periodicity, or a general tournament
spectrum is imported from the common numbers or fields.

## Reproduction and stopping boundary

From the repository root:

```sh
python3 05-knowledge/results/integer32_twoedge_20261113/incoming/check_fano_twists.py
python3 -O 05-knowledge/results/integer32_twoedge_20261113/incoming/check_fano_twists.py
```

The source and both receipts are pinned in SHA256SUMS. Checks stay active
under `-O`. No expensive Hamilton DP was run. The precise remaining dependency
for the table's numerical H and first-arc values is the incoming source's
independent full-start DP, not this small replay.
