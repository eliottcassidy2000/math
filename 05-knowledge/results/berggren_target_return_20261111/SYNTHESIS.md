# Absolute returns and a completely classified infinite metagraph family

**PROVED IN NOTES / FINITE-EXACT, with independent internal reviews.
Chan–Pak OPEN.** Here H(T) is the number of directed Hamilton paths of a
tournament. The target conjecture is that the positive odd values other
than7 and21 all occur. This packet advances explicit constructions and
classifies a substantial family; it does not prove that conjecture or
classify the unrestricted metagraph.

## 1. Two uniform returns at the exact requested target

Set D=2^k, k>=1. The inherited physical construction is

    11D+1 →[+4] 33D+7 →[+2] 99D+23,

where [+c] means the count3H+c. The new R_k source has exactly the last
count and an indegree-one root of defect1. Its unique-start parent is a
fixed C3 funnel with an expandable transitive block. The full module-run
polynomial has no third-run contribution, so the identity holds at every
k, not just the first numerical example.

There is an actual same-count replacement from the inherited endpoint
to R_k. In the best displayed orientation gauge it retains k+5 of the
k+7 old vertices, deleting two and inserting three. It keeps the original
root, a designated port and the two recent newborns. The transitive block
is mapped in reverse order into the converse target; subsequent operations
are mirrored accordingly. No reversal of a retained old arc is hidden.

The next standard two births produce297D+71. A second uniform same-count
replacement gives another defect1 source Q_k, retaining k+8 of k+10 current
vertices. Only one deleted vertex predates that excursion. The next +2
step reaches891D+215. At k1 this is the actual chain

    23 → 73 → 221 → 665 → 1997,

with ready returns at221 and665. All maps, root coordinates, literal row
masks, and independent module-polynomial checks are in
[TARGET_RETURN](reset/TARGET_RETURN.md) and its
[audit](berggren/RESET_AUDIT.md).

This is stronger than the earlier repeatable return whose scalar target
drifted. It still covers two uniform return stages, not an indefinitely
closed source class. The complete ready inverse proves that1997 has no
strong1-ready or2-ready representative at any order; inherited neutral
padding does not create a missing ready value. The finite enumeration is
closed by the cited strong-tournament minimum, not by an arbitrary cutoff.

## 2. Indegree one escapes the branching-defect obstruction

The prior bounded-defect theorem concerns roots with at least two incoming
arcs. Its dyadic-ray sparsity does not include indegree-one roots. Keeping
this distinction produces a new source B with H=1997 and defect3, while
preserving the earlier physical root and its sole incoming neighbor.
The same-count return retains ten of thirteen old vertices. Defects1 and2
are impossible at this count, so3 is the smallest positive defect;
defect0 is not excluded by that statement.

The construction belongs to an explicit infinite family

    H(B_m)=1134·2^m−271,   defect(B_m)=3,   m>=1.

More generally, the first vertex of a two-flank core chart beginning OO
has indegree one and defect H(core), independently of the later cut.
This links the original five-core bank to a hierarchy of persistent
indegree-one interfaces, rather than a bounded list of branching roots.

The particular H1997 source now has both desired branches:

| Target | Actual construction and retained data | Further continuation |
|---|---|---|
| 5993=3·1997+2 | Replace two of twelve old vertices and insert three; retain the same root and sole incoming neighbor; new defect5 | Same-count return to a defect2 source, retaining ten of thirteen; then +4→17983 |
| 5995=3·1997+4 | Start with the standard5997 child; five gap transports and one apex-arc flip retain all fourteen names, peak12215 | Same-count defect1 return retains nine names including the same root and original neighbor; then +2→17987 |

The second branch temporarily loses its one-incoming-arc interface and
restores it at the final same-count return. Its intermediate gap moves
are not individually count-neutral in the full tournament.

The first branch uses S4 as its residual core. An exact all-height inverse
proves that5993 is absent from every two-triangle residual chart, including
arbitrary gaps, and also from every S4→C3 or C3→S4 two-flank chart. It is
present in six normalized S4 charts. Thus changing the core is essential
for these particular menus. The four strong defect2 targets at5993 have
fixed-root/port retention optima10,9,8,8 for this source; this is not an
unrestricted optimum over every possible replacement graph.

The attractive single-module uniformization fails for a precise reason.
The displayed defect3→defect5 map has target

    3684·2^m−1375 = 3H(B_m)+2 + 282(2^m−2).

It preserves all required names for every m but the desired absolute
target only at m1. All84 single-vertex transitive inflations of the six
normalized S4 charts fail the necessary slope3402. This is an excluded
construction method, not a no-return theorem.
See [DEFECT_THREE_ESCAPE](reset/DEFECT_THREE_ESCAPE.md) and
[H15_MENU](reset/H15_MENU.md).

## 3. Why several other repeatability attempts fail

For an indegree-one rooted source retain the triple(h,F,s), with defect
b=F−h. The standard two-birth operation is exactly

    (h,F,s) ↦ (h+2F, 2h+2F−s, F).

Its first output count is3h+2b. Starting at Q_k, defect grows from1 to
192D+47 and then594D+143. The next counts are891D+215 and3057D+739;
the latter overshoots the desired2673D+647 by384D+92. The old H13
translation actuator can realize this difference, but its minimum baseline
already exceeds the entire source count. It cannot be installed as an
induced subgraph. The gap alone is not an available operation.

The older growing cap also fails a stronger test. Although its endpoint
mass stays1 and its feedback number grows, its two-cover coordinate obeys
the characteristic polynomial X^3−6X^2−3X−12. Its ready counts begin
1,5,29,173,1061,6821. No fixed unique-start seed followed indefinitely by
this cap has eventual3H+2 recurrence: the required exponential mode3 is
absent. The conserved endpoint coordinate does not control the cover cost.

An independent Berggren-prefix rewrite is positive. For r>=4,t>=1 and
a strong suffix z of count3·2^(r−3)−1, the source word O^r I^t z has a
same-target3H+2 return obtained by reversing r−1 arcs at the newest vertex.
Every old source arc and every vertex survives. In the correct converse
gauge a subfamily continues

    105·2^t−9 →315·2^t−25 →945·2^t−73,

beginning201→605→1817. A finite higher-core ready return then reaches5453.
The immediate same rewrite cannot recur;1817 is absent from the complete
one-apex spectrum. The entire displayed unbounded-prefix family has only
O(sqrt(X)) presentations below X, even after all eligible suffixes are
included. This explains why this useful operation alone does not supply
global coverage. See [WIDE_RETURN](berggren/WIDE_RETURN.md),
[its size proof](fibres/WIDE_REWRITE_SIZE.md), and
[ALTERNATING_AND_CAP](reset/ALTERNATING_AND_CAP.md).

## 4. A complete all-size classification: triangle chains

Let

    T=TT_g0 → C3 → TT_g1 → ... → C3 → TT_gr,

with r triangles and nonnegative singleton gaps. Let nu be the maximum
matching size in the path of triangles, with an edge exactly when their
intervening gap is zero. Equivalently, sum floor(runlength/2) over all
consecutive-triangle runs. Define C_k(T) as the gcd of H(U) over EVERY
induced extension U obtained by adjoining k arbitrary vertices.

The new, independently reviewed all-k theorem is

    C_k(T)=3^max(r−k−min(k,nu),0),
    L(T)=min{k:C_k(T)=1}=r−nu.

The proof groups coefficients of the full Boolean orientation polynomial
into local path-cover products. Untouched triangles contribute factors3.
Contacts with newborns pay half-units toward the possible losses; gaps
and endpoints pay the remaining local terms, and the only unpaid terms
are bounded by the zero-gap matching. Disjoint H5 and H41 constructors
attain the reverse divisibility. This proves the full filtration, not
merely a finite k2 pattern. The complete C2 coefficient controls include
a fifteen-vertex source covering all2^31 two-birth assignments.

An earlier assertion that an untouched triangle must appear consecutively
in a Hamilton path was false. The corrected argument cyclically relabels
its existing positions, preserving the used old-old edges. The concrete
hostile example and correction lineage are in MISTAKES and the replay.

The family is exactly the class of tournaments whose strong components
are K1 or C3. Its ordinary classes are ordered component words, its
automorphism group has order3^r, and its complete joint enumerator is

    F(z,u,v)=(1+u z^3)/(1−z−u z^4−v u^2 z^6),

where z marks vertices, u triangles, and v the matching number. The
labelled exponential enumerator is F(z,u/3,v). The lifetime enumerator is

    G(z,w)=(1+w z^3)/(1−z−w z^4−w z^6).

For lifetime ell>=1 its layer is

    z^(4ell−1)(1+z^2)^(ell−1)(1+z^3)/(1−z)^(ell+1).

It has eventual polynomial growth of degree ell, leading coefficient
2^ell/ell!, and minimum order4ell−1. The maximum lifetime at order n
within this family is floor((n+1)/4).

Every within-family birth inserts a singleton at a component boundary or
turns two adjacent singletons into a triangle with the newborn. Every
within-family arc flip is3↔111 or the labelled swap of adjacent11.
These rules include their exact multiplicities. A neutral singleton
insertion splitting one triangle run into two odd parts increases lifetime by one.
Gap-transport macros give explicit same-count paths between different
lifetime types, with every intermediate cost retained.

Thus a complete infinite type is now classified by size, composition,
symmetry, internal operations and all external extension divisibilities.
It remains only a subfamily of A000568. A unit extension gcd is not a
target selector: C3→C3 has C_1=1 yet its one-birth spectrum omits29.
Full proofs: [all-k theorem](arithmetic/TRIANGLE_ALL_K_PROOF.md),
[metagraph and enumeration](fibres/TRIANGLE_METAGRAPH.md),
[transport](fibres/GAP_TRANSPORT.md).

## 5. What107,73,17 and30000 actually preserve

The three decorated H9 sources have Berggren A-children107,73,17. Two
are even the same ordinary tournament with different apex charts. The
four homogeneous append coordinates(H,P,e,1) retain the needed terminal
bit; the inherited rank-four lower bound applies to linear homogeneous
representations, not every conceivable nonlinear encoding.

We completely classify all H9 one-apex cut words as I^a v O^b with
v in{OIII,OIOI,OOOI}. Their A-images have closed formulas. No two such
A-image values differ by30000 at any order. This is an actual all-height
statement about the user's source family, rather than a comparison of
three isolated integers. See [H9_FIBRES](fibres/H9_FIBRES.md).

Arithmetically, translation by30000 preserves precisely the prime residue
channels2,3,5, because30000=2^4·3·5^4. The requested examples are

    30017=13·2309,
    30073=17·29·61,
    30107=7·11·17·23.

For prime offsets among the first25 primes, translated primes occur exactly
at11,13,29,47,59,71,89,97. The translated candidate30067 passes trial
division by those25 primes but equals107·281. The full table and exact
Chinese-remainder wheel are in
[the arithmetic note](arithmetic/PROOF.md).

There is exactly one pair of A002496 primes at distance30000:

    547601=740^2+1,   577601=760^2+1.

Factoring (m−n)(m+n)=30000 reduces the unbounded question to fifteen
divisor pairs; every case and primality decision is certified. The
inherited aliquot bridge gives the corresponding A135244 identities
at547602 and577602. These are exact sequence intersections, not an
extension of the earlier global enumeration range or a proof of Landau.

There are also actual rooted tournament sources of counts518511 and
548511 with equal continuation coordinates(E,U,d)=(21,1040031,39).
Their pure H-gap is30000 and has a controlled propagation under common
continuation. Both starting counts are divisible by3. No tournament has
H=0, so a translation from formal zero is not itself a graph constructor.

## 6. The remaining sufficient construction

The useful next object is a family of actual marked sources, closed under
both desired target maps3H+2 and3H+4, with allowed replacements preserving
the count, root and designated port. The family must permit changing
cores and unbounded defect, or provide another way around the established
bounded-family restrictions. The uniform returns above, the finite
defect3 source supporting both target branches, and the exact obstruction lifetime provide
ingredients; they do not yet prove this closure.

The other sufficient route remains surjectivity of the five-core binary
address compiler onto every allowed odd value. The complete triangle
classification helps organize excluded and usable source types, while
the H9 branch collisions show why scalar count, a matching cardinality
of five, or the fixed-root observer alone cannot select the next witness.
All finite universes, scripts, actual maps, failed lifts and independent
audit scopes are preserved for the next attempt.
