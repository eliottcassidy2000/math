# An intrinsic tournament halving decoder and the obstruction to raw tripling

**Status: PROVED** elementary constructions and obstruction below;
**FINITE-EXACT** implementation controls; **OPEN** a structural generator of
all Collatz root certificates. This is one explicitly chosen tournament
encoding, not an assertion that number size selects a unique tournament.
No novelty claim or canon-ID promotion. Date: 2026-09-25.

## Inheritance, portfolio, and live board

Closest proved mechanism: modular substitution and the distinction between
strong components and modules in
[THM-1960, tournaments compose from regular seeds](../../01-canon/theorems/THM-1960-tournaments-compose-from-regular-seeds-the-spectral-substitution-law.md).
The first-doubling seam is inherited from
[THM-371, first doubling for units and pairs](../../01-canon/theorems/THM-371-first-doubling-unit-pair-seam.md).
The present construction proves its own claims; it does not need the
spectral formula or census in THM-1960. The canonical hostile is the
ten-vertex substitution without any pair module in the
[four-core note](collatz_tournament_core_20260925.md).
Corrected near miss: equality of vertex/edge counts did not supply a
halving quotient. Least-used sidecar: the actual guarded arithmetic word,
including its order, rather than the count of odd steps.

Anchor: recover halving from tournament edges. Niche: compatible odd
growth by two. Wildcard: transport this decoder interface to square sums,
Mahler carries, and sparse operation minors. The live board is:

| Concept | Exact object or invariant | Decisive probe |
|---|---|---|
| Odd chain | nested regular cores H_q | extend by two without changing old arcs |
| Doubles seam | graph of two-vertex modules | isolated vertices, edges, or longer paths? |
| Tripling | Q[H,H,H,1] | does it even have a pair module? |
| Root route | guarded inverse-Collatz word | expand each letter and verify its guard |
| Square-sum halving | square-class bit and attachments | divide even labels by two |
| Ordered carry | (length, odd count, carry) | words 01 and 10 have different endpoints |

## 1. The two arithmetic modes and the literal seam

Every positive integer has a unique address

    n = 2^k q = 2^k(2j+1),     k,j >= 0.

The odd chain is q -> q+2; doubling raises k. If S(n)=n+2 and D(n)=2n,
then D S = S^2 D. Thus scale is part of the operation: dividing a
literal +2 edge by two changes it into a +1 edge. On the valuation rows,

    v2(n)=0  => v2(n+2)=0,
    v2(n)=1  => v2(n+2)>=2,
    v2(n)>=2 => v2(n+2)=1.

The last two statements follow from n+2=2(q+1), q odd, and
n+2=2(2^(k-1)q+1), k>=2, respectively. This is an exact three-region
description. It gives an address for every integer; S is not thereby
an inverse-Collatz move.

## 2. A nested odd chain of regular tournaments

A tournament is regular when every vertex has the same out-degree.
Start with the one-vertex H_1. Given regular H_q for odd q, add vertices
a,b. Choose a set S of (q-1)/2 old vertices and orient

    a -> b,
    a -> S -> b,
    b -> (V(H_q) minus S) -> a.

Each old vertex gains one outgoing edge. Vertex a has 1+|S|=(q+1)/2
outgoing edges; b has q-|S|=(q+1)/2. The result H_(q+2) is regular and
contains H_q induced. This constructs every odd order by literal +2.
The implementation takes the first (q-1)/2 labels for S, making the
family deterministic. Different choices are allowed by the proof.

**Pair obstruction.** A regular tournament has no two-vertex module:
if u->v and every outside vertex sees u,v identically, then the
out-degrees of u and v differ by one. This contradicts regularity.
It does not say the tournament has no larger modules.

## 3. Doubling and halving can be exact and intrinsic

Write TT_b for the transitive tournament on b vertices, and define

    E(2^k q) = H_q[TT_(2^k)],                    q odd.       (1)

Here substitution replaces each vertex of H_q with a transitive fiber;
arcs between different fibers are uniform and follow H_q. Vertices are
labeled by (odd-core vertex, position in its fiber). There are no ties.
The orientation is a declared encoding gauge, not an arithmetic order
inferred from n alone.

Define P(E) on the same vertices, with an undirected edge exactly when
the two vertices form a module of E. **PROVED:**

    P(E(2^k q)) is the disjoint union of q paths P_(2^k).   (2)

Inside a transitive fiber, consecutive vertices form a module. A vertex
strictly between two nonconsecutive vertices distinguishes them. For
vertices in two different fibers, the absence of a two-vertex module
in H_q supplies a third core vertex distinguishing those fibers; every
vertex in its fiber also distinguishes the proposed pair. When q=1
there are no cross-fiber pairs. This proves (2), including k=0.

Consequently the user's three regions have a literal structural model:

| Valuation | Pair-module graph | Halving structure |
|---|---|---|
| k=0 | q isolated vertices | no pair to contract |
| k=1 | q disjoint edges | the first complete pairing |
| k>=2 | q disjoint even paths, length at least four vertices | overlapping eligible pairs, unique perfect matching |

An even path has a unique perfect matching: its endpoint forces the
first pair, and induction forces the rest. Thus the matching is
recoverable from the unmarked tournament, without being supplied as
metadata. It consists of consecutive disjoint pairs in each fiber.
Its quotient is

    E(n) / intrinsic matching  ~= E(n/2),         n even.  (3)

Indeed each TT_(2^k) becomes TT_(2^(k-1)); cross-fiber arcs stay uniform.
This proves a universal halving theorem for this encoding, with
canonical quotient **up to isomorphism**. Recovering the original
numerical vertex labels still requires their labeling map. The inverse
operation is uniform substitution by TT_2, since
TT_b[TT_2]=TT_(2b).

The number of eligible pairs is n-q. Their components recover q and
2^k structurally. Regularity of H_q is sufficient; more generally the
same proof works for any odd core with no two-vertex modules. If the
core has a pair module, extra cross-fiber pairs can occur, so the
assumption must not be discarded.

**Root boundary.** This family has E(4)=TT_4. A tournament of order four
partitioned into two uniform pair modules is necessarily TT_4: each
pair and their two-vertex quotient are transitive. Therefore a cyclic
four-core cannot also have this particular complete pair-halving rule.
The numerical Collatz loop at 4 is carried by the arithmetic transition
word below, not claimed to be an internal loop in a tournament.

## 4. The raw three-copy construction cannot satisfy this halving rule

Let H be any regular tournament of odd order q>=3, and let Q be any
four-vertex tournament, with its fourth vertex distinguished. Set

    F_Q(H) = Q[H,H,H,1].

**PROVED:** F_Q(H) has no two-vertex module, for every Q.

If a proposed pair lies in the same H block, regularity supplies an
internal vertex distinguishing it. If its vertices lie in different
blocks, take the one in an H block. That vertex has both an internal
win and an internal loss because q>=3 and H is regular. Its proposed
mate sees the whole block uniformly. One of those internal vertices
therefore distinguishes the pair. This includes a mate in the singleton
block. These cases exhaust all pairs.

Since 3q+1 is even, E(3q+1) has a complete matching of pair modules.
Thus no choice of Q makes the unmodified F_Q(H_q) isomorphic to E(3q+1)
for q>=3. This is an obstruction to simultaneous raw substitution and
this uniform pair-halving decoder, not to all possible tournament
encodings of Collatz. At q=1, F_Q(H_1)=Q and the obstruction need not
hold; a transitive Q gives E(4).

The positive next question is precise: can an explicitly defined
normalization of the three-copy object produce a halving structure,
while retaining enough information to decode the arithmetic route?
Discarding all arcs and rebuilding E(3q+1) from its size would not
answer that question. Arc reversals need their own reconstruction data.

The [complete order-ten repair probe](decoder_pair_repair_20260925.md)
answers a first local version: for H=C_3, the minimum reversals needed
to admit some complete pair-module matching are 3, 6, 7, or 8, depending
on the marked core Q. All 64 cores and 945 matchings per core are
checked. Every optimum retains one pair in each triangle; on this family
the cost has a proved three-choice formula on Q. Imposing the intended
quotient H_5 instead costs 7, 10, 12, 13, or 15 reversals, by an independent
complete check of all 1,451,520 matching/quotient candidates. In one core
class, every optimum now keeps only two internal triangle pairs; the old
optimizer restriction would miss it. Thus the target quotient changes
both the cost and the best contraction. The arithmetic interpretation
of the retained repair data remains OPEN.

## 5. A working arithmetic verifier, with the existence obligation exposed

Use the **ordinary** Collatz map: even x goes to x/2, odd x to 3x+1.
Its inverse steps are

    D(x)=2x,
    O(x)=(x-1)/3,     allowed exactly when x == 4 (mod 6).

A word over {D,O}, read from 4, is a root certificate for its endpoint
when every O guard holds. The guard ensures a positive odd predecessor;
reversing each step then gives a genuine forward route to 4. The script
implements this verifier, rejects unknown symbols and incorrect
endpoints, and independently replays the reversed word forward.

Example: DDODO gives 4 -> 8 -> 16 -> 5 -> 10 -> 3. The root-cycle word
ODD gives 4 -> 1 -> 2 -> 4. There is no unknown branch choice in checking
a supplied finite word.

If n=2^k q and a word W certifies q, then W D^k certifies n. Equivalently,
the forward route n -> ... -> q consists of k forced halvings. Thus
universal ordinary root certification reduces exactly to the odd row.
The script discovers the 1000 odd certificates q<=1999 by bounded
forward iteration, then certifies all 2000 targets n<=2000 by this lift.
An undecided start raises an error; it is never filtered out.

**OPEN:** generate a root word for every odd q from the tournament
structure, without presupposing that its Collatz orbit returns. The
present tournament quotient recovers the dyadic address; the supplied
arithmetic word recovers the route. They are separate pieces of data.
The first theorem does not supply the second existence theorem.

Under the shortcut map (odd step (3x+1)/2), 1 and 2 never visit 4.
The [ordered-carry note](decoder_mahler_catalan_20260925.md) therefore
uses certificates to 4 for n>=3 and handles the terminal cycle separately.
No certificate statement silently changes this clock convention.

## 6. Connections that preserve an actual predicate

The [square/primes note](decoder_prime_square_20260925.md) shows that
halving induced even vertices exchanges square sums with twice-square
sums. It preserves adjacency with one square-class bit, while losing
odd vertices and attachment obligations. Its repaired Q_24 obstruction
and marked-edge decoder at Q_25 show how retaining attachment data can
turn a contraction into an exact reversible construction.

The [Mahler/Catalan note](decoder_mahler_catalan_20260925.md) gives the
ordered-carry endpoint equation 3^a n+R=4*2^L. The clock (L,a) is
insufficient: 01 and 10 have carries 2 and 1. Its direct parity-word
transfer to the Mahler ceiling map sends root 4 to -4/15 and violates
the real safe-tail bound. A decoder must preserve its ordinary domain
as well as its symbolic address.

The [minor note](decoder_minors_20260925.md) finds an exact planarity
threshold at 16 for the explicitly selected sparse graph of literal +2
and doubling. Row-scaled +2 instead gives the planar quadrant grid.
It also separates strongly connected minor contractions from uniform
pair quotients. A compressed path can retain its guarded word; a bare
minor does not reconstruct an arithmetic route.

After these pulls the concept board has a common requirement: recover
the desired operation from the compressed object together with exactly
specified extra data. Here the pairing is intrinsically recoverable;
the missing global data are the odd root certificates. The square case
has a finite marked-edge decoder, while the same-word Mahler transfer
has a proved domain obstruction. These are distinct outcomes, not
assertions that the underlying conjectures are equivalent.

## 7. Reproduction and independent audit

    python 04-computation/experiments/decoder_halving_20260925.py

[Recorded output](decoder_halving_20260925.out). Standard library only.
Explicit universes: 64 nested odd cores through 127; every pair of every
E(n) for 1<=n<=128; halving all even n<=128; all 64 labeled Q for every
odd q=3..15 (448 hostile substitutions); every arithmetic seam through
4096; every positive certificate target through 2000, with a 10000-step
discovery limit and maximum observed certificate length 179.

Independent coordinating audit checked the general proofs, used the
other representative in each contracted pair for n/2<=32, and checked
all 24 labeled regular order-five cores against all 64 Q (1536 further
substitutions). The pair-path proof was independently checked in its
within-fiber and cross-fiber cases. Computations are finite controls;
the universal claims have the elementary proofs above.
