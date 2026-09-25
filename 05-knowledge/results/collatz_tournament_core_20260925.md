# A four-vertex core for tournament tripling: exact counts and certificate obligations

**Status:** PROVED elementary counting, substitution and obstructions;
FINITE-EXACT controls. A universal Collatz certificate remains OPEN.
The edge-orientation rule in the user's proposed construction is not yet
specified; Q below is a parameter, not an inferred user choice.
Session: codex-collatz-bugs-20260925. No novelty or theorem-ID promotion.

## Inheritance and board

The user proposes encoding every integer by tournament structure, with
three copies of an odd-order tournament H and one extra vertex giving
order3A+1. Selecting one vertex from each copy exposes a four-vertex core.
The intended consequence is a provable certificate route to root bug4.

Closest proved mechanism: substitution into a fixed quotient tournament,
as used in [THM-4466, common-edge and substitution classes](../../01-canon/theorems/THM-4466-sharp-tournament-cubic-bound-for-common-edge-and-substitution-classes.md),
section3; here only the elementary uniform cross-block construction is used.
The closest arithmetic four-vertex construction is
[THM-4472, four-vertex reading of 3n plus/minus1](../../01-canon/theorems/THM-4472-four-vertex-reading-of-3n-plus-minus-1.md).
Its vertices are an AM-fair pair and their images, a different construction
from three expanded blocks. Do not identify their cores without a map.
Its corrected near miss: the converse diamonds encode time reversal, not
the positive 3n+1/3n-1 domain swap. See its theorem and the 2026-09-24
corrections in MISTAKES.
Canonical hostiles: a root that is a sink despite a Hamiltonian path;
an even-order tournament without a contractible pair. Least-used sidecar:
the marked block partition, plus a decoder to guarded arithmetic moves.

Anchor: root certificates. Niche: exact four-core substitution. Wildcard:
whether even-order substitutions support intrinsic halving.

| Concept | Preserved information | Missing coordinate / decisive test |
|---|---|---|
| Root bug4 | inverse-Collatz reachability | arithmetic word; no assumed universal membership |
| Tournament order | integer size | orientation rule; arbitrary tournaments share order |
| Four-core substitution | quotient Q | uniform cross-block arcs and marked partition |
| Hamiltonian path | one directed visit per vertex | prescribed root and arithmetic labels |
| Halving | cardinality divided by2 | allowable contraction; pair-module hostile |

Method cards used: "Search the statement before the method" and "Type every
analogy and every implication" in META-PATTERNS. No new card is promoted.

## 1. Edge partition (PROVED)

Use tau(k)=k(k+1)/2 for triangular numbers, to distinguish them from the
Collatz map and the reversal distance f. An order-m tournament has
binom(m,2)=tau(m-1) edges, regardless of their orientations.

Take disjoint blocks H1,H2,H3, each with A vertices, and a new vertex r0.
Choose ri in Hi and write Bi=Hi minus {ri}. The marked core
K={r0,r1,r2,r3} has4 vertices and6 edges. Partition all edges into:

1. Within the three original Hi: 3 binom(A,2).
2. Between distinct Bi: 3(A-1)^2.
3. Core-to-periphery edges not already counted in item1: 9(A-1).
4. Within K: 6.

For item3, each of the 3(A-1) peripheral vertices has four potential core
neighbors; its edge to its own representative was already counted in item1,
leaving three. Thus exactly

    tau(3A)=3 tau(A-1)+3(A-1)^2+9(A-1)+tau(3).

The proposed 3*(3A)=9A term exceeds item3 by9. The minimal witness is
A=1: six actual edges versus fifteen in the proposed count. At A=3 the
correct count is9+12+18+6=45, not54. This is an endpoint-counting repair,
not a failure of the four-core idea.

Equivalently, count only residual interiors:

    binom(3A+1,2)=3 binom(A-1,2)+3(A-1)^2+12(A-1)+6.

Both identities hold for every positive A. Oddness is needed only for
the intended Collatz odd step; counting cannot detect its parity guard.

## 2. When the core is stable (PROVED construction)

Fix an oriented four-vertex tournament Q. Define

    F_Q(H)=Q[H,H,H,{r}],        |F_Q(H)|=3|H|+1.

Replace three Q vertices by copies of H and the fourth by a singleton.
For every quotient arrow i->j, orient ALL edges from block i to block j.
This is tournament substitution. The three copies are actual induced copies,
not merely equal-size sets. Every choice of one representative from each
copy, together with r, induces the same Q. Collapsing the four blocks
gives Q, with the block partition retained as part of the certificate.

This is the precise scope in which the core is independent of A and H.
An arbitrary order3A+1 tournament need not have three isomorphic blocks or
uniform arcs between them. Merely selecting four vertices gives an induced
tournament, not an operation-preserving quotient. The edge count does not
select any of the four isomorphism classes of tournaments of order4.

For a fixed labeling there are2^binom(N,2) tournaments of orderN, since
each unordered pair contributes one orientation bit. On a fixed marked
partition, the family F_Q(H) with three identically labeled copies has
only2^(binom(A,2)+6) choices as H and Q vary, or2^binom(A,2) for fixed Q.
Thus the stable-core construction repeats the input information; its other
orientations are constrained. Extra cross-block orientation bits could carry
additional certificates, but then the uniform quotient requires repair or
those bits must be retained explicitly in its decoder.

Connection contract: input (Q,H,marked copies) -> F_Q(H), inverse quotient
preserves Q and its interblock arrows, forgets all internal orientations;
retain the child tournaments and block addresses to reconstruct the input.
An arithmetic certificate additionally needs numerical labels and a decoder.

## 3. Hamiltonian paths supply no rooted arithmetic certificate (PROVED)

Every tournament has a directed Hamiltonian path. For the existence part,
insert a new vertex v immediately before the first vertex p_i of an existing
path with v->p_i, or append v if none exists. The preceding vertex, if any,
points to v by minimality of i. This proves existence inductively without
assuming Collatz. The stronger parity theorem is
[THM-001, Redei](../../01-canon/theorems/THM-001-redei.md).

In a substitution, a Hamiltonian path of Q can be expanded by inserting a
Hamiltonian path of each block consecutively; uniform cross-block arcs make
every splice valid. This is a positive structural mechanism.

But a designated root may be a sink. Such a tournament still has a
Hamiltonian path, ending at that root, and no directed path leaves the root.
Even a Hamiltonian path beginning at the root would certify tournament
reachability only. Its arcs must decode to valid Collatz transitions before
they can establish f(4,X)=0. Count, orientation, distinguished root and
arithmetic decoding are separate pieces of information.

## 4. A hostile inside the proposed family: halving is extra structure

An unambiguous contraction retaining all external arc directions can
identify a pair {u,v} only when every outside vertex w has the same relation
to u and v. Such a pair is a two-vertex module. Order even does not imply
the existence of any such pair.

**PROVED explicit hostile.** Let H=C3 and Q be a source over C3. Then
F_Q(H) has10=3*3+1 vertices, with three cyclic blocks and a source r.
It has no two-vertex module:

* Two vertices within one C3 are distinguished by its third vertex.
* Vertices in different cyclic blocks are distinguished by a vertex in
  the third block (one of the two blocks beats it, and the other loses).
* The source and a vertex in block i are distinguished by any vertex in
  the block preceding i: it beats the latter and loses to the source.

Thus this valid four-core tripling construction does not support halving
by uniform pair contraction. This refutes that specific automatic-halving
route, not all possible tournament encodings or all halving operations.
Allowing deletion, nonuniform contraction or extra recorded data changes
the operation and requires its own preservation proof.

## 5. The root claim and the missing theorem

For bug identities, the earlier statement f(4,X)=0 means that a path
following inverse-Collatz arrows runs from root4 to X. Universal validity
is equivalent to Collatz, as proved in
[the bug-distance note](collatz_bugs_20260925.md), section2. It is the target
to certify, not an assumption licensed by calling4 the root.

Here the size move A->3A+1 uses the ordinary Collatz odd step C, whereas
the previous note uses shortcut T(A)=(3A+1)/2. Track that clock distinction:
one must also define the even halving operation on tournament encodings.
Under ordinary C the root cycle is4->2->1->4; under shortcut T it is1<->2.
The identity I=3A+1 is exactly the bug whose first odd leaf is A.

A checkable inverse certificate for ordinary C is a word starting at4
in the operations

    D(x)=2x,
    O(x)=(x-1)/3, allowed exactly when x=4 mod6.

Each operation is locally verifiable. The ordinary root cycle is reachable
from4 by these operations, so certificates can include1 and2 without a
separate exception. In the shortcut root-tree convention, they are removed
or contracted as specified in the bug-distance note.

The quotient F_Q(H)->Q has cardinality jump3A+1->4; it is not itself a
Collatz transition. To make it a numerical certificate, give a decoder
that turns it into a finite valid D/O word to the represented integer.
Constructing the tournament by first recording an already-known Collatz
trajectory verifies that trajectory but cannot establish universal coverage.
The structural obligations are therefore an explicit orientation rule,
compatible halving, arithmetic soundness of the decoder, and a proof that
every integer receives a finite certificate. None follows from edge counts.

Board revisit: corrected counts support the core niche; module failure
separates odd construction from halving; the sink separates Hamiltonicity
from rooted reachability; decoding returns the anchor to guarded words.
The original orientation rule remains unspecified, and no tournament
relation is inferred solely from matching cardinalities.

## 6. Exact controls

Run from the repository root:

    python3 04-computation/experiments/collatz_tournament_core_20260925.py

See [output](collatz_tournament_core_20260925.out). Unfiltered universes:
all A=1..25 for counts; all27 representative choices and all45 vertex pairs
in the explicit ten-vertex substitution. Pair enumeration independently
checks the algebraic count; module enumeration checks the proof's actual
contraction consequence. Positive controls retain the quotient and a
Hamiltonian path; hostile controls detect the excess9, failure of pair
halving, and the sink-root failure despite a Hamiltonian path.
