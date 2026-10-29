# Run refinement and closure under actual native complements

2026-10-28. **PROVED, scoped / FINITE-EXACT controls; same-agent audit.**
Local statements C1-C5 belong to this file. No universal smallest carrier for
every possible strong-atom operation is asserted.

## C1. The smallest supplied-run carrier once cuts are legal

Fix an aligned tournament pair (A,B) on n vertices whose agreement arcs form
a directed spanning path forest F, with k components. An input orders those
components as R_1,...,R_k. Observe ordinary unlabelled complement partners of
legal quotient-walk assemblies and the ordered run lengths. Allow a run to
be cut at a specified position, repeatedly, before such an assembly.

**Theorem.** Two such inputs are indistinguishable by this operation menu iff
one common element of G=Aut(A,B) carries their full ordered run words to one
another. Over the fixed pair the minimum is k!/|G| classes.

**Proof.** A common map transports every cut, supplied-run assembly and actual
complement, proving sufficiency. For necessity, cut every original run into
singletons, retaining the declared cut word and original lengths. The child
partner is now A^op. The intrinsic strong C3 context of
[ENDPOINTS E2](../cursor_cm_phase_20261028/ENDPOINTS.md) recovers this entire
partner and the ordered singleton endpoints from its ordinary unlabelled
output. Thus it recovers A and the complete vertex word. The retained lengths
recover the original runs and F. An isomorphism therefore lies in Aut(A,F),
which equals G. The n=1 case is immediate; for n>=2 use the E2 context with
n singleton runs. Finally G acts freely on complete component orderings:
fixing every directed path setwise fixes its vertices pointwise. QED.

This uses the inherited E2 decoder, rather than claiming a new arbitrary-port
reconstruction theorem. The menu explicitly allows named run positions;
unnamed cuts without a recoverable choice require their own observer analysis.

**Concurrent intermediate quotient.** Incoming ENDPOINTS E6 permits a
forest-native complement before assembly, without cutting the runs. It proves
the exact two-sheet group J=Hbar_A intersect Hbar_B, on the same component
names. Thus G<=J<=Hbar_B, and the corresponding state counts lie between the
one-stage endpoint quotient and C1's cut-closed minimum. Maps on the two sheets
need not agree on interiors, which is why their abstract intersection is not
automatically G. The incoming examples all give J=G, but this is not proved
universally. Its strict-strong-child caveat matters here: the eight-vertex B
has a source and sink, so switching it into the source slot would leave that
restricted menu. C2 below instead uses an actual all-strong assembly source.

## C2. One intrinsic cut restores a trit; sixteen vertices realize the loss

The inherited first length-aware loss is at eight vertices (ENDPOINTS E4):

    A = (2,104,99,37,15,145,57,95)
    B = (254,156,152,240,224,198,134,0)
    U = ((0,1,3,5,7),(2),(4),(6))
    V = ((0,1,3,5,7),(4),(6),(2)).

Rows encode out-neighbours as bit masks. A is strong. The initial B is not
required to be strong. Both ordered endpoint/length records are isomorphic,
although no simultaneous source/partner map transports U to V.

Cut the unique nonsingleton run after its first vertex. This choice is
intrinsic to the ordered-run interface. The resulting two covers are

    U' = ((0),(1,3,5,7),(2),(4),(6))
    V' = ((0),(1,3,5,7),(4),(6),(2)).

Their new common partner has trivial automorphism group, so the rooted
endpoint records differ. Across all 24 original run orderings, the 8 old
endpoint cells become 24 cells after this single cut: one whole trit is
necessary. Each of the four internal cut positions separates the displayed
pair, as checked in `ports.json`.

There is also a small **actual native-complement context**. Let E_4 be the
increasing transitive tournament on four vertices with only its first-last
arc reversed. Take S=C3[A,E_4,E_4], with block offsets 0,8,12. All three
children, and S, are strong. Use the two Hamilton paths

    P_U = (0,8,12,1,3,5,7,9,13,2,10,14,4,11,15,6)
    P_V = (0,8,12,1,3,5,7,9,13,4,10,14,6,11,15,2).

Reverse precisely the arcs outside the chosen path. Their partners have

    H(C_PU(S)) = 629677591,
    H(C_PV(S)) = 634905785.

These unequal Hamilton counts certify different ordinary isomorphism classes.
The exact same labelled S is used on both sides. Sixteen is a verified size,
not a claimed minimum. An eighteen-vertex variant is retained in the receipt.

**Failure mechanism.** A target-only endpoint automorphism rotates singleton
runs while moving hidden internal vertices of the source forest. A cut makes
one of those vertices observable. The old one-stage theorem remains correct;
its expressly excluded run-refinement extension is what fails here.

## C3. Actual complement covariance and a closed sufficient carrier

For a Hamilton presentation (T,P), write C_P(T) for native complement. Then

    C_P(C_P(T))=T,
    g C_P(T)=C_(gP)(gT).

Both identities follow arc by arc. Moreover the agreement arcs of T and
C_P(T) are exactly the directed Hamilton path P. Hence the aligned full pair
recovers P, and Aut(T,C_P(T)) is trivial. The ordinary endpoint classes
separately can have large automorphism groups; this does not contradict that
joint rigidity.

A sufficient all-size carrier for general strong-atom composition is the
full aligned relation together with ordered run words, shared vertex identities,
and actual compatible maps. Assemble seams along the quotient walk; reverse
nonpath arcs; then recompute the strong-module decomposition of the resulting
tournament. The new tree is a derived representation, not a presumed invariant
of complement. Re-marking transports the common relation and words together.
This closes the operations by construction. C1 proves minimality only for
the narrower fixed-pair cut/assembly interface, not for this whole menu.

The complement checks in `ports.py` use twelve literal Hamilton presentations
of the order-nine iterated triangle and all 81 source automorphisms. They test
involution, covariance and joint rigidity. This is a finite check of elementary
all-size identities, not an exhaustive native metagraph census.

## C4. Weights must be pushed through the common action

For a finite group action, the transporter kernel is

    K_G(u,v) = #{g in G : g u=v}/|G|.

It is the matrix of the Reynolds averaging operator on the declared tuple
space. An equivariant map f satisfies f_* Reynolds_G = Reynolds_G f_* when
the output action is the induced common action. C3 therefore gives exact
weighted native-complement transport on full presentations. Tensor factors
may be multiplied only for independent child group variables; two ports
using the same child must enter the same joint kernel.

The actual 81/27 hostile is rerun in `ports.json`: equal labelled cyclic
triples coexist with automorphism orders 81 and 27, extension-orbit counts
176 and 512, and transporter 0->4 sizes 27 and 0. It rejects replacing actual
maps by cyclic triples, stabilizer orders, or unrelated marginal kernels.

Next decisive question: for a fixed number of cuts and attachments, how much
of C1's full ordered word can be omitted? The unbounded cut theorem gives the
limiting answer; it does not settle the bounded-depth rate.

## C5. Global converse is not the actual native-complement operation

Incoming [CANONICAL K2](../weighted_crt_clocks_20261028/CANONICAL.md) proves a
one-bit selector for the transformation f(v)->1-f(alpha(v)). This is global
arc reversal, followed by identifying the reversed parent with its original
copy. Its earlier label "native complement" was incorrect and has been
repaired in the theorem, summary, audit and index. Its proof and counts remain
valid for global converse.

Here is a literal separating source. Give parent vertices0,1,2 their C3 arcs
0->1->2->0 and attach vertex3 with 3->0, 1->3, 2->3. The resulting tournament

    S=(2,12,9,1),     P=(3,0,1,2)

is strong and P is Hamilton. Actual C_P(S) retains 0->1 and1->2 but reverses
2->0 into0->2. Its induced parent is transitive. Global converse instead
retains a cyclic parent. Thus the fixed-parent selector cannot implement this
native move: the parent even leaves the recursive family. Four vertices are
a verified witness, not a general minimum claim. `audit.py` checks it.

The strongest survivor is the incoming complete joint-colour canonicalizer
plus its global-converse selector. Closing actual C_P still requires the path
and target decomposition, as in C3. The operation-type repair is logged in
`01-canon/MISTAKES.md` under this packet's slug, without reserving a numeric ID.
