# A closed singleton carrier and its sharp reconstruction boundary

2026-10-28. Identifiers S1--S5 and P1--P3 are local to this file. **PROVED**
means the elementary arguments below, not a global theorem reservation or an
independent-agent review. Computational statements are separately **FINITE-EXACT**.

## S1. Declare the operation language

Use exactly [N4 of the inherited compiler](../strong_attachment_tail_compiler_20261028/THEORY.md):
X_n is the set of tournaments on ordered positions 0,...,n-1 whose consecutive
arcs point forward. Its free bits are the nonconsecutive arcs. The observations
are ordinary isomorphism class q(T) and the blue flag b(T), meaning invariance
under converse followed by reversal of positions (the flag is suppressed at
n<3). Operations are native complement tau, which flips all free bits;
ordered adjoining of a universal source J_L or sink J_R; and **fresh uniform
Hamilton-path re-marking**. Counted fresh re-marking is also covered.

This language has no named-path legality query, selection of a re-marking by
hidden coordinates, interior attachment, or arbitrary joint-port gadget.
Those richer languages can distinguish more states. This qualification is
essential: the old six-vertex quotient-walk counterexample is not refuted.
The target here is behavior of a Hamilton presentation, not reconstruction of
an unlabelled tournament whose class q is already an observation.

Write E(T) for the isomorphism type with its **ordered first and last vertices
marked**. This includes the relation, not just its endpoint degrees or group
order. At n=1 the marks coincide.

## S2. A linear boundary record closes the blue observation — PROVED

For each i, record whether the free forward row (T(i,j): j>=i+2) is all0,
all1, mixed, or empty. Similarly classify each free column
(T(i,j): i<=j-2). Let P_k and S_k say whether the first/last k positions are
fixed by converse/reversal; set them true at k=0,1,2.

The full temporary record consists of row/column trits and these prefix/suffix
bits. Its update formulas are exact. Under tau, exchange0 with1 in every
nonempty constant row/column and leave P,S unchanged. Under J_L:

* prepend an all1 free row (empty if n=1); existing rows are unchanged;
* append a1 to every old column j>=1, prepend an empty column, and keep the
  old first column empty. Appending1 sends all1/empty to all1 and0/mixed to mixed;
* for k>=3, P'_k = [column k-2 is all1 or empty] P_(k-2);
* S'_k=S_k for k<=n, and S'_(n+1)=P'_(n+1).

The J_R formulas are reflected: append an all1 free column, append1 to rows
i<=n-2, and

    S'_k = [row n-k+1 is all1 or empty] S_(k-2),
    P'_k=P_k (k<=n),  P'_(n+1)=S'_(n+1).

For the prefix formula, remove the new universal source and its reflected
last vertex. Blue requires that last vertex to be a sink in this prefix;
after removing the pair, precisely the smaller prefix condition remains.
Consecutive arcs and the arc between the removed endpoints are automatic.
This proves the formulas, including their sufficiency. Reflection proves the
suffix formula. No enumeration hypothesis enters this argument.

We can discard most of that temporary record. Let a be the number of initial
constant free rows before the first mixed row, stopping at n-2, and c the
number of final constant free columns before the first mixed column. Retain

    Psi(T) = (n; u_1,...,u_a; v_1,...,v_c;
                  S_(n-1),...,S_(n-a); P_(n-1),...,P_(n-c); b).

Here u runs from the first row inward and v from the last column inward;
b=P_n, with the observation suppressed separately below order3. This uses at
most4(n-2)+1 Boolean entries plus two lengths for n>=3.

**The compressed update is well-defined at every size.** Under J_L the new
row run is (1,u), and its suffix bits are (b,S_(n-1),...,S_(n-a)). The new
column run is the initial run of1s in v, extended by the newly available
empty old column only if every old free column was1. Every needed new prefix
bit is the displayed formula above: it needs the next retained column/bit,
or is false when the next column was mixed or0, or is the trivial true value
at lengths<=2. Thus no discarded entry is consulted. The new whole blue bit
uses the last column and P_(n-1) by the same rule. J_R is symmetric; tau
complements u,v and leaves the flag lists fixed.

`boundary_profile.py` implements these rules by filling unused entries with
formal mixed/false placeholders, applying the full formula, then truncating.
These placeholders need not be realizable as one tournament; the dependency
argument proves they cannot affect the retained answer.

FINITE-EXACT: all three updates are checked against literal relations for
every presentation through7 (101604 updates). The numbers of distinct Psi
records at n=1,...,7 are1,1,2,7,26,68,154.

## S3. Explicit all-size closed carrier — PROVED

The following is sufficient and closed under the declared operations:

    K(T) = (E(T), E(tau T), Psi(T)).

Tau swaps the two rooted types and applies S2. For J_L, its first component
adjoins a universal source to E(T), marking the new source and old last port.
Its second component adjoins to E(tau T) a vertex that beats **only the old
first port**, marking the same new/old ends. For J_R, adjoin a sink to E(T);
on E(tau T) the new vertex loses **only to the old last port**. These are
well-defined constructions on endpoint-rooted isomorphism classes. S2 supplies
the blue coordinate. The observation q is obtained by forgetting E's ports.

Fresh re-marking has a class-constant row: Aut(T) acts freely on Hamilton
paths, so every presentation of the ordinary class appears with multiplicity
|Aut(T)| and probability |Aut(T)|/H(T). Hence the law of the new K depends
only on q(T), already retained. The same argument covers counted reset.
This proves closure, including any finite adaptive combination of these
operations. It does not reconstruct named transporters or arbitrary attachments.

K is a **sufficient explicit carrier**, not a proved smallest one. Its two
rooted relation types can still be large. The linear bound concerns the extra
blue interface, not the total information in the relations. At order7 its
realized quotient has32272 cells among32768 presentations; this is an upper
bound for the number of minimal behavioral cells. Horizon7 gives the lower
bound24754. We have not shown equality with either bound.

## S4. Reconstruction fails first at seven — PROVED / FINITE-EXACT boundary

Use lexicographic free-edge order (i,j), j>=i+2. The two masks1072 and4176
give outgoing-neighbor bit rows

    A = (66,12,41,17,39,75,30),
    B = (66,20, 9,51,37,71,30).

Both are strong, distinct presentations. Explicit source-to-target maps are

    A -> B:         (0,1,4,2,3,5,6),
    tau A -> tau B: (0,2,4,1,5,3,6).

Both preserve the two endpoint labels. Their first free rows and last free
columns are mixed, and neither presentation is blue. Consequently each has
Psi=(7; empty; empty; empty; empty; false), and S3 proves indistinguishability
for **every future size**, including fresh reset. Direct arc checks verify
the displayed certificates. H(A)=H(B)=127; H(tau A)=H(tau B)=93; all four
automorphism groups are trivial.

There is also a useful direct explanation. In any deterministic context the
original block stays contiguous; its induced relation is T or tau T. Every
new vertex sees its interior uniformly, with an exception possible only at
an original end port. Extend the appropriate endpoint isomorphism by the
identity on all new vertices to preserve ordinary classes.

To see why blue never returns, let L,R be the numbers added on either side
and d=R-L. Reflection on old local coordinates is i -> n-1+d-i. At d=0 blue
would require the original block to be blue. At d>0, old vertex0 reflects
to a right exterior vertex. All reflected pairs of (0,j), j>=2, join that
vertex to older vertices by nonconsecutive arcs. Those arcs were equal when
the exterior vertex was inserted and every later tau flips them together.
Thus blue would force the original first free row constant. At d<0 the last
free column gives the symmetric contradiction. At the first reset the two
ordinary classes agree, so their entire subsequent state laws agree.

The computation rechecks all21845 deterministic contexts with at most seven
insertions, including all interleaved complements and a final complement.
This finite check corroborates the proof; it is not the source of the all-size
quantifier. There are412 paired-root/blue cells with at least two members
satisfying the mixed-boundary criterion, so the obstruction is not isolated.

For every n>=7 there are distinct **strong** presentations with identical
behavior. Apply L^- = tau J_L tau repeatedly to A,B. It retains the old
relation and adds a vertex beating only the old first port. Strongness is
preserved: the new vertex reaches the old strong block through that port,
and any other old vertex reaches the new one. The restriction to the old
block preserves distinctness; behavioral equivalence is preserved by context
composition. Orders7..30 are also checked literally.

Sharpness uses the inherited exact N4 census, replayed here: horizon7 separates
every presentation through5, leaving only five pairs at6; all five separate
by horizon8. Thus **full Hamilton-presentation reconstruction under this exact
menu is REFUTED, with least order7**. Named legality reconstruction and larger
joint-port menus are unaffected.

## S5. Hostile symmetry and phase remain real data

The12-vertex81/27 pair is rechecked independently:67 identical labelled cyclic
triples, H=462915 on both, but automorphism orders81/27, binary attachment-orbit
counts176/512, and transporter0->4 sizes27/0 in this packet's numbering.
These controls concern actual maps, not only orbit totals. They warn against
replacing either E component with cyclic triples, degrees, or a group order.

The inherited [N1/N3 frame and Frobenius results](../strong_attachment_tail_compiler_20261028/THEORY.md)
retain their precise scope: coherent free-frame transports compose along LCA
paths; distinct F64 normal ports admit phases0,2,4, and F256 only0 as tournament
automorphisms. These constraints do not reconstruct missing path alignments.
P2 below gives a fresh literal three-phase control.

## P1. A square of bits is a labelled consecutive pair — PROVED

For an n by n binary matrix M, use its strict upper triangle for a tournament
A on[n], its strict lower triangle for B restricted to[n], and its diagonal
for the arcs *->i of B on[n] union{*}. This is a bijection because

    binom(n,2) + binom(n+1,2) = n^2.

The correct quotient uses **simultaneous relabeling of A and B fixing***.
It is not ordinary matrix conjugation: at n=2, the16 matrices have10
conjugacy orbits but8 synchronized pair orbits. At n=3 these are104 and88.
Separate isomorphism classes also discard the alignment. For fixed classes
A and (B,*), alignments are the double cosets

    Aut(A) \ S_n / Aut(B,*).

The stabilizer of a particular alignment pi is the intersection
Aut(A) intersect pi Aut(B,*) pi^-1. This is the same joint-witness requirement
that protects the81/27 and Paley port controls.

## P2. Recover the classical odd self-converse representation

**CITED classical:** Eplett's [1979 structural theorem, pp.23--24](https://doi.org/10.4153/CMB-1979-004-6)
represents an odd self-converse tournament by a tournament with a loop-allowed
undirected graph on half the remaining vertices. The anti-involution has one
fixed pole. Its choice does not change the representation's isomorphism class:
the anti-involutions are conjugate under the tournament's odd-order
automorphism group. This gives a bijection on unmarked isomorphism classes.

Combine that theorem with the repository's already-audited
[XOR disagreement map](../converse_tournament_count_channels_20261021/session22-converse-review/COROLLARY_REVIEW.md).
Here is the explicit consecutive-pair specialization, using vertices i+,i-,*:

    * -> i+,  i- -> *;
    i+ -> j+ iff A(i,j)=1;
    i- -> j- iff A(i,j)=0;
    i+ -> j- iff A(i,j) XOR B(i,j)=1  (i != j);
    i+ -> i- iff B(*,i)=1.

The involution exchanges + and - and fixes*. The XOR is symmetric under
exchanging i,j, so these rules are consistent. They have the literal inverse
implemented in `pair_square.py`, and are equivariant under simultaneous
relabeling. Eplett then identifies their unmarked output classes exactly.
This is a recovered classical bridge, **not a new structural theorem**.

The bijection of classes need not preserve full automorphism groups. The
paired representation retains the centralizer of the chosen anti-involution.
For Paley7, |Aut|=21 but that centralizer has order3. With pole0 on Z/7,
rho(i)=-i and the surviving map i->2i is its three-phase action. Under the
inherited F8 exponent model it is precisely Frobenius. Choosing the pole
removes a sevenfold owner choice, not the remaining threefold phase.

FINITE-EXACT: all bit squares through n=3, every simultaneous permutation,
both literal inverses, and all self-converse ordinary classes through7.
This checks a structural correspondence, not merely equal numerical lists.

## P3. Count aligned pairs, compatible extensions, and actual root blocks

For an odd cycle partition lambda of n, let c(lambda) be its cycle count,
e(lambda) the number of unordered-pair orbits, and z(lambda) its centralizer
size. A permutation fixes2^e choices of A and2^(e+c) choices of rooted B.
Even cycles fix no tournament. Hence

    D_n = sum_(odd lambda|-n) 2^(2e+c)/z(lambda) = SC(2n+1).

This is exactly the inherited channel F_n(2,1), now given the aligned-pair
interpretation and the explicit classical structural map. The checker reuses
that engine through100. The first counts are1,2,8,88,2752,279968,... .

If B restricted to[n] **equals A**, the situation is different: only the n
root attachment bits remain free, and the count is the pointed tournament
count p_(n+1)=F_n(1,1). That is the domain of root-switch enumeration. Once a
parent representative is fixed, enumerate its binary masks modulo its actual
automorphisms; a designated node by itself cannot encode n arbitrary arcs.

`lazy_root_switch.py` combines incoming exact individualization/refinement
with the inherited complete-root-block theorem and transporter replay. Store
one aligned root row per ordinary child, and for each pointed state retain
its owning block and an anchor port. The anchor's actual isomorphism transports
every record back to that pointed child. Thus full root records, groups and
Qv are recovered on demand, without storing a copy for every pointed state.
This is a concrete implementation of those inherited proofs, not a new
canonization theorem or an asymptotic speed claim.

At9 it generates191536 classes and1716608 pointed states, storing1723824
anchor root records instead of15449472 expanded records. All1763536 retained
maps are checked arc by arc; all lazy expansions through8 and1026 at9 are
recomputed independently. The complete9-bank agrees with the earlier
independent generator. The factored operator satisfies Q^2=9Q and
tr(Q)=1723824. Black/blue adjacency is not generated at this order.

Scalar A000568 and nonempty strong coefficients are separately extended
through350, with exact full-prefix agreement between the incoming degree
buckets and the earlier terminal-contraction implementation. a(350) has17646
decimal digits. These share the same proved Burnside mechanism; agreement is
an alternative implementation check, not an independent proof of Burnside.
They do not reconstruct a metagraph on350 vertices.
