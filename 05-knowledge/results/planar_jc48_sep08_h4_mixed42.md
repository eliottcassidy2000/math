# H4 generators of type (4)(2) are equal

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This is a uniform abstract permutation theorem with a complete finite
pair certificate.  It does not identify an abstract H4 representation
with an actual Keller map.  The final geometric consumer requires its
own proved marked supplier and actual generation.  JC(2) remains OPEN.

## 1. Statement and inherited mechanisms

Let a,b,c,d be permutations of a finite set Omega.  Each has exactly one
4-cycle, exactly one disjoint transposition, and arbitrary other fixed
points.  Suppose

    aba=bab,       bcb=cbc,       cdcdc=dcdcd,
    ac=ca,         ad=da,         bd=db.               (1)

**Theorem.** Then a=b=c=d.  In particular the four generators cannot
generate a transitive action: even their common permutation has the two
distinct nontrivial orbits of sizes 4 and 2.

The ambient degree is arbitrary.  The finite computation below concerns
one pair, whose union of moved supports has size at most twelve.  No
upper bound of twelve is imposed on Omega or on a four-generator action.
The global step is a joint-centralizer argument followed by the proved
[H4 single-cycle equality theorem](planar_jc48_sep08_h4_single_cycle.md).

The closest inherited mechanism is the newly proved
[mixed (3)(3) theorem](planar_jc48_sep08_h4_mixed33.md): a common
centralizer of an ordinary pair cannot move certain joint orbits.  The
[mixed (3)(2) theorem](planar_jc48_sep08_h4_mixed32.md)
also warns that isolating one cycle from an odd-braid pair need not
preserve the odd relation.  Here two minimal six-label hostiles pay that
warning for the present type.  The corrected near miss is to strip a
transposition before proving it is the same transposition for all four
generators.  The least-used sidecar is the restriction of a centralizing
order-four permutation to a six-point transitive orbit: its order-two
possibility would require **three** transpositions, which the declared
cycle type does not possess.

The five live concepts are unique-length cycle blocks, literal odd
relations, joint orbits, semiregular centralizers, and a common factor
that can legitimately be removed.  The map from a pair to its block
intersection matrix retains typed block sizes but loses cyclic order;
the full literal pair enumeration keeps the latter until the exact
relations have been evaluated.  The map from four mixed generators to
four single cycles is used only after the common two-point orbit is
proved.  No quotient is silently assumed to preserve the braid words.

Throughout, X4 and X2 denote the moved 4-block and 2-block of a
permutation x, and X=X4 union X2 its entire moved support.

## 2. The complete pair universe

Fix sigma=(1 2 3 4)(5 6) in S_12.  Every ordered pair of the declared
type on an arbitrary finite set is simultaneously conjugate, on its
union of moved supports, to sigma and one of the partners below,
after padding with fixed points to twelve labels.  Conversely every
partner is tested without inherited filters.

Choose a four-set, one of its six cyclic orders, and a two-set from
its complement.  The exact number of partners is

    binom(12,4)*6*binom(8,2)=83,160
        =binom(12,6)*90.

The source visits each once.  Since the block lengths differ, no
arbitrary block-order gauge is chosen: the rows and columns of the
intersection matrix are always ordered by lengths 4,2.  Write its four
entries in the order

    (|S4 intersect T4|, |S4 intersect T2|,
     |S2 intersect T4|, |S2 intersect T2|).           (2)

The source evaluates both odd relations as literal permutations and,
independently, by nested pointwise letter action.  Commutation is also
checked in both forms.  The complete tables are as follows.

**Ordinary relation sigma tau sigma=tau sigma tau.**

| Typed matrix entries | Number of partners | Nontrivial joint orbit sizes |
| --- | ---: | --- |
| (2,0,0,1) | 480 | (3,6) |
| (2,0,0,2) | 60 | (2,6) |
| (2,2,2,0) | 4 | (6) |
| (4,0,0,1) | 60 | (3,4) |
| (4,0,0,2) | 5 | (2,4) |

There are 609 partners in total.  The last row does not mean equality:
some ordinary pairs have equal typed blocks but different cyclic orders.

**Fifth relation sigma tau sigma tau sigma=tau sigma tau sigma tau.**

| Typed matrix entries | Number of partners | Nontrivial joint orbit sizes |
| --- | ---: | --- |
| (2,0,0,2) | 120 | (2,6) |
| (3,0,0,2) | 24 | (2,5) |
| (3,1,1,1) | 16 | (6) |
| (4,0,0,2) | 1 | (2,4) |

There are 161 partners.  The unique last-row partner is tau=sigma.
Thus a fifth-related pair with the same 4-block is exactly equal.
Also every fifth-related pair has one of the following two alternatives:

    S2=T2,     or     |S4 intersect T4|=3.            (3)

In the second alternative its full matrix is (3,1,1,1), not just an
unlabelled overlap condition.

**Commuting relation.**

| Typed matrix entries | Number of partners | Nontrivial joint orbit sizes |
| --- | ---: | --- |
| (0,0,0,0) | 90 | (2,2,4,4) |
| (0,0,0,2) | 90 | (2,4,4) |
| (4,0,0,0) | 30 | (2,2,4) |
| (4,0,0,2) | 2 | (2,4) |

There are 212 partners.  All tables include orientations and ties
literally.  No relation has been inferred solely from its matrix.

The proof uses only these consequences of the full inventory:

1. Ordinary pairs have 4-block overlap at least two.  If their 4-blocks
   differ, their joint nontrivial orbit sizes are (2,6), (3,6), or (6).
   In the (2,6) case their transpositions are identical.  For every
   ordinary pair, each transposition meets the other's total support.
2. Fifth pairs have 4-block overlap at least two, obey (3), and are
   equal if their 4-blocks coincide.
3. Commuting pairs have no cross-length block intersection, while
   equal-length blocks are equal or disjoint.

These are all-pair statements because the support reduction to twelve
labels is complete.  They are not a sample in the ambient degree.

## 3. Two structural lemmas

The commuting block assertion also has an elementary proof independent
of the table.  A permutation commuting with sigma preserves the unique
4-cycle orbit S4 and the unique 2-cycle orbit S2 setwise.  On S4 it is
a power of the 4-cycle.  The square power would consist of two
transpositions and is impossible for a global permutation of type
(4)(2), which has exactly one transposition.  Hence its restriction is
identity or a 4-cycle on all of S4.  Its own 4-block is therefore equal
to S4 or disjoint from it, and its transposition misses S4.  On S2 its
restriction is identity or precisely that transposition.  The 4-block
cannot meet S2.  This proves all the typed assertions, including absence
of cross-length mixing.

Next let K=<a,b>, and let d commute with a and b.  The permutation d
maps K-orbits to K-orbits of the same size.  In every ordinary pair in
Section 2 the nontrivial orbit sizes are distinct.  Thus d preserves
each nontrivial orbit individually.  Fixed points of K may still be
permuted by d; no conclusion about them is made.

On a transitive K-orbit, any element of its centralizer that fixes one
point fixes all points.  Applying this fact to each power of d shows
that the cyclic group induced by d is semiregular: every orbit of that
cyclic action has the same size, equal to its effective order.  That
order divides four.  In particular:

* On a three-point K-orbit the restriction is identity, since neither
  two nor four divides three.
* On a six-point K-orbit a restriction of order four is impossible.
  A restriction of order two would be three disjoint transpositions.
  Since d has only one transposition and this orbit is preserved,
  that is also impossible.  Thus d fixes the six-point orbit pointwise.

This is the crucial global sidecar.  It works at any ambient degree;
extra fixed points of K outside the displayed orbits do not change it.
On a four-point orbit the restriction could be a 4-cycle, and the proof
does not incorrectly discard that possibility.

There is a further useful consequence of the pair tables.  Suppose
b,c are ordinary-related and c,d fifth-related.  Then d cannot fix the
entire moved support B of b pointwise.  Indeed the ordinary table gives

    |C4 intersect B4|>=2,       C2 intersect B!=empty.

If d fixes B, these imply

    |C4 intersect D4|<=2,       |C2 intersect D2|<=1.

The latter rules out C2=D2, while the former rules out the other
alternative in (3).  This proves the claim with the typed fifth-table
sidecar.  A generic one-third-support bound alone would not pay it.

## 4. The complete four-generator reduction

First suppose A4=B4.  Because a,c commute, their 4-blocks are equal
or disjoint.  The ordinary relation between b,c and A4=B4 forces
C4=A4.  Similarly a,d commute, and the fifth relation between c,d
forces D4=A4.  The equal-4-block row of the fifth table now gives c=d
as permutations, including the transposition.

Since b commutes with d=c and is ordinary-related to c, cancellation
in the ordinary relation gives b=c.  Since a commutes with c=b and
is ordinary-related to b, it also gives a=b.  Thus all four are equal.
This handles every ordinary pair with coincident 4-blocks, including
the non-diagonal orientations in the five-partner ordinary row.

It remains to suppose A4!=B4.  The ordinary table then leaves only

    (2,6),       (3,6),       (6)                     (4)

as joint nontrivial orbit sizes.  The joint-centralizer lemma shows
that d fixes every six-point and three-point orbit.  In the last two
cases of (4), d therefore fixes all of B.  The last paragraph of
Section 3 excludes both cases.

Only (2,6) remains.  Here the exact ordinary row is (2,0,0,2), so

    A2=B2=T,       |A4 intersect B4|=2.               (5)

The common transposition is the two-point joint orbit, while the union
A4 union B4 is the six-point orbit fixed pointwise by d.  We now pay
each step that forces the same transposition for c and d.

1. If C4=A4, then d fixes C4, contrary to the fifth-table requirement
   |C4 intersect D4|>=2.  Thus C4!=A4.
2. Since a,c commute, C4 is disjoint from both A4 and A2=T.  The
   ordinary row (2,2,2,0) for b,c is impossible: that row would put
   both points of B2=T in C4.  Every other ordinary row gives
   C2 intersect B2 nonempty.
3. The commuting a,c pair has its 2-blocks equal or disjoint.  Since
   B2=A2=T, the preceding nonempty intersection forces C2=T.
4. The ordinary b,c relation gives |C4 intersect B4|>=2.  Since d
   fixes B4, it follows that |C4 intersect D4|<=2.  The fifth-table
   alternative (3,1,1,1) is therefore impossible.  The remaining
   alternative forces D2=C2=T.

Thus a,b,c,d all have the **same transposition** on T.  Every one of
them preserves its complement.  Restrict (1) to that complement: the
four resulting permutations are each a single 4-cycle with arbitrary
additional fixed points, and satisfy all six H4 relations.  The proved
[H4 single-cycle equality theorem](planar_jc48_sep08_h4_single_cycle.md)
therefore makes their 4-cycles equal.  Restoring the common transposition
makes a=b=c=d, completing the proof.  In particular no genuinely unequal
4-block case in (4) is realized by a full H4 quadruple.

This is an actual invariant-subset restriction after all four factors
are shown to agree on T.  It is not a claim that cycle extraction is a
homomorphism on arbitrary mixed-type odd pairs.

## 5. Minimal hostiles and the exact geometric consumer

For sigma=(1 2 3 4)(5 6), the following six-label partners are literal
hostiles to premature transposition extraction:

    tau3=(1 5 3 6)(2 4),
    tau5=(1 3 2 5)(4 6).

The pair sigma,tau3 satisfies the ordinary relation, with matrix
(2,2,2,0).  But its isolated transpositions (5 6),(2 4) do not satisfy
that relation.  The pair sigma,tau5 satisfies the fifth relation,
with matrix (3,1,1,1), while (5 6),(4 6) do not satisfy it.  These
hostiles already live on six labels, the minimum possible ambient size
for one permutation of the specified type.  They do not claim to be
full H4 quadruples.  The first failed implication is cycle extraction
from one odd relation.  The repaired form is the common-transposition
restriction established in Section 4.

Conversely the equal quadruple a=b=c=d=sigma is a genuine positive
control for every relation.  It is intransitive, with moved orbits of
sizes four and two, even before extra fixed points are added.  Thus the
theorem excludes transitivity and does not claim that every abstract
representation is trivial.

For an application, suppose a separately proved actual marked supplier
gives four meridians with the relations (1), all of type (4)(2), and
these meridians generate the actual transitive monodromy action.  The
theorem contradicts that transitivity.  In particular it applies to
this cycle-type case of the proved
[mixed-cusp H4 supplier](planar_jc48_sep08_mixed_cusp_braid.md), when
used with its actual meridian and generation hypotheses.  It requires
no full-fixed-point assumption, Euler ledger, retained-sheet count,
or Coxeter involution quotient.  No equality between the actual
complement group and the six-relator presentation is needed.

This closes one additional uniform meridian type.  It does not exclude
other mixed types, construct a global source from an abstract passport,
or settle JC(2).

## 6. Reproduction and finite-control scope

The standalone source imports no previous mathematical implementation.
It enumerates all 83,160 partners, evaluates the literal odd words and
commutators by two separate paths, and records both typed intersection
matrices and joint-orbit sizes.  All table entries are exact counts in
the unfiltered universe.  Its direct joint-orbit computation uses the
permutation edges, not a conclusion inferred from the matrix.

It additionally visits the full ordinary/commuting pair product inside
S_12 and checks every admitted literal common-centralizer action on
its joint orbits.  It checks each two-six propagation case with c
commuting with a, and both explicit split-projection hostiles.  These
are controls for the all-degree proof, not an assertion that a full
H4 quadruple has support at most twelve.  The distinct joint-orbit
sizes and semiregular argument pay that missing ambient information
analytically.

Run from the repository root:

    python3 04-computation/planar_jc48_sep08_h4_mixed42.py
    python3 -O 04-computation/planar_jc48_sep08_h4_mixed42.py

Every gate raises an explicit exception on failure; no Python assert
is used.  The producer's normal, optimized, and frozen output are
byte-identical: **343,032 gates**, 1,263 output bytes.  The exact extra
controls comprise 129,108 common-centralizer trials with 1,892 admitted
edges, 720 fixed-six checks, and 720 two-six transposition-propagation
cases.  Full typed-cell/joint-orbit coupling is checked for every
admitted partner, not only through separate marginal histograms.

Frozen SHA256 pins:

    source (9,443 bytes):
      987f121ef43da7a4ea549c2d4cbd98d7d1c5206a10c0f2504744a2f52b28c9e3
    output (1,263 bytes):
      499481bc6e285848612dff6e35cb00252c45ebb9dbbd435e45b160f180c3e963
    raw labelled semantic trace:
      b911ff8e8201822437bf1c428c07b7445651c9aa37f76b7d980b90728c7a510d

The source, output, and proof are now frozen for independent audit.
The [independent complete audit](planar_jc48_sep08_h4_mixed42_audit.md)
accepts every unbounded-degree step and both343032-gate replays. A
separate support-first/direct-letter/union-find census reproduces all
83160 partners and every coupled typed-intersection/joint-orbit table.
