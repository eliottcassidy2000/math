# Every H4 generator of type (5)(2) is the same permutation

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**

This is an abstract, uniform permutation theorem. Its finite certificate
is a complete pair inventory, not a bound on the ambient degree of a
four-generator representation. The geometric application requires the
separate proved actual marked-meridian supplier and actual transitive
generation. JC(2) remains OPEN.

## 1. Statement, inheritance, and the useful change of mechanism

Let a,b,c,d be permutations of an arbitrary finite set Omega. Suppose
each has one 5-cycle, one disjoint transposition, and otherwise only
fixed points. Assume the six H4 relations

    aba=bab,        bcb=cbc,        cdcdc=dcdcd,
    ac=ca,          ad=da,          bd=db.                 (1)

**Theorem.** Then a=b=c=d. Consequently these four
permutations cannot generate a transitive action: their common
permutation already has distinct nontrivial orbits of sizes five and
two, with any additional fixed labels separate.

The closest proved mechanism is
[h4_mixed42](planar_jc48_sep08_h4_mixed42.md), which retains both
cycle blocks through its odd relations and only removes a factor after
proving a common invariant subset. The
[mixed32](planar_jc48_sep08_h4_mixed32.md) and
[mixed33](planar_jc48_sep08_h4_mixed33.md) results supply warnings about
cycle extraction and block ambiguity. The present unequal block sizes
give a canonical length gauge, but that does not itself make cycle
extraction preserve any Artin relation. Two minimal seven-label
hostiles below explicitly refute that implication.

The corrected near miss is to copy the fixed-prefix proof for (4)(2),
or to assume equal five-point supports imply equal cycles. The new
ordinary pair table gives strict majority, at least three points,
inside a five-cycle. Together with the actual commuting relations this
forces a common five-set much earlier. The final restriction is made
on that proved invariant set and its complement. No inherited
single-cycle equality theorem is required.

The five live concepts are typed cycle supports, literal odd relations,
strict majority intersections, centralizers of a unique cycle, and
restriction to an actual invariant subset. The least-used sidecar is
the distinction between coincident supports and coincident cyclic
orders. The pair-to-matrix map preserves typed overlap counts and
loses cyclic order; the producer evaluates the original words before
taking that map. The four-generator-to-single-block map is used only
after every generator is proved to preserve the same set.

Write X5 and X2 for the 5-cycle support and transposition support of x.
These names are determined by the cycle lengths, with no ordering
choice between equal blocks.

## 2. The exact unfiltered pair universe

Fix

    sigma=(1 2 3 4 5)(6 7)

on fourteen labels. For any pair of the declared type on an arbitrary
finite Omega, its support union has at most fourteen labels. Restrict
to that union, conjugate the first member to sigma, and pad with fixed
points. This gives exactly a pair represented in the universe below.
Both original words and both typed supports are preserved by this map.
Conversely every partner in that universe is tested without filters.

Choose its five-set, put its smallest label first in the cyclic order,
take all 24 orders of the other four labels, and choose an unordered
two-set from the remaining nine labels. Every permutation occurs once:

    binom(14,5)*24*binom(9,2)
      =binom(14,7)*504
      =1,729,728.                                      (2)

The source streams the entire universe, comparing composed permutation
words with independent direct nested letter action for the ordinary,
fifth, and commuting relations on every partner. Its raw digest uses
fourteen image bytes followed by the three relation flags. No pair is
discarded by an overlap or support hypothesis before this test.

Order the matrix entries by lengths 5,2 in both rows and columns:

    (|S5 intersect T5|, |S5 intersect T2|,
     |S2 intersect T5|, |S2 intersect T2|).             (3)

The complete tables, including the coupling to nontrivial joint orbit
sizes, are as follows.

**Ordinary relation.**

| Typed entries | Partners | Nontrivial joint orbit sizes |
| --- | ---: | --- |
| (3,2,2,0) | 10 | (7) |
| (4,0,0,1) | 420 | (3,6) |
| (4,0,0,2) | 35 | (2,6) |
| (5,0,0,1) | 84 | (3,5) |
| (5,0,0,2) | 6 | (2,5) |

Total: 555. In particular every ordinary pair satisfies

    |S5 intersect T5| >= 3.                           (4)

**Fifth relation.**

| Typed entries | Partners | Nontrivial joint orbit sizes |
| --- | ---: | --- |
| (2,1,1,0) | 2,100 | (10) |
| (4,0,0,2) | 105 | (2,6) |
| (4,1,1,1) | 10 | (7) |
| (5,0,0,2) | 6 | (2,5) |

Total: 2,221. In particular every fifth-related pair satisfies

    |S5 intersect T5| >= 2.                           (5)

The last row contains six cyclic orientations; it does not say that
the pair is equal. That distinction is explicitly retained in the
same-support hostile below.

**Commuting relation.**

| Typed entries | Partners | Nontrivial joint orbit sizes |
| --- | ---: | --- |
| (0,0,0,0) | 504 | (2,2,5,5) |
| (0,0,0,2) | 504 | (2,5,5) |
| (5,0,0,0) | 84 | (2,2,5) |
| (5,0,0,2) | 4 | (2,5) |

Total: 1,096. This count is independently recovered from the
centralizer: the two outside-five-block cases each have
binom(7,5)*24=504 possibilities; a common five-block with an outside
transposition has 4*binom(7,2)=84 possibilities; both supports common
give four nonidentity powers of the original 5-cycle.

Only the all-pair bounds (4),(5) and the commuting block assertion
below are needed for the global equality argument. The remaining
inventory and orbit information are retained as independently reusable
sidecars. They are not interpreted as a census of full H4 tuples.

## 3. The elementary commuting sidecar

If x,y commute, y sends each orbit of x to an x-orbit of the same
size. Since x has exactly one orbit of size five, y preserves X5.
Its restriction there commutes with a 5-cycle, hence is a power of
that cycle. Primality of five says that power is either identity or
itself a 5-cycle. Thus

    Y5=X5 or Y5 intersect X5 is empty.                (6)

Its transposition cannot meet X5. The unique two-point orbit X2
is preserved too; y acts there as identity or its own transposition.
Consequently the length-five and length-two blocks do not cross,
and their respective supports are equal or disjoint. The complete
commuting table confirms this argument literally.

One useful but unneeded extension of this sidecar concerns a common
centralizer of an ordinary pair. Its table has distinct nontrivial
joint orbit sizes, so a common centralizer preserves each such orbit.
On a transitive joint orbit a centralizer element fixing one point
fixes all points. Applying this to each power shows that its effective
cyclic action is semiregular. A global permutation of type (5)(2)
has a nontrivial semiregular restriction only on a set of size exactly
five or exactly two: every nontrivial cycle on such a preserved orbit
would need the same length. It therefore fixes all ordinary joint
orbits of sizes three, six or seven. Additional singleton orbits can
still be permuted; nothing is claimed about them. The proof below
does not need this extension.

## 4. Strict majority forces one actual common five-set

Apply (4) to a,b and b,c. If C5 were disjoint from A5, then

    |C5 intersect B5|
       <= |B5 minus A5|
       =5-|A5 intersect B5|
       <=2.

This contradicts the ordinary b,c bound at least three. Since a,c
commute, (6) leaves no intermediate overlap alternative. Therefore

    C5=A5.

Since a,d commute, D5 is equal or disjoint from A5. But (5) for c,d
says that D5 meets C5=A5, so

    D5=A5.

Finally b,d commute, making B5 equal or disjoint from D5=A5. The
ordinary a,b bound excludes disjointness. Thus

    A5=B5=C5=D5=T.                                   (7)

This argument is uniform in Omega. It uses no bound on the union of
all four supports. The support reduction to fourteen labels was used
only to establish the two all-pair facts.

The common set T in (7) is preserved by all four permutations, and so
is its complement. Restrict every actual relation (1) to T. Because
a commutes with c and d, the restrictions c|T and d|T are powers of
the same 5-cycle a|T, hence commute. In a commuting pair x,y, the
relation xyxyx=yxyxy reduces to x^3 y^2=x^2 y^3 and forces x=y.
It follows that

    c|T=d|T.

On the complementary invariant subset, c and d are each exactly one
transposition, with arbitrary further fixed points. For involutions
x,y the fifth relation is equivalent to (xy)^5=1. The product of two
transpositions has order one, two, or three, according as they are
equal, disjoint, or meet in one label. Only order one divides five.
Thus their restrictions to the complement are equal as well, and

    c=d

as full permutations of Omega.

Now b commutes with d=c and satisfies bcb=cbc. Commuting cancellation
in that ordinary relation gives b=c. Then a commutes with c=b and
satisfies aba=bab, giving a=b. This proves a=b=c=d.

This is the promised repair of cycle extraction. Both restrictions
are legitimate because (7) was proved first. The argument never
claims that taking the 5-cycle or transposition component is a
homomorphism on arbitrary odd-related mixed permutations.

## 5. Minimal hostiles, positive controls, and actual scope

On seven labels, the following two partners of
sigma=(1 2 3 4 5)(6 7) are exact hostiles:

    tau3=(1 6 4 7 2)(3 5),
    tau5=(1 4 3 6 2)(5 7).

The first satisfies the ordinary relation with typed entries
(3,2,2,0), but its extracted transposition does not satisfy that
relation with (6 7). The second satisfies the fifth relation with
entries (4,1,1,1), while the extracted transpositions (6 7) and (5 7)
do not satisfy it. Seven is the minimum possible ambient size for
a permutation of this type. These are individual-pair hostiles,
not full H4 representations.

There is a further seven-label warning:

    tau_same=(1 2 4 5 3)(6 7).

It satisfies the fifth relation with sigma and has the same two
typed supports, but is not sigma and does not commute with it.
Therefore same support is not equality. The global commutation
on the common five-set in Section 4 is load-bearing.

Conversely a=b=c=d=sigma satisfies every relation. Its common action
has nontrivial orbits of sizes five and two, so it is intransitive even
before fixed padding. The conclusion is equality, not triviality.

For an actual geometric consumer, assume a separately proved marked
supplier gives four meridian permutations satisfying exactly the
relations (1), all of type (5)(2), and assume those actual meridians
generate the transitive action. The theorem contradicts transitivity.
In particular it can be used on this type in the proved
[mixed-cusp H4 supplier](planar_jc48_sep08_mixed_cusp_braid.md) with
its declared meridian and generation hypotheses. It does not identify
every abstract representation with an actual Keller map. It removes
one moved-support-seven cycle type; no numerical degree-floor upgrade
is asserted here without separately checking the remaining types.

## 6. Exact reproduction and proof boundary

The [standalone source](../../04-computation/planar_jc48_sep08_h4_mixed52.py)
imports no inherited mathematical implementation. It checks the whole
unfiltered pair universe, literal words by two paths, the complete
typed/orbit tables, commuting restrictions, the exact transposition
pair universe on at most four labels, all nonidentity powers on the
five-set, the strict-majority integer inequalities, both extraction
hostiles and the same-support hostile. The finite terminal controls
verify the elementary statements used in the written proof; they do
not substitute for the proof that T is invariant at arbitrary degree.

The raw stream is fixed-width: for each partner in the specified
deterministic construction order, its fourteen image bytes precede
its three Boolean relation flags. This pins every rejected pair as
well as every admitted pair. Orbit sizes are computed only after a
relation is found, since rejected pairs are not part of any claimed
orbit inventory. The source records an additional complete typed-table
digest. No ambient-degree quadruple enumeration is used.

Run from the repository root:

    python3 04-computation/planar_jc48_sep08_h4_mixed52.py
    python3 -O 04-computation/planar_jc48_sep08_h4_mixed52.py

All gates are explicit exception checks, active in both modes. The
normal and optimized producer replays are byte-identical to the frozen
output: **1,749,134 gates, 1,155 output bytes**.

* Source: 8,729 bytes, SHA256
  bb0bb94252ab3e818e26172d2bc224ec1039964fd0e9f290a9cb6057ada301b9.
* Output SHA256:
  6ccb824360d4c5416b88d1fead54503b58e00adea3e412ee0efbccb8c1b5e175.
* Raw fixed-width pair/flags SHA256:
  f8dac4bf55a921332df0c02c5b262c067ac9747e89c807d6d937e7596533d801.
* Complete typed-table SHA256:
  cb623d547de2fe30b050373c9bb5c2acc215156dcc1563451cb189d46ca6bc81.

The [independent full audit](planar_jc48_sep08_h4_mixed52_audit.md)
accepts the uniform majority/common-invariant-set proof and every
source/replay gate. A separate transposition-first enumeration with
union-find joint orbits reconstructs all1,729,728 partners and every
typed table exactly. The theorem is accepted into the proved graph.
