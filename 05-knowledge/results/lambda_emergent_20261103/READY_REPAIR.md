# An induced obstruction, optimal star repair and a commuting ready grid

**PROVED, scoped, with exact native and independent path-count certificates.**
This is a source-changing construction. It does not preserve H497 or prove
that its ternary children return to a ready source of the required count.

## RR1. Minimum starting multiplicity is an induced monotone

For any nonempty induced S in a strong finite tournament T,

    min_v starts_T(v) >= min_w starts_S(w).

Proof. A Hamilton path P on a proper set U can retain its first vertex
while adding an outside x whenever some vertex of U beats x: insert x after
the last vertex of P beating x. The following vertex, if any, is beaten by
x, so the insertion is legal. Strongness supplies an outgoing edge from
every nonempty proper U. Choose an eligible x by a fixed rule depending
only on U, not on P. Repeated insertion extends every path injectively,
because restriction to its original vertex set recovers it.

For a specified start r outside S, fix a shortest path from r to S, first
hitting S at w. Concatenate it with each path of S starting at w and then
extend. This gives at least starts_S(w) distinct spanning paths starting r.
For r in S use the extension directly. The strong ambient hypothesis is
essential; adjoining a universal source gives zero starts at old vertices.

The regular five-vertex tournament R5 has three Hamilton paths from each
vertex, but only two directed Hamilton cycles. Starting multiplicity, not
cycle count, is the useful obstruction. Independent checks cover every
nonempty induced subset of every labelled strong tournament of order five:
16,864 comparisons in [controls.json](controls.json).

## RR2. Every strong one- or two-ready tournament is R5-free

Ready means F-H in {1,2}, with the exact rooted definition in
[READY sections 1--2](../two_entry_funnels_20261103/READY.md).
For indegree one, deletion of the ready root leaves a strong parent with
one or two paths starting at its sole in-neighbour. RR1 excludes R5 in
that parent. An R5 containing the ready root would already give it two
in-neighbours, also impossible.

For indegree two, the complete exceptional language is

    11^k, 11^k 10, 11^k 01, 11^k 10 01, 11^k 01 10,  k>=1.

Here the k initial chain vertices have identical incidences to everything
outside that transitive block. A regular induced tournament cannot contain
two of them: the earlier one has strictly larger induced outdegree.
Thus any R5 would already occur in the corresponding k=1 core. The five
explicit cores, of orders four through six, contain none; their adjacency
rows and complete subset check are in [ready.json](ready.json). This is a
finite proof of the five cases after the all-size block reduction.
Indegree at least three is excluded by the inherited defect inequality.

Consequently **no induced extension of a tournament containing R5 can
acquire a ready root at any vertex**, however large or differently marked.
Both H497 orientations contain R5. Their neutral apex flip exchanges the
R5 copies rather than removing the obstruction. A change to an old edge
or deletion of an old vertex is necessary; additional births alone cannot
solve this source problem.

## RR3. An optimal native star repair

Take a transitive base v0,...,v(m-1) and two apexes a,b with a->b.
Both apexes beat v0 and lose to v(m-1). Let z be the number of missing
a->vi incidences for i<m-1. Reverse exactly these z arcs, leaving every
other incidence unchanged. The result is one-ready at a.

Indeed a now loses only to the last base vertex ell. The parent T-a is
strong, with Hamilton cycle b,v0,...,ell,b. Its only Hamilton path starting
at ell is ell,b,v0,...,v(m-2). The indegree-one ready formula therefore
gives F-H=1. This also proves optimality among edits restricted to the star
at a, with a required ready endpoint: the parent stays strong, so q>s for
every nonempty in-neighbour set. The exceptional indegree-two ready case
requires q=s and is impossible. The final root must have indegree one;
reducing its original indegree z+1 takes at least z arc reversals.

Every reversal has an explicit native implementation. For a missing
i->a, the path

    (i,a,b, all remaining base vertices in increasing order)

is Hamiltonian. Append a helper beating only i. Complement off the path
with helper last, then off the path with helper inserted after i. Their
difference reverses exactly i->a; delete the helper. This costs two native
moves and one helper per arc, retaining every old name and both apexes.
The guard path selects an actual old i; this is a full-source controller,
not an assertion that a forgotten vertex can be re-marked from six counts.

For the two inherited H497 presentations:

|Input|One reversed arc|Ready root profile (h,s,q,F)|NR output H|
|---|---|---|---|
|a->b|v2->a|(395,68,69,396)|1187|
|b->a|v3->b|(311,60,61,312)|935|

There are exactly two R5 copies in each source. The one-edge repairs remove
both, and are minimal since an unchanged source cannot be ready. Literal
rows, both guarded paths, retained names and 170 small-family normalizations
are checked by [ready.py](ready.py). The plus repair changes H by -102;
the minus repair changes H by -186. They are not H-neutral resets.

## RR4. An infinite one-hole family closed under two preparation operations

Let p>=1, q>=0. Give the dominant apex and other apex the words

    w=1^(p+q+1)0,    u=1^p 0 1^q 0,    a->b.

The order is p+q+4. These sources are all one-ready at a and have

    H(p,q)=2^q(2*3^p+6*2^p-3)-1-[q=0].

In particular H(3,2)=395 and H(2,3)=311 are the repaired H497 sources.
Other actual values include H(1,3)=119=7*17, H(2,1)=77=7*11 and
H(3,3)=791=7*113. Their one-ready resets yield 359,233,2375 respectively;
the first two and1187 are prime by exact division. No infinite-prime or
full prime-lane coverage conclusion follows.

For q>=1, H(p,q+1)=2H(p,q)+1. At q=0 the offset is3 instead.
The entire family has rational bivariate generating function

    sum_(p>=1,q>=0) H(p,q) z^p w^q
      = [6z/(1-3z)+12z/(1-2z)-3z/(1-z)]/(1-2w)
        -z/[(1-z)(1-w)]-z/(1-z).

It also retains the full golden observer. Put e=ends_a(1), u2=starts_a(2).
For q>=1,

    e=2^p+1,
    u2=6*3^(p+q)+4*3^p*2^q-2*2^p*3^q+6*2^(p+q)
       -3*2^p-3*2^q-3.

For q=0, e=2^p and u2=10*3^p+2^p-6. Then X=(u2-e)/2,Y=(u2+e)/2.
The repaired golden triples are (395,945,954) and (311,839,844).

Proof: in a fixed order of the apexes, distribute the base vertices among
the gaps before, between and after them, increasing within every gap.
For each gap retain empty, occupied with bad right endpoint, or occupied
with good right endpoint. Summing the two apex orders counts H. For u2
there are two ordered-bin patterns: (a,b) in the first bin and no apex in
the second, or a in the first and b in the second; the gap before a is
forbidden. For e the only pattern is (b,a), with the final gap forbidden.
These give respectively three, three and two available gaps per pattern.

For a common11 letter the gap insertion maps are commuting idempotents.
Their sum A is annihilated by product_(j=0..g)(A-jI), g<=3. The word is
A^p B A^q C, for hole10 and final00. Since p>=1, and separately q>=1,
the zero eigenmode vanishes. Exact values at p,q in {1,2,3} uniquely give
the formulas; the q=0 boundary is fitted in p only. The operator identity
proves all orders. [one_hole.py](one_hole.py) builds each matrix, checks the
identity, recovers the coefficients and compares independent endpoint DP.

## RR5. Two ports to compile; three to support both named updates

A fresh two-port compiler retains (a,b), appends all base vertices in
order, and forgets each newborn. Uniform old-to-new polarity and masks
3,1,3,0 realize the runs11^p,10,11^q,00 respectively.

For an existing source, keep three actual names (a,b,ell). The p-update
inserts a new first base vertex, losing only to a,b. The q-update inserts
one just before ell, beating only ell. In the two-polarity birth convention
both have named mask3; only the uniform polarity changes. They commute up
to exchanging the two fresh labels, and retain readiness. Literal graph
checks, not only equal H values, are saved in [one_hole.json](one_hole.json).

Three names are necessary for these **specified one-birth updates**, with
no retrieval of forgotten vertices: the p-update requires exceptional
names a,b, while the q-update requires ell. At order at least five, choosing
the opposite polarity cannot replace either exceptional set by at most
two names. Any live set supporting either next update must contain all
three. This is a register-liveness lower bound for this instruction set,
not a lower bound for arbitrary native gadgets or all encodings.

The distinction is constructive: a short source template can be compiled
with fewer ports than are needed for arbitrary future edits to its output.
The ready grid still does not contain its ternary children as a proved
closed return family; that remains the next source-level obligation.

## RR6. A commuting repair cube whose peak clock does not commute

On the larger two-apex family with common first1 and last0, let P be the
star normalization, U prepend a common11, and V insert a common11 immediately
before the last00. At the source-template level

    P^2=P,       PU=UP,       PV=VP,       UV=VU.

These identities hold for actual graphs with the indicated fresh-label
transport: P changes only the dominant apex's missing base arcs; U and V
add common11 incidences and do not alter any such missing arc. Thus this
is an action of N^2 times the two-element idempotent monoid. The image of
P is exactly the all-one-dominant ready family considered in RR3.

For fixed base size m and fixed other-apex word, P has2^(m-2) labelled
source presentations in its fibre. Exactly binomial(m-2,z) use z helpers
and2z native moves. Its exact helper-clock generating polynomial is
(1+y)^(m-2). These are labelled word fibres, not unlabelled tournament
isomorphism-class multiplicities.

There is a useful operational failure of naive commutation. For a source
of order n with z>0, all six orders of P,U,V produce the same graph, spend
2z native moves and add two permanent vertices. But putting P last requires
peak order n+3, whereas the other four schedules need only n+2. Its helper
coexists with both new vertices only in the former schedules. The literal
traces in [clock_cube.json](clock_cube.json) demonstrate both cases.

This gives a concrete evaluation-order connection to the lambda-cube cue:
source confluence does not identify resource effects. The retained clock
must say which resource is observed. No claim about lambda normalization
or computational complexity follows just from this commuting diagram.
