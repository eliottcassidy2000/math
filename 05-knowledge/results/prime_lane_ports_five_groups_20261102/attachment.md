# Two named births detect the circulation forgotten by one-step response

**PROVED:** the complete two-birth detector and isomorphism criterion for
oppositely oriented named triangles, the triangle-module one-birth
invisibility theorem, and the smallest nonclosure boundary below.
**INHERITED PROVED:** the quadratic response/cut law and the full order-three
mean-response dimension seven. **FINITE-EXACT:** independent controls of
all labelled parents through order five and all 128 two-birth triangle
attachments. This does not settle the prime or 7p realization lanes.

## Inheritance and the exact target

The closest mechanisms are
[THM-4097, exposed-slot response](../../../01-canon/theorems/THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension.md),
[Session124, canonical weighted-cut law](../session124_tournament_cut_response_20261031/README.md),
and [Session136, positive-cut circulation quotient](../session136_full_tournament_spectrum_20261102/ear_geometry/positive_cut_geometry.md).
The last explicitly identifies opposite named C3 orientations as equal
one-step cut functions. We extend that hostile to a complete two-birth
classification, rather than claiming the first collision as new.

The corrected near miss is that an exact one-step observation should update
itself under arbitrary future insertions. The least-used sidecar is the
**ordered** endpoint pair, with its actual old vertex names. The live board
is `one-step polynomial / ordered two-path cover / invisible circulation /
two-port signature / exact prime target`. The anchor is actual prime/7p
construction; the niche is a smallest exact failure of a proposed state
compression; the wildcard is its determinant detector.

The full fixed-attachment linear carrier is a separate object. Its
order-three dimension seven is already proved in
[the forest-response theorem, §5](../primegame_port_spectrum_20261103/FOREST_RESPONSE_RECURRENCE.md)
using 21 S/B words of lengths at most two. The
[pooled-parity carrier](../session145_endpoint_fence_ramsey_20261103/carrier/POOLED_PARITY_CARRIER.md)
and [Session146's order-ten correction](../session146_carrier_primegame_20261103/carrier10_upper/PROOF.md)
must be read together: the full dimensions at seven through ten are
2668,13028,65908,345106, and the 2520-dimensional order-ten kernel
propagates. No all-order minimality assertion is made here.

## 1. The exact one-step observable and what it forgets

Fix a named tournament T. Let `s_i,e_i` count Hamiltonian vertex orders
starting/ending at i, and let `Q_ij` count ordered nonempty two-path covers
of all old vertices, first path ending at i and second starting at j.
Equivalently, Q counts full old orders with an exposed consecutive i,j slot
and all other adjacencies valid. Q is ordered; it is generally asymmetric.

For an actual new vertex x use `a_i=[i->x]`. Deleting x, including the
possibly invalid old joining slot, gives the inherited identity

    H(T+x)=H(T)+sum_i(e_i-s_i)a_i
              +sum_(i!=j) Q_ij a_i(1-a_j).              (1)

Put `w_ij=(Q_ij+Q_ji)/2`. The endpoint-divergence cancellation of Session124
gives the exact original-orientation cut law

    H(T+x)=H(T)+2 sum_(i->j in T) w_ij (1-a_i)a_j.       (2)

All weights are positive integers and `w_ij>=H(T-{i,j})`. These facts are
inherited, including the orientation-supported lift; they are not a new
weighted-network realization theorem.

For `S={i:x->i}`, write `F(S)=H(T+x_S)` and
`d_i^w=sum_(i->j)w_ij`. The canonical coefficient data for **this one-step
observable** are `(H,w,d^w)`, with `sum_i d_i^w=sum_(i<j)w_ij`. They can be
recovered from the empty, singleton, and two-element cuts:

    H=F(empty),
    d_i^w=(F({i})-H)/2,
    w_ij=(F({i})+F({j})-F({i,j})-H)/2.                 (3)

Thus two named sources have identical one-step responses iff these data
agree. This is the coefficient quotient, not a claim that all formal data
are tournament-realizable or linearly independent on the actual source
set. Arbitrary positive weighted cuts can falsely suggest H7 at order
four: take H3 and triangle weights `(2,1,1)`, whereas every actual
three-vertex source has all pair weights one.

Matrices R having the same cut function can differ by a skew zero-divergence
circulation. The OCF path-network lift and the original-orientation lift
are two specified representatives. Passing to the observable forgets this
circulation; it does not justify deleting the ordered Q coordinates from
a state intended for future updates.

## 2. Triangle modules are invisible to one named birth

**Theorem A.** Let three named old vertices form a cyclic triangle module
in any tournament: each vertex outside the module sees all three uniformly.
Reverse the triangle, leaving every other arc and name unchanged. For
every single fixed new-vertex attachment, the two resulting tournaments
are isomorphic by a permutation fixing every vertex outside that triangle,
including the new vertex.

**Proof.** The three module vertices receive only two possible adjacency
colors to the new vertex, so two receive the same color. Swap those two
vertices. A transposition reverses the cyclic triangle; it preserves all
outside adjacencies by modularity and preserves the new incidences by the
equal-color choice. This gives the stated isomorphism. QED.

For the bare triangle, let `T_+=(2,4,1)` and `T_-=(4,1,2)` in outgoing-row
bitsets on 0,1,2. Then `H=3`, `s=e=(1,1,1)`, and the complete one-step
response is 3 at the two constant cuts and 5 at the other six cuts. Their
ordered Q matrices differ: Q is 2 on the source triangle's arcs and zero
on its reverse arcs. Their symmetric pair weights are all one.

This is named-port ambiguity, not two different ordinary parent classes.
An isomorphism can repair it only if the future ports are transported by
that same isomorphism. A fresh independent renaming or resampling would
change the authorized source and its fixed attachment.

## 3. Complete two-birth detector on the bare triangle

Add two named vertices x,y. Specify all seven new incident arcs by

    a_i=[i->x],  b_i=[i->y]  (i=0,1,2),  c=[x->y].

Let U+ and U- be the five-vertex children of T+ and T- with these **same**
seven bits. Define

    Delta(a,b)=det [1 a_0 b_0; 1 a_1 b_1; 1 a_2 b_2].

**Theorem B.** For every one of the 128 attachments,

    H(U+)-H(U-)
      =2 Delta(a,b) [1-2c + sum_i b_i-sum_i a_i].        (4)

Moreover U+ and U- are isomorphic iff the three two-bit signatures
`(a_i,b_i)` are not all distinct. Distinct signatures can nevertheless
give equal H: this happens exactly when the second factor in (4) vanishes.

**Proof.** If two signatures coincide, their transposition gives the
isomorphism from Theorem A while fixing both new vertices. The determinant
and count difference vanish.

Otherwise the three signatures are three different corners of the binary
square. Permuting old names multiplies both Delta and the difference by
the permutation sign: odd permutations exchange T+ and T-. It therefore
suffices to assign the three signatures in increasing lexicographic order.
There are four missing corners and two values of c. Literal path counting
in these eight fixed tournaments gives:

| Missing corner | c | H(U+) | H(U-) | Delta |
|---|---:|---:|---:|---:|
|(0,0)|0|13|11|1|
|(0,0)|1|11|13|1|
|(0,1)|0|9|9|1|
|(0,1)|1|11|15|1|
|(1,0)|0|11|15|-1|
|(1,0)|1|9|9|-1|
|(1,1)|0|11|13|-1|
|(1,1)|1|13|11|-1|

These entries can be verified by listing the 120 vertex orders, retaining
an order exactly when its four adjacent arcs are valid. The independent
checker also computes each entry by subset DP. The missing corner `(u,v)`
satisfies `u-v=sum b_i-sum a_i`, so the table proves (4).

In six table rows, unequal H already proves nonisomorphism. In the two
equal-H rows, the multisets of Hamilton counts after deleting one vertex
are respectively `{1,3,3,5,5}` and `{1,1,5,5,5}` (possibly exchanged).
These are isomorphism invariants, so those children too are nonisomorphic.
This exhausts all cases and proves the class criterion. QED.

The complete partition of attachment space is:

| Attachments | Ordinary child classes | H |
|---:|---|---|
|80|equal (repeated old port signature)|equal|
|12|different (three distinct signatures)|both 9|
|36|different (three distinct signatures)|different|

The difference histogram is `0:92`, `2:12`, `-2:12`, `4:6`, `-4:6`.
This separates three notions: a port signature distinguishing old vertices,
a resulting ordinary child class, and the scalar target H.

### Minimal nonclosure and an actual prime-target hostile

First add x with old-to-new mask 1. Then add y with old-to-new mask 2,
including bit 3 for x, which is zero here (`y->x`). The exact children are

    U+=(10,20,1,6,13),  H(U+)=13;
    U-=(12,17,2,6,13),  H(U-)=11.                      (5)

Both first children have H5 and are isomorphic as ordinary four-vertex
tournaments, but their named future response data differ. The old sources
had identical **complete** one-step response polynomials. Consequently no
deterministic update from that reduced state and the fixed new cut can
produce the next complete one-step state correctly for all sources.
The full ordered `(s,e,Q)` states do distinguish the initial triangles;
the collision must not be attributed to those full states.

This is the smallest old order at which the reduced-state failure can
occur. At orders zero/one there is only one named tournament; at order
two the two orientations have distinct singleton responses. Hence an
identical complete one-step state there determines the old orientation.

Equation (5) is a source-specific prime realization test: a controller
seeking 11 or 13 cannot erase the old chirality while keeping these actual
port names. It is not a construction of every prime. One may pad the
entire gadget by transitive sources/sinks to propagate the failure to
larger old orders; this preserves both scalar values and vertex transport.

## 4. Four one-birth dimensions versus seven full-response dimensions

For arbitrary rational source laws on the eight labelled tournaments of
order three, the span of all fixed **one-new-vertex** mean responses has
dimension exactly four. Indeed all pair weights equal one, so

    F_T(S)=H(T)+2 sum_(i in S)d_i^+(T)-|S|(|S|-1).

Together with `sum_i d_i^+=3`, the response span is contained in
`span{1,H,d_0^+,d_1^+}`. Conversely the empty cut gives H, singleton cuts
give the degrees, and a singleton/complement pair gives the constant
function because its sum is `2H+4`. These four functions are independent:
on the six transitive orientations H=1 while the degree triples run over
all permutations of `(0,1,2)`, forcing any vanishing linear relation's
degree coefficients to be zero; the cyclic orientations then have H=3,
forcing the remaining coefficients to be zero. The receipt independently
checks exact rank four.

Two new vertices already attain the inherited full fixed-attachment
dimension seven. Here is a compact independent certificate: in the bit
convention of (5), use the seven context pairs

    (0,0),(0,1),(0,2),(0,3),(0,9),(0,10),(1,2).

Each pair is the first birth's three-bit mask followed by the second
birth's four-bit mask. Against old tournament masks 0,...,7 in the
lexicographic pair order `(01,02,12)`, where a bit reverses the increasing
arc, the evaluation matrix is

```text
1 1  3 1 1  3  1  1
3 5  9 9 3  9  5  9
5 3  9 3 9  9  9  5
5 5 13 9 9 13 11 11
1 3  5 5 1  5  3  5
3 1  5 1 5  5  5  3
5 5 13 9 9 11 15 13
```

The first seven columns have determinant `-128`, proving rank at least
seven over Q. Conversely no Hamiltonian vertex order in any fixed
attachment can contain all three old triangle edges. Its old-edge
multilinear polynomial thus has degree at most two, with at most seven
monomials. This gives the all-attachment upper bound seven directly.

Thus two births are necessary and sufficient to attain the full order-three
linear response span. The dimension seven and the existence of length-two
query banks are inherited; this seven-query minor is a smaller explicit
certificate and the four-dimensional one-birth comparison explains the
nonclosure. This does not assert source-law reconstruction: the old
three-edge Walsh character remains invisible to every fixed attachment.
Nor does it claim that seven is an information-theoretic bit cost or a
minimal dimension at every old order.

## 5. Reproduction and remaining realization interface

Run in this directory:

```sh
python3 -I -B attachment_probe.py > attachment_receipt.json
python3 -I -B -O attachment_probe.py > /tmp/attachment-probe-optimized.json
cmp attachment_receipt.json /tmp/attachment-probe-optimized.json
```

The script imports no producer or repository algorithm. It computes Q by
literal exposed-word enumeration and child H by subset DP, checks all
33866 cuts of all1099 labelled parents of orders1..5, all128 two-new
attachments with an additional literal Hamilton count, the exact
eight-case class criterion and deletion decks, the seven-query determinant,
and module-preserving transposition controls. No strongness or isomorphism
filter removes a source. Checks are explicit exceptions, active under `-O`.
The normal and optimized outputs agree byte for byte.

Connection contract: source = a named actual tournament with ordered path
endpoints; target = its complete attachment response, then a specified
prime target; map = fixed arcs to the new vertices, retaining common old
names; preserved = actual tournament realizability and exact H; lost by
the one-step quotient = directed circulation; repair = joint port
signatures/ordered endpoint state, or a proved closed carrier for the exact
operation language. The cheapest decisive test is (5). None of the
formal response span, positivity of cut weights, or a source-law mean
provides an actual source realizing every requested prime or7p.
