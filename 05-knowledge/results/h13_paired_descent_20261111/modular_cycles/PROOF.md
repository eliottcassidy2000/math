# A paired source reset: move the unit defect into an endpoint, then restore it

**2026-11-11. PROVED IN NOTE:** an actual asymmetric paired birth, its
all-address response, exact deletion floor25, and a common-continuation
cost comparison. **FINITE-EXACT:** the explicitly bounded source searches.
**OPEN:** repeatable target selection and Chan–Pak. Inherited modular and
Monster inputs are used for their operation ideas, not reverified or
promoted as tournament theorems.

The concrete new operation is a simultaneous same-count replacement of
the H13 family's211/213 twins. It reduces their order from9 to8, raises
their minimum deletion count from13 to25, and lowers their U coordinate
from347 to269. After every nonempty common O/Z word it produces strictly
smaller counts, by a known amount, while preserving the internal unit gap.

## 1. The recovered modular idea that suggested the operation

The old fifteen-prime story has several resolutions. The most useful here
is the distinction between **rational endpoints** and **a rational edge
joining them**. The recovered level-two supersingular graph has, at the
four primes29,41,59,71, rational j-endpoints and parallel isogeny edges
exchanged by Frobenius. A vertex-only observer loses which lift was used.
The Monster7A normalizer gives a second version: two local actions must
have the same image in one common C2 quotient before their independent
kernel lifts can be paired.

These sources led to an actual construction test. Given two rooted
tournaments with a unit H difference and equal `(E,U,d)`, do not demand
the **same** new-vertex cut on both sides. Search the fibre product of
their two cut sets over the common future sidecar. The two cuts may be
different lifts, just as two compatible local choices need not be equal.

The inheritance is the H13 golden pair from
[the preceding history continuation](../../five_core_web_20261111/history/continuation/PROOF.md).
Its root state is `(H,E,U,d)`: H is the Hamilton path count, E the number
ending at the root, U the weighted root-start two-cover count, and d the
root-deletion count. The common native maps are inherited:

    O(H,E,U,d)=(H,0,2U,d),
    Z(H,E,U,d)=(H+U-E,d,2U-d,d).

The canonical hostile is the sparse image of common O/Z continuation from
one source. The corrected near miss is independent endpoint choices with
no compatible shared lift. The underused sidecar is the **difference
vector of two actual rooted states**.

## 2. An actual birth moves the difference from H to E

Use the old five-vertex core

    Q=(18,8,3,5,14), H(Q)=13.

At left/right addresses0, its two apex cuts6 and3 give rooted states

    (31,7,45,13),       (33,7,45,13).                 (1)

Append y on each side. It loses to the retained root. On the first side
it beats core vertices1,2,3; on the second it beats only core vertex3.
All other comparisons with Q point toward y. The resulting deletion cores,
in order `Q,y`, are

    C-=(50,8,3,5,46,14),        root cut38,
    C+=(50,40,35,5,46,8),       root cut35.            (2)

These are actual tournaments. Both deletion cores have H25. The rooted
states after this **asymmetric** paired birth are exactly

    (73,9,147,25),      (73,7,147,25).                (3)

Thus the original difference `(2,0,0,0)` has become `(0,-2,0,0)`.
The scalar counts have coalesced, but the endpoint difference remains.
One common Z restores the unit H difference:

    (211,25,269,25),    (213,25,269,25).              (4)

One common O instead sends (3) to identical states `(73,0,294,25)`.
Every further O/Z continuation then has identical H on both sides.
This supplies a literal operation for hiding and reactivating a difference
in a coordinate that scalar H forgets.

The finite constants can be computed directly by weighted path enumeration.
In the inherited signature order `(d,k+,k-,i+,i-,o+,o-,j+,j-)`, (2) gives

    C-: (25,279,-25,107,-15,122,-16,43,5),
    C+: (25,279,-25,104,-14,122,-18,44,4).            (5)

In particular the difference of full two-flank outputs is

    H+(A,B)-H-(A,B)=-3A-(A mod2)+2(B mod2).          (6)

The zero-left condition in the next section is therefore essential.

## 3. The entire right-address family, including its parity gate

For arbitrary `B>=0`, define

    M(B)=25g(B)+122B+16(B mod2)+48.                  (7)

Then the actual graphs (2), with left address0 and right address B, have
counts

    M(B),        M(B)+2(B mod2).                     (8)

If `B>0` has b bits and reversal `r_b(B)`, both rooted sidecars are

    E=25,        U=147*2^b-25r_b(B),       d=25.     (9)

There are two independent ways to understand the formula. Substituting
(5) in the inherited compiler proves (7)-(8). Alternatively start from
(3), read B from least significant bit to most significant bit, and apply
0=O,1=Z. The first bit either kills the pure E difference or moves it to
H; all subsequent steps preserve the resulting pure H difference. This
proves the exact factor `B mod2` without extrapolating a bounded census.

Equation(9) follows by the inherited recurrence for U; a nonempty canonical
word ends in Z, so E=d. The same mechanism gives the old H13 formula.

There is also a literal all-B description of the paired birth. Start from
the original H13 pair with right tail B. Let y beat every tail vertex and
have the two prescribed Q cuts above; let the root beat y. Move y before
the tail in the vertex order. The resulting graphs are precisely (2)
with the same right address B. At B1 the birth masks in the old label
order are46 and40; mask42 on the second side gives the same rooted
response. The formula retains the actual cut, not just a desired count.

## 4. A sharp increase in the actual source boundary

**Theorem1.** For every `B>=1`, each graph in (8) has minimum
one-vertex-deletion H exactly25.

Deleting the root gives `C± -> TT_b`, hence25. Deleting a tail vertex
gives the same family at a shorter binary word; its H is at least73.
For a core vertex v, deletion has the form

    d_v g(B)+o_v B+t_v(B mod2)+c_v.

The exact tables, with rows v0 through5, are

|v|C- `(d,o,t,c)`|C+ `(d,o,t,c)`|
|---:|---|---|
|0|`(3,18,2,6)`|`(5,15,3,6)`|
|1|`(11,32,6,18)`|`(13,32,6,20)`|
|2|`(9,26,6,10)`|`(9,41,7,16)`|
|3|`(5,27,5,4)`|`(3,21,3,2)`|
|4|`(13,58,8,24)`|`(11,49,9,12)`|
|5|`(13,32,6,18)`|`(13,32,6,20)`|

Every entry is at least29 at B1. For B>=2, using only g(B)>=1
and discarding parity gives a lower bound at least41. These cover every
possible deletion and prove the theorem.

Consequently these actual graphs cannot be one-birth constructions over
any of the five cores or over an H13 core, under any relabelling or choice
of apex. Such a presentation would exhibit a deletion count at most13.
This is an actual-graph statement; their scalar values may have other,
smaller-parent realizers.

At211/213 the old canonical B4 pair has9 vertices and minimum deletion13.
The new B1 pair has8 vertices and minimum deletion25. Thus shortening
the graph and lowering its parent complexity are different operations.
The actual replacement raises one complexity measure while improving
another, and the next section shows that it also improves every future
O/Z count at a fixed address.

## 5. An all-address cost decrease after a simultaneous same-count reset

The old H13 pair at B4 and the new pair at B1 have the same211/213 totals,
but their common sidecars are respectively

    old: (E,U,d)=(13,347,13),
    new: (E,U,d)=(25,269,25).                         (10)

For any initial rooted state `s=(h,E,U,d)`, the inherited O/Z recurrence
gives the exact scalar response

    F_s(C)=h+d(g(C)-1)+(U-d)C+(d-E)(C mod2).         (11)

This follows by induction on the word, or from the rooted flank compiler.
Applied to the two lower sources in (10), it gives

    F_old(C)=13g(C)+334C+198,
    F_new(C)=25g(C)+244C+186,
    F_old(C)-F_new(C)=90C-12(g(C)-1).                (12)

Both upper sources have exactly2 more paths under every common word.
The inherited elementary bound

    g(C)<=C+1-(C mod2)

therefore gives

    F_old(C)-F_new(C)>=78C+12(C mod2).               (13)

For every `C>0` the new pair has a strictly smaller baseline. This is a
uniform cost comparison of actual continuation programs, not an assumption
that a prescribed target can be placed in either family.

**Theorem2 (general source-comparison rule).** Suppose two actual rooted
sources have the same H and

    d'-d=E'-E=a>=0,      U-U'=b>=0.

Then for every common O/Z address C,

    F_s(C)-F_(s')(C)=(a+b)C-a(g(C)-1)
                      >=bC+a(C mod2).               (14)

Proof: subtract (11). The parity coefficients agree because `d-E=d'-E'`.
Apply the bound on g. QED. Strict decrease follows at C>0 when b>0.
The new concrete reset has a12,b78. Its ready-source defect
`U-H-d+E` drops by78, from136/134 to58/56 for the lower/upper sources.
It does not become1-ready or2-ready at the retained root.

This comparison is the main reusable result. It gives exact inequalities
for replacing a source, including the coordinate that increases and the
coordinate that decreases. Neither H alone nor a shared overlap spectrum
would reveal this order.

## 6. The exact stopping boundary of the native search

At the original31/33 pair, every common retained-root one-birth cut
preserving a pure H gap2 is O, Z, or a dominant neutral vertex. Allowing
different cuts on the two sides yields no additional compatible pair at
that base. At69/71, the asymmetric46/40 and46/42 lifts produce the new
25-source pair. This is why passing from identical cuts to a shared-sidecar
fibre product was useful.

The next cheap probe was pursued. At the equal73/E-different pair (3),
the complete retained-root matching of two128-cut sets finds exactly one
pure-H unit-gap output: common Z, namely211/213. At the new211/213 pair,
the complete matching of two256-cut sets gives only O/Z, neutral-source,
and same-response tail-placement variants. No further deletion-carrier
change appears in that one-birth universe. These are **FINITE-EXACT**
statements with the actual masks in `receipt.json`.

Thus the remaining operation cannot be silently called repeatable.
It may need two births, rerooting, a different source, or a more general
pair of endpoint defects. The other remaining gate is target entrance:
(12)-(14) decrease continuations of this known reset, but do not construct
a witness for an arbitrary missing target. The H13 unit-switch family,
the new25-source family, and a finite union of such fixed one-address
families do not by themselves prove all odd-value coverage.

## 7. Historical recovery and connection contracts

The sources below were read for distinct mechanisms, not as a campaign
to reverify old claims.

|Recovered source|Mechanism used here|Actual target map and its limit|
|---|---|---|
|March16 `07-reflections/supersingularity-triples.md`|Old15-prime and golden/icosahedral compulsion|Led outward to level structures; no number13/25 coincidence used as a theorem|
|September27 `ogg_triangular_chowla_20260927.md`|Changing the retained coordinate turns a nodal-cubic condition into a conic or elliptic condition|Retain the root endpoint rather than only H; no elliptic curve is asserted for the binary image|
|October22 `monster_levels_square_sum_7_21_20261022.md`|Level1/2/6 refinements distinguish objects with the same coarse endpoint|Compare root-response levels of the same211/213 counts; graph/source hierarchy remains explicit|
|October23 `cyclic_attachment_lattices_20261023/MONSTER.md`|Common-C2 fibre product and central-sign projectors|Pair two different birth cuts over exactly the same `(E,U,d)` response; independent cuts are accepted only after their common coordinate matches|
|October31 `session131_reflection_numeration_moonshine_20261031/moonshine/README.md`|Frobenius exchanges parallel level-two edges while fixing both endpoints|The new73 pair agrees in H but differs in E; Z reveals this difference and O erases it|
|July H13 hostile, recovered in the prior five-core history note|Golden triangle-overlap graph with an endpoint-sensitive response|Supplies the exact Q and cuts in (1), rather than another abstract pentad|
|Current S253 integral/affine lattice work|Both H9 charts repair the same affine5-torsion but different targets|We change an actual positive source fibre, with exact cost and deletion data, rather than repeat a linear modular test|

For the modular/Frobenius connection the source is a two-edge lift over
fixed endpoints; the target is a pair of actual rooted attachments over
a common response. The map is the **shared-quotient construction method**,
not a claimed isomorphism between a modular curve and tournament data.
It preserves compatibility of the two chosen lifts. It forgets group,
field and geometric structure; no such structure is used in (2)-(14).
The decisive test was the explicit source-reset search and its literal
counting receipt. The positive result is the source operation itself.

## 8. Reproduction

Run from the repository root:

```bash
python3 -B 05-knowledge/results/h13_paired_descent_20261111/modular_cycles/probe.py
```

The script constructs literal tournaments, computes `(H,E,U,d)` by
subset dynamic programming, derives the finite compiler signatures and
deletion tables, and checks all B0..15 family members. It records every
matching retained-root one-birth pair in the declared small universes and
checks the cost difference for common addresses C0..1023. The all-height
claims follow from the displayed compiler identities and inequalities,
not from enlarging a coverage census. Source hashes and explicit graph
rows for the simultaneous reset are included in `receipt.json`.
