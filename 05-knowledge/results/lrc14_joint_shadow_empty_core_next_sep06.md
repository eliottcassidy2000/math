# Joint residue shadows sharpen the hereditary gcd ceilings

Status: **PROVED RELATIVE TO CITED LOWER-RUNNER LRC + FINITE-EXACT +
INDEPENDENTLY AUDITED.** The analytic exclusions, complete compiler,
realization lemma, and actual boundary rows passed independent review.
LRC(14) remains **OPEN**.
No theorem ID or external priority is claimed.

Every primitive thirteen-speed strict counterexample `V` to LRC(14)
must satisfy the following necessary conditions:

| Size of P subset V | 12 | 11 | 10 | 9 | 8 | 7 |
|---|---:|---:|---:|---:|---:|---:|
| Upper bound on gcd(P) | 1 | 2 | 4 | 9 | **30** | **90** |

For a nonprimitive row multiply every bound by `gcd(V)`. The earlier
`32,96` bounds remain true, but these are stronger. The two new arguments
retain three shadows simultaneously, where independently minimized pair
intersections lost information. Explicit actual body-safe phases at clocks
30 and90 have every lift blocked, even though their full rows are safe at
another phase. Thus a body-phase-uniform lifting argument cannot simply
delete these two largest surviving clocks.

## 1. Inheritance, scope, and the finite input

The closest mechanism is the
[hereditary CRT/Hunter sieve](lrc14_recursive_gcd_empty_core_next_sep06.md),
with its [independent complete compiler](gcd_pair_hunter_audit_empty_core_next_sep06.md).
The exact object was already identified in
[THM-3387, exact sheet-cover atlas](../../01-canon/theorems/THM-3387-exact-cyclic-sheet-cover-atlas-and-q2-gcd-graph.md):
it is the union of bad sets on one labelled fibre. The canonical hostile
is the [clock9 body-safe phase with no surviving lift](gcd_nine_audit_empty_core_next_sep06.md).
The corrected near miss is treating optimal pair intersections as
simultaneously attainable. The least-used sidecar is a small common quotient
on which several masks have distinct-residue shadows.

For a subset with gcd c and d complementary tails w_i, write

```text
g_i=gcd(c,w_i), q_i=c/g_i, k(q)=ceil(q/7).
```

At any body-safe quotient phase, its bad set is contained in a padded mask
on `Z/cZ` pulled back from an affine unit progression of length k(q_i)
on `Z/q_iZ`. Its size is `g_i k(q_i)`. All masks retain the same sheet
labels. This containment, including strict-danger endpoint cases, is
proved in the inherited note and
[THM-4447, composite-clock capacity](../../01-canon/theorems/THM-4447-lrc14-composite-clock-capacity-and-small-clock-reduction.md).

The complete inherited sieve leaves just these maximum-clock rows:

```text
d=5, c=32: g=(1,1,2,4,4),      next largest surviving c=30;
d=6, c=96: g=(1,1,4,4,6,12) or (1,3,4,4,6,12),
                               next largest surviving c=90.
```

These are complete arithmetic candidate lists, not selected speed-height
samples. The bounded universe follows from the inherited pivot inequality
`q<=floor(6d/(7-d))`, the preceding gcd set, and all-divisor absorption.
The cited lower-runner LRC supplies a strictly body-safe phase for every
proper body. An incomplete padded cover then supplies an actual safe lift.

The live concepts are hereditary gcds; labelled mask unions; common
quotient shadows; sharp root-product margins; coefficient-preserving
operations; and actual-versus-relaxed realization. The new projection
preserves forced intersections, loses unit steps outside the quotient,
and needs a CRT multiplicity and a valid forest as its sidecars. The
cheapest decisive test is the three-mask shadow at clock32.

## 2. The clock32 exclusion

The effective orders are `(32,32,16,8,8)`, so padded sizes are
`(5,5,6,8,8)`, totaling32. A full cover would therefore be a partition.
Each order8 mask occupies two residue classes modulo8. The two such masks
must be disjoint, so together they occupy four distinct classes. Each
order32 mask is an affine unit progression of length5; its five residues
modulo8 are distinct. It must meet the four occupied classes. Since the
order8 masks contain every lift of those classes, this is a nonempty
intersection, contradicting the proposed partition.

Equivalently, the union of the two order8 masks and one order32 mask is
at most20; adding the other sizes5 and6 gives a full union of at most31.
The [independent affine-block audit](clock32_audit_empty_core_next_sep06.md)
checks this three-mask inequality exhaustively. Every pair type separately
has minimum intersection zero. A pairwise minimum table alone cannot see
the contradiction.

## 3. The two clock96 exclusions

The full proof and literal controls are in
[the independent clock96 note](clock96_masks_empty_core_next_sep06.md).
Here is a self-contained tree certificate. Name the order8 mask E,
the order16 mask D, and the two order24 masks A,B. The other masks are
X of order96 and Y of order96 or32. The total padded cardinality is
T=102 or103 respectively.

The inherited balanced CRT formula always gives

```text
|X intersect E|>=2, |X intersect D|>=1.
Y order96: |Y intersect E|>=2.
Y order32: |Y intersect A|>=1, |Y intersect B|>=1.
```

The E,D,A,B shadows modulo8 have cardinalities2,3,4,4. All are
distinct-residue shadows because their progression lengths are below8.
If E meets D, their intersection has at least `96/lcm(8,16)=6` points.
If E meets A orB, it has at least `96/lcm(8,24)=4` points. If E meets
none of them, D,A,B lie in the remaining six classes. Consequently D
shares a class with each of A,B; each shared class contributes
`96/lcm(16,24)=2` actual intersection points by CRT.

The following five-edge trees use these simultaneous facts. A zero-credit
edge is still an edge of the tree and makes no unproved intersection claim.
The `EA` case includes `EB` by exchanging A,B.

| Case | Y order96: edges / credit | Y order32: edges / credit |
|---|---|---|
| ED nonempty | ED6, XE2, YE2, DA0, DB0 /10 | ED6, XE2, YA1, YB1, DA0 /10 |
| EA nonempty | EA4, XE2, YE2, XD1, DB0 /9 | EA4, XE2, YA1, YB1, XD1 /9 |
| E disjoint from D,A,B | DA2, DB2, XE2, YE2, XD1 /9 | DA2, DB2, XE2, XD1, YA1 /8 |

For sets indexed by any tree, the union is at most the total sizes minus
the sum of its edge intersections. Pointwise a tree induced on m occupied
vertices has at most m-1 edges when m>=1. The displayed trees therefore
give bounds at most93 for Y of order96 and at most95 for Y of order32.
Neither signature covers96 labels. This is an analytic simultaneous-mask
obstruction; no exhaustive affine search or common-phase approximation
enters the proof.

## 4. Recompiled hierarchy and its exact limits

Re-running the finite hereditary construction with clock32 excluded at
five tails removes32 from that level's allowable child gcds. At six tails
this also deletes `(32;(1,1,2,4,4,32))` and
`(64;(1,1,2,4,4,32))`; clocks32 and64 themselves still survive with other
signatures. The two clock96 exclusions remove two more rows.

The final profile counts at d=1,...,6 are
`1,2,5,19,109,1213`, including the gcd-one signature at each level.
The eight-body gcd set is

```text
{1,2,3,4,5,6,8,9,10,12,15,16,18,24,30}.
```

The seven-body set is the inherited43-element set with96 deleted; its
largest element is90. The exact full signature lists are retained in
[the JSON](lrc14_joint_shadow_empty_core_next_sep06.json).
This compiler imposes only the inherited arithmetic filters and the
three new clock-specific exclusions. It does not classify every possible
simultaneous cover at all the other retained signatures. In particular
1213 is the size of this declared relaxation, not a realizability count.

## 5. Actual body-phase obstructions at the new maxima

Two simple padded partitions explain the boundary. At c30 take order15
blocks `{2,3,4}`, `{7,8,9}`, `{12,13,14}` and order10 blocks
`{0,1}`, `{5,6}`. At c90 take order45 blocks `8..14`, `23..29`,
`38..44`, order30 blocks `3..7`, `18..22`, and order15 block `0..2`.
All masks are pulled back to their respective full clocks. They partition
the labels, with primitive gcd signatures `(2,2,2,3,3)` and
`(2,2,2,3,3,6)`.

More strongly, the following are **literal physical realizations**, with
`V=cC union T`. The tails occur in the same order as the blocks above.

| c | C | body phase y | T |
|---:|---|---|---|
| 30 | `{1,...,8}` | `112/1009` | `(25082,24992,85712,11073,123)` |
| 90 | `{1,...,7}` | `126/1009` | `(542,55082,25292,211773,30513,51126)` |

The body clearances are respectively at least112/1009 and126/1009,
strictly above1/14. The exact strict bad sets

```text
E_w(y)={j in Z/cZ: ||w(y+j)/c||<1/14}
```

are precisely the displayed padded partitions. Thus every lift of the
specified body-safe y is spoiled. Each full row has thirteen distinct
positive speeds, is primitive, and nevertheless is weak-safe at
`x=1/(14c)`. Integer inequalities in the source verify every claim; no
height or denominator minimality is asserted. An exhaustive subset check
gives maxima `(1,2,3,3,30,30)` and `(1,2,3,3,6,90)` at sizes12 through7
for these two rows. Both also pass every stated global subset gcd ceiling.

The analytic realization mechanism is broader. Fix any finite collection
of effective orders q_i dividing c, any affine unit blocks of their maximum
lengths k(q_i), and any nonempty open interval I inside `(0,1)` of body-safe phases.
One can **choose distinct positive tails** and a common rational y in I
whose actual strict bad masks are exactly those prescribed blocks.

To prove this, choose a prime L sufficiently large that `L>7max q_i`
and some `y=p/L`, `0<p<L`, lies in I. Choose a positive unit a_i modulo
q_i inverse to the desired progression step and consider tail

```text
w_i=(c/q_i)(a_i+q_i N_i).
```

On sheet j its phase is `a_i j/q_i+a_i y/q_i+N_i p/L` modulo1.
The offsets giving any prescribed maximum consecutive image block
contain an open interval of length

```text
1/7-(k(q_i)-1)/q_i >= 1/(7q_i).
```

The L choices of N_i modulo L traverse a translated uniform L-grid, so
one lies in that interval. Adding multiples of L to N_i keeps every mask
unchanged while making the tails positive and pairwise distinct and
avoiding the fixed finite body. Their gcds with c remain c/q_i. This also
covers the strict endpoint convention when q_i is divisible by7.

The quantifier matters: the tails are chosen to realize the masks. For
already fixed tails, independently specified phases need not glue. That
obstruction and its complete cochain are retained in
[THM-3402, atomized sheet covers](../../01-canon/theorems/THM-3402-atomized-sheet-covers-and-constructive-cochain-locus.md).
This lemma supplies a hostile for a phase-uniform statement over arbitrary
tails; it does not discard the actual common phase of a given row.

### Seven tails: an actual obstruction at every clock

For every integer `c>=2`, there is a primitive thirteen-speed row
`V=c{1,...,6} union T` whose seven tails are individually coprime to c,
with a strictly body-safe phase having no safe lift, while the same full
row is safe at `x=1/(14c)`. This is an **analytic all-clock corollary**,
independently proof-audited by `certificate_audit`.

Put `k=ceil(c/7)` and prescribe seven unit-step blocks starting at `ik`,
`i=0,...,6`, modulo c. Their unreduced union contains `0,...,7k-1`, so
they cover every residue. The body `{1,...,6}` has a strictly safe open
interval around `1/7`. Choose a prime `L>7c`, a rational y=p/L there,
and residues N_i modulo L as in the realization lemma with a_i=1.
Impose in addition `N_i=1 modulo14` by CRT. Then the tails
`w_i=1+cN_i` retain their prescribed bad masks and satisfy

```text
w_i/(14c) = (c+1)/(14c) modulo1,
1/14 < (c+1)/(14c) <= 3/28 < 13/14.
```

Adding sufficiently large distinct multiples of14L to the N_i preserves
both requirements and separates positive tails. The tails cannot coincide
with body speeds because they are one modulo c. Their gcds with c are
one, so the full row is primitive; the six-body gcd is exactly c. The
body itself is weak-safe at the displayed full-row witness.

Thus the previous seven-tail arithmetic nonfiniteness has an actual
body-phase realization at every clock. It does not show failure of LRC,
or exclude a bound using other phases or additional global information.
The proof is the grid/CRT construction, not an extrapolation from a finite
clock bank; the c30/c90 physical rows above are separate finite controls.

## 6. Verification and next object

The standalone [source](../../04-computation/lrc14_joint_shadow_empty_core_next_sep06.py)
uses integer arithmetic, the analytic inherited finite universe, and a
Prim tree pass. It imports no repository mathematical producer. Its46
explicit gates check the complete excluded-signature lists, hierarchy
counts, relaxed partitions, and the two exact actual rows. The independent
clock32 and96 controls supply the separate geometric audits.

```text
python3 -B 04-computation/lrc14_joint_shadow_empty_core_next_sep06.py --write-profiles /tmp/joint-shadow.json
python3 -B -O 04-computation/lrc14_joint_shadow_empty_core_next_sep06.py
```

Normal and optimized outputs are byte-identical. Frozen raw LF SHA-256:

```text
3f906146953677c5e1734020e97ef82fb801ee66cf5ae7a6c697ce83b8d21245  source
950acfb7073ec93af9372c6ff41a5a010281a5b9f1cf38df5de62b808478fa21  output
935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f  JSON
```

Revisiting the concept board: a common quotient made the pairwise sieve
strictly stronger, while the actual boundary rows show why a further
step needs body-phase selection or inherited subset interactions. The
wildcard's sharp inequality likewise retains integer multiplicities;
the Laurent niche's coefficient lowering retains its actual zero. These
are method connections, not reductions between the theorems. The useful
next LRC object is a retained body component together with its simultaneous
mask word, as in the incoming
[virtual pair-wall construction](overnight9_20260906_lrc_virtual_pair_wall.md),
the incoming [quantitative pair packet](lrc14_quantitative_pair_packet_overnight_hexagon_sep05.md),
or a small joint-shadow obstruction in another residual signature. The
pair packet retains actual dyadic separation and bounds safe measure by
`(3-e)mu(A_delta)/6`; its residue family supplies margin1/12 at unbounded
body and tail heights. This is a conditional geometric supplier, not a
universal consequence of the present gcd restrictions.

Final independent audit by `certificate_audit`: **PASS**, recorded in
[the audit sidecar, Section6](clock32_audit_empty_core_next_sep06.md).
Every full profile list agrees with the frozen independent ancestor after
exactly the stated deletions. The referee checked all six clock96 trees,
the common-rational-phase quantifiers, both actual rows and their subset
maxima; normal/optimized output and regenerated JSON agree byte for byte.
