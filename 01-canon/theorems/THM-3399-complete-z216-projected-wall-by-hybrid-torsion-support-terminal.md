---
id: THM-3399
title: "Complete z1=216 projected wall by hybrid torsion and common-source support terminals"
status: >
  PROVED ANALYTIC MECHANISMS + FINITE-EXACT + INDEPENDENTLY VERIFIED-EXACT
  ROW-195 TWO-HIGH REPAIR + SELF-CONTAINED BOUNDED ROW-195 STATE/ONE-HIGH
  ACCELERATOR. The historical post-THM3361 110-row projected wall is
  redundantly rechecked; after THM-3378, all 109 live rows and all 12 live
  families close and the projected k=3 cap moves 216 ->215. This is a
  necessary projected-screen theorem only, with no physical entry,
  arbitrary-k, rung, or LRC(14) conclusion.
source: root/lrc14-projected-wall/2026-08-14
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-2979-projected-k3-z275-ten-body-status-and-located-torsion-closure
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
  - THM-3351-projected-k3-z216-sixteen-family-translated-two-high-closure
  - THM-3361-projected-k3-L720720-one-high-translated-residue-closure
  - THM-3378-projected-k3-z216-gcd24-L129360-row94-one-high-torsion-closure
  - THM-3391-weighted-common-source-cyclic-support-capacity
related:
  - THM-3385-odd-fibre-doubling-projection-and-half-even-complement-clocks
  - THM-3387-exact-cyclic-sheet-cover-atlas-and-q2-gcd-graph
  - THM-3398-general-finite-mode-sheet-cover-cochain
  - THM-3402-atomized-sheet-covers-and-constructive-cochain-locus
script: 04-computation/lrc14_j7_k3_z216_post_thm3378_complete_wall_hybrid_closure_20260814.py
output: 05-knowledge/results/lrc14_j7_k3_z216_post_thm3378_complete_wall_hybrid_closure_20260814.out
script_sha256: 7964db3080e1329dfdfa66df53c451a653d88e7c411430d76c3cd3b8c239b765
output_sha256: 74b268bb96c436c1c659e743adce5ab99270189e8758466ac5340a700b0be920
semantic_sha256: 125200a407d2d055ff907097cabb9894817cae37f468b55e6a79b1c3d1a77920
audit_script: 04-computation/lrc14_j7_k3_z216_row195_two_high_top_weight_majorant_independent_audit_20260814.py
audit_output: 05-knowledge/results/lrc14_j7_k3_z216_row195_two_high_top_weight_majorant_independent_audit_20260814.out
audit_script_sha256: fccc10392624bbdfc2083993ad51a423e8974c135b9bc635351304a71cb0de74
audit_output_sha256: 55ad1da385d35f5b38fdc7de2d9916f54ee874b806496d4a7b1270cb526ad30c
audit_semantic_sha256: 7a70daaa1f25fc6fb9bfd7469ce3e5a18ba650296119789c81bc70fd32016b9c
accelerator_script: 04-computation/lrc14_j7_k3_z216_row195_three_block_top_h_torsion_accelerator_audit_20260815.py
accelerator_output: 05-knowledge/results/lrc14_j7_k3_z216_row195_three_block_top_h_torsion_accelerator_audit_20260815.out
accelerator_script_sha256: 00e1855dd428466c113d3d559e629baefcfdd9b0ea8ac732069fee3bb4ac74fe
accelerator_output_sha256: b188805a09b484466f4c36395b2b426beba8ad470704d9c28414f92ec7a43f13
accelerator_semantic_sha256: 5ad577c0a388c0021cda9506a8147a25732f36f801e7bcf746b856a5d9ff1070
hash_basis: LF-normalized bytes
---

# THM-3399 -- the projected `k=3,z1=216` wall closes

**PROVED ANALYTIC MECHANISMS + FINITE-EXACT, with independently verified-exact
row-`195` two-high repair and a self-contained bounded row-`195`
state/one-high accelerator.**

## 1. Exact statement and live composition

In the necessary projected `k=3,z1=216` atlas inherited from THM-2941,
THM-3361 leaves a historical queue of `110` wall rows in `12` complete
intrinsic-cost families.  THM-3378 has already closed historical row `94`,
so the current live queue is exactly

```text
historical post-THM3361 indices minus {94},
110 historical rows -> 109 current live rows.                    (1)
```

The integrated exact audit separately and redundantly rechecks all `110`
historical rows through the pinned inherited machinery, including row `94`,
and excludes every one.  It therefore closes all `109` currently live rows,
not a fresh set of `110`.  Composed with the proved THM-3378 baseline, the
exact consequence is

```text
projected ledger       372913 -> 372804,
z1=216 wall rows          109 ->      0,
complete families          12 ->      0,
projected k=3 cap         216 ->    215.                         (2)
```

The historical and live ordered-index digests are respectively

```text
75f7e326e7694af41d6eedc3c520d4f7dd3c38e7fa83b34d7eb5bccc5cda460d,
e0d24b7f8a46cfdaa4b3c2fbe3aca070b4004b7b6a81f1f1083799b8e9138d80, (3)
```

and the exact twelve-family packet has digest

```text
3f9bd7b1ff330b6a804fae4a280265b2b5791d7bd33ddeec00f03f1b377d3921. (4)
```

This is a necessary projected-screen theorem only.  It proves no physical
entry, endpoint origin, owner/current gluing, arbitrary-`k` result, rung, or
LRC(14).

## 2. The exact historical wall

The twelve historical post-THM3361 families are:

```text
rank  gcd        L       cost       historical indices
 1     24    129360    52778880     94,161,174,215,237,263,319,354,369,399,430,443,472
 2     72    229320    58705920     108,123,156,283,406,417,420,435,447
 3     72    504504    69621552     143,159,310,436,448
 4     72    194040    84213360     106,120,154,200,240,281,307,341,372,404,415,418,433,445,475
 5     24    152880    96925920     83,96,163,176,198,217,239,265,321,356,371,380,388,393,401,432,444,456,457,474
 6     24   1681680   117717600     392,473
 7     36   1261260   138738600     201,329,342,373
 8     72    388080   142813440     181,192,229,326,333,363,389,411,439,452,462,466,476
 9     72   1009008   175567392     224,235,367,442,470,479
10     72   2522520   206846640     408,419,446
11     72    458640   221064480     183,195,223,232,327,336,365,384,391,407,413,441,454,461,464,468,478
12     72   5045040   433873440     459,463,477.                    (5)
```

Only rank `1` loses an index when passing from the historical to the live
queue: its live indices are the twelve entries in `(5)` other than `94`.
Every later family is unchanged.

The exact ray/common-status screen over the historical wall partitions

```text
596799 denominator states
  =279934 crude capacity exclusions
   +298387 exact common-status exclusions
   + 18478 residual passports.                                  (6)
```

The per-family `(states,crude,status,residual)` totals are

```text
 1  ( 27892, 13774, 13553,  565)    7  (  1288,   451,   688,  149)
 2  (  6468,  3479,  2672,  317)    8  (110261, 49815, 56348, 4098)
 3  (  4334,  2141,  1847,  346)    9  ( 20198,  8425, 10429, 1344)
 4  ( 19246,  9678,  9031,  537)   10  ( 19051,  7311, 11090,  650)
 5  ( 47122, 24077, 22057,  988)   11  (291553,137923,144782, 8848)
 6  (  4234,  2771,  1324,  139)   12  ( 45152, 20089, 24566,  497).       (7)
```

A separate scalar-only census, using the same pinned prefix and status engine
through a separate execution path, classifies the `110` rows as

```text
43 screen-empty rows,
66 residual rows with strictly positive duplicate-two-high gap,
 1 residual row with nonpositive gap: row 195.                   (8)
```

The integrated terminal recomputes this classification, and its complete
packet digest is frozen at runtime as

```text
c39709a32e08e8c24b7d25259a85a85c2816d09df250a7fdab33feed1943dd89.                                      (9)
```

The weakest strictly positive gap is row `183`, with exact value

```text
16964035/350997016188 > 0.                                    (10)
```

The three heavy-family rows are also positive:

```text
row 459: 522773597341321/196871471610763950,
row 463:   5876658008161/2683097971003450,
row 477:   1185688570401/759552141808460.                       (11)
```

Thus row `195` is the unique obstruction to the tempting universal
positive-gap proof.

## 3. High gate, positive-gap rows, and one-high terminals

THM-2941's strict high gate gives at least one high later drift on every
residual row because the fixed first drift is `216` below the row's high
floor.  Scalar zero-high records are retained as hostile accounting but are
excluded only by this inherited gate.

For each of the `66` positive-gap residual rows, the exact
duplicate-permitting two-high upper bound lies strictly below the necessary
scalar threshold.  The high gate and the gap therefore force exactly one
high drift.  For each resulting one-high case, let `d` be its exact high
denominator and let `C` be the common set of actual complete cells fixed-safe
for the first drift and all low labels.

Two lawful terminals cover all such cases:

1. **Located torsion.**  Two distinct cells of `C` have nonzero residue
   difference of exact effective order `r` with `2<=r<=7` modulo `d`.  Unit
   multiplication preserves exact order, height cancels in the cell
   difference, and a translated strict `d/7` band cannot contain both
   residues.  One common cell survives at every local coordinate.
2. **Translated support cardinality.**  If the located-pair selector does
   not fire, the exact distinct support `C mod d` has size strictly greater
   than `ceil(d/7)`.  Every unit and translated strict band contains at most
   `ceil(d/7)` residue classes, so one residue, hence one actual cell,
   survives.  The known boundary case occurs in the row-`181` terminal.

Row `195` is not expanded through this generic path.  Before mask enumeration,
the bounded accelerator forms the actual common source `C_216` fixed-safe for
the first label only:

```text
|C_216|=119368, first=37800, last=420839,
u32 SHA-256=745925d13ce82fd6f85d5184de86a12a3cb2c370890d6f60ec3d95d354d13cae. (12a)
```

For each of the `78` attained suffix denominators it computes the top-`h`
atom majorant from Section 4.  Over all three-denominator multisets and then
the inherited state screen, the exact typed partitions are

```text
82160 suffix triples = 80687 strict top-h + 1473 failed suffix triples,
695 attained divisor states among those failures = 295 crude + 243 status
                                             + 157 screened residual states,
1565 row-195 residual states = 1408 strict top-h + 157 hostile.          (12b)
```

The weakest strict state margin is `4`, at passport
`(1872,6370,57330,458640)`.  The closest hostile is exact equality, with
suffix triple `(10192,65520,458640)` and majorants
`(22696,31152,65520)`.  On only the `157` hostile states, the pinned
state-local one-high routine yields `599` cases on `156` passports.  All
`599` close by located torsion; no translated-cardinality branch is needed
there.  The qualifying/effective order censuses and weakest witness are

```text
qualifying orders ((2,467),(3,89),(4,7),(6,36)),
effective orders  ((2,469),(3,124),(4,6)),
minimum surplus 1 at passport (2,5733,6370,152880), high d=2,
low labels (237,320), cells 37800 and 37801.                 (12c)
```

The hostile bank also has `109` zero-high records excluded by the strict high
gate.  Its unique passport with no one-high case is `(2,2,3920,6370)`; that
passport has the single denominator-two two-high case routed through the
full two-high audit of Section 5.  The accelerator never constructs the
historical fully expanded one-high tuple, and `1408` states must never be
added to `599` scalar cases.

The final exact mechanism counts, weakest witnesses, and terminal digest are

```text
(('located-torsion', 113727), ('translated-cardinality', 1)),
(1, 94, (2, 2352, 5390, 6468), 2, (260, 275), (2, 2, 2, 1, 9912, 9913, 0, 1, 1, 1, 24414)),
(9, 181, (11, 5390, 7056, 17640), 11, (275, 286), 11, 2),
c39709a32e08e8c24b7d25259a85a85c2816d09df250a7fdab33feed1943dd89.                                      (12d)
```

Every non-heavy compact common-cell array and located-torsion witness other
than row `195` is checked for exact equality with the inherited tuple
implementation.  On row `195`, the integrated verifier instead enumerates the
compact arrays literally over the carrier ranges with the direct cell-clean
predicate and directly rechecks every selected torsion cell.  The exact
accelerator packet has SHA-256

```text
8dcede7f9ee5b92eb9cc45c5b101e5bbd6f0862ca85942d8cc9ff12e73f94b2d. (12e)
```

On the three heavy rows, the compact replay reproduces the frozen inherited
terminal packet digest

```text
9bbac6d9e329647d830097de5a02be273c85b785039c77bd96ee7c4a54a45ed9. (13)
```

The compact representation uses a four-byte unsigned cell array and one
`ceil(d/8)`-byte support bitset at a time; selected witness cells are then
recovered and directly rechecked against the literal cell-clean inequalities.

## 4. Common-source top-`h` atom majorant

The repair of row `195` uses the following named specialization of THM-3391.

Let `C` be a finite common source of total nonnegative weight

```text
W=sum_(c in C) w(c)>0,                                      (14)
```

and let blocker `i` map the same source by `phi_i:C -> Z/d_i Z`.  For every
unit `a` and cyclic translate `t`, define the pullback danger edge

```text
E_i(a,t)={c in C : a*phi_i(c) lies in the translated strict d_i/7 band t}.
                                                               (15a)
```

Aggregate the fibre weights

```text
w_i(r)=sum_(c: phi_i(c)=r) w(c),
h_i=ceil(d_i/7).                                            (15)
```

Let `K_i` be the sum of the `h_i` largest values of `w_i(r)`, padding with
zero fibres if necessary.  Multiplication by a unit permutes residue classes,
and every translated strict-open interval of length `d_i/7` meets at most
`h_i` integer residue classes.  Hence every lawful unit, height, and local
coordinate translate removes weight at most `K_i`.  Therefore

```text
W-sum_i K_i>0                                               (16)
```

leaves one positive-weight element of the **same source `C`** outside every
blocker.  This is the common-source top-`h` atom majorant.

The strict inequality in `(16)` is a sufficient top-`h` terminal.  At its
equality boundary `W=sum_i K_i`, the coarse majorant alone gives no survivor,
but it also does not prove a cover.  Let `lambda_i` be the exact THM-3391
pullback-band capacity of blocker `i`.  Since `lambda_i<=K_i`, any strict
slack `lambda_i<K_i` gives `sum_i lambda_i<W` and the exact common-source
theorem still leaves a survivor.  Only when `lambda_i=K_i` for every `i`,
equivalently `sum_i lambda_i=W` on this top-`h` boundary, does THM-3391's
independent-phase equality criterion apply: a cover exists exactly when
maximum-weight pullback-band edges can be chosen whose positive-weight parts
partition `C`.  With correlated phases, one must then retain and inspect the
exact cover locus.  In particular, top-`h` equality must not be converted
into either an unproved cover or an unproved survivor claim.

There are two useful specializations:

1. On a common modulus `M=lcm(d_1,...,d_m)`, give each represented support
   residue weight one.  A blocker of denominator `d_i` has the crude support
   capacity `ceil(d_i/7)(M/d_i)`.  Support larger than the sum of these
   capacities gives a common residue and therefore a common source cell.
2. On the actual cell source, retain the full multiplicity of every fibre in
   each quotient and compute the `K_i` of `(15)` directly.  This retains
   information destroyed by independently optimized quotient supports and
   is often strictly stronger.

No common `lcm` is mathematically necessary when the actual source and all
quotient maps survive.  The common-modulus lift is merely a cheap sufficient
representation for cases where it closes.

## 5. The unique nonpositive row

The unique row from `(8)` is

```text
row 195: E=(1,5,8,9,13,14), L=458640, high floor=45170,
passport (3920,6370,32760,57330), low label 351,
duplicate-two-high gap=-17651221657/66699959142726.          (17)
```

This is a genuine refutation of the simplified universal-gap route, recorded
as MISTAKE-387.  Exact scalar enumeration leaves `137` two-high cases, no
three-high case, with ordered case digest

```text
5811cb585e35831a9c58b2d0af35ab6db89e4503229e6420259d6ce6d4586c3c. (18)
```

All `137` cases use the same actual common fixed-safe source for labels
`216,351`:

```text
|C|=100776, first cell=37987, last cell=420652,
u32 cell digest=89236b7d22a6be06afecaaaba8f2e0be1346b1cab3c0b720314baa4014d87c36. (19)
```

The exact terminal partition is

```text
76 cases: common-modulus support capacity,
60 cases: actual-cell top-h atom majorant,
 1 case : denominator-two exact projected-measure terminal. (20)
```

For the sixty actual-cell cases, the weakest `(16)` margin is case index
`130`, high denominators `(32760,458640)`:

```text
W=100776, K_1=24800, K_2=65520,
W-K_1-K_2=10456>0.                                        (21)
```

The remaining case has `(d_1,d_2)=(2,2)`.  It lies on the equality boundary
of pointwise capacity and is **not** claimed as a weighted-cardinality
closure.  On the normalized local coordinate `y in [0,1]`, the inherited
exact measure argument writes each high as `z_i=L/2+h_iL` with `h_i>=0`;
one high danger comb removes at most `2/7` of local-coordinate mass from one
fixed-safe cell.  Two remove at most `4/7`, leaving

```text
3/7=39/91>36/91, with exact surplus 3/91.                  (22)
```

Thus every row-`195` two-high scalar packet closes, while the equality case
is kept correctly typed as a measure terminal.

## 6. Exact quotient-loss and positive structural controls

The low labels `216,351` are intrinsically periodic with

```text
P=lcm(3920,6370)=50960=L/9.                               (23)
```

It is tempting to infer that every fixed-clean cell supplies a complete
nine-cell orbit.  That inference is false because the body-safe carrier
`stream.ranges` is not `P`-invariant.  On the actual source `C`, fibre
multiplicities modulo `P` have histogram

```text
((1,760),(2,7944),(3,9848),(4,11626),(5,1616)), maximum 5<9. (24)
```

For source `37987`, the putative orbit has membership mask

```text
(T,F,F,T,T,F,F,F,F).                                      (25)
```

The first missing target `88947=37987+P` remains intrinsically clean for
both `216` and `351`, but lies outside every body-safe carrier range.  This
is the exact missing-coordinate sidecar: quotient periodicity preserves the
clean predicate but destroys carrier membership.

As a positive, non-load-bearing control, case `124` with high denominators
`32760,57330` contains three actual common cells

```text
(198900,46020,122460).                                     (26)
```

Their residues form complete order-three cosets in both quotients, with a
permutation of positions:

```text
mod 32760: (2340,13260,24180)=2340+(0,1,2)10920,
mod 57330: (26910,46020,7800)=7800+(1,2,0)19110.           (27)
```

Every unit permutes the three positions, and a strict `d/7` band is shorter
than the `d/3` separation, so each blocker removes at most one cell and two
blockers leave at least one.  The independent audit exhausts all `6912` and
`12096` units in the two quotients.  This illustrates the located-coset
mechanism related to THM-3398, but no THM-3398 cochain-equality statement is
used in the proof.  THM-3402's atomized constructive cover locus is likewise
an adjacent equality-boundary sidecar, not a load-bearing step.  THM-3398 and
THM-3402 are related, not dependencies.

## 7. Independent two-high audit, bounded one-high accelerator, and runtime gate

The canonical independent row-`195` **two-high** audit does not import a
scratch packet, the scalar scanner, or the integrated verifier.  Starting from
pinned canonical atlas and engine artifacts, it independently reconstructs:

1. the exact atlas row and its `65728=27198+36965+1565` screen;
2. the gap `(17)`, all `137` cases, and the no-three-high result;
3. the actual source `(19)` directly from the body-safe ranges and literal
   cell-clean formula;
4. all common-modulus and top-`h` capacities in `(20)`;
5. the weakest margin `(21)`, the equality typing `(22)`, the carrier hostile
   `(23)--(25)`, and the order-three control `(26)--(27)`.

That audit does not claim an independent replay of row `195`'s one-high cases.
It is load-bearing only for the full `137`-case two-high/no-three-high repair.

Its normal and optimized outputs are byte-identical, with hashes

```text
source   fccc10392624bbdfc2083993ad51a423e8974c135b9bc635351304a71cb0de74,
output   55ad1da385d35f5b38fdc7de2d9916f54ee874b806496d4a7b1270cb526ad30c,
semantic 7a70daaa1f25fc6fb9bfd7469ce3e5a18ba650296119789c81bc70fd32016b9c. (28)
```

Its byte identity is reproduced by

```text
python3 04-computation/lrc14_j7_k3_z216_row195_two_high_top_weight_majorant_independent_audit_20260814.py
python3 -O 04-computation/lrc14_j7_k3_z216_row195_two_high_top_weight_majorant_independent_audit_20260814.py
```

The separate bounded one-high accelerator imports neither a scratch packet
nor the integrated verifier.  It reconstructs the row, screen, exact
`1408/157` state partition, the `599` hostile-state one-high cases, and every
literal compact-cell torsion witness from pinned canonical sources.  It is
therefore independent of the integrated proof artifact and its scratch history,
while deliberately sharing the pinned canonical atlas/status engine rather
than claiming an unrelated reimplementation of that engine.  It bridges the
unique hostile denominator-two case to the independent two-high audit above.
Its ordinary and optimized outputs are byte-identical, with hashes

```text
source   00e1855dd428466c113d3d559e629baefcfdd9b0ea8ac732069fee3bb4ac74fe,
output   b188805a09b484466f4c36395b2b426beba8ad470704d9c28414f92ec7a43f13,
semantic 5ad577c0a388c0021cda9506a8147a25732f36f801e7bcf746b856a5d9ff1070. (28a)
```

Its byte identity is reproduced by

```text
python3 04-computation/lrc14_j7_k3_z216_row195_three_block_top_h_torsion_accelerator_audit_20260815.py
python3 -O 04-computation/lrc14_j7_k3_z216_row195_three_block_top_h_torsion_accelerator_audit_20260815.py
```

The complete row is finally frozen by the integrated packet: its field `35`
is typed only at row `195`, contains the state and case counts separately,
and has digest `(12e)`.  Every other terminal record retains `None` in that
field and every earlier field index is unchanged.

The integrated verifier pins twenty-six direct dependencies, contains no
`assert`-dependent correctness and no floating-point literal, uses explicit
big-endian incremental integers for compact-cell digests, and checks every
selected witness cell by a second literal formula.

The byte-matched final replays freeze the following exact runtime fields:

```text
screen packet SHA-256        748f1a4f9590d0eb72da8a03c3a14d5278a63f26eaba4c99e091f9057b61a0e0,
terminal summary             (67, 18478, 16236, 113728, 17069, 4690, (('common-modulus-support', 76), ('denominator-two-measure', 1), ('located-torsion', 113727), ('top-h-weighted-common-source', 60), ('translated-cardinality', 1)), 113728, 1, 112514, 4232424, 630630, '9bbac6d9e329647d830097de5a02be273c85b785039c77bd96ee7c4a54a45ed9', 137, (('common-modulus-support', 76), ('denominator-two-measure', 1), ('top-h-weighted-common-source', 60)), (10456, 130, (32760, 458640), 24800, 65520), 403104, 1834560, ((195, '5811cb585e35831a9c58b2d0af35ab6db89e4503229e6420259d6ce6d4586c3c', '4cdb5acfda8b7db6e061b54ea2d71102c29a8694a7fce85f6a2c25e317a93ef9', ((351, 100776, 37987, 420652, '89236b7d22a6be06afecaaaba8f2e0be1346b1cab3c0b720314baa4014d87c36'),), ((50960, ((1, 760), (2, 7944), (3, 9848), (4, 11626), (5, 1616)), 5, 37987, 88947, (True, False, False, True, True, False, False, False, False)),)),), ((195, ((119368, 37800, 420839, '745925d13ce82fd6f85d5184de86a12a3cb2c370890d6f60ec3d95d354d13cae'), (78, '6b12fcbf3b4d7415cc4006d72ed57099bd635c09a084dda231e9643dad2c7803'), 'ce68e2ad70faf1bd0d5581ac490d84e345eea97479a42969927431c15bcc5985', (82160, 80687, 1473), (1473, 695, 295, 243, 157), (1565, 1408, 157), '3e3fe58c04ccbbbc1751d48bd199b76782c37589c58625f23767f5f870b8818d', '3ddbcf69c2c1521e1440973cdc46be0e85ca769a6ab68a247c8d00115340e372', 'e997f387c89563b9d87fca78e049f4484f0712a21879a7f4b6088461f704e947', (4, (1872, 6370, 57330, 458640), (1872, 57330, 458640), (18182, 35662, 65520)), ((10192, 65520, 458640), 0, (22696, 31152, 65520)), (599, 156, '935127e640015f44655de13600eb68fe282d7c77fb7c398341d2ef5da4dcbe7e', 'f5d7e673e8394bf35147c3ce8f08e3198fefdc8da9f1295c78af15828a0f7cbe', 109, 'a1111e7538da5d8bebbc938d28db591396fb30f4504c20fe2b4a3023ccbad092', '001859942059725805de93d67634758a40734d8c7e55f3507ab516ddef8b94fa'), (599, 156, 179, 599, ((2, 467), (3, 89), (4, 7), (6, 36)), ((2, 469), (3, 124), (4, 6)), (1, (2, 5733, 6370, 152880), 2, (237, 320), (2, 2, 2, 1, 37800, 37801, 0, 1, 1, 1, 86268)), 352552, 57330), 'a67160cf481123639297f2e121776f6c17e552dbf04dee711bc9eb402a0201ac', (1, '0153f587817922eab7f9012644fe168b000233b3bded211229d3abb934eed194', 0, ('denominator-two-measure', Fraction(3, 7), Fraction(3, 91), '7a70daaa1f25fc6fb9bfd7469ce3e5a18ba650296119789c81bc70fd32016b9c'), (((2, 2, 3920, 6370), (2, 2), 351, Fraction(6457, 31302180)),), (('denominator-two-measure', Fraction(3, 7), Fraction(3, 91)),), (('selfcontained_independent_source', 'fccc10392624bbdfc2083993ad51a423e8974c135b9bc635351304a71cb0de74'), ('selfcontained_independent_output', '55ad1da385d35f5b38fdc7de2d9916f54ee874b806496d4a7b1270cb526ad30c'), ('selfcontained_independent_semantic', '7a70daaa1f25fc6fb9bfd7469ce3e5a18ba650296119789c81bc70fd32016b9c'))), ('2a8584886131c177149055d208f5ccd3b356cc42dcdc84b6533773f6b800b935', '7d8644411a2d78ed40a6c593a7f977979803aac2e0f2427f788ee58da7c9489c', ((1, 65014), (2, 711), (3, 3))), (('selfcontained_accelerator_source', '00e1855dd428466c113d3d559e629baefcfdd9b0ea8ac732069fee3bb4ac74fe'), ('selfcontained_accelerator_output', 'b188805a09b484466f4c36395b2b426beba8ad470704d9c28414f92ec7a43f13'), ('selfcontained_accelerator_semantic', '5ad577c0a388c0021cda9506a8147a25732f36f801e7bcf746b856a5d9ff1070')), ('historical_nonloadbearing', 'a5cdcedc61c07714cacb2164e8793c655fe57ee80fd17ca82149b6710e254036', 314615), 'a353502ccd9bcee61c6981948bbbba9af8a353d17dfc60e25f41c95efdbb6f38')),)),
terminal packet SHA-256      c39709a32e08e8c24b7d25259a85a85c2816d09df250a7fdab33feed1943dd89,
one-high mechanism counts    (('located-torsion', 113727), ('translated-cardinality', 1)),
weakest located torsion      (1, 94, (2, 2352, 5390, 6468), 2, (260, 275), (2, 2, 2, 1, 9912, 9913, 0, 1, 1, 1, 24414)),
weakest translated support   (9, 181, (11, 5390, 7056, 17640), 11, (275, 286), 11, 2),
semantic SHA-256             125200a407d2d055ff907097cabb9894817cae37f468b55e6a79b1c3d1a77920,
final source SHA-256         7964db3080e1329dfdfa66df53c451a653d88e7c411430d76c3cd3b8c239b765,
final output SHA-256         74b268bb96c436c1c659e743adce5ab99270189e8758466ac5340a700b0be920.       (29)
```

The final frozen source passes both commands with byte-identical stdout.  The
two-worker choice changes only scheduling: ordered worker results
are canonically reassembled before any packet or digest is formed.

```text
python3 04-computation/lrc14_j7_k3_z216_post_thm3378_complete_wall_hybrid_closure_20260814.py --processes 2
python3 -O 04-computation/lrc14_j7_k3_z216_post_thm3378_complete_wall_hybrid_closure_20260814.py --processes 2
```

With `(29)` frozen and both replays passed, all historical rows are excluded,
the live queue is empty by `(1)`, and the current-baseline consequence `(2)`
follows.  The first surviving projected wall row at `z1=216` is `none`.

**QED.**
