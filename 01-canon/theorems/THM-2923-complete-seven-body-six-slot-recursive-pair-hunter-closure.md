---
id: THM-2923
title: "Complete seven-body/six-slot recursive pair-Hunter closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY AUDITED.  The marked
  Hunter recursion closes the complete 3,432-root rung E in C([14],7), with
  six external speed slots.  Its sole zero-margin terminal is the doubled AP
  endpoint 2*{1,...,13}, already safe by THM-2907.  This repository-grade
  audit is not external peer review or Lean formalization.  This is not
  unrestricted LRC(14): the at-most-six-in-window / j>=7 loose-spread sector
  remains open.
source: codex-all-root-recursive-pair-hunter-2026-07-29
depends_on:
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2897-partition-cap-tropical-convolution-and-alternating-pair-ladder
  - THM-2904-all-root-ranked-h1-hunter-closure
  - THM-2907-paircap-exception-h4-heavy-link-child-census
  - THM-2913-one-h3-row-pair-hunter-toothpick-closure
  - THM-2915-all-open-centre-exact-child-top-four-closure
related:
  - THM-2892-eight-body-five-slot-heavy-triangle-closure
  - THM-2920-two-h3-row-pair-hunter-recursive-toothpick-closure
  - THM-885-covering-case-decomposition-j56-sweep
verification:
  - 04-computation/lrc14_j6_seven_body_six_slot_recursive_pair_hunter_closure_thm2923.py
  - 04-computation/lrc14_j6_seven_body_six_slot_recursive_pair_hunter_thm2923_stage1.py
  - 04-computation/lrc14_j6_seven_body_six_slot_recursive_pair_hunter_thm2923_stage2.py
  - 04-computation/lrc14_j6_seven_body_six_slot_recursive_pair_hunter_thm2923_composition.py
  - 04-computation/lrc14_j6_seven_body_six_slot_recursive_pair_hunter_thm2923_endpoint_audit.py
  - 04-computation/lrc14_j6_seven_body_six_slot_recursive_pair_hunter_thm2923_inherited_audit.py
  - 05-knowledge/results/lrc14_j6_seven_body_six_slot_recursive_pair_hunter_closure_thm2923.out
  - 05-knowledge/results/lrc14_j6_seven_body_six_slot_recursive_pair_hunter_closure_thm2923.ledger.out
---

# THM-2923 -- complete seven-body/six-slot recursive pair-Hunter closure

**PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY AUDITED.**

The pinned THM-2915 dependency and the complete source/artifact graph replay
cleanly.  An independent hostile implementation and further frozen ordinary
and optimized replays agree exactly.  Here “independently audited” means a
repository-internal mathematical, source, and replay audit, not external peer
review or Lean formalization.  Nothing below claims the unrestricted
fourteen-runner conjecture.

## 1. Exact statement and scope

Fix the seven-in-window rung

```text
E in C({1,...,14},7),       six further distinct speeds at least 15.       (1)
```

THM-2904 and THM-2915 represent this rung by `11,842` marked parent rows and
`51,222` ordered first-centre children.  THM-2915 closes `46,356` children by
the exact global top-four test and leaves `4,866`.  Every one of those failures
is discharged by the following exhaustive marked recursion:

```text
four remaining slots
  exact pair partition                                    1,612
  additional G4 Hunter compatibility                      1,175
  ordered second-centre recursion                         2,079
    children closed by all exact top-three pivots          1,884
    children reaching a failed second pivot                  195

three remaining slots: 228 failed second pivots
  exact pair plus leading singleton                          65
  additional G3 Hunter compatibility                        114
  ordered third-centre recursion                             49
    rows closed by all terminal exact pairs                  48
    endpoint equality rows                                    1.            (2)
```

The unique equality is the family

```text
(2,4,6,8,10,12,14,16,18,20,22,24,26)
  = 2*{1,...,13}.                                            (3)
```

It is an endpoint-safe member of THM-2907's exception bank.  Exact route
composition then gives

```text
parent rows covered by E union P union T2                 11,842
route roots                                                  3,411
proved baseline roots through THM-2913                        351
route/baseline overlap                                        330
complete seven-body union                       3,411+351-330=3,432.          (4)
```

Thus the theorem closes the complete rung `(1)`.  It does **not** treat
families with at most six speeds in `{1,...,14}` and at least seven external
speeds.  That loose-spread sector remains open, so this theorem is not
unrestricted LRC(14).

## 2. The general marked Hunter recursion lemma

Let `R` be a literal residual carrier of mass `h`, let `F` be its forbidden
label sidecar, and suppose `k` speed slots remain, so a hypothetical completion
would choose `k` distinct allowed labels from an infinite pool.  For an allowed
label `w`, write

```text
c(w)=mu(R intersect D_w),
q_1>=q_2>=...>=q_k
```

for the first `k` globally sealed singleton ranks over the whole allowed pool,
and let

```text
B=beta_2(R,F)
```

be the exact largest union of two distinct allowed danger sets.  Define

```text
G_k(a)=a+sum_(j=2)^k min(a,q_j,B-a),
                 0<=a<=min(q_1,B).                         (5)
```

For any putative `k`-cover, order all its labels by nonincreasing singleton
coverage, breaking ties by the fixed global label order, and let the first
coverage be `a`.  Then the `j`-th chosen coverage is at most the global rank
`q_j`.  If `b_j` is label `j`'s contribution outside the first danger set,
then

```text
b_j<=a                    by maximality,
b_j<=q_j                  by the global singleton ranks,
b_j<=B-a                  by the exact pair cap.            (6)
```

Consequently its covered mass is at most `G_k(a)`.  Therefore

```text
max_a G_k(a)<h                                               (7)
```

excludes a `k`-cover.  The objective is continuous, concave and piecewise
affine.  Its exact maximum is attained among the clipped breakpoints

```text
0, min(q_1,B), B/2, q_j, B-q_j       (2<=j<=k).             (8)
```

To obtain the first crossing, the verifier finds the first breakpoint interval
whose endpoint values straddle `h` and solves the affine crossing inside that
interval exactly.  If `(7)` fails, put

```text
lambda=min{a:G_k(a)>=h}.                                    (9)
```

Every hypothetical cover then has maximum singleton at least `lambda`.  When
`lambda>h/7`, THM-735's discrepancy tail makes

```text
H_1(R,F,lambda)={w notin F:c(w)>=lambda}                    (10)
```

finite.  Order this core by decreasing coverage and then increasing label.
Branch on the earliest label among all maximum-coverage labels, subtract its
literal danger set, and enlarge the sidecar by that centre and every earlier
core label:

```text
(R,F,k) -> (R minus D_y, F union {y} union earlier(y), k-1). (11)
```

The earlier-label set is proof-bearing: it makes the branches exhaustive and
disjoint without reusing a maximum already allocated to an earlier branch.
Equation `(11)`, rather than an unmarked residual, is the recursive object.

## 3. Exact pair caps and the four-slot layer

For each of the `4,866` THM-2915 failures, the verifier reconstructs the
nine-fixed-speed child both by sequential subtraction and directly from the
full family.  It retains the inherited prefix, first centre, and all earlier
first centres as `F`.

The exact pair cap is computed from a finite head and an infinite-tail seal.
At horizon `M`, if `q_1(M)` is the finite leading singleton and

```text
c(w)<=h/7+gamma/(M+1)                for every omitted w,   (12)
```

write `tau=h/7+gamma/(M+1)`.  A finite-head label paired with an omitted label
has union at most `q_1(M)+tau`.  At the sealing point the retained finite
winner `B_M` satisfies `B_M>=q_1(M)+tau`; subadditivity also gives
`B_M<=2q_1(M)`, hence `tau<=q_1(M)`.  A pair of two omitted labels therefore
has union at most `2tau<=q_1(M)+tau` as well.  Thus every omitted pair is
bounded by the same finite-tail quantity, and the retained winner is accepted
only when it weakly dominates it.  Equality in the finite singleton-sum
pruning is evaluated, not discarded.

Across all children, the maximum pair horizon is `4,056`, six horizon
extensions occur, `40,753` finite pairs are paid, and the minimum strict
winner-tail gap is

```text
50423/214911136440.                                        (13)
```

Partitioning four labels into two pairs proves that `2B<h` closes `1,612`
children.  The minimum positive margin is

```text
103/40060020.                                              (14)
```

Among the remaining `3,254`, the envelope `G_4<h` closes another `1,175`.
Its minimum positive margin is

```text
4573/1242639216.                                           (15)
```

No equality occurs in either test.  For all `2,079` remaining children,
the exact first-crossing gap is positive; its minimum is

```text
lambda-h/7 >= 792063/74554480.                             (16)
```

The resulting second-centre cores have size histogram

```text
size             1    2    3    4    5   6  7
children        214  525  664  504  137  30  5,           (17)
```

with maximum cutoff `626`.  Their `6,172` ordered pivots are exactly the
cores prescribed by `(10)` and `(11)`.

After fixing a second centre, exactly three labels remain.  The verifier
reconstructs every ten-fixed-speed grandchild, seals the exact global top
three at horizon at most `2,000`, and checks their sum.  Any three-label
completion has covered union at most `q_1+q_2+q_3`, so the strict inequality
`q_1+q_2+q_3<h` closes that pivot; equality or failure is retained for the
next layer.  The minimum tail gap is

```text
1067146889/381167448480.                                   (18)
```

Exactly `5,944` pivots close strictly, with minimum positive margin
`1/697680`.  Hence `1,884` whole first children close and `195` retain at
least one of the `228` failed pivots.

## 4. Three-slot compatibility and the terminal pair

At a failed second pivot, there are three remaining labels.  Computing a pair
statistic does not itself consume a slot; this explicit arity is the repair
mandated by `MISTAKE-324`.

For each twice-subtracted carrier, recompute the exact global pair cap `B` and
the singleton ranks `q_1>=q_2>=q_3`.  The specialization

```text
q_1+B<h                                                   (19)
```

is valid because one member of any three-label completion contributes at most
`q_1` and the other two have union at most `B`; it closes `65` pivots.  The
lawful three-slot envelope `G_3<h` closes another `114`.  Their smallest
positive margins are respectively

```text
10679/986305320,        160063/864913896.                 (20)
```

Neither route has an equality.  The `49` remaining rows again have
`lambda_3>h/7`, with minimum gap

```text
5220695/427756329.                                        (21)
```

Their third-centre cores have histogram

```text
size             1   2   3   4  5
rows             8  20  16   4  1,                       (22)
```

maximum cutoff `441`, and `117` ordered pivots.

Fixing the third centre leaves exactly two labels.  The verifier computes the
exact global pair cap of every eleven-fixed-speed great-grandchild.  A strict
`B<h` closes `116` pivots; its minimum positive margin is

```text
7017779/1460244240.                                       (23)
```

The largest pair horizon at the three-slot/terminal layers is `4,259`, and
the minimum pair-tail gap is

```text
507443/2217543127800.                                     (24)
```

This upgrades `227/228` failed second pivots and `194/195` first children.
The sole remaining terminal has margin exactly zero and identity

```text
E=(2,4,6,8,10,12,14), rank=1, apex=22, prefix=(22),
x=18, y=26, z=16, pair witness=(20,24),                   (25)
```

which reconstructs `(3)`.

## 5. The endpoint equality

At time `t=1/28`, the thirteen speeds in `(3)` have

```text
min_v ||tv||=1/14,
strict-danger owners=empty,
boundary owners={2,26}.                                   (26)
```

After dividing the speeds by two, the residues are exactly
`1/14,2/14,...,13/14`.  More precisely, for every
`0<epsilon<1/364`,

```text
||2(1/28-epsilon)||  = 1/14-2 epsilon,
||26(1/28+epsilon)|| = 1/14-26 epsilon.                  (26a)
```

Thus the minimum immediately on either side is below `1/14`, so `(26)` is an
isolated closed-good endpoint, not a hidden open interval.  The endpoint
auditor asserts both formulas at its exact rational epsilon.  The full family
occurs in THM-2907's independently generated two-family endpoint bank.

Thus the zero terminal is not promoted as a strict pair inequality.  It is
joined to the previously proved endpoint route `E`; this distinction is why
the strict/open convention is part of the theorem rather than formatting.

## 6. Identity-complete route composition

The composition checker does not trust route flags.  It reparses exact
margins and derives every top-four, pair, `G_4`, top-three, pair-plus-singleton,
`G_3`, and terminal-pair decision.  It then performs exact joins on:

```text
THM-2915 children                                      51,222
THM-2915 failed children                                4,866
first recursive rows                                    2,079
second pivots                                            6,172
failed second pivots                                       228
third recursive rows                                        49
third pivots                                                117
parent rows                                             11,842.             (27)
```

Write `P` for the `279` literal branch-closed parent rows, `E` for the `52`
THM-2907 exception rows, and `T2` for the `11,562` open-centre rows completed
by THM-2915 plus the recursion above.  Before applying `E`, exactly one parent
remains: `(25)`'s parent.  It lies in `E minus T2`; the other `51` exception
rows lie in `E intersect T2`.  Therefore

```text
P union E union T2 = all 11,842 parent rows.              (28)
```

The branch-atom census is

```text
E:1, ET:51, HP:215, HQT:1, HT:2874, P:64, QT:37, T:8599. (29)
```

Projecting `(28)` closes `3,411` bodies.  The proved `351`-root baseline
through THM-2913 overlaps these in `330` bodies and supplies the remaining
`21`, proving the exact union `(4)`.  The route-root semantic digest is

```text
75c4e9f27ed81fb08250b6cf808ed8a0db762031246db9358dc287a081a51aa3.  (30)
```

The narrower THM-2913 and THM-2920 implementations are positive controls, not
substitutes for the complete join.  All `42` THM-2913 child identities and all
`367` THM-2920 child identities occur here with identical exact coverages,
margins, and normalized tied ranks.  All `296` THM-2920 roots lie in the
`3,410`-root `P union T2` projection.

## 7. Validity and independent audit

Every proof-bearing carrier is reconstructed both sequentially and directly
from its full fixed family.  The verifier asserts the arities

```text
9 fixed + 4 remaining,
10 fixed + 3 remaining,
11 fixed + 2 remaining,
13 fixed + 0 remaining.                                  (31)
```

All fixed labels are distinct.  Every inherited prefix and earlier-centre
sidecar is disjoint from every newly fixed centre.  Pair witnesses are
re-subtracted and reconstructed as full families.  Infinite singleton and
pair tails are sealed exactly; failed sufficient certificates are never
interpreted as covers.  All comparisons use rational arithmetic.
Pair unions are obtained by literal positioned-carrier subtraction; no
full-circle overlap formula from the corrected THM-594 branch is used
(`MISTAKE-327`).

An independent hostile implementation reconstructed all `4,866` first
children, `6,172` second pivots, `228` failed grandchildren, `49` third cores,
and `117` terminal pairs.  It independently enumerated the finite heads,
recomputed every route bit from exact margins, found `116` positive terminals
and the one equality `(25)`, and reproduced the `11,842`-parent / `3,432`-root
composition.  Its mathematical verdict was **PASS**.  Four fresh frozen
replays—ordinary runs with worker counts `8` and `6`, and optimized runs with
worker counts `7` and `5`, all at distinct seeds—were byte-identical to the
frozen summary and combined ledger.  The repository-internal mathematical,
source-freeze, artifact, and replay audit verdict is **PASS**.

A separate exact reparse audited every invocation of the shared pair-cap
routine: `4,866` first-child calls, `228` grandchild calls, and `117`
third-level calls, hence `5,211` calls in total.  With `tau` the omitted
singleton tail, `q_1` the finite leading singleton, and `B` the retained
finite pair winner, every call satisfies

```text
2 tau <= q_1+tau <= B <= 2q_1.                          (31a)
```

The committed routine now asserts the first two inequalities before accepting
its tail seal.  Across the frozen ledger, the minimum exact slacks in `(31a)`
are

```text
(q_1+tau)-2tau       54542843/11777645880,
B-(q_1+tau)          507443/2217543127800,
2q_1-B               19/1441440.                       (31b)
```

Thus a pair of two omitted labels, a mixed finite/omitted pair, and a retained
finite pair are all covered by the same proof-bearing cap, at every recursive
depth.  This explicitly discharges the last source-level pair-tail audit gate.

## 8. Conditional critical barrier at the next rung

The scope boundary is structural.  At a seven-slot node, suppose

```text
q_7>=h/7,                 B>=2h/7.                        (32)
```

Then every summand in `(5)` equals `h/7` at `a=h/7`, so

```text
G_7(h/7)=h.                                               (33)
```

For `a<h/7`, every one of the seven terms is at most `a`, hence
`G_7(a)<=7a<h`.  Thus the first crossing is exactly

```text
lambda_7=h/7.                                             (34)
```

The strict discrepancy gap needed in `(10)` vanishes, and THM-735 cannot turn
the hostile set into a finite core.  The corresponding six-body census is
scope telemetry awaiting its own audited canonical artifact; it is not a
claim or dependency of this theorem.

Accordingly, on any seven-slot node satisfying `(32)`, blind scalar repetition
has no THM-735-finite core and needs new overlap, compatibility, or relation
structure.  The uncommitted six-body telemetry suggests that this conditional
barrier is uniform, but that claim awaits its own audit and is not used here.

## 9. Locked replay

The top-level verifier hash-pins all five helper programs, regenerates every
stage in a temporary directory, checks the LF-normalized SHA-256 of each
handoff before it is consumed, and emits one combined identity-complete
ledger.  Its constituent semantic digests are

```text
pair layer
  2b0580c8e98e13aeb34adc2f4ec3f8da5be30363adb882cd6794da192ce9476c
second-centre layer
  c9e4a568f17778355b4ce91ee387add8f67dd6dcca257108ec945a2dfe1bd327
grandchild layer
  0b70905f638bb2a548e74a265df76506769f33424989c1385aa170b73c8a302a
third-centre layer
  55c2e237e4b22475d2daa1dd08f8d8df25653952387c62813dd33edc2f66ae67
combined ledger
  a2555163f07f6787e3a789b365ef27c47957ed39a82074ee9067b9bac89f3104
summary
  c2ef8608ac64f068e6cc83481c8b91c126c1146b4808bd21b9e42d61e08d6abb.
```

LF-normalized SHA-256 values at source freeze are

```text
top-level verifier
  2d90e8d34f66d1c624759889a8ef3563aea7e59ec77f5c8696d0641521c453b5
stage-one verifier
  32751b2a5beb789b1657f06d3964cbf24e634251e57f57f299b8f50c647f2103
stage-two verifier
  6f494ee13fd7ecc5a7a9787b75e0cf38c58a3e89b5136470b2332f5a4cece1e5
composition verifier
  7cb3384831151ac3b2dd6d3b38f185d5328ed8c13dfca757621b573a6120600d
endpoint verifier
  4c799992395a064656a2b949b2b81bfd6082e33dce14b4d92611a0c9bb7d05cb
inherited-slice verifier
  cb40a748fecad6659eb5b2a140d2d8d23a966a643bd1a3918b316341644dd78d
output
  993866a4586ead55ff63b65421b558a716ddcf74054ee13226cc9bdcb459884c
combined ledger
  2db51b183502bdb75f0f4ef68e9772566d2564e07c0543298c6de9ed513c7c0c.
```

The reproduction commands are

```bash
python3 04-computation/lrc14_j6_seven_body_six_slot_recursive_pair_hunter_closure_thm2923.py --workers 8 --hash-seed 17
python3 -O 04-computation/lrc14_j6_seven_body_six_slot_recursive_pair_hunter_closure_thm2923.py --workers 7 --hash-seed 2915
```
