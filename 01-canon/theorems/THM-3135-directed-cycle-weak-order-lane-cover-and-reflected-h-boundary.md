---
id: THM-3135
title: "Directed-cycle weak-order lane cover and the reflected H/H2 boundary"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Weak-order edge lanes cover every positive level assignment exactly when
  their directed graph contains a cycle.  With the corrected orientation,
  the two same-pair half lanes repair the H2 reflected tails at caps 7/3 and
  5/2.  For H, every asymptotically eligible standard single-pair half/full
  uniform-envelope lane lies in a DAG, so no family from that class can cover
  all assignments.  The H statement is only an envelope-class obstruction;
  it is not a physical survivor, a global cone closure, or a proof of LRC(14).
audit: >
  The canonical exact probe replays identically under normal and optimized
  Python.  An independent implementation checked the directed-cycle
  equivalence, all 59 and 114 unequal-channel H2 heads, all 90 ordered-pair
  and side choices for H at each cap (including whole-interval lanes), the
  admissible DAG witness, and the three negative complementary limits.
source: root/frontier-synthesis-2026-08-02
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
related:
  - THM-3129-bounded-poset-upset-facet-irredundancy
  - THM-3130-divisor-antichain-totient-response-and-supermodular-witness-loss
  - THM-3132-projected-forest-boundary-parity-and-augmentation
  - MISTAKE-347-reflected-cone-split-tail-orientation
script: 04-computation/lrc_reflected_directed_cycle_lane_boundary_thm3135.py
output: 05-knowledge/results/lrc_reflected_directed_cycle_lane_boundary_thm3135.out
script_sha256: fe007e3aa361b2e4b7abbb980d99b08d1984b4fda2f55d86d661cb251a9d9470
output_sha256: 9eae574cfbe7013b90779b048b05b0632181178dab80a8e41658a7671c41d6e3
independent_script_sha256: 017a03e2f557fdc52c2fd0315f9c91b7e0ddcc8c6f5f0301a2d85474db0645b2
independent_output_sha256: 87a4eaeb45a7fbca5af3aa82c3e1dd822498a92a0a93d5db8f732b2f072b3bf2
hash_basis: LF-normalized bytes
---

# THM-3135 -- directed-cycle weak-order lane cover and the reflected H/H2 boundary

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.**

MISTAKE-347 identified an assignment-level orientation loss in the reflected
split tails of THM-2941.  The correct sidecar is not an unordered collection
of pairwise comparisons: it is the directed graph of the weak inequalities
actually certified.  This gives an exact cover criterion, repairs the `H2`
tail, and isolates why the same move cannot repair `H` inside the standard
single-pair uniform-envelope class.

## 1. The directed-cycle cover criterion

Let `V` be a finite nonempty set and let `G=(V,E)` be a directed graph.  For
each edge define the weak-order lane

```text
U_(u->v)={q in R_(>0)^V : q_u<=q_v}.                         (1)
```

Then

```text
union_(e in E) U_e = R_(>0)^V
       iff G contains a directed cycle.                       (2)
```

Indeed, missing every lane would require `q_u>q_v` on every directed edge.
That is impossible around a directed cycle.  Conversely, if `G` is a DAG,
choose a topological ordering and a positive potential strictly decreasing
along every edge.  The resulting assignment misses every `U_e`.  This also
shows that `(2)` includes equality: equality on even one selected edge is
already covered.

The useful hostile dual is therefore exact:

```text
directed cycle  <=>  assignment cover,
DAG             <=>  a strict-potential adversary.            (3)
```

## 2. Correct typing of a reflected half lane

For an ordered physical pair `(i,j)`, the THM-2941 tail variable is

```text
Q/P=q_j/q_i.                                                  (4)
```

Consequently

```text
low  (i,j), Q/P<=1  gives q_j<=q_i, hence edge j->i;
high (i,j), Q/P>=1  gives q_i<=q_j, hence edge i->j.           (5)
```

MISTAKE-347 changed both the ordered pair and the inequality side, so its two
lanes encoded the same edge twice.  Keeping the pair fixed and using both
halves produces the two-cycle `i->j->i`, which covers all assignments by
`(2)`.  The finite unequal-channel heads below exclude `P=Q`; the inherited
same-level mechanism discharges that equality face.

## 3. Exact H2 repair

Take

```text
H2=(1,3,4,6,8,12),             L=14 lcm(H2)=336.              (6)
```

For the fixed pair `(0,1)`, the existing high lane and the repaired low lane
have the following exact data.  `N` is the uniform numerator in the inherited
THM-2941 tail envelope; `s` is its first positive tail threshold.

| cap | edge | side | `N` | `s` | unequal heads | tail margin at `s` | weakest located head |
|---:|:---:|:---:|---:|---:|---:|---:|---:|
| `7/3` | `0->1` | high | `4` | `7` | `16` | `10714776833746874957/17992296767910983278761` | `5142490788532529/780922437170070195` |
| `7/3` | `1->0` | low | `36/7` | `11` | `43` | `1073139685678090183/6348118491755500615155` | `1611487748805825701/734535644402168025417` |
| `5/2` | `0->1` | high | `4` | `8` | `27` | `201784828118192099/194587023752730087920` | `60797802549184067/10640912096266013135` |
| `5/2` | `1->0` | low | `26/5` | `14` | `87` | `45635336138525012/1043031165396943398095` | `35253803208250917/25848680979062186735` |

Every displayed rational is positive.  The finite head universes are
exhaustive and all `59` heads at cap `7/3` and all `114` heads at cap `5/2`
are positive.  The two rows at each cap form a directed two-cycle.  Hence the
corrected `H2` policy covers every positive level assignment and repairs this
body at both caps.  The old duplicated-edge policy missed `120` of the `720`
strict level orders.

## 4. The H envelope-class obstruction

Now take

```text
H=(1,2,3,4,6,12),              L=14 lcm(H)=168.               (7)
```

At each cap, exhaust all `30` ordered pairs and each of the low, high, and
whole-interval sides.  Retain every lane whose inherited uniform tail
envelope has positive asymptotic limit.  After identifying duplicate
physical inequalities, the complete eligible edge graphs are

```text
cap 2:    0->1->2->3,
cap 7/3:  0->1->2->3,
cap 5/2:  0->1->2.                                      (8)
```

There is no eligible whole-interval lane.  Each graph in `(8)` is a DAG, so
`(2)` proves that no subset of these lanes covers all assignments.  The
single six-distinct assignment

```text
q=(9,8,7,6,10,11)                                         (9)
```

has minimum level `6`, span `11/6<2`, and strictly decreases on every edge
in `(8)`.  It is therefore an admissible hostile witness simultaneously for
all three caps.

The most direct attempt to create a cycle is to add the complementary low
lane on the existing pair `(0,1)`.  Its exact asymptotic envelope limits are

```text
cap 2:   -334872435363311/231417722241366920,
cap 7/3: -39873151682395/12844735947963612,
cap 5/2: -19285497107/4139186882565.                         (10)
```

All are negative, so increasing the tail threshold can never certify that
lane in the displayed uniform envelope.  Thus `H` is the first structural
obstruction to this repair mechanism: the available comparison graph retains
a global potential, while the missing reverse edge is analytically insolvent
in the inherited one-pair bound.

## 5. Forest/cycle operational duality

THM-3132 shows that a forest boundary is the efficient full-rank carrier for
signed difference data (up to the odd augmentation obstruction).  Here the
same acyclicity has the opposite operational meaning: it preserves a scalar
potential and therefore permits a level assignment that defeats every weak
inequality.  A cycle is redundant for linear reconstruction but indispensable
for order-cover contradiction.  Thus “forest versus cycle” is not an
intrinsic quality of a sidecar; it depends on whether the downstream operation
is reconstructing signed differences or forbidding strict descent.

## 6. Verification and scope

Run

```bash
python3 04-computation/lrc_reflected_directed_cycle_lane_boundary_thm3135.py
python3 -O 04-computation/lrc_reflected_directed_cycle_lane_boundary_thm3135.py
```

Both outputs byte-match the stored transcript.  The independent audit is
disjoint code and checks `(2)`, the orientation dictionary `(5)`, all H2
heads, all `90` H lane choices per cap, witness `(9)`, and limits `(10)`.

The promoted claims are exactly the general criterion `(2)`, the two finite
`H2` repairs, and the `H` no-cover result within the standard-half/full-
interval, single-pair uniform-envelope class inherited from THM-2941.  This
does **not** prove that `H` is a physical survivor, rule out a multi-pair or
nonuniform tail design, restore any audit-required global reflected cone,
reduce the `561`-body proved certificate-failure wedge, or prove LRC(14).
