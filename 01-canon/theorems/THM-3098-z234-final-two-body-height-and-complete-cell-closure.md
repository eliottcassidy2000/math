---
id: THM-3098
title: "z234 final two-body height and complete-cell closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Exact height and
  complete-cell certificates close the final two THM-3078 rows.  Composition
  with THM-3094 empties the full pinned z1=234 layer, updates the projected
  k=3 ledger to 374387, and lowers its cap to z1<=233.  It makes no LRC(14)
  claim.
source: codex-thm3094-hostile-audit-2026-08-02
audit: >
  Root independently rederived both height inequalities and every exceptional
  cutoff, checked that all scalar filters are necessary upper relaxations,
  and rebuilt the 1,018 literal pair intersections with pure-Python integer
  bitsets rather than the companion's NumPy path.  This reproduced the empty
  denominator-two box, minimum Bonferroni slack 15470, actual minimum 19908
  at (20174,48540), strict-open endpoint controls, and the atlas composition
  374768-381=374387 with next-layer census 62=45+17.  Fresh normal and
  optimized runs byte-match stored output; LF hashes and documentation pass.
depends_on:
  - THM-3078-z234-direct-farkas-normalization-and-four-two-high-boundary
  - THM-3094-z234-two-high-complete-cell-intersection-closure
  - THM-3071-projected-k3-z236-z235-compositional-descent
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-1166-seven-wall-fano-gcd-discrepancy
  - THM-1094-exact-two-comb-component-theorem
related:
  - THM-3087-z234-two-high-common-torsion-and-finite-height-reduction
script: 04-computation/lrc14_j7_k3_z234_final_two_body_height_complete_cell_closure_thm3098.py
output: 05-knowledge/results/lrc14_j7_k3_z234_final_two_body_height_complete_cell_closure_thm3098.out
bank: 05-knowledge/results/lrc14_j7_k3_z234_final_two_body_mask_bank_thm3098.tsv
script_sha256: 744be96e68469afcf1864fadedfecfa237490b6ec3819f69bac16c52c5b558ed
output_sha256: 722acc446f90a57df8ae84d50eb05d04250b1ffd688fb02e52ac0fcdd4111d07
bank_sha256: a1b3959152cc87b05d352492870803cc23e04fb6b38dd7cbf09d2e9605850144
semantic_sha256: f5427dc8392ed287162f076df9b5e418afcec2a2f3a63a3b94a840c313a5b6bd
hash_basis: LF-normalized bytes
---

# THM-3098 -- z234 final two-body height and complete-cell closure

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Statement

THM-3078 leaves four rows at the two-high boundary of the pinned THM-2941
projected `k=3`, `z_1=234` necessary atlas.  THM-3094 closes two.  The two
remaining rows are

```text
E_A=(1,5,6,9,12,14),
E_B=(2,5,9,11,12,14).                                  (1)
```

Both rows are empty.  More precisely, the complete post-THM-3078 residual
bank in `(1)` consists of

```text
E_A:    16 denominator masks;
E_B: 1,138 denominator masks;
total: 1,154 masks.                                    (2)
```

All three-high packets fail a strict scalar upper bound.  In the exactly
two-high sector, the exact upper envelope leaves `3` cases in `E_A` and `69`
in `E_B`.  A physical height argument makes both high labels finite.  The
`E_A` cases then have no denominator-compatible high label at all.  The
`E_B` cases contain exactly `1,018` ordered literal candidates, and every
one has a common complete safe cell.  Thus no packet remains in `(1)`.

The exact composition consequence is not left implicit.  THM-3078 proves
that the `381` pinned `z_1=234` rows split as `377` closed rows plus its four
boundary rows.  THM-3094 closes two boundary rows and this theorem closes
the other two.  Therefore

```text
all 381 z_1=234 rows are empty,
374768-381=374387,
projected k=3 cap: z_1<=233.                            (3)
```

The companion separately reparses the pinned atlas and confirms that the
next layer `z_1=233` is occupied by `62=45+17` wall/order rows.  No row class
or relaxation remains at `z_1=234`.

## 2. Pinned bank and scalar exhaustion

For a body `E`, let `L` be its period,

```text
H=floor(13L/132)+1,                                    (4)
```

and retain THM-3078's four-denominator mask `d`.  One denominator belongs to
the fixed first label `234`; the other three are the suffix slots.  The bank
is the sorted text projection of exactly the two corresponding THM-3078
screen-checkpoint residuals.  Its bodywise tuple hashes are

```text
E_A, 16 masks:
  fb53e6a6c3e2098560a8971a079aa595a04a95c1fcdff77e3be3a8691bfff6d5;
E_B, 1,138 masks:
  c5c1949cc5a2037ab6af71f8d41d36b6f2161fce0fd710a84315a3a3f4f1218f.
                                                                    (5)
```

Let `R` be the suffix mass required after the first label, let
`Delta(x)` be the inherited exact singleton surplus, and let `u_d` be the
rigorous high-ray supremum for denominator `d`.  The high-ray law is

```text
(x+L)Delta(x+L)=x Delta(x).                             (6)
```

THM-3078's exact unit-ray quotient computes `u_d` without a search horizon;
nonpositive ray amplitudes are safely replaced by upper bound zero.

### Three high labels

For every mask, the companion checks

```text
R-sum_(three suffix slots d) u_d > 0,                   (7)
```

even while granting repeated labels.  The exact minima are

| body | masks checked | minimum gap in `(7)` | denominator witness |
|---|---:|---:|---|
| `E_A` | `16` | `13734121/4337473140` | `(2,980,1764,1960)` |
| `E_B` | `1138` | `3429147350803/807271369268364` | `(1260,10780,16170,21560)` |

Hence every three-high assignment is impossible.  This is required because
THM-3078's failed duplicate-two-high statistic alone did not close that
subsector.

### Exactly two high labels

Choose the lone low suffix slot and enumerate every literal label
`234<x<H` of its denominator.  If

```text
epsilon=R-Delta(x),                                    (8)
```

then a necessary scalar condition for the other two slots is

```text
u_(d_1)+u_(d_2)>=epsilon.                              (9)
```

The exhaustive result is

| body | raw low-slot expansions | survivors of `(9)` | low label | `epsilon` |
|---|---:|---:|---:|---:|
| `E_A` | `2,748` | `3` | `324` | `37/3611790` |
| `E_B` | `1,512,196` | `69` | `243` | `4117/39729690` |

There is no survivor with equal high denominators.  The exact survivor
digests are respectively

```text
1fd3d18d3a4fd592a2391643fe7fdef1ed196e2a237916093e611c3bba76a87f,
ae7d6b1b307992b48ade07543354798f85ea4d4136eb1b1de01b58b959151a77. (10)
```

The positive envelope margins range over

```text
E_A: 31396/122542875 .. 26554/71461845;
E_B: 1063/30830239440 .. 138056623/3349673731404.       (11)
```

Thus `(9)` is only a necessary upper filter, not a claimed physical closure.
The remaining sections retain the literal labels and actual cells.

## 3. A physical two-stage height bound

In either body let `x` be the unique low label in the preceding table and let
`J` be the complete body `1/L` cells which are safe from both `234` and `x`.
Exact weak-endpoint tests give

```text
E_A: L=17640,  H=1738,  x=324, |J|=4446;
E_B: L=194040, H=19111, x=243, |J|=44364.              (12)
```

Choose any cell `I` represented by `J`.  THM-1166 bounds the union of the
three distinct aligned danger combs by `36/91`.  THM-1094, equation `(10)`,
applied to the interval of length `1/L`, gives for every high label `z`

```text
L mu(I intersect D_z) <= 1/7 + 6L/(49z).               (13)
```

If both high labels are at least `N`, coverage of `I` would force

```text
1 <= 36/91 + 2/7 + 12L/(49N).                          (14)
```

The right side is strictly below one when

```text
N>1092L/1421.                                          (15)
```

Therefore the smaller high label is at most the floor in the following
table.  The companion checks every integer from `H` through that weak
endpoint.  A label is ordinary when some cell of `J` is wholly safe from its
strict-open danger comb.

| body | boundary in `(15)` | smaller-label ceiling | nonordinary labels |
|---|---:|---:|---|
| `E_A` | `393120/29` | `13555` | `8820,11760,13230` |
| `E_B` | `4324320/29` | `149114` | `97020,129360,145530` |

For an ordinary smaller label `z`, choose its safe cell.  On that cell only
the other high comb and the three aligned combs remain, so `(13)` forces

```text
w<=13L/49.                                              (16)
```

The nonordinary labels are not loopholes.  Their exact least normalized
loads, attaining cells, far-label cutoffs, and strict margins already at
`w=z` are

| body | `z` | least `L mu(I intersect D_z)` | cell | far cutoff | margin at `w=z` |
|---|---:|---:|---:|---:|---:|
| `E_A` | `8820` | `1/7` | `1350` | `196560/29` | `47/637` |
| `E_A` | `11760` | `3/28` | `1350` | `262080/43` | `435/2548` |
| `E_A` | `13230` | `2/21` | `1591` | `29484/5` | `388/1911` |
| `E_B` | `97020` | `1/7` | `6930` | `2162160/29` | `47/637` |
| `E_B` | `129360` | `3/28` | `6930` | `2882880/43` | `435/2548` |
| `E_B` | `145530` | `2/21` | `6931` | `324324/5` | `388/1911` |

Every displayed far cutoff is strictly below its exceptional `z`.  Since
`z` was chosen as the smaller high label, the other high label satisfies
`w>=z`, and the displayed positive margin contradicts coverage.  Hence the
exceptional case is impossible.  Combining the ordinary and exceptional
branches gives the exact finite contracts

```text
E_A: 1738 <= z,w <= 4680;
E_B: 19111 <= z,w <= 51480.                             (17)
```

This argument is carried on actual complete intervals.  It never substitutes
a scalar safe mass for the projected set.

## 4. Two exact finite exits

### The 16-mask denominator obstruction

The three `E_A` scalar survivors have high-denominator pairs

```text
(2,2520), (2,3528), (2,17640).                          (18)
```

An integer label has denominator `2` relative to `L=17640` only in the
class `8820 mod 17640`.  There is no such label in `[1738,4680]`.  Thus all
three survivor boxes are empty.  The per-survivor count digest is

```text
f73fe2699430f02465a78ae555f36adfeead70a062d1008cbc08f7331f4b11cf. (19)
```

### The 1,138-mask literal common-cell obstruction

For `E_B`, the companion enumerates every ordered denominator-compatible
pair in `(17)`, enforces distinct labels, and applies the exact scalar
inequality using `(6)`.  The finite split is

```text
69 scalar-survivor masks,
7,366,342 raw ordered box pairs,
1,018 ordered literal candidates,
0 empty scalar-survivor masks.                           (20)
```

The `1,018` candidates are `1,018` distinct ordered label pairs, supported
on `342` labels from `19166` through `48540`.  The packet and distinct-pair
digests are

```text
3277b423430c9fe955d78e056942b0e73e0f8ec740f2f240e108451d46ceefe0,
b723c2705b0f25c04e5b07307908ccdfdd666e9ef4501d4baff1cccf41b22aff. (21)
```

For a candidate label `z`, define the literal complete-cell set

```text
A_z={j in J:[j/L,(j+1)/L] is disjoint from D_z}.        (22)
```

Because the danger comb is strict-open and every label in `(17)` is below
`L`, if `r=jz mod L` then membership in `(22)` is exactly

```text
14r>=L,                  14(r+z)<=13L.                  (23)
```

For every one of the `1,018` candidate pairs, exact `int64` evaluation of
`(23)` gives

```text
|A_z intersect A_w|
 >= |A_z|+|A_w|-|J|
 >= 15470 > 0.                                          (24)
```

The minimum Bonferroni slack in `(24)` occurs at `(z,w)=(20174,48540)`.
The actual minimum intersection is larger:

```text
|A_20174 intersect A_48540|=19908,                      (25)
```

with witness cell `6930`.  The complete intersection-record digest is

```text
4c245c14a867a219dea91a7e428fe7de2ba0cce27b1300064ba6fcbf6107fe79. (26)
```

Actual runner labels are distinct, and the finite enumeration enforces
`z!=w`; no duplicated-minimum shortcut is present.

## 5. Lossless projection and the exact direction

Choose the common complete cell supplied by `(24)`.  Its closed interval is
safe from the body, the fixed first drift `234`, the fixed low drift `243`,
and both literal high drifts.  Hence it lies in `S_(E,Z)`.  The map

```text
phi_L(t)=Lt mod 1                                       (27)
```

maps any complete `1/L` interval surjectively onto the circle.  Therefore

```text
P_(E,Z)=phi_L(S_(E,Z))=T.                               (28)
```

The logical direction used from THM-2941 `(25g)` is precisely

```text
literal completion  ==>  P_(E,Z) subset U_A.            (29)
```

It is not the converse direction and is not an inference from an intermediate
mass statistic.  THM-1166 gives, for the three distinct aligned labels,

```text
mu(U_A)<=36/91<1.                                       (30)
```

Equations `(28)--(30)` contradict one another.  This closes all candidates
in `(20)` and therefore both rows in `(1)`.

The weak endpoint inequalities in `(23)` are load-bearing.  At cell `1350`
for `E_A`, body label `14` has zero left-endpoint slack; at cell `6930` for
`E_B`, body label `2` has zero left-endpoint slack.  Replacing the weak
conditions by strict inequalities would silently discard valid complete
cells.

## 6. Composition, evidence, and scope

The companion pins the promoted THM-3078 and THM-3094 source, output, and
semantic hashes.  For `(3)` it also pins THM-3071's source, output, semantic
hash and exact predecessor ledger formula, pins THM-2941's atlas hash, and
reparses all `6,060` structured atlas rows.  It independently recovers

```text
z_1=234: 381 rows = 330 wall + 51 order;
z_1=233:  62 rows =  45 wall + 17 order.                 (31)
```

Thus the decrement in `(3)` is a disjoint set difference in the canonical
projected necessary-row ledger, not an inferred numerical slogan.

Run

```text
python 04-computation/lrc14_j7_k3_z234_final_two_body_height_complete_cell_closure_thm3098.py
python -O 04-computation/lrc14_j7_k3_z234_final_two_body_height_complete_cell_closure_thm3098.py --output <optimized-output>
```

The script contains no truth-bearing Python `assert`.  Normal, optimized,
and stored transcripts are LF-byte-identical and end in
`all_exact_controls=PASS`.

**Scope.**  The theorem empties the complete pinned projected `k=3`,
`z_1=234` layer and proves only the projected cap/ledger consequence `(3)`.
It does not classify physical covers outside this necessary projection, say
anything about arbitrary `k<=1` packets or the final rung, or prove LRC(14).

QED.
