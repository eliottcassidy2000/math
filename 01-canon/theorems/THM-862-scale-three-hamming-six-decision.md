---
id: THM-862
title: Scale-three Hamming-six common-sheet classification and exact metric plan
status: PROVED STRUCTURAL + FINITE-EXACT — the primitive c=3 common-sheet bank has 212 presentations and 1,504 unit contexts, with an exact affine toothpick-code classification and 146,912 cap-admissible first metric edges; the unbounded metric bank remains open
source: codex-2026-07-15-S15/S16 c=3 transport audit
depends_on: [THM-810, THM-815, THM-823, THM-857, THM-859, THM-860, THM-861]
related: [THM-844, THM-847, THM-858, HYP-6820]
verification:
  - 04-computation/lrc13_scale_three_hamming_six_sheet_classification_codex_S16.py
  - 05-knowledge/results/lrc13_scale_three_hamming_six_sheet_classification_codex_S16.out
---

# THM-862 — the scale-three sheet stalk is classified

Let `R subset F_13^*` have six elements and put

```text
A=3([12] minus R) union {w_r:r in R},
w_r=3r (mod 13),                 w_r>0,       w_r!=3r.   (1)
```

Assume that `A` is primitive and `M(A)<=1/13`.  At common scale three the
effective order of a replacement is

```text
D_r=3/gcd(3,w_r) in {1,3}.                              (2)
```

For `D_r=3`, write `e_r=w_r mod 3 in {+1,-1}`.  An order-one unit is
trivial.  THM-810/823's germ argument makes common-sheet coverage a necessary
condition for (1).  This theorem classifies that condition completely and
freezes the exact first layer of the remaining metric recursion.

It does **not** run that recursion to its terminals.

## Theorem

### A. Exact common-sheet bank

After excluding the all-order-one word, which would make every speed in (1)
divisible by three, the complete primitive `c=3` common-sheet bank is

| effective orders | presentations | unit contexts |
|---|---:|---:|
| `1^2 3^4` | 84 | 336 |
| `1 3^5` | 84 | 672 |
| `3^6` | 44 | 496 |
| **total** | **212** | **1,504** |

There are no rows with one, two, or three order-three colours.  In the
all-order-three stratum, `26` presentations have eight unit words and `18`
have sixteen, so

```text
26*8+18*16=496.                                        (3)
```

### B. The signed-provider law and toothpick codes

Put

```text
H={1,5,8,12},
X_-={4,5},                       X_+={8,9}.              (4)
```

At an order-three owner `o`, the owner itself supplies normalized sheet zero.
An off-owner order-three label `r` supplies a nonzero sheet exactly when

```text
r/o in X_- union X_+.                                  (5)
```

It supplies the normalized sign

```text
-e_r  if r/o in X_-,                 +e_r if r/o in X_+.(6)
```

Thus an order-three owner is covered on all three sheets exactly when both
signs occur among (6).  This is one signed not-all-equal provider obligation
per owner.  An order-one colour covers every sheet at its own owner and none
at a distinct replacement owner.

On every surviving presentation the apparently higher-arity obligations
collapse further: the valid unit words form an affine binary code cut out by
disjoint signed pair equations.  Write `e_i e_j=+` for equal signs and
`e_i e_j=-` for opposite signs.  The mixed rows normalize to

```text
C=H:             e_1 e_12=+,   e_5 e_8=+;              (7)
C=H union {2}:   e_1 e_12=+,   e_5 e_8=+,  e_2 free.   (8)
```

For six order-three labels there are exactly five multiplicative orbits:

| normalized `C` | orbit size | affine toothpicks | free labels | words | contexts |
|---|---:|---|---|---:|---:|
| `{1,2,3,4,9,10}` | 12 | `e1e2=e4e9=e3e10=+` | none | 8 | 96 |
| `{1,2,3,5,8,12}` | 12 | `e5e8=e1e12=+` | `2,3` | 16 | 192 |
| `{1,2,5,6,8,11}` | 12 | `e1e6=-`, `e5e8=e2e11=+` | none | 8 | 96 |
| `{1,2,5,8,11,12}` | 6 | `e5e8=e1e12=+` | `2,11` | 16 | 96 |
| `{1,3,4,9,10,12}` | 2 | `e4e9=e3e10=e1e12=+` | none | 8 | 16 |

The last row is the quadratic-residue six-set and its nonresidue mate.  The
pair equations in (7)--(8) and the table are exact descriptions of the full
unit fibres, not only necessary parity tests.  This matching-code form is a
scale-three **toothpick code**: every live sheet stalk is assembled from
disjoint signed two-vertex constraints and free stalks.  The name refers only
to this matching decomposition; it does not revive THM-841's false local
toothpick recurrence for the Farey ladder.

### C. Exact multiplicative and reflection orbits

Let `T_a`, `a in F_13^*`, multiply every replacement label by `a` while
leaving its order and attached unit fixed.  Let

```text
J:e_r -> -e_r                                             (9)
```

on every order-three colour.  Both actions preserve common-sheet coverage and
commute, so the sheet action is `F_13^* x <J>`, of order twenty-four.

The exact orbit census is:

| objects | `1^2 3^4` orbit sizes | `1 3^5` orbit sizes | `3^6` orbit sizes | total orbits |
|---|---|---|---|---:|
| presentations under `F_13^*` | `6x12, 2x6` | `7x12` | `3x12, 1x6, 1x2` | 20 |
| contexts under `F_13^*` | `24x12, 8x6` | `56x12` | `36x12, 10x6, 2x2` | 136 |
| contexts under `F_13^* x <J>` | `12x24, 4x12` | `28x24` | `18x24, 5x12, 1x4` | 68 |

The reflection `J` is free, and no multiplier composed with `J` stabilizes a
context.  Hence it pairs the `136` multiplier-only orbits into `68` sheet
orbits.

These `68` orbits have the following compact canonical schemes.  Unit words
use digits `1,2` in increasing order of the displayed order-three set.

```text
1^2 3^4:
  C={1,5,8,12},
  B in {{2,3},{2,4},{2,6},{2,7},{2,9},{2,11},{4,6},{4,9}},
  e in {1111,1221};                                     16 reps

1 3^5:
  C={1,2,5,8,12},
  B in {3,4,6,7,9,10,11},
  e in {11111,11221,12111,12221};                       28 reps

3^6:
  C={1,2,3,4,9,10},     e in {111111,111221,112112,112222};
  C={1,2,3,5,8,12},     e in {111111,111221,112111,112221,
                                121111,121221,122111,122221};
  C={1,2,5,6,8,11},     e in {111211,112221,121212,122222};
  C={1,2,5,8,11,12},    e in {111111,111121,112211,112221,
                                121121,122221};
  C={1,3,4,9,10,12},    e in {111111,112211}.           24 reps
```

The two `1^2 3^4` label pairs `{2,11}` and `{4,9}` have stabilizer
`{+/-1}`.  In the fourth all-order-three row, the words
`111111,112211,121121,122221` have the same stabilizer; its other two words
have trivial stabilizer.  On the quadratic-residue row, `111111` has the
six-element quadratic-residue stabilizer and `112211` has `{+/-1}`.  All
other displayed context representatives have trivial stabilizer.

Multiplicative inversion is **not** an action on this bank.  It maps back only

```text
86/212 presentations = all 84 rows of type 1^2 3^4 plus 2 QR rows,
352/1504 contexts     = all 336 contexts of that type plus 16 QR contexts.
                                                               (10)
```

It maps no `1 3^5` row back to the bank.  A dihedral multiplicative quotient
would therefore be false.

### D. Sheet symmetry is not metric symmetry

The `68` sheet representatives cannot replace the `1,504` arithmetic
contexts in an exact component recursion.  An explicit all-order-three
context supplies the guardrail.  For

```text
C={1,2,3,4,9,10},                     e=111111,          (11)
```

the least proper packet, its unit reflection, and its label-doubled image are

```text
Q ={1,4,15,16,18,19,21,22,24,25,33,36},
JQ={14,15,17,18,21,24,29,32,33,35,36,38},
T2Q={3,9,19,25,27,28,30,31,33,34,36,37}.               (12)
```

Exact breakpoint evaluation gives

```text
M(Q)=1/4,             M(JQ)=7/26,             M(T2Q)=14/65. (13)
```

The respective witness pairs are

```text
{7/20,13/20},         {1/52,51/52},           {32/65,33/65}. (14)
```

Thus both `T_a` and `J` are only sheet-equivariant.

### E. Exact first metric layer and the remaining run

For a fixed context, the complete proper replacement rays are

```text
D_r=1:  w_r=3r+39k,                     k>=1,
D_r=3:  w_r=u_r+39k,                    k>=0,           (15)
```

where `u_r` is the CRT representative in `[1,38]` satisfying

```text
u_r=3r (mod 13),                    u_r=e_r (mod 3).     (16)
```

Numerically order the six replacements.  At a prefix with `j` replacements,
longest strict-safe component length `L_j`, and `s=6-j` later combs, every
tight continuation obeys THM-815's cap

```text
x_(j+1)<=floor(22s/[13(13-2s)L_j]).                    (17)
```

There are `110` distinct retained-core roots among the `1,504` contexts.  By
stratum there are `84,66,44` distinct roots; these sets overlap across
strata.  Direct exact root-component evaluation gives:

| type | contexts | first `D=1` edges | first `D=3` edges | all first edges | per-context range |
|---|---:|---:|---:|---:|---:|
| `1^2 3^4` | 336 | 10,440 | 22,268 | 32,708 | 53--161 |
| `1 3^5` | 672 | 10,200 | 54,396 | 64,596 | 66--152 |
| `3^6` | 496 | 0 | 49,608 | 49,608 | 68--183 |
| **total** | **1,504** | **20,640** | **126,272** | **146,912** |  |

The exact root real-cap range on this survivor bank is

```text
4752/13 <= B_root <= 1188.                              (18)
```

The `146,912` count and (18) are theorem outputs.  They are not sampled height
cutoffs.

For execution planning only, scale THM-861's frozen `c=2` ratios by this
first layer.  Its depth-two ratio predicts

```text
146912*(641866/6266)=47148908896/3133
                         about 15,049,125 depth-two nodes, (19)
```

while its logical nodes per root predict

```text
1504*(41882982/64)=984,250,077 logical nodes.            (20)
```

Equations (19)--(20) are workload estimates, not `c=3` censuses.  They say to
run a complete depth-two scout before an overnight full tree and to shard the
`1,504` canonical contexts.  A proof run must retain every arithmetic lane,
use THM-857's complete-tooth and streaming-cap certificates, and be checked
by an independent closed-danger reconstruction.

## Proof

### 1. Derivation of the signed-provider law

Fix an order-three provider `r`, owner `o`, and put `a=r/o in F_13^*`.
Let `u` be the CRT class with

```text
u=3r (mod 13),                         u=e_r (mod 3).    (21)
```

On sheet `ell`, coverage asks that the centred residue

```text
z=<u(o^(-1)+13ell)>_(39)                               (22)
```

satisfy `-3<z<=3`.  Modulo thirteen, `z=3a`.  The only possibilities are

| `a` | 1 | 4 | 5 | 8 | 9 |
|---:|---:|---:|---:|---:|---:|
| `z` in `(-3,3]` | 3 | -1 | 2 | -2 | 1 |
| `z mod 3` | 0 | -1 | -1 | +1 | +1 |

The ratio `a=12` would give `z=-3` and is excluded by the left-open
orientation.  No other ratio hits the window.  Reducing (22) modulo three
gives

```text
z=e_r(o^(-1)+ell) (mod 3).
```

Since `e_r^(-1)=e_r`, the normalized owner sheet

```text
y=ell+o^(-1) (mod 3)                                  (23)
```

is precisely zero for the self-provider and (6) for an off-provider.  All
three sheets occur exactly under the signed NAE condition.

At order one, the oriented interval contains only its own nonzero residue.
After lifting to common scale three, that colour fills all three sheets at
its own owner and none at another owner.  This proves the symbolic law.

### 2. Classification and affine collapse

Let `C` be the order-three labels and `B=R minus C` the order-one labels.
If `C` is empty, every core and replacement speed in (1) is divisible by
three, contradicting primitivity.

If `|C|=4`, the four order-three colours must cover their four owners without
help from `B`.  THM-810's equality classification gives

```text
C=aH,                         e_r=e_(-r) on C.           (24)
```

There are three cosets, `C(8,2)=28` choices for `B`, and four unit words,
giving `84` presentations and `336` contexts.

If `|C|=5`, the five order-three colours form THM-823's all-order-three flag:

```text
C=aH union {b},                     b in 2aH.            (25)
```

The two antipodal equations on `aH` remain and `e_b` is free.  Choose the
coset, the forward flag, and the one order-one label outside `C`; this gives

```text
3*4*7=84 presentations,                     84*8=672 contexts. (26)
```

The signed law, evaluated over its finite unit words, rejects `|C|<=3`.
For `|C|=6`, evaluating (5)--(6) over all
`C(12,6)2^6=59,136` set/word pairs gives the five rows in Part B.  Direct CRT
masks and the signed NAE test agree.  Substituting the displayed pair
equations reproduces the complete direct-mask word set on every representative.
Multiplication supplies the stated orbit sizes, proving (3) and Part A.

The verifier checks the direct-mask/NAE identity on all

```text
sum_(k=1)^6 C(12,k)2^k=94,448                            (27)
```

order-three set/word states before decorating by `B`.  It then constructs all
presentations and contexts and independently tests every affine code.

### 3. What the sheet actions preserve

Let `q_o` be the standard integer inverse of `o mod 13`, reduced modulo three.
Equation (23) shows that `T_a` leaves every provider ratio and unit unchanged
but translates the literal mask at owner `o` by

```text
q_o-q_(ao) (mod 3).                                    (28)
```

The unit flip `J` reflects it by

```text
ell -> -ell-2q_o (mod 3).                              (29)
```

These owner-dependent affine relabellings preserve whether the three sheets
are covered, proving the group action in Part C.  Orbit-stabilizer calculation
gives the three exact rows there.  Inversion replaces provider ratios by their
inverses, but the support `{4,5,8,9}` is not inverse-stable; the exact failures
in (10) follow.

Because (28)--(29) depend on the owner, neither is induced by one circle map
on the continuous strict-safe set.  The exact values (13) prove this failure
without relying only on that observation.  The verifier obtains them from
the complete piecewise-linear candidate set with denominators
`2u`, `u+v`, and `|u-v|`.

### 4. Exact recursion plan

The CRT gives (15)--(16).  Distinct labels have distinct residues modulo
thirteen, so every terminal has one numerical replacement order.  Before the
sixth insertion there are at most eleven speeds, so the lower-runner theorem
makes the strict-safe prefix nonempty.  THM-815 then gives (17), exactly as in
THM-857 and THM-861.

At the root, `E(3P)` is the threefold inverse image of `E(P)`, so its component
lengths are those of `E(P)` divided by three.  Enumerating (15) under the exact
root cap gives the table in Part E and the range (18).

There is one legitimate dynamic quotient.  While a numerical prefix contains
only core and order-one speeds, all its speeds remain divisible by three, so
the `Z/3Z` deck and the scale-one quotient are exact.  An order-three insertion
is not `Z/3Z`-deck-invariant, so it invalidates that quotient as a certified
continuation and must materialize literal components.  More generally,
translation by the reciprocal of the current prefix gcd, time reversal
`t->-t`, common dilation, and permutation of already inserted speeds are
genuine metric equivariances.  Multiplicative residue relabelling, unit reflection,
provider-graph isomorphism, and inversion are not.

## Tournament analysis and challenged vertices

Use the order-three labels as vertices and draw the signed edge

```text
r -> o   iff r/o in {4,5,8,9},
sign=- on {4,5},              sign=+ on {8,9}.          (30)
```

For each canonical order-three set, record the number of edges, the number of
unordered pairs carrying zero/one/two arrows, the SCC sizes, and code size:

| `C` | edges | pairs `0/1/2` | SCC sizes | code size |
|---|---:|---|---|---:|
| `{1,5,8,12}` | 8 | `2/0/4` | `4` | 4 |
| `{1,2,5,8,12}` | 10 | `4/2/4` | `4,1` | 8 |
| `{1,2,3,4,9,10}` | 12 | `5/8/2` | `6` | 8 |
| `{1,2,3,5,8,12}` | 14 | `6/4/5` | `4,2` | 16 |
| `{1,2,5,6,8,11}` | 12 | `5/8/2` | `6` | 8 |
| `{1,2,5,8,11,12}` | 12 | `7/4/4` | `4,1,1` | 16 |
| `{1,3,4,9,10,12}` | 12 | `3/12/0` | `6` | 8 |

Every row has equally many positive and negative edges.  None is a tournament:
some pairs are absent or bidirected.  The quadratic-residue row is closest,
with one arrow on twelve pairs and three absent pairs.  Forcing the absent or
double pairs into a tournament destroys the NAE predicate.  Even the coarse
signed graph fingerprint fails to distinguish the third and fifth
all-order-three rows; their labelled signed incidence and affine codes are the
sharper carrier.

Challenge the runner-vertex assumption in stages:

```text
runners
 -> owner-sheet obligations
 -> signed provider incidences
 -> affine toothpick equations
 -> literal strict-safe components and remaining rays.  (31)
```

The matching code preserves the complete common-sheet predicate and destroys
continuous component geometry, lift height, and numerical comb order.  The
faithful object for the next computation is therefore an **orbit bundle**:

```text
(one of 68 signed sheet templates,
 every arithmetic lane in its 1,504-context fibre,
 literal strict-safe components,
 remaining labelled step-39 rays,
 last numerical speed,
 exact shortcut ancestry).                             (32)
```

This is also the Kakeya-needle guardrail.  The provider graph describes how
periodic danger needles meet the three local sheet germs; it does not decide
whether their full continuous combs cover the circle.

## Scope guardrail and reproduction

This theorem proves the complete primitive `c=3` **common-sheet**
classification, the affine unit codes, the exact sheet orbit bank, and the
first metric layer.  It does not prove that any of the `1,504` languages is
strictly loose, count covering terminals, close `c>=4`, or prove global
`n=12` sporadic-branch emptiness.

Reproduce the frozen output with

```bash
python3 \
  04-computation/lrc13_scale_three_hamming_six_sheet_classification_codex_S16.py \
  > /tmp/thm862.out
cmp /tmp/thm862.out \
  05-knowledge/results/lrc13_scale_three_hamming_six_sheet_classification_codex_S16.out

python3 -O \
  04-computation/lrc13_scale_three_hamming_six_sheet_classification_codex_S16.py \
  | cmp - \
  05-knowledge/results/lrc13_scale_three_hamming_six_sheet_classification_codex_S16.out
```

Frozen SHA-256 values are

```text
source  f87b0436c509313d5d90f5f9abe2d06b591e9ab639b5988bf206e4115cd7f39b
output  c8413de89655592b5009aad83596330750d9b5ca9cb407af692fd06f5e353ba8
```
