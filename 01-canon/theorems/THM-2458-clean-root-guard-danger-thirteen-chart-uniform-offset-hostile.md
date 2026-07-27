---
id: THM-2458
title: "Clean-root guard-danger thirteen-chart uniform-offset hostile"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the
  exact fixed-gate clean C_13 cover model, every
  complete atom has at most four roots, so at least four charts are
  needed for a positive uniform owner-minus-replica offset. The 4,430
  ordinary-danger step signatures all have exact rational span
  obstructions. Of 1,949 guard-danger signatures, 215 pass the
  rational span test. The signature (2,5,1,3,5) has an explicit
  positive thirteen-chart atlas; its guard circulant has determinant
  four, so the unique uniform weights are all 1/4 and every chart is
  necessary. Adding one shared guard chart gives a nonconstant
  delta-plus-six-replicas THM-2449 table with exact overlap word
  (0,0,1,1,0,1,1,1,1,0,1,1,0). This is a static fixed-residue
  root-atlas hostile; proved THM-2462 separately realizes all
  thirteen charts on one fixed-speed phase orbit. No full scalar-row
  embedding, linked semantic word, clock, endpoint current, or LRC
  row is constructed.
source: codex-2026-07-26-clean-root-uniform-offset-atlas
depends_on:
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
  - THM-2456-two-root-replica-uniform-offset-boundary
related:
  - THM-2401-common-filter-endpoint-or-first-death-certificate
  - THM-2452-indicator-idempotent-aggregate-endpoint-restoration
  - THM-2456-two-root-replica-uniform-offset-boundary
  - THM-2457-complete-atom-root-cosupport-graph-and-semantic-word-hostile
  - THM-2462-mixed-radix-root-phase-orbit-universality
script: 04-computation/lrc14_clean_root_guard_danger_thirteen_chart_thm2458.py
output: 05-knowledge/results/lrc14_clean_root_guard_danger_thirteen_chart_thm2458.out
script_sha256: 1cabfd915453d43d762a0a1d34e8611f48d20a17825aa4195618532637a4a90d
output_sha256: 1d695057063fbe4b10e00de0d61514365424db69a22c3828f73de238cb2798fb
hash_basis: working-tree bytes (LF)
---

# THM-2458 -- clean-root guard-danger uniform-offset atlas

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2449 identifies the only mixed-zero source/target table after the
positive owner anchor: one owner row plus six identical replicas.
THM-2456 then shows that an invertible two-root convolution turns this
boundary into a uniform owner-minus-replica density offset. Its
two-chart hostile is an abstract Boolean-chart model; it does not
retain the clean root-cover geometry of the current LRC packet.

This theorem installs that geometry exactly. The outcome has two
halves:

```text
ordinary-danger complete atom
  -> every one of 4,430 fixed-step atlases is span-obstructed;

guard-danger complete atom
  -> a fixed-step thirteen-chart hostile survives.               (1)
```

Thus the clean pointwise cover is genuinely stronger than the abstract
two-chart model, but it still does not force a THM-2449 rectangle.

## 1. Exact clean-cover universe

Work on

```text
G=C_13,

K={0,1}.                                                         (2)
```

Here `K` is the factored ordinary two-root danger gate. A **clean
cover** consists of:

- a four-root guard AP

  ```text
  G_(g,d)={g,g+d,g+2d,g+3d},       1<=d<=6;
  ```

- four genuine two-root masks `P_0,...,P_3`;
- every `P_i` disjoint from `K union G_(g,d)`;
- all thirteen roots covered; and
- total incidence fourteen.

Consequently exactly one root has multiplicity two and every other
root has multiplicity one. The four low masks represent the retained
ordinary role `q_*` and the other three ordinary roles. Their labels
are assigned only after the unlabelled cover is constructed.

There are `15,120` distinct unordered mask families:

```text
3,780  with the double root in K intersect G;

11,340 with K and G disjoint and the double root shared by two
       low masks.                                                (3)
```

The deterministic repeated-root recursion inherited by the exact
companion has `18,900` traversal records:

```text
3,780 + 15,120.                                                   (4)
```

The difference is not mathematics. In the disjoint branch the
repeated uncovered root has two indistinguishable copies; the
historical recursion traverses sixty records where canonical
unordered pairing gives forty-five. The theorem's signature atlases
are sets of supports, so this traversal multiplicity disappears.
Equation (4) is retained as a reproducibility census, not mislabeled
as a distinct-cover count.

This exact universe is an additional clean-root model. No claim is
made yet that all of its abstract charts arise from one fixed integer
speed tuple.

## 2. Complete-atom support classification

Choose one low mask as `Q=P_(q_*)` and require that retained `q_*` be
safe. The remaining four root-moving truth bits are:

```text
H, O_1, O_2, O_3,
```

where `H` is the guard and the `O_i` are the other ordinary roles.
The two blocker bits are root-constant on the chart and do not change
the support analysis.

For a truth word `epsilon in {0,1}^4`, with `1` meaning dangerous,
let

```text
C_epsilon
 =Q^c
   intersect
   intersection_(R in {H,O_1,O_2,O_3})
     (R if epsilon_R=1 else R^c).                                (5)
```

This is the complete-atom support **before** the factored gate `K` is
applied. The excess-one cover gives:

1. If `|epsilon|=0`, then

   ```text
   C_epsilon subset K,
   ```

   so these supports cannot positively combine to a uniform vector
   without changing the fixed truth word.

2. If `|epsilon|>=3`, then `C_epsilon` is empty, because no root has
   total incidence three.

3. If `|epsilon|=2`, then `C_epsilon` is either empty or the unique
   double root. It is disjoint from `K`; otherwise the root would have
   incidence at least three.

4. If `|epsilon|=1` and the dangerous role is ordinary, the support
   has size one or two. It has size one exactly when that dangerous
   low mask contains the unique low-pair double root; otherwise both
   of its roots survive. Across the `15,120` canonical mask families,
   with `q_*` and the dangerous ordinary role labelled, the exact
   census is

   ```text
   support size one:    68,040;

   support size two:   113,400.                                  (6)
   ```

5. If `|epsilon|=1` and the dangerous role is the guard, the support
   is the four-root guard AP.

In particular,

```text
|C_epsilon|<=4                                                     (7)
```

for every complete atom. A positive uniform combination must cover
all thirteen roots, so at least

```text
ceil(13/4)=4                                                       (8)
```

charts are necessary. Thus the literal two-chart hostile from
THM-2456 cannot be embedded in this exact clean-cover model.

## 3. Fixed-step signature atlases and exact duals

For a two-root mask, write its unsigned chord step in
`{1,...,6}`. A guard-danger signature is

```text
(d_q,d_H; sort(d_1,d_2,d_3)),                                     (9)
```

and an ordinary-danger signature is

```text
(d_q,d_H,d_D; sort(d_S,d_T)).                                    (10)
```

Translations vary, but every step in (9) or (10) is fixed. Exact
enumeration of the clean-cover universe gives:

```text
guard-danger signatures:       1,949;

ordinary-danger signatures:    4,430.                             (11)
```

For each signature `sigma`, let `M_sigma` be the `13 by N_sigma`
zero-one matrix whose columns are its distinct possible supports
`1_C`. A positive uniform-offset model would require

```text
M_sigma alpha=h 1,              alpha>=0, h>0.                    (12)
```

Every one of the `4,430` ordinary-danger signatures already fails
the weaker rational-span test:

```text
1 notin col_Q(M_sigma).                                          (13)
```

The companion obtains a deterministic exact dual row
`y_sigma in Q^13` with

```text
y_sigma^T M_sigma=0,

y_sigma^T 1 !=0                                                   (14)
```

for every ordinary signature and verifies (14) support by support.
This is stronger than a numerical LP failure and contains no
floating-point tolerance.

Among the `1,949` guard-danger signatures, exactly

```text
215                                                               (15)
```

pass the rational-span test. Equation (15) is not a claim that all
`215` admit a nonnegative solution. The next section gives one exact
positive member, which is all that is needed for the hostile.

## 4. A fixed-signature thirteen-chart atlas

Fix the signature

```text
(d_q,d_H,d_A,d_B,d_C)=(2,5,1,3,5).                               (16)
```

Keep `K={0,1}` and use oriented starts

```text
G_g={g,g+5,g+10,g+15},

Q_q={q,q+2},

A_a={a,a+1},

B_b={b,b+3},

C_c={c,c+5}.                                                      (17)
```

The following thirteen rows are exact clean covers:

| `g` | `q` | `a` | `b` | `c` |
|---:|---:|---:|---:|---:|
| 0 | 9 | 7 | 3 | 12 |
| 1 | 2 | 9 | 5 | 7 |
| 2 | 3 | 10 | 6 | 3 |
| 3 | 2 | 10 | 6 | 7 |
| 4 | 10 | 7 | 2 | 11 |
| 5 | 6 | 3 | 9 | 11 |
| 6 | 2 | 9 | 2 | 7 |
| 7 | 8 | 2 | 2 | 6 |
| 8 | 2 | 6 | 9 | 11 |
| 9 | 8 | 3 | 2 | 7 |
| 10 | 9 | 4 | 3 | 3 |
| 11 | 2 | 6 | 9 | 5 |
| 12 | 8 | 6 | 2 | 11 |

For every row:

- `K,G_g,Q_q,A_a,B_b,C_c` cover `C_13`;
- the incidence multiset is `2,1,...,1`;
- `Q_q,A_a,B_b,C_c` are all disjoint from `K union G_g`; and hence
- the guard-danger, `q_*`-safe, other-ordinary-safe support is exactly
  `G_g`.

As `g` varies, every root belongs to exactly four of the displayed
supports. Therefore

```text
sum_(g=0)^12 (1/4)1_(G_g)=1.                                     (18)
```

This is a positive uniform-offset atlas.

## 5. Thirteen is necessary for this fixed signature

Let

```text
M_(r,g)=1_(r in G_g),                r,g in C_13.                  (19)
```

The matrix is circulant. Exact elimination gives

```text
det M=4.                                                         (20)
```

Equivalently, its multiplier is

```text
1+z^5+z^10+z^15.
```

If it vanished at a thirteenth root `z`, with `z^5!=1`, the geometric
sum would force `(z^5)^4=1`. Since `gcd(4,13)=1`, this would force
`z=1`, where the multiplier is four. Thus it never vanishes.

The unique solution of

```text
M alpha=1
```

is

```text
alpha_g=1/4                 for every g.                            (21)
```

All thirteen weights are nonzero. Hence every positive uniform model
using this fixed signature and these possible guard supports needs all
thirteen charts. The global capacity lower bound (8) and the exact
fixed-signature necessity (21) are different statements.

## 6. The delta-plus-six-replicas hostile

Choose the shared chart

```text
C_*=G_6={3,6,8,11}.                                               (22)
```

For the translated factored gate

```text
K_s={s,s+1},
```

put

```text
w_s=|C_* intersect K_s|.
```

Then

```text
w=(0,0,1,1,0,1,1,1,1,0,1,1,0),          sum_s w_s=8.             (23)
```

Before applying the gate, define the averaged source densities

```text
f_0=(1+1_(C_*))/13,

f_ell=1_(C_*)/13,                    ell=1,...,6.                  (24)
```

They have the uniform difference

```text
f_0-f_ell=1/13.                                                   (25)
```

These are not merely fractional formal vectors. They have an exact
source-exclusive Boolean-chart realization:

- the thirteen owner atlas charts in Section 4, each with weight
  `1/52`;
- seven copies of the shared chart `C_*`, one for each source row,
  each with weight `1/13`; and
- one empty chart of weight `11/52`.

The weights total one. On each nonempty chart only one source row is
present.

Let

```text
(Kf)(s)=f(s)+f(s+1)
```

be the two-root convolution. Equations (23)--(25) give

```text
A_(0,s)=2/13+w_s/13,

A_(ell,s)=w_s/13,                     ell!=0.                     (26)
```

Thus

```text
A_(ell,0)=(2/13)1_(ell=0),

A_(ell,s)-A_(ell,0)-A_(0,s)+A_(0,0)=0                             (27)
```

for every nonzero source row and every target shift. The target word
`w` is nonconstant. Equation (26) is therefore an exact nontrivial
THM-2449 delta-plus-six-replicas table inside the static clean-root
atlas.

## 7. Fixed residues versus one physical phase orbit

The chord steps in (16) have inverse speed residues

```text
(K,q_*,H,A,B,C)=(1,7,8,1,9,8) mod 13.                            (28)
```

Distinct positive integer speeds can have these repeated residues.
Residue typing does not determine the translations in Section 4.

Here is the exact missing phase coordinate. For a positive
`13`-unit speed `u`, put

```text
N_u(y)=floor(uy+1/2),

delta_u(y)=uy-N_u(y) in [-1/2,1/2).                              (29)
```

Let `d=u^(-1) mod 13`. Away from endpoints, an ordinary danger mask
with oriented start `s`,

```text
{s,s+d},
```

occurs in exactly one of the two alternatives

```text
N_u(y)=-us       mod 13,       -1/2<delta_u(y)<-1/14;

N_u(y)=-us-1     mod 13,        1/14<delta_u(y)<1/2.              (30)
```

For a four-root guard mask with start `g`,

```text
{g,g+d,g+2d,g+3d},
```

the corresponding alternatives are

```text
N_u(y)=-ug-1     mod 13,       -1/2<delta_u(y)<-1/7;

N_u(y)=-ug-2     mod 13,        1/7<delta_u(y)<1/2.               (31)
```

Consequently each displayed row `(g;q,a,b,c)` has the machine-readable
nearest-integer residue code

```text
((0,-1),
 (-7q,-7q-1),
 (-8g-1,-8g-2),
 (-a,-a-1),
 (-9b,-9b-1),
 (-8c,-8c-1))                  mod 13,                            (32)
```

where the first/second entry in every ordinary pair uses the
negative/positive delta interval in (30), while the guard pair uses
the intervals in (31). The exact companion prints all thirteen codes.

As a cheap hostile probe, it exhausts the finite lift bank

```text
u_j=bar(u_j)+13t_j,                    t_j in {0,1,2},             (33)
```

requiring the two residue-one speeds to be distinct and the two
residue-eight speeds to be distinct. There are `324` fixed-speed
tuples. Exact rational open-interval intersection finds:

```text
displayed row / fixed-speed tuple incidences: 0.                  (34)
```

The companion also checks every phase interval by evaluating the
physical root mask at an exact rational midpoint. Equation (34) is an
exhaustive statement only for the bounded bank (33). It is an
exploratory negative control, not a proof that no larger fixed-speed
orbit realizes the atlas.

The theorem-level hostile is therefore deliberately scoped as:

```text
fixed residues + finite rational root charts,
```

not

```text
one fixed integer speed tuple + one Haar base orbit.              (35)
```

The later proved THM-2462 gives an explicit mixed-radix fixed-speed
realization of all thirteen rows. Thus even the phase-orbit sidecar is
insufficient. THM-2462 is related context rather than a dependency:
the static atlas and its uniform-offset table are proved here
independently of that later physicalization.

## 8. Exact companion

Run

```text
python 04-computation/lrc14_clean_root_guard_danger_thirteen_chart_thm2458.py
python -O 04-computation/lrc14_clean_root_guard_danger_thirteen_chart_thm2458.py
```

Both transcripts must reproduce

```text
05-knowledge/results/lrc14_clean_root_guard_danger_thirteen_chart_thm2458.out
```

byte for byte. The dependency-free companion uses only integers and
`fractions.Fraction`, with explicit `require` checks. It verifies:

1. the `18,900` traversal-record and `15,120` distinct-cover censuses;
2. every clean-cover incidence and all support impossibilities in
   Section 2;
3. all `4,430` ordinary rational dual certificates;
4. the `1,949/215` guard span census;
5. every row of the thirteen-chart atlas;
6. determinant `4`, the unique `1/4` solution, the overlap word, the
   source-exclusive probability mixture, the anchor, and all mixed
   rectangles;
7. all thirteen nearest-integer phase codes; and
8. the exact bounded fixed-speed probe (33)--(34).

The script contains no Python `assert`, floating point, numerical
linear algebra, or external algebra package.

An independent hostile audit rebuilt the clean-cover universe and
the rational column spaces with a separate basis implementation. It
reproduced the `15,120=3,780+11,340` cover split, the corrected
`68,040/113,400` ordinary one-danger census, the `215/0` feasible
guard/ordinary signature counts, the explicit atlas, determinant
four, the unique `1/4` weights, all phase formulas, and the bounded
`324`-lift zero-incidence control. Normal, optimized, and stored
transcripts agree byte for byte.

## 9. Scope and frontier consequence

This theorem proves that the clean excess-one root cover does not, by
itself, remove the THM-2449 uniform-offset boundary. It also sharpens
the abstract chart question:

```text
two clean charts are impossible;

ordinary-danger one-bit charts are span-obstructed;

one guard-danger signature needs and admits thirteen charts.      (36)
```

It does **not** prove that a scalar-cover packet supplies those
thirteen charts with the required weights. In particular it does not
construct or preserve:

- a single fixed integer speed tuple realizing all translations;
- one common base-parent phase orbit;
- the two blocker truth bits as physical endpoint factors;
- a THM-2305 semantic owner, pure/fork word, or delayed word;
- a source clock, graft, or source skew;
- an endpoint/cospan phase or exact `(m,X)` current.

THM-2462 removes fixed-speed phase compatibility from the debt. What
remains is embedding these six-role charts in one full scalar-cover
row with its valuation/blocker constraints, then retaining a linked
semantic owner/word and canonical endpoint phase. No scalar row is
excluded, the ledger remains `165`, and LRC(14) remains open.
