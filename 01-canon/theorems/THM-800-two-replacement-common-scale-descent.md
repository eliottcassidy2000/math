---
id: THM-800
title: Two-replacement common-scale descent and strict Hamming-two rigidity
status: PROVED (oriented splice-deck descent, using settled lower-runner LRC in the scale-one bound) + FINITE-EXACT (600,756 normalized rows)
source: codex-2026-07-15-S10 (Hamming-two continuation, extracting and repairing the S52/HYP-4103 floor certificate)
depends_on:
  - LRCUpTo13  # only the 10- and 11-speed bounds 1/11 and 1/12 in Part B
related: [THM-724, THM-769, THM-770, THM-795, HYP-4103, HYP-6775, HYP-6800, HYP-6820]
verification:
  - 04-computation/lrc13_hamming_two_common_scale_rigidity_codex_S10.py
  - 05-knowledge/results/lrc13_hamming_two_common_scale_rigidity_codex_S10.out
---

# THM-800 — Two-replacement common-scale descent and strict Hamming-two rigidity

Put

```text
delta=1/13,                         beta=2/25,
M(W)=max_(t in R/Z) min_(w in W)||wt||.
```

## Theorem

Let `c>=1` with `13` not dividing `c`, and let `r,s` be distinct members of
`[12]`.  Choose positive integers `w_r,w_s` such that

```text
w_r=cr (mod 13),      w_s=cs (mod 13),
w_r!=cr,              w_s!=cs,                         (1)
```

and set

```text
A=(c[12]\{cr,cs}) union {w_r,w_s}.                     (2)
```

Then

```text
M(A)>1/13.                                               (3)
```

More precisely, the proof has two parts.

### A. Common-scale descent at exact tightness

For `i in {r,s}` define

```text
D_i=c/gcd(c,w_i).                                       (4)
```

If, contrary to (3), `M(A)=1/13`, then

```text
D_r=D_s=1,
```

or equivalently

```text
c divides w_r  and  c divides w_s.                      (5)
```

This is uniform in `c` and both replacement heights.

### B. Sharp normalized scale-one floor

For every distinct `r,s in [12]` and every `j,k>=1`, put

```text
B_(r,s;j,k)=([12]\{r,s}) union {r+13j,s+13k}.           (6)
```

Then

```text
M(B_(r,s;j,k))>=2/25.                                   (7)
```

The constant is sharp:

```text
B_(4,6;1,1)={1,2,3,5,7,8,9,10,11,12,17,19},
M(B_(4,6;1,1))=2/25,                                    (8)
```

with maximizers `6/25` and `19/25`; at `6/25` the binders are `8` and `17`.

Part A and Part B imply (3).  They do **not** imply the stronger all-scale
bound `M(A)>=2/25`: Part A uses an oriented boundary germ specifically under
the equality assumption `M(A)=1/13`.

## 1. The oriented missing-owner splice deck

Fix an owner `r` and choose `a_r in [12]` with

```text
a_r r=1 (mod 13).
```

The `c` lifted `r`-splices are

```text
t_(r,l)=(a_r+13l)/(13c),             0<=l<c.             (9)
```

At `t_(r,l)`, the deleted runner `cr` would have signed phase `+1/13`.
If `s!=13-r`, the surviving runner `c(13-r)` is the unique core binder, at
signed phase `-1/13`; every other core runner has clearance at least `2/13`.
Moving a sufficiently small distance to the **left** sends that `-1/13` phase
outside the danger arc.  If `s=13-r`, both complementary binders were deleted,
so the core already has clearance at least `2/13` and is safe on both sides.
In either case there is an `epsilon_(r,l)>0` with

```text
(t_(r,l)-epsilon_(r,l),t_(r,l))
  subset {t:min_(v in c([12]\{r,s}))||vt||>1/13}.        (10)
```

For a positive replacement speed, its closed danger tooth covers points
immediately to the left of an endpoint phase `x` exactly when the canonical
signed phase lies in

```text
J=(-1/13,1/13].                                         (11)
```

The asymmetry matters: `+1/13` is eligible, while `-1/13` is not.  Interior
points of the arc remain dangerous on the left; at `+1/13` the left germ
enters the tooth; at `-1/13` it leaves the tooth.

If `A` were tight, the two replacements would cover every germ in (10).
Therefore, at every sheet `l`, at least one of

```text
w_r t_(r,l),       w_s t_(r,l)                           (12)
```

must lie in `J` modulo one.  This is the exact sheet--colour incidence
predicate used below.

## 2. Deck capacity leaves only orders one and two

For a fixed replacement `w`, the phases over the `r`-sheets form

```text
w a_r/(13c) + (w/c)l,             l in Z/cZ.             (13)
```

This is a translated uniform grid of order

```text
D=c/gcd(c,w),                                             (14)
```

each grid point repeated `c/D` times.  A half-open arc of length `2/13`
contains at most

```text
K(D)=ceil(2D/13)                                         (15)
```

points of a uniform `D`-grid.  Hence one replacement is eligible on at most

```text
(c/D)K(D)                                                (16)
```

sheets.

The capacity fractions satisfy

```text
K(1)/1=1,       K(2)/2=1/2,
K(D)/D<=1/3 for every D>=3.                              (17)
```

For `D=3,4,5,6`, this is immediate from `K(D)=1`.  For `D>=7`, use

```text
K(D)/D <= 2/13+1/D <= 2/13+1/7 < 1/3.
```

Two colours must cover all `c` sheets.  If neither has order one, (17) shows
that an order at least three paired with another order at least three covers
at most `2c/3`, while order two paired with order at least three covers at
most `5c/6`.  Thus

```text
some D_i=1,       or       D_r=D_s=2.                    (18)
```

## 3. Order one forces the other order to be one

Suppose `D_r=1`.  Then `w_r=cu` for a positive integer `u`, and (1) gives

```text
u=r (mod 13).                                            (19)
```

On its own `r`-sheets, this replacement has constant phase `+1/13` and is
eligible on every sheet.  On the distinct owner's `s`-sheets its constant
phase is

```text
r s^(-1)/13  (mod 1).                                   (20)
```

Because `r!=s`, this is never `+1/13`.  If `r=13-s`, it is `-1/13`, which is
excluded by the half-open orientation in (11); otherwise its clearance is at
least `2/13`.  Thus the order-one `r`-replacement is eligible on **zero**
`s`-sheets.

The `s`-replacement must consequently cover all `c` `s`-sheets by itself.
But `K(D)<D` for every `D>=2`, so (16) forces `D_s=1`.  The same argument with
`r,s` reversed proves

```text
some D_i=1  implies  D_r=D_s=1.                          (21)
```

This is where an unoriented closed-endpoint argument fails.  For example,

```text
(c,r,s,w_r,w_s)=(2,1,12,28,11)
```

has deck orders `(1,2)`, and its closed endpoint masks cover both missing-owner
sheet families.  The order-one phase at the complementary owner is exactly
`-1/13`; the left-germ carrier correctly rejects it.

## 4. The order-two pair is arithmetically inconsistent

It remains to exclude `D_r=D_s=2`.  Write `c=2d`.  From

`gcd(c,w_i)=d` and `w_i=ci+13h_i`, one obtains

```text
h_i=d e_i,        e_i odd,
w_i=d(2i+13e_i).                                        (22)
```

At an `r`-sheet, the own replacement has phase

```text
1/13 + e_r(a_r+13l)/2.                                  (23)
```

It is eligible exactly on the parity `a_r+l=0 mod 2`.  Put

```text
z=s r^(-1)  (mod 13),       1<=z<=12.                   (24)
```

On the opposite parity, the cross replacement has phase `z/13+1/2`.  Direct
reduction to the signed interval shows

```text
z/13+1/2 lies in (-1/13,1/13]  iff  z in {6,7}.          (25)
```

There is no same-parity cross rescue: `z=1` would mean `s=r`, while `z=12`
gives the excluded left endpoint `-1/13`.

Therefore coverage of the `r`-sheets requires

```text
s r^(-1) in {6,7}.                                      (26)
```

Coverage of the `s`-sheets simultaneously requires the inverse ratio to lie
in the same set.  But

```text
{6^(-1),7^(-1)}={11,2}  (mod 13),
{11,2} intersect {6,7}=empty.                            (27)
```

This contradiction eliminates the second case of (18).  Together with (21),
it proves Part A. ∎

## 5. A periodic-danger measure lemma

For `0<gamma<1/2`, let

```text
D_w(gamma)={t:||wt||<=gamma}.
```

For every interval `I` of length `L`,

```text
meas(I intersect D_w(gamma))
  <= 2 gamma L + 2 gamma/w.                              (28)
```

Indeed, after scaling by `w`, the danger indicator is one-periodic and has
integral `2gamma` on a full period.  Write `wL=N+theta`, `N` integral and
`0<=theta<1`.  The `N` full periods contribute `2gamma N`; the remaining arc
contributes at most `2gamma`.  Divide by `w` to obtain (28).

## 6. Analytic reduction of the scale-one chart

Fix `r,s,j,k` in (6), and write

```text
P=[12]\{r,s},          u=r+13j,          v=s+13k.
```

The settled ten-speed Lonely Runner bound gives a point where every member of
`P` has clearance at least `1/11`.  Since `max(P)<=12`, the Lipschitz lemma
gives a `beta`-safe interval of radius

```text
rho_10=(1/11-beta)/12=1/1100.                            (29)
```

If the two replacement combs covered this interval, (28) would force

```text
1/u+1/v >= rho_10(1-4beta)/beta = 17/2200.              (30)
```

But

```text
2/259 < 17/2200.                                        (31)
```

Consequently at least one replacement is at most `258`.  Call the smaller
one `w_b` and the larger one `w_a`.

Now adjoin `w_b` to `P`.  This eleven-speed core has a point of clearance at
least `1/12`.  Its maximum speed is `w_b`, so it has a `beta`-safe interval of
radius

```text
rho_11=(1/12-beta)/w_b=1/(300w_b).                       (32)
```

If the last replacement covered that connected interval, the interval would
have to lie in one of its teeth.  Comparing the interval length with a
`beta`-tooth gives

```text
1/(150w_b) <= 4/(25w_a),
w_a <= 24w_b.                                           (33)
```

Thus every possible counterexample to (7) lies in the finite box

```text
w_b<=258,              w_b<=w_a<=24w_b.                 (34)
```

## 7. Exact closure of the finite box

The exact replay labels the smaller replacement as `w_b=s+13k`, then ranges
over the other label `r!=s` and every `w_a=r+13j` satisfying (34).  This gives
exactly

```text
600,756 rows.                                            (35)
```

Every row is discharged by one of two integer-exact certificates.

1. In `393,962` rows, some `m in {2,...,12}` divides no speed.  Then `t=1/m`
   has clearance at least `1/m>=1/12>2/25`.

2. Each of the other `206,794` rows has a coprime rational witness `a/q`
   satisfying

   ```text
   25 min_(w in B)|aw|_q >= 2q.                         (36)
   ```

   The verifier checks the explicit denominator list

   ```text
   {13h:2<=h<=20}
     union {q:8<=q<=40, 13 does not divide q}
     union {25,50,w_a+w_b,|w_a-w_b|,w_b+1,w_a+1}.
   ```

There are zero failures.  The rational-certificate denominator histogram is

```text
{25:1, 26:5810, 39:53153, 52:61310, 65:63693,
 78:5767, 91:14531, 104:748, 117:876, 130:41,
 143:704, 156:78, 169:82}.                              (37)
```

This proves (7).  For the row in (8), an independent exact piecewise-linear
maximin scan checks every self-cusp denominator `2w` and every pair-crossing
denominator `u+v` and `|u-v|`.  It returns exactly `2/25`, with maximizers
`6/25,19/25`, proving sharpness. ∎

## 8. Completion of strict Hamming-two rigidity

Every packet in (2) is a complete nonzero residue transversal modulo `13`, so
at each nonzero thirteenth its clearance is exactly `1/13`.  Hence

```text
M(A)>=1/13.                                              (38)
```

If equality held, Part A would give `w_r=cu_r,w_s=cu_s`.  Because `c` is a
unit modulo `13`, positivity and (1) imply

```text
u_r=r+13j,       u_s=s+13k,       j,k>=1.               (39)
```

Thus `A=cB_(r,s;j,k)`.  Multiplication by `c` is onto on the circle, so

```text
M(cB)=M(B).                                              (40)
```

Part B gives `M(A)>=2/25>1/13`, contradicting equality.  Combining this with
(38) proves (3). ∎

## Tournament Analysis and assumption challenge

The theorem challenges the assumption that the useful vertices must be
runners or even unoriented endpoint sheets.

- **Exact vertices:** the `c` lifted splice sheets for each missing owner and
  the two replacement colours.  Their incidence is membership in the oriented
  half-open phase arc (11).  This preserves the exact predicate “both colours
  cover every left core-safe germ.”
- **Residue quotient:** in the order-two chart it preserves the cross-parity
  rule `{6,7}` and exposes the inverse obstruction (27).  It destroys deck
  multiplicity, order, and germ width, so it is not faithful before the
  capacity reduction.
- **Rejected bare carriers:** runner vertices lose the simultaneous `c` sheet
  obligations; fixed circle sections and section boundaries lose missing-owner
  identity; gap and cover-arc vertices lose the gcd deck; Fourier modes retain
  average density but lose the half-open boundary bit; matroid circuits retain
  dependence but not metric tooth orientation.  Proof-obligation vertices are
  the smallest carrier among these that keeps the theorem's predicate.

For diagnostic Tournament Analysis on the twelve nonzero residue owners, the
pair observable is whether an order-two replacement at one owner covers the
opposite parity at the other.  A live pair is oriented provider-to-recipient;
silent pairs use increasing residue as the tie gauge.  Reversing the live-pair
switch gives

```text
live edges                         24
edge flips                         24
score histogram                    {2:2,4:3,5:1,6:1,7:3,9:2}
directed triangles                 40
SCC sizes                          [12]
Hamiltonian-path count             125,355.
```

One tie-gauge Hamiltonian path is

```text
10 -> 6 -> 7 -> 1 -> 12 -> 2 -> 11 -> 9 -> 5 -> 8 -> 3 -> 4.
```

The switched tournament has the same coarse fingerprint and a different
Hamiltonian path.  This cyclic shadow correctly advertises nontrivial residue
interaction, but the proof is the stronger labelled assertion that the two
opposite cross obligations cannot both be live.  The tournament alone forgets
that conjunction.

## Reproduction

```bash
python3 04-computation/lrc13_hamming_two_common_scale_rigidity_codex_S10.py
```

The replay uses only Python integers and `fractions.Fraction`; floating point
is absent.  It checks the complete order-two residue/parity table, the
complementary closed-endpoint guardrail, all `600,756` scale-one rows, the sharp
maximin row, and the tournament fingerprints.
