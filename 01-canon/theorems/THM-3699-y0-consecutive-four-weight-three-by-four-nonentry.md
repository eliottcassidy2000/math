---
id: THM-3699
title: "Consecutive four-weight three-by-four nonentry and a two-weight adjunction gate"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the y=0
  collision ring, no Darboux pair can have reduced grading supports of sizes
  3 and 4 when both supports lie in any common four-consecutive-weight window.
  The synchronized arm law forces that window to be {-2,-1,0,1}.  Three
  missing-weight cases there collapse on a singleton bucket.  In the fourth,
  same-sign endpoint commutation
  upgrades the weight -2 coefficient to carry two arm factors, so every
  scalar address vanishes on the arm divisor.  Consequently the seven-space
  of THM-3698 remains Darboux-empty after adjoining one homogeneous weight
  -1 function and one homogeneous weight-zero function.  Arbitrary 3x4
  support words, higher-dimensional adjunctions, and JC(2) remain OPEN.
source: jc-quartic-c3-construction / 2026-08-22
audit: >
  PASS.  The independent hostile audit checked all four bucket tables, both
  zero/nonzero singleton derivative gates, the strict-opposite product
  derivative, the endpoint logarithmic derivative and common-square upgrade,
  the common arm factor in all three scalar addresses, the arm-law reduction
  from every translated consecutive window to the unique window containing
  both -2 and 1, and the exact scope of the one-dimensional Darboux shear.
  Normal and optimized companion runs byte-match the stored transcript.
depends_on:
  - THM-3695-y0-collision-ring-danielewski-embedding-and-seven-piece-floor
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
  - THM-3698-y0-collision-seven-function-pluecker-compression-no-go
related:
  - THM-3592-universal-exponent-two-three-by-three-weight-darboux-nonentry
  - THM-3693-y0-collision-ring-two-by-three-weight-darboux-no-go
script: 04-computation/jacobian_y0_consecutive_four_weight_three_by_four_thm3699.py
output: 05-knowledge/results/jacobian_y0_consecutive_four_weight_three_by_four_thm3699.out
script_sha256: 645799c066f84ec768a53ee5bd912ed16715ee15e5df76d58f0f6e05edec9e21
output_sha256: b6a2ad08cd285ab929c27b93223ef98c7f5f7ff8a9e4de93f20756d66497219a
hash_basis: LF-normalized bytes
---

# THM-3699 -- the conductor-natural consecutive `3 x 4` word is empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This closes one
natural all-degree slice of the first support cell left by THM-3695.  The
proof is structural rather than a bounded coefficient search.

All rings are over `C`.  Retain the graded embedding of THM-3695/3696

```text
h=1-b^2,
D=C[b,c,e]/(c^2e-h),                  wt(b,c,e)=(0,1,-2),
R=C[e,ce,bc] subset D.                                      (1)
```

For a homogeneous element write

```text
F_r=c^r f(b).
```

The ambient Danielewski bracket is

```text
{c^r f,c^s g}=c^(r+s+1)[s f'g-r fg'].                    (2)
```

The coefficient modules needed below are the following consequences of the
complete table in THM-3696:

```text
M_1=bR_0,
M_0=R_0,
M_-1 subset hC[b],
M_-2 subset hC[b].                                      (3)
```

Here a retained weight-zero piece is nonconstant.  Scalar weight-zero pieces
are removed before support is counted, since they are bracket-invisible.

## 1. Every translated consecutive window reduces to one case

For `t in Z`, put

```text
W_t={t,t+1,t+2,t+3}.                                   (4a)
```

Suppose that the reduced supports have sizes three and four and are both
contained in `W_t`.  At either arm `b=+-1`, THM-3696 says every scalar
address vanishes except one with weights `(-2,1)` or `(1,-2)`.  Since the
full bracket is the unit, such an address must occur.  Thus `W_t` contains
both `-2` and `1`.  They are distance three apart, so the only
four-consecutive integer window containing both is

```text
t=-2,                         W_t={-2,-1,0,1}.          (4b)
```

It remains to close this forced window.

## 2. Forced-window bucket grammar

Let `P,Q in R` satisfy `{P,Q}=1`, and suppose that, after possibly exchanging
the outputs,

```text
supp(P)=S,                  supp(Q)=W,
W={-2,-1,0,1},              S subset W, |S|=3.          (4)
```

A pair of weights `(r,s)` contributes to output weight `r+s+1`.  Every
nonzero output-weight bucket must therefore vanish, while the weight-zero
bucket must equal one.  There are exactly four possibilities for `S`.

## 3. Three singleton collapses

### 3.1 `S={-2,-1,0}`

Output weight `2` is the singleton `(0,1)`, so

```text
{P_0,Q_1}=c^2 p_0' q_1=0.                              (5)
```

The active coefficient `q_1` is nonzero.  Hence `p_0'=0`, making `P_0`
scalar and therefore not retained.  This contradicts `(4)`.

### 3.2 `S={-2,-1,1}`

Output weight `2` is now the singleton `(1,0)`, and

```text
{P_1,Q_0}=-c^2 p_1 q_0'=0.                             (6)
```

Thus the retained `Q_0` is scalar, again a contradiction.

### 3.3 `S={-2,0,1}`

The output-weight-one bucket is

```text
{P_0,Q_0}+{P_1,Q_-1}=0.                                (7)
```

The first bracket vanishes identically because both entries are functions of
`b`.  Formula `(2)` turns the second equation into

```text
0={c p,c^-1 q}=-c(pq)'.                                (8)
```

Thus `pq` is constant.  But `(3)` gives `b|p` and `h|q`; the nonzero product
is divisible by the nonconstant polynomial `bh`.  This is impossible.

## 4. Endpoint commutation erases the last scalar address

It remains to take

```text
S={-1,0,1}.                                             (9)
```

The lowest output bucket is the singleton `(-1,-2)`.  If its coefficient
polynomials are `p,q`, equation `(2)` says

```text
0=-2p'q+pq'=p^3(q/p^2)'.                              (10)
```

Hence `q=lambda p^2` in `C(b)` for a nonzero scalar `lambda`.  Since `p,q`
are nonzero polynomials, this is a polynomial identity.  Membership
`p in M_-1` gives `h|p`, and therefore

```text
h^2|q.                                                  (11)
```

There are precisely three scalar addresses:

```text
(-1,0),                    (0,-1),                    (1,-2). (12)
```

The first is divisible by `h` because its weight `-1` coefficient is; the
second is divisible by `h` for the same reason on the other output.  For the
last address, write `q=h^2t`.  Formula `(2)` gives

```text
{c p_1,c^-2 q}=c^0[-2p_1'q-p_1q'],                    (13)
```

and both `q` and `q'` are divisible by `h`.  Thus every summand in the entire
scalar bucket is divisible by `h`.  Their sum cannot be the unit `1`.

Sections 2--3 prove

```text
supp(P),supp(Q) subset W_t for some t in Z,
{|supp(P)|,|supp(Q)|}={3,4}
                 ==> {P,Q} notin C*.                  (14)
```

The same proof applies after exchanging the outputs; reversing the bracket
only changes its sign.

## 5. The first two-new-weight enlargement of the seven-space

Let `E` be the seven-dimensional space of THM-3698.  Its weight profile is

```text
dim E_-2=4,                 dim E_1=3,                 (15)
```

with no other weights.  Choose arbitrary homogeneous functions

```text
F_-1 in R_-1,              F_0 in R_0,                (16)
```

and put

```text
E^+=E+span_C(F_-1,F_0).                                (17)
```

Then `E^+` contains no Darboux pair.  Indeed, THM-3695 forces at least three
retained weights in each output and at least seven in total.  Since every
element of `E^+` is supported in `W`, the only possible profiles are

```text
3 x 4,                    4 x 3,                    4 x 4. (18)
```

The first two are excluded by `(14)`.  In a `4 x 4` pair, the two weight
`-1` components are nonzero multiples of the single vector `F_-1`.  A
Darboux shear

```text
(P,Q) -> (P,Q-lambda P)                                (19)
```

cancels the weight `-1` component of the second output while preserving
`{P,Q}=1`.  The new pair has at most `4 x 3` support.  THM-3695 rules out a
drop below three pieces, and `(14)` rules out exactly `4 x 3`.

Only the one-dimensionality of one new weight space is needed for the shear;
the symmetric weight-zero shear is a second check.  This corollary is about
the particular common span `(17)`.  It does **not** exclude a general `4 x 4`
pair whose components at weights `-1` and `0` may be linearly independent.

## 6. Exact frontier

The result closes every common four-consecutive-weight word, with the
conductor-natural window `W_{-2}` supplying the only nontrivial case.  It also
closes the first nine-dimensional enlargement suggested by THM-3698.  It does not close
arbitrary `3 x 4` placements with other weights.  By the repaired
every-line/source-jet argument in THM-3695, `2 x 5` is not a live collision-
ring cell: the unique first cell is `3 x 4` up to exchanging outputs.  The
next useful search is therefore another `3 x 4` support word, with the
synchronized three-branch scalar law of THM-3696 retained as a sidecar.

## Reproduction

```bash
python3 -B 04-computation/jacobian_y0_consecutive_four_weight_three_by_four_thm3699.py
python3 -B -O 04-computation/jacobian_y0_consecutive_four_weight_three_by_four_thm3699.py
```

Both commands must agree byte for byte with the frozen transcript.  **QED.**
