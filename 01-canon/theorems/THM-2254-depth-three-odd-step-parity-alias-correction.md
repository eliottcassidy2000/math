---
id: THM-2254
title: "Depth-three odd-step parity-alias correction"
status: >
  PROVED + VERIFIED-EXACT + NEGATIVE AUDIT / CORRECTION. For the
  three-state annular partition underlying the all-equal (3,4,5) branch,
  exact 13-root counting gives P=(-I+14 Pi)/13: both centered directions
  have eigenvalue -1/13. Swapping the two equal-mass inner states gives a
  second stochastic matrix JP with the same 169-step square but the wrong
  sign on the antisymmetric state. The formerly reported u=H and u=6H
  whole-clause closures used JP. Under the correct P their exact bounds are
  28460/199927 and 734515/5198102, both above 961/6930. The correction
  supplies the exact lost midpoint Z/2 bridge coordinate but no exclusion:
  the all-equal branch and LRC(14) remain open.
source: codex-sol-2026-07-25-depth-three-parity-alias-audit
depends_on:
  - THM-2250-depth-three-pair-incidence-partition-reduction
related:
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
  - THM-2232-same-core-signed-eigen-markov-dual-exclusion
  - THM-2243-composite-union-transfer-dual-depth-three-profile-exclusion
script: 04-computation/lrc14_depth_three_parity_alias_correction_thm2254.py
output: 05-knowledge/results/lrc14_depth_three_parity_alias_correction_thm2254.out
script_sha256: 315144dc6abc55e030f7181cf4acf5b8bdc9212f49ca5e3394853bcdc5566409
output_sha256: a823391c809438e0afae6e221ba96d5234d27aa3409dcab8673cf7f0ddf2ffd6
hash_basis: working-tree bytes (LF)
---

# THM-2254 -- the odd-step parity alias

THM-2250 leaves one depth-three scalar branch:

```text
(lambda_0,lambda_1,lambda_2)=(3,4,5),
u_0=u_1=u_2=u.                                      (1)
```

A subsequent working calculation appeared to close the hostile relations
`u=H` and `u=6H` by recovering a unit-13 transition from the already known
unit-169 transition. That recovery chose the wrong stochastic square root.
This theorem records the exact correction, identifies the coordinate lost
under squaring, and evaluates both controls with the correct root.

The support reduction and unequal-core exclusions in THM-2250 are
unaffected. Only the attempted all-equal continuation is retracted.

## 1. Exact three-state process

On `T=R/Z`, put

```text
D_1={z:||z||<1/14},       E_1={z:||z||<1/7},

A=D_1,       B=E_1\D_1,       C=E_1^c.              (2)
```

Null endpoints may be assigned arbitrarily. For Haar-uniform `x`, define

```text
z_t=13^t u x,       Z_t=state(z_t) in {A,B,C}.       (3)
```

The stationary law is

```text
pi=(1,1,5)/7.                                        (4)
```

Every boundary of the partition in (2), and every pullback of such a
boundary under multiplication by `13`, lies on the `1/182` grid. One
midpoint in each open grid cell therefore gives an exact finite
enumeration. The current/next cell counts are

```text
             next A   next B   next C
current A       2        4        20
current B       4        2        20
current C      20       20        90.               (5)
```

Dividing each row by its mass gives the forward transition

```text
       [ 1   2  10 ]
P=1/13[ 2   1  10 ].                                (6)
       [ 2   2   9 ]
```

This is an exact word law, not a one-step approximation. For a point in
future state `A`, `B`, or `C`, respectively, its thirteen roots have state
counts

```text
(1,2,10),       (2,1,10),       (2,2,9).             (7)
```

The counts are pointwise constant on each future state. Backward root
induction therefore proves, for every finite word,

```text
Pr(Z_0=i_0,...,Z_m=i_m)
 =pi_(i_0) product_(t<m) P_(i_t,i_(t+1)).            (8)
```

In particular `pi P=pi`, and the chain is reversible.

## 2. The two square roots and the lost sign

Let `Pi` be the rank-one stochastic projector whose every row is `pi`, and
put `rho=-1/13`. Directly from (6),

```text
P=rho I+(1-rho)Pi=(-I+14Pi)/13.                     (9)
```

Thus every centered observable, not only one selected parity observable,
lies in the same two-dimensional eigenspace:

```text
pi f=0       implies       Pf=-(1/13)f.             (10)
```

Let `J` exchange `A` and `B` and fix `C`. Since `pi` is `J`-invariant,
`J` commutes with `P`. The parity-twisted root

```text
        [ 2   1  10 ]
Q=JP=1/13[ 1   2  10 ]                              (11)
        [ 2   2   9 ]
```

is also stochastic. Nevertheless

```text
          [ 25  24  120 ]
P^2=Q^2=1/169[ 24  25  120 ].                       (12)
          [ 24  24  121 ]
```

For the centered basis

```text
phi=(1,-1,0),       psi=(5,5,-2),                   (13)
```

the one-step actions are

```text
P phi=-(1/13)phi,       Q phi=+(1/13)phi,
P psi=-(1/13)psi,       Q psi=-(1/13)psi.           (14)
```

Squaring erases exactly the sign of the antisymmetric `A-B` channel. This
is the missing `Z/2` sheet: unit-169 endpoint data alone cannot decide
between `P` and `Q`.

The already proved binary danger law selects `P` immediately. Lumping
`A` against `B union C`, the correct transition has

```text
Pr(A_next | A_current)=1/13,
Pr(A_next | current not A)=2/13.                    (15)
```

The twisted root instead gives `2/13` and `11/78`. Thus `Q` conflicts with
the exact one-core law before any carrier optimization is performed.

## 3. The exact midpoint bridge retained by the Z/2 sidecar

The coordinate lost by the 169-step quotient can be written without
metaphor. Given even-time endpoints `Z_0=i`, `Z_2=j`, the exact midpoint
bridge is

```text
K_(ij)(k)=P_(ik)P_(kj)/(P^2)_(ij).                  (16)
```

Its signed `A-B` expectation is

```text
                 endpoint j
             A       B       C
endpoint A  -3/25     0     -1/12
i        B    0      3/25    1/12
         C  -1/12    1/12      0.                  (17)
```

Replacing `P` by `Q` negates every entry of (17), while leaving the
endpoint kernel (12) unchanged. The corresponding correct conditional
probability that the midpoint is dangerous is

```text
                 endpoint j
             A       B       C
endpoint A   1/25   1/12    1/12
i        B   1/12   4/25    1/6
         C   1/12   1/6    20/121.                 (18)
```

Equations (17)--(18) are a concrete refinement of every 169-fibre count
that uses an odd-time danger atom: retain the midpoint annular label, or
equivalently its antisymmetric sign `phi`. They also mark the limit of this
refinement. The exact `P`-word evaluation below already retains the full
three-state midpoint, so the same `Z/2` sidecar cannot sharpen it again.

## 4. The all-equal whole-clause cylinder

Write

```text
X_t=1_{Z_t=A},
W_t=X_t OR X_(t+1) OR X_(t+2).                      (19)
```

For the all-equal profile (1), THM-2250 proves

```text
support(R_+) subset {W_1=W_3=1}
```

and, using the centered nonnegative charge

```text
C=W_2+(1/169)(1-W_0),
```

reduces the positive residual mass `p` to

```text
p<=B_P(g)
  :=E_P[g(Z_0)W_1W_3{W_2+(1/169)(1-W_0)}],          (20)
```

whenever the full phase-fibre conditional expectation satisfies

```text
E[1_(C_H) | z_0]=g(state(z_0))                      (20a)
```

for the relation under study. Merely knowing the conditional expectation
given the coarse state would not suffice, because the future word is a
function of the full phase `z_0`.

There are two exact hostile controls:

```text
u=H:       g=(0,0,1);

u=6H:      g=(5/6,5/6,4/6).                         (21)
```

For the second line, condition on `z=ux=6Hx`. Among the six roots for
`Hx`, five lie in `C_H` when `z` is in `A` or `B`, and four lie in `C_H`
when `z` is in `C`. The count changes only at the state boundaries, so
(21) is pointwise exact off a null set.

## 5. Corrected hostile values

Enumerating the six-state words in (20) using (8) gives

```text
B_P(u=H)
 =28460/199927
 =961/6930 + 728279/197927730
 >961/6930,                                         (22)

B_P(u=6H)
 =734515/5198102
 =961/6930 + 3386176/1286530245
 >961/6930.                                         (23)
```

Neither correct bound closes its hostile relation.

The formerly reported values are exactly the evaluations under the
spurious root `Q`:

```text
B_Q(u=H)
 =354710/2599051
 =961/6930 - 5649673/2573060490,

B_Q(u=6H)
 =1057220/7797153
 =961/6930 - 7929973/2573060490.                    (24)
```

The corrections from (24) to (22)--(23) are

```text
15270/2599051,       89105/15594306,                (25)
```

respectively. Both cross the target in the wrong direction for exclusion.
The session-level claim that these two relations were closed is therefore
retracted.

For a further hostile control, deleting the guard entirely gives

```text
B_P(1)=514705/2599051.
```

Even the diagnostic decorrelated replacement by the guard mean `5/7`
gives

```text
(5/7)B_P(1)
 =2573525/18193357
 =961/6930 + 50101739/18011423430.                  (26)
```

Thus an independence heuristic cannot repair this cylinder estimate.
Equation (26) is a stopping control, not an independence assertion about
the arithmetic process.

## 6. Preserved/lost-coordinate ledger and next interface

```text
source:
  the exact unit-13 roots of the annular three-state partition;

169-step quotient preserves:
  the stationary law, both even-time endpoints, P^2, and the magnitude
  1/169 of every centered two-step eigenvalue;

169-step quotient destroys:
  the sign of the antisymmetric A-B midpoint channel;

repair map:
  retain phi(Z_(2r+1)), or the full bridge K in (16), between consecutive
  169 checkpoints;

predicate sharpened:
  every cylinder using an odd-time danger atom X_(2r+1);

limit:
  the exact P-word law already contains this bridge and still gives
  (22)--(23), so no second use of the same parity sheet can close them.
                                                               (27)
```

The next information missing from (20) is not another square-root choice.
The exact residual has

```text
R_+
 =1_(C_H) product_(i=1)^5 (1-1_(D_(q_i))).          (28)
```

Inequality (20) discarded all five unit-danger punctures in (28). A
strict all-equal improvement must therefore retain some interaction of
those punctures with the whole-clause carrier, or introduce another
coordinate not measurable by the exact three-state word. No such
improvement is claimed here.

## 7. Exact audit and scope

Run

```bash
python3 04-computation/lrc14_depth_three_parity_alias_correction_thm2254.py
python3 -O 04-computation/lrc14_depth_three_parity_alias_correction_thm2254.py
```

The companion uses only exact rational arithmetic. It verifies:

```text
the 182-cell forward count and pointwise thirteen-root reverse law;
stationarity, reversibility, and the rank-one spectral decomposition;
P^2=(JP)^2 and both centered one-step eigen-signs;
the exact midpoint signed and danger bridge tables;
the binary danger-lump hostile check;
the u=H and u=6H full-phase fibre guard weights;
all correct, twisted, unguarded, and decorrelated cylinder values;
an independent 5,198,102-cell direct integration of the three correct
cylinders, agreeing with the Markov-word evaluation;
ordinary/optimized transcript identity.             (29)
```

This is a correction and a negative structural theorem. It excludes no
new coefficient relation or valuation profile. After THM-2250, the
all-equal `(3,4,5)` branch remains open, and LRC(14) remains open.
