---
id: THM-2184
title: "Two-scale tail continuation profile and effective endpoint-grid stability"
status: >
  PROVED + VERIFIED-EXACT. Let a finite core E have all endpoints on the
  1/L grid and let the inserted tail be W_i=NLc_i+r_i. Its exact safe
  measure differs by at most 5||r||_1/(2NL) from a fixed rational two-torus
  continuation profile which retains the complete joint normalized-tail
  law. The profile is exactly decidable by a finite rational arrangement.
  Zero residue recovers THM-2182's exact product, proportional residue gives
  an effective near-aligned product theorem, and one tail recovers the
  qualitative O(1/W) scale of THM-2174. More strongly, every profile with at
  most six tail speeds has the pointwise floor 1-m/7. Hence every arbitrary-
  residue seven-core/six-tail affine ray closes effectively; the zero-profile
  residual for a positive-measure core can occur only with at least seven
  tail constraints and is not an LRC 7+6 obstruction. For seven-element
  cores inside {1,...,13}, THM-2166 gives the uniform near-aligned cone
  NL>=89||r||_1.
source: codex-2026-07-24-relation-carry-spectrum
depends_on:
  - THM-2182-endpoint-grid-product-and-tail-overlap-sidecar
  - THM-2166-hybrid-core-smoothing-low-carry-crossing
related:
  - THM-2162-signed-endpoint-cocycle-and-bv-component-split
  - THM-2174-endpoint-phase-scale-obstruction
script: 04-computation/lrc14_two_scale_tail_continuation_thm2184.py
output: 05-knowledge/results/lrc14_two_scale_tail_continuation_thm2184.out
script_sha256: 246d28ecdb9bc816091b3c1a76f6458bdd594794c1000cb72fa452d426cb730f
output_sha256: 28132bb89bc6fb836e8af3099d3afecc091137da8d4fac0f9462af275bc2706f
hash_basis: working-tree bytes (LF)
---

# THM-2184 -- two-scale tail continuation

At radius `1/14`, put

```text
G={x in R/Z:||x||>=1/14},
G_E={t:et in G for every e in E},
chi=1_G.                                              (1)
```

THM-2182 proves an exact product when a complete tail is scaled by an
endpoint-grid multiple. THM-2174 proves that a fixed nonzero residue retains
an `O(1/W)` endpoint current. The theorem below unifies these facts for an
arbitrary finite tail without discarding any of its internal overlaps.

## 1. The continuation profile and theorem

Let `E` be a finite set of positive integers, and choose `L>=1` such that

```text
14e divides L                  for every e in E.      (2)
```

The empty core is allowed. Fix positive integers

```text
c=(c_1,...,c_m)                                      (3)
```

and arbitrary integers `r=(r_1,...,r_m)`. For `N>=1`, put

```text
Q=NL,
W_i(N)=Qc_i+r_i,                                     (4)
```

and assume the displayed speeds are positive. Repetitions cause no problem
for the measure identity; LRC applications may separately require
distinctness.

Define the joint two-scale continuation profile

```text
Phi_(c,r)(t)
 =integral_(R/Z) product_(i=1)^m
      chi(c_i x+r_i t) dx,                           (5)

P_(E;c,r)
 =integral_(G_E) Phi_(c,r)(t)dt.                     (6)
```

> **Two-scale continuation theorem.** With
>
> ```text
> R=||r||_1=sum_i |r_i|,
> ```
>
> one has
>
> ```text
> |measure(G_(E union {W_i(N)}))
>       -P_(E;c,r)|
> <=5R/(2NL).                                         (7)
> ```

The target in (6) is finite and exact. On the two-torus with coordinates
`(t,x)`, its boundary arrangement consists of the rational lines

```text
et=+/-1/14 mod 1,
c_i x+r_i t=+/-1/14 mod 1.                           (8)
```

Their finitely many cells have rational vertices and rational areas.
Therefore `P_(E;c,r)` is rational and can be decided by exact arrangement
arithmetic. No limiting numerical oracle is hidden in (7).

## 2. Midpoint cell disintegration

Condition (2) makes `1_(G_E)` constant almost everywhere on every
`1/L` grid cell, hence on every `1/Q` cell. Let

```text
I_j=[j/Q,(j+1)/Q),
tau_j=(j+1/2)/Q,
A_j=the almost-everywhere value of 1_(G_E) on I_j.   (9)
```

Write `t=(j+u)/Q` on `I_j`. Since every `c_i j` is an integer,

```text
W_i(N)t
 =c_i u+r_i(j+u)/Q                 mod 1.             (10)
```

Thus the exact safe measure is

```text
M_N
 =1/Q sum_(j=0)^(Q-1) A_j
    integral_0^1 product_i
      chi(c_i u+r_i(j+u)/Q)du.                       (11)
```

Freeze only the slow variable at the cell midpoint:

```text
S_N=1/Q sum_j A_j Phi_(c,r)(tau_j).                  (12)
```

For a fixed factor in (11), the frozen phase is

```text
y_i(u)=c_i u+r_i tau_j
```

and the omitted displacement is

```text
r_i(u-1/2)/Q,
```

of magnitude at most `|r_i|/(2Q)`. If the two interval indicators differ,
`y_i(u)` lies within that distance of one of the two endpoints of `G`.
Because `c_i!=0`, the map `u|->c_i u+r_i tau_j` preserves Haar measure.
The two endpoint neighborhoods therefore have total measure at most

```text
2|r_i|/Q.                                            (13)
```

Product telescoping and (13) give, on each cell,

```text
integral_0^1
 |product_i chi(c_i u+r_i(j+u)/Q)
  -product_i chi(c_i u+r_i tau_j)|du
 <=2R/Q.                                              (14)
```

After the outer factor `1/Q` and the sum over at most `Q` cells,

```text
|M_N-S_N|<=2R/Q.                                     (15)
```

It remains to compare the midpoint sum with the profile integral. A circular
interval and its translate by `h` have symmetric-difference measure at most
`2||h||`. Another product telescope therefore gives

```text
|Phi_(c,r)(t)-Phi_(c,r)(s)|
 <=2R|t-s|                                           (16)
```

whenever `t,s` are represented in one grid cell. Midpoint quadrature on a
cell now costs at most

```text
2R integral_(-1/(2Q))^(1/(2Q)) |z|dz
 =R/(2Q^2).                                          (17)
```

Since `G_E` is a union of selected cells up to endpoints,

```text
|S_N-P_(E;c,r)|<=R/(2Q).                             (18)
```

Equations (15) and (18), with `Q=NL`, prove (7). QED.

## 3. Exact specializations and effective exits

### Zero residue

If `r=0`, then

```text
Phi_(c,0)(t)=measure(G_C),
C={c_1,...,c_m},                                     (19)
```

and the error in (7) is zero. Hence

```text
measure(G_(E union NLC))
 =measure(G_E)measure(G_C),                          (20)
```

which recovers THM-2182.

### Proportional residue

If `r_i=k c_i` for one integer `k`, translation `x|->x+kt` in (5) again
gives the constant profile (19). Therefore

```text
|measure(G_(E union (NL+k)C))
   -measure(G_E)measure(G_C)|
 <=5|k| sum_i c_i/(2NL).                             (21)
```

This is an effective near-aligned endpoint-grid theorem. In the `7+6`
case, if `E={e_1<...<e_7}`, the elementary bounds used in THM-2182 give

```text
measure(G_E)measure(G_C)>=1/(49e_2).                 (22)
```

Consequently every positive, distinct row `E union (NL+k)C` is strictly
safe once

```text
N>245 e_2 |k| (sum_i c_i)/(2L).                      (23)
```

The aligned case `k=0` is exact for every `N`.

### Arbitrary residue

For every fixed `t`, each individual danger set in the `x` fiber of (5) has
Haar measure `1/7`: multiplication by the positive integer `c_i`, followed by
translation by `r_i t`, preserves Haar measure. The union bound therefore
gives the pointwise profile floor

```text
Phi_(c,r)(t)>=1-m/7                    whenever m<=7. (24a)
```

In particular, when `m=6` and `E={e_1<...<e_7}`, (13) of THM-2182 gives

```text
P_(E;c,r)>=measure(G_E)/7>=1/(49e_2).                (24b)
```

This bound is independent of the normalized tail and of its residues. Thus
every positive, distinct affine row

```text
E union {NLc_i+r_i:1<=i<=6}
```

is strictly safe once

```text
N>245e_2||r||_1/(2L).                                (24c)
```

This is the arbitrary-residue strengthening of (23), not merely a
proportional-tail result.

For the literal defect-six core bank

```text
E subset {1,...,13},                 |E|=7,
```

THM-2166's exact sweep sharpens the crude second-speed floor to

```text
measure(G_E)>=45107/229320.                           (24d)
```

For six tails, (24a) and (6) therefore give

```text
P_(E;c,r)>=45107/1605240.                             (24e)
```

Combining this with (7), every positive distinct row on such a ray is
strictly safe whenever, with `Q=NL` and `R=||r||_1`,

```text
Q>(4013100/45107)R.                                   (24f)
```

Since `4013100/45107<89`, the clean integer cone

```text
Q>=89R                                                (24g)
```

is sufficient when `R>=1`; for `R=0` the exact product is already
positive. This is uniform over all `1716` seven-element cores in the bank.

More generally, for any tail length, if

```text
P=P_(E;c,r)>0,
```

then (7) proves strict safety for every

```text
N>5||r||_1/(2LP).                                    (24)
```

Thus every fixed core, normalized tail, and residue vector has an effective
finite terminal unless its rational two-torus profile is zero. The zero
profile is the honest residual: a fixed rank-at-most-two strip arrangement,
not an unbounded speed family. When `measure(G_E)>0`, (24a) shows that it
requires at least seven tail constraints and therefore cannot occur in the
LRC `7+6` application.

For `m=1`, Haar invariance makes `Phi=6/7`; (7) supplies the same qualitative
`O(1/W)` scale as THM-2174. THM-2174's signed endpoint formula is much sharper
in that specialization and detects the exact numerator `C_r`.

## 4. Exact controls and transfer boundary

For a genuinely nonconstant profile, take

```text
E=empty,       L=14,
c=(1,1),       r=(-1,1).                             (25)
```

The character map

```text
(t,x)|->(x-t,x+t)
```

is a surjective torus endomorphism. Moreover

```text
Phi_(c,r)(0)=6/7,
```

whereas Haar pushforward under that endomorphism gives

```text
P_(empty;c,r)=(6/7)^2=36/49.                         (26)
```

Thus the profile is genuinely nonconstant, while its average is still exact.
The actual pairs are `(14N-1,14N+1)`. Exact interval arithmetic at
`N=1,2,5,10` gives respective deviations

```text
-2/3185, -2/12789, -2/80017, -2/320117,             (27)
```

all within (7).

For a nonempty proportional control, take

```text
E=(1,2,3,4,5,6,8),       L=1680,
C=(3,12),                k=1.                        (28)
```

Here

```text
measure(G_E)=27/70,
measure(G_C)=3/4,
P=81/280.                                             (29)
```

At `N=1,2,5`, the exact deviations of
`measure(G_(E union (NL+1)C))` from `P` are

```text
47/1412040, 47/2823240, 47/7056840,                 (30)
```

again inside (21). With `k=0`, the exact measure is `81/280`, as (20)
requires. The companion computes all one-dimensional measures by rational
danger-interval union and independently checks the displayed arithmetic.

The typed ledger is:

```text
source:       endpoint-aligned core plus affine tail ray;
map:          midpoint disintegration on NL cells;
preserved:    full joint normalized-tail law Phi;
lost:         O(||r||_1/(NL)) slow phase drift;
finite exit:  exact rational two-torus profile P;
residual:     only P=0, including possible measure-zero weak witnesses.     (31)
```

This theorem does not prove that every LRC row lies on one affine tail ray,
does not classify the zero-profile two-torus arrangements, and does not by
itself prove LRC(14). QED.
