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
  residue seven-core/six-tail affine ray, and an explicit l1 tube around the
  entire endpoint grid, closes effectively. The normalized tail and residue
  may vary inside that tube. The zero-profile residual for a positive-measure
  core can occur only with at least seven tail constraints and is not an LRC
  7+6 obstruction. For seven-element cores inside {1,...,13}, THM-2166 gives
  the stronger uniform cone NL>=89||r||_1. The proof works unchanged for any
  rational-grid measurable core. On THM-2168's scalar 5+3 residual this gives
  a uniform three-tail tube and a formally valid explicit congruence
  criterion. Combining the all-depth invoice with the full six-coefficient
  endpoint grid and THM-2192's unique-deepest law shows that this particular
  canonical full-grid instance never fires; compressed residual grids are
  not excluded. This is a route no-go, not a closure of scalar 5+3.
source: codex-2026-07-24-relation-carry-spectrum
depends_on:
  - THM-1166-seven-wall-fano-gcd-discrepancy
  - THM-2137-deep-scalar-tail-boundary-complexity
  - THM-2166-hybrid-core-smoothing-low-carry-crossing
  - THM-2168-three-target-second-depth-majorization
  - THM-2182-endpoint-grid-product-and-tail-overlap-sidecar
  - THM-2192-scalar-five-plus-three-root-sheet-chord-invoice
related:
  - THM-2162-signed-endpoint-cocycle-and-bv-component-split
  - THM-2174-endpoint-phase-scale-obstruction
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
script: 04-computation/lrc14_two_scale_tail_continuation_thm2184.py
output: 05-knowledge/results/lrc14_two_scale_tail_continuation_thm2184.out
script_sha256: 246d28ecdb9bc816091b3c1a76f6458bdd594794c1000cb72fa452d426cb730f
output_sha256: 28132bb89bc6fb836e8af3099d3afecc091137da8d4fac0f9462af275bc2706f
independent_script: 04-computation/lrc14_scalar_endpoint_grid_no_go_thm2184.py
independent_output: 05-knowledge/results/lrc14_scalar_endpoint_grid_no_go_thm2184.out
independent_script_sha256: 43d914f9214ab3f9e9f5f979da2d79b88238b826efb7346376995d075d139aff
independent_output_sha256: 417b2288811fc9bd9c62bc36509d9e937fca7961722fe026c020f6afde91ea58
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

In fact no fixed-ray hypothesis is needed for the quantitative tube. Combining
(7) and (24b) gives, uniformly in the positive integer vector `c` and the
integer vector `r`,

```text
measure(G_(E union {NLc_i+r_i}))
 >=1/(49e_2)-5||r||_1/(2NL).                         (24h)
```

Consequently the row is strictly safe whenever

```text
||r||_1<2NL/(245e_2).                                (24i)
```

Both `c` and `r` may vary with `N` here. In particular, every sequence of
positive, distinct six-tails

```text
W_i(N)=NLc_i(N)+r_i(N),       ||r(N)||_1=o(N),
```

is eventually safe beside the fixed seven-core. This is an explicit
endpoint-grid tube, not just a union of fixed affine rays. The theorem does
not assert that an arbitrary tail admits such a simultaneous approximation.

### Rational-grid cores and the scalar `5+3` terminal

Nothing in the proof of (7) uses the special form of `G_E` beyond the
cell-constancy in (9). Let `A` be any measurable subset of the circle whose
indicator is constant almost everywhere on each `1/J` cell, and define

```text
P_(A;c,r)=integral_A Phi_(c,r)(t)dt.                 (GC.1)
```

Repeating Sections 1--2 with `1_A` in place of `1_(G_E)` gives

```text
|measure(A intersection intersection_i {t:W_i(N)t in G})
      -P_(A;c,r)|
 <=5||r||_1/(2NJ),       W_i(N)=NJc_i+r_i.           (GC.2)
```

Thus the continuation theorem is really a rational-grid core theorem; an
ordinary lonely set is only its main application.

Apply this to the scalar survivor in Section 9 of THM-2168. In that notation
put

```text
A=C_H minus union_(i=1)^5 D_(q_i),
C_H={t:||Ht||>1/7},                                  (SC.1)
```

and choose `J` divisible by `7H` and every `14q_i`. Then `1_A` is
`1/J`-cell-constant almost everywhere, while THM-2137 gives

```text
measure(A)>=961/6930.                                (SC.2)
```

For any three-tail profile, the pointwise union bound in the `x` fibre gives
`Phi>=4/7`. Hence

```text
P_(A;c,r)>=1922/24255,                               (SC.3)

measure(A intersection intersection_(j=1)^3
        {t:(NJc_j+r_j)t in G})>0
whenever ||r||_1<3844NJ/121275.                      (SC.4)
```

There is a sharper terminal for the original three blockers. Write them as

```text
13s_j=Mw_j,             M=13^d,
M=nJ+k,                 n>=1,       k in Z,          (SC.5)
```

where `d` is their minimum common 13-adic depth and the `w_j` are distinct.
Taking `c_j=w_j` and `r_j=kw_j` in (GC.2) recovers the original
speeds exactly. The residue is proportional, so translation in `x` makes
the profile constant. The sharp three-comb theorem of THM-1166 supplies

```text
Phi_(w,kw)=measure(G_{w_1,w_2,w_3})>=55/91,
P_(A;w,kw)>=961/11466.                               (SC.6)
```

Consequently the scalar `5+3` covering branch is impossible whenever

```text
|k| sum_(j=1)^3 w_j < 961nJ/28665.                   (SC.7)
```

This is an intrinsically small-quotient obstruction. Indeed, `J` is even
while `M=13^d` is odd, so `k` is always nonzero. If
`X=max(H,q_1,...,q_5)`, then `X>=B/6`, `J>=7X`, and the THM-2168 invoice
`B>=(12493/35640)M` gives

```text
M/J <=213840/87451<2.446.                            (SC.8)
```

Since `sum_j w_j>=6`, for `n>=3` the left side of (SC.7) is at least
`6(n-213840/87451)J`, which is greater than its right side. Thus only the
neighboring grid scales `n=1,2` remain as apparent candidates under this
weak maximum-coordinate estimate; choosing `n=1` also represents `M<J`
through a negative residue.

There is, however, a decisive obstruction to using (SC.7) with the
**canonical full six-coefficient endpoint grid** of the actual THM-2168
residual.  The guard and five unit coefficients there are pairwise distinct.
Write

```text
J=14Q.
```

The grid divisibilities imply `q_i|Q`.  They also imply `H|2Q`, and the
oddness of `H` improves this to `H|Q`.  Thus

```text
a_0=Q/H,                 a_i=Q/q_i,  i=1,...,5,
```

are six distinct positive integers.  Ordering them and using the sixth
harmonic number gives the sharp divisor envelope

```text
B/Q=sum_(i=0)^5 1/a_i
   <=1+1/2+1/3+1/4+1/5+1/6
    =49/20,
B<=7J/40.                                             (SC.9)
```

This envelope is arithmetically sharp even with an odd guard and
thirteen-unit coefficients: take `Q=60`, `H=15`, and
`(q_i)=(60,30,20,12,10)`.  Combining (SC.9) with THM-2168's all-depth
invoice yields

```text
J/M >=(40/7)(12493/35640)
     =12493/6237
     =2+19/6237,
M/J <=6237/12493<1/2.                                (SC.10)
```

THM-2192 makes the normalized-tail obstruction still farther away.  Order
the actual blocker valuations as

```text
lambda_1<=lambda_2<lambda_3,          d=lambda_1.
```

If `lambda_2=d`, then two distinct `w_j` are thirteen-units and the third
is divisible by thirteen, so

```text
sum_j w_j>=1+2+13=16.
```

If `lambda_2>d`, the last two normalized coefficients are divisible by
`13` and `13^2`, respectively, and the stronger floor is

```text
sum_j w_j>=1+13+169=183.                             (SC.11)
```

In particular, (SC.10) makes `k=M-nJ` negative for every `n>=1`, and the
two valuation cases above give

```text
|k| sum_j w_j/J
 >=16(n-6237/12493).
```

The difference between this lower bound and the right side of (SC.7),
after division by `J`, is increasing in `n`; already at `n=1` it is

```text
16(1-6237/12493)-961/28665
 =219788159/27547065>0.                              (SC.12)
```

Hence neither `n=1` nor `n=2` -- indeed, no `n>=1` -- can satisfy (SC.7)
for any `J` satisfying the full-grid divisibilities chosen after (SC.1).
Moreover `J` is a multiple of fourteen, so

```text
k=13^d-nJ ==(-1)^d (mod 14);
```

`k` is odd, and since it is negative, `|k|` is `13` modulo fourteen for
even `d` and `1` modulo fourteen for odd `d`.  The signed quotient remainder
is therefore fully controlled, but in the direction opposite to the
continuation tube.

The criterion (SC.7) remains formally valid for an abstract rational-grid
core.  What (SC.9)--(SC.12) prove is a precise no-go for this canonical
full-grid application: a grid required to resolve each of the six shallow
coefficients is already more than twice the minimum deep scale.  They do not
exclude a smaller grid arising from boundary masking or cancellation in the
actual residual set.  The remaining scalar target must therefore use either
such compressed endpoint information or data discarded by the proportional
continuation profile, such as THM-2192's root-sheet ownership, signed overlap,
or matching-support constraints.  This does not exclude the scalar branch.
The independent companion checks the sharp divisor row, both valuation
minima, every rational constant in (SC.10)--(SC.12), the monotone `n=1`
endpoint, and the signed mod-fourteen remainder law with checks active under
optimized Python.

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

THM-735 is the essential prior comparison. Its simultaneous multi-peel
already proves that every fixed positive-measure body with at most six
sufficiently remote tail speeds is eventually safe, with no alignment or
separation required among the tail speeds. THM-2184 does not supersede that
far-cone theorem. It instead gives the exact joint continuation law, a
component-count-free grid error, and the explicit endpoint-grid tube (24i).
Thus the two results are complementary descriptions of the same terminal
region: THM-735 is alignment-free at sufficient distance, while THM-2184 is
phase-sensitive and exact near the endpoint grid.

The typed ledger is:

```text
source:       endpoint-aligned core plus grid-decomposed tail;
map:          midpoint disintegration on NL cells;
preserved:    full joint normalized-tail law Phi;
lost:         O(||r||_1/(NL)) slow phase drift;
finite exit:  exact rational two-torus profile P;
residual:     only P=0, including possible measure-zero weak witnesses.     (31)
```

This theorem does not force every LRC row into the THM-735 far cone or the
endpoint-grid tube, does not classify the zero-profile two-torus
arrangements, and does not by itself prove LRC(14). QED.
