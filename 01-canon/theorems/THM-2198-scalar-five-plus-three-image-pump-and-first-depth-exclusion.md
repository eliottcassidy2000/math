---
id: THM-2198
title: "Scalar five-plus-three image pump and first-depth exclusion"
status: >
  PROVED + VERIFIED-EXACT. In the scalar five-unit/three-deep survivor, the
  mixed residual expands under multiplication by thirteen by a factor at
  least 13/10 on its first image. Until the least blocker valuation is
  exhausted, every image remains inside the correspondingly divided
  three-comb carrier. If the least valuation is unique, removing its unit
  owner leaves measure at least 2593/69300 which propagates to the second
  valuation. Independently, the deepest blocker has a positive-measure
  private transition word: repeated-minimum profiles force owner word
  (0,0,1), while a strict valuation ladder forces (0,0,1) and then
  (*,0,1). Sharp two-/three-comb caps classify equality as the ratios
  (1,13) and (1,12,13). An exhaustive exact 13^3 annulus argument eliminates
  the entire first profile (lambda_1,lambda_2,lambda_3)=(1,1,2), with a
  uniform conditional-capacity margin of 86. The scalar 5+3 branch remains
  open at all other profiles, so this is not a proof of LRC(14).
source: codex-2026-07-24-scalar-five-plus-three-transition
depends_on:
  - THM-1166
  - THM-2137
  - THM-2138
  - THM-2168
  - THM-2192
related:
  - THM-2184
  - THM-2196
  - THM-2197
script: 04-computation/lrc14_scalar_five_plus_three_image_pump_thm2198.py
output: 05-knowledge/results/lrc14_scalar_five_plus_three_image_pump_thm2198.out
script_sha256: 5129dff3175bb5c433b395c1f6be6e36b21fb86798510968a3e188df69c6c9b4
output_sha256: e96d69f6e06e05a61397c3a6d43ece2b8317a748b9ba52cdca13a811c95c2124
hash_basis: working-tree bytes (LF)
---

# THM-2198 -- the scalar `5+3` image pump

Use the scalar notation of THM-2168 and THM-2192:

```text
C_H={t:||Ht||>1/7},
D_a={t:||at||<1/14},
A_0=C_H minus union_(i=1)^5 D_(q_i),                 (1)
```

where `H,q_1,...,q_5` are positive thirteen-units, `H` is odd, and the
`q_i` are distinct. Let the three actual blockers, also distinct as
inherited from the LRC row, be

```text
b_j=13s_j,
lambda_j=nu_13(b_j),          lambda_1<=lambda_2<lambda_3, (2)
```

after relabelling. Thus the unique-deepest conclusion of THM-2192 has
already been used. Suppose, outside a null set,

```text
A_0 subset D_(b_1) union D_(b_2) union D_(b_3).      (3)
```

Put `T(t)=13t` on `R/Z`. The following sections give the image pump,
owner transition, and first exact discharge.

## 1. Root-image expansion

For a Borel set `A subset C_H`, define

```text
n_A(y)=#{x in A:T(x)=y}.                              (4)
```

The thirteen-root guard count in THM-2192 gives

```text
#{x:T(x)=y, x in C_H}
 =10 if y in E_H,
 =9  if y in C_H,                                    (5)
```

away from endpoints. Consequently

```text
0<=n_A(y)<=10,
{y:n_A(y)>0}=T(A)                                    (6)
```

up to finite endpoint conventions. Here `T(A)` is Borel: partitioning the
circle into the thirteen standard sheets writes it as a finite union of
Borel images under sheetwise Borel isomorphisms. Haar disintegration gives

```text
measure(A)=(1/13) integral n_A(y)dy
          <=(10/13)measure(T(A)).                     (7)
```

Therefore

```text
measure(T(A))>=(13/10)measure(A).                     (8)
```

For every Borel set `B`, multiplication by thirteen also obeys

```text
measure(T(B))>=measure(B).                            (9)
```

Indeed `B subset T^(-1)(T(B))`, while Haar invariance gives
`measure(T^(-1)(E))=measure(E)`. We also use below that the Lipschitz map
`T` sends null sets to null sets.

THM-2137 gives the sharp five-unit residual floor

```text
measure(A_0)>=961/6930.                              (10)
```

Combining (8)--(10), and writing `A_r=T^r(A_0)`, yields

```text
measure(A_r)>=12493/69300              for every r>=1. (11)
```

The number in (11) is not a union-bound heuristic. It is the residual mass
in (10), multiplied by the exact largest proportion `10/13` of a
thirteen-root fibre which can lie on the guard-danger side.

In the rooted sheet coordinates of THM-2197, the set

```text
Z(y)=A_0 intersection T^(-1)({y})                   (12)
```

is precisely the safe-sheet deficiency, and `n_(A_0)(y)=|Z(y)|`.
Thus `A_1` is its nonempty projection. The new ingredient is the occupancy
sidecar `|Z|`; merely recording whether `Z` is empty loses (7).

## 2. Division support and the first quantitative owner transition

If `b=13^r c`, then

```text
x in D_b     iff     T^r(x) in D_c.                  (13)
```

Multiplication by `T^r` sends null sets to null sets. Applying (13) to (3)
therefore gives, for every `0<=r<=lambda_1`,

```text
A_r subset union_(j=1)^3 D_(b_j/13^r)                (14)
```

almost everywhere.

Assume first that the least valuation is unique:

```text
lambda_1<lambda_2<lambda_3.                           (15)
```

Put

```text
a_1=b_1/13^lambda_1,
R_1=A_(lambda_1) minus D_(a_1).                       (16)
```

The coefficient `a_1` is a thirteen-unit, so its danger comb has measure
`1/7`. Equations (11), (14), and (16) give the exact positive remainder

```text
measure(R_1)
 >=12493/69300-1/7
 =2593/69300.                                        (17)
```

Moreover, for every `0<=u<=lambda_2-lambda_1`,

```text
T^u(R_1)
 subset D_(b_2/13^(lambda_1+u))
       union D_(b_3/13^(lambda_1+u)),                 (18)

measure(T^u(R_1))>=2593/69300.                        (19)
```

At the second valuation put

```text
a_2=b_2/13^lambda_2,
c_3=b_3/13^lambda_2,
R_2=T^(lambda_2-lambda_1)(R_1).                       (20)
```

Then `a_2` is a thirteen-unit, `13|c_3`, and

```text
2593/69300<=measure(R_2)
             <=measure(D_(a_2) union D_(c_3))
             <=25/91.                                (21)
```

The last inequality is the sharp pair-union cap obtained from the
THM-1166 intersection floor. Equality in the last inequality is possible
only when

```text
c_3=13a_2.                                           (22)
```

Thus saturation of the two-comb carrier forces

```text
lambda_3=lambda_2+1
```

and equality of the remaining unit parts. This is an equality-carrier
classification, not a claim that every survivor saturates (21).

## 3. A deepest-private transition word

The quantitative remainder (17) is not the only transition forced by the
cover. THM-2138 proves that five unit masks and any two positive-valuation
masks cannot cover `C_H`. Hence

```text
P_3=A_0 minus (D_(b_1) union D_(b_2))                (23)
```

contains a nonempty open set. The full cover (3) makes `b_3` its unique
deep owner:

```text
P_3 subset D_(b_3)                                   (24)
```

almost everywhere. Its images have positive measure by (9). Equation (13)
also preserves every owner bit whose coefficient can still be divided:

```text
T^r(P_3) subset D_(b_3/13^r)
                      for 0<=r<=lambda_3,             (25)

T^r(P_3) intersection D_(b_j/13^r)=empty
                      for j=1,2 and r<=lambda_j,     (26)
```

up to null boundaries.

This gives the exact minimal owner words.

### Repeated least valuation

If

```text
lambda_1=lambda_2<lambda_3,                           (27)
```

then at the first transition `r=lambda_1`,

```text
T^r(P_3)
 subset D_(b_3/13^r) minus (D_(a_1) union D_(a_2)),  (28)
```

with positive measure. In blocker order `(1,2,3)`, the forced word is

```text
(0,0,1).                                             (29)
```

At the same level (14) and THM-1166 give

```text
12493/69300<=measure(A_r)<=36/91.                    (30)
```

If the upper bound is attained, the three normalized blockers are a scaled
permutation of

```text
(1,12,13).                                           (31)
```

Because exactly two normalized blockers are thirteen-units, (31) means

```text
{a_1,a_2}={g,12g},       b_3/13^r=13g               (32)
```

for one thirteen-unit `g`. In particular the valuation gap is one.

### Strict valuation ladder

Under (15), equation (26) first gives

```text
(0,0,1) at r=lambda_1.                               (33)
```

At `r=lambda_2` the first blocker can no longer be divided integrally, so
its old zero bit is honestly lost. The surviving word is exactly

```text
(*,0,1).                                             (34)
```

Indeed `T^lambda_1(P_3) subset R_1`, and hence
`T^lambda_2(P_3) subset R_2 minus D_(a_2)`. The latter set has positive
measure, although (17) does not give it a uniform mass floor. This
distinguishes:

```text
quantitative transition:
    R_2 has mass at least 2593/69300 in the last two carriers;

topological owner transition:
    a positive piece of R_2 is privately owned by the deepest carrier. (35)
```

The star in (34) is not a cosmetic wildcard. Multiplication beyond
`lambda_1` has destroyed the integral division identity for `b_1`.

## 4. Exact exclusion of the first valuation profile

The first unique-deepest profile is

```text
(lambda_1,lambda_2,lambda_3)=(1,1,2).                (36)
```

It is empty.

Put

```text
N=13^3=2197,
U_N={z mod N:13 does not divide z and 7||z||_N>N}.   (37)
```

Before normalization, the primitive guard-danger annulus is

```text
U_N(H)={z mod N:13 does not divide z and 7||Hz||_N>N}.
```

Because `H` is a thirteen-unit, the coordinate bijection
`z' = Hz mod N` sends `U_N(H)` onto the canonical set `U_N` in (37).
In the new coordinate every terminal coefficient `a` becomes
`aH^(-1) mod N`; this preserves unit/sign classes and all three
thirteen-valuations. The five unit masks are indexed, up to sign, by the

```text
phi(13^3)/2=1014
```

unit classes modulo `13^3`.

Write the two depth-one blockers as `13a,13b`. Their masks on `U_N` depend
only on the sign classes of `a,b` modulo `13^2`; there are

```text
phi(13^2)/2=78
```

such classes. Repetition is allowed, because distinct integer blockers may
have the same sign class at this modulus. Thus there are exactly

```text
78*79/2=3081                                           (38)
```

unordered shallow pairs with repetition.

The depth-two blocker is invisible on `U_N`: for a unit numerator `z`,

```text
13^2 c z/13^3=cz/13
```

is a nonzero thirteenth root, whose circle norm is at least
`1/13>1/14`.

Fix one shallow pair and let `R_(a,b)` be the subset of `U_N` safe from
both shallow blockers. For each of the `1014` unit sign classes `q`, put

```text
m_q(a,b)=#{z in R_(a,b):14||qz||_N<N}.               (39)
```

Order these `1014` integers decreasingly:

```text
m_(1)>=...>=m_(1014).                                 (40)
```

Any five actual unit coefficients use at most five distinct sign-class
masks. Repeating a class does not enlarge their union. After duplicates are
removed, the elementary union bound therefore gives the
conditional-capacity upper bound

```text
#(R_(a,b) intersection union_(i=1)^5 unit_mask_i)
 <=m_(1)+...+m_(5).                                  (41)
```

This is why five times the largest individual mask would be the wrong
quantity: it spends the same residue class five times even though its union
is idempotent.

The exact audit enumerates the full universes in (38)--(40) and proves

```text
#R_(a,b)-sum_(i=1)^5 m_(i)>=86                        (42)
```

for every shallow pair, including repetitions. The unique worst row is

```text
(a,b)=(14,46),
#R=1046,
(m_(1),...,m_(5))=(204,200,190,186,180),             (43)
```

with unit labels `(183,799,599,1000,1007)`, so the conditional cap is
`960`.

Thus at least one primitive guard-danger torsion residue is safe from all
eight terminal masks. No endpoint ambiguity can hide it: `13^3` is coprime
to both seven and fourteen, so every guard and terminal inequality at such
a residue is strict. The uncovered residue thickens to an uncovered open
interval, contradicting the almost-everywhere cover. This proves (36).

The companion performs direct depth-one versus reduced-mask checks for all
`78` classes, verifies that all twelve nonzero depth-two unit parts give the
empty mask, checks all `3081*1014` conditional incidences, and freezes the
full table digest

```text
d31e5c874d8b5893ff33fc35095c18dcdf865c7f8296285c1b3441a8d8d679d9.
```

Normal and optimized Python outputs are byte-identical.

## 5. Carry-chart and continuation boundary

THM-2197 proves that static mask addition factors through the Boolean
deficiency semilattice. Equations (7), (18), and (25) explain why the
all-depth transition does not: deletion and multiplication require
occupancy, blocker labels, and the coefficient winding which says when
division by thirteen is still integral.

THM-2196 gives every zero-Haar row a finite semilinear relation/carry chart.
The present owner word factors through such a chart only after retaining
the **full carry vector and resulting integer coefficients**. A stabilized
chart or its finite residue alone does not determine the valuations in (2),
the sign classes in (39), or the event at which the star in (34) appears.
Consequently THM-2196's qualitative projective finiteness does not turn
(42) into an effective audit of all remaining profiles; coefficient winding
is still the required sidecar.

The exact connection ledger is

```text
source:       mixed scalar residual A_0;
map:          repeated multiplication by thirteen;
preserved:    residual mass, while support owners divide integrally;
lost:         a blocker bit immediately after its valuation is exhausted;
sidecar:      rooted occupancy, labelled owner incidence, valuation,
              and full coefficient/carry winding;
decisive test:
              the (1,1,2) annulus, where all 3081 shallow pairs fail. (44)
```

## 6. Scope

The theorem:

1. supplies the first quantitative multiplication-by-thirteen pump;
2. separates repeated-minimum and strict-ladder valuation profiles;
3. forces a positive deepest-private transition through the second
   valuation;
4. identifies the sharp `(1,13)` and `(1,12,13)` equality carriers; and
5. eliminates the complete first profile `(1,1,2)`.

It does not eliminate `(1,1,lambda_3)` for `lambda_3>=3`, any strict ladder,
or any profile with larger common depth. It gives no uniform lower measure
for the private piece after the second unit owner is removed, no upper bound
on the unit boundary height, and no proof of LRC(14). QED.
