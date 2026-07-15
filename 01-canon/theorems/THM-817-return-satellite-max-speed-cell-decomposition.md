---
id: THM-817
title: Return satellites are max-speed cells and admit an adaptive erosion selector
status: PROVED (exact max-speed-cell decomposition and adaptive selector) + VERIFIED (independent Fraction interval replay and sharp linear-satellite family)
source: codex-2026-07-15-S10 return-satellite continuation
depends_on:
  - THM-803
  - THM-807
related:
  - THM-772
  - THM-789
  - THM-797
  - HYP-6820
verification:
  - 04-computation/lrc13_return_satellites_cell_classifier_codex_S10.py
  - 05-knowledge/results/lrc13_return_satellites_cell_classifier_codex_S10.out
---

# THM-817 — Return-satellite max-speed-cell decomposition

Use the two-sheet notation

```text
gamma=2/143,
R_U={d:max_(u in U)||ud||<gamma},
B=max(U),
L=sum_(u in U)u.
```

Every component of `closure(R_U)` is one explicitly computable interval
inside a tooth of the maximum speed.  Consequently its component count
`N_R` is at most `B`, and the THM-803 erosion selector admits the adaptive
bound

```text
|Sigma(U;x,y)| <= 2c_E N_R+2W-2g
                <= 20B^2+22B-2g,
```

where `c_E` is the number of deep-set components, `W=max(x,y)`, and
`g=gcd(x,y)`.  This replaces the former universal `200B^2+22B` estimate
and specializes to THM-807 when `N_R=1`.

The bound `N_R<=B` is of the correct order under the present arithmetic
hypotheses: an explicit primitive, divisor-complete, exact signed-complement
family has `N_R=3+1440n`.  Thus signed-wall support, divisor completeness,
the parity bridge, and the named scalar taxes do not imply connected,
bounded-component, or sublinear return geometry.

## 1. Exact max-speed cell theorem

For `k in Z/BZ`, let `[uk]_B` be the balanced residue in
`(-B/2,B/2]`.  Define the open scaled cell

```text
I_k=(-gamma,gamma)
    intersect intersection_(u in U)
      ((-B gamma-[uk]_B)/u, (B gamma-[uk]_B)/u).          (1)
```

Then

```text
R_U = disjoint_union_(k:I_k nonempty) (k+I_k)/B,         (2)

closure(R_U)
    = disjoint_union_(k:I_k nonempty) (k+closure(I_k))/B. (3)
```

The unions are circular; using balanced `k` puts every represented interval
inside `[-1/2,1/2]`.  A zero-width intersection of the closed inequalities is
correctly omitted: `R_U` is open, so an isolated closed-feasible point is not
in its closure.

### Proof

The return set of the largest speed `B` is the disjoint union of the `B`
teeth

```text
d=(k+e)/B,                    |e|<gamma.                 (4)
```

Fix one tooth.  The `u`-constraint is

```text
|uk+ue-mB|<B gamma                                  (5)
```

for some integer `m`.  If (5) is feasible, then

```text
|uk-mB| < B gamma+u gamma <=2B gamma=4B/143<B/2.     (6)
```

Thus `m` is unique and `uk-mB=[uk]_B`.  Solving (5) for `e` gives exactly the
`u`-interval in (1).  An intersection of intervals is an interval, so no
largest-speed tooth can split into several return components.  Distinct
`B`-teeth are separated by a positive gap.  This proves (2), and taking the
closure of the finite disjoint union proves (3).

This also classifies all endpoint owners.  The left owners are the speeds
attaining the maximum lower endpoint in (1); the right owners attain the
minimum upper endpoint.  No endpoint reconstruction from a pairwise sum of
anonymous arcs is necessary.

## 2. Component and label bounds

Let

```text
N_R(U)=#{k:I_k is nonempty}.                              (7)
```

The cell theorem immediately gives

```text
N_R(U)<=B.                                                (8)
```

There is also a cheap gcd label sieve.  A surviving label satisfies, for
every `u`,

```text
|[uk]_B| <= B gamma+u gamma <=4B/143.                    (9)
```

Put `g_u=gcd(u,B)`.  Multiplication by `u` modulo `B` takes each image residue
with multiplicity `g_u`, and its image residues are the multiples of `g_u`.
Consequently

```text
N_R(U)
 <= min_(u in U) g_u(2 floor(4B/(143g_u))+1).            (10)
```

Formula (10) deliberately includes boundary labels and is therefore a safe,
slightly conservative bound.  The exact test is still (1).

The zero cell is always

```text
closure(I_0)/B=[-gamma/B,gamma/B],                       (11)
```

recovering THM-807's central component.  Cells occur in `+/-k` pairs.  If
`U` is primitive, some speed is odd, so the self-antipodal `k=B/2` tooth is
impossible: its odd-speed phase remains at least `1/2-gamma>gamma`.  Hence

```text
N_R(U)=1+2*(number of satellite pairs)                   (12)
```

for primitive cores.

## 3. Adaptive compression of THM-803's selector

Let `c_E` be the number of closed circular components of

```text
E_U={t:phi_U(t)>=1/11}.
```

By the cell theorem, `closure(R_U)` has `N_R` components.  Each sum of one
`E_U` component and one return cell is one circular arc or the whole circle.
Thus

```text
K_U=E_U+closure(R_U)
```

is represented by at most `c_E*N_R` arcs, not `c_E^2` arcs.  With
`a=(x+y)/2`, `b=|x-y|/2`, `W=max(x,y)`, and `g=gcd(x,y)`, the exact
endpoint/cusp proof of THM-803 therefore gives

```text
|Sigma(U;x,y)| <= 2c_E N_R+2W-2g.                       (13)
```

Since `c_E<=L` and `N_R<=B`,

```text
|Sigma| <=2LB+2W-2g.                                    (14)
```

Under THM-772, `L<=10B` and `W<=11B`, so the former

```text
200B^2+22B
```

bound improves uniformly to

```text
20B^2+22B-2g.                                            (15)
```

More importantly, (13) is adaptive.  It recovers THM-807's linear selector
when `N_R=1`, and pays only for the satellite cells that actually survive the
ten simultaneous residue intervals (1).  Constructing those cells costs
`O(B|U|)` exact rational comparisons.

## 4. Exact signed-complement family with linearly many satellites

For every `n>=0`, put

```text
B_n=506+360360n,
U_n=(1,2,3,4,7,B_n-6,B_n-3,B_n-2,B_n-1,B_n).            (16)
```

The increment `360360=lcm(1,...,13)` preserves every relevant residue and
divisibility condition.  Directly:

```text
gcd(U_n)=1,
every m=2,...,12 divides a member of U_n,
no member is divisible by 13,
{u mod 13:u in U_n}=F_13^* minus {5,8},
P(U_n)={1,2,3,4,5,6}.                                   (17)
```

Thus these are primitive, divisor-complete, exact signed-complement cores
with the mandatory parity bridge.

### Exact return components

The speed `1` first confines return times to `|d|<gamma`.  On the positive
side, speed `7` sharpens this to

```text
0<=d<gamma/7.                                            (18)
```

Write `B_nd=k+e`, with `|e|<gamma`.  For the five high speeds
`B_n-r`, `r in {0,1,2,3,6}`, no phase can wrap again inside (18), and their
constraints reduce to

```text
6d-gamma < e < gamma.                                    (19)
```

Hence the positive `k`-cell is exactly

```text
((k-gamma)/(B_n-6),
 min((k+gamma)/B_n,gamma/7)).                            (20)
```

It is nonempty exactly when

```text
k < gamma(B_n+1)/7 = 2(B_n+1)/1001.                     (21)
```

The other possible endpoint comparison is weaker for every `B_n>=506`.
Now

```text
2(B_n+1)/1001 = 720n+1014/1001,                         (22)
```

so there are exactly `720n+1` positive satellites and the same number of
negative satellites.  Together with (11),

```text
N_R(U_n)=3+1440n.                                        (23)
```

This is a linear lower bound in `B_n`.  Therefore the exact signed complement,
the parity-twisted support gate, divisor completeness, and primitivity cannot
imply connected returns, boundedly many satellites, or any universal
`o(B)` return-component bound.

The family also passes the named scalar taxes for `(x,y)=(13,5)`.  The metric
pins and determinant inequality are immediate for `B_n>=506`.  Since adding
speeds can only decrease the maximin and

```text
M({1,2,3,4,7})=1/5,
```

the required metric bound is exact: the consecutive-speed pigeonhole lemma
gives the upper bound from `{1,2,3,4}`, while `t=1/5` gives distance at least
`1/5` for all five displayed speeds.  Therefore its
`rho=(M(U_n)-1/13)/B_n` satisfies

```text
1/(13*5)+2rho
 <=1/65+16/(65B_n)
 <=36/845=2/(13*13)+2/(13*5).                            (24)
```

This is a method limit for any attempt to infer connectedness from the current
arithmetic and scalar gates.  It is not a tight twelve-speed counterexample.
Indeed, using only the central return component of `U_0`, put

```text
s=479/616 in E_(U_0),         d=1/36179 in closure(R_(U_0)),
t=s+d=1575487/2026024.
```

The active equalities are

```text
phi_(U_0)(s)=1/11,            max_u ||ud||=2/143,
Q_(9,4)(t)=226661/2026024 < 11/13,
11/13-Q_(9,4)(t)=1487667/2026024.                         (25)
```

The last margin is strict, so moving `d` slightly into the open central return
cell still violates the erosion containment.  Thus the family isolates a
return-geometry method limit while remaining decisively outside the tight
branch.

### Smallest row

For `n=0`,

```text
U_0=(1,2,3,4,7,500,503,504,505,506)
```

and the exact closed return components are

```text
[-2/1001,-141/71500]       endpoint owners (7,500),
[-1/36179,1/36179]         endpoint owners (506,506),
[141/71500,2/1001]         endpoint owners (500,7).       (26)
```

Here `L=2535`, `N_R=3`, and `(x,y)=(13,5)`.  The coarse THM-803 arc bound
drops from

```text
2L^2+26 = 12,852,476
```

to the adaptive

```text
2L N_R+26-2 = 15,234.                                   (27)
```

## 5. Tournament telemetry and challenged vertices

For (26), take the three return components—not the ten runners—as vertices.
The pair observable is centre order on the signed circle cut at `-1/2`; switch
to component width, with signed cell label as the tie Hamiltonian path.  The
orders are

```text
(-1,0,1),                      (-1,1,0),                  (28)
```

with one edge flip.  Both tournaments are transitive: score histogram
`(0,1,2)`, no directed cycle, singleton SCCs, and one Hamiltonian path.

This is telemetry only.  The negative satellite hands off from endpoint owner
`7` to owner `500`, while the positive satellite hands off from `500` to `7`.
The reciprocal owner incidence vanishes in either bare tournament.  Runner
vertices lose cell labels and widths; denominator vertices lose endpoint
owners; unsigned component vertices identify the two reciprocal satellites;
fixed circle sections lose the maximum-speed tooth action; Fourier modes see
return density but not whether an intersection cell has positive width.

The predicate-preserving carrier is the signed max-speed cell

```text
(k, I_k, left owners, right owners),                     (29)
```

followed by its incidence with every closed component of `E_U`.  This is the
smallest quotient encountered here that preserves both disconnected-return
geometry and the erosion containment predicate.

## 6. Exact replay

The replay cross-checks the cell classifier against an independent literal
open-interval intersection on 213 deterministic primitive rows, reconstructs
the connected THM-803/807 benchmark rows, verifies (17)--(29), and records the
component-tournament telemetry.  Its digest is

```text
286777a892bb16a926fc2dffcb08c8832987da6682cfd3cfed0331cf9a5795f6.
```

No floating point or sampled-circle verdict enters the classifier.
