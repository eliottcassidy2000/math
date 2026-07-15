---
id: THM-820
title: Scale-one Hamming-five finite collar reduction
status: PROVED (uniform two-box reduction) + FINITE-EXACT (792 height-one H5 rows, 924 height-one H6 rows, and 201,376 five-cycle band rows)
source: codex-2026-07-15-S10 radius-five continuation
depends_on:
  - LRCUpTo13  # only the settled seven-speed 1/8 and eight-speed 1/9 bounds
  - THM-806    # universal collar, half-open handoff rule, and no 2-/3-cycle lemmas
  - THM-815    # exceptional four-cycle classification and spread at most four
related: [THM-795, THM-800, THM-804, THM-810, THM-816, HYP-6820]
verification:
  - 04-computation/lrc13_scale_one_hamming_five_finite_reduction_codex_S10.py
  - 05-knowledge/results/lrc13_scale_one_hamming_five_finite_reduction_codex_S10.out
---

# THM-820 — Scale-one Hamming-five finite collar reduction

Put

```text
delta=1/13,                         [12]={1,...,12},
M(W)=max_(t in R/Z) min_(w in W)||wt||.
```

## Theorem

Let `R` be a five-element subset of `[12]`, put `P=[12] minus R`, and
choose a positive lift height `h_r>=1` for every `r in R`.  Write

```text
B=P union {u_r:r in R},             u_r=r+13h_r,          (1)
```

and sort the five replacements numerically as

```text
x<v<w<y<z.                                                   (2)
```

If `B` is tight at `delta`, then the following uniform finite reduction
holds.  With `m=max(P)<=12`,

```text
x<=floor(1456m/45)<=388.                                    (3)
```

Moreover exactly one of the following two search boxes contains the row.

### A. Recursive doubling box

```text
v<=2x,            w<=2v,            y<=2w,            z<=2y. (4)
```

In particular the five numerical caps are

```text
(x,v,w,y,z)<=(388,776,1552,3104,6208).                     (5)
```

### B. Exceptional top-four box

If `v>2x`, the top four labels form one multiplicative translate of

```text
{1,2,4,8} in F_13^*.                                       (6)
```

Their owner-handoff graph is THM-815's exceptional four-cycle, and

```text
x<286m/15,                         hence x<=228,
v<=floor(819x/40),
v<=floor((7/2)/(15/(104m)-1/x)) when the denominator is positive,
v<=1986,                    max{v,w,y,z}<=4v<=7944.         (7)
```

Thus the entire proper scale-one Hamming-five chart is reduced uniformly to
two explicit finite boxes.  This theorem is a **finite reduction**, not a
claim that every arbitrary-height row in those boxes has already been
rejected.

### C. Exact height-one slice

For `h_r=1` on all five labels, every one of the `C(12,5)=792` packets is
loose.  The exact minimum is

```text
min M(B)=1/12,                                             (8)
```

attained uniquely for

```text
R=(1,3,5,7,9),
B=2[11] union {11}.                                       (9)
```

Its complete maximizer set is

```text
{1,5,7,17,19,23}/24.                                     (10)
```

### D. The next height-one seam

Among all `C(12,6)=924` height-one six-replacement packets, the minimum is
`1/13` and is again unique:

```text
R=(1,3,5,7,9,11),                  B=2[12].               (11)
```

This is a doubled AP, not a new non-AP residual.  It is a method boundary:
radius five approaches a scale-changing AP face, and at radius six the face
is reached exactly.

## 1. Complete residues and the periodic danger estimate

Packet (1) contains exactly one representative of every nonzero residue
modulo `13`.  At every nonzero thirteenth it therefore has clearance exactly
`1/13`, so

```text
M(B)>=1/13.                                                (12)
```

The only case requiring exclusion is equality.

For a positive speed `u`, let

```text
D_u(delta)={t:||ut||<=delta}.
```

THM-806's periodic-danger estimate says that on any interval `I` of length
`L`,

```text
meas(I intersect D_u(delta))<=2L/13+2/(13u).              (13)
```

It is obtained by scaling `I` by `u`, counting complete periods, and bounding
the two end fragments by one danger tooth.  All reciprocal constants below
come directly from (13).

## 2. The seven-core interval gives the preliminary least bound

The retained core `P` has seven speeds.  The settled seven-speed Lonely
Runner bound gives a point where every `p in P` has clearance at least
`1/8`.  Since `m=max(P)`, its open `delta`-safe Lipschitz interval has length

```text
L_7=2(1/8-1/13)/m=5/(52m).                               (14)
```

If `B` were tight, the five replacement danger combs would cover this
interval.  Summing (13) five times gives

```text
L_7 <= 10L_7/13 + (2/13)sum_(r in R)1/u_r,
sum_(r in R)1/u_r >= (3/2)L_7 = 15/(104m).               (15)
```

Every replacement is at least `x`, hence the left side is at most `5/x`.
Consequently

```text
x <= (5)(104m)/15 = 104m/3 <= 416,                       (16)
```

This is the preliminary uniform cap.  It already makes the least coordinate
finite without any collar assumption.  Section 4 uses the exact cycle bank to
sharpen `416` to the theorem-facing `388` in (3).

## 3. The exact collar relation

At the inverse thirteenth belonging to an owner label `r`, THM-806 supplies
the universal left core-safe collar of length `1/156`.  This proof uses only
that retained speeds lie in `[12]`; it remains valid for the present
seven-speed core.

If `u_r>24`, the left endpoint of the owner's aligned danger tooth lies
strictly inside the collar.  Tightness then forces a different replacement
to cover the left germ.  Write

```text
q -> r       when provider q covers owner r,
lambda=u_q/u_r,
z=q r^(-1) in F_13^*.                                    (17)
```

The exact half-open condition is

```text
-1 < z-2lambda-13j <= 1               for some j in Z.  (18)
```

Putting `k=z-13j`, its ratio band is

```text
(k-1)/2 <= lambda < (k+1)/2.                            (19)
```

In particular, for a subunit ratio,

```text
lambda<1  =>  z=2 and 1/2<=lambda<1.                    (20)
```

Thus a provider less than half the owner can never discharge that owner.
THM-806 also proves that (18) admits no directed two-cycle and no directed
three-cycle.  Both proofs are local to the alleged cycle, so they remain
valid inside an induced subgraph on any subset of the five owners.

## 4. The collar dichotomy

Assume first that `v>2x`.  Since every proper lift is at least `14`,

```text
v>28>24.                                                  (21)
```

All four owners `v,w,y,z` therefore have collar obligations.  By (20), the
anchor `x` is too small to provide to any of them.  Hence the **induced graph
on the top four owners** has positive indegree at every vertex.  It contains
a directed cycle; the no-two- and no-three-cycle lemmas force a directed
four-cycle.  Only at this point, with all four induced owners known to exceed
`24`, do we apply THM-815.  Its four-cycle classification gives exactly (6),
the centre word `(2,2,2,5)`, and speed spread at most four.

Now suppose `v<=2x`.  If `w>2v`, then `w,y,z>24`, and neither `x` nor `v`
can provide to the induced top three.  Positive indegree there would force a
forbidden two- or three-cycle.  Thus `w<=2v`.  If `y>2w`, the induced top two
would have to provide to each other, a forbidden two-cycle.  Thus `y<=2w`.
Finally, if `z>2y`, every other replacement is below half of `z`, so `z` has
no provider.  This proves `z<=2y`, hence (4).

It remains to sharpen the least-coordinate cap used in (5).  If `x<=24`,
then (3) is automatic: a seven-element subset of `[12]` has `m>=7`, and
`floor(1456*7/45)>24`.  Suppose instead that `x>24`.  Every one of the five
owners now has positive collar indegree.  A minimal directed cycle has length
four or five.

- In length four, THM-815 gives the centre word `(2,2,2,5)`.  One speed is at
  least twice another and hence at least `2x`, so

  ```text
  sum_(r in R)1/u_r <= 4/x+1/(2x)=9/(2x).                (21a)
  ```

- In length five, the exact bank in Section 6 shows that every feasible word
  contains a centre at least four.  Its provider/owner ratio is at least
  `3/2`, so one speed is at least `3x/2` and

  ```text
  sum_(r in R)1/u_r <= 4/x+2/(3x)=14/(3x).               (21b)
  ```

The five-cycle bound is weaker and therefore uniform.  Combining (21b) with
the lower bound (15) gives

```text
15/(104m) <= 14/(3x),
x<=floor(1456m/45)<=floor(1456*12/45)=388.               (21c)
```

This proves (3), and recursive doubling gives the sharpened caps (5).

The distinction from radius four is structural.  One small anchor no longer
forces immediate doubling: four remaining large owners can support the
exceptional four-cycle.  The next section supplies the missing metric input.

## 5. Reciprocal mass bounds the exceptional cycle

The exceptional four-cycle already improves the anchor cap before we adjoin
anything.  With `v` the least top speed, its centre-five edge gives
`sum_top 1/u<=7/(2v)`.  Since this branch assumes the strict inequality
`v>2x`,

```text
sum_(r in R)1/u_r = 1/x+sum_top 1/u
                   < 1/x+7/(4x)=11/(4x).                (21d)
```

Combining (21d) with (15) preserves strictness:

```text
15/(104m) < 11/(4x),
x < 286m/15 <= 3432/15=228.8,
hence the integer x<=228.                                  (21e)
```

Now adjoin the anchor to the retained core:

```text
Q=P union {x}.                                             (22)
```

This is an eight-speed set.  Because `x>=14>max(P)`, its maximum speed is
exactly `x`.  The settled eight-speed Lonely Runner bound supplies clearance
at least `1/9`, hence a `delta`-safe interval of length

```text
L_8=2(1/9-1/13)/x=8/(117x).                              (23)
```

The four top danger combs must cover it.  Four applications of (13) give

```text
L_8 <= 8L_8/13 + (2/13)sum_top 1/u,
sum_top 1/u >= (5/2)L_8 = 20/(117x).                    (24)
```

THM-815's centre-five edge has one speed at least twice another.  With `v`
the least of the top four, one reciprocal is therefore at most `1/(2v)` and
the other three are at most `1/v`.  Thus

```text
sum_top 1/u <= 7/(2v).                                   (25)
```

Combining (24)--(25),

```text
v <= (7/2)(117x/20)=819x/40.                             (26)
```

There is a second top-only lower bound already inside the seven-core estimate.
Subtracting the anchor reciprocal from (15) gives

```text
sum_top 1/u >= 15/(104m)-1/x.                            (26a)
```

When the right side is positive, combine (26a) with (25):

```text
v <= (7/2)/(15/(104m)-1/x).                              (26b)
```

Thus the exact row-dependent cap is the minimum of (26) and (26b), omitting
(26b) when its denominator is nonpositive.  Here `7<=m<=12` and, by (21e),
`14<=x<=228`.  The replay takes the exact rational floor of both caps on this
`6*215`-row superset.  The maximum of their minimum is

```text
v<=1986,                                                   (26c)
```

uniquely at `(m,x)=(12,97)`, where the two floored caps are `(1986,2046)`.
At `x=98` the second cap has already fallen to `1928`.  With the cycle spread
at most four, `max_top<=4v<=7944`.

There is also a short monotonic check behind the exact loop.  The second cap
is weakest at `m=12`; where positive, its unrounded value is

```text
4368x/(15x-1248),
```

which decreases with `x`, while `819x/40` increases.  Their real crossover is
`x=4384/45`, strictly between `97` and `98`.  Hence rows through `97` are
bounded by the increasing cap, rows from `98` onward by the decreasing cap,
and the two displayed floor calculations prove the uniform maximum `1986`.
The replay's rational loop independently mirrors this argument.  This proves
(7).

This is the metric step that the bare collar graph cannot see.

## 6. Spanning five-cycle classification

The new possible cycle length is small enough to classify completely.  This
is the input used in the sharpened least-coordinate bound (21b).  On a
directed five-cycle,

```text
product lambda_i=1,              product z_i=1 mod 13.   (27)
```

Every `lambda_i>=1/2`; hence every `lambda_i<=16`.  Consecutive labels are
distinct, so an integer centre cannot be `0` or `1` modulo `13`.  The complete
centre bank is therefore

```text
{2,...,12,15,...,25,28,...,33}.                          (28)
```

For the five sorted centres, multiplying (19) gives the exact product filters
(necessary for the labelled cycle, and sufficient for choosing positive real
ratios of product one before enforcing distinct intermediate labels)

```text
product(k_i-1)<=32<product(k_i+1),
product k_i=1 mod 13.                                    (29)
```

The exact enumeration of all

```text
C(28+5-1,5)=C(32,5)=201,376                             (30)
```

nondecreasing multisets leaves six numerical words:

```text
(2,2,2,2,9),   (2,2,2,2,22),   (2,2,2,3,6),
(2,2,2,4,11),  (2,2,3,3,4),    (2,2,3,5,5).             (31)
```

The last word is impossible for distinct labels: its lower product is already
`32`, so every ratio must equal its lower endpoint; in particular its
centre-three edge has ratio one, forcing equal adjacent speeds.  The other
five types are genuinely realized by chordless all-large handoff cycles.
Up to cyclic rotation they give `16` band words and `15` multiplicative
label-set orbits.  The replay prints an exact realization of every type.

## 7. Exact height-one census and the doubled-AP seam

For a finite speed set, `min_u ||ut||` is piecewise linear.  A positive local
maximum occurs at a self-cusp or where two signed branches meet.  It is
therefore enough to scan exact rational candidates with denominator

```text
2u,                   u+v,                   or |u-v|.   (32)
```

The replay uses integer residues at every such candidate.  On the `792`
height-one Hamming-five rows, there are `66` distinct maximum values.  The
lowest strata are

```text
value       1/12    1/11    3/31    1/10
row count       1       2        5        3.             (33)
```

The unique minimum packet (9) has a transparent recursive description.
Its subset `2[11]` is a scaled eleven-term AP and has maximum `1/12`.  Its
eight AP equality clocks are

```text
{1,5,7,11,13,17,19,23}/24.                              (34)
```

The additional speed `11` destroys exactly `11/24` and `13/24`, leaving
(10).  Thus its maximum is exactly `1/12`; the full census proves uniqueness.

Lifting the next odd label replaces `11` by `24` and turns the packet into
`2[12]`.  The complete height-one radius-six census verifies that this is the
unique minimum row and that its maximum is exactly `1/13`, with maximizers

```text
{a/26:1<=a<=25, a!=13}.                                  (35)
```

Consequently the shallow charts are not a monotone sequence of unrelated
finite boxes.  They are approaching a scale-changing AP face: radius five is
the codimension-one predecessor `2[11] union {11}`, and radius six completes
the doubled AP.

## 8. A persistent top-four SCC refutes naive doubling

Let

```text
R={1,2,3,4,8},                         u_3=16,
c_N=1+13N,
(u_1,u_2,u_4,u_8)=c_N(79,54,30,34).                    (36)
```

Multiplication by `c_N` preserves every top label modulo `13` and every top
speed ratio.  For every `N>=0`, the induced top graph is exactly

```text
1 -> 8 -> 4 -> 2 -> 1                                  (37)
```

and the full live relation has SCC sizes `(4,1)`.  The edges from the top SCC
to the anchor vary with `N mod 8`: replacing `N` by `N+8` changes the scaled
phase by a multiple of `13`, and the replay checks the eight exact cases.
There is never an edge back from the anchor: for anchor-provider label `3`,
the subunit residue quotients toward owners `1,2,4` are `3,8,4`, rather than
the required `2`, while toward owner `8` the quotient is `2` but
`16/(34c_N)<1/2`.  Thus the top-four scale tends to infinity while the least
replacement stays `16` and the collar SCC remains unchanged.

This family does not consist of tight packets: at `N=0`, for example,
`M=5/21`.  It is a method-limit family.  It cleanly refutes the naive radius-
five extension of THM-815's recursive doubling argument.  The reciprocal
component step (23)--(26), which uses full tightness rather than only local
collar obligations, is essential.

## 9. Tournament Analysis and assumption challenge

At the collar stage the pair observable is the antisymmetric left-handoff
difference from (18).  A live edge is oriented provider-to-owner.  Silent
pairs can be completed by increasing label or by decreasing label; deleting
one edge from a live directed cycle gives a shared tie Hamiltonian path.

The spanning cycle examples with centre words

```text
(2,2,2,2,9)       and       (2,2,2,2,22)                (38)
```

have the same ordered labels `(1,7,10,5,9)`, exactly the same five live edges,
and the same completed tournament fingerprints.  In the increasing-label
gauge the score histogram is `{1:1,2:3,3:1}`, there are four directed
triangles, one SCC of size five, and thirteen Hamiltonian paths; the reverse
gauge has the same score and triangles and eleven Hamiltonian paths.  Yet the
closing ratio lies respectively near `9/2` and `22/2`, and the full packet
maxima are respectively `21/109` and `107/559`.

The persistent family (36) is even sharper.  Its live SCC sizes stay `(4,1)`;
both tie gauges have score histogram `{1:2,2:1,3:2}`, three directed triangles,
one completed SCC, and nine Hamiltonian paths, while the metric scale is
unbounded.

The theorem-bearing objects are therefore not bare runners or bare arcs.
There are two typed carriers at different depths.

1. **Owner-exit / provider-tooth incidence.**  A local object retains the
   owner endpoint, provider tooth, residue quotient `z`, integer centre `k`,
   speed ratio, and half-open boundary flag.  Its tournament shadow is useful
   for finding cycles but forgets reciprocal scale.
2. **Dynamic component / remaining-tooth incidence.**  Put
   `U={u_r:r in R}`.  For a set `S` of already adjoined replacements, define

   ```text
   E_S={t:min_(q in P union S)||qt||>1/13}.               (39)
   ```

   Retain every connected component of `E_S`, its exact endpoints and width,
   and which teeth of each remaining speed meet or contain it.  Passing from
   `S` to `S union {u}` is exact interval subtraction.  The final LRC
   predicate is precisely whether `E_U` is empty, so this alternating
   component-obligation carrier is predicate-preserving.

The assumption challenge is explicit.

- Runner vertices preserve pair labels but destroy simultaneous component
  coverage.
- Residue vertices preserve the packet orbit but destroy centre `k` and
  speed ratio, as (38) demonstrates.
- Gaps or fixed circle sections lose moving tooth scale unless exact endpoints
  are attached.
- Bare wall events preserve order but lose owner identity and component width.
- Fourier modes preserve average danger mass but not the connected component
  needed for exact containment.
- Bare proof obligations forget whether their witness intervals are mutually
  compatible.

The correct recursive object is a typed incidence tree whose tournament is a
controlled projection, not a tournament on runners pretending to be the
predicate itself.

## Reproduction and certificate digests

```bash
python3 04-computation/lrc13_scale_one_hamming_five_finite_reduction_codex_S10.py
```

The height-one census hashes, for every row in lexicographic missing-label
order, the labels, reduced maximum, number of maximizers, and every reduced
maximizer as big-endian unsigned 64-bit fields.  The digests are

```text
Hamming five  323913c07710f8321d29d07c9dfdfb4d13107549063587584b4ee0ceab923dc9
Hamming six   c41767013f5403a00004609ba31b165896a9a94416e0693393aecd9915737ccc
```

The replay byte-matches the stored output, whose SHA-256 is

```text
fb5f88213b591b69a6f6565b91776f9ab00ffa03fa6df35ebebb5a4264268eab.
```

The result was initially headed for THM-819, but live `main` assigned THM-819
to the primitive harmonic law while this computation was running.  It was
renumbered to THM-820 before publication; no THM-819 radius-five artifact was
created.
