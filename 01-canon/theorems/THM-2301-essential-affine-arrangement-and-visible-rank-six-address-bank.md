---
id: THM-2301
title: "Essential affine arrangements and the visible rank-six address bank"
status: >
  PROVED + VERIFIED-EXACT. The sharp complement bound for m essential
  affine hyperplanes in F_q^d is (q-m+d-1)(q-1)^(d-1) in the stated
  range. The actual degree-526 Fourier survivor relations for THM-2298's
  comparison already span rational rank at least six. Combining THM-2275, THM-2295,
  and THM-2298 with centered p-adic saturation gives, for p=7 or 13,
  six mod-p-independent rows of heights at most
  20,196,196,206,309,526. In the no-dark-column branch this yields
  at depth n at least 9*12^5*13^(6n-6) unanchored and
  9*12^4*13^(5n-5) normalized thirteen-adic all-unit addresses, or
  3*6^5*7^(6n-6) and 3*6^4*7^(5n-5) septimal addresses. Intrinsic
  bright/dark alternatives and explicit prescribed-c_1 packet heights
  are proved. Septimal darkness is literal deletion from the complete
  degree-526 survivor support; septimal bright addresses are genuine
  nonzero interval Fourier terms, but their signed aggregate can cancel.
  No ancestry incidence, scalar profile exclusion, or proof of LRC(14)
  follows.
source: codex-2026-07-25-rank-six-kakeya-lift
depends_on:
  - THM-2284-thirteen-adic-anchored-rank-three-plucker-lift
  - THM-2275-mixed-scalar-relation-and-guard-blocker-crossing
  - THM-2295-weighted-basis-safe-floor-and-scalar-rank-five-harvest
  - THM-2298-weighted-rank-five-facet-deficit-and-scalar-rank-six-harvest
related:
  - THM-2282-thirteen-adic-saturation-and-unit-anchored-minor
  - THM-2294-anchored-plucker-tournament-and-kakeya-address-bank
script: 04-computation/lrc14_visible_rank_six_kakeya_lift_thm2301.py
output: 05-knowledge/results/lrc14_visible_rank_six_kakeya_lift_thm2301.out
script_sha256: f79a10e3184b35d0a78209d943462fadfb468761e6023fd93c2ad9b2bf3879f9
output_sha256: 3e8dee251a65e8219539421d1557504613df5c921b1c149af29ec3c079092305
hash_basis: working-tree bytes (LF)
---

# THM-2301 -- essential arrangements and visible rank-six addresses

**PROVED + VERIFIED-EXACT.**

THM-2294's anchored rank-three packet leaves at least `72` all-unit
addresses over `F_13`. The count is one instance of a sharp theorem for
essential affine arrangements. THM-2275 supplies a height-`20` seed,
THM-2295 supplies rational rank five at height `196`, and THM-2298 crosses
the rank-five facet to rank six at height `526`. These mechanisms compose,
but only after repairing two losses:

```text
rational rank need not survive reduction modulo p;
an all-unit address need not have a sign-definite Fourier contribution.
```

The first loss is repaired by THM-2284's centered saturation operator. The
second remains the exact stopping boundary. At seven, however, the address
bank does at least land inside the nonzero Fourier support of every interval
factor.

## 1. Sharp essential affine arrangements

Let `q` be a prime power and let

```text
H_1,...,H_m
```

be affine hyperplanes in `F_q^d`. Assume that their normals span the dual
space and that

```text
d<=m<=d+q-1.                                         (A1)
```

Then

```text
|F_q^d minus union_(i=1)^m H_i|
 >=(q-m+d-1)(q-1)^(d-1).                             (A2)
```

### Proof

Choose `d` hyperplanes with independent normals. Their common intersection
is one point. Translate it to the origin and make an invertible linear
coordinate change. The chosen hyperplanes become

```text
x_1=0,...,x_d=0.
```

Their complement is the torus `(F_q^*)^d`, of size `(q-1)^d`. Every
remaining affine hyperplane meets this torus in at most `(q-1)^(d-1)`
points: fix all but a variable with nonzero coefficient, and there is at
most one possible value of the last variable. A union bound inside the
torus gives

```text
(q-1)^d-(m-d)(q-1)^(d-1)
 =(q-m+d-1)(q-1)^(d-1).
```

Coincident hyperplanes only make the union smaller, so no deduplication
hypothesis is hidden here. QED.

The bound is sharp. Take the `d` coordinate hyperplanes and `m-d` distinct
nonzero translates

```text
x_1=c_1,...,x_1=c_(m-d).                             (A3)
```

The first coordinate then has `q-m+d-1` available nonzero values and every
other coordinate has `q-1`.

For at most eight affine hyperplanes in dimension `r-1` over `F_13`, (A2)
specializes to

```text
(r+3)12^(r-2),                    2<=r<=9.            (A4)
```

Thus the rank-three count is `6*12=72`, the rank-five count is
`8*12^3=13824`, and the rank-six count is

```text
9*12^4=186624.                                       (A5)
```

At `q=7,r=6`, the corresponding anchored count is

```text
3*6^4=3888.                                          (A6)
```

## 2. Column packets: unanchored and normalized counts

Let `R` be an `r by 9` matrix over `F_q` of row rank `r`, with columns

```text
V_0,...,V_8 in F_q^r.
```

Call coordinate `i` **dark** if `V_i=0`.

If no coordinate is dark, the nine equations

```text
lambda.V_i=0,                    0<=i<=8,             (B1)
```

form an essential central arrangement in the `r`-dimensional coefficient
space: the normals are the columns, and the columns span because `R` has
row rank `r`. There may be fewer than nine distinct hyperplanes, which only
improves (A2). Therefore the number of unanchored coefficient vectors
giving all-unit words is at least

```text
(q+r-10)(q-1)^(r-1).                                 (B2)
```

Now prescribe an anchor `a`. In the affine slice

```text
A_a={lambda:lambda.V_a=1},                           (B3)
```

the other eight zero equations are empty when `V_i` is a nonzero multiple
of `V_a`, and otherwise are affine hyperplanes. Their normals are the
images of `V_i` in

```text
F_q^r/F_q V_a.
```

Those images span the quotient. Applying (A2) in dimension `r-1` gives at
least

```text
(q+r-10)(q-1)^(r-2)                                  (B4)
```

normalized all-unit words. The unanchored all-unit set is exactly the
disjoint union of its `q-1` anchor-value slices; scaling each slice to
(B3) explains the factor `q-1` between (B2) and (B4).

Both constants are sharp for abstract rank-six packets. Take

```text
V_a=e_1,
V_1=e_2,...,V_5=e_6,
V_(5+j)=-j e_1+e_2,                j=1,2,3.          (B5)
```

On the normalized slice `(1,x_1,...,x_5)`, the eight forbidden
hyperplanes are

```text
x_1=0,1,2,3,       x_2=0,...,x_5=0.                 (B6)
```

They leave `(q-4)(q-1)^4` points. Without normalization, the same packet
leaves `(q-4)(q-1)^5`. These are (B4) and (B2) for `r=6`.
This is an abstract sharpness model; no positive-kernel or cover claim is
made for it.

## 3. The degree-526 visible survivor span already has rank six

For a live scalar cover, write

```text
w=(w_0,...,w_8) in Z^9,              w_i!=0,
Lambda(w)={m in Z^9:m.w=0}.                          (C1)
```

Let

```text
Xi_526={k in Z:
        |k|<=526 and (k=0 or 7 does not divide k)}
```

and define

```text
W_vis=span_Q(Lambda(w) intersection Xi_526^9),
L_vis=W_vis intersection Z^9.                        (C2)
```

Then

```text
dim_Q W_vis>=6.                                      (C3)
```

### Proof

Use the `N=264` squared-Fejer product `Q_264` from THM-2298. Its
one-variable factors approximate the concentric translated safe intervals
of lengths `5/7` and `6/7`. At every frequency `|k|<=526`, the Jackson
multiplier is positive. For `k!=0`, the remaining interval coefficient
vanishes exactly when `7|k`. Hence the exact Fourier support of `Q_264` is

```text
Xi_526^9.                                            (C4)
```

Suppose `s=dim_Q W_vis<=5`, and put `L=L_vis`. The lattice `L` is saturated,
has rank `s`, and lies in `Lambda(w)`. Every coordinate character is
nonzero on its annihilator torus `K_L`: if `e_i` belonged to `L`, then
`e_i.w=w_i=0`, contrary to (C1).

The finite line and torus survivor sets agree exactly:

```text
Xi_526^9 intersection Lambda(w) subset L
```

by the definition of `W_vis`, while

```text
Xi_526^9 intersection L
 subset Xi_526^9 intersection Lambda(w)
```

because `L subset Lambda(w)`. Thus the complete averages of `Q_264` on the
scalar line and on `K_L` agree.

The null scalar safe event bounds the line average by `9 eta_264`.
THM-2298's weighted rank-five facet theorem, with THM-2295 for lower ranks,
bounds the `K_L` average below by

```text
900/31213-9 eta_264.
```

Equality of averages would force

```text
900/31213<=18 eta_264,
```

contradicting THM-2298's exact positive `N=264` margin. Therefore (C3)
holds. QED.

This visible refinement makes septimal darkness exact. Every generator in
`Lambda(w) intersection Xi_526^9` has, in each coordinate, either a literal
zero or a seven-unit. Consequently, for every coordinate `i`,

```text
x_i=0 mod 7 for every x in L_vis

iff m_i=0 for every m in Lambda(w) intersection Xi_526^9

iff W_vis subset {x:x_i=0}.                          (C5)
```

The nontrivial forward implication holds because a visible generator whose
coordinate is divisible by seven must have that coordinate equal to zero.
Thus the dark alternative at seven is a literal support-deletion branch
for the complete degree-`526` survivor set, not just a congruence label.

## 4. Centered saturation and its exact height invoice

We use a general consequence of THM-2284's pivot-extension lemma. Let
`Lambda subset Z^N` be saturated and suppose

```text
dim_Q span_Q(Lambda intersection [-H,H]^N)>=s.       (D1)
```

For any odd prime `p`, there are rows `rho_1,...,rho_s` in `Lambda`
which are independent over `F_p` and satisfy

```text
B_1=H,
B_(k+1)=max(H,ceil((B_1+...+B_k)/2)),
||rho_k||_infinity<=B_k.                             (D2)
```

For the first row, repeatedly divide a nonzero bounded row by `p` while
its reduction vanishes. Saturation keeps the quotients in `Lambda`. Given
`k` rows, choose a bounded row outside their rational span and apply the
centered pivot-extension lemma. The repeated cancel-and-divide operation
preserves rational novelty, eventually extends the mod-`p` rank, and has
the height bound in (D2). This proves the assertion by induction.

There are two useful applications.

### The smallest universal packet

THM-2275 supplies a nonzero relation of height at most `20`. Divide by its
content; the primitive row has nonzero reduction at every prime. Starting
with this row, use four rationally new height-`196` rows from THM-2295 and
then one rationally new height-`526` row from THM-2298. At each step apply
the centered pivot extension. For either `p=7` or `p=13`, the successive
bounds are

```text
(20,196,196,206,309,526),
sum=1453.                                            (D3)
```

Indeed the first five partial sums are `20,216,412,618,927`; the next
bounds are respectively

```text
max(196,10), max(196,108), max(196,206),
max(196,309), max(526,ceil(927/2)).
```

This is the packet used for the sharp address heights in Section 5.

### The visible-carrier packet

For intrinsic statements about `L_vis`, seed the construction with
THM-2275's height-`462` crossing relation. Every nonzero coefficient of
that row is prime to seven. Dividing by its content preserves this property,
so the primitive row belongs to the visible generating set in (C2) and has
nonzero reduction modulo both seven and thirteen. Use rationally new
height-`526` visible generators thereafter. The bounds become

```text
(462,526,526,757,1136,1704),
sum=5111.                                            (D4)
```

Thus the visible carrier itself has a bounded mod-`p` rank-six packet.
Neither construction retains six rows at its original rational height;
the displayed growth is a uniform upper bound, not a claim of necessity.

## 5. All-depth rank-six address banks

Let `Vtilde_i in Z^6` be the columns of the packet in Section 4 and let
`V_i` denote their reductions modulo `p`.

The mod-`p` row rank supplies a `6 by 6` minor which is a `p`-unit.
Consequently the coefficient map

```text
(Z/p^n Z)^6 -> (Z/p^n Z)^9,
lambda |-> lambda R
```

is injective at every depth. The counts below are therefore counts of
distinct relation addresses, not merely coefficient presentations.

If some column is packet-dark, every packet word has coefficient divisible
by `p` there, so the packet has no all-unit word. If no column is dark,
Section 2 gives the following depth-one counts:

```text
p=13:
  unanchored >=9*12^5=2239488,
  anchor-normalized >=9*12^4=186624;

p=7:
  unanchored >=3*6^5=23328,
  anchor-normalized >=3*6^4=3888.                   (E1)
```

For `n>=1`, reduction of the unanchored coefficient space modulo `p` has
fibres of size `p^(6n-6)`. Reduction of a normalized slice

```text
A_(a,n)={lambda in (Z/p^n Z)^6:
         lambda.Vtilde_a=1 mod p^n}                  (E2)
```

has fibres of size `p^(5n-5)`. A coefficient nonzero modulo `p` remains a
unit at every lift. Therefore the all-depth counts are

```text
p=13:
  unanchored >=9*12^5*13^(6n-6),
  normalized >=9*12^4*13^(5n-5);

p=7:
  unanchored >=3*6^5*7^(6n-6),
  normalized >=3*6^4*7^(5n-5).                     (E3)
```

Choose centered representatives of the six combination coefficients:

```text
|lambda_j|<=(p^n-1)/2.
```

By (D3), every address in (E3) has an exact integer relation
representative of height at most

```text
1453(p^n-1)/2.                                       (E4)
```

The normalized count may be anchored at any packet-bright coordinate.

There is a useful labelled restriction on every thirteen-adic unit pivot.
Modulo thirteen, the three blocker scalar values vanish, while the six
guard/unit scalar values are nonzero. For every six-row relation packet,
the `6 by 6` submatrix on those six guard/unit columns annihilates their
nonzero scalar-weight vector. It is therefore singular. Consequently every
unit `6 by 6` pivot contains at least one blocker column.

The exact possible pivot atlas has

```text
1 blocker + 5 guard/units:  C(3,1)C(6,5)=18,
2 blockers + 4 guard/units: C(3,2)C(6,4)=45,
3 blockers + 3 guard/units: C(3,3)C(6,3)=20.        (E5)
```

There are `83` possible labelled pivot sets in total. Their three-column
complements have types `2+1`, `1+2`, and `0+3`, respectively. At least one
of the `83` pivots is a unit because the packet has row rank six. The atlas
does not prescribe which blocker occurs, and in particular does not force
the anchor `c_1`.

## 6. Packet-independent bright/dark alternatives

The preceding dark branch depends on the selected six rows. A larger
height makes it intrinsic to `L_vis`.

Continue (D2) through the full rank

```text
6<=rank(L_vis)<=8.
```

A mod-`p` basis of `L_vis`, continuing the visible recurrence (D4), has
height bounds dominated by

```text
462,526,526,757,1136,1704,2556,3834,
sum<=11501.                                          (F1)
```

### Thirteen

If no coordinate functional vanishes on `L_vis/13L_vis`, the union of its
nine coordinate kernels has at most

```text
9*13^(rank(L_vis)-1)<13^rank(L_vis)                  (F2)
```

points. Choose a residue outside their union and lift it using centered
basis coefficients in `[-6,6]`. This produces one all-unit relation row
of height at most

```text
6*11501=69006.                                       (F3)
```

Replace one basis row by this row and select five of the other basis rows.
The five largest bounds in (F1) sum to

```text
3834+2556+1704+1136+757=9987.
```

The resulting bright rank-six packet has row-height sum at most

```text
69006+9987=78993.                                    (F4)
```

Because it contains an all-unit row, every coordinate is packet-bright.
It may be normalized at any prescribed coordinate, including `c_1`, and
its address representatives have height at most

```text
78993(13^n-1)/2.                                     (F5)
```

Otherwise a coordinate is intrinsically dark modulo thirteen on `L_vis`.
This is a genuine carrier congruence but need not be literal support
deletion, because a visible coefficient can be a nonzero multiple of
thirteen.

### Seven

Suppose no coordinate is intrinsically dark on `L_vis/7L_vis`. Choose a
six-dimensional subspace contained in none of the nine coordinate kernels.
If `r=rank(L_vis)`, the proportion of six-subspaces contained in one fixed
hyperplane is

```text
[r-1 choose 6]_7/[r choose 6]_7
 =(7^(r-6)-1)/(7^r-1)
 <7^(-6).                                            (F6)
```

The union bound over nine kernels is strictly less than one, so the desired
subspace exists. Express one of its bases in the basis (F1), with centered
coefficients in `[-3,3]`. Each of the six resulting rows has height at most
`3*11501=34503`, and their total height is at most

```text
6*34503=207018.                                      (F7)
```

This bright packet may again be normalized at `c_1`; its representatives
have height at most

```text
207018(7^n-1)/2.                                     (F8)
```

If the bright case fails, (C5) says that one coordinate is literally zero
on every degree-`526` visible survivor and on all of `W_vis`. This is the
sharp support-dark alternative.

## 7. What the septimal bank does and does not prove

Let `J_0` be the translated guard-safe interval of length `5/7`, and let
`J_1,...,J_8` be the translated ordinary safe intervals of length `6/7`.
For a translated interval `J` of length `ell` and a nonzero integer `k`,

```text
Fourier(1_J)(k)
 =phase(J,k)*sin(pi*k*ell)/(pi*k).                   (G1)
```

For `ell=5/7` or `6/7`, this coefficient is nonzero exactly when `7` does
not divide `k`. Hence every all-seven-unit exact relation address

```text
m=(m_0,...,m_8) in Lambda(w)
```

from (E3) satisfies

```text
product_(i=0)^8 Fourier(1_(J_i))(m_i)!=0.            (G2)
```

Since `m.w=0`, its torus character restricts to the constant character on
the scalar line `t |-> wt`. Thus the septimal bank consists of genuine
nonzero Fourier terms surviving line restriction. The relations produced
after saturation can have height above `526`, so this statement concerns
the exact interval product, not necessarily the finite `Q_264` support.

This does not prove positive safe mass. All relation frequencies restrict
to the same zero line frequency; their complex phases and sine signs can
cancel in the constant-term aggregate. The finite-field count also forgets
amplitude:

```text
|Fourier(1_J)(k)|=|sin(pi*k*ell)|/(pi|k|),
```

so denominators can make deep-address contributions very small. A positive
Jackson kernel need not reach the address height in (E4). What remains is a
sign-definite subsum, a unique exposed frequency, or an owner-labelled
noncancellation theorem.

The thirteen-adic bank is weaker analytically: being a thirteen-unit does
not prevent a coefficient from being divisible by seven, so an address can
land on a zero interval Fourier factor. Likewise, `chi_7(-1)=-1` permits
the tournament-valued quadratic-character orientation from THM-2294, but
that orientation does not retain the phases or amplitudes in (G1).

## 8. Fixed-section lift and exact stopping boundary

THM-2298's fixed original-row section doubles only the guard coefficient:

```text
(x_H,x_rest) |->(2x_H,x_rest).                       (H1)
```

Multiplication by two is invertible modulo both seven and thirteen. The
lift preserves mod-`p` row rank, bright/dark columns, address counts, and
the septimal support statement. It doubles every row-height bound. Thus
the basic address height (E4) becomes

```text
1453(p^n-1),                                         (H2)
```

while (F5) and (F8) become respectively

```text
78993(13^n-1),
207018(7^n-1).                                       (H3)
```

The connection and loss ledger is:

```text
source:
  THM-2275's primitive and visible crossing seeds, THM-2295/2298's
  bounded rank ladder, THM-2298's actual Jackson support, THM-2284's
  centered p-adic pivot extension, and the interval Fourier zero set at
  multiples of seven;

target:
  visible rational rank six, bounded mod-seven and mod-thirteen rank-six
  packets, sharp all-depth affine address banks, the rank-six blocker
  pivot atlas, and exact septimal support-dark versus bright alternatives;

map:
  restrict the rank argument to actual finite Fourier survivors, saturate
  by centered cancel-and-divide, then avoid the essential coefficient
  hyperplane arrangement;

preserved:
  exact integral relation equations, all nine labels, mod-p row rank,
  all p-adic address lifts, the degree-526 septimal support, and the
  fixed-section lift;

destroyed:
  the original six short rows, the height-526 cutoff after saturation,
  coefficient signs and amplitudes, Fourier phases, root character,
  owner chronology, and ancestry incidence;

sharp hostile controls:
  the rank-six parallel-pencil packet (B5), intrinsic carrier darkness,
  signed cancellation among nonzero septimal survivors, and amplitude
  decay at deep addresses;

needed sidecar:
  discharge the support-dark deletion branch or force a sign-definite
  visible aggregate, then identify it with the marked ancestry channel.
```

No scalar profile is excluded. LRC(14) remains open.

## 9. Exact reproduction

Run

```bash
python3 04-computation/lrc14_visible_rank_six_kakeya_lift_thm2301.py
python3 -O 04-computation/lrc14_visible_rank_six_kakeya_lift_thm2301.py
```

Both executions must match

```text
05-knowledge/results/lrc14_visible_rank_six_kakeya_lift_thm2301.out
```

byte-for-byte after LF normalization. The companion checks the saturation
recurrence, the general arrangement formula, the sharp rank-three,
rank-five, and rank-six constants, every normalized address in the sharp
rank-six models over `F_7` and `F_13`, all-depth fibre arithmetic, the
intrinsic Grassmannian union bound, and all scalar/original height ledgers.
QED.
