---
id: THM-2413
title: "Prime-index affine drift and twin-center summand--multiplicand weld"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED, with explicitly
  separated CITED consequences.
  A rational line through prime-index points is exactly a level set of
  H_(a,b)(n)=b p_n-a n, so A373813 is a minimum affine-drift-level-set
  cover. For slope two, H_(2,1) has one initial descent and is
  nondecreasing from n=2 onward; its plateau edges are exactly twin-prime
  pairs. The centers A014574 are therefore an additive fixed-gap pattern
  cut out by multiplicative atoms. The unique repeated plateau is the
  startup triple 3,5,7, and every later slope-two twin line has only two
  prime points. The reciprocal-center identity is exact; convergence uses
  the cited Brun theorem. Equal-parent deletion gives additive defect
  x-1 and multiplicative defect tau(x)-2, while prime valuations make the
  latter's nonuniform additive coordinate system exact.
source: codex-2026-07-26-prime-index-drift
depends_on:
  - THM-362-natural-operation-graph-shadows
related:
  - THM-361-product-sum-defect-normal-form
  - HYP-1966-pair-lens-twin-prime-slices
  - HYP-1994-twin-goldbach-necklace
external: >
  OEIS A014574 and A373813; Kominers--Mrazovic--Pomerance--Sole,
  "Lines in the Prime Number Graph" (2026),
  https://math.dartmouth.edu/~carlp/lineprime060126.pdf.
script: 04-computation/prime_index_drift_twin_center_thm2413.py
output: 05-knowledge/results/prime_index_drift_twin_center_thm2413.out
script_sha256: c72f49d17e314df356f8d493abdaef32264d7cc0fcc1519fba0fbb70bfd66037
output_sha256: 4fd91fd2d31d270e02c1a75228048d52c1d4968cd4d5fc5436e039f40992d2de
secondary_script: 04-computation/prime_drift_twin_center_weld_thm2413.py
secondary_output: 05-knowledge/results/prime_drift_twin_center_weld_thm2413.out
secondary_script_sha256: 35d0a623b87ae5f0fa5e651c7c956b2c9343b5e8f14862dd535c8aa3911edc48
secondary_output_sha256: bc612c8b6419f3752744d318d1880378230c39abdade7f29497ca606dab72287
secondary_scope: independent normal-run evidence; its Python assertions are not optimization-safe
hash_basis: working-tree bytes (LF)
---

# THM-2413 -- prime drift and the twin-center operation weld

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED, WITH CITED
CONSEQUENCES SEPARATED.**

Let `p_n` be the `n`-th prime and let

```text
P_n=(n,p_n).
```

The central object is not the unlabelled graph of the points. It is the
family of integer-valued affine drifts

```text
H_(a,b)(n)=b p_n-a n,          gcd(a,b)=1, b>0.      (1)
```

Different slopes are different gauges on the same sequence. Twin-prime
centers are exactly the zero-increment edges of the slope-two gauge.

## 1. A373813 is an affine-drift level-set cover

Every line containing two prime-index points has rational slope `a/b` in
lowest terms and an equation

```text
b y-a x=c
```

with integers `a,b,c`, `b>0`. Therefore

```text
P_n lies on that line
iff b p_n-a n=c
iff H_(a,b)(n)=c.                                   (2)
```

Conversely, every level set in (2) is the intersection of the prime plot
with that rational line. A line used only for one point may be replaced by
the horizontal line `y=p_n`, so singleton cover pieces also have the form
(2).

It follows that the minimum number `L(N)` of straight lines covering

```text
P_1,...,P_N
```

is exactly the minimum number of **truncated fibres**

```text
F_(a,b,c)^(N)={n<=N:H_(a,b)(n)=c}
```

whose union is `{1,...,N}`. Equivalently one may select global level sets
whose union contains the first `N` indices; later prime points hit by those
lines are irrelevant. This is the sequence A373813.

This is a set cover, not necessarily a partition: affine patches may overlap
and their index sets need not be intervals.

For `i<j<k`, three points are collinear exactly when

```text
(p_j-p_i)/(j-i)=(p_k-p_j)/(k-j),                    (3)
```

or equivalently when the second divided difference vanishes. At equally
spaced indices this becomes

```text
Delta_h^2 p_i=p_(i+2h)-2p_(i+h)+p_i=0.              (4)
```

Thus A373813 is a sparsest cover by setwise degree-at-most-one Newton
divided-difference patches (Gregory--Newton on equally spaced subfamilies).
It measures vanishing second difference at selected index sets; it is not
the multiplicand graph itself.

## 2. The slope-two drift

Put

```text
h_n=H_(2,1)(n)=p_n-2n.                              (5)
```

Its forward difference is

```text
Delta h_n=(p_(n+1)-p_n)-2.                          (6)
```

The initial values are

```text
h_1,h_2,h_3,h_4,h_5,h_6=0,-1,-1,-1,1,1.            (7)
```

The gap `p_2-p_1=1` creates the unique initial descent. For every `n>=2`,
both primes are odd and their gap is at least two, so

```text
Delta h_n>=0.                                       (8)
```

Moreover,

```text
Delta h_n=0
iff p_(n+1)-p_n=2
iff (p_n,p_(n+1)) is a twin-prime pair.              (9)
```

The twin primes are therefore exactly the plateau edges of a monotone
integer drift after the one startup defect.

## 3. Why the startup defect is the only slope-two triple

If one slope-two line contains prime-index points at three indices
`i<j<k` with `i>=2`, then

```text
h_i=h_j=h_k.
```

By monotonicity, `h` is constant at every intermediate index. In particular
there are two consecutive prime gaps equal to two, so some

```text
p, p+2, p+4
```

are all prime. One of these three numbers is divisible by three. Hence the
only possibility is

```text
(p,p+2,p+4)=(3,5,7).                                (10)
```

It occurs at indices `(2,3,4)` on

```text
y=2x-1.                                             (11)
```

The omitted index `1` cannot create another fibre: `h_1=0`, whereas
`h_n=p_n-2n` is odd for every `n>=2`. Thus the line of slope two through
`P_1` contains no later prime-index point.

Consequently:

```text
maximum slope-two occupancy = 3;

the unique three-point slope-two line is (11);

every later slope-two line carrying a twin edge has exactly two
prime-index points.                                  (12)
```

This is a precise propagation of the irregular beginning: the exceptional
gap from `2` to `3`, followed by the overlapping twin pairs `(3,5)` and
`(5,7)`, creates the only repeated slope-two plateau. It does not imply
that all later prime irregularity is caused by addition.

## 4. A014574 as an additive pattern cut out by multiplicative atoms

Let

```text
c=(p_n+p_(n+1))/2=p_n+1
```

on a plateau edge. Then

```text
p_n=c-1,                     p_(n+1)=c+1.            (13)
```

THM-362 identifies primes as the edge labels/atoms in the transitive
reduction of divisibility. Hence the twin-center set has the intrinsic
description

```text
C_twin={c>=3: c-1 and c+1 are multiplicative atoms}. (14)
```

The additive relation supplies the fixed gap

```text
(c-1)+2=c+1
```

and the midpoint, while multiplicative atomicity filters the two endpoints.
This is exactly A014574:

```text
4,6,12,18,30,42,60,72,102,108,138,150,... .         (15)
```

Each center has the simultaneous weld

```text
(c-1)+(c+1)=2c,                                     (16)

(c-1)(c+1)=c^2-1.                                   (17)
```

For `c>=4`, this has the equivalent divisor-lattice form

```text
c is a twin center
iff Div(c^2-1)={1,c-1,c+1,c^2-1}.                  (17a)
```

The forward direction is the product of two distinct primes. Conversely,
an integer with four divisors is either a product of two distinct primes or
a prime cube. In the prime-cube case, two middle divisors at distance two
force `p^2-p=2`, hence `p=2` and `c=3`, excluded by `c>=4`. Thus the interval
in (17a) is a Boolean divisor diamond exactly at the twin centers.

The sum and product also have

```text
S=2c,             P=c^2-1,             S^2-4P=4.    (17b)
```

The discriminant reconstructs the adjacent pair `c-1,c+1` from `(S,P)`;
it does **not** by itself prove primality. The atomic/Boolean-diamond
condition in (17a) is the missing coordinate.

Thus the same centered pair feeds a labelled summand edge and a labelled
multiplicand edge. The simple shadows destroy the common center and the gap
label, so the operation cospan is essential.

For every pair except `(3,5)`, the lower prime exceeds three. Reduction
modulo six then forces the endpoints to be `6m-1,6m+1`. Therefore

```text
c=4, or 6 divides c.                                (18)
```

The centers `4` and `6` are the two edges of the unique triple (10).
The center set in the open HYP-1994 twin-Goldbach necklace is therefore
precisely A014574; this identification does not prove that hypothesis.

## 5. Reciprocal support identity

For each twin center `c`,

```text
1/2(1/(c-1)+1/(c+1))

=c/(c^2-1)

=1/c+1/(c(c^2-1)).                                  (19)
```

Let

```text
B_2=sum_(twin pairs (p,p+2)) (1/p+1/(p+2))
```

be Brun's twin-prime constant, with overlapping pairs counted separately.
Summing (19) gives the exact identity

```text
sum_(c in C_twin) 1/c

=B_2/2-sum_(c in C_twin) 1/(c(c^2-1)).              (20)
```

The correction series converges absolutely by comparison with
`sum c^(-3)`. The finiteness of `B_2` is Brun's theorem, so the convergence
of the reciprocal-center support in (20) is a **CITED consequence**, not a
new proof of Brun's theorem.

This realizes the “reciprocal of a sequence as a harmonic subset” lens
without conflating value sums and index sums.

## 6. Why twin centers do not control the full line cover

Every isolated pair of points determines a line. Therefore a later
twin-prime edge, whose slope-two line has occupancy only two by (12), gives
no automatic saving in a minimum line cover.

The prime plot also has richer non-twin affine patches. For example, the
eight points with indices

```text
6,7,10,12,13,16,18,21
```

lie on

```text
y=4x-11.                                            (21)
```

Thus A014574 is an exact slope-two subcarrier of A373813, but it cannot by
itself explain the global sublinear behavior of `L(N)`.

The 2026 paper *Lines in the Prime Number Graph* proves the **CITED**
statements

```text
L(N)=o(N);

#{n<=N: L(n)>L(n-1)}=L(N);

sum_(L(n)>L(n-1)) 1/p_n < infinity.                 (22)
```

The support in the last sum is the set of “awkward primes.” It is a second,
different convergent harmonic subset: twin-center support comes from
slope-two plateaus, while awkward-prime support comes from jumps of the
global minimum cover.

## 7. Summand and multiplicand shadows

THM-362 proves

```text
x -> z by some positive summand   iff x<z,

x -> z by some nonunit factor    iff x properly divides z.         (23)
```

Both shadows are transitive:

```text
z=x+a, w=z+b  => w=x+(a+b);

z=xa,  w=zb   => w=x(ab).
```

But they are structurally different. The first is a total order with one
successor in its transitive reduction; the second is a partial order whose
cover edges have prime quotients. Forgetting the witness `a` or factor `a`
turns the additive graph into a nearly featureless order and hides the prime
labels of multiplication.

There is an exact way to expose the difference. Delete the edge whose two
parents are equal:

```text
x R_+^neq z       iff x<z and z!=2x,

x R_x^neq z       iff x properly divides z and z!=x^2.             (23a)
```

At a fixed root `x>=2`, count two-step paths whose two displayed edges
survive but whose composite direct edge is the deleted equal-parent edge.
Additively these are exactly the ordered positive solutions

```text
c+d=x,
```

so there are `x-1`. Multiplicatively they are exactly the ordered solutions

```text
cd=x,                     c,d>=2,
```

so there are `tau(x)-2`. Consequently

```text
x is prime
iff tau(x)-2=0
iff R_x^neq has no rooted two-step transitivity defect.             (23b)
```

This does not say addition causes primes. The same diagonal deletion
produces a uniform additive defect and a divisor-sensitive multiplicative
defect; primality is the zero set of the latter.

Unique factorization makes the “nonuniform base” precise:

```text
nu:N_(>=1) -> direct_sum_p N,

nu(ab)=nu(a)+nu(b),                 nu(a)=(nu_p(a))_p.              (23c)
```

Multiplication is coordinatewise addition on infinitely many prime axes.
Its unweighted Hasse rank is `Omega(a)=sum_p nu_p(a)`, while its ordinary
metric is `log a=sum_p nu_p(a)log p`; the axes have nonuniform weights.

Ordinary addition is nonlinear in these coordinates. If

```text
a=p^r u,      b=p^s v,             p does not divide uv,
```

then

```text
nu_p(a+b)=min(r,s),                         r!=s,

nu_p(a+b)=r+nu_p(u+v),                      r=s.                    (23d)
```

The equal-valuation wall is exactly where carries/cancellation live. This,
rather than a causal slogan, is the rigorous propagation of additive
boundary irregularity into the prime-coordinate description of
multiplication.

Multiplication can be unfolded as a family of additive rays

```text
0,x,2x,...,yx=xy,                                   (24)
```

but (24) has a different step and length for every `(x,y)`. There is no
single “uniform base” recovered by the unlabelled additive shadow. Equations
(14), (16), and (17) are the lawful intersection: an additive local pattern
selects two vertices, multiplicative atomicity decides whether they form a
twin pair, and the center is the needed sidecar.

## 8. Equality and failure boundaries

1. Monotonicity of `h_n` begins at `n=2`, not `n=1`; the descent in (7) is
   load-bearing.
2. A line cover may use singletons. The horizontal replacement in Section 1
   is required for the exact drift-cover equivalence.
3. Vanishing `Delta_h^2` is imposed only within a selected affine patch.
   The patches need not be consecutive or disjoint.
4. A twin edge is not a three-point collinearity certificate. Apart from
   `(3,5,7)`, its line has only its two defining points.
5. The reciprocal identity (19) is elementary. Convergence of the twin-prime
   half uses Brun's theorem.
6. The finite A373813 prefix below is a control, not an asymptotic proof.

## 9. Exact companion

The dependency-free exact companion:

- builds the first `17,984` primes;
- checks `141,600` rational-line/drift incidences using a determinant path;
- verifies `804` second-difference/collinearity controls;
- checks the slope-two derivative law through the first `10,000` primes;
- verifies `1,270` twin centers, their additive/multiplicative welds, and
  the mod-six boundary;
- finds the unique slope-two triple `(3,5,7)`;
- checks `200` exact reciprocal identities;
- retains the eight-point slope-four line as a hostile to slope-two
  sufficiency; and
- independently solves the minimum line-cover problem through `N=27`,
  reproducing

  ```text
  1,1,2,2,2,3,3,3,3,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,6,6,7.
  ```

Run

```bash
python3 04-computation/prime_index_drift_twin_center_thm2413.py
python3 -O 04-computation/prime_index_drift_twin_center_thm2413.py
```

Both outputs must byte-match the stored transcript. Every executable check
raises explicitly under optimized Python.

## 10. Independent hostile audit

The audit reconstructed every affine-drift fibre and slope-two boundary
without using the primary incidence routine. A separate SciPy/HiGHS set-cover
formulation reproduced the exact A373813 prefix through `N=27`; direct
divisor enumeration reproduced the Boolean-diamond and diagonal-deletion
defect laws; and exact rational arithmetic reproduced the center-reciprocal
identity. It retained the startup descent, `3,5,7`, discriminant-four
composite, and slope-four eight-point hostiles.

The primary normal and optimized transcripts byte-match the stored output.
The secondary operation-weld program independently matches its stored
normal-run transcript; because it uses Python `assert`, it is not evidence
for optimized execution. Brun convergence and the 2026 prime-line
asymptotics remain `CITED`, not computationally promoted.
