# Critical seven-tail grids force a bounded six-component scale

**Status: PROVED grid count and scale reduction, with
[independent analytic and exact audit accepted](third_20260906_grid_audit.md).**
The direct conclusion is weak full-row safety from an explicit grid count.
Combining it with inherited subset-gcd caps and an actual inert decoder edge
gives a universal finite upper bound on the smaller component's physical
scale in a hypothetical balanced failure. LRC(14) remains OPEN.

## Inheritance and connection

[THM-2060 / sharp CRT tail-coset saturation](../../01-canon/theorems/THM-2060-crt-tail-coset-saturation.md)
and [THM-2064 / multitail sheet capacity](../../01-canon/theorems/THM-2064-multitail-sheet-capacity-and-dyadic-seam.md)
already give the exact marginal ceiling count. The new step retains the
overlaps those marginal counts discard at the critical seven-tail budget.
The [current joint-shadow caps](lrc14_joint_shadow_empty_core_next_sep06.md)
supply necessary gcd bounds in a strict failure. The full decoder domain
and ratio atlas are inherited from
[THM-3818 / inert pair packet](../../01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md).
Proper six-component LRC is the same CITED input used in the existing entry
theorems. No new literature claim is made.

The hostile is a selected grid killed by a tail divisible by its scale;
coprimality cannot be silently imported from the global primitive gcd.
The corrected near miss is treating a contained interval intersection as a
triangular overlap: the corrected microscopic part of
[THM-739 / exact Bernoulli overlap](../../01-canon/theorems/THM-739-pairwise-coprime-bad-overlap-exact-bernoulli-closed-form.md)
retains its smaller-width cap. Its main Bernoulli identity remains valid.
The least-used sidecar is the **common origin of all danger sets before
projecting to their separate sheet orders**.

The board is: actual translated grid; individual sheet gcds; hereditary
divisor capacities; overlap multiplicity; primitive ratio components; and
the resulting physical scale. No tournament is imposed on symmetric danger
incidences. A failed sufficient count is not an unsafe row.

## 1. A direct count for every translated grid

For positive integers `u_1<...<u_b`, integer `t>=1`, and real alpha, let

```text
X_alpha={(alpha+j)/t mod1: j=0,...,t-1},
D_u={x: ||ux||<1/14},       d_i=gcd(t,u_i),
B(t,U)=sum_i d_i ceil(t/(7d_i)).
```

THM-2060 gives `|X_alpha intersect D_(u_i)|<=d_i ceil(t/(7d_i))`.
At a grid point with danger multiplicity r, replacing the sum of indicators
by their union saves exactly `max(r-1,0)`. Therefore any lower bound C on
that total multiplicity excess gives

```text
number of weak-safe points >= t-B(t,U)+C.             (1)
```

Two alternative overlap credits are useful. They must not simply be added
without proving their combined pointwise multiplicity bound.

**Origin credit.** For `i>=2`, the open circle interval
`||x||<1/(14u_i)` lies in the first i danger sets. These intervals are
nested. Their indicator sum is at most the danger multiplicity excess.
An open interval of length ell contains at least `ceil(t ell)-1` points
of every translated t-grid. Hence

```text
C_0(t,U)=sum_(i=2)^b (ceil(t/(7u_i))-1)               (2)
```

is a valid credit, uniformly in alpha. This retains every small-label layer,
not only a selected pair.

**Pair credit.** Choose two labels `u=d p`, `v=d q`, with coprime `p<q`.
Let `mu(p,q)` and `J(p,q)` be the measure and number of open components of
`D_p intersect D_q`, and put `e=gcd(t,d)`. Multiplication by d sends the
t-grid to a translated `(t/e)`-grid with multiplicity e. Counting each open
component gives

```text
|X_alpha intersect D_u intersect D_v|
 >= e[ceil(t mu/e)-J] >= t mu-eJ.                    (3)
```

The first lower bound follows because the sum of ceilings is at least the
ceiling of the sum. A negative lower bound is harmless and can be replaced
by zero. This overlap is one valid multiplicity-excess credit in (1).
More generally, a forest on the b labels permits summing its pair credits:
the induced forest on r dangerous vertices has at most r-1 edges.

Now take a primitive thirteen-speed row in the form

```text
tV union gU, |V|=6, |U|=7,
gcd(V)=gcd(U)=gcd(t,g)=1.
```

For any supplied V-safe phase eta, the t lifts `(eta+j)/t` preserve its
clearance. In the U clock their images form exactly a translated t-grid
because gcd(t,g)=1. Thus any strictly positive lower bound in (1) gives a
common physical phase of clearance at least1/14. This implication does not
require decoder equality. It is a pointwise existence argument, not an
average safe-mass assertion. Strict clearance is not inferred from weak
grid endpoints.

## 2. Hereditary caps bound the critical ceiling excess

In a hypothetical strict failure, the gcd of a subset of sizes
7,8,9,10,11,12 is at most90,30,9,4,2,1 respectively. Apply this to all six
physical V labels together with k chosen physical U labels. Their gcd is
exactly the gcd of the corresponding `d_i=gcd(t,u_i)`. Consequently

```text
gcd(d_i: i in S) <= C_k,
(C_1,...,C_7)=(90,30,9,4,2,1,1).                    (4)
```

At seven tails the critical excess is an integer

```text
E(t,U)=B(t,U)-t >=0.
```

For any divisor h of some d_i, (4) bounds the number of d_i divisible by h.
In particular the multiplicity of an individual value d is at most
7 for d=1, 5 for d=2, 4 for d=3,4, 3 for d=5,...,9,
2 for d=10,...,30, and 1 for d=31,...,90.
Using only these individual multiplicities gives a valid relaxation; the
complete subset caps remain part of the original domain.

If `7` does not divide t, put `tau=t mod7`. Since every d_i divides t,

```text
7E=sum_i d_i*((-tau*d_i^(-1)) mod7).                (5)
```

Taking the seven largest weights in the finite multiplicity-relaxed bag
gives the following sharp upper bounds for the abstract subset-cap problem:

| tau | E bound | One attaining seven-value gcd row |
|---:|---:|---|
| 1 | 429 | 85,78,88,71,81,64,74 |
| 2 | 426 | 86,79,72,85,78,65,71 |
| 3 | 438 | 87,80,89,73,82,66,75 |
| 4 | 435 | 88,81,74,86,67,79,90 |
| 5 | 447 | 89,82,90,75,83,68,76 |
| 6 | 445 | 90,83,76,87,69,80,62 |

Every displayed row satisfies all127 nonempty subset conditions (4).
Taking t to be a suitable multiple of the row's lcm realizes the indicated
residue tau and all its sheet gcds. This establishes sharpness for that
abstract gcd-cap object, not for actual inert connectivity, the full
joint-shadow word, or strict failure.

If `h=v_7(t)>=1`, only d_i with `v_7(d_i)=h` contribute to E. For h=1,
(4) permits at most three such labels; the three largest possible
contributions give bounds151,134,152,135,146,122 according to
`(t/7) mod7`. These are upper bounds; no joint attainment is claimed.
For h=2, only d_i=49 is possible and it occurs at most once, so E<=42.
For h>=3 no d_i<=90 carries the full valuation, so E=0. Therefore

```text
E(t,U)<=447                                             (6)
```

in every strict failure. If U contains1, fixing one d_i=1 improves the
nonseptimal bounds to377,376,385,384,393,392, and hence uniformly E<=393.
For example (2) then supplies weak safety whenever `t>2758*u_2`.
Without a unit, the corresponding simple condition is `t>3136*u_2`.

## 3. The actual ratio atlas gives a universal scale ceiling

The corrected exact primitive-pair geometry gives

```text
J(p,q)=2 ceil((p+q)/14)-1,
mu(p,q)=[p+sum_(k=1)^(ceil((p+q)/14)-1)
                      min(2p,p+q-14k)]/(7pq).      (7)
```

Equivalently mu is the inherited Bernoulli formula. A useful coarse control
is the universal minimum `mu>=1/91`: for `pq>=27` it follows from
`mu>=1/49-1/(4pq)`, and the36 coprime pairs with `pq<=26` are exact controls,
with equality at `(1,13)`. This coarse minimum is not advertised as sharp
on the smaller inert atlas.

Suppose U contains a pair whose primitive coefficient sum is at most356
and has only inert primes, each with exponent at most two. This hypothesis
is supplied, in particular, by every actual THM-3818 balanced decoder entry:
its larger primitive component is connected by such edges. The present
grid implication itself needs neither decoder equality nor a physical box.
Select any such pair. Its pair e in (3) satisfies `e<=30` by the eight-subset cap.
Equations (1), (3), and (6) guarantee weak safety if

```text
t mu(p,q)>447+30J(p,q).                              (8)
```

The complete actual atlas consists of5,855 coprime pairs over94 allowed
sums. Exact rational maximization gives

```text
max (447+30J)/mu = 6019965/62,
unique maximizing pair (p,q)=(5,348),
J=51, mu=62/3045.                                   (9)
```

Thus every such primitive balanced row with **t>=97097** is weakly safe.
Equivalently, a hypothetical strict failure in this branch has

```text
1<=t<=97096,       1<=g<=90.                         (10)
```

This is a uniform bound on the six-component physical scale, independent
of both primitive component maxima and of the presence of a unit or a
parity phase. It is not a finite enumeration of all remaining component
shapes, and it does not assert every small scale is feasible or unsafe.
The proof uses the actual connected ratio edge, full physical phase clock,
and necessary gcd caps. No arbitrary rank-eleven matrix is substituted.

## 4. Verification and next exact test

The [independent audit](third_20260906_grid_audit.md) reconstructs the
clipped overlap geometry, the multiplicity relaxation, and the full inert
atlas, including both boundary conventions. The standalone
[audit source](../../04-computation/third_20260906_grid_audit.py) and
[frozen output](third_20260906_grid_audit.out) record **303,923 always-active
gates**; normal and optimized output agree byte for byte. The source checks
all5,855 actual ratio pairs, all94 sums, literal arc sweeps and translated
grid cells, every abstract attaining row, and the repaired clipping hostiles.
The analytic proof supplies the all-parameter statements; the exhaustive
finite atlas supplies its exact maximum.

The finite scale bound makes a new refinement possible: retain the divisibility
`d_i|t` in the finite cost bag rather than allowing all d up to90. The same
applies to e. This can only strengthen the safety test. The [completed full-word refinement](third_20260906_grid_refined.md)
now records8,202 remaining scales, all at most16,704. A surviving scale
is still an uncertified quotient object, not an unsafe row.
