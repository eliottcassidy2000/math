---
id: THM-2310
title: "Quantitative beta-shift gluing of positive handoffs"
status: >
  PROVED CANDIDATE UNDER INDEPENDENT AUDIT. For the integer beta-shift, two
  positive finite unions of intervals meet after r iterates whenever
  12 beta^r mu(A)mu(B)>Var(1_A)Var(1_B). Consequently every finite path of
  positive finite-interval transition carriers can be realized, in order, by
  one orbit after finite explicitly computable waits. In particular a pure
  directed 2-cycle or 3-cycle in THM-2305 gives a single-orbit closed label
  itinerary, although not a periodic orbit. Thus disjoint incoming/outgoing
  owner subsets are not the residual obstruction at support level. The waits
  alter the prescribed clocks and do not supply relative Fourier phase, so
  LRC(14) remains open.
source: codex-2026-07-25-quantitative-handoff-gluing
depends_on:
  - THM-2303-terminal-component-phase-current-and-defect-rank
  - THM-2305-canonical-blocker-word-handoff-hypergraph
related:
  - THM-2075-safe-child-homeomorphism-and-wall-word-conjugacy
  - THM-2077-terminal-kakeya-needle-and-recursive-quarter-escape
  - THM-2288-shallow-owner-bv-mixing-and-delayed-blocker-handoff
  - THM-2296-prescribed-expiration-return-or-bounded-ancestry-resonance
---

# THM-2310 -- quantitative beta-shift gluing of positive handoffs

**PROVED CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2305 isolates the honest terminal relation. A pure word gives a
time-oriented edge between exclusive owners, while the simultaneous word gives
a two-head hyperedge. Even if pure edges form a directed label cycle, their
positive arrival and departure subsets inside the intermediate owner can be
disjoint.

At the level of support, that is not a permanent obstruction. The expanding
circle map is quantitatively mixing, and finite scalar carrier sets have
bounded variation. Positive handoffs can therefore be connected after bounded
waiting times. What the wait does *not* preserve is exactly what THM-2303 says
is missing: the prescribed component phase current.

## 1. A sharp-form bounded-variation mixing estimate

Let

```text
T_beta(x)=beta x mod 1,                    beta>=2,  (1)
```

and let `A,B` be positive-measure finite unions of intervals in the circle.
Endpoints are understood modulo null sets. Write

```text
V_A=Var(1_A),          V_B=Var(1_B).                 (2)
```

For an interval union whose essential boundary has `J` points, `V=J`.

> **BV beta-mixing lemma.** For every integer `r>=0`,
>
> ```text
> |measure(A intersection T_beta^(-r)B)
>      -measure(A)measure(B)|
> <=V_A V_B/(12 beta^r).                             (3)
> ```
>
> In particular,
>
> ```text
> 12 beta^r measure(A)measure(B)>V_A V_B             (4)
> ```
>
> forces `measure(A intersection T_beta^(-r)B)>0`.

### Proof

Use the Fourier convention

```text
f_hat(n)=integral_T f(x)exp(-2*pi*i*n*x)dx.          (5)
```

Integration by parts for a circle BV function gives

```text
|f_hat(n)|<=Var(f)/(2*pi*|n|),          n!=0.        (6)
```

In `L^2`,

```text
1_B(T_beta^r x)
 =sum_(n in Z) 1_B_hat(n)exp(2*pi*i*n*beta^r*x).
```

The nonconstant correlation series is absolutely convergent by (6), and

```text
measure(A intersection T_beta^(-r)B)
 =sum_(n in Z) 1_A_hat(-beta^r n)1_B_hat(n).         (7)
```

The `n=0` term is `measure(A)measure(B)`. Therefore

```text
|sum_(n!=0) 1_A_hat(-beta^r n)1_B_hat(n)|
 <=V_A V_B/(4*pi^2 beta^r)
      sum_(n!=0)1/n^2
 =V_A V_B/(12 beta^r),                              (8)
```

because `sum_(n!=0)n^(-2)=pi^2/3`. This proves (3), and
(4) makes its lower bound strictly positive. QED.

The least `r` satisfying (4) is an explicit wait. The constant `1/12` is the
direct product of the two universal BV coefficient bounds; no claim is made
that it is the optimal mixing constant for indicator sets.

## 2. Transition carriers and arrival support

A **finite-interval transition carrier**

```text
(F,k): X -> Y                                          (9)
```

consists of a positive-measure finite union of intervals `F`, an integer
`k>=0`, and the containment

```text
F subset X intersection T_beta^(-k)Y                 (10)
```

up to null sets.

Let `P_beta` be the normalized Perron operator. Define the essential arrival
support

```text
Arr(F,k)={y:P_beta^k 1_F(y)>0}.                       (11)
```

It is a finite union of intervals contained in `Y`. Every positive value of
`P_beta^k 1_F` is at least `beta^(-k)`, so

```text
measure(Arr(F,k))>=measure(F)>0.                     (12)
```

Every jump of `P_beta^k1_F` is an image of a jump of `1_F`. Positivity can
change only at those jumps. Hence

```text
Var(1_(Arr(F,k)))<=Var(1_F).                         (13)
```

These two elementary facts retain enough information to glue carriers.

## 3. Two positive handoffs always compose after a bounded wait

Suppose

```text
(F,k): X -> Y,
(G,l): Y -> Z                                       (14)
```

are positive finite-interval transition carriers. Put

```text
A=Arr(F,k),       p=measure(F),       q=measure(G).
```

Choose any `r>=0` satisfying

```text
12 beta^r p q
 >Var(1_F)Var(1_G).                                 (15)
```

Equations (3), (12), and (13) give

```text
B=A intersection T_beta^(-r)G,
measure(B)>0.                                       (16)
```

Lift this overlap back to the first carrier:

```text
F' =F intersection T_beta^(-k)B.                    (17)
```

The transfer identity and the pointwise lower bound on the arrival density
give

```text
measure(F')
 =integral_B P_beta^k1_F
 >=beta^(-k)measure(B)>0.                            (18)
```

Every `x in F'` obeys

```text
x in X,
T_beta^k x in Y,
T_beta^(k+r)x in G subset Y,
T_beta^(k+r+l)x in Z.                               (19)
```

Thus the two handoffs occur along one orbit, separated by the bounded wait
`r`. Notice that the orbit may leave `Y` and return to it during the wait;
only the named checkpoints in (19) are asserted.

## 4. Every finite positive path has a single-orbit realization

Let

```text
(F_i,k_i):X_(i-1)->X_i,              1<=i<=m,        (20)
```

be a finite sequence of positive finite-interval carriers. Repeatedly apply
Section 3.

More explicitly, after a partial path there is a positive finite-interval
set of reachable points in `X_i`. Mix that set with `F_(i+1)`, restrict to
the positive overlap, and lift through the preceding finite-to-one map as in
(17)--(18). Images, preimages, and finite intersections of interval unions
remain finite interval unions. Induction therefore yields nonnegative waits

```text
r_1,...,r_(m-1)                                    (21)
```

and a positive-measure subset of `F_1` whose points traverse every carrier in
the prescribed order.

At each induction step the least wait satisfying (4) for the current
reachable set and the next carrier is computable from their interval
endpoints. The reachable mass can shrink, so this statement does not replace
those recursive data by a fictitious uniform short wait.

If the labels in (20) form a directed cycle, append `F_1` as the final target
set. The same argument gives a positive set of initial points whose orbit
traverses every cycle edge and later lands back in `F_1`. This is a
**single-orbit closed label itinerary**. It is not a claim that the endpoint
equals the starting point, that the orbit is periodic, or that the itinerary
can be repeated indefinitely.

## 5. Application to the LRC blocker-word graph

Use THM-2305's scalar notation with

```text
T=T_13.
```

A positive pure word from owner `j` to owner `t` is exactly the carrier

```text
F_(j,t)
 =E_j intersection T^(-k_j)E_t,       k_j=lambda_j+1. (22)
```

It has mass `rho_(j,{t})>0`. All sets are Boolean combinations of finitely
many scalar comb intervals, so they satisfy the hypotheses above.

If the conditional pure-edge graph in THM-2305 contains

```text
j <-> t
```

or

```text
1 -> 2 -> 3 -> 1,                                   (23)
```

Sections 3--4 turn (23) into a single-orbit closed label itinerary after
finite waits. Thus the possible disjointness of an arrival subset and the
next departure subset is no longer an obstruction to **support-level**
composition.

There is also a direct row-dependent wait bound. If `S` is the scalar speed
sum and an edge has clock `k`, then

```text
Var(1_(F_(j,t)))<=2S(1+13^k).                       (24)
```

Indeed `E_j` has at most `2S` boundary points, while the preimage under
`T^k` of the at-most-`2S` boundary points of `E_t` has at most
`2S*13^k` points. For two consecutive pure carriers of masses `p,q` and
clocks `k,l`, it is enough to choose `r` with

```text
12*13^r*p*q
 >4S^2(1+13^k)(1+13^l).                             (25)
```

This bound is enormous but finite under the existing global scalar ceiling.
It is a structural certificate, not a practical row elimination.

The simultaneous word

```text
j -> {a,b}
```

remains a genuine hyperedge. Mixing cannot choose one head without losing the
exact terminal word, so the fork branch of THM-2305 is untouched.

## 6. The word-overlap lens and its quantitative role

Zakharov's cited word-overlap theorem gives a symbolic companion. For
length-`n` word families `A_n,B_n` over a finite alphabet, absence of every
directed suffix--prefix overlap implies

```text
mu(A_n)mu(B_n)
 <=(1/n)(n/(n+1))^(n+1).                            (26)
```

For base-`beta` cylinder unions, a suffix--prefix overlap of length `j` is
exactly an intersection after the shift `r=n-j`. Thus (26) is another
finite-wait criterion.

Actual LRC interval endpoints are generally not beta-adic cylinders. If
`A^(n)` is the union of depth-`n` cylinders contained in `A`, then

```text
measure(A^(n))>=measure(A)-V_A/beta^n,              (27)
```

and similarly for `B`. Consequently Zakharov's theorem applies lawfully once

```text
(measure(A)-V_A/beta^n)_+
(measure(B)-V_B/beta^n)_+
 >(1/n)(n/(n+1))^(n+1).                             (28)
```

Since the left side tends to a positive product and the right side tends to
zero, some finite `n` always works. For interval unions, however, the Fourier
estimate (3) gives the much shorter logarithmic-in-`1/(mu(A)mu(B))` wait.
The word theorem contributes the correct suffix--prefix representation and
first-hit perspective; it is not needed as an external dependency of the BV
proof.

## 7. Exact stopping boundary

The gluing map preserves

```text
positive measure,
source and target labels at named checkpoints,
the order of a finite pure-edge path,
and finite interval complexity.                    (29)
```

It destroys or fails to control

```text
the original consecutive expiration clocks,
the root sheet and address,
the exact common Fourier index,
the continuous terminal base phase,
and every phase-tree transport ratio.               (30)
```

Inserting a wait `r` replaces a prescribed spectral channel by a
`13^r`-shifted divisibility problem. Support mixing therefore cannot repair
THM-2303's rank-one phase defect or THM-2302's missing unit-colored marked
degree.

The sharpened frontier is:

```text
pure-edge support graph:
  composable after bounded waits;

fork branch:
  still needs a two-blocker intersection mechanism;

signed Fourier proof:
  still needs phase transport at the prescribed clock.                (31)
```

No owner is forced onto the pure branch, the conditional three-owner
hypothesis of THM-2305 is not supplied, no scalar profile is excluded, and
LRC(14) remains open. QED.
