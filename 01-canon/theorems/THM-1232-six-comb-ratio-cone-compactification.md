---
id: THM-1232
title: The six-comb slow-gap ratio cone is compact — functional dual drift bounds d_3/c<84/5, and three/two/one-tooth component spans propagate this to d_4/c<2254/15, d_5/c<10857/8, and d_6/c<40747
status: PROVED analytic consequence of THM-1198 plus exact interval-component geometry.  Every putative six-comb cover lies in one universal normalized ratio box.  The carrier c and its phase orbit remain unbounded, so this is not universal six-comb noncoverage or LRC(14)
source: codex-2026-07-18-S78 functional-tail compactification
depends_on: [THM-1198, THM-1196]
related: [THM-1176, THM-1178, THM-1199]
script: 04-computation/lrc14_six_comb_ratio_compactification_thm1232.py
output: 05-knowledge/results/lrc14_six_comb_ratio_compactification_thm1232.out
---

# THM-1232 — every relative speed is uniformly bounded

Let a complete `c`-slow gap be covered by six strict danger combs with

```text
c<d_1<d_2<d_3<d_4<d_5<d_6.                           (1)
```

Then

```text
d_3/c < 84/5,
d_4/c < 2254/15,
d_5/c < 10857/8,
d_6/c < 40747.                                        (2)
```

Thus the last-speed coordinate cannot escape even by dragging its entire
prefix with it.  This is stronger than a fixed-prefix tail bound: all six
ratios lie in one compact box independent of `c` and the slow-gap address.

## 1. The functional drift bounds the third speed

Put

```text
x_i=d_i/c,                  L_i=6x_i/7.               (3)
```

THM-1198 proves that every six-cover obeys

```text
sum_i Pbar(L_i)>=1,                                    (4)
```

where `Pbar(L)<=7/36` everywhere and, for `L>3`,

```text
Pbar(L)=1/7+1/(7L)=1/7+1/(6x).                        (5)
```

If `x_3<=7/2`, the first inequality in (2) is immediate.  Otherwise all four
indices `i>=3` lie on the tail branch (5).  Order and the strict inequalities
in (1) give

```text
1 <= sum_i Pbar(L_i)
  < 2(7/36)+4(1/7+1/(6x_3))
  =121/126+2/(3x_3).                                  (6)
```

The strict sign is available twice: `d_1,d_2>c` exclude the unique global
equality slope `L=6/7`, and `d_4,d_5,d_6>d_3` strictly improve three tail
terms.  Since `1-121/126=5/126`, (6) gives

```text
5/126<2/(3x_3),             hence x_3<84/5.           (7)
```

This is the missing scale-free consequence of the functional `H`-drift.
The scalar harmonic inequality alone cannot produce (7).

## 2. A three-comb component span

Write `w_i=1/(7d_i)` for one danger-tooth width.  THM-1196 proves that a
connected component of `D_(d_5) union D_(d_6)` has span at most

```text
w_5+2w_6.                                             (8)
```

This is smaller than `3w_4`, while the safe gap between consecutive `d_4`
teeth is `6w_4`.  Therefore a component of the faster two-comb union cannot
bridge two distinct `d_4` teeth.  It follows that every connected component
of

```text
D_(d_4) union D_(d_5) union D_(d_6)                   (9)
```

contains at most one `d_4` tooth.  Components of the faster union meeting
that tooth extend by at most their full span on each side, so (8) gives

```text
span(9) <=w_4+2(w_5+2w_6)<7w_4=1/d_4.               (10)
```

After deleting the first three combs, THM-1198 leaves dual mass strictly
larger than

```text
1-3(7/36)=5/12.                                       (11)
```

Since its density is at most `7/6`, the normalized survivor length is
strictly larger than `5/14`; in the physical slow gap this is

```text
|E_3|>15/(49c).                                       (12)
```

The first three combs cut the slow gap into at most

```text
C_3<=1+sum_(i=1)^3 ceil((6d_i+c)/(7c))               (13)
```

components.  Under the six-cover, each lies in one component (9), so
(10)--(13) imply

```text
d_4/c<(49/15)C_3.                                    (14)
```

This is the `j=3` continuation of THM-1196's `j=4` flood tail.  The same
simple recursion stops one rung earlier: a three-comb faster component can
approach the `6w_3` gap closely enough that four remaining combs need not
contain only one `d_3` tooth.

## 3. Exact propagation through the ceiling ledger

From (7), for `i<=3`,

```text
(6d_i+c)/(7c)<509/35<15,
C_3<=1+3*15=46.                                      (15)
```

Equation (14) therefore gives

```text
d_4/c<(49/15)*46=2254/15.                            (16)
```

THM-1196's two-comb component-span law is

```text
d_5/c<(21/8)C_4,
C_4<=1+sum_(i=1)^4 ceil((6d_i+c)/(7c)).               (17)
```

Using (16), every ceiling in (17) is at most `129`, because

```text
(6d_i+c)/(7c)<4513/35<129.
```

Thus

```text
C_4<=517,                 d_5/c<(21/8)*517=10857/8. (18)
```

Finally THM-1196's last-tooth law gives

```text
d_6/c<7C_5,
C_5<=1+sum_(i=1)^5 ceil((6d_i+c)/(7c)).               (19)
```

Equation (18) makes every ceiling in (19) at most `1164`, since

```text
(6d_i+c)/(7c)<32575/28<1164.
```

Consequently

```text
C_5<=5821,                  d_6/c<7*5821=40747,      (20)
```

completing (2).

## 4. What compactness changes

The six-comb slow-gap problem no longer has a projective ratio tail.  Any
remaining sequence of covers with `c -> infinity` has, after passing to a
subsequence, convergent ratios and convergent phase-orbit data inside a fixed
compact torus box.  What remains is not to control arbitrarily separated
teeth, but to exclude (or classify) limiting real covers and then obtain a
uniform arithmetic margin for the cyclic phase orbits.

This does **not** make the problem finite in `c`.  The residues
`k(d_i mod c)/c`, common divisors, and endpoint-owner words can keep changing
at every scale.  A finite enumeration that silently drops those phase data
would not follow from (2).

## 5. Tournament and carrier audit

The runner-order tournament, oriented by `sign(d_j-d_i)`, is transitive with
score histogram `0,1,2,3,4,5`, no directed cycles, singleton SCCs, and one
Hamiltonian path.  It retains order but none of the argument above.

The faithful object is a two-layer obligation DAG:

```text
six phase-free load chambers -> third-ratio gate (84/5)
                              -> j=3/j=4/j=5 component spans
                              -> exact ceiling propagation.               (21)
```

We challenged runners, arcs, gaps, residues, beat points, endpoint events,
Fourier modes, Fano lines, and proof obligations as vertices.  The useful
vertices here are proof obligations carrying both a metric constant and a
prefix depth.  Quotienting to runners destroys the nonlinear envelope;
quotienting to load chambers destroys component chronology.  Both layers are
necessary.

## 6. Reproducibility and frontier

The companion exact-rational script checks every constant and strict ceiling
transition in (6)--(20) under ordinary and optimized Python.  It is an
arithmetic consumer, not a substitute for THM-1198's arrangement theorem or
the interval-component proofs.

Normal and optimized outputs are byte-identical to the stored ledger.  Frozen
SHA-256 hashes are

```text
source  67896239d77767ae7111f7a7fbe165c33fdfcd47adaeeb9c1f617a9e489528bd
output  18abfa432c5f67c631c5a9fb74776ff7f1bc20788255f0f57f98f781ebf41836
```

The result leaves three honest tasks: classify the limiting real six-covers
inside (2), recover a positive uniform margin on the arithmetic phase orbits,
or compose the compact box with the existing finite-carrier BV certificates.
LRC(14) remains open.
