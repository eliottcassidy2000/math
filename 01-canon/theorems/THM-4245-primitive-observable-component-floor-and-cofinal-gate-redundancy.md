---
id: THM-4245
title: "Primitive-observable danger-overlap component floor and cofinal gate redundancy"
status: >
  PROVED ANALYTIC ROUTE OBSTRUCTION + VERIFIED-EXACT + INDEPENDENTLY
  HOSTILE-AUDITED. For every distinct actual pair a,b, the THM-4233
  primitive-observable gate forces a+b-gcd(a,b)>=2930, hence a+b>=2931 and
  max(a,b)>=1466. Therefore the entire scalar gate on two fixed-pool
  outsiders is already covered by THM-4231 and cannot remove a finite
  residual edge. The gate is nonvacuous globally; failure is not Haar danger,
  and LRC(14) remains OPEN.
source: root/cross-frontier-bridge/2026-08-26
depends_on:
  - THM-4233-pair-specific-primitive-observable-oscillation-haar-charts
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
related:
  - THM-4228-common-gcd-two-outsider-periodic-observable-haar-ray
  - THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer
  - THM-4238-forty-small-label-uniform-r590-haar-tail-closure
  - THM-4242-fixed-fifty-direct-r590-tail-and-twenty-three-label-chart
primary_script: 04-computation/lrc14_primitive_gate_danger_overlap_height_grid_audit_thm4245.py
primary_output: 05-knowledge/results/lrc14_primitive_gate_danger_overlap_height_grid_audit_thm4245.out
independent_audit_script: 04-computation/lrc14_primitive_gate_danger_overlap_height_teeth_audit_thm4245.py
independent_audit_output: 05-knowledge/results/lrc14_primitive_gate_danger_overlap_height_teeth_audit_thm4245.out
primary_script_sha256: 32aa2d64d8c181c936b5677d778f5f4e54428add0e6ae075ce22b212a1fd7f0e
primary_output_sha256: b30d8a34962ef1944b89b9a97b3b55717f266360293c48565272563a784340fa
independent_audit_script_sha256: 0ddcceff34cb261e702ae51e209ca20144a4ea6460e5a4bad7626087f9ab93b8
independent_audit_output_sha256: 1878be50ae7f913797e37a172fffba91966305c2f5d6d22d47e305a4b728bb7b
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT AFTER STRENGTHENING AND ONE TOPOLOGY REPAIR. Hostile review
  sharpened the danger-overlap and component bounds, then caught that isolated
  safe points prevent equality of safe/danger component counts. The repaired
  exit-boundary injection proves the needed inequality with unchanged
  constants. Two independent residual reconstructions, each under normal and
  optimized Python, check all 115,429 primitive ratios, the nonvacuous
  (3713,5149) gate control, zero residual gate passes, and exact hashes.
---

# THM-4245 -- primitive-observable danger-overlap component floor and cofinal gate redundancy

**PROVED ANALYTIC ROUTE OBSTRUCTION + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED; LRC(14) REMAINS OPEN.**

## 1. Definitions and theorem

For a positive integer `n`, put (endpoint choices are null and immaterial)

```text
D_n={x in R/Z: ||nx||<1/14},
G_n=(R/Z)\D_n={x in R/Z: ||nx||>=1/14}.                 (1)
```

For distinct positive integers `u,v`, retain exactly the THM-4233 joint-pair
observable

```text
A_(u,v)=G_u intersect G_v,
beta_(u,v)=mu(A_(u,v)),                                 (2)

H_(u,v)(t)=integral_0^t(1_(A_(u,v))(s)-beta_(u,v)) ds,
omega_(u,v)=max H_(u,v)-min H_(u,v).                   (3)
```

The value of `H` depends on a choice of origin, while `omega` does not.  Let

```text
T=1650/28710227=(1650/8281)/3467.                       (4)
```

This is the exact THM-4233 depth-eight transfer threshold.  Given distinct
positive integers `a,b`, set

```text
g=gcd(a,b),       u=a/g,       v=b/g.                  (5)
```

> **Theorem (height obstruction).** If the THM-4233 pair gate
>
> ```text
> beta_(u,v)>=66/91,          omega_(u,v)/g<=T          (6)
> ```
>
> holds, then
>
> ```text
> a+b-g>=2930,          a+b>=2931,          max(a,b)>=1466. (7)
> ```

The quantifiers are universal over every distinct positive pair `a,b`; the
statement is not restricted to the displayed THM-4233 families or to the
fixed pool.

## 2. Component-oscillation lemma

Let `A` be a nonempty proper finite union of `c` positive-length circular
intervals, let `beta=mu(A)`, and put

```text
H(t)=integral_0^t(1_A(s)-beta) ds,
omega=max H-min H.                                     (8)
```

Then

```text
c omega >= beta(1-beta).                               (9)
```

Indeed, choose the origin in the positive-length complement of `A`; this only
adds a constant and a translate to `H`, so it preserves `omega`.  Write the
component lengths of `A` as `ell_1,...,ell_c`.  On component `i`, `H` has
slope `1-beta`, hence rises by

```text
(1-beta) ell_i.                                        (10)
```

Every such rise is at most the global range `omega`.  Summing `(10)` and
using `sum_i ell_i=beta` gives

```text
beta(1-beta)=sum_i (1-beta)ell_i <= c omega,            (11)
```

which proves `(9)`.  This is a lower bound; THM-4228's familiar
`omega<=beta(1-beta)` is the complementary upper bound.

## 3. Primitive danger-overlap and component bounds

The following estimates are universal for every two distinct coprime positive
integers `u,v`:

```text
mu(D_u intersect D_v)<=1/14,
beta_(u,v)<=11/14,
c_(u,v)<=u+v-1.                                        (12)
```

Relabel the pair so that `v=max(u,v)>=2`.  The `v` danger components of
`D_v` have the almost-everywhere disjoint parametrization

```text
x=(j+y)/v,       j in Z/vZ,       -1/14<y<1/14,         (13)
```

with `dx=dy/v`.  Fix `y`.  The condition that `(j+y)/v` also belongs to
`D_u` is

```text
||u(j+y)/v||<1/14.                                     (14)
```

Because `gcd(u,v)=1`, multiplication by `u` permutes `j mod v`.  Thus, as
`j` varies, the values in `(14)` are the `v` equally spaced points

```text
{k/v+(u/v)y mod 1 : k in Z/vZ}.                        (15)
```

The circular danger arc `||t||<1/14` has length `1/7`.  Any open circular arc
of that length contains at most `ceil(v/7)` points of a shifted `1/v` lattice:
after cutting at a point outside the finite lattice, `m` contained points span
`m-1` spacings strictly shorter than the arc, which gives
`m<=ceil(v/7)`.  Moreover,

```text
ceil(v/7)<=v/2                                         (16)
```

for every `v>=2`: check `2<=v<=7` directly, while for `v>=8`,
`ceil(v/7)<=v/7+6/7<=v/2`.  If `N(y)` denotes the number of indices in
`(14)`, Fubini applied to `(13)` now gives

```text
delta_(u,v):=mu(D_u intersect D_v)
 = (1/v) integral_(-1/14)^(1/14) N(y) dy
 <= (1/v)(1/7)(v/2)=1/14.                              (17)
```

Since every `D_n` has measure `1/7`, inclusion-exclusion yields

```text
beta_(u,v)
 =1-mu(D_u union D_v)
 =1-2/7+delta_(u,v)
 =5/7+delta_(u,v)<=11/14.                              (18)
```

The upper bound is attained by the primitive pair `(1,2)`: the zero danger
component of `D_2` is `(-1/28,1/28)`, contained in the zero component
`(-1/14,1/14)` of `D_1`, while the other `D_2` component is disjoint from
`D_1`.  Hence `delta_(1,2)=1/14` and `beta_(1,2)=11/14`.

Finally, `D_u` and `D_v` are unions of `u` and `v` positive-length circular
danger components, and their two zero components overlap on a neighborhood of
zero.  Their union therefore has at most `u+v-1` positive-length components.
Whenever the gate density holds, `A_(u,v)` is nonempty and proper.  Traversing
the circle, every positive-length component of `A_(u,v)` begins at the exit
boundary of a distinct positive-length component of `D_u union D_v`.
Isolated safe points can occur when danger components meet, but they do not
count toward `c_(u,v)`.  Hence

```text
c_(u,v)<=#components(D_u union D_v)<=u+v-1,             (18a)
```

which proves the last estimate in `(12)`.

## 4. Proof of the height obstruction

Under the density half of `(6)`, `(18)` gives

```text
66/91 <= beta_(u,v) <= 11/14.                          (19)
```

The function `x(1-x)` is decreasing on `[1/2,1]`, so

```text
beta_(u,v)(1-beta_(u,v))
 >=(11/14)(3/14)=33/196.                               (20)
```

Apply the component-oscillation lemma `(9)` and the component bound `(12)`:

```text
omega_(u,v)
 >= beta_(u,v)(1-beta_(u,v))/c_(u,v)
 >= 33/[196(u+v-1)].                                   (21)
```

Since `a+b-g=g(u+v-1)`, division by `g` yields

```text
omega_(u,v)/g >= 33/[196(a+b-g)].                      (22)
```

Combining `(22)` with the oscillation half of `(6)` gives

```text
a+b-g >= 33/(196T)
       =33*28710227/(196*1650)
       =585923/200
       =2929 + 123/200.                                (23)
```

The left side is an integer, hence `a+b-g>=2930`.  Since `g>=1`, this implies
`a+b>=2931` and then `max(a,b)>=ceil(2931/2)=1466`.  This proves `(7)`.
**QED.**

If one discards both overlap refinements and uses only
`beta<=6/7`, `c<=u+v`, the same argument gives the older intermediate bound
`a+b>=2131`.  That estimate is superseded here and is not claimed sharp.

The obstruction does not make the THM-4233 gate globally vacuous.  The fixed
primitive pair `(u,v)=(3713,5149)` with `g=1` has the exact values

```text
beta=98322360/133827659,
beta-66/91=16387038/1739759567>0,

omega=277071798/4823550313337,
T-omega=82934716896/2826229070241355051>0.              (23a)
```

Thus both halves of `(6)` hold, consistently with
`a+b-g=3713+5149-1=8861>=2930`.  Both independent repository paths recompute
this certificate; it is the nonvacuity control, whereas `(1,2)` is the sharp
`beta<=11/14` equality control.

## 5. Exact consequence for the live pair graph

THM-4231 proves every fixed-pool outsider pair with `max(a,b)>=770`.  The
theorem above shows that every pair which can satisfy the existing THM-4233
gate already has `max(a,b)>=1466`.  Therefore:

> **Cofinal incompatibility corollary.** Relative to the fixed pool and the
> current deck threshold `T`, the entire THM-4233 pair-observable gate on two
> distinct outsiders is coverage-subsumed by THM-4231.  It cannot remove an
> edge from the finite THM-4231 residual, before or after THM-4238/4242.

In particular, THM-4231/4238/4242 leave exactly `181,126` edges, all with
maximum endpoint at most `769`.  Since the endpoints are distinct and `g>=1`,

```text
a+b-g<=768+769-1=1536.                                 (24)
```

Equation `(22)` gives the uniform strict margin

```text
omega/g >= 33/(196*1536)=11/100352,

11/100352 - 1650/28710227
 =3065953/58798544896>0.                               (25)
```

No live edge can pass `(6)`.  This conclusion uses only the height cap and
the analytic lemma; the exhaustive audits below are hostile checks and locate
the actual closest edge rather than carrying the proof.

## 6. Closest exact hostile in the live residual

Two independent exhaustive paths agree that the closest post-THM-4242 edge
to the oscillation gate is

```text
(a,b)=(466,699)=233(2,3).                              (26)
```

For the primitive pair `(2,3)`, use the grid `N=84`.  The joint-safe
components are

```text
[3,26], [30,39], [45,54], [58,81]                     (27)
```

in grid ticks.  Hence

```text
c=4,       beta=64/84=16/21,
min H=-67/1764,       max H=67/1764,
omega=67/882.                                          (28)
```

After the common dilation,

```text
omega/g=67/205506,

(omega/g)/T=39256841/6920100
             =5.6728719238...>1.                       (29)
```

This is a hostile to any hope that the current gate only narrowly misses the
finite graph.  Gate failure is not Haar danger: `(466,699)` is not asserted
unsafe, and no body-level conclusion is drawn from `(29)`.

## 7. Inheritance pass and connection contract

The closest inherited mechanism is THM-4170's exact `(mass,components)`
discrepancy estimate, used by THM-4231/4238/4240/4242 and spliced with literal
finite heads.  The canonical hostile is THM-4207: two lawful newcomer
marginals need not compose.  The corrected near miss is MISTAKE-520's attempt
to create a new subset layer after full-pool heredity had already closed it;
MISTAKE-524 separately forbids using a reserved forward theorem.  The
least-used sidecar tested here is THM-4233's centered primitive of the
**joint** two-comb observable.  THM-4242's new gcd-fibre tariff was also
inspected; it is a different sufficient, overlap-losing sidecar.

```text
source:       pair (a,b)=g(u,v) and A_(u,v)=G_u intersect G_v
target:       THM-4233 deck gate beta>=66/91, omega/g<=T
map:          joint-safe components -> their H-increments -> omega lower bound
preserved:    common dilation, joint density, component count, exact threshold
destroyed:    component addresses and detailed endpoint phase after summing
sidecar:      danger overlap <=1/14; c_(u,v)<=u+v-1; beta_(u,v)<=11/14
hostile:      (466,699)=233(2,3), the exact closest live edge
decisive test: compare 33/[196(a+b-g)] with T.          (30)
```

The result redirects the current finite frontier: do not rescan the existing
THM-4233 scalar gate on the `181,126` edges.  A useful pair-specific successor
must enlarge the deck tolerance, use component addresses/phase cancellation
beyond the scalar range, or add a different sidecar.  The proved
THM-4242 twenty-three-label chart alone cannot remove another graph edge: it
covers `binom(23,9)=817,190` selected bodies, whereas an edge is removed only
after all `binom(30,9)=14,307,150` pool bodies are closed.  Small-ray/full-body
analytic-literal splices and exact gcd-fibre overlap remain logically live.

## 8. Exact audits and reproduction

The two scripts reconstruct the post-THM-4242 graph directly from the frozen
THM-4231 `K(q)` table and apply all `48+32+36` proved edge closures.  They use
explicit `require` checks, not Python optimization-sensitive assertions.

The primary path sorts the exact common endpoint grid and classifies every
open cell at its midpoint.  The independent path never builds that
arrangement: it intersects the two already ordered safe-tooth lists and
integrates the rational primitive at component endpoints.  Their residual
construction also differs: AST literal extraction plus ordered filtering
versus regex table isolation plus set difference.

Reproduce from the repository root:

```bash
python3 04-computation/lrc14_primitive_gate_danger_overlap_height_grid_audit_thm4245.py \
  --repo .
python3 -O 04-computation/lrc14_primitive_gate_danger_overlap_height_grid_audit_thm4245.py \
  --repo .

python3 04-computation/lrc14_primitive_gate_danger_overlap_height_teeth_audit_thm4245.py \
  --repo .
python3 -O 04-computation/lrc14_primitive_gate_danger_overlap_height_teeth_audit_thm4245.py \
  --repo .
```

Each normal/optimized pair must be byte-identical.

The grid streams must also byte-match
`05-knowledge/results/lrc14_primitive_gate_danger_overlap_height_grid_audit_thm4245.out`,
and the teeth streams must byte-match
`05-knowledge/results/lrc14_primitive_gate_danger_overlap_height_teeth_audit_thm4245.out`.

Both paths verify:

```text
post-THM-4242 edges                    181,126
residual FNV                 bdf59726990a6c92
residual SHA-256             c0e2fe1c...b9be6bf
distinct primitive ratios                 115,429
density-upper-bound violations                  0
density-gate failures                            0
oscillation-gate passes                           0
closest edge                         (466,699)=233(2,3). (31)
```

Over the full `115,429`-element primitive-ratio universe induced by the live
residual, each script explicitly checks `beta<=11/14` and `c<=u+v-1`.  They
also reproduce the THM-4233 `(1,13)` hostile, the `(1,2)` equality control
`beta=11/14`, and the nonvacuous `(3713,5149)` gate certificate in `(23a)`
exactly.  Artifact hashes are frozen in the theorem frontmatter.

## 9. Duplication and scope firewalls

A current-tree exact search for

```text
585923/200 | integer minimum 2930 | danger-overlap component floor
```

and the corresponding `c omega >= beta(1-beta)` route found no earlier
proved theorem, hypothesis, script, or result duplicating this statement.
THM-4228 contains the opposite upper estimate
`omega<=beta(1-beta)`; it does not imply `(9)` or the height threshold.

This theorem does **not**:

1. assert that any failed THM-4233 gate is unsafe;
2. improve the exact `181,126` edge count;
3. lower THM-4231's cutoff `770`;
4. prove that `590` is a minimal fixed-fifty or small-ray boundary;
5. turn `chi_50>=23` into a universal full-pool edge closure;
6. supply physical entry of an arbitrary thirteen-speed row; or
7. prove LRC(14).
