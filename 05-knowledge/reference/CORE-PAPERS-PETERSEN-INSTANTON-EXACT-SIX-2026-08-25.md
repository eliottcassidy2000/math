# Petersen minors, instanton passages, and exact period six: source ledger

> **Audit date: 2026-08-25.** This file separates three primary-source
> imports from elementary finite checks. Pintér and Shapiro are recent
> external preprints; neither is promoted to internally `PROVED`. Stoll's
> published rational-six-cycle conclusion is conditional. The elementary
> order and squaring-map calculation is routed to
> [THM-4139](../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md),
> which was subsequently promoted to `PROVED / FINITE-EXACT / VERIFIED` with
> an independent audit. This ledger is source metadata and scope control, not
> a second theorem record.

## Status at a glance

| item | exact status here | imported content |
|---|---|---|
| Pintér, arXiv:2607.22267v2 | **CITED / EXTERNAL PREPRINT v2**; appendix models **VERIFIED-EXACT** independently | a `P/e`-minor-free bridgeless multigraph has a nowhere-zero `4`-flow |
| Shapiro, arXiv:2608.23342v1 | **CITED / EXTERNAL UNREFEREED PREPRINT v1**; statement and proof interface audited | fixed-mean Poisson law for instanton passage counts in a double well |
| Stoll, arXiv:0803.2836v2 | **PUBLISHED / CONDITIONAL** for rational-six-cycle nonexistence; Theorem 6 is unconditional as an implication from `rank J(Q)=3` | arithmetic-dynamical curve and Chabauty reduction for quadratic six-cycles |
| `2^6-1` and `z -> z^2` calculation | **ELEMENTARY / VERIFIED-EXACT** here and in promoted THM-4139 | no primitive prime at exponent six, but exact order six at modulus `9`, and exact-period-six torsion strata |

## 1. Petersen contraction and nowhere-zero four-flows

### Primary source and exact import

- József Pintér,
  [*Nowhere-zero 4-flows in graphs excluding a proper minor of the Petersen
  graph*](https://arxiv.org/abs/2607.22267v2), arXiv:2607.22267v2;
  [official HTML](https://arxiv.org/html/2607.22267v2).
- Put `P` for the Petersen graph and `Q=P/e` for the contraction of any edge
  of `P`. Theorem 1.1 states that every finite bridgeless undirected
  `Q`-minor-free multigraph admits a nowhere-zero `Z_2^2`-flow, hence a
  nowhere-zero `4`-flow.
- Together with the Thomas--Thomson `P-e` theorem, Corollary 1.2 says that a
  finite bridgeless graph with no nowhere-zero `4`-flow contains **every
  proper minor** of `P`, in particular both `P-e` and `P/e`.
- For a loopless cubic graph, a nowhere-zero `Z_2^2`-flow is equivalent to a
  proper three-edge-coloring. Thus every snark in the usual finite bridgeless
  cubic sense contains both proper Petersen minors. This consequence is
  exact; it is not a classification of snarks.

The proof route is a minimal flow obstruction, reduced using imported
Thomas--Thomson lemmas to a three-connected, minimum-degree-three,
girth-at-least-five, almost-four-connected graph. The planar branch uses the
Four Color Theorem. The nonplanar branch uses the Triplex/Petersen/
Dodecahedron/Basket atlas and Norin--Thomas jump/cross extensions. Those
imported structural theorems were not re-proved in this audit.

### Independent finite check

An independent exact Python replay reconstructed `Q`, all seven appendix
hosts, the stated deletions, branch sets, and the fourteen `Q`-adjacency
witness edges in each model. It
asserted pairwise disjoint nonempty branch sets, connectivity of every branch
set, validity of retained host edges, and realization of every adjacency of
`Q`. All checks passed:

```text
host       original edges  retained edges  branch-set sizes
P                    15              15     2,1,1,1,1,1,1,1,1
Triplex              18              17     appendix model
Basket               18              17     appendix model
D(0,3)-jump          31              23     appendix model
D(0,4)-jump          31              22     appendix model
D(0,5)-jump          31              21     appendix model
D-cross              32              23     appendix model
models verified: 7
```

The SHA-256 of the canonical **summary payload** consisting of the sorted
edge set of `Q` and each host's retained-edge/branch-size record is

```text
e48fc97585a3ba2138935439914bd94957ab9a022f8e1da012e5d2c528934f5c
```

This is a semantic-output checksum, not a source-file hash and not a hash of
the full branch-set input. The assertions did inspect the full reconstructed
models. The replay does **not** certify the Thomas--Thomson lemmas, the
Norin--Thomas theorem, or the global minimal-counterexample proof.

### Controls and transfer boundary

- **Positive obstruction control:** `P` has no nowhere-zero `Z_2^2`-flow and
  contains `Q`.
- **Hostile converse control:** `K_10` contains `Q` as a minor but has a
  nowhere-zero `Z_2^2`-flow. Indeed its standard one-factorization has nine
  perfect matchings; assign three matchings to each of the three nonzero
  elements of `Z_2^2`. Each vertex then sees three edges of each value, whose
  sum is zero. Thus containing `Q` is necessary for a flow obstruction, not
  sufficient.
- A minor certificate forgets the embedding of its branch sets. It also
  forgets the ordered boundary support and extension multiplicities used in
  the repo's snark four-pole calculus. Those data require a separate
  boundary-state sidecar.
- Do not call snarks “primes” merely because every one contains these
  Petersen minors. Minor containment is neither unique factorization nor a
  decomposition operation; the correction recorded in
  [MISTAKE-507](../../01-canon/MISTAKES.md) applies.

## 2. Double-well instantons and a Poisson passage law

### Primary source and exact import

- Jacob Shapiro,
  [*Instantons in a Double-Well are Poisson Distributed*](https://arxiv.org/abs/2608.23342),
  arXiv:2608.23342v1, submitted 2026-08-24;
  [official HTML](https://arxiv.org/html/2608.23342).

The paper studies

```text
H_lambda = -Delta + lambda^2 V
```

on `R^n`. Its hypotheses include an even `C^2`, nonnegative, confining
potential with exactly two nondegenerate minima `+d,-d`, a two-level gap
`E_2-E_1 >= c lambda`, and the stated one-well gap and localization
hypotheses. Passages are alternating hitting times between shrinking balls of
radius `lambda^(-1/4)` around the wells, under a localized
Feynman--Kac-weighted Brownian bridge.

For the paper's negative hopping coefficient

```text
rho_lambda = -lim_(beta -> infinity) Z_(lambda,1)/(beta Z_(lambda,0)),
```

Theorem 2.1 states that for fixed `N>0`, independent of `lambda`, and each
fixed `k`, the passage count at time `N/|rho_lambda|` tends to

```text
exp(-N) N^k/k!.
```

It also obtains

```text
-lambda^(-1) log |rho_lambda| -> S(d,-d),
```

and Corollary 2.2 gives `E_1-E_0=2|rho_lambda|(1+o(1))`.

The mechanism is more portable than the metaphor: the normalized passage
kernel is split into a time-independent rank-one ground-state term plus a
defect. Integrated-mass, first-moment, and pointwise-in-time defect estimates
vanish. Lemma 4.2 gives the fixed-`k` sector ratio `N^k/k!`; Lemma 4.3 supplies
a summable uniform majorant of the form `C_N 2^(-k)`, permitting normalization
and exchange of limit and sum.

The abstract probability gate is elementary and theorem-ready: if
nonnegative weights have fixed-`k` ratios tending to `N^k/k!` **and** admit a
summable uniform majorant, their normalized laws converge to `Poisson(N)`;
on the countable state space, discrete Scheffe then gives total-variation
convergence. The domination hypothesis is essential. Fixed-`k` convergence
alone fails after adding mass `exp(L^2)` at an index `k=L` escaping to
infinity.

### Scope firewall

- The theorem fixes `N`; this audit found no uniform-in-`N` conclusion.
- It is a count-law theorem. Do not silently upgrade it to convergence of the
  full passage point process.
- A deterministic exact-period cycle has no probability measure, rare-event
  scaling, or escaping-sector estimate. The number `63` therefore has no
  Poisson interpretation here.
- Two special target fibers in a planar Jacobian candidate are not stochastic
  wells. A lawful transfer would first have to construct a positive path
  measure or operator, a rank-one leading kernel, a spectral gap, and the
  three defect estimates. Without those typed objects, “rank one plus defect”
  is only a proof-interface analogy.

## 3. Stoll's quadratic six-cycle reduction

### Primary source and exact status

- Michael Stoll,
  [*Rational 6-cycles under iteration of quadratic polynomials*](https://arxiv.org/abs/0803.2836),
  arXiv:0803.2836v2;
  [printed-version PDF](https://arxiv.org/pdf/0803.2836v2). Published in
  *LMS Journal of Computation and Mathematics* **11** (2008), 367--380.

The headline nonexistence result is **conditional**, not an unconditional
resolution. Theorem 7 assumes that the `L`-series of
`J=Jac(X_0^dyn(6))` extends to an entire function and satisfies the standard
functional equation, and assumes weak BSD for `J`; under those hypotheses
there is no rational exact length-six cycle for `x -> x^2+c`.

The unconditional intermediate Theorem 6 says: **if** `rank J(Q)=3`, then
`X_0^dyn(6)` has only the ten rational points listed in the paper, and hence
`X_1^dyn(6)` has only cusps and no rational exact six-cycle. The rank-three
divisor subgroup and Chabauty computation prove the rational-point statement
for points mapping into its saturation; the analytic/BSD input is used to
bound the global rank by three.

The dynamical curve `X_1^dyn(6)` has the order-six automorphism induced by
iteration; its quotient `X_0^dyn(6)` has genus four. A rational point of the
quotient can represent a Galois-stable six-cycle without any individual
point of that cycle being rational. This distinction is essential for the
paper's torsion examples.

### Exact torsion bridge, and its limits

Stoll records quotient points whose cycles contain:

```text
c=0:   zeta_9;
c=-2:  zeta_13 + zeta_13^(-1), or zeta_21 + zeta_21^(-1).
```

The first is the fixed torus map `z -> z^2`, where `ord_9(2)=6`. For the
other two, the typed bridge is the Chebyshev semiconjugacy

```text
pi(z)=z+z^(-1),       pi(z^2)=pi(z)^2-2.
```

Consequently the exact period of `pi(zeta_m)` under `x -> x^2-2` is the least
positive `k` with `2^k = +1` or `-1 (mod m)`. For `m=21`, the proper divisors
`1,2,3` of six give residues `2,4,8`, none congruent to `+/-1`, while
`2^6=1 (mod 21)`. For `m=13`, those proper divisors again fail and
`2^6=-1 (mod 13)`. Both traces therefore have exact period six, even though
`zeta_13` itself has period twelve under squaring.

These are algebraic, Galois-stable cycles, not rational individual
six-cycles. Stoll's theorem neither proves Zsigmondy's theorem nor identifies
all complex exact-period-six points of `z -> z^2`.

## 4. The `2^6-1` exception and the squaring map

For a prime `p`, being a primitive prime divisor of `2^n-1` is equivalent to
`ord_p(2)=n`. At `n=6`,

```text
2^6-1 = 63 = 3^2 * 7,
ord_3(2)=2,             ord_7(2)=3.
```

Thus exponent six has no primitive **prime** divisor, the familiar
Zsigmondy exception. Nevertheless

```text
ord_9(2)=6, because 2^3=-1 (mod 9).
```

It is accurate to call `9` a primitive prime-power modulus at exponent six;
it is inaccurate to call `9` a primitive prime divisor. More generally LTE
gives

```text
v_3(2^(2*3^j)-1)=1+j,
ord_(3^a)(2)=2*3^(a-1)  for a>=1.
```

For `f(z)=z^2` on `G_m`, a primitive `m`-th root has exact period
`ord_m(2)`. Mobius inversion gives the exact-period-six polynomial

```text
D_6(z)
 = ((z^63-1)(z-1))/((z^7-1)(z^3-1))
 = Phi_9(z) Phi_21(z) Phi_63(z).
```

Hence the exact-period-six conductors are `9,21,63`, and there are

```text
phi(9)+phi(21)+phi(63)=6+12+36=54
```

exact-period-six points. The `63` solutions of `f^6(z)=z` decompose by exact
period as `1+2+6+54=63`.

The hostile control is `ord_21(2)=lcm(ord_3(2),ord_7(2))=6`. Exact period six
therefore has two incomparable minimal mechanisms: prime-power depth at `9`
and CRT order assembly at `21`; `63` is their join, not the unique primitive
carrier. A Boolean exact-period observer cannot distinguish conductors `21`
and `63`. The lost sidecar is the conductor, equivalently its prime-adic
valuation vector.

An independent integer calculation checked the factorization, the orders for
all divisors of `63`, the three cyclotomic factors, and the exact-point counts.
Its canonical semantic-output SHA-256 is

```text
32c19eb63315d92d6f2dcb1859a3aa8b445d2aae8b517772831b00c3f117e4b3
```

The checked order table was

```text
m          1   3   7   9   21  63
ord_m(2)   1   2   3   6    6   6
```

The convention in the `m=1` display is the fixed identity stratum, not the
usual multiplicative-order definition modulo one. As a positive control,
`31` is a primitive prime divisor of `2^5-1` and `ord_31(2)=5`.

The rigorous typed map is

```text
source:    modulus m with gcd(m,2)=1
target:    primitive torsion stratum mu_m^prim under f(z)=z^2
preserved: m | 2^n-1  iff f^n fixes that stratum;
           ord_m(2)=n iff its points have exact period n
lost by radical/prime-support quotient: conductor depth and CRT assembly
sidecar:   (v_p(m))_p, or the full conductor m
```

## 5. Cross-frontier firewalls

1. **Planar Jacobian conjecture.** In the current planar-JC route, the
   occurrence of `21` as a function-field degree and `63=3*21` as a pole
   divisor degree is differently typed from torsion conductors and fixed-point
   counts. There is no canonical squaring action or preserved Keller
   predicate connecting them. The target-fiber multiplicities
   `1,2,2,3,7,3,3` also rule out the naive connected cyclic-Galois
   degree-`21` reading; the relevant seam is separately excluded in
   [THM-4130](../../01-canon/theorems/THM-4130-theta-only-extremal-seam-critical-monodromy-obstruction.md).
   Equality of the numerals is not a mathematical map.
2. **The word “Jacobian.”** Stoll's `Jac(X_0^dyn(6))` is the abelian variety
   attached to a genus-four curve. It is not the Jacobian determinant in the
   planar polynomial-map conjecture.
3. **Snarks.** Pintér supplies a necessary proper-minor carrier for every
   cubic flow obstruction. It supplies neither boundary-coloring
   multiplicities nor a factorization theory. Those require the ordered
   four-pole sidecar.
4. **Poisson language.** Shapiro's law is powered by a normalized stochastic
   path measure, rare-event time scaling, a rank-one kernel, and a uniform
   summable tail bound. Exact-period-six torsion and deterministic snark
   colorings have none of these by default. No Poisson claim transfers without
   constructing all four objects.
5. **Arithmetic dynamics.** Stoll concerns `Q`-rational cycles in the family
   `x^2+c`; the elementary factorization above counts complex torsion points
   for the single map `z^2`. The Chebyshev quotient at `c=-2` additionally
   identifies `z` with `z^(-1)`, changing the period test from `2^k=1` to
   `2^k=+/-1` modulo the conductor.

## Reproduction boundary

The two SHA-256 values above certify canonical finite **outputs**, not the
external papers or their source files. The Petersen replay is independent of
the appendix's prose but covers only the seven explicit minor models. The
arithmetic replay is exhaustive for divisors of `63` and the displayed
cyclotomic factorization. No computation here verifies Shapiro's analytic
estimates, Stoll's Chabauty/rank calculations, or Pintér's imported structural
graph theory. Those claims retain exactly the external statuses shown in the
opening table.
