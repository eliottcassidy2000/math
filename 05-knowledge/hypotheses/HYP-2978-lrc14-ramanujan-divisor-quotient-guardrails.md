---
id: HYP-2978
title: LRC14 Ramanujan-divisor quotient guardrails
status: PROOF-INTERFACE / quotient-admissibility proof lane with finite collision audit; not a proof
source: codex-2026-06-24-S161
script: 04-computation/lrc14_ramanujan_divisor_quotient_guardrails_codex_s161.py
result: 05-knowledge/results/lrc14_ramanujan_divisor_quotient_guardrails_codex_s161.out
related:
  - HYP-2979
  - HYP-2977
  - HYP-2976
  - HYP-2975
  - HYP-2974
  - HYP-2973
  - HYP-2972
  - HYP-2971
  - HYP-2970
  - HYP-2969
  - HYP-2968
  - HYP-2963
  - HYP-2956
  - HYP-2946
  - HYP-2938
  - HYP-2887
  - HYP-2886
  - HYP-2264
  - THM-406
  - THM-572
  - OPEN-Q-108
artifacts:
  - 04-computation/lrc14_ramanujan_divisor_quotient_guardrails_codex_s161.py
  - 05-knowledge/results/lrc14_ramanujan_divisor_quotient_guardrails_codex_s161.out
  - 07-reflections/lrc14-ramanujan-divisor-quotient-guardrails-codex-s161.md
---

# HYP-2978: LRC14 Ramanujan-Divisor Quotient Guardrails

This hypothesis is the quotient-admissibility lane requested on 2026-06-24.
The guiding principle is:

```text
A quotient is admissible for an LRC14 proof only if it preserves the predicate
needed by the next implication, or records an explicit certificate explaining
what information was intentionally forgotten.
```

Equivalently, each quotient must declare:

```text
preserved predicate
forgotten labels
compensating transform or side-channel
defect certificate when a forgotten label is load-bearing
```

External seed:

```text
sigma_k(n), Dirichlet convolution, multiplicativity, Jordan/Dedekind unit
capacity, Ramanujan sums c_q(n) as primitive-root power sums, divisor-summatory
hyperbolic-simplex counts, Lambert/Eisenstein coefficient pushes, and Euler
pentagonal recurrences for sigma.
```

Internal seed:

```text
every failed scalar LRC quotient forgot a label: irreducible cores,
unital-design incidence, Faulhaber odd moments, Pollock defect degree,
unit-distance norm layers, tiling/solid units, endpoint owners,
C27/K33 state debt, or harmonic dual data.
```

The immediate LRC14 hook is HYP-2974's divisor-curried Fourier coefficient:

```text
hat F_S(k) = sum_{v in S, v|k} sin(pi*(k/v)/7)/(pi*(k/v)).
```

Mode `k` sees both the divisor fiber `v|k` and the quotient `k/v`.
Ramanujan sums `c_q(n)`, the sums of nth powers of primitive q-th roots, are
therefore candidate exact-period unit characters: they retain primitive
residue-period data after averaging over units instead of collapsing to a bare
divisor count.  HYP-2979 is the companion retained-packet route: exact-period
Ramanujan projectors for q-ladders, endpoint sums, primitive phase packets, and
shifted danger-count autocorrelations.

## Computation

Script:

```text
04-computation/lrc14_ramanujan_divisor_quotient_guardrails_codex_s161.py
```

Stored output:

```text
05-knowledge/results/lrc14_ramanujan_divisor_quotient_guardrails_codex_s161.out
```

The script verifies through `n<=80` the identities:

```text
sigma_0=tau
sum_{d|n} phi(d)=n
phi=mu*id
psi=id*|mu|
sum_{d|n} J_2(d)=n^2
c_q(n)=sum_{d|gcd(q,n)} d*mu(q/d)
```

It then audits named LRC14 rows:

```text
AP, GW 12->24, residue liar 12->26, near 12->36, petals 10->20
and 13->26, P10+GW, P10+K33, covering 12->84/168, covering 6->98.
```

Main numerical readout:

```text
qdiv_only route-mixing collisions                 1
open_state_only route-mixing collisions           1
mod14_residue_multiset route-mixing collisions    1
ramanujan_14_profile route-mixing collisions      1
unit_counts_14_27_41 route-mixing collisions      2
divisor_lcm_scalars route-mixing collisions       1
guarded_packet_signature route-mixing collisions  0
```

The collisions are the point.  Coarse divisor data, `qdiv`, open/zero-open
status, mod-14 residues, and a single Ramanujan shell are all useful features,
but they are not admissible proof quotients by themselves.  They identify rows
with different proof routes: AP/GW boundary atoms, q-witness rows, K33/state
lifts, unit-petal rows, and covering-moment rows.

For example, the AP, GW, residue-liar, near-K33, petal, and P10-splice rows
mostly share the same lcm scalar data:

```text
lcm = 2^3*3^2*5*7*11*13
tau(lcm)=192
phi(lcm)/lcm=192/1001
psi(lcm)/lcm=2304/715
```

Many also share the same `c_14` primitive-shell profile:

```text
((-6,1), (-1,6), (1,6)).
```

Thus lcm-divisor data and the `q=14` primitive shell cannot replace exact
q/Farey/Haar/C27/K33 labels.  Conversely, adding q-threshold, Haar open state,
Ramanujan profiles at `14,27,41`, unit witness counts, and the current route
label separates the named audit bank with zero route-mixing collisions.  This
does not prove LRC14; it states the smallest currently honest packet type.

## Tournament Analysis

Vertices are quotient channels rather than runners:

```text
raw_divisor_counts
squarefree_psi_support
totient_jordan_unit_capacity
gcd_strata
ramanujan_primitive_shell
exact_period_packet
labelled_lrc_packet_sheaf
```

Pairwise observable:

```text
A -> B iff A retains more of the declared proof payload:
divisor capacity, squarefree support, unit capacity, gcd strata,
primitive phase, exact q/Farey packet, and endpoint/state/dual labels.
```

Stored fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3_cycles=0
Hamiltonian path:
labelled_lrc_packet_sheaf >
exact_period_packet >
ramanujan_primitive_shell >
gcd_strata >
totient_jordan_unit_capacity >
squarefree_psi_support >
raw_divisor_counts
```

The readout is not that Ramanujan shells are complete.  It is that they are the
first reasonable arithmetic side-channel above scalar divisor counts, while
endpoint owners, q/Farey data, and state-lift labels remain load-bearing.

Assumption challenge: considered vertices included runners, residues, divisors,
unit twists, gcd strata, primitive Fourier modes, endpoint owners, C27/K33
packets, and proof obligations.  This audit chooses quotient channels.  The
preserved predicate is whether a channel retains enough payload to support the
next LRC implication without silently mixing AP/GW, q-witness, K33, petal, and
covering routes.  It intentionally destroys raw runner identity when that
identity is not part of the next proof predicate.

## Web Pass And Source Trail

The divisor-function page makes the useful warning explicit.  `sigma_a(n)` is
both a divisor-fiber count/sum and a convolutional coefficient: its Dirichlet
series is `zeta(s) zeta(s-a)`, its Lambert series re-sums the same divisor
fibers into harmonic/modular coefficients, and for `k>0` it has a Ramanujan
expansion

```text
sigma_k(n) = zeta(k+1) n^k sum_m c_m(n)/m^(k+1).
```

So `sigma`, `tau`, and `psi` are not proof predicates by themselves.  They are
pushforwards of divisor fibers, useful only if the kernel of the pushforward is
declared.

One-hop reading sharpened the algebra:

- Dirichlet convolution says arithmetic functions form a local Dirichlet ring;
  multiplicative functions form a subgroup under convolution, and identities
  such as `d=1*1`, `sigma_k=Id_k*1`, `phi*1=Id`, and `mu*1=epsilon` are packet
  laws rather than loose analogies.
- Jordan's totient `J_k` counts primitive `k`-tuple capacity and satisfies
  `sum_{d|n} J_k(d)=n^k`, so it is the higher-rank exact-period version of the
  `phi` packet used in HYP-2899.
- The divisor summatory function counts lattice points under `jk<=x`; the
  divisor problem is exactly a boundary-defect problem after quotienting the
  product lattice by the hyperbola.
- Ramanujan sums retain the primitive unit shell and have orthogonality over
  lcm periods.  Murty's survey records the Carmichael-style shifted
  orthogonality `M(c_r(n)c_s(n+h))=c_r(h)` when `r=s` and `0` otherwise, which
  directly suggests shifted packet autocorrelations for LRC endpoint sums.
- The supercharacter reading of Ramanujan sums says the same thing in
  representation language: quotient by unit-group orbits is legitimate only
  because the orbit character is retained.

External pages read in this pass:

```text
https://en.wikipedia.org/wiki/Divisor_function
https://en.wikipedia.org/wiki/Multiplicative_function
https://en.wikipedia.org/wiki/Dirichlet_convolution
https://en.wikipedia.org/wiki/Ramanujan_sum
https://en.wikipedia.org/wiki/Jordan%27s_totient_function
https://en.wikipedia.org/wiki/Divisor_summatory_function
https://en.wikipedia.org/wiki/Lambert_series
https://en.wikipedia.org/wiki/Eisenstein_series
https://arxiv.org/abs/1201.1060
https://mast.queensu.ca/~murty/HRJ-2013.pdf
```

## Guardrail Theorem Target

Let `P` be the predicate needed by the next LRC implication, and let
`pi:X->Q` be a proposed quotient.  The quotient is admissible only if every
fiber of `pi` satisfies at least one of:

```text
P is constant on the fiber;
the forgotten coordinate is reconstructible from retained packet data;
the forgotten coordinate is killed by a proved orthogonality/dual certificate;
the fiber is assigned to an explicit residual bucket with endpoint/state labels.
```

Equivalently, every coordinate forgotten by an LRC14 proof quotient must be:

```text
invariant under the proof step;
reconstructible from retained q/Farey, unit, gcd, and Ramanujan data;
annihilated by a proved orthogonality or dual certificate; or
placed in an explicit residual bucket with endpoint/state-lift labels.
```

Expected falsifier shape:

```text
two rows with the same scalar divisor signature but different LRC route.
```

The S161 audit already exhibits such pairs for multiple scalar channels, so
those channels are demoted to features.  Any final quotient must include
Ramanujan packet labels or exact endpoint-owner/Farey data.

This is the common shape behind the older lessons:

- irreducibility/coefficient quotients need Newton slope, fixed-divisor, or
  residue-depth side channels before sign/count shadows are trusted;
- unital prompts split geometric unitals, algebraic identity preservation, and
  unit groups instead of using the same word for incompatible kernels;
- Faulhaber moment shadows preserve positivity but not compatibility/integrality
  packets;
- Pollock dyadic residues fail singly, while carry-pair defect lifts recover the
  load-bearing predicate;
- unit-distance scalar edge counts need spine/bulk/direction/ear labels;
- tiling and regular-solid analogies work only when curvature, norm recursion,
  and local vertex-figure defects are retained;
- HYP-2899's `Div(D) x B_r` product lattice is the exact denominator/far-packet
  version of the same rule.

Thus the fundamental object is not a connection between two scalars.  It is a
connection whose kernel has been named and discharged.

## Tournament-Adjacent Proof Angle

Use Tournament Analysis on quotient channels or proof obligations, not on
runners.  The proof-facing pairwise observable is:

```text
Does channel A separate every route collision that channel B mixes, while
preserving the LRC predicate: strict witness, AP/GW boundary equality, or named
state-lift debt?
```

Ties are broken by the Hamiltonian path:

```text
exact-period packet
> endpoint-owner label
> dual certificate
> residual bucket
```

The S161 named-row audit gives the first transitive channel tournament.  The
proof-facing version should be stronger:

```text
Every post-THM-571 Moon-core packet maps into the labelled LRC packet sheaf.
If a coarser quotient is used, its kernel must be contained in the AP/GW
equality atoms plus packets killed by Ramanujan orthogonality, Toeplitz PSD,
moment duals, taut-current potentials, or HYP-2908/THM-572 state lifts.
```

A hypothetical LRC14 counterexample is then not just "not seen" by current
features.  It must lie in the intersection of the kernels of q/Farey witness
ladders, Ramanujan primitive shells, Haar-open fronts, Toeplitz PSD tests,
danger-count and multiplicity moment duals, endpoint taut currents, C27/K33
owner labels, and covering boundary-moment packets.  HYP-2978's job is to make
that intersection exact and then empty it, or name the first genuinely new
residual.

Concrete route toward LRC14:

```text
Use divisor/totient/Jordan/psi data as capacity ledgers.
Use Ramanujan c_q profiles as primitive-shell Fourier guards.
Then hand off to HYP-2974/HYP-2977 harmonic duals or HYP-2965/HYP-2969
endpoint packets before claiming that a qdiv>14 zero-open packet cannot exist.
```

## Inquiry Bucket

- Turn regular A-functions into "declared divisor-subset" quotient metadata.
- Test divisor-summatory hyperbolic-simplex defects against Pollock and tiling
  carrier packets.
- Use Euler pentagonal sigma recurrences as sparse boundary moments in the
  HYP-2974/HYP-2977 harmonic stack.
- Treat totient as a finite Fourier transform of gcd and compare with
  Ramanujan exact-period packets.
- Use Busche-Ramanujan identities as non-coprime repair terms for any product
  quotient that tries to forget shared prime support.
- Extend the named-row audit to HYP-2963 packet families and then to the
  finite-check endpoint, looking for the first quotient whose route-mixing
  count remains zero without including the route label itself.
- Test HYP-2979's shifted Ramanujan autocorrelations against HYP-2973 count
  moments and HYP-2974 Toeplitz/Fejer failures.
