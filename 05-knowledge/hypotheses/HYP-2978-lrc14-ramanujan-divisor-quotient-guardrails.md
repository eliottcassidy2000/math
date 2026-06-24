---
id: HYP-2978
title: LRC14 Ramanujan-divisor quotient guardrails
status: INQUIRY / quotient-admissibility proof lane with finite collision audit; not a proof
source: codex-2026-06-24-S161
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
---

# HYP-2978: LRC14 Ramanujan-Divisor Quotient Guardrails

This hypothesis reserves the divisor/Ramanujan quotient-admissibility lane
requested on 2026-06-24.  The guiding principle is:

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

The intended external seed is the divisor-function neighborhood:
`sigma_k(n)`, Dirichlet convolution, multiplicativity, Ramanujan sums
`c_q(n)` as primitive-root power sums, and the bridge from divisibility data to
Fourier/cyclotomic packets.  The intended internal seed is the repeated repo
lesson that scalar quotients are useful only after labelled fibers are retained:
irreducible cores, unital designs, Faulhaber moment positivity, Pollock degree
jumps, unit-distance norm layers, tiling/solid analogies, and the current LRC14
dual stack.

The immediate LRC14 hook is HYP-2974's divisor-curried Fourier coefficient:

```text
hat F_S(k) = sum_{v in S, v|k} sin(pi*(k/v)/7)/(pi*(k/v)).
```

Mode `k` sees both the divisor fiber `v|k` and the quotient `k/v`.
Ramanujan sums `c_q(n)`, the sums of nth powers of primitive q-th roots, are
therefore candidate exact-period unit characters: they retain primitive
residue-period data after averaging over units instead of collapsing to a bare
divisor count.

HYP-2979 is the companion retained-packet route: exact-period Ramanujan
projectors for q-ladders, endpoint sums, and primitive unit phase packets.

Finite audit added by S161:

```text
qdiv_only route-mixing collisions                 1
open_state_only route-mixing collisions           1
mod14_residue_multiset route-mixing collisions    1
ramanujan_14_profile route-mixing collisions      1
unit_counts_14_27_41 route-mixing collisions      2
divisor_lcm_scalars route-mixing collisions       1
guarded_packet_signature route-mixing collisions  0
```

The named-row table shows the main warning sharply.  AP, GW, the residue liar
`12->26`, near/K33 `12->36`, petals, and `P10+K33` all share the same coarse
lcm divisor scalars; several also share the same `c_14` profile.  Those
channels mix AP/GW boundary, q-witness, K33, petal, and covering proof routes.
The over-labelled packet signature avoids route mixing precisely because it
keeps q/Farey, Haar-open status, endpoint/state labels, and the relevant dual
certificate exits attached.

Tournament Analysis over quotient channels is transitive:

```text
labelled_lrc_packet_sheaf
  > exact_period_packet
  > ramanujan_primitive_shell
  > gcd_strata
  > totient_jordan_unit_capacity
  > squarefree_psi_support
  > raw_divisor_counts
```

with `score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}` and no directed 3-cycles.
The readout is not that Ramanujan shells are complete.  It is that they are the
first reasonable arithmetic side-channel above scalar divisor counts, while
endpoint owners and state-lift labels remain load-bearing.

## Web pass and source trail

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

## Guardrail theorem target

The abstract lesson across the repo can now be stated as a quotient-kernel
rule.  Let `P` be the predicate needed by the next LRC implication, and let
`pi:X->Q` be a proposed quotient.  The quotient is admissible only if every
fiber of `pi` satisfies at least one of:

```text
P is constant on the fiber;
the forgotten coordinate is reconstructible from retained packet data;
the forgotten coordinate is killed by a proved orthogonality/dual certificate;
the fiber is assigned to an explicit residual bucket with endpoint/state labels.
```

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

## Tournament-adjacent proof angle

Use Tournament Analysis on quotient channels or proof obligations, not on
runners.  The pairwise observable is:

```text
Does channel A separate every route collision that channel B mixes, while
preserving the LRC predicate: strict witness, AP/GW boundary equality, or named
state-lift debt?
```

Ties are broken by the Hamiltonian path

```text
exact-period packet
> endpoint-owner label
> dual certificate
> residual bucket
```

The S161 named-row audit already gives the first transitive channel tournament.
The proof-facing version should be stronger:

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

Remaining work:

1. Extend the collision audit from named rows to packet families in the
   HYP-2963 bank.
2. Turn the quotient-kernel rule into a theorem-facing lemma over the current
   labelled packet sheaf.
3. Test HYP-2979's shifted Ramanujan autocorrelations against HYP-2973 count
   moments and HYP-2974 Toeplitz/Fejer failures.
4. Stress-test this admissibility criterion against new packet-family mixed
   fibers:

```text
Any quotient used to rule out an LRC14 counterexample must retain enough
divisor/cyclotomic phase data to distinguish AP/GW boundary atoms, positive
Toeplitz/Ramanujan exits, and K33/state-lift debts.
```

Expected falsifier shape: two rows with the same scalar divisor signature but
different LRC route.  The S161 audit already exhibits such pairs for multiple
scalar channels, so those channels are demoted to features.  Any final quotient
must include Ramanujan packet labels or exact endpoint-owner/Farey data.

Broader theorem target:

```text
Every post-THM-571 Moon-core packet either
  has a quotient whose preserved data forces a known dual certificate,
  or exhibits a Ramanujan/divisor defect showing exactly which label
  the quotient illegally forgot.
```
