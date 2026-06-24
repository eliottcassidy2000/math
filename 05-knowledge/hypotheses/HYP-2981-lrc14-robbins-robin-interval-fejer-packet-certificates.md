---
id: HYP-2981
title: LRC14 Fejer interval certificates and Robbins/Robin quotient-bridge guardrails
status: PROOF-INTERFACE / packet-anchored interval-certificate scaffold and precision blueprint; not a proof
source: codex-2026-06-24-S162
related:
  - HYP-2979
  - HYP-2978
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
  - HYP-2966
  - HYP-2965
  - HYP-2964
  - HYP-2963
  - HYP-2956
  - HYP-2908
  - THM-572
  - OPEN-Q-108
artifacts:
  - 04-computation/lrc14_packet_fejer_interval_scaffold_codex_s162.py
  - 05-knowledge/results/lrc14_packet_fejer_interval_scaffold_codex_s162.out
  - 07-reflections/lrc14-robin-robbins-fejer-interval-scaffold-codex-s162.md
  - 04-computation/lrc14_fejer_interval_packet_certificates_codex_s162.py
  - 05-knowledge/results/lrc14_fejer_interval_packet_certificates_codex_s162.out
  - 07-reflections/lrc14-robbins-robin-interval-fejer-packet-certificates-codex-s162.md
---

# HYP-2981: LRC14 Fejer Interval Certificates and Robbins/Robin Guardrails

This hypothesis sharpens the immediate proof target after HYP-2974/S157:

```text
prove no non-AP/GW labelled packet fiber can remain Toeplitz-PSD-blind once
the floating Fejer evaluations are replaced by rigorous interval certificates.
```

It integrates the HYP-2981/T1065 reservation and the later packet scaffold.
The spelling collision is useful but should not be blurred: the graph theorem
is Robbins, while the divisor-function/RH criterion is Robin.  A third
side-reading, associated in the repo with Neville Robbins style
divisor/partition fibers, reinforces the same point: divisor functions and
partition products become safe only after the fixed fiber has been declared.
The shared rule is:

```text
connection = map + declared kernel + certificate for every forgotten bridge.
```

The motivating finite evidence is strong.  HYP-2974's full HYP-2963-bank scan
checks `21913` rows.  AP and Goddyn-Wong are the only zero-safe rows, and all
`21911` positive rows have an explicit Fejer-vector PSD violation by degree
`<=280`.  The hardest named packets are exactly the ones the current proof
geography predicts:

```text
P10+GW                         degree 280
K33 / 12->36                   degree 159
covering 12->168               degree  63 but small margin
two drop(12,13)->add(14,29)    degree  41 and smallest full-bank margin
```

The remaining gap is not conceptual negativity.  It is formalization:
the signs are currently floating or prototype interval trigonometric
evaluations.  They need to become interval-arithmetic certificates attached to
labelled packet families `P(S)`.

## Spectrum Lesson

The user's spectrum point changes how to prioritize proof effort.  The family

```text
M_k = {1,2,...,k-1,2k}
```

has lonely value `G=2/(2k+1)`, so the second LRC-spectrum point satisfies

```text
1/(k+1) < sigma_2(k) <= 2/(2k+1).
```

The gap from the tight AP value is only

```text
2/(2k+1) - 1/(k+1) = 1/((k+1)(2k+1)).
```

So any lift-depth lemma that must resolve the second spectral point is forced
into depth on the order of `(k+1)(2k+1)`.  For `k=13` runners below the LRC14
threshold this is already a large finite debt.  The Fejer/Toeplitz packet
route is attractive because it does not wait for all runners to return to a
total reset clock.  It certifies a labelled packet fiber by a local harmonic
contradiction to `F_S=C_S-1>=0`.

## Certificate Shape

Both S162 computations use the Fejer quadratic form

```text
Q_d(x)=6/7 + 2 sum_{k=1}^d (1-k/(d+1)) cos(2*pi*k*x)
             sum_{v|k, v in S} sin(pi*(k/v)/7)/(pi*(k/v)).
```

Each term is a divisor-curried atom labelled by `(k,v,m=k/v)`.  For a rigorous
certificate, it is enough to outward-enclose the rational center `x`, the
root-of-unity cosine, the algebraic sine `sin(pi*m/7)`, and the reciprocal of
`pi`, then show the resulting interval for `Q_d(x)` stays negative.

The companion scaffold
`04-computation/lrc14_packet_fejer_interval_scaffold_codex_s162.py` gives the
candidate certificate object.  It imports exact safe components and packet
keys, uses rational interval arithmetic with a hard-coded rational enclosure
of `pi`, and emits interval upper endpoints `<0` for five selected rows:

```text
near/K33 12->36
P10+GW
covering 12->168
two drop(12,13)->add(14,29)
single swap 6->63
```

This is not yet the production certificate emitter.  The current rational
`pi` enclosure and Taylor interval engine must be replaced by a formally
sourced backend, and the final pass should group the HYP-2963 bank by packet
family rather than treating `21911` rows individually.

## Precision Blueprint

The budget script
`04-computation/lrc14_fejer_interval_packet_certificates_codex_s162.py`
recomputes selected hard S157 certificates, expands each certificate into its
divisor-curried atom bank, and gives a deliberately conservative precision
budget.  If the atom L1 is `L` and the negative margin is `eta`, relative atom
intervals of size roughly `2^-p` with

```text
L * 2^-p <= eta/8
```

leave slack for outward rounding.  On the selected hard rows, this yields a
surprisingly finite-looking burden:

```text
P10+GW                         degree 280, atoms 862, bits 17
K33 / 12->36                   degree 159, atoms 490, bits 19
covering 12->168               degree  63, atoms 190, bits 22
two drop(12,13)->add(14,29)    degree  41, atoms 122, bits 27
```

The smallest full-bank floating margin is the last row:

```text
Q_41(347/4312) = -3.36091433e-7.
```

This suggests the exact interval problem is not numerically huge.  The real
work is quotient discipline: the interval certificate must be attached to the
right labelled packet fiber, not to an unstructured list of rows.

## Robbins Versus Robin

Robbins' graph theorem says a connected graph has a strong orientation exactly
when it has no bridge.  The proof uses ear decompositions: each new ear can be
oriented without losing strong connectivity.  The LRC14 certificate analogue is
a no-bridge assembly rule:

```text
exact rational center
  -> divisor atom bank
  -> trig interval enclosure
  -> signed margin
  -> labelled packet fiber
  -> route handoff
```

If any arrow in that dependency graph is quotient-forgotten, the proof can get
stuck even though the floating scalar looked good.

Robin's number-theoretic theorem is the cautionary scalar shadow: a single
inequality involving `sigma(n)` is equivalent to the Riemann hypothesis for
`n>=5041`.  That makes `sigma` powerful, but not safe as a quotient in this
project.  The divisor-function page records that divisor sums are
multiplicative pushforwards of prime-power fibers and that their Dirichlet
series, Lambert series, Eisenstein coefficients, and Ramanujan expansions all
re-express the same divisor data.  Therefore `tau`, `sigma`, and `psi` are
allowed only after their kernel is controlled by Ramanujan modes, endpoint
owners, Toeplitz intervals, or a residual packet label.

External source trail:

```text
https://en.wikipedia.org/wiki/Divisor_function
https://en.wikipedia.org/wiki/Ramanujan_sum
https://en.wikipedia.org/wiki/Dirichlet_convolution
https://en.wikipedia.org/wiki/Robbins%27_theorem
https://mathworld.wolfram.com/RobinsTheorem.html
```

Useful facts imported from the web pass:

- divisor functions are multiplicative divisor-fiber sums and have Dirichlet
  series such as `sum sigma_a(n)/n^s = zeta(s) zeta(s-a)`;
- Ramanujan sums are sums of powers of primitive roots of unity and give
  integer primitive-period traces;
- Dirichlet convolution is the product law on arithmetic-function packets;
- Robbins' theorem is the no-bridge strong-orientation criterion;
- Robin's theorem makes the sigma-function scalar shadow deep enough to be
  dangerous if used without its missing divisor-fiber labels.

## Guardrail Theorem Target

The admissible interval-certificate theorem should read:

```text
For every post-THM-571 labelled packet family P(S), one of the following holds:

1. P is the AP/GW boundary atom, with zero strict mass and AP/GW endpoint
   zero-sum current.
2. P admits a familywise exact interval Fejer certificate:
   Q_d(x) < 0 for a rational center x, degree d, and divisor-curried atom bank.
3. P routes to a Ramanujan/Toeplitz handoff: a primitive exact-period packet
   reduces the interval work to a late q packet such as q in {25,27,41}.
4. P creates the HYP-2908/THM-572 tournament-conflict state-lift debt.
```

This is stronger than "the bank has floating negative values" and weaker than
"one bounded degree closes all LRC14."  It is the right middle object:
finite, labelled, and compatible with the repo's repeated quotient lesson.

## Tournament Analysis

Tournament vertices are proof obligations and certificate carriers, not
runners:

```text
labelled_packet_interval_certificate
> Robbins_no_bridge_assembly
> endpoint_owner_packet_anchor
> Fejer_divisor_atom_bank
> Ramanujan_exact_period_projector
> Dirichlet_convolution_packet_law
> floating_Fejer_shadow
> Robin_sigma_scalar_shadow
> raw_divisor_count
```

The pairwise observable is retention of the LRC cover predicate, quotient
kernel control, exactness path, phase information, packet labels, and
anti-scalarization.  The S162 tournament is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3_cycles=0
hamiltonian_paths=1
```

This says the current frontier is not another raw runner tournament.  It is a
certificate-assembly tournament whose top vertex is a labelled interval proof
object.

## Next Work

1. Replace the scaffold's hard-coded `pi` enclosure with a formally cited
   interval backend such as arb/MPFI-style ball arithmetic or a Lean-imported
   rational certificate.
2. Replace per-row certificates by packet-family templates: AP/GW, K33,
   unit-petal, covering, few-apex, and selected two-swap families.
3. Use Ramanujan exact-period projectors from HYP-2979 to pre-split late
   packets before interval work.
4. Prove a Robbins-style no-bridge lemma for admissible proof quotients:
   every forgotten coordinate is either reconstructible, annihilated by a dual
   certificate, or explicitly routed to a named residual bucket.
