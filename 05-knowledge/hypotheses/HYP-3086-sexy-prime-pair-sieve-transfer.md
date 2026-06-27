---
id: HYP-3086
title: Sexy prime conjecture as the h=3 pair-sieve ray and what the LRC/Tournament Analysis stack transfers
status: SYNTHESIS / proof-interface map; not a proof
source: codex-2026-06-27-S254
tangent: T1168
related:
  - HYP-1964
  - HYP-1965
  - HYP-1966
  - HYP-2982
  - HYP-3032
  - HYP-3083
  - HYP-3084
  - THM-503
  - THM-518
  - THM-573
---

# HYP-3086: Sexy prime pair sieve transfer

## Claim

The repo's existing prime-pair lens already places the sexy prime conjecture in
the right coordinate system: sexy prime pairs are the fixed half-gap ray

```text
(p, p+6) = (m-3, m+3), so h=3 and gap=2h=6.
```

The transferable contribution from the LRC/Tournament Analysis stack is not a
proof of infinitely many such pairs. It is a finite-address proof interface:
keep the midpoint, gap ray, local residue ledger, sieve weight, prime-power
sidecars, and remaining proof debt visible until a genuine lower-bound
prime-pair theorem has paid for them.

## Local pair-sieve geometry

For a fixed half-gap `h`, the pair `(m-h, m+h)` is killed modulo an odd prime
`q` when

```text
m = h mod q or m = -h mod q.
```

For the twin ray `h=1`, these are two bad midpoint residues for every odd
prime. For the sexy-prime ray `h=3`, the `q=3` obstruction collapses:

```text
m = 3 mod 3 and m = -3 mod 3 both mean m = 0 mod 3.
```

That collapse changes the Hardy-Littlewood fixed-gap chord by

```text
(3-1)/(3-2) = 2
```

relative to twins. This exactly matches the older HYP-1966/HYP-1965 pair-lens
readout: gap `6` should occur at about twice the twin-pair frequency after the
same global scaling. The new S254 scout through `10^6` gives:

```text
h=1 gap=2 count=8169
h=3 gap=6 count=16386
count(gap 6)/count(gap 2)=2.006
predicted chord ratio=2.000
```

This closes the local admissibility and singular-series-shape part of the
story. It does not close infinitude.

## What transfers from the repo

1. Pair-first ontology. HYP-1965 and HYP-1966 say a prime-pair proof should not
   begin with isolated prime endpoints. The natural cell is
   `(midpoint, gap, local residue ledger, singular multiplier, carry debt)`.

2. Controlled forgetting. The LRC finite-address discipline says a quotient may
   forget coordinates only after naming the predicate preserved and the
   destroyed data. For sexy primes the preserved predicate is:

   ```text
   infinitely many m with m-3 and m+3 both prime.
   ```

   Destroyed data include endpoint identities, residue channels modulo local
   primes, prime-power channels, squarefree blindness, and exceptional
   distribution moduli.

3. Analytic-sieve guardrails. HYP-2982/HYP-3032 warn that `mu^2/phi`, large
   sieve clocks, and related weights are capacity meters unless a proof also
   records what they cannot see. For sexy primes that means squarefree weights
   cannot be allowed to hide prime-power and repeated-prime obstructions, and
   upper-bound sieve success cannot be confused with simultaneous primality.

4. Side-channel jackknife. The obstruction-carrier atlas style is useful here:
   delete or retain one local prime, one gap ray, or one weight sidecar and
   measure whether the gap-6 lower-bound predicate survives. This localizes
   which part of the proof is still doing unpaid work.

5. Recent LRC14 level-7 and gK8 work transfers only methodologically. Incoming
   THM-573 sharpens the LRC covering route to a level-7 lift sieve and reduces
   the residual to rows with at most six multiples of `7`; incoming HYP-3084
   refutes the literal Clebsch-design covariance shortcut and replaces it with
   a reflection-symmetric pairwise co-emptiness matrix whose Perron mode should
   be bounded through a `3x3` symmetric block. These are strong examples of
   local gates plus finite residual ledgers plus low-order certificates. A
   sexy-prime proof would need an analogous ledger, but for actual prime-pair
   distribution rather than LRC covering rows.

## What does not transfer

The LRC singular-series material is not the Hardy-Littlewood prime-pair Euler
product. THM-503 and THM-518 emphasize that the LRC singular object is an
archimedean or geometric singular integral, while fixed prime gaps are governed
by local congruence obstructions and an Euler product. Treating these as the
same object would be an unsafe quotient.

Likewise, the existing LRC proof infrastructure does not break the parity
barrier. It does not by itself prove that two sifted endpoints are both prime,
nor does it supply the required level of distribution for a lower-bound sieve
on the exact `h=3` ray. Maynard/GPY-style bounded-gap machinery gives bounded
gaps under deep distributional input, but it does not force the exact gap `6`
without a stronger fixed-gap prime-pair theorem.

The incoming level-7 and Perron-covariance refinements make this distinction
sharper, not weaker. THM-573 is a pigeonhole lift over seven phases for LRC
covering rows, while HYP-3084 is a low-dimensional moment/eigenvalue route for
sector co-emptiness. Neither supplies equidistribution of primes in arithmetic
progressions or a parity-breaking prime-pair detector.

## Remaining proof obligations

The local sexy-prime packet is now clear, but the proof still needs:

1. A lower-bound sieve specialized to the fixed ray `h=3`.
2. A mechanism that breaks parity and detects two primes, not merely
   almost-primes or unsieved survivors.
3. Uniform distribution of primes in arithmetic progressions at a strength
   sufficient for the fixed-gap lower bound, or a replacement structural input.
4. A finite-address ledger for every quotient used in the argument: local
   residue channel, prime-power sidecar, weight type, exceptional modulus,
   parity debt, distribution debt, and terminal exit.

## Tournament Analysis

Candidate vertex sets considered: prime endpoints, prime pairs, midpoints,
fixed gap rays, local residue channels, sieve weights, almost-prime witnesses,
exceptional moduli, and proof obligations.

Chosen vertices for the useful tournament: proof obligations, fixed-gap pair
rays, and residue/weight channels. This preserves the exact gap-6 infinitude
predicate while avoiding the false assumption that raw primes or raw pair
counts are the right proof atoms.

Pairwise observable:

```text
local obstruction collapse,
Hardy-Littlewood chord,
observed fixed-gap row count normalized to twins,
squarefree blindness,
prime-power sidecar retention,
parity debt,
distribution debt.
```

Switch: orient `A -> B` when `A` preserves the exact `h=3` predicate with
fewer unnamed residue, weight, parity, or distribution debts than `B`.

Tie Hamiltonian path:

```text
h3_lower_bound_sieve
> parity_breaking_prime_detection
> distribution_in_AP_input
> local_residue_ledger
> squarefree_primepower_sidecars
> pair_surface_counts
> raw_gap_count
```

Assumption challenge: the vertices should not be runners, arcs, or individual
prime endpoints by default. A quotient to gap counts preserves the visible
frequency signal but destroys the midpoint residue ledger and the sieve-weight
debt; a quotient to local residues preserves admissibility but destroys global
distribution. The proof object has to keep both until a lower-bound theorem
discharges them.

## Next artifact

Build a `sexy_prime_pair_ledger` with fields:

```text
midpoint_m
gap=6
local_bad_residue_word
surviving_primorial_class
hardy_littlewood_chord
sieve_weight_type
almost_prime_status
parity_debt
distribution_modulus
exceptional_modulus_flag
prime_power_sidecar
carry_debt
terminal_exit
```

The test for progress is strict: the ledger should say exactly which proof
obligation remains unpaid after local admissibility and the gap-6 singular
series multiplier are accounted for.
