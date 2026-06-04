# HYP-2215: Side-Channel Jackknife Is the Shared Transfer Move

**Status:** OPEN synthesis / method hypothesis.

## Claim

The current threads around `H=21`, LRC `n=14`, unit-distance `n=21`,
A000568, Schanuel/pi-e, twin primes, and Goldbach share one proof shape:

```text
raw scalar quotient + hidden witness + retained side channels.
```

The transferable method is not raw scalar matching.  It is a side-channel
jackknife: remove or compress one retained channel at a time and measure which
obstruction mass is released.

## Evidence

S639 records seven carrier rows:

- `H=21`: scalar Hamiltonian-path count; side channels are strongness, `c3`,
  SCC, OCF, and `Phi3`/self-converse structure.
- LRC `n=14`: scalar `C=27=3^3`; side channels are unit/nonunit shells, gcd
  strata, owner, carry, lift, and self-converse round classes.
- unit-distance `n=21`: scalar `21` vertices / `57` unit edges; side channels
  are `20`-edge unit spine, `37`-edge centered-hex bulk, direction support,
  traceability, and embedding.
- A000568: scalar isomorphism counts; side channels are Burnside cycle types,
  automorphisms, self-converse fixed layers, q-deformations, and transporters.
- Schanuel/pi-e: scalars `S=e+pi`, `P=e*pi`; side channels are trace, norm,
  discriminant, `log(pi)`, and transverse power sums.
- twin primes: scalar gap `2`; side channels are local residues, singular
  series, parity barrier, level of distribution, and almost-prime shadows.
- Goldbach: scalar `N=p+q`; side channels are local residues, circle-method
  major/minor arcs, singular series, and representation counts.

The transfer-method tournament in S639 is transitive and ranks

```text
side_channel_jackknife
> spine_bulk_decomposition
> finite_window_extinction
> local_prime_ledger
> transverse_shadow_fallout
> Burnside_transporters
> quotient_collapse_test
> Schanuel_style_completion
> raw_scalar_numerology.
```

## New Applications

1. For twin primes and Goldbach, run LRC-style channel jackknives over local
   moduli and residue classes before any global density heuristic.
2. For unit-distance constructions, search scalar twins: same edge count, but
   different spine/bulk/direction/traceability packets.
3. For A000568, extend hard terms by fixed, merged, q-shadow, and transporter
   companions before chasing the next term directly.
4. For famous conjectures, ask what the finite resolved face is.  H=21 closes
   by retaining strongness and `c3<=10`; pi/e gets transverse fallout by
   retaining `log(pi)` and power-sum shadows.

## Method

Given a hard conjecture:

1. Name the raw scalar quotient.
2. Name the hidden witness.
3. List side channels.
4. Jackknife the side channels.
5. Look for scalar twins with different retained packets.
6. Promote the first non-redundant retained packet to the next proof object.

**See also:** `04-computation/obstruction_carrier_atlas_s639.py`,
`05-knowledge/results/obstruction_carrier_atlas_s639.out`,
`07-reflections/obstruction-carrier-atlas-s639.md`, HYP-2214, HYP-2213,
HYP-2212, HYP-2211, HYP-2209, HYP-2206, HYP-2200.
