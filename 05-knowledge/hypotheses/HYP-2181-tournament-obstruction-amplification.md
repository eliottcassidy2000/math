---
id: HYP-2181
status: OPEN synthesis from S615; shared carrier structure supported, exact shared 1.014 tournament exponent unproved
source: user-2026-06-03; codex-2026-06-03-S615
tags: [tournaments, forbidden-H, LRC, unit-distance, Collatz, amplification, 1.014, coimage, side-channels]
---

# HYP-2181: forbidden tournament values are carrier obstructions that propagate into LRC, unit distance, and Collatz

The hard part shared by LRC, the unit-distance problem, Collatz, and the
tournament `H` spectrum is not a raw graph or a raw quotient.  It is a retained
carrier with a small forbidden arithmetic face.  In tournaments, the visible
face is the permanent gap

```text
H(T) != 7, 21.
```

In LRC, the corresponding numbers reappear as shell moduli and residue-depth
obstructions: `C=2n-1=21` at `n=11`, `C=27=3^3` at `n=14`, and the new HYP-2177
law says the doubling sporadic is governed by

```text
gcd(n-2, 2n-1) = gcd(3, 2n-1).
```

So the numbers `7`, `21`, and `27` should be treated as obstruction carriers,
not isolated coincidences.

## What is supported

The repo's solved tournament face gives the template.  The OCF packages
Hamiltonian-path counting as an independence-polynomial evaluation at `2`.
The values `7 = Phi_3(2)` and `21 = 3 Phi_3(2)` are the first two permanent
forbidden values, and existing notes now treat them as the only known permanent
gaps.
The rebased HYP-2180 mechanism thread gives the cleanest tournament explanation
so far: `H` is multiplicative over strong components, so achievable values come
from the semigroup of strong-tournament `H` values; `7` is prime and not strong,
while `21=3*7` inherits the missing factor.  HYP-2179 supplies a complementary
finite bridge: `7` and `21` remain absent in exhaustive/sampled tournament
spectra and are excluded by round LRC worry-set tournaments.  The S615
transfer-atlas artifact `HYP-2178-tournament-gap-transfer-exponent.md` recasts
the same fact as a coimage side-channel certificate.

The LRC side has become sharply compatible with this.  S613 shows that proof
burden follows `C=2n-1` residue depth: `C=23` is clean, but `C=21=3*7` grows
the survivor burden.  HYP-2177 then proves the doubling-sporadic mod-3 law on
the pinch lattice by the Euclidean identity above.  Thus `3` controls the
doubling face, while `7` appears as the first squarefree companion in `C=21`.
Incoming S612b adds the single-swap classification: the classic sporadics
are reflection-swaps, while the infinite `V*` family is the unique AP swap whose
mirror remains in the AP.  This strengthens the carrier/coimage reading because
the visible collapsed row is not enough; the retained mirror channel determines
which sporadic mechanism persists.

The unit-distance side supplies the amplification model.  Sawin's explicit
lower bound gives more than `n^1.014` unit-distance pairs infinitely often, and
the construction works by using a high-dimensional arithmetic carrier
(number fields with many small-norm primes) before projecting to the plane.
This mirrors the repo's coimage lesson: the visible planar graph is not the
master object.

The Collatz side supplies the residual-gap model.  A nontrivial cycle would
force a two-block equation involving `2^E - 3^k`; density and drift see the
bulk, while the hard residue is a low-rank arithmetic gap.

## What is not yet supported

S615 did not find a public or in-repo theorem saying that tournament
Hamiltonian-path growth has exponent `1.014`.  Classical maximum-`H` tournament
results instead say `max H(T)` is of order `n!/2^(n-1)`, with Alon's polynomial
upper slack and random-regular constant-factor lower asymptotics.  Strong
tournament minimum-`H` results have base `5^(1/3)`, again not `1.014`.

Therefore the phrase "shared 1.014 exponent" should currently be read as a
falsifiable research target:

1. find a constrained tournament carrier whose amplification exponent is
   genuinely `1.014`, or
2. demote the equality and keep the real shared structure: arithmetic carriers
   whose visible coimages miss a thin amplification or obstruction channel.

## Next proof targets

1. Build a carrier tournament for the unit-distance CM/class-field construction:
   vertices should be split primes, norm-one generators, projection coordinates,
   or proof obligations, not planar points.
2. Revisit the `H=7/21` proofs as a reusable obstruction schema: identify the
   exact low-order OCF packet that cannot be realized by a tournament conflict
   graph, then ask which LRC shell ledgers realize the same packet.
3. Compare `C=21=3*7` and `C=27=3^3` as two different propagation modes:
   squarefree product burden versus prime-power carry depth.
4. Search for an exponent only after fixing the carrier.  Candidate observables:
   number of compatible proof obligations, density of nontrivial SCCs in the
   carrier tournament, CM split-prime incidence growth, or the decay rate of
   forbidden OCF packets.

## Artifacts

- `04-computation/tournament_obstruction_amplification_s615.py`
- `05-knowledge/results/tournament_obstruction_amplification_s615.out`
- `07-reflections/tournament-obstruction-amplification-s615.md`
