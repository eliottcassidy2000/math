---
title: LRC14 Hurwitz-Markov-Pell cannonball carrier
status: SYNTHESIS / finite arithmetic sidecar scout; not a proof of LRC14
session: codex-2026-06-27-S243
script: 04-computation/lrc14_hurwitz_markov_pell_cannonball_s243.py
result: 05-knowledge/results/lrc14_hurwitz_markov_pell_cannonball_s243.out
tags: [lrc14, hurwitz, markov, pell, cannonball, diophantine-approximation, controlled-forgetting, tournament-analysis]
---

# HYP-3075: LRC14 Hurwitz-Markov-Pell Cannonball Carrier

This pass merges the user's Hurwitz theorem, Markov numbers, cannonball
problem, and Pell-number prompt into the current LRC14 proof-carrier ledger.
The result is not a direct LRC proof route.  It is a warning about how
Diophantine boundary data must be carried through a quotient.

## Scout

```text
04-computation/lrc14_hurwitz_markov_pell_cannonball_s243.py
05-knowledge/results/lrc14_hurwitz_markov_pell_cannonball_s243.out
```

The script generates Markov triples by Vieta mutation through coordinate
`10000`, Pell numbers through the same range, square-pyramidal/cannonball
solutions through `n<=100000`, the first endpoint-wall Pell families from the
triangular tower archive, and a proof-carrier Tournament Analysis.

Main exact readouts:

```text
Markov-Pell intersection:
  1, 2, 5, 29, 169, 985, 5741

Fixed-2 Markov-Pell branch:
  (1, 2, 5)
  (2, 5, 29)
  (2, 29, 169)
  (2, 169, 985)
  (2, 985, 5741)

Cannonball solutions through n<=100000:
  sum_{i<=1} i^2 = 1^2
  sum_{i<=24} i^2 = 70^2
```

The important bridge is:

```text
Pell P5 = 29
Pell P6 = 70
Pell P7 = 169
P5 * P7 - P6^2 = 1
70^2 = 4900
29 * 169 = 4901
```

So the nontrivial cannonball square root `70` is a Pell carry between two
neighboring Markov-Pell wall numbers `29` and `169`.  This is exactly the kind
of signal that controlled forgetting should not scalarize: the visible square is
not the proof object; the proof object is the quadratic-unit address, carry
residue, and endpoint wall.

## Hurwitz / Markov Reading

Hurwitz's theorem says every irrational has infinitely many approximants
`p/q` with error below the `1/(sqrt(5) q^2)` threshold, and the constant is
sharp at the golden-ratio class.  Markov's refinement stratifies the early
Lagrange spectrum by Markov numbers:

```text
L_m = sqrt(9 - 4/m^2)
```

For the Markov-Pell branch, the corresponding approximation coefficient
`1/L_m` descends rapidly toward `1/3`.  This is useful LRC language, but it is
facing the wrong direction if used naively.  Classical Hurwitz/Markov finds
exceptionally good approximation times; LRC needs an anti-Bohr witness time
that avoids every forbidden integer neighborhood.

The transfer is therefore:

```text
best-approximation wall
  -> exceptional continued-fraction / height / quadratic-unit sidecar
  -> anti-Bohr endpoint ledger with named destroyed coordinates
```

This extends HYP-3062's Roth-Minkowski lattice fence.  Roth/Minkowski says do
not use volume or exponent pressure after deleting lattice/height data.
Hurwitz/Markov adds: do not use best-approximant constants after deleting the
continued-fraction period, Markov depth, or quadratic-unit branch that made the
constant sharp.

## Cannonball / Pell Reading

The triangular tower archive already showed that many endpoint coincidences are
Pell walls:

```text
m(m+1) = 2n(n+1)
m^2 = 2n^2 + 2n + 1
m^2 = n(2n+1)
```

Those are infinite quadratic-unit endpoint families.  The cannonball problem is
different: `1^2+...+24^2=70^2` is a global square-pyramidal coincidence, and
the known nontrivial solution behaves like an elliptic stop rather than a
simple infinite Pell wall.

The synthesis is still valuable:

```text
Pell endpoint wall:
  infinite shell-address/carry family

Cannonball square:
  rare scalar square whose root is a Pell carry between Markov-Pell walls

LRC14 packet:
  visible blocked/open token only becomes exact after endpoint owner, shell
  address, carry residue, and route/certificate sidecars are restored
```

That is the same lesson as HYP-2456's Beatty-Pell crossover word: the visible
token is a carry-decorated projection, not the native proof coordinate.

## Proposed HYP-2963 Fields

```text
hurwitz_markov_approximant_class
lagrange_markov_depth
continued_fraction_period_word
markov_pell_fixed_coordinate
pell_wall_unit
pell_cassini_gap
cannonball_square_pyramid_gate
endpoint_shell_address
quadratic_carry_residue
preserved_anti_bohr_predicate
destroyed_owner_or_scale_coordinate
required_sidecar_or_exit
```

These fields should only appear as sidecars inside packet rows.  They are not a
new scalar certificate.

## Proof Angle A: Markov-Depth Sidecar

The first possible proof angle is to use Markov depth as a finite sidecar for
the existing q=7 resonance wall.

HYP-2745 already writes the relevant discrepancy as three legs:

```text
G_P(p,q) = [2AB(P-A)(P-B) + 2C(P-C)] / P
A = ||p||_P
B = ||q||_P
C = ||pq||_P
```

HYP-2753 says the LRC14 full-residue wall is not just the two visible legs; the
third leg is the hidden coordinate.  Markov theory gives this a better
arithmetic name: it is a best-approximant/continued-fraction depth sidecar,
not a slope-only quotient.  Future q=7 residual ledgers should therefore
record `lagrange_markov_depth` and `three_leg_residue` only when endpoint owner,
route, exact scale, and certificate exits are still present.

## Proof Angle B: Cannonball/Pell Wall Ledger

The second possible proof angle is a wall-ledger theorem:

```text
Whenever a visible scalar token is explained by a Pell/cannonball coincidence,
replace the token by shell address + quadratic carry + endpoint atom before
running status/route purity.
```

This is directly applicable to the proposed Q27 address ledger from HYP-2456
and to the HYP-3074 route-state closure median interface.  A denominator or
blocked token is not proof-bearing until the carry coordinate explains whether
the row is an endpoint atom, a neighboring open row, a deletion target, or a
named residual debt.

## Tournament Analysis

Candidate vertex sets considered included runners, speed gaps, fixed circle
sections, section boundaries, wall-crossing events, residues, cover arcs,
Fourier modes, matroid circuits, proof obligations, Markov triples, Pell wall
atoms, cannonball square-pyramid events, and proof-carrier ledgers.

Chosen vertices are proof carriers and arithmetic sidecar types, not runners.
Pairwise observable: retained critical LRC axes, total retained payload, and
destroyed owner/route/anti-Bohr coordinates.  Switch/gauge: orient toward the
carrier preserving more critical LRC axes; tie by total payload and fewer
destroyed coordinates.

Stored fingerprint:

```text
score_order:
  labelled_lrc_packet_ledger
  route_state_closure_median
  cross_carrier_portfolio
  beatty_pell_endpoint_wall
  markov_three_leg_resonance
  markov_pell_fixed_two_branch
  hurwitz_threshold
  cannonball_square_pyramid_gate

score_hist = {7:1, 6:1, 5:1, 4:1, 3:1, 2:1, 1:1, 0:1}
directed_3cycles = 0
hamiltonian_path_count = 1
```

Assumption challenge: the tempting vertices are Markov numbers or Pell numbers.
That destroys the LRC predicate.  The quotient preserves only arithmetic wall
type.  To preserve the LRC predicate it must be pulled back to packet rows with
endpoint owner, exact scale, route, legal exit, and certificate sidecars.

## Conclusion

Hurwitz/Markov/Pell/cannonball data should enter the LRC14 program as an
exceptional-wall address system.  It is especially useful for naming when a
visible scalar is an endpoint wall or carry projection.  It is not a direct
proof quotient because it forgets the anti-Bohr predicate, endpoint owner, and
route/certificate payload needed to prove LRC14.
