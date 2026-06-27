---
id: HYP-3075
title: LRC14 Hurwitz-Markov-Pell modular cusp carrier
status: SYNTHESIS / finite arithmetic and modular-cusp sidecar scout; not a proof of LRC14
session: codex-2026-06-27-S243
updated: codex-2026-06-27-S245
tangents: [T1158, T1161]
techniques: [LTI-223, LTI-226]
tournament_techniques: [LTT-121, LTT-124]
scripts:
  - 04-computation/lrc14_hurwitz_markov_pell_cannonball_s243.py
  - 04-computation/lrc14_modular_cusp_pochhammer_hurwitz_s245.py
results:
  - 05-knowledge/results/lrc14_hurwitz_markov_pell_cannonball_s243.out
  - 05-knowledge/results/lrc14_modular_cusp_pochhammer_hurwitz_s245.out
tags: [lrc14, hurwitz, markov, pell, cannonball, q-pochhammer, modular-functions, cusp, diophantine-approximation, controlled-forgetting, tournament-analysis]
related:
  - HYP-3083
  - HYP-3077
  - HYP-3076
  - HYP-3074
  - HYP-3073
  - HYP-3072
  - HYP-3071
  - HYP-3062
  - HYP-3060
  - HYP-3058
  - HYP-3009
  - HYP-2963
  - HYP-2745
  - HYP-2753
  - HYP-2456
  - HYP-2454
  - HYP-2627
  - HYP-2617
  - HYP-2614
  - THM-538
  - THM-572
  - OPEN-Q-108
---

# HYP-3075: LRC14 Hurwitz-Markov-Pell Modular Cusp Carrier

This hypothesis now has two compatible layers.

The S243 layer merges Hurwitz's theorem, Markov numbers, Pell walls, and the
cannonball problem into the LRC14 proof-carrier ledger.  The result is not a
direct proof route.  It is a warning about how Diophantine boundary data must
be carried through a quotient.

The S245 layer adds the user's q-Pochhammer and full-modular-group cue.  A
modular function for the full modular group is invariant under modular
transformations and has a meromorphic expansion at the cusp `i infinity`; with
`q=exp(2*pi*i*z)`, its Laurent q-expansion has only finitely many negative
power terms.  In LRC language, those finitely many negative terms are named
cusp debt.  They must be retained before any q-series tail is treated as a
theorem carrier.

The shared principle is finite address before infinite tail:

```text
rare scalar coincidence
  -> retain its arithmetic/cusp address
  -> then use the forced infinite tail or mutation orbit
```

## Scouts

```text
04-computation/lrc14_hurwitz_markov_pell_cannonball_s243.py
05-knowledge/results/lrc14_hurwitz_markov_pell_cannonball_s243.out

04-computation/lrc14_modular_cusp_pochhammer_hurwitz_s245.py
05-knowledge/results/lrc14_modular_cusp_pochhammer_hurwitz_s245.out
```

The S243 script generates Markov triples by Vieta mutation through coordinate
`10000`, Pell numbers through the same range, square-pyramidal/cannonball
solutions through `n<=100000`, the first endpoint-wall Pell families from the
triangular tower archive, and a proof-carrier Tournament Analysis.

Main exact S243 readouts:

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
of signal that controlled forgetting should not scalarize: the visible square
is not the proof object; the proof object is the quadratic-unit address, carry
residue, and endpoint wall.

The S245 script packages three address/tail mechanisms:

```text
(q;q)_infty nonzero terms through q^32:
(0,1), (1,-1), (2,-1), (5,1), (7,1), (12,-1), (15,-1), (22,1), (26,1)

1/(q;q)_infty = 1,1,2,3,5,7,11,15,22,30,42,56,77,...

Delta/q=(q;q)_infty^24 begins:
1, -24, 252, -1472, 4830, -6048, -16744, 84480

j(q) = q^-1 + 744 + 196884 q + 21493760 q^2 + ...
j_principal_part = { -1: 1 }
j^2_principal_part = { -2: 1, -1: 1488 }
```

For the Hurwitz surface

```text
w^2 + x^2 + y^2 + z^2 = w*x*y*z
```

the same grammar appears in Diophantine form.  From seed `(2,2,2,2)`, the Vieta
move

```text
x_i -> product(other_three) - x_i
```

preserves the equation.  Through max entry `500`, the scout sees

```text
(2,2,2,2)
(2,2,2,6)
(2,2,6,22)
(2,2,22,82)
(2,6,22,262)
(2,2,82,306)
```

Every row is checked on the equation.  Again, the theorem-facing datum is not
the raw product/sum equality; it is the seed plus legal mutation word.

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

## Modular Cusp / q-Pochhammer Reading

The q-Pochhammer cue sharpens the sidecar rule.  The product
`(q;q)_infty=product_{n>=1}(1-q^n)` is sparse by the generalized pentagonal
theorem, while its reciprocal is the partition function tail.  The eta/Delta
lane then turns that tail into a finite-principal-part statement:

```text
Delta = q*(q;q)_infty^24
1/Delta has principal part q^-1
j has principal part q^-1
```

For the full modular group, meromorphicity at the cusp means there is only a
finite negative-power budget.  The infinite positive q-tail can be forced by
modularity only after the finite principal part has been named.  This is the
modular analogue of a Vieta orbit from a finite Hurwitz seed.

LRC translation:

```text
finite principal part / seed present
  -> q-Pochhammer, eta/Delta, j-tail, Hurwitz, Markov, or Pell tail may be used

finite principal part / seed absent
  -> raw q-series or raw Diophantine equality is only a scalar shadow
```

This directly complements HYP-3076.  HYP-3076 says equal sixth-power scalars
must declare native support-six status, padding status, owner gcd, and residue
phase.  HYP-3075 says q-series and Hurwitz/Pell scalars must declare finite
cusp principal part, seed, mutation, and arithmetic-address data.

## Proposed HYP-2963 Fields

Add the following fields to any HYP-2963 packet, modular/analytic proof route,
or Diophantine sidecar that wants to use q-series, Hurwitz/Markov/Pell data, or
rare scalar coincidence data:

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
modular_cusp_principal_part_order
finite_negative_power_budget
principal_part_coeff_vector
q_pochhammer_tail_signature
eta_delta_denominator_lane
j_rational_function_address
hurwitz_vieta_seed
hurwitz_mutation_depth
cusp_tail_discharge_route
preserved_anti_bohr_predicate
destroyed_owner_or_scale_coordinate
destroyed_cusp_or_owner_coordinate
required_sidecar_or_exit
```

These fields should only appear as sidecars inside packet rows.  They are not a
new scalar certificate.  HYP-3077's median-hull scheduler gives the natural
next interface: finite principal parts and seed/mutation words are sidecar
coordinates whose majority centers can become legal scheduler centers, but a
typed scheduler center is still nonterminal until it names a specific discharge
atom or residual debt.

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
and to the HYP-3074/HYP-3077 route-state median interfaces.  A denominator or
blocked token is not proof-bearing until the carry coordinate explains whether
the row is an endpoint atom, a neighboring open row, a deletion target, a
scheduler center, or a named residual debt.

## Proof Angle C: Modular Principal-Part Ledger

The modular analogue is a finite cusp-debt theorem:

```text
Whenever an analytic packet uses a modular or q-Pochhammer tail, first replace
the scalar tail by principal part + eta/Delta lane + exact packet owner.
```

This does not prove LRC14 by modular functions.  It says that if a proof route
uses partition tails, eta-products, Dedekind/cotangent sums, or modular
invariance language, the finite negative-power budget must be part of the row.
Only then can the tail be tested by cycle-image generation, median sidecar
closure, AP/GW boundary stopping, or THM-572/F7 debt.

## Tournament Analysis

Candidate vertex sets considered included runners, speed gaps, fixed circle
sections, section boundaries, wall-crossing events, residues, cover arcs,
Fourier modes, matroid circuits, proof obligations, Markov triples, Pell wall
atoms, cannonball square-pyramid events, q-coefficients, modular functions,
Hurwitz quadruples, finite principal parts, seed/mutation addresses, and
proof-carrier ledgers.

Chosen vertices are proof carriers and sidecar types, not runners,
q-coefficients, or raw Hurwitz quadruples.  The quotient preserves whether an
arithmetic or modular tail is legally controlled and deliberately forgets
individual coefficients unless they are retained as principal-part debt.

S243 pairwise observable: retained critical LRC axes, total retained payload,
and destroyed owner/route/anti-Bohr coordinates.  Switch/gauge: orient toward
the carrier preserving more critical LRC axes; tie by total payload and fewer
destroyed coordinates.

S243 stored fingerprint:

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

S245 pairwise observable:

```text
retained_cusp_address
finite_principal_part_order
tail_control
hurwitz_mutation_legality
preserved_lrc_predicate
```

S245 switch/gauge:

```text
A -> B when A retains the finite cusp/arithmetic address needed before B's
scalar tail can be used.
```

S245 transitive path:

```text
labelled_lrc_packet_sheaf
> modular_cusp_principal_part
> full_modular_group_invariance_gate
> q_pochhammer_eta_tail
> hurwitz_vieta_seed_orbit
> continued_fraction_markov_address
> pell_wall_unit_address
> raw_q_series_coefficients
> raw_hurwitz_scalar
```

S245 fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
```

Assumption challenge: the tempting vertices are Markov numbers, Pell numbers,
q-coefficients, modular functions, Hurwitz quadruples, or runners.  Those
destroy the LRC predicate.  To preserve the LRC predicate, the arithmetic or
modular address must be pulled back to packet rows with endpoint owner, exact
scale, route, legal exit, finite principal part or seed, and certificate
sidecars.

## Next Pull

Run this sidecar on an actual HYP-2963 residual family with analytic,
q-series, or Diophantine language attached.  A row is proof-safe only if it
records:

```text
exact_M_or_Farey_address
modular_cusp_principal_part_order
principal_part_coeff_vector
q_pochhammer_tail_signature
hurwitz_vieta_seed_or_mutation_depth
continued_fraction_or_pell_address
cycle_class_image_status
route_state_closure_status
median_center_kind
cusp_tail_discharge_route
```

The desired lemma shape is: once AP/GW boundary atoms, finite principal parts,
and named Hurwitz/Pell addresses are separated, every remaining modular or
Diophantine tail either descends through a legal packet family, is generated by
a certificate/cycle image, lands in a legal median scheduler center with a
specific refinement, or becomes explicit THM-572/F7 debt.

## Conclusion

Hurwitz/Markov/Pell/cannonball and q-Pochhammer/modular-function data should
enter the LRC14 program as exceptional-wall address systems.  They are useful
for naming when a visible scalar is an endpoint wall, carry projection, finite
principal part, or seed-generated tail.  They are not direct proof quotients
because they forget the anti-Bohr predicate, endpoint owner, cusp debt, and
route/certificate payload needed to prove LRC14.
