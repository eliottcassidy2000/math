---
id: HYP-3058
title: LRC14 Fermat-Catalan hyperbolic reciprocal bound
status: SYNTHESIS / proof-interface analogy; not a proof
source: codex-2026-06-26-S222
tangent: T1140
related:
  - HYP-3055
  - HYP-3054
  - HYP-3043
  - HYP-3039
  - HYP-3012
  - HYP-3009
  - HYP-3003
  - HYP-3002
  - HYP-2998
  - HYP-2963
  - HYP-2945
  - HYP-2937
  - HYP-2934
  - HYP-2928
  - THM-572
  - OPEN-Q-108
---

# HYP-3058: LRC14 Fermat-Catalan Hyperbolic Reciprocal Bound

## Claim

The Fermat-Catalan condition

```text
1/a + 1/b + 1/c < 1
```

should be imported into the LRC14 proof stack as a hyperbolic
reciprocal-curvature sidecar, not as a scalar shortcut and not as a black-box
number-theory theorem.

For a triple of packet orders, clocks, route moduli, or compression depths

```text
tau = (a,b,c),
```

define

```text
chi_orb(tau) = 1/a + 1/b + 1/c - 1,
kappa(tau)   = 1 - 1/a - 1/b - 1/c.
```

Then `kappa(tau)>0` is the same abstract threshold that appears in
hyperbolic triangle groups and three-cone-point orbifolds.  In LRC14 terms,
it marks a three-lane quotient whose lost payload has negative-curvature debt:
the quotient should not be flattened into a Euclidean scalar until the
observer-extension/cut payload is retained, reconstructed, annihilated,
descended, stopped at AP/GW boundary structure, or routed to named residual
debt.

## Where The Bound Appears

The same inequality is a standard boundary in several mathematical settings.

| Setting | Object | Meaning of `1/a+1/b+1/c<1` |
|---|---|---|
| Fermat-Catalan / generalized Fermat | Primitive solutions to `x^a+y^b=z^c` | The hyperbolic exponent regime where primitive solutions are conjecturally finite. |
| Triangle groups | Triangles with angles `pi/a`, `pi/b`, `pi/c` | The triangle is hyperbolic; area is `pi*kappa`. |
| Orbifolds | Sphere with three cone points `(a,b,c)` | Orbifold Euler characteristic is negative: `chi_orb<0`. |
| Coxeter / tiling trichotomy | Reflection groups and regular tessellation signatures | Spherical if the sum is `>1`, Euclidean if `=1`, hyperbolic if `<1`. |
| Seifert / Brieskorn settings | Base orbifolds `S^2(a,b,c)` and Brieskorn links | The base has hyperbolic orbifold type; the link/singularity inherits finite-vs-hyperbolic structure. |
| Singularity thresholds | Brieskorn-Pham forms such as `x^a+y^b+z^c` | The reciprocal sum is a log-canonical-style threshold crossing. |

The useful abstraction is not "these domains solve LRC14."  The useful
abstraction is that a three-direction compression has a curvature sign, and
the sign says whether local quotient flattening is safe.

## LRC14 Transfer

HYP-3009 already treated Fermat-Catalan power-lift data as a packet field:
perfect powers on the unit-excess chain, such as `p=2` giving
`q=27=3^3`, are no-lift guardrails rather than proof certificates.  HYP-3054
and HYP-3055 supply the controlled-forgetting wrapper: before a quotient is
used by the next operation, name the observer-extension payload it would
forget.

This hypothesis adds one packet field:

```text
hyperbolic_triple_signature:
  triple_orders: (a,b,c)
  reciprocal_sum: 1/a + 1/b + 1/c
  curvature_margin: 1 - reciprocal_sum
  orbifold_euler: reciprocal_sum - 1
  regime: spherical | euclidean | hyperbolic
  triangle_signature: Delta(a,b,c)
  discharge_route: AP/GW boundary | q-witness | C27 petal | K33 state-lift |
                   Fejer/Toeplitz certificate | named F7 debt
```

Candidate LRC readings of `(a,b,c)` include:

- exact denominator/order side, such as q/Farey binding scale;
- automatic or lacunary language order, such as fibbinary/Moser/Hurwitz
  state depth;
- route-incidence order, such as C27, K33, or covering branch arity;
- observer-extension order, such as endpoint-owner cut, deletion fiber, and
  rectangle/hourglass cycle residue;
- certificate order, such as primitive-period deck, Fejer degree, and
  state-lift obligation.

The field is theorem-facing only after the three orders are declared.  A raw
triple of attractive numbers is numerology until it names the preserved LRC
predicate, the destroyed coordinate, and the discharge route.

## The `(2,3,7)` Resonance

The smallest hyperbolic triangle signature is

```text
1/2 + 1/3 + 1/7 = 41/42 < 1,
kappa(2,3,7)=1/42.
```

This is suggestive for LRC14 because the existing packet stack repeatedly
touches:

```text
14 = 2*7,
q=27 = 3^3,
3/41 near-miss scale,
C27 petal/two-block branch,
K33 state-lift branch,
AP/GW boundary atoms.
```

The proposed use is conservative.  `(2,3,7)` should be a route heuristic for
first hyperbolic debt, not a proof.  The actual LRC proof still needs exact
`M`, endpoint owners, safe topology, route labels, harmonic certificates, and
state-lift sidecars.

In the language of HYP-3055, `(2,3,7)` is another controlled-forgetting alarm:
when a three-lane packet has curvature margin `1/42`, do not collapse it to a
single scalar until the missing observer/cut payload is accounted for.

## Tournament Analysis

Vertices are proof carriers and quotient obligations, not runners:

```text
labelled_packet_sheaf
hyperbolic_reciprocal_signature
observer_extension_cut_payload
triangle_orbifold_guard
fermat_catalan_power_guard
exact_M_Farey_node
C27_petal_two_block_route
K33_state_lift_route
automaton_gap_language
raw_exponent_numerology
```

Pairwise observable: for carriers `A,B`, compare which one preserves more of
the LRC predicate tuple

```text
(boundary/open status, exact scale, endpoint owner, route schedulability,
 certificate handoff, named residual debt).
```

Binary gauge:

```text
A -> B
```

when `A` keeps every field `B` keeps and at least one additional field needed
to discharge a route/status-changing quotient fiber.  Ties are broken by the
Hamiltonian path displayed below.

Conservative retention path:

```text
labelled_packet_sheaf >
hyperbolic_reciprocal_signature >
observer_extension_cut_payload >
triangle_orbifold_guard >
fermat_catalan_power_guard >
exact_M_Farey_node >
C27_petal_two_block_route >
K33_state_lift_route >
automaton_gap_language >
raw_exponent_numerology.
```

Assumption challenge: this pass explicitly considered runners, exponent
triples, residues, exact denominator clocks, packet fibers, triangle cone
points, endpoint owners, cover arcs, Fourier/Fejer certificates, and proof
obligations.  The selected vertices are proof carriers because the LRC
predicate is not preserved by exponent data alone.

Preserved predicate: whether a compressed three-lane packet has legal
LRC14 discharge after the next observer/cut/certificate operation.

Destroyed information: exact identity of runners, endpoint-owner transfer,
safe interval topology, deletion fiber, rectangle/hourglass residue,
route label, and certificate payload unless retained elsewhere.

## Next Pull

Add the following fields to the HYP-2963 packet-ledger vocabulary when a
triple of meaningful orders is available:

```text
hyperbolic_triple_signature
reciprocal_sum
reciprocal_curvature_margin
orbifold_euler_sign
triangle_group_signature
fermat_catalan_power_guard
first_hyperbolic_debt
hyperbolic_debt_discharge_route
```

Then audit whether known residuals cluster into:

```text
spherical/easy capacity
euclidean boundary atom
hyperbolic finite-debt packet
```

without losing exact `M`, endpoint-owner, Haar/Fejer, primitive-period, or
state-lift sidecars.
