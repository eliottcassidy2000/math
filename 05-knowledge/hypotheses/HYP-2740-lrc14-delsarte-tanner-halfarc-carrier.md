---
id: HYP-2740
title: LRC14 Delsarte duals have a Tanner carrier, but the useful invariant is half-arc sign rigidity
status: OPEN guardrail and proof-order target
source: codex-2026-06-21; renumbered from a transient HYP-2738 collision after incoming mainline work
depends_on:
  - THM-534
  - HYP-2726
  - HYP-2737
related:
  - HYP-2723
  - HYP-2735
  - HYP-2638
  - HYP-2217
  - HYP-2205
---

# HYP-2740: Delsarte/Tanner/Half-Arc Carrier

## Claim

The THM-534 Delsarte dual polynomials

```text
g(t) = sum_r y_r binom(t,r),  t=0..6,
```

do define an exact Tanner-style bipartite carrier:

```text
check nodes = factorial moment rows r with y_r != 0,
variable nodes = missed-depth atoms t=0..6,
edge r--t exists iff r <= t, with weight y_r binom(t,r).
```

This carrier preserves the exact dominance predicate `g(t) >= 1[t=0]`, hence
the Delsarte/LP certificate.  But it is not an LDPC/Tanner-expander mechanism:
it is a dense nested Ferrers graph with many 4-cycles.  The prompt's Tanner
graph inspiration is therefore useful as a carrier audit, not as a sparse
local proof import.

The useful residue is instead a Doyle-Holt-like half-arc fact.  The sign data
is rigid: under all side-preserving automorphisms of the undirected
moment-depth carrier, no edge orbit mixes positive and negative weights.  The
same undirected support can be viewed without signs, but its two orientations
are not interchangeable by carrier symmetry.  This is the finite certificate
analogue of an edge-transitive object whose arcs split into two types.

## Exact Evidence

Script:

```text
04-computation/lrc14_delsarte_tanner_carrier_codex_20260621.py
```

Stored output:

```text
05-knowledge/results/lrc14_delsarte_tanner_carrier_codex_20260621.out
```

The three binding duals have exact scaled forms:

```text
K8:  g scale 10 = [10,0,0,1,0,0,10]
     y scale 10 = [10,-10,10,-9,6,0,0]
     Kraw scale 80 = [5,0,2,0,3,0,0]

K9:  g scale 18 = [18,5,0,0,2,3,0]
     y scale 18 = [18,-13,8,-3,0,0,0]
     Kraw scale 144 = [12,2,4,3,0,0,0]

K11: g scale 6 = [6,3,1,0,0,1,3]
     y scale 6 = [6,-3,1,0,0,0,0]
     Kraw scale 24 = [3,1,1,0,0,0,0]
```

Carrier diagnostics:

```text
K8:  edges=25, density=5/7,  girth=4, 4-cycles=65,
     automorphisms=6,   edge-orbits=15, mixed-sign-orbits=0.

K9:  edges=22, density=11/14, girth=4, 4-cycles=53,
     automorphisms=24,  edge-orbits=10, mixed-sign-orbits=0.

K11: edges=18, density=6/7,  girth=4, 4-cycles=35,
     automorphisms=120, edge-orbits=6,  mixed-sign-orbits=0.
```

Single-edge sign flips preserve dominance only for the negative edge family:

```text
K8:  10/25 feasible flips
K9:  10/22 feasible flips
K11:  6/18 feasible flips
```

Whole-row sign flips likewise show that negative moment rows can be flipped
without breaking the weak inequality, while positive rows are certificate
critical.  This separates "dominance feasibility" from "the intended minimal
Delsarte certificate": arbitrary flips are not a proof route, but their
asymmetry identifies the signed orientation carried by the certificate.

## External-Source Reading

The web prompts fit the following division of labor.

- Delsarte LP certificates live in association schemes and dual positivity,
  not in sparse local checks.  Samorodnitsky's Delsarte-LP optimum result is a
  warning that the LP shadow can be structurally limited if not tied to
  realizability/generator data.
- New dual LP solutions via Laplacian/association-scheme methods suggest that
  "graph energy" is a natural carrier for the certificate, but not necessarily
  an LDPC Tanner graph.
- Li's uniqueness and parity phenomena for Delsarte optima suggest the next
  real target here: explain the transition from K8 even-only Krawtchouk support
  to the mixed K9/K11 supports by a puncture/extend parity lemma.
- The antipodal/unit-distance Delsarte method confirms that unit-distance
  transfer should be spectral/Fourier-positive rather than sparse-check local.
- The Tanner graph question for quantum LDPC codes is best read as a search
  for a local carrier behind a global extremal function.  The audit says the
  local carrier exists here but is too dense and chain-like to carry expansion.
- The Doyle-Holt graph analogy survives only at the sign-orbit level: support
  symmetry and oriented/sign symmetry are different invariants.

## Proof Order

This hypothesis changes the next target from "find a Tanner theorem" to:

```text
prove a puncture/extend parity lemma for the K8/K9/K11 Delsarte duals,
then splice it to HYP-2737's generated row-slice balance.
```

More concretely:

1. Treat K8's even-only support `(K0,K2,K4)` as the even half-tiling side of
   the atom quotient.
2. Treat K9/K11's mixed low-degree supports as punctured or extended versions
   of that even certificate.
3. Show that the sign-rigid binomial carrier is compatible with the
   HYP-2737 row word after quotienting by missed depth.
4. Only after that, scalarize to the Delsarte cap.

## Root-Parity Addendum

Follow-up script:

```text
04-computation/lrc14_delsarte_root_parity_atlas_codex_20260621.py
```

Stored output:

```text
05-knowledge/results/lrc14_delsarte_root_parity_atlas_codex_20260621.out
```

This tests the puncture/extend target directly.  For each degree `d`, enumerate
all normalized root-polynomial certificates

```text
g_R(t) = prod_{a in R}(t-a) / prod_{a in R}(0-a),  R subset {1,...,6},
```

that satisfy both dominance `g_R(t) >= 1[t=0]` and Krawtchouk nonnegativity.
Abstract Delsarte feasibility is not unique:

```text
degree 2: 3 feasible root certificates
degree 3: 3 feasible root certificates
degree 4: 5 feasible root certificates
```

The generated AP occupancy law then selects exactly the THM-534 roots:

```text
k=8:       degree 4 best roots = (1,2,4,5), known_best=True
k=9,10:    degree 3 best roots = (2,3,6),   known_best=True
k=11..13:  degree 2 best roots = (3,4),     known_best=True
```

Every candidate tournament under the observable "lower `L_y` on the generated
AP row is better" is transitive, with one Hamiltonian path.  Thus the known
Delsarte duals are not selected by abstract LP positivity alone; they are
selected by the actual LRC/generated-depth occupancy distribution.  This is the
precise finite form of the proof-order warning:

```text
generated word -> depth occupancy -> Delsarte parity certificate.
```

## Tournament Analysis

Vertices are proof lenses, not runners or arcs:

```text
global_delsarte_lp
> krawtchouk_parity_puncture
> halfarc_sign_orbit
> tanner_local_carrier
> unit_distance_delsarte_transfer
> modular_belyi_address.
```

Pairwise observable:

```text
(preserves_bound, exact_audit, LRC_transfer, sign_structure, locality)
```

Switch/gauge: keep generated atom words before scalarizing; Tanner locality is
a tie-break, not the main invariant.  The computed tournament is transitive:
`directed_3cycles=0`, `Hamiltonian_paths=1`.

## Assumption Challenge

The useful vertices are not runners, tournament arcs, unit-distance points, or
Tanner checks alone.  The chosen quotient uses missed-depth atoms and moment
rows because it preserves the exact Delsarte dominance predicate.  It destroys
the actual speed-set relation lattice and the HYP-2737 generated odometer
word, so it cannot replace those proof obligations.  Its job is to classify
which signs and parities a later generated-word proof is allowed to use.
