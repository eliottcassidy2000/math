---
source: codex-2026-05-31-S450
status: exploratory synthesis with web/repo search
tags:
  - lonely-runner
  - analogies
  - anti-bohr
  - distance-graphs
  - zonotopes
  - endpoint-incidence
---

# LRC Analogy Atlas: What Is Actually Isomorphic?

This session searched the repo and the web for problems analogous to the
Lonely Runner Conjecture.  The useful answer is not "many problems involve
circles."  The useful answer is:

```text
LRC is a one-parameter subgroup trying to avoid coordinate neighborhoods,
and the hard object is the finite protected boundary of the bad cover.
```

This is a companion to the upstream S430 lens-atlas session.  HYP-1900 names a
labelled incidence core; this S450 pass sharpens the proposed invariant to the
protected anti-Bohr boundary and asks whether every serious analogue preserves
positive gaps, unprotected endpoints, and endpoint debt export.

The exact formulations preserve that object.  The looser analogies preserve
only the visual, volume, or coloring flavor.

## Current Web Status

As of this session, the frontier has moved.  Rosenfeld proved the eight-runner
case by computer-assisted finite checking:

- https://arxiv.org/abs/2509.14111

Trakulthongchai refined the approach with a sieve and obtained nine and ten:

- https://arxiv.org/abs/2511.22427

Sungkawichai and Trakulthongchai then pushed the same finite-checking/sieving
line to `k in {10,11,12}` in the Wills stationary-runner notation, i.e. the
titles' eleven, twelve, and thirteen runners:

- https://arxiv.org/abs/2604.23906

So, if these current preprints are accepted as the working frontier, the repo's
old fixation on `n=14` was not a random numerological itch.  It is exactly the
next first-open denominator after the 2026 computational progress.

## Exact Or Near-Exact Models

### 1. Endpoint arc cover

This is the repo's native model after THM-357.

For speeds `V={v_1,...,v_k}`, with `n=k+1`,

```text
F(V)=union_v {t in R/Z : ||v t|| < 1/n}.
```

The conjecture says `F(V)` never covers the whole circle.  THM-357 sharpens
that into:

```text
positive complement interval,
or unprotected boundary endpoint,
or full open cover with every endpoint protected.
```

The third case is the only dangerous one.  This is the finite protected
boundary hypergraph.

### 2. Distance graph regular circular coloring

For the integer distance graph `G(Z,D)`, a multiplier coloring

```text
c_t(x)=t*x mod 1
```

is proper at circular threshold `1/n` exactly when

```text
||d*t|| >= 1/n for all d in D.
```

Thus this is not just analogous.  It is the same multiplier search.  The
important qualifier is "regular" or "multiplier": arbitrary graph colorings
are too flexible and forget the one-parameter subgroup.

Web anchor:

- https://arxiv.org/abs/0710.4495

### 3. Subtorus/cube avoidance

The speed vector `v` gives a line or one-dimensional subtorus

```text
t -> t*v in (R/Z)^k.
```

LRC asks whether that line intersects the safe cube

```text
[1/n,1-1/n]^k modulo Z^k.
```

This is the same inequality, but shifted from a cover of the parameter circle
to a hitting problem in the speed torus.

Web anchors:

- https://arxiv.org/abs/2304.01462
- https://arxiv.org/abs/2605.27941

### 4. View obstruction

The view-obstruction formulation asks whether a ray through a periodic field
of obstacles can avoid the obstacles.  It is exact at the correct parameter
translation, but operationally it tends to emphasize line of sight rather than
the finite endpoint-protection graph.

Web anchors:

- https://www.sciencedirect.com/science/article/pii/S0095895697917706
- https://www.quantamagazine.org/new-strides-made-on-deceptively-simple-lonely-runner-problem-20260306/

### 5. Zonotope covering radius

Henze-Malikiosis and later work project the cube/line arrangement into a
lattice zonotope covering-radius question.  This is one of the strongest
external formulations for the repo because it turns endpoint debt into
covering debt.

The Cambridge finite-checking paper explicitly recalls the zonotopal
reformulation and states that a time exists for a velocity vector exactly when
the associated zonotope satisfies the relevant covering condition.

Web anchors:

- https://arxiv.org/abs/1609.01939
- https://www.cambridge.org/core/journals/forum-of-mathematics-sigma/article/linearly-exponential-checking-is-enough-for-the-lonely-runner-conjecture-and-some-of-its-variants/A51A991DE89B8C9C2E2FF13FBD4501DA
- https://arxiv.org/abs/2603.24784

### 6. Lonely Runner spectra

Giri-Kravitz define spectra using one-dimensional subtori of `(R/Z)^n` and
their `L_infinity` distance from the center.  Their result that accumulation
points of `S(n)` are `S(n-1)` resonates hard with the repo's recursive
tournament work:

```text
delete runner / drop coordinate / take spectral accumulation
```

That is not the endpoint graph itself, but it is a real moduli-space shadow of
the same problem.

Web anchor:

- https://arxiv.org/abs/2304.01462

## Good Analogies That Are Not Isomorphisms

### Nowhere-zero flow compression

Bienia-Goddyn-Gvozdjak-Sebo-Tarsi used LRC in connection with nowhere-zero
flows.  This is important, but it is more application than isomorphism.  The
flow problem shares modular avoidance and value compression, but its native
boundary object is a graph-flow conservation constraint, not the LRC endpoint
hypergraph.

Web anchor:

- https://www.sciencedirect.com/science/article/pii/S0095895697917706

### Danzer/dense forest and general visibility problems

These share the "line of sight through obstacles" smell.  They are useful for
intuition about hitting all long line segments or rays, but the quantifiers and
dimension shift.  They do not automatically preserve the one-dimensional
subgroup or endpoint-protection incidence.

### Covering systems, residue covers, and Erdos-Straus style ledgers

The repo's base-42/Erdos-Straus material is analogous at the quotient-cover
level:

```text
identity families cover residue classes;
LRC intervals cover multiplier residues.
```

This is useful because both ask what survives a first quotient cover.  It is
not an exact LRC model because Erdos-Straus has algebraic equality boundaries,
not circular-arc endpoints.

### Graph reconstruction and union-closed sets

These are farther away but share the repo's residue principle:

```text
do not study only the projection; study the boundary residue it forgets.
```

Graph reconstruction has deletion decks; union-closed sets have coordinate
frequency shadows.  They are not LRC, but they help name what we care about:
the obstruction living at the edge of a projection.

### Zeckendorf and Ostrowski normal forms

Zeckendorf decomposition belongs in this thread as a normal-form model, not as
an exact LRC formulation.  The repo's earlier bridge is:

```text
Zeckendorf = independent sets in P_infty at fugacity x=1
Tournament OCF = independent sets in Omega(T) at fugacity x=2
```

For LRC, the useful imported structure is the no-adjacent-carry rule.  If the
endpoint-protection hypergraph can be quotiented to a path-like debt automaton,
then exposed endpoints may have a unique Zeckendorf-style normal form: repairs
are carries, adjacent carries are illegal, and a proof certificate is a
nonzero normal form that cannot be annihilated.

There is also an anti-Bohr version: Ostrowski numeration for irrational
rotations encodes return times by continued-fraction denominators.  For the
golden slope, all continued-fraction digits are `1`, so the Ostrowski expansion
is exactly Zeckendorf.  The golden rotation is therefore the clean model case
for anti-Bohr boundary recursion with no adjacent carries.

See `07-reflections/lrc-zeckendorf-bridge-s451.md` and HYP-1902.

## Exact Finite Check From S450

`04-computation/lrc_analogy_atlas_s450.py` reruns small exact boundary checks:

```text
initial n=14:
  boundary_only, forbidden=1, unprotected=6, first_unprotected=1/14

n14 seven-ladder:
  positive_gap, forbidden=142/143, max_gap/th=5/924,
  unprotected=84, first_unprotected=9/98

n14 S380 gate ladder:
  positive_gap, forbidden=142/143, max_gap/th=5/1848,
  unprotected=168, first_unprotected=15/196
```

The key lesson survives every formulation: near-counterexamples get close in
volume/gap width, but they still have exposed boundary debt.  A true
counterexample would need to erase both the real gap and the protected-boundary
defect.

## What Underlying Structure Matters?

The structure that matters is not:

```text
speeds as integers,
volume alone,
chromatic number alone,
or a visual obstacle metaphor alone.
```

The structure that matters is:

```text
(one-parameter subgroup) + (forbidden neighborhoods) + (finite boundary incidence).
```

Different exact models move this triple around:

```text
circle cover:     parameter-side bad arcs
distance graph:   failed multiplier colorings
torus:            line missing a safe cube
zonotope:         covering-radius debt
spectrum:         moduli of 1D subtori and deletion accumulation
IP:               rows=endpoints, columns=speeds
Bruhat-Tits:      endpoint denominators as p-adic frontier mass
```

That is the universality class.

## New Hypothesis

HYP-1901 says the next proof object should be a functor from every exact LRC model to the
same finite protected-boundary incidence category.

Objects:

```text
critical boundary states
```

Morphisms:

```text
protection/coverage by a generator
```

A model is useful to the extent that it preserves:

```text
positive gap,
unprotected endpoint,
endpoint debt export under gates,
and deletion/coordinate recursion.
```

This suggests a practical test for future analogies: if an analogy cannot say
what an unprotected endpoint becomes, it is probably poetic but not a proof
engine.

## Next Experiments

1. Build a converter from endpoint hypergraphs to distance-graph failed
   multiplier colorings and back.
2. Add a zonotope-facing ledger: endpoint debt should become a covering-radius
   row deficit or a facet-normal debt.
3. Compare the S450 endpoint debt table with the S420 integer-program rows and
   ask whether each uncovered row has a zonotope facet interpretation.
4. Use the spectra recursion `accumulation(S(n))=S(n-1)` as the clean model for
   the repo's "delete runner / add two / double column" recursion.
5. Treat `n=14` as the first post-2026 frontier and search for a certificate
   that can be read simultaneously as endpoint peel, IP dual, and zonotope
   covering debt.
6. Build an endpoint-layer graph for known near-counterexamples and test
   whether exposed debt has a Zeckendorf/Ostrowski no-adjacent normal form.
