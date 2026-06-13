# Residue Calculus Feedback Loop

**Session:** codex-2026-05-30
**Layer:** synthesis / hypotheses, not canon
**Threads sampled:** projection defects, H=63 complete Omega, good-cut interval gas,
bucket transport, endpoint transfer, Paley codes, beta constraint propagation.

The common question behind several unrelated-looking files is:

```text
After a natural projection forgets structure, what residue remains?
```

The answer is different in each thread, but the shape is similar enough to
be useful.  A projection can kill almost everything, preserve a parity shadow,
hide mass inside fibers, or expose a quotient-boundary flow.  The thing that
matters is often not the projected invariant itself, but the defect between
the upstairs object and its shadow.

## Sampled Threads

### 1. Deletion / old-coordinate projection

The H=63 n=8 classes are exact projection kills.  All odd cycles share one
core vertex, deleting it leaves a transitive tournament, and the conflict graph
is complete:

```text
Omega(T) = K31
H(T) = I(K31,2) = 1 + 2*31 = 63
```

The THM-025 n=9 real-rootedness counterexample is the adjacent near-kill:

```text
94 odd cycles -> delete vertex 3 -> 2 odd cycles
rho_3 = 92/94 = 0.979
residue alpha(T-v) = [1,2,1]
```

That two-cycle residue is small but structurally poisonous: it preserves the
independent triple responsible for the non-real-rooted independence polynomial.

### 2. Good-cut height / interval projection

Good cuts begin as a fixed-base-path statistic on the staircase tiling grid.
THM-349 collapses the raw bucket count to a one-dimensional interval-cover gas
on the cut path:

```text
B_N(x) = B_{N-1}(x) + sum_{L=2..N} c_L x^L B_{N-L-1}(x)
```

Then THM-353 identifies the top good-cut bucket with strong connectivity in the
concrete staircase tournament family.  So the same statistic has two readings:

```text
interval-cover partition function
graded absence of base-path cut obstructions
```

The residue here is not a leftover cycle family; it is the set of missing cut
obstructions and the transport mass that crosses between good-cut buckets.

### 3. Quotient transport / bucket balance

THM-350 and THM-351 turn bucket transport into a finite conservation law:

```text
2*internalLineCount_b + crossHalf_b = |fiber_b| * |moves|
```

This is the transport version of projection defect.  A bucket map forgets the
actual tiling, a move system perturbs the upstairs cube, and the residue is the
escaping half-line mass.  Internal mass is paired by a fixed-point-free
involution; cross mass is what the quotient cannot hide.

This is why bucket balance feels like a checksum.  The formula does not tell us
which neighboring buckets are important, but it tells us how much leakage must
exist before any interpretation of spine/ribs/sea geometry is allowed.

### 4. Paley / interval / code shadow

Paley T7 and interval T7 are deliberately awkward for a single invariant:

```text
same score sequence: regular
same odd-cycle support count: 36
different directed-cycle multiplicity: 80 vs 59
different disjoint pairs: 7 vs 14
different even-graph projection: 14 edges of degree 4 vs 7 edges of degree 2
```

The projection to supports says they are the same size of shadow; the fibers
over the supports say otherwise.  The projection to GF(2) code data makes the
same split visible in another language: the Paley adjacency gives the simplex
or Golay-side code at p=7 and p=23, while p=11 is full-rank/trivial in the
raw parity-check view.

So Paley is not merely "more symmetric."  It is flatter across several
projections at once: deletion losses are uniform, scores are uniform, and its
code shadow is classical when the rank defect is nontrivial.

### 5. Homology / constraint propagation

The beta and metallic-mean thread is less settled, but it belongs in the same
map.  The constraint complex asks whether a cycle survives after quotienting
by boundaries.  The metallic-mean notes did not find a direct spectral
quantity, but the useful object may be a propagation residue:

```text
constraint rows upstairs
dependency/boundary quotient downstairs
rank defect or beta_1 as surviving residue
```

This is the homological cousin of deletion loss and transport leakage.

## Working Dictionary

| Thread | Upstairs object | Projection | Residue |
|---|---|---|---|
| OCF / Omega | directed odd cycles | delete vertex or support set | kept cycles, support multiplicity, alpha-vector |
| H=63 | complete odd-cycle family | delete core vertex | zero residue |
| THM-025 | near-complete odd-cycle family | delete high-loss vertex | two-cycle disjoint residue |
| Good cuts | staircase tilings | cut-path interval union | missing cuts, bucket height |
| Bucket balance | Boolean cube half-lines | quotient bucket map | crossHalf leakage |
| Endpoint transfer | endpoint-extension cube | child quotient class | parity boundary / SC indicator |
| Paley codes | QR tournament matrix | GF(2) row/null space | rank defect, weight distribution |
| Path homology | chain complex | boundary / deletion quotient | beta/rank defect |

## Hypotheses Generated

### HYP-1781: Projection-Transport Residue Calculus

For tournament computations on fixed-path tilings, a combined residue vector
should explain more anomalies than any single invariant:

```text
R(T) = (
  deletion-loss profile,
  support-multiplicity profile,
  Omega alpha-vector,
  even-graph projection where defined,
  quotient-transport leakage by selected mask families,
  endpoint-transfer parity row
)
```

Expected use: prefilter real-root failures, H-gap exceptions, and path-homology
anomalies before running full expensive invariants.

### HYP-1782: Single-Core Target Spectrum

For complete-core tournaments obtained by inserting one core vertex into a
transitive tournament, the weighted signature value `r_core(s)` controls which
complete-Omega Hamiltonian counts can occur:

```text
H = 1 + 2*r_core(s)
```

The observed absence of `r_core=3` and `r_core=10` through `m=40` suggests the
complete-core mechanism cannot realize H=7 or H=21, while `r_core=31` first
appears at `m=7`, exactly where H=63 unlocks in n=8.

### HYP-1784: Flat-Versus-Localized Residue Duality

Extremal structured tournaments separate into two residue regimes:

```text
Paley-like: flat residue across deletion/code/score projections.
Core-like: localized residue concentrated at one kill or near-kill vertex.
```

The two regimes can both maximize or expose special H-values, but for opposite
reasons.  Paley spreads cycle mass evenly; H=63 concentrates it so hard that
Omega becomes complete.

## Feedback Loop For The Next Session

1. Build a small `residue_vector(T)` utility by composing existing
   `projection_defect_bridge_s12.py`, `goodcut_transport_excess_s15.py`, and
   `tournament_tda.py` feature blocks.  The first seed table now lives in
   `04-computation/residue_vector_codex_2026_05_30.py`.
2. Test whether known THM-025-style real-root failures have high
   `max_v rho_v` or a distinctive small `alpha(T-v)` residue.
3. Compare Paley, interval, near-transitive, and H-maximizer families by the
   two axes:

```text
localization = max_v loss_v / |Omega|
mobility     = normalized quotient-transport crossHalf
```

4. Try to prove the first non-computational lemma in HYP-1782: characterize
   the image of the single-core weighted signature function modulo small
   targets, especially the permanent absence of 3 and 10.

The productive move is to stop asking whether these threads are "really the
same theorem."  They are not.  They are the same diagnostic shape: a projection,
a hidden fiber, and a residue that remembers what the projection forgot.
