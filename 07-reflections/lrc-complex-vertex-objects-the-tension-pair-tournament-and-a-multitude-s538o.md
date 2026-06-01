---
source: oracle-2026-06-01-S538o
status: synthesis + computation (complex/precise vertex objects for the LRC tournament; the tension-pair tournament; honest gap negative; a posed multitude)
tags:
  - lonely-runner
  - tournament-vertices
  - tension
  - cocycle
  - difference-set
  - three-gap
  - oriented-matroid
  - multitude
---

# Complex, Precise Vertices for the LRC Tournament: the Tension-Pair Tournament, an Honest Gap Negative, and a Posed Multitude

Pushing the "what are the vertices?" question into more intricate, precise objects.
Two prior lessons constrain the search: **restriction comes from retained metric**
(S535), and **the speeds live in a flow/tension duality** (S537). So the fertile
complex vertices carry *derived* metric/algebraic data with built-in **consistency
laws**. One genuinely new construction lands cleanly; one hoped-for one fails
honestly; and a multitude is posed precisely.

## The star: vertices = pairs, carrying difference-speeds (a tension)

Let the vertices be the `C(n,2)` **pairs** (edges of `K_n`, observer included), and
decorate pair `{i,j}` with its **difference-speed** `w_{ij} = v_i - v_j`. The
half-turn comparator on the difference-phases `frac(w_{ij} t)` makes a digraph on the
pairs — a **second-order LRC on the difference set** `{w_{ij}}`.

The key structure is a **consistency law**: the `w_{ij}` are a **tension /
coboundary**, obeying the cocycle identity
```
w_{ij} + w_{jk} + w_{ki} = 0   on every triangle      (verified True).
```
So the `C(n,2)` difference-speeds are **not free**: they are determined by the `n-1`
potentials (speeds), leaving `C(n-1,2)` cocycle constraints. This restricts the
realizable pair-structures sharply (`lrc_complex_vertex_objects_s538.py`):

- **n=5 (10 pair-vertices):** realizable **labeled** pair-tournament types =
  **24,535** out of `2^{45} ≈ 3.5·10^{13}` — a vanishing fraction (`~7·10^{-10}`).
  The cocycle (additive consistency) is the restriction.
- **n=4 (6 pair-vertices):** the structure is a tension-valued **digraph with ties**
  (two pairs with equal `w` always tie), so it is *not* a plain tournament — the
  realized `74` digraph-types are not comparable to `A000568(6)=56`. The **ties are
  meaningful**: coincident difference-speeds are exactly the *additive coincidences*
  of the speed set (its Sidon-defect). So the pair tournament reads the **additive
  combinatorics** of `{v_i}` that the runner tournament cannot see.

The observer's loneliness lives inside this object as a marked sub-family: the
observer-pairs `{0,i}` carry `w_{0i} = -v_i`, and loneliness = these are all far from
`0`. So the pair/tension tournament *contains* the original LRC on its observer-star
and *adds* the difference-set (Sidon/cocycle) structure on the rest — a strictly
richer, strongly-restricted lift.

## The honest negative: gaps are not three-gap-rigid for multi-speed

A natural complex object: vertices = the `n` **gaps** between consecutive points,
with the apex (largest gap) = loneliness (S530). One hopes for **three-gap /
Steinhaus rigidity** (≤ 3 distinct gap lengths). Computed: the distinct-gap-length
histogram runs up to `n` distinct lengths (e.g. at `n=6`, most samples have all `6`
lengths distinct), and the quantized gap-multisets are *not* a small set. **The
three-gap theorem is a single-rotation phenomenon and does not transfer to the
multi-speed (multi-rotation) setting.** So gaps are a *weaker* vertex object than
sectors (S536) — honest negative, worth recording so it is not re-attempted.

## The multitude (posed precisely, by their consistency law)

Each is "which iso-classes are exhibitable," restricted by a *different* law:

- **Harmonic / spectral:** vertices = frequencies `m ∈ {1,…,n-1}`; edges by the phase
  of the character `ĝ(m)`; the **dual of the S537 flow** picture. Restriction: the
  Fourier support / character relations.
- **Arrangement cells:** vertices = cells of the **combined braid `{x_i=x_j}` +
  threshold `{x_i=±1/n}` arrangement** on `T^{n-1}`; each cell = `(cyclic order,
  loneliness pattern)`; the orbit is a closed walk. Restriction: the orbit visits a
  thin slice (S521o). The most complex; the genuine home of LRC.
- **Incidence cells (sector × runner):** vertices = occupied `(sector, runner)`
  cells (S536 ⊗ runners); a bipartite lift. Restriction: doubly (occupancy + order).
- **Matroid flats:** vertices = flats of the **resonance matroid** `M_v` (S537);
  edges by closure/containment; LRC reads the connectivity of `M_v`. Restriction:
  representability over `Q`.
- **Time-frequency (Gabor) cells:** vertices = `(sector, harmonic)` pairs — the joint
  **space ⊗ frequency** lift unifying S536 (space) and S537 (frequency); an
  "uncertainty" tournament whose realizable cells obey a discrete **uncertainty
  principle** (a cell can't be sharp in both). The deepest unification.
- **Wiring-diagram events:** vertices = crossing events (adjacent transpositions);
  realizable = **stretchable** allowable sequences (S535 MAP-wire). Restriction:
  stretchability.

## The organizing principle

> **A complex vertex object restricts the realizable iso-classes exactly to the
> extent it carries a CONSISTENCY LAW** — cocycle (pairs/tension), occupancy
> (sectors), Fourier support (harmonics), stretchability (wiring), representability
> (matroid flats), or an uncertainty bound (Gabor cells). The pair/tension tournament
> is the cleanest new instance: the cocycle is a hard algebraic constraint, giving a
> `~10^{-9}` realizable fraction at `n=5`, and it makes the difference-set (additive)
> structure of the speeds a *first-class* part of the tournament — the structure the
> plain runner tournament is blind to.

## Verdict / next
- New star: the **tension-pair tournament** (vertices = pairs, cocycle-constrained,
  a second-order LRC on the difference set); strongly restricted; sees the Sidon
  structure via ties.
- Honest negative: gaps are *not* three-gap-rigid for multi-speed.
- Multitude posed by consistency law; the **time-frequency (Gabor) cells** and the
  **arrangement cells** are the most promising complex frontiers.
- Concrete next: (1) make the pair-tournament a genuine tournament via a Sidon-aware
  tiebreak and recompute `R`; (2) build the `(sector,harmonic)` Gabor tournament and
  test its uncertainty-type restriction; (3) the arrangement-cell walk as the exact
  LRC object.

## Artifacts
```
04-computation/lrc_complex_vertex_objects_s538.py
05-knowledge/results/lrc_complex_vertex_objects_s538.out
```
Related: S535 (mapping spectrum / metric restriction), S536 (sectors/DFT), S537
(flows/tension), S530 (apex=largest gap), S521o (arrangement/permutohedron),
S533 (channels = matroid/character).
