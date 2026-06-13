# Support-Residue Calculus

**Session:** codex-2026-05-30-support-residue

The random walk through the repo kept returning to the same shape:

```text
choose supports -> project/forget -> measure the residue that survives.
```

This is visible in several threads that looked unrelated.

## The Good-Cut Thread

The good-cut count first looked like a base-path coordinate: count cuts crossed
by upward tiles. THM-349 then showed the raw bucket count is a one-dimensional
interval-cover gas. THM-350 showed quotient transport is governed by finite
half-line conservation and involution pairing.

The missing interpretation is now THM-354:

```text
goodCutCount(T, P) = n - #SCC(T)
```

for every Hamiltonian path `P`. So the coordinate costume falls away. A bad cut
is not a tiling accident; it is exactly a boundary between consecutive strong
components in the condensation order. The top bucket is strong connectivity,
the bottom bucket is transitive condensation, and the intermediate buckets are
the SCC-defect ladder.

This explains why `g` descends to `G_n/Z_2`: isomorphism preserves SCC count,
and complement reverses the condensation order without changing its components.

## The Omega Thread

The H=63 classes and the THM-025 real-root counterexample are two opposite
residue profiles of the odd-cycle support family.

For the n=8 H=63 classes:

```text
delete core vertex -> zero odd cycles
Omega = K31
H = 1 + 2*31
```

This is exact projection kill. Every odd-cycle certificate uses the same core
vertex, and the residue after deleting it is empty.

For THM-025:

```text
delete vertex 3 -> two old cycles remain
I(Omega,x) = 1 + 94x + 10x^2 + x^3
```

This is near-kill with a dangerous residue. Almost everything dies, but the
surviving cycles support just enough disjointness to create the unique
independent triple and break real-rootedness.

The useful variable is not only `H` or `alpha_1`. It is the projection residue:
which supports die, which remain, and what independence geometry the remainder
has.

## The Paley-Interval Thread

Paley versus Interval is the same story with the projection changed. Paley and
Interval can share visible support shadows while differing in multiplicity and
disjointness. Paley has more cycle multiplicity over the same support shadow;
Interval creates co-occurrence gradients and disjointness excess.

THM-143 makes this one-dimensional again: in the Interval tournament, cycle
co-occurrence through a distance `d` is linear in `d`. The slope

```text
binom(m-2, k-3)
```

is a support-residue amplifier. Longer cycles convert a small spatial gradient
into higher `alpha_j` mass, eventually overtaking Paley's multiplicity
advantage.

## The Random-Sample Threads

The deliberately noisy sample reinforced the same pattern.

- THM-053 says diagonal transfer-matrix entries are signed position residues:
  after forgetting the full path, only even-minus-odd position bias remains.
- THM-060 says the blue skeleton's bipartition is the unpaired residue of
  consecutive triples after the GS flip pairs away the rest.
- THM-094 says mod 2 kills all tournament-dependent lower Taylor terms; the
  surviving top coefficient is Redei parity.
- THM-164 says Walsh degree projection reveals a support tower: scores see the
  3-cycle layer, degree 4 sees the first non-score cycle residue, and higher
  degrees see finer packing.
- THM-145 proposes that spectral elementary symmetric functions determine the
  GLMY omega profile for circulants: topology as the residue of Fourier
  magnitude data.

These are not the same theorem, but they share an algebra:

```text
projection + cancellation + surviving support geometry = invariant.
```

## A Working Dictionary

| Project object | Supports | Projection | Residue |
|---|---|---|---|
| good-cut count | path cuts | SCC condensation | `n - #SCC` |
| H via OCF | odd cycles | conflict graph IP at x=2 | independent cycle packets |
| H=63 single-core | odd cycles | delete core vertex | empty residue |
| THM-025 | odd cycles | delete near-core vertex | two-cycle residue with alpha triple |
| Paley/Interval | cycle supports | support shadow / Fourier magnitude | multiplicity vs disjointness |
| transfer diagonal | Hamiltonian paths | vertex position parity | even-minus-odd bias |
| bucket balance | tiling half-lines | quotient bucket map | paired internal lines plus escape |
| path homology | allowed paths | boundary image projection | ghost or genuine old cycle |

## New Hypotheses

1. **Residue-first real-root filter.** Non-real-rooted `I(Omega,x)` should be
   enriched among tournaments with high deletion-loss fraction and a small
   nonempty residue whose independence profile is top-heavy. Alpha counts alone
   are too compressed.

2. **Good-cut transport is SCC transport.** Every `Delta g` transport statement
   should be rephrased as a statement about changes in the SCC condensation
   count. This should remove the path-coordinate dependence from HYP-1770 and
   HYP-1777 entirely.

3. **Interval gradients are residue amplifiers.** THM-143's linear
   co-occurrence slope should predict not just `alpha_2` excess but the first
   Walsh degree where Interval overtakes Paley. The slope profile is the
   one-dimensional residue that higher OCF weights amplify.

4. **Ghost-cycle failures are near-kill residues.** HYP-408-style failures
   should correlate with old-projection residues that are small but not empty,
   analogous to THM-025 rather than H=63.

The feedback loop is clear: when a new invariant appears coordinate-dependent,
look for the quotient whose components it actually counts; when a theorem
fails rarely, look for an almost-killed support family with a small structured
residue.

**S3 continuation:** `07-reflections/residue-feedback-loop.md` sharpens this
into HYP-1780: first obstructions may be stratified by residue rank rather
than by raw support size.

**S355 continuation:** `07-reflections/residue-calculus-theses.md` and
HYP-1785 turn the rank slogan into an implemented feature layer.  The new
`tournament_tda.py` block records SCC residue features plus max-loss
deletion-residue alpha/rank data.  The first probe separates:

- H=63 single-core exact kills: `rank_res=0`;
- THM-025 near-kill: `keep=2`, `rank_res=2`, `I(R,2)=9`;
- Paley/Interval broad fiber residues: rank 2 but not small near-kills.
