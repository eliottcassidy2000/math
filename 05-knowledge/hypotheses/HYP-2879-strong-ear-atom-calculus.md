---
id: HYP-2879
status: SUPPORTED by exact iso-class audit n<=8; open as all-n strong-ear reducibility / cofinite atom calculus
source: codex-2026-06-22-S99
tags: [tournaments, h-spectrum, strong-components, ears, insertion, forbidden-H, finite-basis, lrc14, residues, tournament-analysis]
related:
  - HYP-2877
  - HYP-2878
  - HYP-+2878
  - HYP-+2879
  - HYP-2874
  - HYP-2872
  - HYP-+2876
  - HYP-2875
  - HYP-2873
  - THM-115
  - THM-520
---

# HYP-2879: strong H-atoms grow by labelled ears with an exact cut polynomial

The strong-component atom ledger from HYP-2877 becomes recursive once the
allowed contraction is reversed.  The right growth move is not an arbitrary
even-graph edge operation; it is a one-vertex strong ear.

Let `T` be a tournament and add a new vertex `x`.  Write `sig[v]=1` if
`x -> v`, and `sig[v]=0` if `v -> x`.  Define:

- `start_T[b]`: Hamiltonian paths of `T` starting at `b`;
- `end_T[a]`: Hamiltonian paths of `T` ending at `a`;
- `Q_T[a,b]`: old-vertex permutations in which `a` is immediately followed by
  `b` and every other adjacent step is a valid edge of `T`.

Then

```text
H(T+x) =
    sum_{b: sig[b]=1} start_T[b]
  + sum_{a: sig[a]=0} end_T[a]
  + sum_{a: sig[a]=0, b: sig[b]=1} Q_T[a,b].
```

This is the exact insertion analogue of the LRC character-sum formula.  It has
linear boundary terms plus a cut/resonance matrix `Q_T`.  The new vertex is a
finite residue/cut certificate, not a scalar perturbation.

If `T` is strong, every nonconstant signature gives a strong child: `x` has an
out-neighbor in `T` and an in-neighbor from `T`; strongness of `T` then lets
`x` reach all vertices and be reached from all vertices.  Thus the strong-ear
constraint is exactly `sig` nonconstant.

## Exact audit through n=8

Script `04-computation/tournament_strong_ear_atoms_codex_s99.py` imports the
validated non-isomorphic tournament tower from the S5/S6 strong-spectrum work
and audits every nonconstant strong ear through `7 -> 8`.

- Counts match A000568 through `n=8`.
- Insertion formula failures: `0`.
- Strongness failures for nonconstant ears from strong parents: `0`.
- Strong-deletion reducibility is exact through `n=8`: every strong class on
  `n=4..8` has at least two vertices whose deletion remains strong.
- Ear-generated H spectra equal the full strong H spectra for every transition
  `3->4`, `4->5`, `5->6`, `6->7`, and `7->8`.

The largest transition is the useful stress test:

```text
7 -> 8:
  strong-ear cuts checked: 44478
  ear H-values: 297
  target strong spectrum: 297
  missing target values: []
```

## Finite basis signal

The S34 LRC rational-witness route says a finite basis of denominators can
certify many covering sets.  The strong-ear analogue is a finite basis of
parent atoms or cut weights.

For `7 -> 8`, a greedy parent-atom basis of `12` parents covers all `297`
strong `n=8` H-values.  The first few parents already cover most of the bank:

```text
parent171:H67  covers 96 new values
parent128:H131 covers 62
parent17:H27   covers 52
parent360:H131 covers 32
parent9:H33    covers 19
```

The cut-weight basis is sharper:

```text
w=3 covers 295/297 strong n=8 values;
w=1 supplies the two missing values.
```

Equivalently, balanced ears (`w=3` or, by complement, `w=4`) miss only
`H=49` and `H=75`.  Those are boundary-ear values.  This is the exact strong
atom analogue of "one denominator can resonate to zero, so use a small finite
basis."

Mac-mini S36's HYP-+2879 adds a useful disambiguation.  In this note `w` means
the ear cut weight, while HYP-+2879 uses "weight" for the number of strong
component atoms in the Moon product.  The labels should not be identified.
What transfers is sharper: the two balanced-ear misses `49` and `75` are
verified single irreducible strong atoms on `n=7`, so the boundary cut
`w=1` branch is exposing a single-atom obstruction rather than a product
branch.

## Proof program

The all-n version should split into three lemmas.

1. **Strong-ear reducibility:** every strong tournament on `n>=4` has a
   vertex whose deletion remains strong, ideally at least two such vertices.
   The audit verifies this through `n=8`; it should be a Moon/ear-decomposition
   theorem for strong tournaments.
2. **Ear cut formula:** the displayed `start/end/Q` formula is elementary and
   formalizable by partitioning Hamiltonian paths according to the position of
   the new vertex.
3. **Cofinite atom generation:** use strong-min lower bounds for the small
   forbidden boundary and use balanced-ear cut polynomials to show direct
   strong-core re-entry.  The observed `w=3` near-surjectivity at `7->8`
   suggests balanced ears are the generic/high-main-term branch; exceptional
   low values route to boundary ears.

This reframes `{7,21}` as absent solutions of recursive cut polynomials before
the Busch/Moon strong-min boundary lets the spectrum re-enter.  They are not
minor-closed shadows of the even-graph quotient.

## LRC / residue transfer

S34's rational witness count has the form

```text
N(S,D) = main term + resonance error.
```

HYP-2879 gives the tournament analogue:

```text
H(T+x) = boundary main terms + Q-cut resonance.
```

The productive proof state is therefore the labelled boundary carrier:
denominator basis plus resonances on the LRC side, and `(start,end,Q)` plus cut
signature on the tournament side.  This also explains why runner deletion is
the wrong minor operation for LRC: the preserved object lives on residues/cuts,
not on raw subsets.

S35's HYP-+2878 makes this bridge sharper.  Its LRC atom is a residue-covering
component: a single strongly resonant obstruction can force `N(S,D)=0` at a
prime denominator.  This hypothesis is the tournament-side analogue: a single
strong parent plus one labelled ear can force or miss an `H` value through the
`Q` matrix.  The shared obstruction is not the raw vertex set.  It is a finite
boundary carrier plus a resonance ledger.  In the corrected S35 route,
covering at one prime is a structured interval event of order `~10%`, not a
Poisson-tiny event; the proof pressure comes from CRT over-determination across
a many-prime finite basis.  The tournament analogue is that one cut weight can
miss by resonance, while several cut weights or parent atoms break the aligned
miss.

KPS S31e's HYP-2878 gives the concrete odd-cycle target for this transfer.
The even-graph metagraph `E_7` has exactly `1496` chordless `C_5` holes and
`196 = 14^2` chordless `C_7` holes.  Since this ear calculus routes `{7,21}` to
cut-polynomial failures rather than even-graph minors, the next single-object
test is not minor containment.  It is whether the exposed-slot graph encoded by
`Q_T` has canonical odd cycles corresponding to the E7 `C_5` pentagon
obstruction and the apex-prime `C_7` heptagon family.  The S31e enrichment
narrows the target further: the E7 heptagons are concentrated on `34/54`
classes, complement-paired, and centered on `9`-edge classes rather than the
naive `7`-edge cycle.  Thus the `Q_T` test should look for complement-symmetric
central exposed-slot classes, not merely any chordless heptagon.

## Artifacts

- `04-computation/tournament_strong_ear_atoms_codex_s99.py`
- `05-knowledge/results/tournament_strong_ear_atoms_codex_s99.out`
- `07-reflections/strong-ear-atoms-and-finite-bases-codex-s99.md`
