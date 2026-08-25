---
id: HYP-2879
status: >
  PARTIALLY PROVED by THM-4097. The exact insertion formula and all-order
  strong-ear reducibility are proved (the latter from Moon 1966, Theorem 2);
  the cofinite atom calculus, all-order finite cut-basis law, and interval
  tiling remain OPEN. The older iso-class audit through n<=8 remains an exact
  finite control.
source: codex-2026-06-22-S99
tags: [tournaments, h-spectrum, strong-components, ears, insertion, forbidden-H, finite-basis, lrc14, residues, tournament-analysis]
related:
  - HYP-2877
  - HYP-2878
  - HYP-2881
  - HYP-+2878
  - HYP-+2879
  - HYP-2874
  - HYP-2872
  - HYP-+2876
  - HYP-2875
  - HYP-2873
  - THM-115
  - THM-520
  - THM-4097
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

THM-4097 resolves the first two parts of the earlier three-lemma program and
sharpens the remaining carrier.

1. **Strong-ear reducibility -- PROVED from a cited theorem.** Moon 1966,
   Theorem 2 gives `s(n,n-1)=2`: every strong tournament of order at least four
   has at least two induced strong parents of order one less. Thus every strong
   tournament is recursively a nonconstant one-vertex ear from `C_3`.
2. **Ear cut formula -- PROVED.** THM-4097 proves the displayed
   `start/end/Q` formula, compiles `Q` by subset path DP, and refines the
   response to an integral symmetric cut weight `w` plus zero-sum orientation
   field `h`. The exact `623/735` complement pair proves that the unoriented
   `(w,{S,V\S})` quotient is lossy.
3. **Cofinite atom generation -- OPEN:** use strong-min lower bounds for the small
   forbidden boundary and use balanced-ear cut polynomials to show direct
   strong-core re-entry.  The observed `w=3` near-surjectivity at `7->8`
   is followed at `8->9` by the exact two-weight basis `{3,4}`: weight four
   covers `1478/1482` values and misses `89,93,105,125`, all supplied by weight
   three. This remains finite-order evidence, not an all-order theorem.

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
`196 = 14^2` chordless `C_7` holes.  S37/HYP-+2880 verifies the support-level
identity `directed C5 = H=7 K3` in cycle space.  HYP-2881 adds the quotient
guardrail: under the fixed-path path-fundamental cycle-space map, the THM-200
directed pentagon support is a single E7 class, while an E7 C5 metagraph hole
is a five-class quotient cycle.  The useful target for the exposed-slot graph
encoded by `Q_T` is therefore incidence profile: do the `Q_T` odd cycles mark
the same obstruction layer, especially the
`k3_forces_pentagon` classes that hit many E7 C5 holes without defining them?

The S31e enrichment still narrows the apex-prime side: the E7 heptagons are
concentrated on `34/54` classes, complement-paired, and centered on `9`-edge
classes rather than the naive `7`-edge cycle.  Thus the `Q_T` test should look
for complement-symmetric central exposed-slot classes and C5/C7 incidence
profiles, not merely any chordless heptagon.

## Artifacts

- `04-computation/tournament_strong_ear_atoms_codex_s99.py`
- `05-knowledge/results/tournament_strong_ear_atoms_codex_s99.out`
- `07-reflections/strong-ear-atoms-and-finite-bases-codex-s99.md`
- `01-canon/theorems/THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension.md`
- `04-computation/tournament_order9_strong_ear_spectrum_thm4097.py`
- `04-computation/order_nine_strong_ear_cut_field_thm4097_independent_audit.py`
