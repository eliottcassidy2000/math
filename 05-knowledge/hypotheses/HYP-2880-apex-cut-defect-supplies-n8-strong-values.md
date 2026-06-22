---
id: HYP-2880
status: SUPPORTED by exact rooted n=8 apex-extension audit; open as a constructive strong-atom operation
source: codex-2026-06-22-S99
tags: [tournaments, tiling-model, apex-tile, h-spectrum, strong-atoms, n8, balanced-cuts, tournament-analysis]
related:
  - HYP-2879
  - HYP-2877
  - HYP-2878
  - HYP-2872
  - HYP-2008
  - THM-333
  - THM-520
results:
  - 05-knowledge/results/tournament_apex_cut_weight_ledger_codex_s99.out
---

# HYP-2880: apex-cut defect supplies n=8 strong values

Adding the fixed-path apex vertex gives a useful constructive operation for
strong `H` values, but the proof carrier is not raw cut balance.  It is the
insertion-weight profile of the apex over old Hamiltonian paths.

Start with a rooted `n=7` tournament on vertices `0..6` with fixed path
`0 -> 1 -> ... -> 6`.  Add vertex `7` at the end of the fixed path, so `6 -> 7`
is fixed.  The six free apex-row bits decide which vertices in `{0,...,5}` are
beaten by `7`.  Let

```text
w = #{i in {0,...,5} : 7 -> i}.
```

For every Hamiltonian path `P` of the old tournament, the apex contributes the
number of legal insertion slots:

```text
before P[0]      if 7 -> P[0],
after P[-1]      if P[-1] -> 7,
between a,b      if a -> 7 -> b.
```

Thus

```text
H(T + apex) = sum_{old Hamiltonian paths P} insertion_weight_apex(P).
```

This formula is exact and is the right rooted refinement of the user's
balanced-cut lead.

## Exact audit

The S99 script audits all `2^15 = 32768` rooted `n=7` bases and all `64` apex
cuts, hence all `2^21 = 2097152` rooted `n=8` rows represented as apex
extensions.

For strong `n=8` extensions by raw apex weight:

```text
w=1: strong_rows=191394, H_count=269, min=2, max=351, has49=True, has75=True
w=3: strong_rows=645112, H_count=463, min=2, max=513, has49=True, has75=True
```

Balanced raw `w=3` cuts are indeed broad suppliers in the rooted ledger: in
`[1,100]` they miss only

```text
1,5,7,20,21.
```

But raw rooted `w` does not isolate the `49/75` clue: `w=1` also supplies both
values, and raw rooted `w=3` can supply both as well.  This does not contradict
HYP-2879.  HYP-2879 works at the non-isomorphic strong-ear target-spectrum
level, where balanced cut weight `w=3` covers `295/297` strong `n=8` values and
misses exactly `49,75`, while adding boundary weight `w=1` covers all.  The
rooted fixed-path ledger here deliberately sits before that quotient; it shows
that raw cut cardinality washes out the exceptional mechanism.  The sharper
structure is visible only in the insertion distributions.

## The single-unbalanced apex defect

The single-unbalanced witnesses are nearly constant-weight over old
Hamiltonian paths:

```text
H8=49: base_H=25, insertion distribution {1:1, 2:24}
       so 49 = 2*25 - 1.

H8=75: base_H=39, insertion distribution {1:3, 2:36}
       so 75 = 2*39 - 3.
```

Both use strong `n=7` bases.  This is the plausible apex-tile signal.  A
single unbalanced apex row almost doubles the old strong atom and then
subtracts a small odd source-sink defect.  That is exactly what the S530 apex
tile suggests: the source-sink tile is a boundary term that closes the path
into a cycle, so it should appear as a one-row defect rather than as an
ordinary balanced interior cut.

The post-pull KPS S31f scripts sharpen the caveat.  In the fixed-base tiling
model, the literal longest tile `(n-1,0)` is not required for `H=49` or `H=75`
at `n=7`: `apex_tile_H49_75_kps.py` finds `apex_NOT` witnesses for both.  So
the word "apex" in this hypothesis means the boundary row / insertion-defect
operation, not the literal longest tile.  The same pull strengthens the atom
side: `strong_atom_7adic_kps.py` confirms that `49=7^2` and `75` are
irreducible `m=7` strong atoms, while `7` and `21` remain forbidden and
7-multiple atoms first appear at the apex index `m=7`.

Balanced `w=3` witnesses have a different profile.  Examples for the same
values use distributions such as

```text
H8=49: {2:4, 3:11, 4:2}
H8=75: {2:2, 3:13, 4:8}
```

These are centered around the balanced weight `3`; they are useful suppliers,
but less diagnostic of the apex boundary mechanism.

## Corrected proof target

The constructive strong-atom operation should be stated as:

```text
old strong atom H
  -> apex extension with near-constant insertion weight c
  -> new strong value c*H - small apex defect.
```

For the `w=1` witnesses above, `c=2` and the defects are `1` and `3`.

This gives a concrete way to search for cofinite strong-atom realization after
HYP-2877 and HYP-2879.  Instead of asking whether arbitrary `n=8` rows realize
a missing value, ask whether a labelled old atom admits an apex cut whose
insertion profile has controlled defect and survives the strong-ear quotient.
That is a preservation-labelled operation, not a lossy even-graph contraction.

## Assumption challenge

Candidate vertices considered: runners, free arcs, raw cut size, balanced cuts,
unbalanced cuts, insertion slots, strong components, old `H` atoms, apex tile
orientation, even-graph holes, and proof obligations.

Chosen quotient: apex insertion-weight profile.  It preserves `H` exactly under
apex extension and destroys most old-row geometry.  Raw balanced/unbalanced cut
size is too coarse: it sees that rooted `w=3` supplies many values, but it
misses the quotient-level fact from HYP-2879 that balanced strong ears still
leave exactly the `49,75` boundary holes.  The near-doubling minus odd-defect
law is the information retained by the single `w=1` apex row.

## Next tests

1. Classify all near-constant apex insertion profiles `c*H-d` for `n=7 -> n=8`.
2. Determine whether every large odd missing strong value has a predecessor
   atom and an apex defect `d` in a finite allowed set.
3. Compare the odd defects with E_7 odd holes from HYP-2878: the `d=1,3`
   defects may be the low odd-cycle boundary, while the balanced `w=3` layer is
   the interior supplier.
4. Formalize the insertion formula as the constructive counterpart to
   `H(T)=prod H(C_i)` from HYP-2877.
