---
id: HYP-2032
status: SUPPORTED
source: codex-2026-06-01-S542
related:
  - HYP-2022
  - HYP-2024
  - HYP-2025
  - HYP-2026
  - HYP-2027
  - HYP-2028
  - HYP-2029
  - HYP-2031
  - HYP-2008
  - THM-382
  - THM-383
  - THM-384
---

# HYP-2032: finite Gabor trienerments turn observer tie-freeness into a marked zero column

**Claim.** The useful Gabor lift of the LRC sector model is not another global
realizability question.  It is a marked finite time-frequency trienerment:

```text
vertices = (two-sector window a, harmonic m)
G_c(a,m) = c_a zeta^(ma) + c_{a+1} zeta^(m(a+1))
```

where `c` is the fixed `n`-sector occupancy vector.  The observer danger
window is `a=n-1`, covering sectors `n-1` and `0`.  Therefore:

```text
c_0=c_{n-1}=0    iff    G_c(n-1,m)=0 for every harmonic m.
```

Thus runner-space LRC deletes all observer near-ties, while Gabor-space LRC
creates a marked zero column, a large unresolved tie block among observer
time-frequency atoms.

**Evidence.** `lrc_gabor_zero_column_s542b.py` constructs the finite Gabor
phase trienerment on `(a,m)` atoms.  Edges are directed by phase half-turn
order; ties mean zero/unresolved or same atom-angle cell.  Classes are
canonicalized under the fixed circular scaffold rather than arbitrary
`S_{n^2}` relabeling.

Global sector-vector image:

```text
n=3: raw=6   raw dihedral=2   Gabor scaffold classes=5
n=4: raw=20  raw dihedral=4   Gabor scaffold classes=18
n=5: raw=70  raw dihedral=10  Gabor scaffold classes=44
n=6: raw=252 raw dihedral=26  Gabor scaffold classes=178
```

Good-sector image:

```text
n=3: good raw=1   good Gabor classes=1
n=4: good raw=4   good Gabor classes=4
n=5: good raw=15  good Gabor classes=11
n=6: good raw=56  good Gabor classes=37
```

In every audited `n`, good vectors are exactly the vectors whose observer
Gabor window has `n` zero atoms.  Some non-good even cases have a single zero
atom from harmonic cancellation, so the full zero column is the right marked
target, not mere existence of a zero coefficient.

Bounded fixed-clock menus:

```text
n=4, B<=10: classes seen 18/18, good 4/4, open_good=108, wall_only=1
n=5, B<=8:  classes seen 44/44, good 11/11, open_good=67, wall_only=2
n=6, B<=7:  classes seen 135/178, good 22/37, open_good=20, wall_only=1
```

The AP row is wall-only in all three scans, so the open Gabor zero-column
target still needs the THM-383 compactified boundary layer.

**Trienerment duality.** The runner trienerment and the Gabor trienerment use
ties in opposite local directions:

```text
runner trienerment:
  observer lonely  = observer tie-degree 0

Gabor atom trienerment:
  observer lonely  = observer window has n zero atoms,
                     hence many unresolved/tie edges
```

The bounded samples show this numerically.  Average Gabor tie counts are
higher at observer-tie-free chambers than at observer-tied chambers:

```text
n=4: 73.3 vs 58.4
n=5: 148.4 vs 97.9
n=6: 304.3 vs 252.3
```

**Interpretation.** HYP-2028 says every raw sector-vector is globally
realizable, so this finite Gabor lift cannot prove LRC by global class
existence.  Its leverage is as a fixed-clock menu/fiber restriction:

```text
For every primitive clock, does its Gabor trienerment menu hit
the marked zero-column family, or a compactified wall analogue?
```

This is the concrete version of the upstream Gabor/trienerment proposal.  It
also matches HYP-2027's consistency-law principle: the vertex set becomes more
meaningful because the Gabor atoms preserve the sector/harmonic scaffold and
carry a zero-column law, not because arbitrary atom classes are intrinsically
important.

**Concurrent integration.** Upstream HYP-2031 probes an unanchored
`(sector,harmonic)` Gabor-cell trienerment and finds real uncertainty/tie-axis
structure but mixed good/bad fingerprints.  HYP-2032 is the observer-anchored
complement: use a two-sector window so the LRC predicate itself becomes a
marked all-zero observer column.  The lesson is consistent: unanchored Gabor
support is structure, not a certificate; observer anchoring is what turns the
Gabor angle into a target map.

**Predictions.**

1. The equivalence `good face <=> observer Gabor zero column` holds for all
   `n` for the two-sector finite Gabor window.
2. The marked zero-column Gabor class family is pure for open LRC cells under
   the fixed circular scaffold.
3. AP/regular-polygon rows require compactified Gabor wall atoms, likely by
   half-weight endpoint windows or wall-owner labels.
4. Exact cyclotomic phase comparison will preserve the observed class counts
   while removing the current floating atom-angle tolerance.
5. The next useful isomorphism group is the finite Heisenberg/symplectic
   scaffold preserving adjacency windows and harmonic sign, not full
   unmarked trienerment isomorphism.

**Files.** `04-computation/lrc_gabor_zero_column_s542b.py`;
`05-knowledge/results/lrc_gabor_zero_column_s542b.out`;
`07-reflections/lrc-gabor-zero-column-trienerment-s542b.md`.
