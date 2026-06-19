# HYP-2650: LRC Invariant Separation

**Status:** OPEN / structural invariant separation scout  
**Date:** 2026-06-19  
**Session:** codex-2026-06-19-S32  
**Tangent:** T897  
**Script:** `04-computation/lrc_invariant_separation_codex_s32.py`  
**Result:** `05-knowledge/results/lrc_invariant_separation_codex_s32.out`

## Claim

The LRC structures currently visible in the repo are determined only after a
finite address is retained on the wall or crossing arrangement.  Scalar
summaries are useful routing labels, but they are not determinants.

Two concrete forms of the same principle appeared in the S32 scout:

```text
max-min M(S):       optimal clearance word = (q, folded residues, active runners, crossing source)
LRC14 sector model: measured state word   = wall atom mass + missed seventh-sector set
```

This complements HYP-2648.  HYP-2648 identifies the measured cyclic state word
as the LRC14 carrier; HYP-2650 asks which weaker invariants are already
separating, and which ones provably mix target fibers.

## Exact Max-Min Separation

The script enumerated `1743` primitive speed sets from `(k,B)=(5,13),(6,11)`
and grouped them by candidate invariants.  Coarse speed summaries all have mixed
fibers for the exact value `M(S)`:

```text
sumset_excess: mixed M fibers 9, largest mixed fiber 548
fold_count:    mixed M fibers 7, largest mixed fiber 603
fold_profile:  mixed M fibers 113, largest mixed fiber 228
gap_pattern:   mixed M fibers 441, largest mixed fiber 9
```

Even the optimal denominator `q` and `(q, active-count)` mix values.  The first
tested finite word that separated `M` in this bank was

```text
(q, j, active runner values)
```

and the slightly coarser clearance word `(q, sorted folded residues)` also had
no mixed fibers.  This is not claimed as a theorem, but it gives a clear target:
for exact max-min problems, the determinant is the local crossing word on the
THM-524 envelope, not the raw speed row.

AP-floor examples show why.  `S=(1,2,3,4,5)` has `M=1/6`, `q=6`, `j=1`,
folds `(1,1,2,2,3)`, active runners `(1,5)`, and sources
`peak:3, sum:1+5, sum:2+4`.  Nearby rows can share additive labels while
changing the active crossing source and exact `M`.

## LRC14 Sector Separation

The script also enumerated all `1287` primitive `k=9` offset rows
`E={0}+8-subsets of [1,13]` and compared summaries against `L_y`.  Again,
coarse summaries mix target fibers:

```text
sumset_excess:        mixed L_y fibers 11, largest mixed fiber 278
fold_count:           mixed L_y fibers 11, largest mixed fiber 409
fold_profile:         mixed L_y fibers 130, largest mixed fiber 6
transition_signature: mixed L_y fibers 54, largest mixed fiber 3
```

But the missed-count histogram, state-mass table, and full state word had no
mixed `L_y` fibers in this bank.  The histogram separates `L_y` because `L_y` is
literally a valuation of missed-count masses.  HYP-2648's full measured state
word then carries the richer structure that the scalar `L_y` discards:
transition complexity, sector bias, signed transport, fold targets, and
coimage phase.

The top rows illustrate the quotient ladder:

```text
AP9              L_y=26083/52920, states=32, entropy=3.73574, transitions=((1,68),(2,8),(3,2))
near-AP defect   L_y=38681/79380, states=30, entropy=3.71220, transitions=((1,77),(2,8),(3,1))
```

Both are low-complexity measured words, but their exact state masses and
transition words differ.  Larger state support and higher transition entropy
are candidate certificates for routing rows away from the dangerous near-AP
template.

## Proof Target

The common object is an addressed wall/crossing sheaf:

```text
atom mass
  -> wall/crossing address
  -> missed-sector or clearance state
  -> active runner / fold target / coimage phase
  -> scalar valuation
```

For LRC14, a plausible proof skeleton is:

1. classify high-`L_y` measured state-word templates;
2. prove the AP and one-step near-AP template are the only low-complexity
   contenders near the cap;
3. route every other row by a retained address: Freiman small-excess,
   relation-covered GAP slack, signed coimage cancellation, or far-element
   plateau contraction;
4. only then take scalar valuations.

The negative form is equally useful: any candidate invariant with mixed target
fibers cannot be the determinant unless an address is reattached.

## Tournament Analysis

The tournament vertices are proof carriers rather than runners:

```text
addressed_wall_sheaf
optimal_clearance_word
measured_state_word
missed_histogram
additive_freiman_labels
fold_profile
raw_speeds
```

The observable is information loss under quotienting.  Orient `A -> B` when
`B` is a quotient of `A` that can mix target fibers after scalarization.  The
computed Hamiltonian path is:

```text
addressed wall sheaf
> optimal clearance word
> measured state word
> missed histogram
> additive Freiman labels
> fold profile
> raw speeds
```

Alternate vertex sets considered: runners, gaps, fixed seventh-sectors, wall
atoms, section boundaries, pair-crossing events, residue/coimage classes,
Fourier modes, matroid circuits, and proof obligations.  The quotient that best
preserves the LRC predicate in this scout uses wall atoms and addressed
crossing/state labels.  It destroys raw runner identity unless active-runner,
fold-target, or coimage labels are attached.

Challenged assumption: the determinant of LRC structure is not a tournament on
runners or arcs.  The useful tournament is a comparison of quotient maps from a
finite addressed wall sheaf to scalar invariants.

## Next Work

1. Prove or refute clearance-word separation beyond the small exact bank.
2. Search for rows with identical missed-count histogram but different signed
   transport, fold target, or coimage phase.
3. Define a canonical addressed sheaf object that simultaneously contains the
   THM-524 crossing source and the HYP-2648 state word.
4. Use the sheaf quotient tournament as an invariant-ranking harness in future
   LRC14 proof scripts.
