# HYP-2648: LRC14 State-Word Invariants

**Status:** OPEN / structural invariant scout  
**Date:** 2026-06-19  
**Session:** codex-2026-06-19-S31  
**Tangent:** T895  
**Script:** `04-computation/lrc14_state_word_invariants_codex_s31.py`  
**Result:** `05-knowledge/results/lrc14_state_word_invariants_codex_s31.out`

## Claim

For a row `E`, the structural object behind the current LRC14 invariants is the
measured cyclic word

```text
W(E) = {(I, |I|, M_E(I)) : I is a wall atom, M_E(I) subset {1,...,6}}
```

where `M_E(I)` is the set of inner seventh-sectors missed by `{e*x : e in E}`
on the atom `I`.  The important scalar invariants are valuations or quotients of
this word, and AP-to-defect comparisons are couplings of two such words on a
common wall refinement.

This is the next layer after HYP-2647: signed wall transport is not only an
addressed matrix; it is a transport between measured state words.

## What It Preserves And Destroys

The state word preserves the exact LRC predicate seen by the sector-cover model:
which inner sectors are hit or missed on every wall atom, with exact atom mass.
It also preserves cyclic adjacency of state changes.

It destroys runner identities inside a state unless they are added back as an
address.  That is why useful refinements attach fold target, moving endpoint
sector, relation/coset class, or proof-obligation labels to the atoms before
scalarizing.

## Quotient Ladder

The same finite object gives the existing invariants by projection:

- `p_t = meas{|M_E|=t}` and `L_y(E)=sum_t p_t g_k(t)` are the missed-count histogram.
- HYP-2644's far-element plateau is the coarser quotient `p0 + p1/7`.
- HYP-2643's fold target profile is an additive address attached to selected state transitions.
- HYP-2647's signed wall transport is a valuation of the common-refinement coupling `W(A,D)`.
- HYP-2646's mod-7 signed reciprocal factorization is the Fourier/coimage dual of retaining finite address before taking signed sums.
- HYP-2640's relation rank is capacity/preimage data, not the coefficient; the coefficient is state address plus signed phase.

## Computed Evidence

The S31 scout verifies that the state-word coupling reproduces the known HYP-2642/HYP-2647 number for
`AP9=(0,1,2,3,4,5,6,7,8)` and `D=(0,1,2,3,4,5,6,7,9)`:

```text
signed D-AP = -887/158760
positive    = 9749/158760
negative    = 2659/39690
neutral     = 4393/5880
support     = 76 state-pairs
```

The neutral mass is about three quarters of the interval.  The signed drop is
therefore not a bulk-density effect; it is the residue left after a large
addressed cancellation.

Single-row signatures show why scalar features are too lossy:

```text
AP9:              32 states, H=3.73574, p0+p1/7=1229/2744, folds=12
D top+1:          30 states, H=3.71220, p0+p1/7=773/1764, folds=12
M9 spectrum 2k:   43 states, H=3.64936, p0+p1/7=59597/123480, folds=12
sporadic 1347:    43 states, H=4.78764, p0+p1/7=81929/452760
KPS pocket:       50 states, H=4.69887, p0+p1/7=1446451/8376060
```

AP9, the near-AP defect, and the doubled-top spectrum row all have the same
visible fold count, but their state support, transition word, single-sector
bias, and transport couplings differ.  The sporadic and KPS pockets have much
larger state support/transition complexity; they are not AP-like after the
state-word quotient even when additive relation ledgers are rich.

## Proof Target

Define a finite addressed state category:

```text
state atom
  -> missed-set label
  -> missed count
  -> cyclic transition type
  -> fold target
  -> moving-sector transport address
  -> mod-7 coimage / reciprocal phase
  -> signed valuation
```

Then prove the high-`L_y` rows are forced into a small set of state-word
templates.  The k=9 endpoint defect should be the maximal non-AP template inside
the near-AP transport class, and all other templates should route to one of:

1. Freiman small-excess finite certificate (HYP-2638).
2. Relation-covered but non-AP GAP slack (HYP-2639).
3. Signed coimage/reciprocal cancellation tail (HYP-2646/HYP-2636).
4. Far-element plateau contraction (HYP-2644).

## Tournament Analysis

Vertices are proof quotients, not runners or arcs:

```text
measured state word
signed state transport
missed-set histogram
p0/p1 plateau
fold-target profile
mod-7 reciprocal phase
raw relation rank
raw runners
```

The observable is deterministic information retention for the LRC bad-set
measure.  Orient `A -> B` when `B` is a quotient or signed valuation of `A`.
This gives the Hamiltonian path

```text
state word -> transport -> histogram -> plateau -> fold profile
-> mod-7 phase -> relation rank -> runners
```

Challenged assumption: the useful tournament vertices do not have to be runners,
arcs, or residue classes.  Here they are finite addresses and proof obligations
that remember exactly the amount of LRC predicate information needed before the
final scalar is taken.

## Next Work

1. Add fold target and moving-sector labels directly to the S31 state-pair table.
2. Search for state-word template collisions: rows with identical histogram and
   fold profile but different `L_y` or transport sign.
3. Prove a state-word compression lemma for high `L_y`: low entropy plus low
   transition complexity should imply near-AP or Freiman small-excess.
4. Translate the state-word coupling into the HYP-2646 reciprocal language by
   Fourier transforming state indicators instead of raw runner indicators.
