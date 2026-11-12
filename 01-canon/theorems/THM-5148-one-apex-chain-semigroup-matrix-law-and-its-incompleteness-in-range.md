---
id: THM-5148
title: "Block form of the O/Z/J root-change grammar, and the coverage of two of its sub-languages. (A), (B), (D) KNOWN, re-derived in block form: joining a one-apex block with apex word w at its apex acts on the rooted state (A, W, d) = (H, o + e, o - e) by M_w = [[g(P), Q+1, -[Q even]], [P+1, 2^n, 0], [[P even], 0, -[n = 0]]]; with a transitive chain of i >= 1 vertices inserted, by the 2 x 2 matrix M(p, i) = [[g(p), 2^i (R(p)+1)], [p+1, 2^(i+L(p))]]; the family obtained is the O/Z/J grammar (J = new apex, O / Z = spine letter 1 / 0). (C) FINITE-EXACT: the sub-language in which every root change is preceded by a sink (a chain inserted at every join) misses 1214 odd numbers below 2^14 at full closure. (E) FINITE-EXACT: the sub-language with at most L spine letters per root change misses 1355, 127, 49 odd numbers below 2^17 for L = 3, 4, 5; for L = 5 (63 blocks) no odd number of (10391, 2^19] is missed. The full grammar misses exactly 29 odd numbers, the largest 2647: KNOWN (HYP-9494 (D), to 10^9), re-computed here to 2^18."
status: >
  (A), (B), (D): KNOWN, re-derived in block form. The transport of the rooted triple through a block by an integer
  matrix and the semigroup law are the universal seam transfer (seven_return_native_recursions_20261029/SEAM.md
  (1)-(7); THM-4875 (B)); for this family the letter form is the O/Z/J grammar (THM-4904, Setting; Sessions
  197-200); the absorption of leading letters 0 and of trailing letters 1 is Session205's JZ = OJ
  (THM-4888 (D) 1); the join law is THM-4866 (B). New here: only the closed entries of a whole one-apex block.
  Brute force: 1470 of 1470 chains twice (author); every stage of 600 random and 3,234 exhaustive chains, 87,825
  joins, and 797,161 O/Z/J words against the matrix products (audit).
  The first item of (E) in the pre-audit draft (29 misses below 2^18) is KNOWN: Session197 C1_COMPLETE.md (to
  5,001), THM-4884 (2b) (to 50,001), HYP-9494 (D) (to 10^9, kind-pasteur-2026-11-07-S14); re-computed here to 2^18
  and by the audit to 2^16.
  (C), (E): FINITE-EXACT (exhaustive closures: every value of the sub-language below the bound is generated),
  pre-registered in three stages (commits 7479352f36, 870ad0c5d7, d67741ce21). (C) was recomputed by the audit with
  independent programs without normal forms (identical miss lists); the numbers of (E) below 2^16 and 2^17 were
  reproduced by the audit as a by-product, the statements of (E) were not audited. Not found elsewhere in the
  repository by one search each (audit for (C); author for (E): the values 10391, 130247 occur nowhere else); a
  search, not a clearance. Related finitely generated statement: HYP-9640 (B).
  Predictions: eleven; C1, C2, D1, D3, E1 CONFIRMED; C5 CONFIRMED for tuples only (the pre-audit draft's 'mean
  multiplicity 82.9' counted each chain many times: 14.3 distinct chains per odd number); C3, C4, D2, E2, E3
  FAILED; D2, D3 were registered against a fact already in canon.
  AUDIT: independent, own code, PASS WITH CORRECTIONS
  (05-knowledge/results/kps_s25_euler_field_bridge_law_20261111/audit/AUDIT.md; M1-M4 concern this file; all
  applied; MISTAKE-768).
  The Chan-Pak conjecture remains OPEN.
session: kind-pasteur-2026-11-11-S25
source: 05-knowledge/results/kps_s25_euler_field_bridge_law_20261111/README.md
scripts:
  - 05-knowledge/results/kps_s25_euler_field_bridge_law_20261111/scripts/chain_semigroup.py, direct_joins.py, finite_blocks.py, finite_blocks_trend.py
related:
  - 01-canon/theorems/THM-4904-root-change-word-fibers-carmichael-chamber-parity-law-sqrt-density.md
  - 01-canon/theorems/THM-4888-finite-oz-banks-mod3-law-involution-rigidity-residue-confinement.md
  - 01-canon/theorems/THM-4884-definability-lenses-korselt-congruence-of-spectrum-polynomial-strong-multiplicity-census-minimal-rank-sublanguages.md
  - 01-canon/theorems/THM-4875-old-mark-digit-lower-bound-lorentzian-gadget-adjoint-and-one-apex-joins-realise-every-odd-number-to-2-32-except-7-and-21.md
  - 01-canon/theorems/THM-4866-single-cycle-floor-bottom-band-uniquely-hamiltonian-one-arc-lorentzian-joins-and-crevasses.md
  - 01-canon/theorems/THM-5070-port-chain-transfer-for-odd-cycle-packing-counts-and-ternary-coverage-to-a-million.md
  - 05-knowledge/hypotheses/HYP-9494-root-change-carmichael-hierarchy-two-root-change-coverage-c1-depth-three-adic-factor.md
  - 05-knowledge/hypotheses/HYP-9640-chan-pak-as-a-matrix-semigroup-orbit-and-one-defect-digraph-spectra.md
  - 05-knowledge/hypotheses/HYP-9677-block-semigroup-conjecture-apex-path-tournaments-reach-every-large-odd-number.md
---

Notation. `H` = number of Hamiltonian paths; for a tournament `X` and a vertex `u`:
`W+_X(u) = sum over directed paths P starting at u of H(X - P)`, `W-_X(u)` the same for paths ending at `u`,
`tau(Z) = sum over all directed paths P of Z, the empty one included, of H(Z - P)`; `H` of the empty tournament is
1. `o`, `e`: the parts of `W+` with an odd, even number of vertices in `P`. `g(0) = 1`, `g(2v+1) = g(v)`,
`g(2v) = g(v) + v + (v mod 2)`; `L(p)` = bit length; `R(p)` = reversal of the complement of `p` on `L(p)` bits.
A one-apex tournament has a vertex (the apex) whose removal leaves a transitive tournament (the spine); its apex
word has letter 1 where the apex beats the spine vertex. Proofs, credits and data: part B of the source note.

## (A) Join with an inserted chain (KNOWN, re-derived)

`T = X => C => Y`, `C` transitive with `i >= 1` vertices, all arcs from earlier to later blocks except `v -> u`:

    H(T) = H(X) H(Y) + 2^i W+_X(u) W-_Y(v),        W+_T(v) = H(X) W+_Y(v) + 2^i W+_X(u) tau(Y - v).

This is THM-4866 (B) applied to `(X => C, u)`, whose rooted triple is `(H(X), 2^i W+_X(u), 0)`, plus the third row
of the seam matrix of `C => Y`.

## (B) Pure chains of one-apex blocks (block form)

For a one-apex tournament rooted at its apex with apex word of value `p` without leading letter 0: `H = g(p)`,
`W+ = p + 1`, `W- = R(p) + 1`, `tau = 2^(L(p))` (with `k` leading letters 0: `W- = 2^k (R(p) + 1)`,
`tau = 2^(k + L(p))`). With chains of `i_1, ..., i_k >= 1` vertices between the blocks,

    H = top-left entry of M(p_k, i_k) ... M(p_1, i_1) M(p_0, 1),     M(p, i) = [[ g(p), 2^i (R(p)+1) ], [ p+1, 2^(i+L(p)) ]].

In O/Z/J terms these are the words in which every `J` is immediately preceded by `O`.

## (C) The pure chains are not complete in range (FINITE-EXACT)

- Depth at most 2: 11.06 per cent of the odd numbers of `[2^17, 2^18)` missed; depth at most 3: 5.26 per cent
  (3446 of 65536), 9190 below `2^18`, with 14.3 distinct chains per odd number (671,731 of depth 2, 1,199,555 of
  depth 3).
- Full closure: 414 odd numbers missed below `2^12`, 1214 below `2^14`. First misses: 7, 13, 21, 29, 35, 41, 43,
  49, 55, 59.

## (D) Direct joins (KNOWN, re-derived): the block matrices

Let `(X, u)` have `(A, W, d) = (H(X), o + e, o - e)` and let `Y` be one-apex with apex `v` and apex word `w` of
length `n`, `P` = value of the letters 1 (first spine vertex = top bit), `Q` = value of the letters 0 (first spine
vertex = unit bit). The join `X => Y` with the arc between `u` and `v` reversed, rooted at `v`, has

    A' = g(P) A + (Q + 1) W - [Q even] d,      W' = (P + 1) A + 2^n W,      d' = [P even] A - [n = 0] d,

i.e. `(A', W', d')^T = M_w (A, W, d)^T` with `M_w = [[g(P), Q+1, -[Q even]], [P+1, 2^n, 0], [[P even], 0, -[n = 0]]]`.
The values of the O/Z/J grammar are the first coordinates of the vectors `M_(w_k)...M_(w_0)(1, 0, 0)^T`
(`(D, E, U) = (A, d, W)`; `J` = the block with the empty word, `O`, `Z` = spine letters 1, 0, the rooted edge = the
block with word `1`). After the first block `0 <= d <= A <= W` (KNOWN, Session200), and `A' >= W`.

## (E) At most `L` spine letters per root change (FINITE-EXACT)

- Words of length at most `L`, below `2^17`: `L = 3` (15 blocks): 1355 misses, largest 130247; `L = 4` (31
  blocks): 127 misses, largest 130247; `L = 5` (63 blocks): 49 misses, largest 10391.
- `L = 5` below `2^19` (57,008,668 live states): the same 49 misses; no miss in `(10391, 2^19]`.
- Below `2^16`: `L = 1`: 28,170 misses; `L = 2`: 10,078.
- For comparison (KNOWN): all words miss exactly 29 odd numbers, the largest 2647 (HYP-9494 (D), to `10^9`).
