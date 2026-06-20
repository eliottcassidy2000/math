# HYP-2675: LRC(14) wide-branch sector ridge

**Status:** supported computationally / proof still open (codex-s47,
2026-06-20).

The KPS comfortable-margin route reframes the remaining sector proof obligation
as a direct wide-row bound:

> for primitive `E` with `0 in E`, `|E|=k=8..12`, and `span(E)>14`, prove
> `p0(E)=meas(S7(E)) <= cap_k`.

This hypothesis reserves the next scout around that obligation.  The intended
test is not a standalone `Delta_w` constant.  It records exact `p0`, missed
sector distributions, additive-energy/Freiman fingerprints, squarefree divisor
profiles, shell-1 occupancy, state-word entropy, and a row-risk Tournament
Analysis for wide rows.

## Exact findings

Stored output:
`05-knowledge/results/lrc14_wide_branch_ridge_codex_s47.out`.

The first correction is terminological: the sector condition `span(E)>14` is
not one proof branch.  It splits into:

- **boundary one-far rows:** `span(E)>14` but second-largest `<=14`;
- **true-wide-base rows:** second-largest `>14`.

In the exact finite boxes:

- `k=9`, `{0}+8` choices from `[1,20]`: `125970` raw rows, `122922`
  primitive rows with `span>14`.  The all-`span>14` leader is boundary,
  `E=(0,2,4,6,8,10,12,14,15)`, with
  `p0=437/1176` and cap margin `20627/168168`.  The true-wide leader is
  `E=(0,4,6,8,10,12,14,15,16)`, with `p0=321/980` and margin
  `11681/70070`.
- `k=8`, bound `18`: boundary leader
  `(0,3,5,7,9,11,13,15)` has margin `130133/840840`; true-wide leader
  `(0,3,6,9,12,14,15,18)` has margin `2819/17640`.
- `k=10`, bound `16`: boundary leader
  `(0,2,4,6,7,8,10,12,14,16)` has margin `1303/9555`; true-wide leader
  `(0,2,4,6,8,10,12,14,15,16)` has margin `1175/7644`.

Named KPS/codex rows match the comfortable-margin reading:

- KPS spread-base row `(0,3,5,7,9,11,13,14,15)` has
  `p0=56593/210210`, margin `94609/420420`.
- HYP-2671 dyadic row `(0,1,2,4,8,12,16,20,24)` has `p0=29/112`,
  margin `3769/16016`.
- The best structured three-cluster probe
  `(0,1,2,10,11,12,20,21,22)` has `p0=2813/10780`, margin `8174/35035`.

## Claim being tested

The only rows that can approach the direct cap retain a low-rank
near-consecutive/GAP scaffold.  Once the state-word support spreads across
genuinely wide additive structure, `p0` drops with comfortable margin.  The old
large-`Delta_w` resonances are compensated by small plateau.

The proof target is now two-stage:

1. a **boundary collar lemma** for `second-largest <=14`, where the leaders are
   AP-like one-far rows;
2. a **true-wide sector-cover deficit lemma** for `second-largest >14`, likely
   expressed in Freiman/additive-energy and measured state-word terms.

## Tournament Analysis

The scout uses row/proof-obligation vertices rather than runner or arc vertices.
The pairwise observable is the exact risk ratio `p0/cap_k`.  The stored
tournament is transitive on the top twelve vertices (`0` directed 3-cycles);
the Hamiltonian path starts with the bounded/near-consecutive rows, not the
true-wide rows.  This preserves the LRC predicate `p0(E)<=cap_k` while
destroying runner-label data after the S3 sector quotient.

## Artifacts

- `04-computation/lrc14_wide_branch_ridge_codex_s47.py`
- `05-knowledge/results/lrc14_wide_branch_ridge_codex_s47.out`

**Open:** no theorem is claimed yet.  The missing proof is still the analytic or
finite-to-infinite statement turning these exact ridge fingerprints into
`span>14 => p0(E)<=cap_k`.
