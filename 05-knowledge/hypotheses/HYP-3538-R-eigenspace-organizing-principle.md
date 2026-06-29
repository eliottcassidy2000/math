---
id: HYP-3538
title: The even/Brouwer (+) odd/Borsuk-Ulam two-index split IS the +-1 eigenspace split of the reversal R -- TESTED -- the cap's pairwise co-emptiness matrix M is R-symmetric with the PERRON (bulk, SOS-provable) mode in the R-EVEN eigenspace and a nonzero R-ODD eigenspace carrying the signed (negatively-contributing) obstruction; the witness index f lives on the R-fixed half-system (THM-583). The organizing principle is SUPPORTED as a real structural feature of the cap.
status: SUPPORTED (Task 2 verified for the binding configs); the precise claim "R-odd IS the binding deviation" needs the S75e-SOS gap localization. Not a proof of LRC.
source: mac-mini-2026-06-29-S8
related:
  - THM-583   # f (witness) via the half-system -- the R-odd/fixed witness datum
  - THM-582   # the two-index synthesis (Redei odd / lonely even)
  - THM-581   # floor is even-category (existence)
  - HYP-3085  # kps reflection-Perron co-emptiness M (the Perron mode)
  - HYP-3242  # cap = Euler char of the cover nerve
  - HYP-3244  # half-tiling lift->quotient (must not drop the eps=-1 part)
script: 04-computation/lrc_cap_R_eigenspace_obstruction_macmini_20260629.py
result: 05-knowledge/results/lrc_cap_R_eigenspace_obstruction_macmini_20260629.out
---

# HYP-3538 -- the two-index split IS the R +-1 eigenspace split (tested)

## The organizing principle (user)
If `eps` is the obstruction class, the LRC "signed certificate" = the R-ODD (`eps=-1`) eigenspace of
the cap; the R-symmetric half data on `(0,1/2)` is the SOS/Brouwer-provable piece.  The even/Brouwer
(+) odd/Borsuk-Ulam split (mac-mini S78 / THM-582) is then literally the `+-1` eigenspace split of the
reversal `R`.  The half-tiling lift->quotient (HYP-3244) drops the `eps=-1` eigenspace -- the
coordinate the half-compression must NOT discard is the R-odd part.

## TASK 2 test (binding configs): SUPPORTED
The pairwise sector co-emptiness matrix `M` (6x6, inner sectors 1..6) is `R`-symmetric, where `R` is
the time reflection `t->-t` acting on sectors as `(1 5)(2 4)` fixing `3,6`.  R-even eigenspace dim 4,
R-odd dim 2 = `span(e1-e5, e2-e4)`.  Decomposing `M = M_even (+) M_odd`:

| config | cap(cover) | Perron eig | Perron in | R-odd block eigs | tr(M_odd) | R-odd -> S2 |
|--------|-----------|-----------|-----------|------------------|-----------|-------------|
| consec_8 | 0.327 | 0.867 | **R-EVEN** | 0.172, 0.096 | 0.268 | -0.134 |
| consec_9 | 0.416 | 0.747 | **R-EVEN** | 0.144, 0.088 | 0.232 | -0.116 |
| {1,5,7,8,9} | 0 | 1.412 | **R-EVEN** | 0.291, 0.245 | 0.536 | -0.268 |
| {1,11,12,13} | 0 | 1.842 | **R-EVEN** | 0.277, 0.245 | 0.521 | -0.261 |

Findings: (1) the **Perron/bulk mode is R-EVEN in every binding config** -- the SOS-provable bulk is
the R-symmetric part, as the principle predicts.  (2) `M` has a **nonzero R-ODD eigenspace** (2-dim)
that contributes NEGATIVELY to `S2` (`-tr(M_odd)/2`) -- the signed residual sitting beyond the even
bulk.  (3) The cap's inclusion-exclusion splits into R-fixed-subset (positive) + R-paired-subset
(negative) parts.  So the `+-1` eigenspace split is a REAL structural feature of the cap, with the
R-even = provable bulk and R-odd = the signed obstruction.

## Refinement of the two-index picture
This sharpens THM-581/THM-582.  The split is per-OBLIGATION:
- **Floor / existence** (lonely MEASURE, sigma-invariant as a function of `t`): purely R-EVEN
  (THM-581) -- a lonely tournament has a source, not self-converse, so no R-odd witness.
- **Cap / concentration** (the co-emptiness MATRIX `M`, the gK8 bound): genuinely carries an R-ODD
  eigenspace (this test) -- the signed obstruction the SOS even-bulk misses.
So the user's principle holds on the CAP side; the floor side is already pure-even.  Both are the
`eps=-1` content of `R`: on the cap it is `M_odd`; on the witness it is `f` (THM-583, the half-system
count).

## Open / next
Confirm "R-odd IS the binding deviation" by localizing the S75e Fejer-Bochner (cosine = R-even) SOS
gap: show the SOS minorant captures `M_even` and the residual is exactly `M_odd`.  Then the LRC cap
proof = [R-even part by S75e SOS] + [R-odd part by the Borsuk-Ulam/odd-degree certificate], the
literal `+-1` eigenspace split.  And ensure any half-compression retains `M_odd` (HYP-3244 warning).
