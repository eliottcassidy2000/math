---
id: THM-956
title: THE BLOCK-STRUCTURE REDUCTION — residual families whose sorted speeds split into blocks of size ≤ 6 with junction ratios ≥ G₀(k_src, k_tgt) (exact proven table, worst case 2030 at 6→6) are LONELY by THM-955 + window_tail_glue with NO search; the entry-length law l_k = max(2, 4k/(7−k)) (the k ≥ 4 naive-entry failure caught and fixed); hence THE DENSE CORE IS EXACTLY the single ≥7-comparable-block families, where the 7-wall pair-overlap crumb is the one missing elementary ingredient; plus the floor-transfer lemma (width δ ⟹ ≥ ⌈δD⌉ − 1 grid points: item 3's continuous→discrete bridge)
status: constants derived exactly (G₀ table in the out file); end-to-end verification 60/60 towers with exact rational witnesses at proven junctions; floor transfer verified exact
source: opus-2026-07-17-S339 (owner: work the dense core, the 7-wall pair-floor, and the a-priori floors)
depends_on: [THM-955 (cluster gap widths), LRCLacunaryNest.window_tail_glue/norm_glue (the formal glue), THM-936 (accessible taxonomy), the formalization picture items 1-3]
scripts: 04-computation/block_reduction_G0_opus_S339.py -> 05-knowledge/results/block_reduction_G0_opus_S339.out
---

# THM-956 — the block-structure reduction

**Theorem.** Sort the 13 speeds; split at every junction with ratio at
least the G₀ table value (source-block size to target-block size; exact
rationals: 1→1 needs 3.1, 6→6 needs 2030; worst case over block shapes at
the residual's own compression cap 13). If every resulting block has ≤ 6
speeds, the family is LONELY: THM-955 certifies a safe window inside each
block's running window — with the ENTRY-LENGTH LAW l_k = max(2, 4k/(7−k))
(the naive l = 2 glue entry makes the width bound NEGATIVE for k ≥ 4; the
glue must demand size-dependent entries) — and window_tail_glue jumps each
junction. No search, no within-block structure assumptions beyond the
ratio cap. Verified end-to-end: 60/60 random towers, exact witnesses.

**The dense core, exactly.** The residual minus this theorem = families
with a block of ≥ 7 comparable speeds. At k = 7 the union bound saturates
(7·(1/7) = 1): loneliness needs a quantified PAIR-OVERLAP crumb inside the
window — the 7-wall pair-floor, now provably THE unique remaining
ingredient of the elementary route. Open items 1 and 2 are one item.

**The floor transfer (item 3's bridge).** A safe width δ contains
≥ ⌈δD⌉ − 1 points of every 1/D-grid (exact; verified): THM-955's
continuous width floors and the census funnel's discrete liveCount floors
are the same pigeonhole seen from the two sides.

**Lean note.** The composition is already ~90% formal: window_tail_glue
and norm_glue are sorry-free; the cluster-gap draft supplies the widths;
the G₀ arithmetic is norm_num-checkable rationals.
