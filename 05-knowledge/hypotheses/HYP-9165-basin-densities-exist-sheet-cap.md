---
id: HYP-9165
title: "Basin densities exist on both sheets of 3x+1, with dens B({1}) = 0.3269, dens B({5,7,10}) = 0.3248, dens B({17,...,91}) = 0.3484 for 3x-1 and dens B(5) = 0.938 for 3x+1; the smallest 3x-1 basin is the five-cycle's, so the sheet cap on cycle-uniform sheet-blind lower bounds is the limit of the {5,7,10} basin density, >= 0.3250 on the evidence to 2^32 (still rising), not 1/3"
status: >
  OPEN. FINITE-EXACT support: counting densities on [1, 2^29] (3x-1) and [1, 2^30]
  (3x+1), stable to three decimals across dyadic ranges (THM-4517). Existence of the
  natural density of any Collatz basin is open (it implies, for the plus sheet, that
  almost every orbit reaches 1 if the trunk-entry densities e_i sum to 1). A proof
  would fix the exact constant that no sheet-blind cycle-uniform basin bound can exceed.
  Drift test to 2^32 (dyadic rows): the {1} basin is flat at 0.32675 from 2^28 on; the
  {5,7,10} basin still rises 5e-5 per doubling (0.3250474 at [2^31,2^32)) and the
  {17,...} basin falls at the same rate (0.3482025); the sheet cap is therefore between
  0.3250 and the limit of the rising basin, and the first gap (1.7e-3 at 2^32) is still
  closing. Next test: 2^36 with a 2-bit sieve (16 GB) or a residue-class extrapolation.
source: collatz-necklace-20260929 session (mac-mini), 2026-09-29
related:
  - 01-canon/theorems/THM-4517-no-root-uniform-positive-density-of-collatz-basins.md
  - 05-knowledge/results/collatz_necklace_20260929_fair_splits_power_clocks_basins.md
  - 05-knowledge/results/mazur_positive_density_20260928.md
---

# HYP-9165 -- basin densities exist; the sheet cap is 0.3248

See THM-4517 section 3 for the tables. The hypothesis is the existence of the limits and the strict ordering `dens B({5,7,10}) < dens B({1}) < dens B({17,...})` on the minus sheet. Its consequence for the barrier atlas: a sheet-blind argument giving every positive cycle basin lower density `>= c` is impossible for `c > 0.3248` (and for `c > 1/3` unconditionally, THM-4517(2)).
