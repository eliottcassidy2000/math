---
id: HYP-9772
title: "Single-kernel conjecture: every odd N >= 36835 equals H(A) + F2(A) = 3H(A) + P2(A) for some tournament A, i.e. is the Hamiltonian-path count of the double cone a -> A -> c -> a; fixed-kernel line coverage of the pair spectrum is a one-dimensional statement, true through order 10 exactly and across the order 10 -> 11 -> 12 junctions to 1,238,097"
status: >
  OPEN (CONJECTURE), klein-2026-11-15-lines.
  FINITE-EXACT: order 10 alone covers every odd N in [36835, 241521] (all 9,733,056 classes). Exact data alone (the
  census of orders <= 10 plus the exact reducible pairs of orders 11-14) miss no odd multiple of 3 up to 1.4e6.
  WITNESSED (explicit tournaments from deterministic tilted sampling, spot-checked by an independent Python
  recount): every odd N in [36835, 1238097].
  PROVED: the order-parity obstruction for odd readers (from THM-5172 (C)), and the derived-kernel reduction.
  NOT CLAIMED: Chan-Pak (OPEN); any local limit theorem.
session: klein-2026-11-15-lines
date: 2026-11-15
source: 05-knowledge/results/chanpak_fixed_kernel_line_coverage_20261115/LINES.md
scripts:
  - 05-knowledge/results/chanpak_fixed_kernel_line_coverage_20261115/paircensus.c, tiltsample.c, verify_witnesses.py
  - 05-knowledge/results/chanpak_fixed_kernel_line_coverage_20261115/linecov.py, multiplicity.py, junction11.py, reducible.py, union_ext.py
related:
  - wave 35 of chanpak_waves_20261113 (pair-spectrum fullness reduction), THM-5222 and THM-5172 (C) (cyclic join),
    THM-5127 (B) (tilted ensemble), HYP-9692 (band fullness, the 2D version), HYP-9722 (finite carrier conjecture)
---

# HYP-9772 - the single-kernel conjecture

**Conjecture.** For every odd N ≥ 36,835 there is a tournament A with

    N = 3H(A) + P₂(A) = H(A) + F₂(A),

where F₂ = P_A(2) counts ordered 2-slot path covers, empty slots allowed. Equivalently, by the cyclic-join
identity, N = H(a → A → c → a): a one-vertex reader a, the core A, and a connector c.

**Why it is one-dimensional.** Each pair (H, P₂) lies on exactly one line 3x + y = N. So covering every line is the
surjectivity of one integer statistic, which needs about log N pairs per line. Full two-dimensional occupancy of the
pair band needs about N. At order 10 the bulk lines already carry 4–18 pairs each.

**The order-parity obstruction (PROVED).** One-vertex readers have odd order, so by THM-5172 (C) odd-order cores
give only N ≡ 0 (mod 3). The classes N ≢ 0 (mod 3) therefore come from even orders, and windows two orders apart
must overlap.

**Derived kernels (PROVED reduction).** For the core A ⇒ TT_k, of order n + k, the K₁ statistic equals
H(A) + 2^k F₂(A), because (H, F₂) is multiplicative under ordinal sums.
- So the conjecture follows from bulk coverage of the five statistics H + 2^k F₂, k = 0, ..., 4, at every order
  n ≥ n₀, with window ratios above 2.
- The classes N ≢ 0 (mod 3) use n + k even.
- Each of these is a one-dimensional local-limit statement for the tilted ensemble of THM-5127 (B).

**Evidence.**
- *Order 10.* Every odd N in [36,835, 241,521], with at least 4 pairs per bulk line.
- *Orders 11 and 12.* Explicit witnesses from tilted sampling; the hard even-order junction saturates (505 → 7
  empty lines over 4.8M targeted samples).
- *Union.* Exact reducible pairs at orders 11–14 together with the samples cover every odd N in
  [36,835, 1,238,097].
- *Residues.* All residues mod 3, 5, 7, 9, 16 and 25 occur.

The analogous statement for the two-vertex reader (N = 5H + 3P₂) holds on [188,269, 1,597,157].

**Local limit law, 2026-11-15 (klein-2026-11-15-llt, [LLT.md](../results/chanpak_generic_lines_llt_20261115/LLT.md)).**
The exact per-line class counts of 3H + P₂ fluctuate like √(D/μ) around a smooth density, with D ≈ 2–4.
- The relative fluctuation is 1.54, 0.23 and 0.17 at orders 8, 9 and 10, and the worst bulk line carries 0.25–0.37
  of the smooth density.
- The singular series is flat: residue classes modulo 8–32, 3–27, 5, 7 and 24 are equal within 1%.
- Unit changes of S are vanishingly rare among local moves, so no single-move switching proof exists.

The analytic input this conjecture needs is a ratio local limit law, of L² or bilinear type.
