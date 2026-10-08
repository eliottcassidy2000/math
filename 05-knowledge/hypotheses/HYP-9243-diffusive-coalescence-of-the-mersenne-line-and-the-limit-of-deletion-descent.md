---
id: HYP-9243
title: "Diffusive coalescence of the Mersenne line: below a generic odd exponent E, the block of smaller exponents whose deletion chains are absorbed by Terras time T has extent S(T) of order sqrt(T) (heavy-tailed), because the relation level between two clusters is a recurrent walk; truncated at the orbit horizon sigma_T(M_K) ~ 8.64 K this leaves about sqrt(K) final clusters, which is the orphan law HYP-9242 (fraction of order 1/sqrt(K)), and implies that deletion certificates alone cannot ground a generic giant Mersenne number"
status: >
  OPEN. NUMERICAL: coalescence profiles below six random odd exponents in [1e10, 1e11] (deletions D <= 2048, T <= 2^20):
  upper median S(T) = 10 at T = 4096, 136 at 2^17, 304 at 2^20 (median S(T)/sqrt(T) 0.38 and 0.30); 1-4 live clusters among 2048
  deletions at 2^20; the t = 23 fan's lower block has extent 1910 at T = 2.74e6 (sqrt(T) = 1655); exit-depth tail for odd giant
  exponents P(no exit, D <= 64, by W) = 0.54, 0.38, 0.26, 0.16, 0.092, 0.047 at W = 128..4096. HEURISTIC: the step to the orphan
  law and to the grounding limit. HYP-9242's exact constant (fraction x sqrt(K) = 0.89-1.14) lies inside the predicted band
  0.76-1.1 for c = 0.30-0.45; c is poorly determined (six profiles).
source: opus-2026-10-07-S21, 05-knowledge/results/mersenne_line_barriers_20261007.md (section 4)
related:
  - 01-canon/theorems/THM-4605-mersenne-exits-by-pair-chain-absorption-the-t23-fan-escapes-at-depth-1.9e7.md
  - 05-knowledge/hypotheses/HYP-9242-orphan-law-deletion-orphans-decay-like-inverse-square-root-of-log-n.md
  - 01-canon/theorems/THM-4581-haar-coalescence-affinely-related-collatz-orbits-merge-almost-surely.md
scripts:
  - 04-computation/experiments/mersenne_line_coalescence_20261007.py (+ .out) (S), (L), (V)
  - 04-computation/experiments/mersenne_line_barriers_20261007.py (+ .out) (F)
---

# HYP-9243 — diffusive coalescence of the Mersenne line

## Statement

* For a generic odd exponent E, let `S_E(T)` be the largest S such that every deletion chain `D ≤ S` (THM-4605's source-reference chains) is absorbed into the source by Terras time T.
* Then `S_E(T) ≍ √T` in distribution, with a heavy-tailed, jump-like profile.
* Consequently, among odd `K ≤ X` about `√X` exponents are bottoms of final clusters at their orbit horizon. This is the orphan law HYP-9242.
* A deletion descent from a generic giant `M_E` cannot reach the verified range.

## Mechanism

* Two clusters are related by `y = 3^k x + e`, with k their relative odd-step count.
* At flips k moves toward or away from 0 with equal probability (THM-4581 (1b)). So k performs a recurrent walk with no drift, the debt walk of S19/S20.
* Exponents D apart start at `k = −D` and coalesce on time scale about `D²`.

## Evidence

See the status.
* The t = 23 fan's escape (THM-4605) is one realised level-1 merge: depth 1.9·10⁷.
* Its level-2 merge has depth beyond 1.3·10⁸, with the walk spread to −12297, about √n.

## What would refute it

* A source whose absorbed extent grows linearly in T over two decades.
* Or an orphan fraction decaying faster than any power of 1/√K beyond K = 12800.
