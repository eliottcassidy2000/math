---
id: THM-5265
title: "The open half of the converse (digit 0 => near-join count) is a density statement. (A) Extended converse: NJ_ext = {h(T) + 2^k B_t(T) : k >= 1} union {h(T) + B_t(T) : e_t(T) = 0} is exactly the set of counts of tournaments with a non-sink digit-0 vertex and equals {counts} . {counts of strong tournaments with a digit-0 vertex}; NJ_ext minus NJ consists of midpoint values h(A)h(D) + B_t(A)P_D(2) (D strong) only; the odd numbers missing from NJ_ext are 1, 7, 21 and the sterile-only counts (11, 13, 47, 67 below 125); HYP-9453 is equivalent to (every odd v >= 23 other than 47, 67 is the count of a strong tournament with a sole exit at a dead-end port) and (every midpoint value is a near-join value). (B) Same-core obstruction: the midpoint values 15 (D = C_3) and 37 (D the strong order-4 tournament) are not h(D'') + 2^k P_{D''}(2) for any tournament D'' of any order and any k >= 1; identities P_D(2) = h(D) + sum_x B_x(D) and B_w(D + w) = h(D) + sum_{x in N+(w)} B_x(D). (C) Reach: near-join values from all rooted tournaments of orders <= 9 realise every odd number <= 398,129 except 1, 7, 11, 13, 21, 47, 67 (first miss 398,131); the dead-end appends from bases <= 9 add nothing. (D) Above the reach, midpoint values (cores <= 7, strong sources <= 8) are missed by the bases-<= 9 family at the random-odd rate (ratios 1.01, 0.99 in the first two dyadic windows): no rewriting family protects them; D + w rewrites 36 of 43 small sources, sole-exit w 13 of 43, induced-inside-S bases 38 of 62, halved partners 53 of 397. (E) For a prime p >= 5, 63 | 2^(p-1) - 1 iff 21 | 2^(p-1) - 1 iff 9 | 2^(p-1) - 1 iff p = 1 (mod 6), and 7 | 2^(p-1) - 1 iff p = 1 (mod 3): the Mandelbrot hyperbolic-component counts 2^(p-1) - 1 divisible by both holes 7 and 21 are exactly the prime periods p = 1 (mod 6). (F) The ratio of self-converse tournament classes to converse pairs is 2, 4, 0.545, 0.478, 0.0525, 0.0292 at orders 4-9 (A002785 = 2, 2, 8, 12, 88, 176, 2752), maximal at order 5; all strong tournaments of order <= 5 are self-converse."
status: >
  PROVED (klein-2026-11-16-conv2): (A) from THM-5264 (A), (D), (E) and THM-5263 (B), (I3); (B) by the bound
  P_{D''}(2) >= (n+1) h(D'') and the classification of tournaments with count 1, 3, 5 (one strong component C_3 or the
  strong order-4 tournament, Busch's floors); the identities by splitting pairs at the start of the second path; (E) by
  the orders of 2 modulo 7, 9, 21, 63. FINITE-EXACT: (C) (nj9.c, all 198,948 classes of orders 1-9, every root, values
  <= 2^23), the family censuses in (D), (F) (labelg canonical forms to order 9). NUMERICAL: the random-rate statement in
  (D). Not independently audited. CITED: THM-5263, THM-5264, THM-5260 (B), THM-4922, HYP-9453, HYP-2080, Busch's strong
  floors, chanpak_waves_20261113 waves 54-56 (tail formula, midpoint reading, Zsigmondy/Mandelbrot).
session: klein-2026-11-16-conv2
date: 2026-11-16
source: 05-knowledge/results/converse_midpoint_reach_threads_20261116/CONVERSE2.md
scripts:
  - 05-knowledge/results/converse_midpoint_reach_threads_20261116/scripts/nj9.c (+ receipts/nj9.out)
  - 05-knowledge/results/converse_midpoint_reach_threads_20261116/scripts/da_test.py (+ receipts/da_test.out)
  - 05-knowledge/results/converse_midpoint_reach_threads_20261116/scripts/hyp_w.py (+ receipts/hyp_w.out)
  - 05-knowledge/results/converse_midpoint_reach_threads_20261116/scripts/inside.py (+ receipts/inside.out)
  - 05-knowledge/results/converse_midpoint_reach_threads_20261116/scripts/halved.py (+ receipts/halved.out)
  - 05-knowledge/results/converse_midpoint_reach_threads_20261116/scripts/sc_ratio.py (+ receipts/sc_ratio.out)
  - 05-knowledge/results/converse_midpoint_reach_threads_20261116/scripts/primes_near.py (+ receipts/primes_near.out)
related:
  - 01-canon/theorems/THM-5264-the-out-degree-1-law-for-digit-0-vertices-is-proved-at-every-order-two-path-lemma-unique-digit-0-vertex-in-a-strong-tournament-and-the-converse-digit-0-implies-a-near-join-count-below-22228.md
  - 01-canon/theorems/THM-5263-digit-0-vertices-backward-junctions-sole-exits-and-the-near-join-the-primes-11-13-47-67-are-permanent-near-join-exceptions.md
  - 01-canon/theorems/THM-4922-first-creators-of-hamiltonian-path-counts-converse-fibers-endpoint-genealogy.md
  - 05-knowledge/hypotheses/HYP-9453-the-near-join-sweep-toward-chan-pak-completeness-density-of-b-sets-two-strong-realisers-for-every-attained-odd-count-at-least-15.md
  - 05-knowledge/results/chanpak_waves_20261113/README.md
  - 05-knowledge/results/converse_midpoint_reach_threads_20261116/THREADS.md
---
# THM-5265 — the open half of the converse is a density statement

**(A)** A non-sink digit-0 vertex is a sole exit at a dead-end port (THM-5264 (A)), so its count is `h(S − z) + B_{t_0}(S − z)`
with `e_{t_0} = 0`; the strong component of `z` is a strong digit-0 tournament and the other components multiply. Bases
with a sink give `k ≥ 1` near-joins (THM-5263 (I3)), bases without a sink give midpoint values (THM-5264 (D)). The
sterile-only counts are the counts of no dead-end append (THM-5264 (F)); 7, 21 are no counts; 1 has no strong realiser. ∎

**(B)** `P_{D″}(2)` counts ordered pairs of paths covering `D″`; grouping by the start of the second path gives
`P = h + Σ_x B_x ≥ (n+1)h`. So `v = h + 2^k P ≥ 5h` and `h ≤ v/5`. Counts multiply over strong components and strong
tournaments of order ≥ 4 (resp. ≥ 5) have count ≥ 5 (resp. ≥ 9), so `h = 1, 3, 5` force `D″` transitive, padded `C_3`,
padded strong-4, with `P = 2^n, 12·2^j, 32·2^j`; `h = 7` is impossible. For `v = 15`: `14 = 2^m` or `12 = 12·2^m` need `m = 0`
(`k = 0`). For `v = 37`: `36`, `34 = 12·2^m`, `32 = 32·2^m` likewise. ∎ `B_w(D + w)`: the pair `(P_1, P_2)` with `P_2 = (w)`
contributes `h(D)`; otherwise `P_2 = (w, x, …)` with `x ∈ N⁺(w)` and the rest is counted by `B_x(D)`. ∎

**(C), (D), (F)** are censuses (CONVERSE2.md §3–4, THREADS.md §1); the random-rate comparison uses the odd-number gap density
of the same bitmap. **(E)** `ord_7(2) = 3`, `ord_9(2) = ord_{21}(2) = ord_{63}(2) = 6`. ∎

The owner's claim that the open half "can be proven" is not available as a rewriting: (B) and (D) show that any proof must
be arithmetic (a density statement about `B`-sets, HYP-9453 §2), and the equivalence in (A) is its exact form. Chan–Pak OPEN.
