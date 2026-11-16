---
id: THM-5253
title: "(A) AUDIT SOUND: the sharp one-vertex growth lemma for span-3 ribbons, H(T) <= H(T+x) <= 39 H(T) (SHARP_GROWTH.md, tail deletion), with the per-letter sharp constants 1, 9, 17, 19, 33, 33, 39, 31 (SHARP_EIGHT.md), is correct; the exhaustive maxima of H(W l)/H(W) over all 8^n literal span-3 words of length n <= 9 reproduce the eight constants letter by letter, with max Z_l/H = C_l - 1; hence HYP-9803 (b), Z_l <= 38 H, is PROVED and the age-structured potential proposed for it is withdrawn as a proof device. (B) PROVED: age-structured caps -- every profile coordinate of the 119-state transfer satisfies x_s(T) <= kappa(s) H(T), where kappa(s) is the exact rational obtained by iterated tail deletion over the last six old sites and depends only on the head and tail ages in s; max kappa = 30 (two-fragment states <= 68/3), so every raw transfer coordinate is at most 30 H(T) (THM-5252 (B) had 26,425); kappa is sharp for exactly 16 states. (C) FINITE-EXACT: the observed maxima of x_s/H over words of length <= 9 are 14, 14, 7 for two-, three-, four-fragment states; three caps grow between lengths 7 and 9 (13/9 -> 25/17, 8/5 -> 44/27, 7/5 -> 25/17), so not every supremum is a finite-length maximum. (D) REFUTED (pre-registered): no history-refined box cone of history length 0..3 with data caps is letter-invariant; value iteration stalls on dead directions (K = 0) or diverges (K = 1). (E) Classical remark: the Riesel covering set {3,5,7,13,17,241} of 509203 is the set of primitive prime divisors of 2^d - 1 over d | 24, six-pronged because 2^6 - 1 = 63 is the Zsigmondy exception."
status: >
  AUDIT (A): klein-2026-11-15-age audited proofs written by other sessions on 2026-11-15/16 (chanpak_target_ribbon_audit_20261116
  /SHARP_GROWTH.md, chanpak_eight_braid_20261116/SHARP_EIGHT.md); both checkers replayed; corroborated by exhaustive data over
  134,217,728 words of length 9 (packet AUDIT.md). PROVED (B): elementary (packet CAPS.md section 2), not independently audited;
  the kappa table is exact by computer (profile_caps.py, exact fractions, 4,096 local orientations per state). FINITE-EXACT (C)
  and REFUTED (D) with pre-registration (packet PREREGISTERED.md, P1-P5; P2 and P5 refuted, P1 and P3 confirmed, P4 split).
  (E) classical (Zsigmondy/Bang), verified by sympy (clocks_check.py). CITED: THM-002 (odd-cycle packing identity), THM-4880
  (the profile transfer), THM-5243 (sealed gadgets, 39 = 5 * 2^3 - 1), THM-5252 (B), CLOCKS.md of
  logarithmic_repairs_covers_20261116.
session: klein-2026-11-15-age
source: 05-knowledge/results/age_structured_caps_audit_20261115/README.md
related:
  - 01-canon/theorems/THM-5243-literal-span-3-growth-decomposition-and-exact-bounded-growth-below-h-100-two-apex-tapes-attain-every-odd-number-below-2-15-except-7-and-21-and-the-seam-product-decoder-verified.md
  - 05-knowledge/hypotheses/HYP-9803-two-apex-completeness-with-overlap-at-most-three-every-odd-number-except-7-and-21-is-the-path-count-of-a-transitive-spine-with-two-apexes-and-the-sharp-bounded-growth-lemma-z-at-most-38h.md
  - 05-knowledge/hypotheses/HYP-9802-literal-span-3-certificates-need-a-bounded-growth-lemma-and-rank-adaptive-memory-short-cheap-cycles-are-cut-within-nine-letters.md
  - 05-knowledge/results/chanpak_target_ribbon_audit_20261116/SHARP_GROWTH.md
  - 05-knowledge/results/chanpak_eight_braid_20261116/SHARP_EIGHT.md
  - 05-knowledge/results/logarithmic_repairs_covers_20261116/CLOCKS.md
---
# THM-5253 — the sharp growth lemma audited, age-structured caps, and the refuted cones

Notation as in THM-4880 and THM-5243: a literal span-3 word builds a tournament one vertex at a time, each new vertex
beating a subset of the last three vertices; the profile state records, for each path fragment of a cover, the age of
its head (0, 1, 2, or −1 once it has left the window unfed) and of its tail (0, 1, 2, or −1 once free); `x_s(T)` counts
the covers with profile `s`, `H(T) = Σ_{complete s} x_s(T)`.

**(A) Audit.** The tail-deletion proof of `H(T+x) ≤ 39·H(T)` is sound (packet AUDIT.md §1, seven steps checked): paths
ending at `x` inject into Hamiltonian paths of `T`; a path `P₁, x, P₂` with `P₂` starting at `b ∈ N⁺(x)` is recovered,
for fixed `b` and fixed deleted block `S` (the maximal terminal block of `P₁` beaten by `b`), from a Hamiltonian path of
`T − S` and a Hamiltonian path of `T[S]`, and `H(T−S)·H(T[S]) ≤ H(T)` by the odd-cycle packing identity; the admissible
`(b, S)` number at most 39 after the one impossible subset `{n−2, n−6}` is removed. The eight per-letter constants
follow from the weighted version over the last six old sites. Independent corroboration: over all words of length
`n ≤ 9` the maximum of `H(Wℓ)/H(W)` for the letters 0–7 is exactly 1, 9, 17, 19, 33, 33, 39, 31 and the maximum of
`Z_ℓ/H` exactly 0, 8, 16, 18, 32, 32, 38, 30. The maximisers (THM-5243 (2), the sealed gadgets) contain the proof's
fixture as the `W′ = ∅` case.

**(B) Age-structured caps.** For a profile with `k ≥ 2` fragments, take as base the fragment with the dead head (if any)
and merge the others in a fixed order by the same deletion: delete from the current path its maximal terminal block
beaten by the next head, then append the next fragment. The deleted blocks `S_2, …, S_k` are disjoint subsets of the
out-neighbourhoods of window vertices, so they lie among the last six old sites, and the cover is recovered from the
final Hamiltonian path of `T − ∪S_i` and the ordered blocks. Hence `x_s(T) ≤ H(T)·Σ_{(S_i)} Π_i e_i(S_i)/H(T[S_i])`,
where `e_i(S_i)` counts the orders of the block ending at the tail of the path merged so far; maximising over the
`2^12` orientations of the back arcs among the last six sites gives `κ(s)`. Table (packet `data/kappa_table.txt`):
`κ = 1` on the 14 complete states; two-fragment states `≤ 68/3`; three- and four-fragment states `≤ 30`. Every raw
coordinate is therefore at most `30·H(T)`. ∎

**(C) Observed caps.** Over all words of length ≤ 9 the largest `x_s/H` is 14 (states `((−1,−1),(2,0))` and
`((−1,−1),(1,0),(2,−1))`); 16 states attain their `κ`; three states' maxima grew from length 7 to length 9
(`((−1,−1),(0,−1),(2,2))`: 13/9 → 25/17; `((−1,−1),(0,−1),(1,1))`: 8/5 → 44/27; `((−1,−1),(0,−1),(1,−1))`: 7/5 → 25/17),
their proved caps being 24, 12, 12.

**(D) Cones.** The box cone `{x ≥ 0 : f_j·x ≤ cap_j(key)·H(x)}` over 138 features, with caps the exact data maxima
grouped by the last `K` letters, is not letter-invariant for `K = 0, 1, 2, 3` (one linear programme per key, letter and
feature); raising the violated caps and iterating stalls (`K = 0`: directions whose image has no complete state) or
diverges (`K = 1`: the `Z` cap passes 49,000 in ten rounds). The deletion proof works where cones fail because it uses
`H(T−S)·H(T[S]) ≤ H(T)`, a fact about sub-tournaments that the state vector does not encode.

**(E) The six-prong wheel.** `2^24 − 1 = 3²·5·7·13·17·241`; the primitive prime divisors of `2^d − 1` for `d | 24` are
3, 7, 5, (none), 17, 13, 241 for `d = 2, 3, 4, 6, 8, 12, 24` and none for `d = 1`; the empty level `d = 6` is
`63 = 3²·7`, the Zsigmondy exception. The Riesel covering of `509203·2^e − 1` over `e mod 24` and the Sierpiński covering
of `78557·2^e + 1` over `e mod 36` were replayed. No tournament statement is implied (packet THREADS.md, T4–T5).

**What is new.** The audit with its exhaustive corroboration, the cap theorem (B) with constant 30, the growing-caps
phenomenon (C), the pre-registered refutation (D), and the Zsigmondy reading (E). Known and used: THM-002, THM-4880,
THM-5243, THM-5252, SHARP_GROWTH.md, SHARP_EIGHT.md, CLOCKS.md.
