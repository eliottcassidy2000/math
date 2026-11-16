---
id: HYP-9803
title: "(a) Two-apex completeness: every odd N other than 7 and 21 is the Hamiltonian-path count of a transitive spine with two apexes, the second beating at most three spine vertices (so every realizable conflict graph needed for Chan-Pak is a clique with a second apex at overlap at most three); (b) the sharp bounded-growth lemma for literal span-3 words: the number Z_l(W) of two-path covers bridged by the next vertex is at most 38 H(W), equivalently H(W l) <= 39 H(W), with equality at every H that some word attains (W' 0 0 0 7 6 0 followed by {2,3})"
status: >
  OPEN (CONJECTURE), klein-2026-11-15-ribbon. (a) FINITE-EXACT: all 2^23 two-apex tapes with spine 11 attain every odd
  N < 2^16 but 7, 21 (first miss 77,533; the one-apex family alone misses half); with the second apex beating <= 3 spine
  vertices, spine 12 attains every odd N < 2^16 but 7, 21 (first miss 117,345); overlap 1 gives 56 %, overlap 2 gives
  98 % below 2^14, and only 0.24 % of the values below 2^16 need overlap 3; orientation A -> B suffices below 2^16;
  spine 13 was pre-registered (packet section 5).
  (b) FINITE-EXACT: the maxima 1, 3, 5, 11, 23, 39, 39, 39 at lengths 0-7 (maximisers 0^(n-3) 7 6 0, l = {2,3}, with
  H(W) = 1 and Z = 38; sealing by three empty letters reproduces the ratio at every H), the lemma exact (max ratio 39, equality at every attained H) for every word with H(W) <= 40 via the closed set
  of 16,794 states, and PROVED the decomposition
  H(Wl) = e + s + Z_l with e, s <= H(W). Two certificate routes are PROVED to fail (annihilator LP infeasible; coordinate-cap
  cone not invariant), so a proof must use the ages of fragment heads and tails. (b) is the lemma HYP-9802 (ii) needs; (a) is
  the closure lemma UNIFYING.md section 3 asked for, in its strongest finite form.
session: klein-2026-11-15-ribbon
date: 2026-11-15
source: 05-knowledge/results/chanpak_span3_growth_two_apex_seam_20261115/README.md
related:
  - 05-knowledge/hypotheses/HYP-9802-literal-span-3-certificates-need-a-bounded-growth-lemma-and-rank-adaptive-memory-short-cheap-cycles-are-cut-within-nine-letters.md
  - 05-knowledge/results/chanpak_apex_ocf_literal_cones_20261115/UNIFYING.md
  - 05-knowledge/results/chanpak_signed_seam_bank_20261116/PRODUCT.md
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md
---
# HYP-9803 — two-apex completeness, and the sharp bounded-growth lemma

**(a) Two-apex completeness.** Let the spine be the transitive tournament on `0 < 1 < … < m−1`, and add apexes `A`, `B`
with `A → B`, `A` beating the spine vertices in a set `S_A`, `B` those in `S_B`, `|S_B| ≤ 3`. Conjecture: for every odd
`N ∉ {7, 21}` some such tournament has exactly `N` Hamiltonian paths. In the odd-cycle conflict-graph frame (THM-002,
Grinberg–Stanley) the count is `1 + 2(a + b + c) + 4d`, so the statement is that cliques with one further apex at
overlap at most three realise every value the independence polynomial at 2 must take. Evidence: packet §3 (every odd
`N < 2^16` except 7, 21 with spine ≤ 12 and overlap ≤ 3, and with spine ≤ 11 over all tapes); the overlap needed is 0 for
1,605 values, 1 for 19,467, 2 for 11,615, 3 for 79 below `2^16`. Tests: spines 12 and 13 (pre-registered Q4–Q6), and an explicit digital formula for the count when
`|S_B| ≤ 3` (the one-apex count is `1 + p(u)` with `p` the inversion-pair statistic of the apex packet), which would
turn (a) into a statement about three digital functions.

**(b) The sharp bounded-growth lemma.** For a literal span-3 word `W` and letter `ℓ`, `H(Wℓ) = e_ℓ + s_ℓ + Z_ℓ` where
`e_ℓ, s_ℓ ≤ H(W)` count the paths to which the new vertex is appended or prepended and `Z_ℓ(W)` the two-path covers it
bridges (packet §1.1, PROVED). Conjecture: `Z_ℓ(W) ≤ 38·H(W)` for all `W, ℓ`; hence `H(Wℓ) ≤ 39·H(W)` (HYP-9802 (ii)'s lemma with its conjectured
constant). Equality holds for `W = W' 0 0 0 7 6 0` and `ℓ = {2, 3}` for every word `W'` (three empty letters seal `W'`, so
the ratio 39 of the gadget `7 6 0` is reproduced at `H = H(W')`): the constant is sharp at unboundedly large `H`.
Evidence: all words of length ≤ 7 (and ≤ 9 if Q1–Q2 pass), and the closed set of 16,794 states with `H ≤ 40` (every word
with `H(W) ≤ 40`, of any length), where the maximum ratio is exactly 39 at every attained `H`. What a proof
needs: a bound on two-fragment weights by complete weights that uses the ages of the fragments' heads and tails (a young
head must be fed within three letters or the state dies unless it is prepended; a blocked tail frees itself within three
letters), since the polyhedral relaxations by span or by coordinate caps are both refuted (packet §1.4).

**Why they matter.** (a) is the non-branching route to Chan–Pak completeness proposed in UNIFYING.md: one closure
property of realizable conflict graphs instead of more families. (b) removes the readout cap from letter-level cone
certificates for literal span-3 growth (HYP-9742), which together with a rank potential (HYP-9802 (iii)) is the known
scheme for a certificate above the floor `5^(1/3)`.
