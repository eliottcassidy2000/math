---
id: THM-4521
title: "Coherence, braiding and energy transfer for the Collatz parity graphs and the Syracuse level operator: (1) every periodic parity word closes a cycle in the parity graph G_sigma at every level k, so the level-1 inventory {loop at 0 (density 0, vanishing), the 2-cycle (density 1/2, periodic/contracting), loop at 1 (density 1, expanding)} persists at every level and fixes rho_min = 0, rho_max = 1 -- the base cases of THM-4474 cohere upward exactly as the triangle and pentagon axioms do in Mac Lane's theorem, and this is why no residue-class certificate of either kind exists at any level; (2) negation is the half-turn k -> k + L/2 of the unit cycle Z/L_n (since 2^(L/2) = -1 mod 3^n), it conjugates the Fourier coefficients, and it intertwines the level operator with its complex conjugate, R T_n = conj(T_n) R (the braiding/hexagon compatibility of the sheet involution with the recursion); (3) exact energy transfer: ||T_n g||^2 = sum_xi |Ghat(xi)|^2 E_g(xi) with E_g the energy spectrum of omega_n g, so 1/9 <= ||T_n g||^2/||g||^2 <= 1 with the Parseval value 1/3 exactly when the twisted spectrum is uncorrelated with the geometric symbol, and a ridge is a level whose twisted vector has low-frequency mass on the unit cycle; (4) the 3-node and 5-node induced-motif census of G_sigma is identical for 3x+1 and 3x-1 (sheet-blind), by THM-4474(E)"
status: >
  PROVED (all four items are short: (1) the periodic point x_w in Z_2 of a
  periodic word has residues mod 2^k that follow the word, so they close a walk
  of G_sigma; the three level-1 cycles are the words 0, 10, 1; (2) 2 has order
  L = 2 3^(n-1) mod 3^n and 2^(L/2) = -1, so omega_n(k + L/2) = conj omega_n(k)
  and the half-turn R commutes with the circulant; (3) Parseval for the circulant
  convolution; (4) negation is an isomorphism G_sigma -> G_(nu sigma) with
  nu(+) = -). VERIFIED numerically (the session note): the energy identity to six digits at
  every level n <= 12 (covariance +0.0002..+0.0006 for n >= 4, i.e. S20's slowly
  rising Fourier mass); intertwining to 1e-15; closed walks for the five known
  cycle words and the rational points 1/5, 19/5 at k <= 8; motif censuses at
  k = 4, 5 identical for +/- (20/200 and 20/193 types), random strategies differ. NOT independently audited. Nothing here is a
  Collatz step: (1) is the mechanism behind THM-4474's "no certificate", (2)-(3)
  are exact identities on the unit cycle (the circulant of THM-4520), (4) types
  the tournament/motif approach of the owner's prompt as sheet-blind.
  Dictionary of sixes (PROVED, elementary): in the 3n+1 form the in-degree of n
  is 2 iff n = 4 mod 6 (in T-form iff n = 2 mod 3); the multiples of 3 form the
  transient (vanishing) class; the odd non-multiples of 3 are the two classes
  +-1 mod 6 = (Z/6)^x, which the sheet involution swaps, and every prime > 3
  lies in them; the level-2 unit cycle (Z/9)^x = Z/6 is the hexagon on which
  the level-2 operator has characteristic polynomial 63 lambda^6 + lambda^3 - 1
  (THM-4520 with L = 6; 63 = 2^6 - 1 = 9 x 7 is ord_9(2) = 6, the mod-6
  session's row factorization).
source: collatz-necklace-20260929 session (mac-mini), part 4 (2026-09-30); owner's prompt: Mac Lane coherence (triangle, pentagon), tournament 3-/5-vertex classification, the trichotomy {vanishing, periodic, expanding} of Gervacio-Maehara-Ramos, The Pentagon Graph Operator (arXiv:2604.18984), hexagon identities, the number 6
depends_on:
  - 01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md (parity graphs, cycle densities, negation isomorphism (E), Banach periodic points)
  - 01-canon/theorems/THM-4520-collatz-level-operator-is-a-gauss-twisted-circulant-with-spectrum-on-the-half-circle.md (the level operator and its symbol)
related:
  - 01-canon/theorems/THM-4484-free-and-sporadic-cycles.md (the densities 0, 1/2, 1, 2/3, 7/11 of the five integer cycles)
  - 05-knowledge/results/collatz_mod6_20260917_synthesis.md (63 = 2^6 - 1 = 7 x 9, the six-state subshift, multiples of 3 transient)
  - 05-knowledge/results/collatz_coherence_20260930_trichotomy_hexagon_motifs.md (the session note)
note: 05-knowledge/results/collatz_coherence_20260930_trichotomy_hexagon_motifs.md
scripts: 04-computation/experiments/collatz_coherence_20260930_{coherence, periodic_walks}.py
outputs: 05-knowledge/results/collatz_coherence_20260930_{coherence, periodic_walks}.out
hash_basis: raw LF bytes (hash in the note)
audit: NOT independently audited; each item is a few lines from cited canon.
---

# THM-4521 -- coherence of the base cycles, the half-turn braiding, exact energy transfer, motif sheet-blindness

**PROVED + VERIFIED; not independently audited.** Session note: [collatz_coherence_20260930_trichotomy_hexagon_motifs](../../05-knowledge/results/collatz_coherence_20260930_trichotomy_hexagon_motifs.md).

## 1. The base cycles cohere upward

Let `G_sigma^{(k)}` be the level-`k` parity graph of THM-4474 (nodes `Z/2^k`, two edges from `s` to the lifts of `T_sigma(s) mod 2^{k-1}`). A periodic parity word `w` of length `p` with `a` ones has the Banach periodic point `x_w = c_w/(2^p - 3^a) in Z_2` (THM-4471/4474), whose `T`-orbit has parity word `w` repeated. The residues of that orbit mod `2^k` therefore form a closed walk of length `p` in `G_sigma^{(k)}` whose nodes have the parities of `w`, for **every** `k`. So the cycle inventory of `G^{(1)}` — the loop at `0` (word `0`, density `0`: *vanishing*), the 2-cycle `0 <-> 1` (word `10`, density `1/2`: *periodic*, the cycle `{1,2}`), and the loop at `1` (word `1`, density `1`: *expanding*, the cycle `{-1}`) — persists at every level, and with it `rho_min = 0` and `rho_max = 1` for `sigma = +` (THM-4474's statement that Collatz has no residue certificate of either kind at any level). This is a coherence phenomenon in Mac Lane's sense: the consistency (here, the obstruction) enforced on the smallest configurations propagates to all larger ones with no new input. The trichotomy `{vanishing, periodic, expanding}` of the pentagon graph operator (Gervacio–Maehara–Ramos: every graph is pentagon-vanishing, -periodic or -expanding, by pigeonhole on bounded orders) is, for Collatz, the density trichotomy of THM-4484's integer cycles (`0`; `1/2`; `1, 2/3, 7/11`) and of THM-4474's strategies (class (i): all cycles contracting; class (ii): a closed class with all cycles expanding; the rest undecided, with Collatz undecided at every level).

## 2. The half-turn is the sheet involution and braids the recursion

On the unit cycle `Z/L`, `L = L_n = 2 3^{n-1}`, `2^{L/2} = -1 mod 3^n`, so `u -> -u` is `k -> k + L/2` and `omega_n(k + L/2) = e(-2^k/3^n) = conj omega_n(k)`. Hence, with `R` the half-turn: (i) `mu_hat_n(2^{k+L/2}) = mu_hat_n(-2^k) = conj mu_hat_n(2^k)` — the sheet involution `n -> -n` is complex conjugation of the Fourier coefficients, which is why sheet-blind statements are exactly the conjugation-invariant ones; (ii) `R T_n = conj(T_n) R` where `conj(T_n)` has the conjugate gauge: the involution intertwines the level operator with its conjugate, and the recursion respects it (`f_n = conj f_n circ R` at every level). This is the hexagon-type compatibility the owner asked about: the braiding (the half-turn) commutes with the tensor-like structure (the level operator) up to conjugation, exactly.

## 3. Exact energy transfer

For `g` on `Z/L` and `T_n g = G * (omega_n g)`, Parseval gives `||T_n g||^2 = sum_xi |Ghat(xi)|^2 E_g(xi)` with `E_g(xi) = |(omega_n g)^(xi)|^2/L` (`sum_xi E_g = ||g||^2`) and `|Ghat(xi)|^2 = 1/|2 e(xi/L) - 1|^2 in [1/9, 1]`, of mean `1/3`. So `1/9 <= ||T_n g||^2/||g||^2 <= 1`, with equality to the Parseval value `1/3` iff the twisted energy spectrum is uncorrelated with the symbol; the excess is the covariance `sum_xi (|Ghat(xi)|^2 - 1/3) E_g(xi)`, positive exactly when `omega_n g` has more than its share of low-frequency mass on the unit cycle. A ridge (S22 §2d, my part 2) is a level whose twisted vector `omega_n f_{n-1}` is locally slowly varying along the cycle; the first step from the constant vector has zero covariance because `omega_n` itself has a flat spectrum (THM-4520(4)). For the actual cocycle the covariance is small at every level computed (see the note), which is the Parseval-rate observation of S20 restated as an exact identity.

## 4. Motif censuses are sheet-blind

THM-4474(E): negation `x -> -x` is an isomorphism `G_sigma -> G_{nu sigma}` with `(nu sigma)(m) = -sigma(-m)`; for the constant strategies `nu(+) = -`. Hence every isomorphism-invariant statistic of the parity graphs — in particular the counts of induced 3-node and 5-node subdigraphs by isomorphism type, the "tournament classification" of the owner's prompt — coincides for `3x+1` and `3x-1` at every level, while random strategies differ. The motif approach is typed sheet-blind (SHEET control), and cannot see the sign of the intercept.
