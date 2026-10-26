---
id: THM-4803
title: "Intransitivity is the relation modulo its modular tree: the cyclic-triple hypergraph C3(T) of a tournament determines T up to exactly N(T) = 2^p * prod(k_N!) labelled choices (dualize each of the p prime nodes, permute the k_N children of each linear node), so 2 <= N(T) <= n! with the primes at 2 and the transitive tournament at n!; the strong components are the connected components of C3(T); H, |Aut| and the modular tree are C3-invariants while the scores are not; the number of intransitivity structures is c(n) = 1,1,2,3,9,32,247,3328,93612,4774417 (unlabelled) and 1,1,2,11,208,8563,680426,101747483,28904318404,15824793243139 (labelled), obeying a species equation over the prime quotients with their automorphism-plus-anti-automorphism cycle indices; A000568 through n = 10 is recovered from the prime cycle-index census (23 types through n = 9; 7,906,161 primes at n = 10), extending A280772 by 136655, 7906161 and showing A003505 erroneous; the intransitivity observable loses E[log2 N] = 1.22 bits at n = 10 and one bit in the limit"
status: >
  CITED + PROVED + FINITE-EXACT (klein, 2026-10-26; not independently audited). The fibre description and formula are
  the published theorems of Boussaïri–Ille–Lopez–Thomassé (Discrete Math. 277 (2004) 29–43: a prime tournament is
  determined up to duality by its 3-cycles) and Boussaïri–Chergui–Ille–Zaidi (arXiv:1805.04607, Corollary 57: the number
  of realizations of a realizable 3-uniform hypergraph is 2^{#prime nodes} times the product of the factorials of the
  arities of the linear nodes); CITED, and VERIFIED-EXACT here on every isomorphism class for n <= 7 by exhaustive
  enumeration of the same-cyclic-triples tournaments (c3class.c, 0 mismatches on 530 classes) and through the exact
  counts for n <= 10. PROVED here (elementary, in the source note): the easy direction, the extremes 2 and n!, Lemma C
  (strong components = hypergraph components), the C3-invariance of H and |Aut|, and E[log2 N(T)] -> 1 under the
  uniform labelled law. FINITE-EXACT: the counts c(n), c_s(n), l(n), the fibre spectra, tree-shape tables, c3 level
  widths (n = 9, 10), Slater layers (n = 8, 9), |Aut| strata (n = 10), the prime census (A280772 a(9), a(10)), each by
  two routes (nauty canonical forms of the hypergraphs versus Pólya substitution over the prime cycle-index census)
  agreeing for every n <= 10. The erratum on A003505 (197, 4008 versus 203, 3785 at n = 7, 8) rests on this census and
  on THM-1960's independent n <= 6 values.
session: klein-2026-10-26 (codex/session-tournament-twist-unify)
source: 05-knowledge/results/a000568_intransitivity_families_20261026/INTRANSITIVITY.md (Theorem I, Lemma C, §3–§5) and FAMILIES.md
scripts:
  - 05-knowledge/results/a000568_intransitivity_families_20261026/scripts/tourn_enum.c (nauty census; receipts/census-n3-n9.txt, census-n10.txt)
  - 05-knowledge/results/a000568_intransitivity_families_20261026/scripts/c3class.c (brute-force fibres; receipts/c3class-check.txt)
  - 05-knowledge/results/a000568_intransitivity_families_20261026/scripts/species.py, species10.py (+ species-normal.txt / species-optimized.txt, species10-*)
  - 05-knowledge/results/a000568_intransitivity_families_20261026/scripts/analyze.py (+ analyze-normal.txt / analyze-optimized.txt)
related:
  - 01-canon/theorems/THM-1960-tournaments-compose-from-regular-seeds-the-spectral-substitution-law.md
  - 01-canon/theorems/THM-1805-the-vandermonde-is-a-signed-tournament-sum-intransitivity-cancels.md
  - 01-canon/theorems/THM-1420-no-linear-tournament-invariants-and-the-chirality-that-survives.md
  - 01-canon/theorems/THM-4739-kirchhoff-and-hodge-laws-of-the-tournament-metagraph-intransitivity-currents-second-order-determination-and-the-consecutive-start-port-ledgers.md
  - 01-canon/theorems/THM-462-the-cubic-spectrum-of-tournaments-is-gap-free.md
  - 01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md
  - 05-knowledge/results/session56_rmq_lca_connections_20261025/tournament-hierarchy/REPORT.md
  - 05-knowledge/hypotheses/HYP-9351-second-order-determination-of-the-metagraph-and-the-remaining-consecutive-port-patterns.md
---

# THM-4803 — Intransitivity is the relation modulo its modular tree

Let `T` be a tournament on `V`, `|V| = n`, and `C3(T)` the 3-uniform hypergraph of its cyclic
triples (its intransitivity structure). Let `𝒯(T)` be the modular decomposition tree (Gallai):
linear nodes with `k_N` children (strong components of the node's quotient), prime nodes with a
prime quotient; `p(T)` = number of prime nodes.

1. **Fibre (CITED, verified).** The tournaments `T'` on `V` with `C3(T') = C3(T)` are exactly those
   obtained by permuting the children of each linear node and dualizing or not each prime node;
   these moves act freely, so their number is `N(T) = 2^{p(T)} ∏_N k_N!`. For `n >= 3`:
   `2 <= N(T) <= n!`, `N(T) = 2` iff `T` is prime, `N(T) = n!` iff `T` is transitive.
2. **Lemma C (PROVED).** The connected components of `C3(T)` are the strong components of `T`.
3. **Invariants (PROVED).** `c3`, the multiset of strong-component types, the modular tree with
   prime quotients up to duality, the Hamiltonian-path count `H`, and `|Aut(T)|` are functions of
   `C3(T)`; the score sequence is not. Chirality pairs `(T, T*)` are `C3`-equivalent, which is why
   every reversal-invariant invariant ties on THM-931's wild pairs. Adding the score function
   leaves only balanced-prime dualizations undetermined (`1, 2, 5, 14, 53` classes, `n = 3..7`).
4. **Counts (FINITE-EXACT, species PROVED).** With `c_s` the strong intransitivity classes,
   `c = MSET_{>=1}(c_s)` and `c_s = X + Σ_{[Q]_± prime/duality} Z_{Aut±(Q)}[c]`; evaluated from the
   prime cycle-index census and matching nauty canonical forms for `n <= 10`:
   `c(n) = 1, 1, 2, 3, 9, 32, 247, 3328, 93612, 4774417`, `c_s(n) = 1, 0, 1, 1, 6, 22, 214, 3074,
   90255, 4680547`, labelled `l(n) = Σ_{[T]} (n!/|Aut T|)/N(T) = 1, 1, 2, 11, 208, 8563, 680426,
   101747483, 28904318404, 15824793243139`. None is in the OEIS.
5. **A000568 from families (FINITE-EXACT).** Session56's `A = 1/(1 - S)`, `S = x + Σ_Q Z_{Aut Q}[A - 1]`
   evaluated with the census (23 cycle-index types through `n = 9`; `7,906,161` primes at
   `n = 10`, `5,525` self-converse) returns `A000568(n)` for `n = 1..10`; A280772 (simple
   tournaments) continues `136655, 7906161`; A003505's `203, 3785` are wrong. Exact tables by root
   type, strong-component count, tree shape, fibre `N`, `c3` level (`n = 9, 10`), Slater index
   (`n = 8, 9`) and `|Aut|` (`n = 10`) are in FAMILIES.md.
6. **Information (PROVED + exact).** As a Hartley fibre (MISTAKE-231's frame) the intransitivity
   observable loses `log₂ N(T) ∈ [1, log₂ n!]` bits; under the uniform law `E[log₂ N] = 2.19, 2.97,
   3.04, 2.60, 2.06, 1.64, 1.38, 1.22` for `n = 3..10` and `-> 1`, since the expected number of
   non-trivial modules `Σ_m C(n,m) 2^{-(m-1)(n-m)}` times `n log₂ n` tends to 0; the prime
   fraction is `0.815` at `n = 10`. Per `c3` level the mean loss decreases monotonically and the
   prime fraction increases; primes need `c3 >= n - 2` (Moon) and attain it for `5 <= n <= 10`.

Proofs, tables and receipts: the source packet. The hard direction of (1) is cited, not re-proved;
everything else is proved or exhaustively checked as stated.
