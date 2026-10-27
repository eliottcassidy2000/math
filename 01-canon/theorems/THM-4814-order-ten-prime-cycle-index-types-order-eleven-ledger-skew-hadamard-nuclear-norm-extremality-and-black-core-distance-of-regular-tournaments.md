---
id: THM-4814
title: "Order-10 prime cycle-index types and the order-11 ledger; skew-Hadamard tournaments are the nuclear-norm extremals and the spectral exception bound; regular tournaments are at least m native moves from a singleton (exactly m black moves at orders 7, 9, and for P_11); doubly regular tournaments make A_{q-1} a Z[(1+sqrt(-q))/2]-module"
status: >
  FINITE-EXACT (klein, 2026-10-27; not independently audited): (A) the nine automorphism / anti-automorphism
  cycle-index types of the 7,906,161 prime tournaments of order 10 (registry-free nauty census in eight parts;
  totals agree with THM-4803's run), and the self-converse decomposition identity reproducing the self-converse
  census 1, 1, 2, 2, 8, 12, 88, 176, 2752, 8784 at every n <= 10. PROVED: the two cycle lemmas of (A); the
  self-converse decomposition identity; the nuclear-norm bounds and their equality cases in (C); the linear-parse
  bound, the exception lower bound and the comparison with THM-4810's score bound; the move lower bound and the
  endpoint/parity rule of (D); the module structure of (E). FINITE-EXACT given CITED inputs: the order-11 ledger (B)
  (OEIS A000568(11) = 903,753,248 and A002785(11) = 279,968 as inputs; A051337(11) and A054946(11) as agreeing
  cross-checks); the nuclear-norm maximizer census n <= 9; the black-move witnesses (every regular class at n = 7
  and n = 9; P_11 and twelve random regular 11-tournaments; the regular 5-tournament has no black move). CITED: the
  concurrent dyadic-module packet's all-tree bound nu <= (2b-3) n ceil(log2 n) (live worktree, not on main at
  filing), the paid-ROOT edge identity and Eisenstein rows, THM-4810 / session74 (native moves, score bound),
  THM-4730 (pure-blue alphabet), THM-1370 / THM-4790 / THM-4805 (gap pairs). OPEN: HYP-9422.
  FORWARD (2026-10-27): the order-11 types and the order-12, 13, 14 ledgers (P(11) confirmed by the labelled-count
  identity, HYP-9422 (4) resolved) are THM-4816 / THM-4817.
session: klein-2026-10-27 (codex/session-overnight-20261027)
source: 05-knowledge/results/order10_prime_types_n11_ledger_spectral_blackcore_20261027/README.md (§1–§5)
scripts:
  - 05-knowledge/results/order10_prime_types_n11_ledger_spectral_blackcore_20261027/scripts/tourn_types.c (census; data/census_n10_part*.txt, data/prime_types_n*.txt)
  - 05-knowledge/results/order10_prime_types_n11_ledger_spectral_blackcore_20261027/scripts/species11.py, species11b.py (+ receipts/species11-*, species11b-*)
  - 05-knowledge/results/order10_prime_types_n11_ledger_spectral_blackcore_20261027/scripts/nuclear.py, nucmax.py, nucmax9sv.py, parse_bound.py (+ receipts/nuclear-*, nucmax-*)
  - 05-knowledge/results/order10_prime_types_n11_ledger_spectral_blackcore_20261027/scripts/blackcore.py (+ receipts/blackcore7-*, blackcore9-*, blackcore11paley-*)
  - 05-knowledge/results/order10_prime_types_n11_ledger_spectral_blackcore_20261027/scripts/misc_checks.py (+ receipts/misc_checks-*)
related:
  - 01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md
  - 01-canon/theorems/THM-4810-h-spectrum-atlas-strong-fibres-exact-low-spectra-complement-parity-invariant-score-monotonicity-paley-budgets.md
  - 01-canon/theorems/THM-4730-pure-blue-strong-alphabet.md
  - 01-canon/theorems/THM-1960-tournaments-compose-from-regular-seeds-the-spectral-substitution-law.md
  - 01-canon/theorems/THM-1875-transitive-skew-char-is-the-binomial-reflection-tangent-revival.md
  - 01-canon/theorems/THM-4793-fermat-numbers-under-collatz-valuation-two-descent-exact-merges-pepin-affine-cycle-type-and-the-five-ladder.md
  - 05-knowledge/hypotheses/HYP-9413-minimum-intransitivity-chain-family-fibonacci-primes-and-the-intransitivity-count-asymptotics.md
  - 05-knowledge/hypotheses/HYP-9422-regular-tournaments-exactly-m-black-moves-from-a-singleton-plus-strong-core-nuclear-norm-maximizers-off-the-hadamard-orders-and-the-order-11-types.md
  - 05-knowledge/results/paid_root_cycle_budget_20261026/README.md
  - 05-knowledge/results/session74_eleven_mechanisms_20261027/session74-paley/REPORT.md
---

# THM-4814 — Order-10 prime types, the order-11 ledger, nuclear-norm extremals, black-core distance, Paley root lattices

Notation: `T` a tournament on `n` vertices, `S_T` its skew adjacency (`+1` for `i -> j`), `ν(T) = ‖S_T‖_*`
the sum of its singular values; a tournament is prime when it has no nontrivial module; native move
`τ_Q(T) = T^op xor E(Q)` for a directed Hamiltonian path `Q`, blue iff the positional reversal of `Q`
is an anti-automorphism (THM-4810, session74, THM-4732); `A(n)`, `S(n)`, `P(n)`, `SC(n)` the numbers of
isomorphism classes of all / strong / prime / self-converse `n`-tournaments; `c(n)`, `c_s(n)`, `l(n)`
the intransitivity (cyclic-triple hypergraph) counts of THM-4803.

**(A) Order-10 types (FINITE-EXACT; lemmas PROVED).** The 7,906,161 prime tournaments of order 10
fall into nine cycle-index types: rigid 7,899,762 (not self-converse) and 5,504 (self-converse,
anti-automorphism type `2^5`); `|Aut| = 3` with types `1^4 3^2` (684 non-SC; 8 SC with anti types
`2^2 6^1 ×2, 2^5 ×1`; 8 SC with `2^5 ×3`) and `1^1 3^3` (184, non-SC); `|Aut| = 5` with type `5^2`
(6 non-SC; 3 SC with anti `10^1 ×4, 2^5 ×1`; 2 SC with `2^5 ×5`). *Cycle lemma:* an automorphism of
prime order `p` of a prime tournament on more than `p` vertices has at least two `p`-cycles (a lone
cycle would be a module, since fixed vertices relate uniformly to it). *Anti-automorphism lemma:* an
anti-automorphism has at most one fixed point and no odd cycle of length `>= 3`. These determine the
shapes of every row of the `n <= 10` type table. *Self-converse decomposition:*
`SC(n) = SC_nonstrong(n) + SC_strong_nonprime(n) + SC_prime(n)`, where the non-strong term counts
dual-palindromic sequences of strong classes and the strong non-prime term is, per self-converse prime
quotient `Q`, `2·#(Aut±(Q)-orbits) − #(Aut(Q)-orbits)` of module assignments, an odd `ψ`-cycle carrying a
self-converse content `sc(x^l)` and an even one a free content `a(x^l)`; it returns the census values
at every `n <= 10`.

**(B) Order-11 ledger (FINITE-EXACT given CITED inputs).** From the cycle indices of the orders
`<= 10` alone: 19,288,658 non-strong classes, 87,588,161 strong classes with a nontrivial module,
3,453,658,999,907,072 labelled strong tournaments with a nontrivial module, 2,768 + 17,893 self-converse
non-prime classes, labelled strong `s_lab(11) = 35,262,092,417,856,512` (= A054946). With
`A(11) = 903,753,248` (A000568) and `SC(11) = 279,968` (A002785): `S(11) = 884,464,590` (= A051337),
`P(11) = 796,876,429` (A280772, new term), 259,307 self-converse primes, `c_s(11) = 442,336,617`,
`c(11) = 447,114,461`, `P_lab(11) = 31,808,433,417,949,440`, `l(11) = 16,901,802,400,214,078`. The species
identity cannot give `A(11)` itself: the order-`n` primes are their own `x^n` term (HYP-9413 (d) is
resolved in this form).

**(C) Nuclear-norm extremals and the spectral exception bound (PROVED; census FINITE-EXACT).**
`‖S_T‖_F² = n(n−1)`, so `ν(T) <= (n−1)√n` for odd `n`, with equality iff `S Sᵀ = nI − wwᵀ` for a `±1`
vector `w ∈ ker S`, iff `T` is a one-vertex deletion of a skew-Hadamard tournament of order `n+1`
(`n ≡ 3 mod 4`; doubly regular tournaments are the case `w = ±1`); and `ν(T) <= n√(n−1)` for even `n`
with equality iff `S Sᵀ = (n−1)I`. Maximizers: two at `n = 4`, six at `n = 7` (`P_7` and the other
deletions), four at `n = 8` (the skew-Hadamard tournaments, among them the source over `P_7`), none
attaining the bound at `n = 5, 6, 9` (216 classes with singular values `√7^6, √15^2, 0` at `n = 9`).
`ν(TT_n) = Σ_k |cot((2k−1)π/2n)| = (2/π) n (ln n + γ + ln(4/π)) + o(n)`. Linear parse:
`ν(TT_k[B_1..B_k]) <= ν(TT_n) + Σ_i (ν(TT_{m_i}) + ν(B_i))`. If `T, T'` differ on `e` pairs then
`ν(T) − ν(T') <= min(4e, √(8en))`, so a doubly regular `T` needs `e >= q²/8 − (q^{3/2}/π) ln q (1+o(1))`
exceptions against any fixed bank, linear or (by the cited all-tree bound) tree-shaped. For linear
parses THM-4810's score bound gives the stronger `e >= (q²−1)/8 − q(b−1)/2`; it is void for regular
tree parses, where the spectral bound stands: the two obstructions have different scopes.

**(D) Black-core distance (lower bound PROVED; witnesses FINITE-EXACT).** For regular `T` on
`n = 2m+1` vertices and any native word after which a vertex has score `0` or `2m`, the word has length
`>= m` (`s'(v) = n−1−s(v) + [start] − [end]` changes `|s(v) − m|` by at most one); on a geodesic the
singleton-to-be is an endpoint of every path, start when below and end when above the middle score,
and `|s(a) − m| ≡ depth (mod 2)`. All-black words of length exactly `m` ending at a singleton plus a
strong core exist for every regular class at `n = 7` (black geodesic words over all singletons:
30,716; 41,046; 32,928) and `n = 9` (15 classes), and for `P_11` and twelve random regular
11-tournaments; the regular 5-tournament is the pure-blue atom `T5` (THM-4730) with no black move at
all (distance 2 with blue moves). Hence the black distance to a singleton plus a strong core is exactly
`m` at `n = 7, 9` and for the tested order-11 tournaments, and infinite at `n = 5`.

**(E) Doubly regular root lattices (PROVED).** If `S² = J − qI` and `S1 = 0` (doubly regular, e.g.
Paley `P_q`), then `θ = (I+S)/2` maps the sum-zero lattice `Λ = Zⁿ ∩ 1^⊥ ≅ A_{q−1}` to itself
(`((I+S)v)_i = 2v_i + 2Σ_{j∈N⁺(i)} v_j`) and satisfies `θ² − θ + (q+1)/4 = 0` on `Λ`; so `Λ` is a module
of rank `(q−1)/2` over `Z[(1+√−q)/2]` on which `S` acts as `√−q`. For `q = 3` this is the Eisenstein
statement `S = ±(1+2ω)` on `A_2 = Z[ω]`, `(1+2ω)² = −3` (paid-ROOT A1, dyadic-module L5–L6, CITED).
Checked for `q = 3, 7, 11, 19, 23`.

Proofs and data: the packet README. Scope: (B) rests on two OEIS inputs; (C)'s all-tree bound is a
concurrent session's; (D) is a witness census beyond the proved lower bound; nothing here is audited.

**Update 2026-10-27 (klein-2026-10-27, THM-4818).** (A)'s anti-automorphism lemma sharpens to the 2 mod 4 law: every cycle of an anti-automorphism has length `≡ 2 (mod 4)` (a pair swapped by an even power is contradictory), so the only anti-types at order 11 are `2^5 1`, `6 2^2 1`, `10 1`; the lemma as stated here is correct but not sharp. (B)'s order-11 types are now computed by symmetry-first enumeration (7 s; 12 types; two exact validations) and the order-12 ledger follows (packet `order11_prime_types_n12_ledger_dipole_census_20261027`; the theorem-level ledger is THM-4816/4817, reserved by mac-mini-2026-10-27-dipole). (D): all 1,223 regular classes of order 11 have all-black 5-move witnesses. (E): the Eisenstein unit `1+2ω` is also the qutrit overlap `3·(F*DF)_{jk}` of the three mutually unbiased bases of `C^6` (external note `external_preprints_mub6_local_smoothing_eth_merge_20261027.md`, section 1).
