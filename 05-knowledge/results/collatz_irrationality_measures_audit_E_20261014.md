# Audit E (2026-10-14) — THM-4634, THM-4633, HYP-9263, lane file sections 3–4, scripts (session mac-mini-2026-10-14-irrmeasure)

**Auditor:** independent adversarial audit (Claude, worktree `math-wt-irrmeasure`), 2026-10-14. Files audited under their
current ids (the session's THM-4632/HYP-9262 were renumbered THM-4634/HYP-9263; opus-S31's THM-4632/HYP-9262 are out of scope).
**Own computations:** `04-computation/experiments/irrmeasure_20261014/audit_E_{convergents,parity,farey,khinchin}.py` (+ `.out`).
**Sources read:** Hercher, arXiv:2201.00406v3 (full text) and the JIS corrigendum; Zudilin, arXiv:math/0404523 (full text);
Lagarias 1985 survey (HTML, ergodic-theory section); arXiv:2610.10192 abstract; Barina 2025 (J. Supercomput. 81, DOI
10.1007/s11227-025-07337-0) via publisher/VUT pages; Lagarias annotated bibliography II (entry 101); Wikipedia "Collatz conjecture".

## Verdicts

| file | verdict |
|---|---|
| `01-canon/theorems/THM-4634-…parity-number-dictionary…` | **PASS WITH CORRECTIONS** (two false sentences, one proof gap, small slips) |
| `01-canon/theorems/THM-4633-…admissible-minimum-budget…` | **PASS WITH CORRECTIONS** (the Hercher clause and the "first survivor"/"mechanism" claims are wrong; statements 1–4, 5-first-clause, 6-formula are right) |
| `05-knowledge/hypotheses/HYP-9263-…` | **PASS WITH CORRECTIONS** (one false numerical evidence claim; one misrounded literature value) |
| `05-knowledge/results/collatz_irrationality_measures_20261014.md` §3–4 | **PASS WITH CORRECTIONS** (inherits the above; obligation 2 is not achievable as phrased) |
| `log23_convergents.py/.out` | PASS (all numbers reproduced; scope: convergents only) |
| `parity_number_checks.py/.out` | PASS WITH CORRECTION (label in (d)) |
| `log23_khinchin.py/.out` | PASS (numbers reproduced; it does not compute the `10^-3` statistic quoted in HYP-9263) |

## THM-4634 — corrections

1. **Statement 3, last sentence.** Quote: "The converse implication … is Lagarias's periodicity conjecture (OPEN; it contains the
   Collatz conjecture)." Wrong. The Periodicity Conjecture (`Q_∞^{-1}(Q ∩ Z_2) = Q ∩ Z_2`) says every rational 2-adic integer has an
   eventually periodic orbit; it allows nontrivial integer cycles. Lagarias's survey (ergodic-theory section) states only: "the truth
   of the Periodicity Conjecture implies the truth of the Divergent Trajectories Conjecture". Corrected: "(OPEN; it contains the
   no-divergent-trajectories part of the Collatz conjecture, not the cycle part)". Statement 5's last sentence is consistent only
   because of its parenthesis "(with the denominator form of 4)" — keep it.
2. **Statement 6 and proof (6).** Quote: "`1/5 → 4/5 → 2/5 → 1/5` (word `100`, `ρ = 1/7`)". `ρ(1/5) = 0.(100)_2 = 4/7`; `1/7 = 0.(001)_2`
   is `ρ(4/5)` (word `001`, `x_0 = 2^2/(8-3) = 4/5`); `ρ(2/5) = 2/7`. Same mislabel in `parity_number_checks.py` (d). Corrected:
   "(word `100`, `ρ(1/5) = 4/7`; the doubling 3-cycle `{4/7, 1/7, 2/7}`)".
3. **Proof of (3), converse.** Quote: "given a word, the displayed `x_0` has the word as its period word (direct substitution; all 247
   words … checked exactly)". Substitution only shows `x_0` is the fixed point of the affine composition `F_w(x) = (3^k x + S)/2^P`;
   it does not show that `x_0` actually follows the parities `w`. Missing (classical) step: by Terras the set `C_w` of 2-adic integers
   whose first `P` parities are `w` is one residue class mod `2^P`; on `C_w`, `T^P = F_w` maps into `Z_2`, so `C_w ⊆ {x : 3^k x + S ≡ 0
   (mod 2^P)}`, itself one class; hence equality, and `x_0` (which satisfies `3^k x_0 + S = 2^P x_0`) lies in `C_w`, so `T^P x_0 = F_w(x_0) = x_0`.
   The 247-word check is fine as FINITE-EXACT evidence, not as the proof.
4. **Statement 5.** Quote: "partial quotients `a_(n+1) ≍ q_n^(m-2)` along a subsequence". The stratum `μ = m` is
   `limsup log a_(n+1)/log q_n = m − 2`, i.e. `a_(n+1) = q_n^(m−2+o(1))` along a subsequence and not larger; not `≍`. Also
   `dim_H {μ = m} = 2/m` needs Jarník's exact-order theorem (Jarník 1929/1931), not only the `{μ ≥ τ}` version; say so (the citation covers it).
5. **Reading, Fisher paragraph.** Quote: "`1 - H(log_3 2) ≈ 0.049` bits". `log_3 2 = 0.630930`, `H = 0.949956`, `1 − H = 0.05004` → "`≈ 0.050`".
6. **Sharpening (no error).** For every `x ∈ Z_2` whose orbit reaches `{1, 2}`, the lowest-terms denominator of `ρ(x)` is exactly `3·2^L`
   with `L` the hitting time of `{1, 2}`: the last prefix digit is the parity of `T^(L-1)x ∈ {4, 1/3}`, so the numerator `3a + c` is odd.
   Verified for all `n ≤ 100000` (`audit_E_parity.py` (3)). The script comment "prefix may end in the periodic pattern early" is vacuous.
7. **Trivial.** Status says "Jarnik (1931)"; the 2/τ paper is Mat. Sb. 36 (1929) as in `related`.

Verified correct (own computation or argument): Terras bijection `Z/2^k → {0,1}^k` for `k ≤ 16`; `Q` is a 2-adic isometry
(`|Q(x)−Q(y)|_2 = |x−y|_2` on 2000 random 40-digit pairs; in general because `x ≡ y (mod 2^k)` iff the first `k` parities agree);
the non-injectivity set of `ρ` is exactly the orbits reaching `0` or `−1` (eventually constant parity; `Q` injective), e.g.
`ρ(−2) = ρ(−1/3) = 1/2` with `T(−2) = −1`, `T(−1/3) = 0`; the cycle formula with `d_i` counted from `x_0` on `{1,2}` (word `10`),
`{−1}` (`1`), `{−10,−7,−5}` (from `x_0 = −5`, word `110`), the 11-cycle from `−17` (word `11110111000`, `k = 7`, `2363/(2^11−3^7) = −17`),
`1/5` (word `100`); statement 4 (`D_3` = lowest-terms denominators `3·2^j` = expansions ending in `(10)^∞`/`(01)^∞`; Collatz ⟺ all
`ρ(n) ∈ D_3`); statement 5's transport: `ρ_*(Haar) = Lebesgue` (Q is a bijection mod `2^k` for every `k`, the coin measure codes to
Lebesgue), and dimension is preserved because the coding map is 1-Lipschitz and the preimage of an interval of length
`ℓ ∈ [2^(−k−1), 2^(−k))` lies in at most two cylinders of diameter `2^(−k) ≤ 2ℓ`, so `H^s(ρ^(−1)B) ≤ 2^(1+s) H^s(B)`; `S_τ` is `G_δ`;
`S_1 ∖ S_2 = ρ^(−1)(Q)` with the convention `μ(rational) = 1`; the density remark (`T^N n ≥ 1` and the odd elements of a
non-eventually-periodic integer orbit tend to infinity, so `Σ ln(1 + 1/(3n_i)) = o(k)` and `liminf k/N ≥ log_3 2`).

## THM-4633 — corrections

1. **Hercher clause (title, status, statement 5 second clause, proof (5) last sentence, Reading; lane §3.2).** Quote: "**at least
   `137,528,045,312`** (Hercher, Corollary 29; CITED), which is `q_23`, the first surviving convergent denominator — the exclusion ledger
   is the mechanism behind Hercher's number" and "Hercher's Corollary 29 reaches `q_23` itself by controlling the intermediate fractions
   between consecutive convergents (his `m`-cycle machinery with the `2^71` input)". Wrong on three counts.
   * (a) The plain inequality does **not** exclude the intermediate fraction `s = (p_21+p_22)/(q_21+q_22) = 114208327604/72057431991` at
     `2^71`: `s − θ = 1.1034·10^-22 < 1/(3·2^71 ln 2) = 2.037·10^-22`, budget `m_max(s) = 4.358·10^21 = 3780.4·2^60 > 2^71`. The Stern–Brocot
     enumeration of **all** fractions with `k < q_23` in the window (`audit_E_farey.py`) finds exactly one at `2^71` (`s`; ten at `2^68`),
     so the ledger — which is Eliahou's method, and Eliahou's method includes intermediate fractions — gives `k ≥ 72,057,431,991`
     (T-period `K ≥ 114,208,327,604`, Collatz length `≥ 186,265,759,595`), not `q_23`. These are precisely the figures in Hercher's own
     text: introduction "the current `K > 7.2·10^10`" (Hercher–Puchert), and Remark 28 "using previous methods one would have to show
     `δ < (K+L)/K < δ + 1.1032·10^-22`, and therefore `X_0 ≥ 3781·2^60`" (ours: `1.1034·10^-22`, `3780.4·2^60`).
   * (b) Hercher's Corollary 29 reads: "If `X_0 ≥ 1536·2^60 = 3·2^69` then every nontrivial cycle contains at least `K > 1.375·10^11` odd
     numbers." It is conditional, its paper-time input was `X_0 = 704·2^60` (Barina's website, 2023; Definition 4), and its mechanism is
     not the ledger: Theorem 27 sharpens the product inequality by the factor `3/4` (`Σ T(n_i) ≤ K·(3/4)/X_0`, hence
     `(K+L)/K < δ + 1/(4 X_0 ln 2)`; threshold `2836·2^60`, ours `2835.3·2^60`), and a five-week C++ computation (Lemma 26 recursion with
     residue classes modulo powers of two) lowers the threshold to `1536·2^60`. The JIS corrigendum concerns the proof of Theorem 21
     (p. 14), not Theorem 27 or Corollary 29.
   * (c) Barina 2025 (`2^71 = 2048·2^60 ≥ 1536·2^60`) satisfies the hypothesis, so `k > 1.375·10^11` is now unconditional; and since `s` is
     the only admissible `k < q_23` at `2^71`, the next admissible denominator is `q_23`. So "`k ≥ 137,528,045,312`" is a correct
     **consequence** of Hercher Cor. 29 + Barina 2025 + the window enumeration (and agrees with Wikipedia's 2025 "217976794617
     (355504839929 without shortcut)" = `p_23`, `p_23 + q_23`), but Hercher's literal statement is `K > 1.375·10^11`, and the ledger is not
     its mechanism. Corrected sentence: "at least `137,528,045,312 = q_23` (Hercher 2023, Corollary 29, conditional on `X_0 ≥ 3·2^69`,
     now met by Barina 2025; the plain ledger leaves the semiconvergent `114208327604/72057431991` alive until `X_0 ≥ 3781·2^60 ≈ 2^71.9`,
     and Hercher excludes it by Theorem 27 (factor 3/4) plus a computer-assisted residue-class analysis; given that exclusion, the
     enumeration of the window shows the next admissible denominator is `q_23`)". The title's "the first survivor is
     217976794617/137528045312" must read "the first surviving convergent".
2. **Statement 2, asymptotics.** Quote: "`m_max(p_j/q_j) = a_(j+1) q_j^2/(3 ln 2) · (1 + O(1/a_(j+1) q_j))`". False. Exact:
   `m_max(p_j/q_j) = q_j (q_(j+1) + q_j/θ_(j+2))/(3 ln 2)` with `θ_(j+2)` the complete quotient; relative correction
   `(q_(j−1)/q_j + 1/θ_(j+2))/a_(j+1) = O(1/a_(j+1))`, not `O(1/(a_(j+1) q_j))`. Numerically: `j = 1`: `1.159` vs `0.481` (ratio 2.41);
   `j = 23`: `5.10·10^22` vs `4.55·10^22` (ratio 1.12). Proof (2) drops `q_(j−1)`.
3. **Statement 6, "`μ_0 = 5.117`" row.** Quote: "the exponent that a three-term linear-form bound `|p + q_1 ln 2 + q_2 ln 3| > H^(-4.117)`
   would give for `θ`, reported in the literature search". The number is (a misrounding of) the two-term bound for `ln 3` itself:
   `μ(ln 3) ≤ 5.125` (Salikhov 2007), `≤ 5.1163051` (Wu–Wang, J. Number Theory 142 (2014)), `≤ 5.116201` (Bondareva–Luchin–Salikhov);
   it concerns `|q ln 3 − p|` and says nothing about `|q ln 3 − p ln 2|`. No three-term bound with exponent 4.117 exists in the sources
   checked. Delete the row (or relabel it "`μ(ln 3) ≤ 5.1163`, irrelevant to `θ`"). The script tabulates `4.117` (→ `185,635`) while the
   theorem reports `5.117` (→ `1.7·10^4`): inconsistent.
4. **Statement 6, Rhin rows.** Zudilin (math/0404523, Theorem 3) states Rhin 1987 as "the irrationality exponent of any nonzero
   `γ ∈ Q log 2 + Q log 3` satisfies `μ(γ) < 8.616`" — stated per `γ`. Whether Rhin's theorem is a uniform linear-form bound
   `|x_1 ln 2 + x_2 ln 3| > H^(−7.616)` (which would transfer to `θ` with exponent `8.616`, so that "a different statement, quoted for
   scale" would be wrong) could not be checked (Rhin's paper inaccessible). Simons–de Weger 2005 do use Rhin's result (Lagarias
   bibliography II, entry 101) together with Laurent–Mignotte–Nesterenko (their NAW 2004 survey), but the exponent "`13.3`" (`μ_0 = 14.3`)
   appears in no accessible source: keep both rows CITED-from-memory and say so in the status.
5. **"NEW here" (status).** The budget `m_max` and the exclusion ledger are Eliahou's 1993 method with the roles inverted (fix the
   fraction, solve for the largest admissible minimum); Hercher's Remark 28 performs exactly this inversion. Eliahou's method also
   handles intermediate fractions, which the ledger omits. The self-contained Legendre route (`k > 4.95·10^10`) is a weakening of
   Eliahou's bound that avoids the intermediate-fraction theory; label it as such, not as new.
6. **Minor.** HYP-9263/lane: "first surviving denominator `1.375·10^11 ≈ sqrt(2^71)·3`": the ledger's first surviving denominator is
   `7.2·10^10 ≈ 1.5·sqrt(2^71)`; Hercher's is `≈ 2.8·sqrt(2^71)`. Lane §7 obligation 2 ("state Hercher's intermediate-fraction step in the
   budget language … self-contained in-repo") is not achievable: in budget language `s` survives until `X_0 ≥ 3781·2^60`; the bound
   needs Theorem 27 and Hercher's computation.

Verified correct: statement 1 (identity, both inequalities, equality only for the trivial cycle); statement 2's inequality for common
multiples; statement 3 (Legendre's hypothesis `|θ − p/q| < 1/(2q^2)` on the reduced fraction; `1/(2k^2) ≤ 1/(2q^2)` handles
non-reduced `K/k`; odd-index convergents are the ones above `θ`); the table (recomputed with mpmath at 300 digits and sympy's exact
continued fraction: all `p_j, q_j`, signs, budgets, statuses at `2^71` agree; `a_0..a_39` max `55 = a_14`); the threshold
`sqrt((3/2)·2^71 ln 2) = 4.9548·10^10`; the self-contained contradiction in (5) is complete (any `k ≤` threshold has `K/k` reducing to
`p_j/q_j` with `q_j ≤ k`, hence `j ≤ 21`, all budgets `≤ 2.17·10^20 < 2^71`); statement 6's formula and the arithmetic of all rows
(`7.0·10^10`, `1.7·10^7`, `1.7·10^4`, `33`, `329`); Eliahou's generators `301994, 17087915, 85137581 = p_13, p_15, p_16`; THM-4591's
`301,993` is `p_13 − 1`; Barina 2025 exists (J. Supercomput. 81, 2025, DOI 10.1007/s11227-025-07337-0; `2^71`) — article number not checked.

## HYP-9263 — corrections

1. **Evidence.** Quote: "`log a_(n+1)/log q_n < 10^(-3)` beyond `n = 1000`". False: 311 indices `n ∈ (1000, 20000)` exceed `10^-3`; the
   maximum is `0.00453` at `n = 1179` (`a_1180 = 565`); the last exceedance is `n = 5113` (`a_5114 = 997`). `log23_khinchin.py` never computed
   this (it prints the global maximum, `1.0` at `n = 2`). Corrected: "`< 5·10^-3` for `n > 1000` and `< 10^-3` for `5113 < n < 20000`".
2. **Status.** "refinements to about 5.117 reported" → `5.1163051` (Wu–Wang 2014), `5.116201` (Bondareva–Luchin–Salikhov).
3. Verified: `μ(log 2) ≤ 3.57455` (Marcovecchio 2009); `μ(log 3) ≤ 5.125` (Salikhov 2007); `μ(ζ(3)) ≤ 5.513891` (Rhin–Viola 2001);
   `μ(π) ≤ 7.103205…` (Zeilberger–Zudilin 2020); Ball–Rivoal; Zudilin's one-of-four; Khinchin; Roth. arXiv:2610.10192 (Liu–Jiang–Zhang,
   v1 7 Oct 2026, v2 8 Oct 2026): the abstract claims `μ(log r) = 2` for every positive rational `r ≠ 1` (separated-weight determinant
   argument; interpolation theorem for logarithmic jets; blow-up with Serre vanishing), with consequences for nonzero rational linear
   combinations of such logarithms, rational affine changes and logarithms of positive rationals raised to integer powers; it does not
   mention ratios of logarithms or `log 3/log 2`. The HYP's description is accurate. Khinchin statistics reproduced at independent
   precision (22,000 and 23,000 digits agree on all 20,000 quotients): geometric mean `2.672405`, `log q_n/n = 1.181941`, `a_18287 = 24,477`,
   ten largest below 25,000.

## Lane file §3–4

Inherits items THM-4633 1–6 and THM-4634 1–2 (§3.1 "all 247 period words … rational cycles" fine; §3.2's "Hercher's Corollary 29 gives
`≥ 137,528,045,312 = q_23` — the same convergent" must be reworded as in THM-4633 item 1). §4: the Liu–Jiang–Zhang row is accurate; the
essay correction is right (`a_(n+1) = a_n^c` gives `μ = c + 1`, so squaring gives 3, verified numerically; Cahen's constant `μ = 3`,
Duverney–Shiokawa 2020); the Borel-class remarks are trivially right. Rows on hexagonmath 2610.00200, the openai/math primitive-roots
preprint and the Kaczynski repository were not checked (typed CITED/ANALOGY there; no claim depends on them).

## Typing

Nothing numerical-only is typed PROVED. THM-4634's PROVED-by-assembly is acceptable once items 1–3 are fixed. THM-4633: statements
1–3, 6 PROVED, 4 FINITE-EXACT, 5-first-clause PROVED modulo CITED Barina — all stand; the CITED second clause of 5 is mis-described
(item 1) and "NEW here" overclaims (item 5). HYP-9263 OPEN is right; its Evidence has one false sentence (item 1).

## Scripts (this audit)

* `audit_E_convergents.py/.out` — continued fraction two ways (mpmath 300 dps, sympy exact), convergents, budgets, exact convergent
  formula and its `O(1/a_(j+1))` correction, all intermediate fractions above `θ` with budgets and status at `2^68`, `3·2^69`, `2^71`,
  first surviving fraction per `M_0`, Legendre threshold, `μ_0` table.
* `audit_E_farey.py/.out` — Stern–Brocot enumeration of every fraction with `k < q_23` in `(θ, θ + 1/(3 M_0 ln 2))` for
  `M_0 ∈ {2^68, 3·2^69, 2^71, 3781·2^60, 2^72}`; Hercher's Remark 28 numbers reproduced.
* `audit_E_parity.py/.out` — cycle formula on the five named cycles, `ρ` values, denominators and `j = L` for `n ≤ 100000`,
  non-injectivity example, Terras bijection `k ≤ 16`, isometry sample, KL number, density identity on `n = 27`, exponent of `a_(n+1) = a_n^2`.
* `audit_E_khinchin.py/.out` — 20,000 partial quotients at 22,000 and 23,000 digits, Khinchin/Lévy statistics, exceedances of `10^-3`.

## Summary

The dictionary theorem (THM-4634) is right in substance: the Terras bijection, the isometry, the rationality/periodicity equivalence,
the cycle formula (checked on all five named cycles), the `3·2^j` form of the Collatz conjecture and the Jarník transport all hold; fix
the "contains the Collatz conjecture" sentence (it contains only the no-divergence part), the `ρ = 1/7` mislabel, and the proof gap for
"every word is a period word". The budget theorem (THM-4633) is right where it is self-contained (identity, budget, Legendre route,
`k > 4.95·10^10`, the table) but wrong about Hercher: at `2^71` the plain ledger leaves the semiconvergent `114208327604/72057431991`
alive (budget `3780·2^60 ≈ 2^71.9`, exactly Hercher's Remark 28 figure), so the ledger yields `k ≥ 72,057,431,991`; the bound
`k ≥ 137,528,045,312` is Hercher's conditional Corollary 29 (`X_0 ≥ 3·2^69`, Theorem 27's factor 3/4 plus a five-week computation) made
unconditional by Barina 2025, not "the exclusion ledger". The asymptotic error term in statement 2 is wrong (`O(1/a_(j+1))`), the
`5.117` row is a misread of `μ(ln 3) ≤ 5.1163`, and the `8.616`/`13.3` rows remain unverified. HYP-9263's literature is accurate (the
Liu–Jiang–Zhang abstract claims `μ(log r) = 2` for positive rationals and rational linear combinations, nothing on ratios) but its
"`< 10^-3` beyond `n = 1000`" claim is false (311 exceedances, max `0.0045`). The "NEW here" of THM-4633 is Eliahou 1993 restated
(Hercher's Remark 28 already inverts the budget); the only self-contained novelty is the weaker Legendre-route bound.
