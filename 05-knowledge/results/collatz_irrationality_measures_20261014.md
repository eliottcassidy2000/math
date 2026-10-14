# Collatz through irrationality measures: the parity-number dictionary and the convergent budget

**Session:** mac-mini-2026-10-14-irrmeasure (worktree `math-wt-irrmeasure`). **Owner prompt (2026-10-14):** view the
positive integers as irrationality measures; in the accelerated Collatz map `1 → 2 → 1` is the only abnormality of the
connected tree; hypothesis `μ = 2` for all odd zeta values and all logarithms (arXiv:2610.10192); sources to consider:
the Kaczynski works repository, "think Fisher information" (statistics aims for Dirac distributions), the automorphisms
of symmetric and alternating groups, the owner's Ramanujan-prime repository, the 7/8 zero-free half-plane paper
(hexagonmath 2610.00200), and the openai/math primitive-roots preprint; plus a pasted essay drawing the analogy.
**Status:** PROVED (THM-4632 dictionary, THM-4633 budget theorem with a self-contained `4.95·10^10` odd-element bound)
+ FINITE-EXACT (checks, convergent ledger) + CITED (Lagarias, Terras, Jarník, Eliahou, Hercher, Barina) + OPEN
(HYP-9262). Audit pending at the first checkpoint (section 8).

## 1. Inheritance pass

* **Closest proved mechanisms.** THM-4591 (integer `T`-cycles to period 301,993 and Ellison's shapes: the convergent
  `301994/190537` is its period bound); THM-4484/4516 (free and sporadic cycles, Fermat–Catalan clocks); THM-4623
  (order dictionary, Koopman torsion); THM-4523 (rigidity); THM-4629/4631 and HYP-9260/9261 (opus S29/S30, 2026-10-13:
  convergents and semiconvergents of `θ = log_2 3` as resonance frequencies, a Möbius dictionary of parity words, the
  Farey census of `θ`, Artin/Wieferich multiplicities of the `px+1` transfer operators).
* **Canonical hostile examples.** The negative and rational cycles (`{-1}`, `{-10,-7,-5}`, the 11-cycle through `-17`,
  `1/5 → 4/5 → 2/5`): every statement about "the only cycle" must be stated on `N`, not on `Z_2`.
* **Corrected near miss.** MISTAKE-602 (a semiconvergent called a convergent): THM-4633 uses Legendre's theorem, which
  yields convergents only below the threshold `k^2 ≤ (3/2) m ln 2`, and cites Hercher for the intermediate-fraction
  analysis above it.
* **Least-used sidecar.** The parity-vector map `Q` as a 2-adic *isometry* (Terras's bijection mod `2^k`), which
  transports real Diophantine strata to `Z_2` without loss.

## 2. Concept board

1. `ρ(x) = Σ v_i 2^(-i-1)`: the parity vector of `x ∈ Z_2` read in binary; `ρ(Tx) = {2ρ(x)}`.
2. `μ(ρ(x))`: the irrationality exponent of the parity number; `1` = eventually periodic, `2` = generic, `m` = Jarník stratum.
3. `D_3 = {a/(3·2^j)}`: the parity numbers of orbits reaching `{1, 2}`; `{2/3, 1/3}` = the trivial cycle.
4. The budget `m_max(p/q) = q/(3(p ln 2 - q ln 3))`: the largest minimum a cycle with ratio `K/k = p/q` can have.
5. The verified limit `M_0 = 2^71` as the price every convergent must pay.
6. Divergence = irrational parity number with odd-digit density `≥ log_3 2`; its "information cost" `1 - H(log_3 2)`.

## 3. Results

### 3.1 THM-4632 — the parity-number dictionary (PROVED; KNOWN ingredients assembled)

`ρ` conjugates `T` to doubling; `ρ(x) ∈ Q` iff the orbit is eventually periodic, which forces `x ∈ Q ∩ Z_2` with
`x_0 = Σ 3^(k-1-i) 2^(d_i)/(2^P - 3^k)` (Lagarias); the converse is the periodicity conjecture. `ρ(1) = 2/3`, `ρ(2) = 1/3`,
and **Collatz ⟺ every `ρ(n)` has lowest-terms denominator `3·2^j`**; nontrivial cycles are `ρ`-values with odd denominator
part `≠ 3` dividing `2^P - 1`; divergence is irrationality. The strata `{μ(ρ(x)) ≥ τ}` have 2-adic Hausdorff dimension
`2/τ` (Jarník–Besicovitch through the isometry `Q`), `μ = 2` Haar-a.e., and `μ = 1` is the countable set of eventually
periodic points: the integers `m ≥ 2` label strata of "constructed continued fractions" (`a_(n+1) ≍ q_n^(m-2)`), and the
conjecture places `N` in the stratum `μ = 1`. Checks: Terras bijection to `k = 14`; denominators `3·2^j` for `n ≤ 2000`;
247 period words of length `≤ 7` give rational cycles, the integer ones being `{-1}, {1,2}, {-10,-7,-5}` (and `{0}`).

### 3.2 THM-4633 — cycles as rational approximations with a budget (PROVED + FINITE-EXACT + CITED)

`2^K/3^k = Π(1 + 1/(3n_i))` gives `0 < K/k - θ < 1/(3 m ln 2)`; a convergent `p_j/q_j > θ` admits only cycles with minimum
below `m_max ≍ a_(j+1) q_j^2/(3 ln 2)`. With `m ≥ 2^71` the convergents `2/1` (budget 1.16: the trivial cycle, `m = 1`),
`8/5`, `65/41`, `485/306`, `24727/15601`, `125743/79335`, `301994/190537`, `17087915/10781274`, `272500658/171928773`,
`630138897/397573379`, `10439860591/6586818670` are excluded and `217976794617/137528045312` survives; Legendre's theorem
then gives, self-contained, that every nontrivial cycle has more than `4.95·10^10` odd elements, and Hercher's Corollary 29
gives `≥ 137,528,045,312 = q_23` — the same convergent. An effective exponent `|θ - p/q| > c q^(-μ_0)` gives
`k > (3 c M_0 ln 2)^(1/μ_0)`: `7.0·10^10` for `μ_0 = 2`, `1.7·10^7` for `3`, `1.7·10^4` for `5.117`, `329` for `8.616`; the
actual continued fraction (partial quotients `≤ 55` through index 40) already behaves like `μ = 2` in the verified window.

### 3.3 HYP-9262 — the owner's `μ = 2` hypothesis (OPEN)

Recorded with the literature: Liu–Jiang–Zhang (arXiv:2610.10192, 7–8 Oct 2026, unrefereed) claim `μ(log r) = 2` for
positive rationals `r` — not the ratio `θ`; known bounds `μ(log 2) ≤ 3.57455`, `μ(log 3) ≤ 5.125`, `μ(ζ(3)) ≤ 5.513891`,
`μ(π) ≤ 7.1032`; `ζ(5)` not known irrational. Collatz reading: `μ(θ) = 2` would force `k ≳ sqrt(M_0)` odd elements,
exactly the scale the ledger attains; it says nothing about divergence.

## 4. The linked sources, typed

| source | what it is | connection found | type |
|---|---|---|---|
| arXiv:2610.10192 (Liu–Jiang–Zhang) | `μ(log r) = 2` for `r ∈ Q_(>0)`, determinant/interpolation method | supports HYP-9262(b) for `log 2`, `log 3`; not for `θ` | CITED, unverified claim |
| Kaczynski works (micah541 repository, Lean 4) | boundary functions; sets of curvilinear convergence are `F_σδ`; Baire classes; Wedderburn; Sutcliffe's digit-reversal conjecture | the Collatz convergence set in `Z_2` is a countable `F_σ`, the strata `{μ ≥ τ}` are `G_δ`: "sets of convergence are low Borel" is the shared pattern; no mechanism transfers | ANALOGY (Borel class) |
| Fisher information / Dirac | zero-variance limit = exact knowledge | rational parity number = deterministic itinerary; `μ = 2` = law of large numbers of the Gauss map; divergence costs `1 - H(log_3 2)` bits per step (THM-4632 reading) | DICTIONARY, no theorem |
| Automorphisms of `S_n`, `A_n` (`S_6` exceptional) | a unique exception in a family | rhetorical parallel to "the trivial cycle is the only abnormality"; the dictionary's version is "`3` is the only odd denominator" | ANALOGY only |
| owner's Ramanujan-prime repository | `N_R = 38,358,837,683` without RH, from the OpenAI `7/8` zero-free half-plane and Bhattacharjee's correction of Büthe's `E_2` | none to Collatz; records the dependency chain on the unrefereed zero-free region | CITED, no connection |
| hexagonmath 2610.00200 | `L_F(s, r)` zero-free for `Re s > 7/8` for all Weil-group representations, extending the OpenAI cubic-theta method | effective Chebotarev would sharpen THM-4631's Artin/Wieferich multiplicities statistically; nothing used | CITED, unverified |
| openai/math primitive roots (4 Oct 2026) | claims Artin's primitive-root conjecture (infinitude, `≥ c_a x/(log x)^2` primes in `(x, 2x)`) for every admissible base, via "the Hecke zero-free argument" | base 2: infinitely many `p` with fast-mode multiplicity `[(Z/p)^×:⟨2⟩] = 1` in THM-4631, and infinitely many candidate rigid `px+1` maps for THM-4523's criterion (which needs mod `p^2`, i.e. non-Wieferich too) | CITED, unverified; conditional remark |
| pasted essay | analogy integers ↔ measures | its "`μ = 4` by `a_(n+1) = a_n^2`" is wrong (`a_(n+1) = a_n^c` gives `μ = c + 1`, so squaring gives `3`; Cahen's constant `μ = 3` is right, Duverney–Shiokawa 2020); its "engineered precision" intuition is THM-4633's budget | CORRECTED |

## 5. Computations and reproduction

    cd 04-computation/experiments/irrmeasure_20261014
    python3 log23_convergents.py 71     # convergents of log2 3, budgets, 2^71 exclusion, Legendre threshold, mu_0 table
    python3 parity_number_checks.py     # Terras bijection, denominators 3*2^j, rational cycles from words, doubling orbits

Hostile controls: the ledger's first surviving denominator equals Hercher's published odd-element bound; the integer
cycles found among short words are exactly those of THM-4591 with period `≤ 7`; the trivial cycle saturates the product
identity (`4/3 = 1 + 1/3`).

## 6. What was not done

No new cycle bound beyond Hercher's; no finite irrationality-exponent bound for the parity number of a hypothetical
divergent trajectory (good approximations need not come from the binary expansion); no use of the zero-free-region
papers (their status is unrefereed and they do not bear on `θ`).

## 7. Next obligations

1. Extend the continued fraction of `θ` far beyond index 40 and record the Khinchin statistics of its partial quotients
   (cheap; tests HYP-9262 numerically and extends the ledger).
2. State Hercher's intermediate-fraction step in the budget language (which fractions between `p_21/q_21` and `p_23/q_23`
   lie within `1/(3·2^71 ln 2)` of `θ`), to make the `137,528,045,312` bound self-contained in-repo.
3. Lagarias's periodicity conjecture in the dictionary: which rational 2-adic integers of small height have rational
   parity numbers (a census of `ρ` on `Q ∩ Z_2` by height).

## 8. Audit

Pending at the first checkpoint; to be filled when the independent audit returns.

## 9. META-PATTERNS candidate

**Owner analogy → conjugacy dictionary.** Trigger: an analogy between a dynamical problem and a classification
by exponents. Move: find the conjugacy that makes the analogy literal (here the parity vector read as a binary
expansion, with the doubling map), name each failure mode by its Diophantine type, and attach the known
quantitative theorem (the convergent budget). Counterindication: when the analogy's objects are not conjugate
(the essay's "measure of an integer"), do not invent an invariant; the dictionary must be a theorem.
