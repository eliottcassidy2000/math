---
id: THM-4481
title: "An entropy law for sign flips: every stationary law of every q n +- 1 sign strategy has 1 - h(pi(odd)) <= merge entropy <= pi(R u R*) <= pi(odd); hence rho_max >= 0.2270922 always, no q n +- 1 with odd q >= 23 has a bounded-lookahead-provable sign strategy at any level, and provable 5n+-1 strategies need constant flip mass (the Haar question for 5n+1 reduces to concentration)"
status: >
  PROVED + INDEPENDENTLY AUDITED; FINITE-EXACT and VERIFIED data.
  In the strategy cube of THM-4474 for odd q (T_sigma(n) = n/2 or
  (qn + sigma(n mod 2^k))/2), let R be the flipped residues (sigma = -)
  and s* = s - 2q^(-1) mod 2^k the partner of a flip s: the unflipped odd
  residue whose target pair coincides with that of s. Let pi be any
  stationary law of the uniform-lift chain on G_sigma, and p = pi(odd).
  (1) Gain identity: p = 1/2 - sum_(s in R) pi(s) g(s), where g(s) is the
  number of odd steps the flip at s removes from the next k-1 forced
  steps. It also holds with any theta in place of 1/2.
  (2) Entropy/merge law: 1 - h(p) <= sum over merge pairs of
  (pi(s)+pi(s*)) h(pi(s)/(pi(s)+pi(s*))) <= pi(R u R*) <= p. The chain
  makes 1 bit per step, and the parity sequence can carry only h(p) bits;
  information is destroyed only where a flip and its partner merge.
  (3) Consequences.
  (a) p + h(p) >= 1, so rho_max >= pi(odd) >= p0 = 0.2270922 for every
  strategy of every odd q at every level.
  (b) Class (i) is empty at every level for every odd q >= 23, since
  log_23 2 = 0.2211 < p0.
  (c) Every closed class of a class-(i) 5n+-1 strategy has
  pi(R u R*) > 1 - h(log_5 2) = 0.01391 and pi(R) > 0.00126.
  (d) The Haar distance of a class-(i) strategy is at least
  0.01391/max_(R u R*) (dpi/dU). So 5n+1 can lie in the Haar closure of
  class (i) only if the invariant densities of provable strategies
  concentrate without bound on their flips.
  (e) For q = 3 the law is void (log_3 2 > 1/2). This is why THM-4479's
  flips can be exponentially sparse.
  FINITE-EXACT: 44 <= delta_8(5) <= 50 (RC2 lower bound; certified upper
  set). The least rho_max over level-k strategies is 1/2 for k <= 6 and
  3/7 for k = 7, 8 for q = 5, 7, 9, 11 (q = 7: through k = 9), so 7n+1,
  9n+1 and 11n+1 have no provable strategy at those levels.
  VERIFIED: certified class-(i) 5n+1 flip sets for k <= 19, with Haar
  fraction 0.391 (k=8) down to 0.180 (k=19).
  OPEN: whether 5n+1's Haar distance tends to 0 (HYP-9141); whether any
  q in 7..21 is ever provable. Collatz itself is untouched.
  UPDATE 2026-09-26 (THM-4486): the min-max density is a game value.
  Every strategy has a cycle of density >= log_(q+1) 2, which beats p0 for
  q <= 19. For 5n+1 the min-max density is exactly 2/5 for k >= 15 (the
  sporadic cycle 1,3,8,4,2). Stationary-law (entropy) floors cannot exceed
  1/3. No provable sign strategy exists for 7n+-1 at k <= 22, nor for
  q = 9..21 at k <= 18.
source: collatz-procgen-20260922 session, drift lane (2026-09-26); audited and promoted by the session orchestrator 2026-09-26
depends_on:
  - 01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md
related:
  - 01-canon/theorems/THM-4479-strategy-cube-distance-to-provability-sharp.md (q = 3, where the law is void)
  - 01-canon/theorems/THM-4480-peak-discounted-provability-price.md (the same constant 1 - H(log_q 2) as the exponent of a vanishing price)
  - 05-knowledge/hypotheses/HYP-9141-5n-plus-1-haar-closure-concentration.md
note: 05-knowledge/results/procgen_drift_20260926_positive_drift_provability.md
scripts: 04-computation/experiments/procgen_drift_20260926_{lib,ihs,local,minmax,run}.py and procgen_drift_20260926_engine.c (sha256 in the note, section 10)
script_audit: 04-computation/experiments/procgen_drift_20260926_orchestrator_check.py
output: 05-knowledge/results/procgen_drift_20260926.out
output_audit: 05-knowledge/results/procgen_drift_20260926_orchestrator_check.out
output_sha256: b72eaa28c95dbb642ba25dfed9eac09cb5ce3139ccd1d67ef397397e9fd59f83
hash_basis: raw bytes
audit: >
  The orchestrator checked, line by line and found sound:
  * the proof of Theorem 1 (k-step mixing of the base chain, and the
    perturbation expansion pi = pi Q_0^k + sum_(j<k) (pi Delta) Q_0^j
    evaluated on the odd nodes);
  * the proof of Theorem 2: predecessor structure (a pair has two odd
    predecessors iff s_- is flipped and its partner is not); entropy rate
    1 bit per step; the chain rule run backwards; the conditional entropy
    localized at merge pairs;
  * Corollaries 2.1-2.3, including the Jensen bound
    (a+b) h(a/(a+b)) <= a log_2((a+b)/a) + a log_2 e.
  Independent code (procgen_drift_20260926_orchestrator_check.py, written
  without reading the lane's scripts) confirms:
  * Theorems 1 and 2 on every closed class of 790 strategies
    (q = 3, 5, 7, 9, 23; k = 3..9);
  * p0 = 0.2270922, with log_23 2 < p0 < log_21 2;
  * no class-(i) strategy exists for q = 7 or q = 23 at k = 2, 3, 4
    (exhaustive).
  The lane's full pipeline was re-run (DRIFT_KMAX=19; 1235 s, 280 MB). Its
  output is identical to the committed .out except for timing fields; this
  includes the RC2 lower bound delta_8(5) >= 44 and every certificate to
  k = 19. Script hashes match the note's section 10.
---

# THM-4481 -- an entropy law for sign flips

**PROVED + INDEPENDENTLY AUDITED.** Full note: [procgen_drift_20260926_positive_drift_provability](../../05-knowledge/results/procgen_drift_20260926_positive_drift_provability.md).

## 1. Setting

* **Strategies.** For odd `q` and level `k`, a sign strategy `sigma` on the odd residues mod `2^k` defines `T_sigma(n) = n/2` (even) or `(qn + sigma(n mod 2^k))/2` (odd).
* **Parity graph.** `G_sigma` has nodes `Z/2^k` and edges `s ->` the two lifts of `T_sigma(s) mod 2^(k-1)`.
* **Chain.** The uniform-lift chain takes each edge with probability `1/2`; an adversary appends a fair parity bit.
* **Flips and partners.** `R` is the set of flipped residues. The partner of a flip `s` is `s* = s - 2q^(-1) mod 2^k`. It satisfies `(q s* + 1)/2 ≡ (qs - 1)/2 (mod 2^(k-1))`, so `s` and an unflipped `s*` send the orbit to the same pair. That is a **merge**.
* **Class (i)** (THM-4474 A): every cycle of `G_sigma` has odd density `< log_q 2`. On every closed class, `pi(odd) <= rho_max` (THM-4474 C).

## 2. The two identities

**Theorem 1 (gain).** For every stationary law `pi` of the chain whose appended bit is Bernoulli(`theta`),
`pi(odd) = theta - sum_(s in R) pi(s) g(s)`, where `g(s) = |w(s)_(1..k-1)| - |Phi_(k-1)((qs-1)/2)|`.

*Proof.*
* In parity-word coordinates the unflipped chain `Q_0` shifts left and appends a Bernoulli(`theta`) letter, so `Q_0^k(s, ·)` is the Bernoulli product law for every `s`.
* Write `Q_sigma = Q_0 + Delta`, where `Delta` is supported on the rows of `R`.
* Iterating `pi = pi Q_0 + pi Delta` gives `pi = pi Q_0^k + sum_(j<k) (pi Delta) Q_0^j`.
* Evaluate on the odd nodes: the first letter after `j` shifts is the `j`-th letter of the forced word. ∎

**Theorem 2 (entropy/merge law).** With `p = pi(odd)`:
`1 - h(p) <= sum_(merge pairs) (pi(s)+pi(s*)) h(pi(s)/(pi(s)+pi(s*))) <= pi(R ∪ R*) <= p`.

*Proof.*
* **Predecessors.** Each target pair has one even predecessor. It has two odd predecessors exactly at a merge, i.e. when the flip `s` is in `R` and `s*` is not.
* **Entropy rates.** The stationary node chain has entropy rate exactly 1 bit per step (two distinct successors, each with probability `1/2`), while the parity sequence has rate at most `h(p)`.
* **Chain rule backwards.** `H(Y_0^n | X_0^n) <= k + n H(Y_0 | Y_1, X_0)`. So `1 <= h(p) + H(Y_0 | Y_1, X_0)`.
* **Localization.** `H(Y_0 | Y_1, X_0)` is nonzero only at merges, where it equals the middle sum.
* **Last two inequalities.** `h <= 1`, and the nodes `s, s*` over all merge pairs are distinct odd nodes. ∎

## 3. Consequences

* **Universal floor.** `p + h(p) >= 1` forces `p >= p0 = 0.2270922`, so `rho_max >= p0` for every strategy.
* **Large multipliers.** Class (i) is empty at every level for every odd `q >= 23`: `log_23 2 = 0.2211 < p0`, while `log_21 2 = 0.2277 > p0`.
* **Constant flip mass for 5n±1.** On every closed class of a class-(i) strategy, `pi(R ∪ R*) > 1 - h(log_5 2) = 0.01391` and `pi(R) > 0.00126`. The analogous constants are `0.0605` for `q = 7` and `0.1006` for `q = 9`.
* **The Haar question reduces to concentration.**
  * `Haar(R) >= U(R ∪ R*) > 0.01391 / max_(R ∪ R*) f` and `Haar(R) > 0.01391^2/||f||_2^2`, where `f = dpi/dU`.
  * So 5n+1 is in the Haar closure of class (i) only if the invariant densities of provable strategies blow up on their flips. The certified sets show `max f <= 2.4` for `k <= 14`, i.e. no concentration.
* **Three drift regimes.**
  1. **`q = 3`.** Negative drift, `log_3 2 > 1/2`: the law is void, and the flip distance tends to 0 exponentially (THM-4479).
  2. **`5 <= q <= 21`.** Positive drift with `log_q 2 > p0`: the flips need constant stationary mass. For 5n+1 the Haar distance decreases polynomially on the data (about `3.1/k` to `3.4/k` for `k <= 19`); its limit is OPEN. For `q = 7, 9, 11` no provable strategy exists at the computed levels.
  3. **`q >= 23`.** No provable sign strategy at any level.
* **The same constant twice.** `1 - H(log_q 2)` appears here as a constant stationary mass. In THM-4480 it is the exponent of a vanishing price, because arbitrary edits can see height and residue classes cannot.
