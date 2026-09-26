---
id: THM-4480
title: "The price of fixed-horizon provability is the peak-discounted undecided density: rho^peak_L/M*_L <= eps_L(q) <= rho^peak_L for every odd q, so the exponent is 1 - H(log_q 2) for every multiplier (0.013911 for 5n+1, whose undecided density does not tend to 0), and for q = 3 the price is rho_L exp(-Theta(L^(1/3))), answering THM-4478's P1 negatively for arbitrary edits"
status: >
  PROVED + INDEPENDENTLY AUDITED (Theorems 1-3); PROVED modulo one CITED
  theorem (Theorem 4, Mogul'skii 1974); FINITE-EXACT; EMPIRICAL (pairing).
  Let q be odd, T(n) = n/2 or (qn+1)/2, and c = log_q 2. For a parity word
  u put w_j(u) = q^(e_j)/2^j. Bad_L is the set of words with w_j > 1 for
  1 <= j <= L, and rho_L = |Bad_L|/2^L. Let eps_L(q) be the least density of
  edits E(G) = {v : G(v) != T(v)} over maps G under which every n >= 2
  falls below itself within L steps (THM-4478's setting).
  (1) rho^peak_L/M*_L <= eps_L(q) <= rho^peak_L, where
  rho^peak_L = 2^-L sum over u in Bad_L of 1/max_(0<=j<L) w_j(u), and
  M*_L <= L(L+1)(2L-2+3q)/(6q). The upper bound sends each undecided
  source to 1 at the highest point of its first L steps. The lower bound
  is a single-scale integer-capacity count. It also holds for the flipped
  pair indices of the pairing family.
  (2) rho^peak_L <= (q/2) 2^(-(1-H(c))L) for every L, and
  eps_L(q) >= 2^(-(1-H(c))L - O(L^(1/3))). So -(1/L) log_2 eps_L(q) ->
  1 - H(log_q 2): 0.050044 (q=3), 0.013911 (q=5), 0.060510 (q=7),
  0.100619 (q=9). This holds for every odd q, whatever the drift sign. For
  q = 5 the undecided density tends to 0.176 and does not tend to 0.
  (3) 2^((1-H)L) rho^peak_L = exp(-Theta(L^(1/3))), with elementary
  constants. Hence for q = 3, rho^peak_L/rho_L = exp(-Theta(L^(1/3))) and
  eps_L(3) = o(rho_L/L^A) for every A: THM-4478's question P1 has a
  negative answer for arbitrary edits. The periodic strategy cube is
  polynomially sharp (THM-4479), because residue classes cannot see
  height.
  (4) [modulo Mogul'skii] ln rho^peak_L = -(1-H(c)) L ln 2 -
  kappa_q L^(1/3)(1+o(1)), where kappa_q = (3/2)(pi^2 c(1-c))^(1/3)
  (ln z_q)^(2/3) and z_q = q min(1, c/(1-c)); kappa_3 = 2.1076,
  kappa_5 = 2.4360, kappa_7 = 2.4104.
  UPDATE 2026-09-26: for q = 3 the sharp second-order term holds with an
  elementary O(log L) error, no Mogul'skii needed:
  ln rho^peak_L = -(1-H)L ln 2 - kappa_3 L^(1/3) + O(log L). This was proved
  independently by the crossroads223 bridge note and by the robin lane
  (procgen_robin_20260926_robin_inequality.md, Corollary 5; sine sub- and
  supersolutions of the letter walk in a strip).
  EMPIRICAL: in the pairing family a partner-isolated peak catch tracks
  2 rho^peak within a factor 1.2-1.5 for L = 8..32, and every n <= 10^6
  descends. Whether delta_L (pairing) <= poly(L) rho^peak_L is OPEN.
  These are fixed-horizon modification prices. Collatz is OPEN.
source: collatz-procgen-20260922 session, peak lane (2026-09-26). It proves the incoming-digest proposal "catch at the peak" and closes the orchestrator's open q = 5 exponent (procgen_wave13_20260926_orchestrator_findings.md, section 2). Audited and promoted by the session orchestrator 2026-09-26.
depends_on:
  - 01-canon/theorems/THM-4478-critical-tube-affine-capacity-sharp-provability-price.md
related:
  - 01-canon/theorems/THM-4475-price-of-provable-descent-tends-to-zero.md
  - 01-canon/theorems/THM-4479-strategy-cube-distance-to-provability-sharp.md (periodic contrast)
  - 01-canon/theorems/THM-4477-cauchy-schwarz-price-bound-and-am-fair-criticality.md
  - 05-knowledge/results/procgen_wave13_20260926_orchestrator_findings.md
note: 05-knowledge/results/procgen_peak_20260926_peak_discounted_price.md
scripts: 04-computation/experiments/procgen_peak_20260926_{lib,integers,pairing,run}.py (sha256 in the note, section 10; re-checked unchanged)
script_audit: 04-computation/experiments/procgen_peak_20260926_orchestrator_check.py
output: 05-knowledge/results/procgen_peak_20260926.out
output_audit: 05-knowledge/results/procgen_peak_20260926_orchestrator_check.out
output_sha256: 9144dfecbbb30e9ee74a550db827393323ba04ead6c02b5bb8650029805ff679
hash_basis: raw bytes
audit: >
  The orchestrator read Theorem 1 (upper and lower bound), Theorem 2 (the
  majorant and the Chernoff tilt z = qc/(1-c) >= 1, which is
  c(q+1) >= 1), Lemma C and Theorem 3 line by line and found them sound.
  Theorem 4 relies on the cited Mogul'skii lemma, whose statement the lane
  took from Gantert-Hu-Shi 2011, Lemma 2.1. Independent code
  (procgen_peak_20260926_orchestrator_check.py, written without reading
  the lane's scripts) confirms:
  * an exact DP over (e_j, argmax) equals direct enumeration for q = 3, 5, 7;
  * rho_L and rho^peak_L equal the note's values (q = 3, L = 8..64;
    q = 5, L = 16);
  * Theorem 2's chain holds for q = 3, 5, 7, 9 and L <= 48;
  * the peak construction makes every n <= 3*10^5 descend within L
    (q = 3, 5, 7), within the counting bound of Theorem 1;
  * M*_L satisfies its closed-form bound;
  * rho_L(5) stays near 0.19 while rho^peak_L(5) decays.
  The lane's full pipeline was re-run (272.5 s, 369 MB). Its output is
  identical to the committed .out apart from the [time] lines. Script
  hashes are unchanged from the note's section 10.
---

# THM-4480 -- the peak-discounted price of provability

**PROVED + INDEPENDENTLY AUDITED** (Theorems 1–3); Theorem 4 modulo a CITED lemma. Full note: [procgen_peak_20260926_peak_discounted_price](../../05-knowledge/results/procgen_peak_20260926_peak_discounted_price.md).

## 1. Statement

Let `q >= 3` be odd, `T(n) = n/2` or `(qn+1)/2`, and `c = log_q 2`. For a word `u`, write `e_j(u)` for the number of ones among its first `j` letters and `w_j(u) = q^(e_j)/2^j`.
* `Bad_L = {u : w_j(u) > 1 for 1 <= j <= L}`, and `rho_L = |Bad_L|/2^L`.
* `w*(u) = max_(0<=j<=L-1) w_j(u)`, and `rho^peak_L = 2^-L sum_(u in Bad_L) 1/w*(u)`.
* `eps_L(q)` is the least density of `E(G) = {v : G(v) != T(v)}` over all maps `G` with `G^j(n) < n` for some `1 <= j <= L`, for every `n >= 2`.

**Theorem 1.** `rho^peak_L / M*_L <= eps_L(q) <= rho^peak_L`, where `M*_L = sum_(k<L) sum_(e=floor(kc)+1)^(k) (floor(e/q)+1) <= L(L+1)(2L-2+3q)/(6q)` (with the `k = 0` term equal to 1).

* **Upper bound (catch at the peak).**
  * Let `P(n)` be the first maximiser of `T^j(n)` over `0 <= j <= L-1`, for every `n` that does not descend within `L` steps. Put `G = 1` on `E = {P(n)}` and `G = T` elsewhere.
  * Every source either descends by itself or meets `E` no later than its peak time `<= L-1`.
  * On the class of a bad word, `T^j(n) = w_j(n + h_j)` with `h_j >= 0`. So `P(n) >= w*(u) n`, and hence `|E ∩ [1,Y]| <= Y rho^peak_L + |Bad_L| + |X_L|`, where the exceptional set `X_L` is finite.
* **Lower bound (single-scale capacity).**
  * Take the sources `n <= Y/w*(u) - (L-1)/q` in bad classes; there are at least `Y rho^peak_L - 3|Bad_L|` of them.
  * Each must meet an edit `v = T^k(n) <= Y` with `k <= L-1`.
  * For fixed `(v, k, e)` the sources lie in an interval of length `e/q`, because the affine offset satisfies `0 <= h_k <= e/q` (THM-4478 §2).
  * So `v` serves at most `M*_L` sources, and `|E(G) ∩ [1,Y]| >= (Y rho^peak_L - 3|Bad_L|)/M*_L`. No slope cap is needed. ∎

**Theorem 2.** For every `L`, `rho^peak_L <= (q/2) sum_(q^e > 2^L) C(L,e) q^-e <= (q/2) 2^(-(1-H(c))L)`.
* The first inequality holds because `w_L <= (q/2) w*` on bad words.
* The second is Chernoff with tilt `z = qc/(1-c)`. It needs `z >= 1`, i.e. `c(q+1) >= 1`, i.e. `(q+1) ln 2 >= ln q`, which is true for every `q`.
* With Theorems 1 and 3, `-(1/L) log_2 eps_L(q) -> 1 - H(log_q 2)` for every odd `q`. ∎

**Theorem 3.** `2^((1-H)L) rho^peak_L = exp(-Theta(L^(1/3)))`.
* Change measure: under `Q`, letters are iid Bernoulli(`c`), and `2^-L = Q(u) 2^(-(1-H)L) lambda^(S_L)` with `lambda = (1-c)/c` and `S_j = e_j - jc`. Then `rho^peak_L = 2^(-(1-H)L) E_Q[lambda^(S_L) q^-M ; bad]`, where `M = max_(j<L) S_j`.
* A bad word keeps the centred walk in `(0, M+1)`. Blockwise Paley–Zygmund gives confinement cost `exp(-Omega(L/M^2))` against the discount `z^-M`, and optimising `M ~ L^(1/3)` gives the upper bound.
* An explicit confined event gives the lower bound.
* For `q = 3` (`lambda < 1`), `rho_L = 2^(-eta L)` up to `poly(L)`, so `rho^peak_L/rho_L = exp(-Theta(L^(1/3)))`. ∎

**Theorem 4 (modulo Mogul'skii's small-deviation lemma, CITED).** The sharp constant is `kappa_q = min_theta [theta ln z_q + pi^2 c(1-c)/(2 theta^2)]`, with `z_q = q min(1, c/(1-c))`. The fits give `2.095` (q=3, prefactor fixed) against `kappa_3 = 2.1076`.

## 2. Meaning

* **Drift sign does not matter for arbitrary edits.**
  * The exponent of the arbitrary-edit price is the entropy deficit `1 - H(log_q 2)` of the critical density, for every multiplier.
  * For `q = 3` it coincides with the undecided density's exponent (THM-4478).
  * For `q = 5` it does not: most undecided orbits escape upward and are caught cheaply at their peaks. Only the critical band, orbits staying near their starting height for `L` steps, has to be paid in full.
* **Height is what periodic modifications cannot see.**
  * In the periodic strategy cube (THM-4479) the price is polynomially sharp: `N_k <= delta_k <= |Bad_k|`.
  * Arbitrary edits gain `exp(-Theta(L^(1/3)))` for `q = 3`, and the whole exponent for `q = 5`. An edit at height `v` costs density weight `~1/v`.
  * This is the archimedean-versus-2-adic split again.
* **THM-4478's refinements.**
  * Its error term `O(sqrt(L log L))` improves to `Theta(L^(1/3))`.
  * Its question P1 is answered negatively for arbitrary edits.
  * For the pairing family `rho^peak_L/M*_L <= delta_L <= 2 rho_L` is proved, and the data track `rho^peak`.
* **Collatz.** Nothing here bears on Collatz itself. These are fixed-horizon modification prices.
