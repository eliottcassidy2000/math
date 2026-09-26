---
id: THM-4477
title: "Cauchy-Schwarz lower bound for the price of provable descent, delta_L >= rho_L^2 / M2(L) >= 2^(-(0.1445+o(1))L); no distribution-only argument beats exponent 2(1-h); the square loss is forced by q = 3's AM-fairness g(2) = (1+q)/4 = 1"
status: >
  PROVED + INDEPENDENTLY AUDITED (Theorems CS, T and H); PROVED
  computer-assisted (Theorem M).
  Let W_L(v) = sum_(k<L) sum_(T^k n = v) 3^(a(n,k))/2^k be the hub weight of
  THM-4475 §3 and M2(L) = E_Haar[W_L^2] (W_L is a function of v mod 3^(L-1)).
  (CS) Every L-step-provable member of the pairing family has upper flip
  density >= rho_L^2 / M2(L).
  (T) The backward transfer operator Lf(v) = f(2v)/2 + (3/2)[v = 2 (mod 3)]
  f((2v-1)/3) satisfies
  <Lf, Lg> = <f,g> + (1/4)(<f o tau, g> + <f, g o tau>), with tau(u) = 3u+1.
  For qx+1 the diagonal coefficient is g_q(2) = (1+q)/4, which equals 1
  iff q = 3: THM-4470's AM-fairness.
  (M) ||Z_k||^2 <= 253.82 * 1.0312629^(k-14). This is a floating-point
  super-eigenvector certificate on classes mod 3^15 with a 10^-9 margin.
  Independently, an exact-rational certificate at 3^10 gives
  theta_10 <= 1.0491582. So delta_L >= 2^(-(0.1445+o(1))L) (exponent 0.1693
  from the exact certificate alone), improving THM-4475's 0.7737.
  (H) Every Holder exponent p != 2 gives a worse exponent, and any argument
  using only the law of W_L gives at best O(rho_L^2/L). So 2(1-h) = 0.1001
  is the moment-method limit. The square is exactly the Cauchy-Schwarz
  (AM-HM) loss, which abstracts away which orbit meets which hub. It is
  forced because the cascade's second root kappa = 2 of g = 1 holds iff
  E_fwd[w] = (1+q)/4 = 1, iff q = 3.
  OPEN: M2(L) = poly(L) (HYP-9139, which would give 0.1001 within this
  moment method). UPDATE 2026-09-26: the sharp exponent 1-h (HYP-9137) is
  PROVED by THM-4478, using a growth-band restriction and actual integer
  incidence, beyond the distribution-only scope of this theorem.
source: collatz-procgen-20260922 session, Cauchy-Schwarz lane (2026-09-25), from the coordinator's observation g_q(2) = (1+q)/4 and the owner's thesis that Cauchy-Schwarz / AM-GM marks where structure is abstracted away; audited and promoted by the session orchestrator 2026-09-25
depends_on:
  - 01-canon/theorems/THM-4475-price-of-provable-descent-tends-to-zero.md
  - 01-canon/theorems/THM-4470-collatz-pairing-ladder-am-fair-and-defect-blind.md
related:
  - 05-knowledge/hypotheses/HYP-9137-sharp-provability-price-exponent.md
  - 05-knowledge/hypotheses/HYP-9139-second-moment-of-hub-weight-polynomial.md
  - 05-knowledge/results/collatz_procgen_20260924_inverse_tree_mod192.md (Props. 9-10: g(s) = 2^-s + (3/2)^s/3)
script: 04-computation/experiments/procgen_cauchy_20260925_majorant.py
script_audit: 04-computation/experiments/procgen_cauchy_20260925_orchestrator_check.py
output: 05-knowledge/results/procgen_cauchy_20260925.out
output_audit: 05-knowledge/results/procgen_cauchy_20260925_orchestrator_check.out
script_sha256: 4d604d93f7a843b9b73bf5f911af239128fffa7d26f7e04cd3a72a23928bcb93
script_audit_sha256: d6ba66e8e267899d8066b359a4b3d1fc8e9e4a5859e03cbedc252ba032887723
output_sha256: 09a6ea552f9a22f5ec0ca0ccfd4cf3be25fe93055bd219eb09d3341bec56b69c
output_audit_sha256: 8eff1060daa8580b37aef6a46b31a0a747a3eaca274a74f6bc644bcb16ed40e2
hash_basis: raw LF bytes
audit: >
  The orchestrator re-derived the following by hand:
  * Theorem CS: THM-4475 §3 plus one Cauchy-Schwarz step, with Lemma 1 (the
    depth-k tree depends on v mod 3^k) making log-averages Haar averages;
  * Theorem T's inner-product identity: the Haar Jacobian of the odd
    branch is 1/3, and 2v = 3w + 1 under w = (2v-1)/3;
  * the majorant recursion Phi.
  Independent code confirms:
  * the exact M2(8) = 501975/4096;
  * the converged Perron roots of Phi (theta_1 = 1.34307,
    theta_10 = 1.04911; the lane certified 1.04981);
  * a rigorous exact-rational certificate theta_10 <= 1.0491582 on the
    non-zero classes. The multiple-of-3 classes evolve exactly with ratio
    1/4.
  An earlier orchestrator attempt that floored those classes produced a
  weaker 1.0540; that was an artifact of the audit code, not a fault of
  the lane. The lane's pipeline was re-run: identical output apart from
  timing lines. The r = 15 certificate is floating point with a 10^-9
  inflation (per-entry rounding error below 10^-14), not interval
  arithmetic. The exact certificate at r = 10 independently proves an
  exponent of at most 0.1694.
---

# THM-4477 -- the Cauchy–Schwarz price bound and the AM-fair criticality of q = 3

**PROVED + INDEPENDENTLY AUDITED.** Full note: [procgen_cauchy_20260925_cauchy_schwarz_criticality](../../05-knowledge/results/procgen_cauchy_20260925_cauchy_schwarz_criticality.md).

## 1. The bound

* **The setting (THM-4475).** Every `L`-step-provable member of the pairing family must flip a pair on each Collatz-undecided orbit (`n in Bad_L`) within `L` steps. A flip at `v` serves logarithmic weight at most `W_L(v)/v`.
* **Cauchy–Schwarz.** It gives `rho_L <= (delta_bar · M2(L))^(1/2)`.
* **The weight's dependence.** The depth-`k` backward tree of `v`, with its weights, is a function of `v mod 3^k`. So logarithmic averages over integers are Haar averages.

## 2. Where `g(2) = 1` enters, and where it does not

* **Transfer operator.** `Z_k := sum_(T^k n = v) 3^a/2^k = L^k 1`, the density of the forward walk `X -> X/2` or `(3X+1)/2` on `Z_3` from Haar measure.
* **The identity.** `<Lf, Lg> = <f,g> + (1/4)(<f∘tau, g> + <f, g∘tau>)`, with `tau(u) = 3u+1`. The diagonal coefficient `1/4 + 3/4 = g(2)` is neutral exactly because `q = 3`; for general `q` it is `(1+q)/4`.
* **Why the second moment still grows.** `||Z_(k+1)||^2 = ||Z_k||^2 + gamma_k/2`. The growth sits entirely in the cross term `gamma_k = <Z_k∘tau, Z_k>`. It arises because the two children of a node, `3x+1` and `x`, are deterministically linked rather than independent.
* **Data.** Exact values: `||Z_15||^2 = 6.913` (increments about 0.357; independent siblings would give 0.5), and `M2(16) = 762.55`.

## 3. The certified growth rate (Theorem M)

* **The majorant.** Class masses `u_k(c) = ∫_c Z_k^2` for `c mod 3^r` satisfy `u_(k+1) <= Phi(u_k)`, where
  `Phi(u)(c) = u(2c)/4 + [c = 2 (mod 3)] ( (3/4) U(Ec) + (sqrt3/2) sqrt(u(2c) U(Ec)) )`,
  and `U` sums the three lifts of a class mod `3^(r-1)`.
* **The certificate.** A super-eigenvector with `Phi(e) <= theta e` bounds the growth.
* **Values.**
  * `theta_15 = 1.0312629`: floating point, with margin.
  * `theta_10 <= 1.0491582`: exact rationals, orchestrator.
* **Result.** `delta_L >= 2^(-(2eta + log2 theta)L - O(log L))`, with `eta = 1 - h = 0.050044`. That is `0.1445` (via `theta_15`) and `0.1693` (via `theta_10`).

## 4. The moment-method limit (Theorem H)

* **Hölder.** For every `p != 2`, `E[Z_k^p] >= g(p)^k` makes the Hölder exponent strictly worse than `2 eta`.
* **Distribution-only arguments.** An adversary placing flips on the heaviest classes needs only `O(rho_L^2/L)` density. So no argument using only the law of `W_L` beats `2 eta`.
* **Why the square.** The tail exponent of the backward cascade is the non-trivial root `kappa` of `g_q(s) = 2^-s (1 + q^(s-1)) = 1`. The best moment loss is `rho -> rho^(kappa/(kappa-1))`, and `kappa = 2` iff `E_fwd[w] = (1+q)/4 = 1` iff `q = 3`.
  * The equation is the same as THM-4470's pair-sum AM-fairness.
  * The quadratic loss is precisely the owner's point that Cauchy–Schwarz (here AM–HM on the hub weight along bad orbits) abstracts away structure: which orbit meets which hub.
  * The data show the abstraction is costly. Bad orbits carry about `2–2.6` times the average hub weight, where Cauchy–Schwarz must allow `1/rho_L`, which is about `31` at `L = 16`.

## 5. Beyond moments (for HYP-9137)

**Current resolution, 2026-09-26:** [THM-4478](THM-4478-critical-tube-affine-capacity-sharp-provability-price.md)
proves HYP-9137 via critical growth bands and affine integer intervals.
The historical harmonic route below remains interesting, but is not needed
for the sharp exponent. HYP-9139 itself remains OPEN.

* **The harmonic bound (lane Theorem P, PROVED).** `delta_bar >= H_L = E[1_Bad / max_(j<L) W_L(X_j)] >= max(rho_L^2/M2(L), rho_L/max W_L)`.
* **Data.** Numerically `H_L ≈ 3.2 rho_L/M2(L)` for `L <= 16`, i.e. only a polynomial loss against `rho_L`.
* **What HYP-9137 needs.** It follows if typical bad orbits are hub-biased by at most `2^(o(L))`. That is quenched independence of 2-adic badness and 3-adic hub weight. It holds in the i.i.d. mean-field model and is OPEN for the real tree.

## 6. Controls

* **SHEET.** `Z^-_k(v) = Z_k(-v)`, so every bound here is identical on the `3n−1` sheet.
* **DRIFT.** For `5x+1`, `||Z^(5)_k||^2 >= (3/2)^k` because `g_5(2) = 3/2`. Both the Cauchy–Schwarz bound and the harmonic bound tend to 0, although the observed `5x+1` price stays near 0.64 of an undecided density of about 0.176. That floor is a multi-flip phenomenon.
