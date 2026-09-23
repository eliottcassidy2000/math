---
id: THM-4467
title: "Uniform Polya-capacity gap for AMM 12592: every exactly fair critical-run extractor has C >= 11/8"
status: >
  PROVED (classical inputs CITED: Polya 1928, Fatou 1906) + certificates
  FINITE-EXACT + INDEPENDENTLY AUDITED. No deterministic exactly fair
  unknown-bias coin extractor has a pathwise deadline T(n) <= (1+gamma)n + D
  for any constant D when Lambda(gamma) = log R(W_gamma, 0) > 0. The
  certificates give Lambda(7/20) > 0, by pure interval arithmetic, and
  Lambda(3/8) > 0, by a subordination polynomial whose boundary is checked on
  2^21 FFT nodes under an explicit a-priori floating-point error allowance.
  Hence the uniform constant satisfies C* >= 27/20 (interval-rigorous) and
  C* >= 11/8 (with the FFT caveat). The method's numerically validated reach
  is C* >= 1.3775246. This is the first uniform positive gap. It generalizes
  THM-3342, and together with the balanced-block constructions it gives
  C* in [11/8, 1 + log_5(phi^2)].
source: collatz-procgen-20260923 AMM lane (mac-mini), audited and promoted by the session orchestrator 2026-09-23
depends_on:
  - THM-2966-spine-normal-form-for-critical-run-fair-extractors
related:
  - THM-3342-sublinear-deadline-excess-is-impossible-for-fair-critical-run-extractors
  - THM-3009-archimedean-floor-for-balanced-block-extractors
  - THM-3027-capacity-threshold-is-log-sqrt5-phi
script: 04-computation/experiments/amm12592_procgen_20260923_polya_capacity.py
output: 05-knowledge/results/amm12592_procgen_20260923.out
script_sha256: f01a4fd77cfa606e83818c4b5aaff1057e9968b85872d6d4e97b8e0271d1fd3d
output_sha256: f3c4528a51adbb62503fa4982ff525b0a756d0aecd0c3e018261ee2217a92048
certificate: 05-knowledge/results/amm12592_procgen_20260923_subordination_certificate.json
certificate_sha256: 5eccbc0c85794a35297aceb75bb2b3df51f19a71bc9035983e283bc7cf5a5fc5
hash_basis: raw LF bytes
audit: >
  Orchestrator re-derived continuation, gluing (Lemma A1 geometry), the
  Catalan quotient, the capacity form of Polya's theorem and the Fatou-Gauss
  endgame by hand. Independent code
  (04-computation/experiments/amm12592_procgen_20260923_orchestrator_audit.py)
  re-verified the 7/20 two-disk certificate (containment and closed-form
  conformal radius, with its c = 0 sanity value). The stored 3/8 certificate
  was re-run with --verify-only.
---

# THM-4467 -- a uniform capacity gap for time-limited fair extraction

**PROVED + FINITE-EXACT certificates + INDEPENDENTLY AUDITED.** The full
exposition, tables and failure anatomy are in
[amm12592_procgen_20260923_uniform_frontier](../../05-knowledge/results/amm12592_procgen_20260923_uniform_frontier.md) §3.

## 1. Inheritance

* **Closest proved mechanism:** THM-3342. There, Carlson's theorem gives
  per-extractor nonattainment of `n + o(n)`, with no uniform gap.
* **Canonical hostile example:** the balanced-block class. There THM-3009 and
  Long's 2026-09-20 draft pin the constant `C_* = 1 + log_5(phi^2) = 1.59799`.
* **Least-used sidecar:** the *two-point* integrality of the spine functions,
  at `p = 0` and at `p = 1`. Folding by `p -> 1-p` turns it into one-point
  integrality on a quotient domain.

## 2. Statement

Use the spine normal form (THM-2966): `F(p) = sum_m p^m q W_m(p)`, `G`
likewise, with `F(p) + G(1-p) = 1/2`, and `d_m = T(m) - m - 1`. Define

```text
S(p) = |p| + |1-p|,
U_g = { min(|p|,|1-p|) S(p)^g < 1 },
W_g = { p(1-p) : p in U_g },
Lambda(g) = log R(W_g, 0).
```

**Theorem.** If `Lambda(gamma) > 0`, no exactly fair extractor satisfies
`T(n) <= (1+gamma) n + D` for all `n`. Since `Lambda` is nonincreasing, this
gives `C* >= 1 + gamma`.

**Certified instances.**
* `Lambda(7/20) >= log R(w(D(13/200, 39/50) ∪ D(187/200, 39/50)), 0) = 0.0059354 > 0`,
  by interval arithmetic. So `C* >= 27/20`.
* `Lambda(3/8) >= log |P'(0)| = log 1.00071974 > 0`, by a subordination
  polynomial `P` with `P(∂D) ⊂ W_(3/8)`. The boundary is checked on `2^21`
  nodes with explicit perturbation bounds. So `C* >= 11/8`.

## 3. Proof (summary)

1. `|W_m(p)| <= S(p)^(d_m)` (since `0 <= w_(m,k) <= C(d_m,k)`). So `F` converges
   on `Omega_0 = {|p| S^gamma < 1}`, and `G` on the mirror domain.
2. The fairness identity holds on the connected overlap. So `F` and `G`
   continue to the simply connected `U_gamma`.
3. `Delta = F - G` is `(p -> 1-p)`-even and `Sigma = F + G - 1/2` is odd.
   Through `z(w) = sum Cat(n-1) w^n`, both `Delta` and `Sigma/(1-2z)` descend to
   integer power series (up to the constant `-1/2`) that are analytic on
   `W_gamma`.
4. Pólya (1928): an integer power series analytic on a domain of conformal
   radius `> 1` at `0` is rational. So `F` and `G` are rational.
5. Fatou: `F = A/B` with `B(0) = 1`, and `G = C/D` with `D(0) = 1`. Then
   `C(1-p)/D(1-p) = (B - 2A)/(2B)` in lowest terms. Primitivity forces
   `D(1-p) = +-B`, so `B in 2Z[p]`, which contradicts `B(0) = 1`. ∎

## 4. Scope, boundary and what it supersedes

* **Scope.** The theorem holds for every deterministic exactly fair extractor,
  with no balancing or symmetry assumption.
* **Boundary.**
  * `Lambda(C_* - 1) < 0` and `Lambda(1) < 0`, where constructions exist.
    The obstruction correctly does not fire there.
  * The method's reach is `gamma_1 = 0.3775246`, from two independent Laplace
    solvers.
* **What it supersedes.** It replaces the ledger status "no positive uniform
  lower gap". THM-3342 is its `gamma -> 0` shadow.
* **Open.**
  * Whether `C* < 3/2` (HYP-9129).
  * Whether super-blocks prove `C* < C_*` (HYP-9128). There is FINITE-EXACT
    block-level evidence at `53/34`, `83/53` and `157/100`, and numerical
    asymptotics of about `1.570`.

## 5. Update (2026-09-23): extended to gamma = 0.377, and the upper side

* **Extension.** A second subordination certificate at `gamma = 377/1000`
  gives `C* >= 1.377`, with the same FFT caveat as the `3/8` certificate:
  * the certificate polynomial has `b_1 = 1.0001609`;
  * the boundary is checked on `2^23` nodes, with a maximum certified bound
    of `0.99999023 < 1`.

  The orchestrator re-verified it with `--verify-only`. See
  `05-knowledge/results/amm12592_procgen_20260923_hyp9128_proof.md` §8.
* **Upper side.** THM-4468 proves `C* <= 159/100`.
* **Proved window:** `1.377 <= C* <= 1.59`.
* **The owner's structured-majorant question** (same note, §9) has a
  negative answer. No termwise majorant can beat `S^(d_m)` exponentially,
  because fairness-kernel moves saturate single levels. The loss is
  cross-level.
