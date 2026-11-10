---
id: THM-5127
title: "The heat flow of the Hamiltonian-path count and the mod-4 law of a vertex cube. (A) The W-polynomial is the noise semigroup of the arc cube applied to H: (N_rho H)(T) = sum_k a_k(T) p^k q^(n-1-k) and W_T(r) = (2r)^(n-1) (N_{1/(2r)} H)(T) for |r| >= 1/2; in the coordinate w = p(1-p) the noised count is sum_j gamma_j(T) w^j with integer gamma_j(T), gamma_0(T) = H(T) and sum_j gamma_j(T) 4^(-j) = n!/2^(n-1). (B) Tilted ensemble: if each pair i < j is oriented i -> j independently with probability p, then E_p[H_n] = sum_j gamma_{n,j} (p(1-p))^j with gamma_{n,j} the gamma-vector of the Eulerian polynomial; for n >= 3 it decreases strictly and continuously from n!/2^(n-1) at p = 1/2 to 1 at p = 1. (C) For a tournament R on m vertices the pair term mu_R(i,j) of the cut quadratic is even and >= 2 (KNOWN); X -> (H(R + v_X) - H(R))/2 mod 2 is a quadratic form over F_2 with polar form beta_R = mu_R/2 mod 2; the mod-4 imbalance of H over the 2^m out-sets is 0 or +-2^(m - rank(beta_R)/2); hence the mod-4 discrepancy of H over labelled tournaments of order m+1 is at most the mean of 2^(-rank(beta_R)/2) over order m. (C4) |E_X e(theta H(R + v_X))|^2 <= E_d prod_{i in U} |cos(pi theta sum_{j in V'} mu_R(i,j) d_j)|. FINITE-EXACT: beta_R on all classes of orders 2-9 and all labelled tournaments of orders 2-7 (rank statistics within 1.2 points of a uniformly random alternating matrix at orders 8, 9, and within 0.08 points at order 10 in a pre-registered run); class counts at orders 9, 10; the full deletion window [3h, 5h) for 237 <= h <= 1023 and [3h, 9h) for 365 <= h <= 1023 with bases of order <= 10. Chan-Pak OPEN"
status: >
  PROVED: (A), (B), (C) (2)-(4), (C4), and the implications Proposition 5.1 and Proposition 6.1 of the source note.
  KNOWN, credited:
  - forward-edge numbers, palindromy and the W-polynomial (THM-059, THM-062); the transitive case
    a_k(TT_n) = A(n,k), i.e. THM-059's master polynomial F_(n-1);
  - the algebra of (A) for every n: THM-4757 section 5 (sum_L (n-L)! N_L t^L = (1-t)^(n-1) F(T,(1+t)/(1-t)),
    [r^(n-1-L)] W(r) = 2^(-L)(n-L)! N_L) with THM-4761 (1) (N_L/(n)_L is the degree-L Walsh chaos); THM-071 covers
    n = 5, 7 only. (A) is their probabilistic reading; the gamma-coordinates are the change of variable
    rho^2 = 1 - 4w;
  - gamma-nonnegativity of the Eulerian polynomial (Foata-Schutzenberger 1970, CITED);
  - the cut quadratic of a one-vertex ear and its submodularity (THM-4097 section 2, THM-4114); its form (4.2) and the
    evenness of its coefficients (THM-4097 eq. (6a), THM-4118 eqs. (4), (9)); strict positivity
    mu_R(i,j) >= 2 H(R - {i,j}) (Session 136 positive_cut_geometry.md eq. (3); THM-5122). So (C) (1) is KNOWN; the
    insertion argument of the source note is a second proof of mu >= 2 that does not use the odd-cycle formula;
  - the character sum of a quadratic form over F_2 (Dickson, CITED);
  - the window induction as a route (klein, HYP-9453 (2), for the near-join); the forced lower boundary
    n >= (3h+1)/2 (Session214, HYP-9469; THM-4876 (2), (7));
  - the Gaussian regime at p = 1/2 (THM-4761, CONDITIONAL on Moran's CITED theorem; THM-589
    metagraph-H-variance-closed-form-nud-succession);
  - the Bourgain-Kontorovich analogy for word semigroups (Session205, opus S82, opus S90).
  NEW here: the noise/tilt reading and the gamma form, positivity reading and strict monotonicity of Theorem B; the F_2
  quadratic form of a vertex cube and its rank bound for H mod 4 ((C) (2)-(4)); (C4) (opus S90's note has an exact
  ear-cube Weyl identity found independently the same day); the censuses listed under FINITE-EXACT.
  FINITE-EXACT: sections 4.1, 5 and 6 of the source note (programs and outputs in the packet); pre-registered P1-P2
  (order-10 rank distribution) CONFIRMED.
  NUMERICAL: the bridging tilt p_n and its asymptotic.
  NOT CLAIMED: the Chan-Pak conjecture, positive density of the spectrum, or any part of HYP-9659.
  AUDIT: PASS WITH CORRECTIONS (05-knowledge/results/tilted_ensemble_local_limit_20261110/audit/AUDIT.md; 5 MAJOR:
  M1 (C) (1) and (4.2) are prior canon, M2 the all-n algebra of (A) is THM-4757/4761, M3 the hypothesis (B) did not
  imply its stated consequence, M4 Proposition 6.1's hypothesis mis-worded, M5 a statement misattributed to THM-4876;
  22 minor; all applied; every theorem and every number reproduced with independent code; MISTAKE-753).
session: mac-mini-2026-11-10-ensemble
source: 05-knowledge/results/tilted_ensemble_local_limit_20261110/README.md
related:
  - HYP-9659
  - THM-059
  - THM-062
  - THM-071
  - THM-589 (metagraph-H-variance-closed-form-nud-succession)
    - THM-4097
  - THM-4114
  - THM-4118
  - THM-4757
  - THM-4761
  - THM-4756
  - THM-5122

  - THM-4876
  - THM-4885
  - THM-5124
  - HYP-9453
  - THM-1370 (h-spectrum-omits-7-21-all-n)
verification:
  - 05-knowledge/results/tilted_ensemble_local_limit_20261110/verify_cut_mod4.py
  - 05-knowledge/results/tilted_ensemble_local_limit_20261110/mallows.py
  - 05-knowledge/results/tilted_ensemble_local_limit_20261110/betarank.c
  - 05-knowledge/results/tilted_ensemble_local_limit_20261110/pairs.c
  - 05-knowledge/results/tilted_ensemble_local_limit_20261110/analyze_window.py
  - 05-knowledge/results/tilted_ensemble_local_limit_20261110/hist.c
  - 05-knowledge/results/tilted_ensemble_local_limit_20261110/analyze_hist.py
  - 05-knowledge/results/tilted_ensemble_local_limit_20261110/out/
---

# THM-5127 - the heat flow of H and the mod-4 law of a vertex cube

Full statements, proofs and tables: [source note](../../05-knowledge/results/tilted_ensemble_local_limit_20261110/README.md).
The Chan-Pak conjecture (every odd number except 7 and 21 is a Hamiltonian path count) is **OPEN**.

## Notation

- `H(T)`: number of Hamiltonian paths of a tournament `T` on `n` vertices. `a_k(T)`: number of permutations with
  exactly `k` forward steps (THM-062). `W_T(r) = sum_k a_k(T) (r + 1/2)^k (r - 1/2)^(n-1-k)`.
- `N_rho`: reverse every arc independently with probability `q = (1 - rho)/2`; `p = 1 - q`.
- `P_p`: each pair `i < j` oriented `i -> j` with probability `p`, independently.
- `R + v_X`: add a vertex `v` with out-neighbourhood `X` to `R`.

## (A) Noise identity (PROVED)

1. `(N_rho H)(T) = sum_k a_k(T) p^k q^(n-1-k)`.
2. `W_T(r) = (2r)^(n-1) (N_{1/(2r)} H)(T)` for `|r| >= 1/2`.
3. `sum_k a_k(T) x^k y^(n-1-k) = sum_j gamma_j(T) (xy)^j (x+y)^(n-1-2j)` with unique integers `gamma_j(T)`;
   `(N_rho H)(T) = sum_j gamma_j(T) w^j`, `w = pq`; `gamma_0(T) = H(T)`; `sum_j gamma_j(T) 4^(-j) = n!/2^(n-1)`.

Mechanism: a permutation uses `n-1` distinct pairs, so its steps survive independently. Equality boundary: `|r| < 1/2`
is not a probability.

## (B) Tilted mean (PROVED)

`E_p[H_n] = sum_k A(n,k) p^(n-1-k) q^k = sum_j gamma_{n,j} (p(1-p))^j`, with `gamma_{n,j} >= 0` the gamma-vector of
the Eulerian polynomial (`gamma_{n,0} = 1`, `gamma_{n,1} = 2^n - 2n`). For `n >= 3` the mean is strictly decreasing
and continuous in `p` on `[1/2, 1]`, from `n!/2^(n-1)` to `1`. Checked exactly on all labelled tournaments for
`n <= 7`: gamma-vectors `(1,2)`, `(1,8)`, `(1,22,16)`, `(1,52,136)`, `(1,114,720,272)`.

## (C) Mod-4 law (PROVED)

Cut quadratic (KNOWN): `H(R + v_X) = H(R) + sum_{i in X} a_i - sum_{{i,j} subset X} mu_R(i,j)`, where
`a_i = H(R + v_{i}) - H(R)` (THM-4118 eq. (4)) and `mu_R(i,j)` counts the ordered two-path covers of `R` with one path ending at one of
`i, j` and the other starting at the other. Equivalently `mu_R(i,j)` is the number of Hamiltonian paths of `R`
through the arc between `i` and `j` plus the number of Hamiltonian paths of the arc-reversed tournament through the
reversed arc.

1. (KNOWN: THM-4097 eq. (6a), THM-4118 eq. (9); Session 136; THM-5122.) `mu_R(i,j)` is even and at least
   `2 H(R - {i,j}) >= 2`.
2. `q_R(X) = (H(R + v_X) - H(R))/2 mod 2` is a quadratic form on `F_2^{V(R)}` with polar form
   `beta_R(i,j) = mu_R(i,j)/2 mod 2`.
3. `#{q_R = 0} - #{q_R = 1}` is `0` or `+-2^(m - rank(beta_R)/2)`; it is `0` iff `q_R` is not identically zero on the
   radical of `beta_R`.
4. For `T` uniform on `m+1` labelled vertices,
   `|P(H = 1 mod 4) - P(H = 3 mod 4)| <= E_R[2^(-rank(beta_R)/2)]`, `R` uniform on `m` labelled vertices.

(C4) For disjoint `U, V' subset V(R)`, real `theta`, and `d_j` independent in `{-1, 0, 1}` with probabilities
`1/4, 1/2, 1/4`: `|E_X e(theta H(R + v_X))|^2 <= E_d prod_{i in U} |cos(pi theta sum_{j in V'} mu_R(i,j) d_j)|`.

Failure boundary: the count is an exact quadratic only in the arcs at one new vertex.

## FINITE-EXACT

- `beta_R` on all 198,947 tournament classes of orders 2-9 and all labelled tournaments of orders 2-7: no parity or prediction failure; rank distribution at order 8: `2, 205, 3795, 2878` classes of rank `2, 4, 6, 8`; at order 9:
  `8, 436, 32535, 158557`; at order 10 (pre-registered, ranks only, 9,733,056 classes): `7, 668, 210562, 5444657,
  4077162` of rank `2, 4, 6, 8, 10`. Uniformly random alternating matrices: `0.4205 / 0.5585 / 0.0210` (order 8, ranks
  `8 / 6 / lower`) against `0.4183 / 0.5516 / 0.0301`; `0.8394 / 0.1593 / 0.0013` (order 9) against `0.8278 / 0.1699 / 0.0023`; `0.4197 / 0.5591 / 0.0213` (order 10) against
  `0.4189 / 0.5594 / 0.0217`. (Census fractions are over classes.)
- Labelled mean imbalance at `m = 2..7`: `0.5, 0.5, 0.40625, 0.25, 0.12793, 0.08835`; so the mod-4 discrepancy at
  order 8 is at most `0.0884`. Exact discrepancies at orders 3-7: `1/2, 1/2, 13/32, 1/4, 853/8192`; at orders 8, 9, 10 (audit, class
  weights): `1361/32768`, `63103/8388608`, `650471/268435456` (0.0415, 0.00752, 0.00242).
- Class counts: at order 10 every odd `h <= 14683` except 7, 21 occurs; in the window `[3963, 10192]` (class mean
  +- sd) every odd `h` has at least 1,504 classes; residue classes mod 3, 4, 5, 7, 8 are balanced to 1.2% (the only effect seen is an excess of multiples of 3);
  the local variance of the class counts is about twice Poisson.
- Deletion pairs with bases of order `<= 10` and all out-sets: for every odd `h` in `[237, 1023]` every odd `n` in
  `[3h, 5h)` is a one-vertex extension of a realiser of `h`; for `h` in `[365, 1023]` the same holds on `[3h, 9h)`.
  Below 237 exactly 48 values of `h` have an incomplete window `[3h, 5h)`. No single base of order `<= 10` covers
  more than 44% of `[3h, 5h)` for any `h >= 237` (audit).

## Consequences that are proved

- Proposition 5.1: the tilted local limit law HYP-9659 (A) implies that all sufficiently large odd numbers are
  Hamiltonian path counts.
- Proposition 6.1: the window law HYP-9659 (C) (bases of any order, every odd `h >= 237` that is a count) implies the
  Chan-Pak conjecture.
