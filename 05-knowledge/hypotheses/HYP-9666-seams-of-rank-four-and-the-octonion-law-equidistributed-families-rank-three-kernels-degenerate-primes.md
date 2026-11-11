---
id: HYP-9666
title: "Seams of rank four and the octonion law: (A) is there a family of rooted tournaments whose rank-4 fusion states are equidistributed modulo odd primes, so that its counts carry the factor (1 - p^-4) sigma_3(p^v)/p^(3v) of E8; (B) rank-3 kernels (three visits) and the factor with sigma_2; (C) the law of a kernel at the primes of det G is a function of its Smith form; (D) module states at 3"
status: >
  OPEN ((A), (B) questions with a stated prediction and falsifier; (C) a conjecture that is elementary and not
  written out; (D) a PROVED consequence of THM-4756 recorded as a constraint). PRE-AUDIT DRAFT (2026-11-11).
  Nothing here is tested numerically. The only proved local law of a rank-4 form in the repository (THM-4916,
  two root changes) is NOT the uniform one.
session: opus-2026-11-11-S92
date: 2026-11-11
source: 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/README.md
related:
  - 01-canon/theorems/THM-5135-* ((1), (3), (4))
  - 01-canon/theorems/THM-4916-root-change-local-factor-identities-confinement-multipliers-product-family.md
  - 05-knowledge/results/session264_constructive_20261111/fusion/PROOF.md (fusion law, Theorem F1)
  - 01-canon/theorems/THM-5133-flip-and-spin-kernel-signature-sign-tournaments-of-twisted-algebras-and-the-pg32-extremal-theorems.md (visit matrices; PRE-AUDIT)
  - 01-canon/theorems/THM-4756-forward-step-polynomial-mod-3-is-h-times-eulerian.md
---

# HYP-9666 - seams of rank four and the octonion law

By THM-5135 (1), (3) a unimodular seam of rank 4 with uniform independent states has, at every prime, the local
factor `(1 - p^-4) sigma_3(p^v)/p^(3v)`, the local data of the lattice E8 (integral octonions, `240 sigma_3(n)`).
The join calculus has two unimodular rank-4 forms (THM-5135 (4)): the K-form of THM-4916 and the corner of
Session264's fusion law. Whether any family of tournaments realises the law is open.

**(A) An equidistributed rank-4 family.** Question: is there a family of rooted tournaments, closed under a
finite set of operations, whose fusion states `(h, H - h, U - h, d - h)` (Session264) are equidistributed modulo
every odd `p` on the cells of the family, as one-apex states are (THM-5128 (B2))? Prediction if yes: the fusion
counts `H(A * B)` have limiting class factor `1 + p^-3 - p^-4` at the class 0 modulo `p` and
`(1 - p^-4) sigma_3(p^v)/p^(3v)` at prime powers; the bias at 3 would be `2/81`, at 5 `4/625`. Falsifier: a family
with proved equidistribution and a different factor. Known obstruction: for the root-change states of THM-4916
the factor at 0 is `1 + (3p - 2)(p - 1)/(2 p^3)`, of order `3/(2p)`, because of the confinement cascade; and in
every family the states are odd in their first two coordinates sums (`H`, `h` odd), so the law cannot hold at 2.

**(B) Rank 3.** Kernels with three visits (THM-5133 (A)) give forms of rank 3. Prediction for uniform states:
factor `(1 - p^-3) sigma_2(p^v)/p^(2v)`; by THM-5135 (2) this is the law of six squares only at `p = 1 (mod 4)`.
No kernel of rank 3 has been examined.

**(C) Degenerate primes.** Conjecture: for `B(u, v) = u^T G v` with uniform states modulo `p^k`, the class law
depends only on the Smith normal form of `G` over `Z_p` (elementary divisors `p^(e_1), ..., p^(e_r)`), and is the
additive convolution of the laws of the scaled tables `p^(e_i) a_i b_i`. Checked instances (THM-5135 (4)): Session261's
kernel at 7 and 47 (elementary divisors `1, p`: the rank-1 law), `[[11, 6], [6, 3]]` at 3.

**(D) Module states at 3 (constraint, PROVED from THM-4756).** For unrooted module states `(H, Pi_2)` one has
`Pi_2(X) = (-1)^|X| H(X) (mod 3)` (THM-4756, as quoted in THM-5133 (H)). Hence the Gaussian kernel of THM-5133 (B),
`H(K[A, B]) = Pi_2(A) Pi_2(B) + H(A) H(B)`, is `((-1)^(|A| + |B|) + 1) H(A) H(B)` modulo 3: divisible by 3 whenever
`|A| + |B|` is odd. So module kernels do not have the uniform law at 3 and no series `sigma(n)/n` can be asserted
for them at that prime. The rooted one-apex states of F1 have no such constraint.
