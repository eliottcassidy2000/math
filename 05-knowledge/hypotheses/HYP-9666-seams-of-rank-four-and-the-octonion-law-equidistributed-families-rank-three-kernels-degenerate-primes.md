---
id: HYP-9666
title: "Seams of rank four and the octonion law: (A) is there a family of rooted tournaments whose left and right rank-4 fusion states are equidistributed modulo odd prime powers, so that its counts carry at the odd primes the factor (1 - p^-4) sigma_3(p^v)/p^(3v) of E8 (at 2 the law fails in every family); (B) rank-3 kernels (three visits) and the factor with sigma_2, for module states only at the primes p >= 5 not dividing det N; (C) the law of a kernel at the primes of det G is a function of its Smith form: PROVED by the audit, now THM-5135 (4)(b); (D) module states at 3: every module kernel is a multiple of H(A) H(B) modulo 3"
status: >
  OPEN ((A), (B): questions with a stated prediction; (A) says what would settle it negatively; (B) concerns, for
  module states, only the primes p >= 5 not dividing det N. (C) is no longer open: PROVED by the audit, elementary,
  now THM-5135 (4)(b). (D): a PROVED consequence of THM-4756, recorded as a constraint and extended by the audit
  to every module kernel); audited 2026-11-11 with THM-5134 / THM-5135 (PASS WITH CORRECTIONS, MISTAKE-763).
  Report: 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/AUDIT_REPORT.md.
  The session ran no computation for (A), (B), (D). Apart from the datum of HYP-9660 (D) quoted in (A), every
  check in this file is the audit's (receipts under scripts; the last one is the verification of the corrections,
  with independent code). FINITE-EXACT: the fusion law (3) of Session264 on the 5625 pairs of labelled pointed
  tournaments (root = vertex 0) of orders at most 4; the parities of the fusion states of every labelled pointed
  tournament of order at most 7, and the absence of affine relations among them modulo 3, 5, 7 at each of the
  orders 4 to 7; kernels of rank 3 at orders 4 and 5; the collapse modulo 3 of every two-port kernel of order
  3, 4, 5 to a multiple of H(A) H(B), for modules A, B of order at most 4 with K[A, B] of order at most 10.
  NUMERICAL (exploratory): the class-0 factors of H(A * B) modulo 3, 5, 7 over all labelled pointed tournaments
  of one order, orders 4 to 7; they do not decide (A).
  The only proved local law attached to a rank-4 form in the repository is that of the K-form on the root-change
  words of THM-4916 (at most two root changes, letter-count model; one word, hence one free state, per value; its
  by-value limit is OPEN, HYP-9506 (A), HYP-9660 (A)). It is NOT the uniform one: its factor at the class 0 is
  neither the value 1 + p^-3 - p^-4 of two independent uniform states nor the value 1 of one uniform state under
  a non-zero linear functional.
session: opus-2026-11-11-S92
date: 2026-11-11
source: 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/README.md
scripts:
  - 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/fusion_states.py (+ .out; audit, for (A), the fusion law (3), the parities, the residue vectors and the affine relations of the fusion states)
  - 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/fusion_law_small.py (+ .out; audit, for (A), exploratory class factors of H(A * B))
  - 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/rank_laws.py (+ .out; audit, the laws of rank 2, 3, 4 for uniform states, the share 17/32 at 2, THM-4916 against the uniform values)
  - 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/kernel_mod3.py (+ .out; audit, for (B) and (D), every two-port kernel of order 3, 4, 5 modulo 3)
  - 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/hypotheses/pi2_mod3.py (+ .out; audit, for (D), the congruence of THM-4756 and the Gaussian kernel)
  - 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/proof-5135/a4_kernels.py (+ .out; audit, section D the two Gram matrices of rank 4, section F THM-4916 against the uniform values)
  - 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/proof-5135/a6_identity_rank3.py (+ .out; audit, section V, for (B), kernels of rank 3)
  - 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/verify-HYP-9666/verify_hyp9666.py (+ .out; verification of the corrections, independent code, every number of (A), (B), (D) except the datum of HYP-9660 (D))
related:
  - 01-canon/theorems/THM-5135-* ((1), (2), (3), (4)(a); its (4)(b) is the Smith-form law, formerly (C) here; (4)(c) the two forms of rank 4)
  - 01-canon/theorems/THM-5134-* (F1, the seam of rank 2 with equidistributed states; audited together)
  - 01-canon/theorems/THM-5128-f1-by-value-local-densities-first-order-class-law.md ((A) class law of a gluing; (B2) cells)
  - 01-canon/theorems/THM-4916-root-change-local-factor-identities-confinement-multipliers-product-family.md (letter-count model; the K-form is (1d))
  - 05-knowledge/results/root_change_local_factors_20261107/PROOF.md (source of THM-4916; sections 1 (1d), 2.1, 3 item 2, one free state)
  - 05-knowledge/hypotheses/HYP-9506-* ((A) census limit of the two-root-change family, OPEN)
  - 05-knowledge/hypotheses/HYP-9660-* ((A) by-value limit of two root changes, conjecture; (D) the two-apex grammar G(2,0))
  - 05-knowledge/results/session264_constructive_20261111/fusion/PROOF.md (fusion law, Theorem F1)
  - 01-canon/theorems/THM-5133-flip-and-spin-kernel-signature-sign-tournaments-of-twisted-algebras-and-the-pg32-extremal-theorems.md (visit matrices (A), Gaussian kernel (B), the congruence quoted in (H); audited 2026-11-11, kind-pasteur, PASS WITH CORRECTIONS, MISTAKE-762; the statements used here are unchanged by that audit)
  - 01-canon/theorems/THM-4756-forward-step-polynomial-mod-3-is-h-times-eulerian.md
  - 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/AUDIT_REPORT.md (independent audit; findings F010, F026, F027, F032, F039, F041, F055, F074, F103 to F107, F110 mention this file)
---

# HYP-9666 - seams of rank four and the octonion law

By THM-5135 (1), (3) a unimodular seam of rank 4 with uniform independent states has, at every prime, the local
factor `(1 - p^-4) sigma_3(p^v)/p^(3v)`, the local data of the lattice E8 (integral octonions, `240 sigma_3(n)`).
Two forms of rank 4 occur in the join calculus (THM-5135 (4)(c)), both of determinant `-1` in the coordinates in
which they are displayed there: the K-form (Session199; recorded in THM-4916 (1d); which lattice its states fill
was not examined by the audit), and the corner of Session264's fusion law. The corner is unimodular in the
coordinates `(h, H - h, U - h, d - h)`, but the states do not fill the last two: they satisfy `U = d (mod 2)`, a
sublattice of index 2, which their pairs `(U - h, d - h)` span (`(1, -1)` for the arc rooted at its tail, `(2, 0)`
for the cyclic triangle). In the coordinates `(a_o, a_e)` of that sublattice (the path masses of (A)) the corner
has determinant `-4`, so it is unimodular at the odd primes only. For uniform states both forms have the rank-4
law (audit, direct enumeration modulo 2, 3, 4, 5, 9; audit/proof-5135/a4_kernels.out [D]), so a failure of the law
is a property of the family of states. Whether any family of tournaments realises the law is open.

Receipts named `audit/...` are in `05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/`.

**(A) An equidistributed rank-4 family.** The fusion law (Session264, Theorem F1, eq. (3); primes for `B`),
`H(A * B) = h h' + (H - h) h' + h (H' - h') + (U - h)(V' - h') + (d - h)(e' - h')`, pairs the left state of `A`
with the right state of `B`. Question: is there a family of rooted tournaments, closed under a finite set of
operations, whose left states `(h, H - h, U - h, d - h)` and right states `(h, H - h, V - h, e - h)` are
equidistributed modulo every odd `p` and every power of it on the cells of the family, as one-apex states are
(THM-5128 (B2))? Both sides are needed: by THM-5128 (A), if only one side is uniform modulo `p`, the excess of the
class 0 over the mean of the other classes is the mass of the other side at the zero state. Prediction if yes:
the fusion counts `H(A * B)` have limiting class factor `1 + p^-3 - p^-4` at the class 0 modulo `p` and
`(1 - p^-4) sigma_3(p^v)/p^(3v)` at prime powers; the bias at 3 would be `2/81`, at 5 `4/625`.

What would settle it negatively: a congruence satisfied modulo an odd prime by the fusion states of every pointed
tournament, as (D) at 3 for module states. None exists among affine relations modulo 3, 5, 7: at each of the
orders 4, 5, 6, 7 the left states of all labelled pointed tournaments satisfy no affine relation modulo 3, 5 or 7,
nor do the right states, and at order 7 they take all 81 residue vectors modulo 3, 591 of the 625 modulo 5 and
1384 of the 2401 modulo 7 (FINITE-EXACT, audit/hypotheses/fusion_states.out).

Known obstruction (letter-count model). THM-4916 computes the class law of the values
`K(t.alpha, t.beta) = F(t.beta alpha^dagger)` of the K-form over the words with at most two root changes (`M_2`
weights; one word, hence one free state, per value): the factor at the class 0 is `1 + (3p - 2)(p - 1)/(2 p^3)`
(`34/27` at 3, `151/125` at 5), an excess of order `3/(2p)`, because of the confinement cascade. This is neither
the value `1 + p^-3 - p^-4` of two independent uniform states nor the value 1 of one uniform state under a
non-zero linear functional: the states reached with few root changes are confined. The by-value limit of that
family is OPEN (HYP-9506 (A) and HYP-9660 (A) conjecture different values).

The prime 2. The law fails there in every family: `H(A * B)` is odd (Redei), while the uniform rank-4 law gives
the class 0 modulo 2 the share `17/32`. In the states: `h = H(A - z)` and `H` are path counts, hence odd, so the
first coordinate is odd and the second even; and `U = d (mod 2)`, so the last two coordinates have the same
parity (likewise `V = e (mod 2)` in the right state). In the path-mass coordinates of Session264's proof,
`a_o = (U + d)/2 - h` and `a_e = (U - d)/2` (the masses of the non-trivial root-start paths with an odd, an even
number of vertices; `b_o`, `b_e` likewise from `V`, `e` for the root-end paths), the last two terms of the law are
`2 (a_o b_o + a_e b_e)`. (FINITE-EXACT, audit/hypotheses/fusion_states.out: the law (3) on the 5625 pairs of
labelled pointed tournaments (root = vertex 0) of orders at most 4, the parities for every labelled pointed
tournament of order at most 7; the share `17/32`: audit/hypotheses/rank_laws.out.)

Existing data. The closest data on a family of higher rank are those of the two-apex grammar `G(2,0)`
(HYP-9660 (D), FINITE-EXACT, all `4^n` words, `n <= 16`): the class factors of `H` are within `0.0004` of 1 at
`p = 5, 7, 11, 13` at `n = 16`, which is below the rank-4 bias `4/625` at 5; that `G(2,0)` is a seam of this kind
is not shown there. NUMERICAL (audit, exploratory; audit/hypotheses/fusion_law_small.out): for `A`, `B`
independent and uniform among the labelled pointed tournaments of order `n`, the class-0 factor of `H(A * B)` is
`1.31, 1.30, 1.27, 1.16` modulo 3 for `n = 4, 5, 6, 7` (uniform value `1.0247`), `1.26, 1.18, 1.13, 1.05` modulo 5
(`1.0064`) and `0.25, 0.71, 0.92, 0.993` modulo 7 (`1.0025`). These data do not decide (A).

**(B) Rank 3.** Kernels with three visits (THM-5133 (A)) give forms of rank 3: when no port is visited more than
three times, `H(K[A, B])` is a bilinear form in the states `(P_1, P_2, P_3)` of `A` and of `B` (`P_i` the number of
ordered covers by `i` non-empty paths, `P_1 = H`, `P_2 = Pi_2 - 2 H`; THM-5133), of rank 3 when its `3 x 3` visit
matrix `N` has `det N != 0`. Prediction for uniform states: factor `(1 - p^-3) sigma_2(p^v)/p^(2v)` at the primes
`p` not dividing `det N` (THM-5135 (1), (4)(a)); by THM-5135 (2) this is the law of six squares only at
`p = 1 (mod 4)`.
For module states `(P_1, P_2, P_3)` the prediction can concern only primes `p >= 5` not dividing `det N`:
`P_2 = 2 u_2`, `P_3 = 6 u_3` (`u_i` the number of unordered covers by `i` paths), and modulo 3 every kernel of
THM-5133 (A), whatever its number of visits, is `c(|A| mod 2, |B| mod 2) H(A) H(B)` (THM-4756; audit: all kernels
of order at most 5; see (D)).
Kernels of rank 3 exist (FINITE-EXACT, audit/proof-5135/a6_identity_rank3.out [V]; the session had examined none).
Order 4: of the 768 kernels (labelled tournament, ordered port pair), 48 have a port visited three times; in all
of them both ports are, and `|det N| = 32`, for instance `N = [[5, 3, 0], [3, 7, 2], [0, 2, 2]]` (arcs `0 -> 3`,
`1 -> 0`, `2 -> 0`, `2 -> 1`, `3 -> 1`, `3 -> 2`, ports `(0, 3)`); so for uniform states the rank-3 law holds at
every odd prime. Order 5: of the 20,480 kernels, 3,840 have both ports visited three times and none four times,
all with `det N != 0`, `|det N|` in `{6, 14, 32, 128, 192, 228, 242, 248}`. Neither the session nor the audit
examined the distribution of the module states `(P_1, P_2, P_3)` modulo a prime `p >= 5`.

**(C) Degenerate primes.** PROVED by the audit, now THM-5135 (4)(b): the law depends only on the Smith form of `G`.

**(D) Module states at 3 (constraint, PROVED from THM-4756).** For unrooted module states `(H, Pi_2)` one has
`Pi_2(X) = (-1)^|X| H(X) (mod 3)` (THM-4756, as quoted in THM-5133 (H)). Hence the Gaussian kernel of THM-5133 (B),
`H(K[A, B]) = Pi_2(A) Pi_2(B) + H(A) H(B)`, is `((-1)^(|A| + |B|) + 1) H(A) H(B)` modulo 3: divisible by 3 whenever
`|A| + |B|` is odd. By the same congruences (modulo 3, `P_2 = ((-1)^n + 1) H` and `P_m = 0` for `m >= 3`, `n` the
order; THM-4756, steps 2-3) in the expansion of THM-5133 (A), every module kernel, whatever its number of visits,
is a multiple of `H(A) H(B)` modulo 3: `H(K[A, B]) = c(|A| mod 2, |B| mod 2) H(A) H(B) (mod 3)`, explicitly
`c = N_11 + N_12 s_B + N_21 s_A + N_22 s_A s_B` with `s_X = (-1)^|X| + 1` and `N` the visit matrix (the constant was
written out and checked by the verification of the audit corrections on all kernels of orders 3, 4, 5:
audit/verify-HYP-9666/verify_hyp9666.out [K], 0 violations). So module kernels
do not have the uniform law at 3 and no series `sigma_1(n)/n` can be asserted for them at that prime. The rooted
one-apex states of F1 have no such constraint.
(FINITE-EXACT, audit. audit/hypotheses/kernel_mod3.out: every labelled kernel of order 3, 4, 5 with every ordered
port pair, 48 + 768 + 20,480 kernels, `A` and `B` over the 8 isomorphism classes of orders 1 to 4 with `K[A, B]` of
order at most 10: 0 violations. audit/hypotheses/pi2_mod3.out: `Pi_2 = (-1)^n H (mod 3)` for every labelled
tournament of order at most 7; the Gaussian law and its reduction modulo 3 for one kernel of order 4 on 2129 pairs
`(A, B)` of labelled tournaments with `|A|, |B| <= 4` (a sample of 600 at orders 4 and 4), of which the 1188 with
`|A| + |B|` odd all have `3 | H(K[A, B])`.)
