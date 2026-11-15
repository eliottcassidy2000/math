---
id: THM-5172
title: "Hosts and cores. (A) The invertible root-change semigroup contains no polynomial family at a fixed determinant, hence no orbit of a unipotent (Heisenberg) subgroup: on reachable states H(J^3 s) = 3H + 3L - E >= 3H, so a fixed block skeleton yields O((log X)^(r+1)) counts up to X. (B) Profile congruences: for every tournament of order n, F_2 = (-1)^n H (mod 3), F_2 = 0 (mod 4) for n >= 2, F_2 = 0 (mod 8) for even n >= 4. (C) The cyclic join is a Z[sqrt(-2)] inner product: 3 H(A -> B -> c -> A) = Re(alpha_A conj(alpha_B)), alpha = (f - 3h) + f sqrt(-2); hence N = 2 h_A h_B (mod 3), or N = 0 (mod 3) when both cores have odd order, 3N = u_A u_B (mod 32), and a core joined to itself gives 3N = u^2 + 2 f^2 and N = 3 (mod 8)"
status: >
  PROVED (elementary). (A): invariants by induction on the letters; the J^3 identity is exact; the counting is
  direct. (B): from the odd-cycle form of the path-cover polynomial (THM-4900; its even coefficients w_j).
  (C): a rewriting of the reflection identity 3N - 2 f_A f_B = (f_A - 3h_A)(f_B - 3h_B) of
  chanpak_absolute_reflection_20261114 (KNOWN, credited; that note's identity and the cyclic join are not new
  here), followed by the congruences of (B).
  FINITE-EXACT checks:
  - (A): the invariants and the J^3 bound on all 7,478,943,733 reachable states of the full grammar with
    H <= 2^18 (scripts/invariants.c); 45, 211, 528, 946, 1486 words J^a0 JOZ J^a1 with count <= 10^3, 10^6, 10^9,
    10^12, 10^15 (scripts/polylog_skeleton.py).
  - (B): all 9,733,056 classes of order 10 and all classes of orders <= 9 (gentourng census); zero violations.
  - (C): 100,000 random core pairs (scripts/join_check.py also checks the cyclic-join count itself by brute force
    on 300 random joins).
  NOT CLAIMED: any coverage statement; Chan-Pak (OPEN); the band-fullness reduction, which is HYP-9692.
session: mac-mini-2026-11-14-hosts
source: 05-knowledge/results/chanpak_waves_20261113/README.md#wave-25-hosts-and-cores--no-heisenberg-orbit-the-cyclic-join-in-zsqrt-2-and-band-fullness
scripts:
  - 05-knowledge/results/chanpak_waves_20261113/scripts/invariants.c
  - 05-knowledge/results/chanpak_waves_20261113/scripts/polylog_skeleton.py
  - 05-knowledge/results/chanpak_waves_20261113/scripts/join_check.py
  - 05-knowledge/results/chanpak_waves_20261113/scripts/prof.c, profhf.c (profile census from gentourng)
related:
  - THM-5162 (the semigroup), HYP-9682, HYP-9692 (band fullness), THM-4900 (odd-cycle form of P_T(q)),
    chanpak_absolute_reflection_20261114 (cyclic join, reflection identity, dyadic windows)
---

# Hosts and cores

Notation:
- P_T(q) is the number of ordered q-slot path covers of T, empty slots allowed.
- h = P_T(1) = H(T) and f = P_T(2) = F_2(T).
- u = f - 3h.
- The cyclic join of nonempty tournaments A, B is the tournament A -> B -> c -> A: all arcs from A to B, from B to
  a new vertex c, and from c to A.

## (A) No Heisenberg orbit in the invertible semigroup

Letters act on (H, E, L, h) from (1, 1, 1, 1) as in THM-5162.

**Invariants.** On every reachable state,

    0 <= E <= H,     1 <= h <= H,     L >= h,     L >= E.

*Proof (induction on letters).*
- **O** gives (H, 0, 2L, h). It keeps every inequality: 2L >= L >= h.
- **Z** gives (H + L - E, h, 2L - h, h). It keeps them: L - E >= 0, so H does not decrease; E' = h <= H <= H';
  and L' = 2L - h >= h = E' because L >= h.
- **J** gives (H + L - E, H - E, H + L, H). It keeps them: H' >= H >= E' >= 0, h' = H <= H', and
  L' = H + L >= H = h'. Also L' >= H >= H - E = E'. ∎

So H never decreases along a word. Moreover, J^3 sends H to

    H(J^3 s) = 3H + 3L - E >= 3H + 2L >= 3H.

Hence a run of a root changes multiplies the count by at least 3^floor(a/3).

**Proposition.** Fix a block skeleton: the non-J blocks B_1, ..., B_r and a source. The words
J^(a_0) B_1 J^(a_1) ... B_r J^(a_r) with count <= X number at most 3^(r+1) C(log_3 X + r + 1, r + 1), which is
O((log X)^(r+1)). The invertible elements of a fixed determinant 2^K have finitely many skeletons, so they
produce O_K((log X)^(K+1)) counts up to X.

**Corollary (no unipotent family).** Let P, Q be in GL_3(Q), let s be a fixed source state, and let
u_t = exp(tX) with X ≠ 0 nilpotent. Then P u_t Q lies in the semigroup for no infinite arithmetic progression of t.

*Proof.* All such products share one determinant. They are pairwise distinct, so they come from distinct words.
Their counts e_1 P u_t Q s are a polynomial in t of degree at most 2.
- If that polynomial is nonconstant, a progression gives ≫ T distinct words with counts at most O(T^2),
  contradicting the polylogarithmic bound.
- If it is constant, infinitely many words share one count, again a contradiction.

The same holds for any orbit of a lattice in the Heisenberg group of THM-5162 (B). ∎

So the Bourgain–Kontorovich mechanism for Apollonian packings, in which an arithmetic stabilizer contributes
polynomial families, has no counterpart inside this semigroup. Polynomial freedom has to come from elsewhere. The
host/core constructions supply it: two full profile coordinates per core.

## (B) Profile congruences

By THM-4900,

    P_T(q) = q^n + sum_{j >= 1} w_j q^(n-2j),

where each w_j is a nonnegative even integer (weights 2^|K| of odd-cycle packings). So P_T(-q) = (-1)^n P_T(q).
Hence:
- f = P(2) ≡ P(-1) = (-1)^n P(1) = (-1)^n h (mod 3);
- for n >= 2, f ≡ 0 (mod 4), because 2^n ≡ 0 and each w_j 2^(n-2j) has n - 2j >= 1 and w_j even;
- for even n >= 4, f ≡ 0 (mod 8), because then n - 2j >= 2 for every j. ∎

## (C) The cyclic join as an inner product in Z[sqrt(-2)]

Put alpha_T = u_T + f_T sqrt(-2). The reflection identity of chanpak_absolute_reflection_20261114,
3N - 2 f_A f_B = u_A u_B for N = H(A -> B -> c -> A), reads

    3N = u_A u_B + 2 f_A f_B = Re(alpha_A conj(alpha_B)).

- **Geometry.** The first core's ratio h/f fixes the angle of alpha_A. Once N and A are fixed, the second core's
  vector lies on the lattice line Re(alpha_A conj(beta)) = 3N, perpendicular to alpha_A. This is the owner's
  "bounded arithmetic progression": in coordinates, a f_B - c h_B = N with a = f_A - h_A and c = f_A - 3h_A, cut
  off by h_B >= 1 and f_B >= (n_B + 1) h_B.
- **Mod 3.** With f ≡ eps h (eps = (-1)^n), N ≡ h_A h_B ((eps_A - 1)(eps_B - 1) + 2) (mod 3). That is, N ≡ 2 h_A h_B
  if a core has even order, and N ≡ 0 if both cores have odd order. A prime target other than 3 therefore needs a
  core of even order.
- **Mod 32.** Both f are divisible by 4, so 2 f_A f_B ≡ 0 (mod 32) and 3N ≡ u_A u_B (mod 32).
- **Diagonal joins.** A -> A -> c -> A gives 3N = u^2 + 2f^2, a norm from Z[sqrt(-2)]. Here u is odd, so u^2 ≡ 1
  (mod 8) and N ≡ 3 (mod 8). Examples: 3 (singleton), 11 (TT_2), 99 (C_3), 779 (the (5, 32) core).
