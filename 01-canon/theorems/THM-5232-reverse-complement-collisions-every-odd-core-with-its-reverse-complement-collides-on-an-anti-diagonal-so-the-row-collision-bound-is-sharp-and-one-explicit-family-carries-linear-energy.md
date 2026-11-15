---
id: THM-5232
title: "Reverse-complement collisions: in a carrier row G(u) = cu + kp(u) + phi[u odd] with k | c 2^(j-1), every odd a-bit core x with top bit 1 collides with its reverse-complement rc_a(x) = 2^a - 1 - rev_a(x) along the anti-diagonal w + r = (1 + c/k) 2^j - 1, because p(rc x) = p(x) and x - rc x = rev x - rev rc x; hence N(l, w) >= (c/2k) 2^l / w on a positive-log-proportion set of prefixes in every dyadic block (the shape of HYP-9752 is sharp), and this one explicit family has 2^(-1 - ceil(log2(1 + c/k))) N members below N (44% of all collisions of the 3-cycle row at L = 18)"
status: >
  PROVED (Lemma 1 by an inversion-pair symmetry; the theorem from the two-ended lattice-line equation, which is
  re-derived here and is equation (1) of chanpak_cone_collision_placement_20261115/TWISTED_LATTICES.md, credited).
  FINITE-EXACT replay:
  - Lemma 1 for all a <= 16;
  - for four rows at L = 18: every family member is a literal collision, the member count equals the closed form,
    and the energy shares are as stated;
  - the sharpness constant c/(2k).
  NUMERICAL: among palindromic-difference cores, equal p occurs exactly for the reverse-complement cores (plus two
  sporadic cores at odd a), for a <= 20.
  NOT CLAIMED: the upper bound of HYP-9752 (OPEN); positive density; Chan-Pak (OPEN).
session: klein-2026-11-15-rc
date: 2026-11-15
source: 05-knowledge/results/chanpak_reverse_complement_collisions_20261115/RC_COLLISIONS.md
scripts:
  - 05-knowledge/results/chanpak_reverse_complement_collisions_20261115/verify_rc_collisions.py
  - 05-knowledge/results/chanpak_reverse_complement_collisions_20261115/coll_classes.py, hd_family.py, pal_struct.py
related:
  - HYP-9752 (row collision lemma; this proves its lower-bound shape), THM-5192 (A) (carrier form),
    chanpak_dyadic_collision_ports_20261115 (split identity, signature rigidity),
    chanpak_cone_collision_placement_20261115/TWISTED_LATTICES.md (two-ended equation, single-core rays)
---

# THM-5232 - reverse-complement collisions

Let p(0) = 0, p(2v + 1) = p(v), p(2v) = p(v) + v + (v mod 2), and G(u) = cu + kp(u) + φ(u mod 2) with c, k ≥ 1. Write
rev_a for bit reversal in an a-bit window and rc_a(x) = 2^a − 1 − rev_a(x).

**Lemma 1.** For every a ≥ 1 and 0 ≤ x < 2^a:
- p(rc_a x) = p(x);
- x − rc_a x = rev_a x − rev_a(rc_a x).

*Proof.* The recursion unrolls to

    p(x) = Σ_(i<j, b_i=0, b_j=1) 2^(j−i−1) + #{i : b_i = 0, b_(i+1) = 1}.

Reverse-complement maps b to b′_i = 1 − b_(a−1−i). That maps each inversion pair (i, j) to the pair
(a−1−j, a−1−i), which has the same gap, and each pattern 01 to a pattern 01, so p is invariant. For the second
identity, rev_a(rc_a x) = 2^a − 1 − x, so both sides equal x + rev_a x − 2^a + 1. ∎

**Theorem.** Let m = 1 + c/k. Take a ≥ 2, j ≥ 1 with k | c 2^(j−1), an odd x ∈ [2^(a−1), 2^a) and r ∈ [0, 2^j). Put
y = rc_a x, t = rev_j r and w = m 2^j − 1 − r. Then

    G(2^j(2^a w + x) + t) = G(2^j(2^a w + y) + t),

and distinct (a, j, x, r) give distinct ordered pairs, with top split level l = a + j and prefix w.

*Proof.* The two-ended split at the common suffix t and the common prefix w gives the collision equation

    d w + δ r + [w odd] − ε[r even] = ((c/k) 2^j + 2^j − 1) δ + p(x) − p(y),

with δ = x − y, d = rev_a x − rev_a y and ε = (x mod 2) − (y mod 2) = 1. By Lemma 1, d = δ and p(x) = p(y), so the
equation reads δ(w + r − m 2^j + 1) + [w odd] − [r even] = 0. Our w makes the first term zero. Since m 2^j is even,
w + r is odd, so the parities of w and r differ and the second term is zero. ∎

## Consequences

1. **Sharpness of HYP-9752.** Suppose w ∈ [(m−1)2^j, m 2^j − 1] and a = l − j ≥ 2. Then every one of the 2^(a−2)
   cores of span a contributes a collision at (l, w). So

       N(l, w) ≥ 2^(l−j−1) ≥ (c/(2k)) 2^l / w        (ordered pairs).

   These prefixes meet every dyadic block in an interval of logarithmic length log₂(m/(m−1)). Hence the bound
   N(l, w) ≤ C 2^l/(1 + w) is sharp in shape, and C ≥ c/(2k): this gives 3 for the 3-cycle row and 1 for 2U + p.
2. **An explicit family of linear size.** If k | c, the members below N = 2^L number

       Σ_(j≥1) Σ_(a≥2) 2^(a−2) clip(2^(L−j−a) − (m−1) 2^j, 0, 2^j),

   which is about 2^(−1−⌈log₂ m⌉) N. At L = 18 this is 44.1% of all collisions of 18U + 3p + 2[odd], 19.7% for
   2U + p and 17.9% for 12U + 3p. It also explains why HYP-9752's per-level totals grow linearly in L − l.
3. **Arithmetic dichotomy.** The family needs odd(k) | c. When that fails, anti-diagonal collisions survive only
   at the single core difference δ(x) = k/gcd(k, c 2^j). The row 189U + 11p has no family member through L = 18,
   and its energy is 50–250 times smaller than that of the physical rows. Physical carrier rows have k | c.

## Interpretation (heuristic)

In a carrier tournament (THM-5192 (A)), the converse turns the X-tape encoding (most-significant-first, mark =
digit) into the Y-tape encoding (least-significant-first, mark = complemented digit), so rc is the digital trace
of H(T) = H(T^op). The theorem says that converting one middle tape block in this way, at a matched prefix, keeps
the count.

Proofs of HYP-9752 must therefore subtract the explicit lattice-ray families, and bound only the generic
single-point lines.
