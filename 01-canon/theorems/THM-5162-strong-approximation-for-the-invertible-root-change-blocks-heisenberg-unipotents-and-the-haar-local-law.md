---
id: THM-5162
title: "Strong approximation for the invertible root-change blocks. (A) Shape: every invertible block J O v Z acts on the rooted triple (H, E, lam) by [[b - Q, -1, b], [1, 0, 0], [2^n - R, 0, 2^n]], where b - 1 is the binary number of the Z-positions, R its n-bit reversal and Q their pair energy; J w is invertible iff w is empty or w = O v Z. (B) Two blocks of equal length differ by a unipotent element of one Heisenberg group. (C) The invertible blocks generate {g in GL_3(F_p) : det g in <2>} modulo every odd prime p, and {g in GL_3(Z/q) : det g in <2>} modulo every odd squarefree q. (D) Haar local law: from any p-primitive state, long random invertible words give counts distributed mod p as the first coordinate of a uniform nonzero vector of F_p^3"
status: >
  PROVED (elementary; exact computations recorded below). (A), (B), (D): proofs complete in this file. (C): proof
  complete for every odd prime p >= 5 and every odd squarefree q whose prime factors are >= 5; p = 3 rests on the
  exhaustive search of wave 23 (the seven blocks generate all 11,232 elements of GL_3(F_3)), so (C) holds for every
  odd squarefree q.
  FINITE-EXACT checks: the closed forms of (A) for all 2047 blocks with |w| <= 12; the invertibility criterion for
  all 2047 blocks J w with |w| <= 10; the unipotent shape of (B) for all 155 same-length pairs with |w| <= 6; the
  witness determinant 3/4 and the commutator in (C).
  NUMERICAL (by-value census, not part of the proof): the law (D) and the source law (E) on [2^12, 2^24); see
  HYP-9682.
  NOT CLAIMED: uniformity of the mixing rate in the modulus (expansion; would follow from Salehi Golsefidy-Varju
  applied to the commutator subgroup, not written out here); by-value (count-ordered) versions of (D); density one;
  Chan-Pak (OPEN).
  KNOWN and used: the O/Z/J grammar and its block matrices (THM-4904, THM-5148); JZ = OJ (THM-4888 (D) 1); the rank-one
  resets need a singular block (THM-5120); Zariski density of the seven blocks in GL_3 and the exhaustive images mod 3
  and mod 5 (wave 23 of chanpak_waves_20261113). The Heisenberg-unipotent argument, the shape formula with the
  reversal R and the surjection at every odd prime are NEW here (no strong-approximation statement in canon before
  this file).
session: mac-mini-2026-11-14-density
source: 05-knowledge/results/chanpak_waves_20261113/README.md#wave-24-the-invertible-root-change-orbit--reset-free-haar-distributed-and-complete-beyond-6223
scripts:
  - 05-knowledge/results/chanpak_waves_20261113/scripts/shape_lemma.py
  - 05-knowledge/results/chanpak_waves_20261113/scripts/unipotent.py
  - 05-knowledge/results/chanpak_waves_20261113/scripts/modp_image.py (wave 23; p = 3, 5)
related:
  - HYP-9682 (the invertible root-change conjecture this theorem supplies the local half of)
  - THM-5120, THM-5148, THM-4904, THM-5128, THM-5134, HYP-9660, HYP-9677, HYP-9494
---

# Strong approximation for the invertible root-change blocks

Letters act on the root state (H, E, L, h), starting from (1, 1, 1, 1):

    O: (H, 0, 2L, h)      Z: (H + L - E, h, 2L - h, h)      J: (H + L - E, H - E, H + L, H).

Every word builds a rooted tournament whose Hamiltonian-path count is H (THM-4904, THM-5120). A *block* is J
followed by a maximal O/Z run w. Since J ignores h, a block acts linearly on (H, E, lam), lam = L (THM-5148).

## (A) Shape of a block

**Lemma A.** Let w = O v Z with v in {O, Z}^(n-2), and let S be the set of Z-positions of w (1-based; 1 is not in S,
n is). Put

    P = sum_{s in S} 2^(s-1),   b = 1 + P,   R = sum_{s in S} 2^(n-s),
    Q = sum_{t in S} ( sum_{s in S, s < t} 2^(t-1-s) + [t-1 in S] ).

Then J w acts on (H, E, lam) by

    M_w = [[b - Q, -1, b], [1, 0, 0], [2^n - R, 0, 2^n]],   det M_w = 2^n.

As v ranges over {O, Z}^(n-2), b ranges over the odd numbers of (2^(n-1), 2^n), each once, and R = rev_n(P) is the
reversal of the n-bit pattern of P. The block J alone is [[1, -1, 1], [1, -1, 0], [1, 0, 1]], det 1. A block J w is
invertible iff w is empty or w = O v Z.

*Proof.* J sends (H, E, lam) to (H1, E1, L1, h1) = (H - E + lam, H - E, H + lam, H). Along the run, O doubles L and
sets E = 0, and Z adds L - E to H, then sets E = h1 and L to 2L - h1. Hence before position t,

    L_(t-1) = 2^(t-1) L1 - h1 sum_{s in S, s < t} 2^(t-1-s)   and   E_(t-1) = h1 [t-1 in S]   (t >= 2).

Summing the increments over t in S gives H' = H1 + P L1 - Q h1 = (1 + P - Q) H - E + (1 + P) lam. The last letter is
Z, so E' = h1 = H. Finally lam' = L_n = 2^n L1 - R h1 = (2^n - R) H + 2^n lam. Expanding the determinant along the
second row gives 2^n.

The b-values form the stated range because P is any even number whose bit n - 1 is set.

For invertibility:
- If w ends in O, the second row of M_w vanishes.
- If w starts with Z, then J Z sends (H, E, lam) to (H + 2 lam, H, H + 2 lam), which does not depend on E, so M_w
  is singular.

Every nonempty w either starts with Z, or starts and ends with O, or has the form O v Z. ∎

Checked exactly: the closed forms for all 2047 blocks with |w| <= 12, and the invertibility criterion for all
2047 blocks with |w| <= 10 (scripts/shape_lemma.py).

So the last block reads a count as

    N = b (H + lam) - E - Q H,     and the next E is the current H.

The pair (H, E) is a shift register of consecutive counts, as in a negative continued fraction
q_(k+1) = a q_k - q_(k-1), with lam as a carry. The perturbation Q is the binary pair energy of the one-apex count
(THM-5093). The pair (P, R) is a point of THM-5141's reversal graph.

## (B) Heisenberg unipotents

**Lemma B.** Let A = M_w and B = M_(w') with |w| = |w'| = n. Write a, b, c for the entries (1,1), (1,3), (3,1). Then

    A^(-1) B = [[1, 0, 0], [alpha, 1, beta], [gamma, 0, 1]],
    beta = b_A - b_B,   gamma = (c_B - c_A)/2^n,   alpha = a_A - a_B + b_A gamma.

In the ordered basis (e_H, e_lam, e_E), this matrix is lower unitriangular. So every such quotient lies in the
unipotent radical U of the Borel subgroup stabilizing the flag F: <e_E> ⊂ <e_lam, e_E>. Its image in
U/[U, U] ≅ G_a^2 is (gamma, beta). Commutators lie in the centre:

    [u1, u2] = I + (beta1 gamma2 - beta2 gamma1) E_(E,H).

*Proof.* Use A^(-1) = [[0, 1, 0], [-1, a - bc/d, b/d], [0, -c/d, 1/d]] with d = 2^n, and multiply. ∎

Checked exactly for all 155 same-length pairs with |w| <= 6 (scripts/unipotent.py).

**Corollary B' (Hammersley digits).** Fix the anchor A_n = M_(O^(n-1) Z). Its entries are
P = 2^(n-1), Q = 0, R = 1, so A_n = [[2^(n-1)+1, -1, 2^(n-1)+1], [1, 0, 0], [2^n - 1, 0, 2^n]].

Every block of length n is then A_n u, with u in U. Write the block as w = O v Z, and let k in [0, 2^(n-2)) be v
read as a binary number, Z = 1, with the first letter of v as the lowest bit. Then the image of u in U/[U, U] is

    (gamma, beta) = (-2 rev_(n-2)(k)/2^n, -2k).

So the 2^(n-2) blocks of length n are one anchor times a Heisenberg "digit" whose abelianized coordinates form the
Hammersley point set {(k, rev_(n-2)(k))}. This is a (0, n-2, 2)-net in base 2: every dyadic box of the right area
holds exactly one point.

*Proof.* By Lemma A, P = 2^(n-1) + 2k and R = 1 + 2 rev_(n-2)(k). Then apply Lemma B. ∎ Checked for n = 5, 6, 7.

The invertible orbit is thus a walk whose step at scale 2^n is a fixed anchor times an equidistributed Heisenberg
digit. The reversal graph of THM-5141 is exactly this digit set.

## (C) Strong approximation

**Theorem.** Let Gamma be generated by J, JOZ, JOOZ, JOZZ, JOOOZ, JOZOZ, JOZZZ, or by any set of invertible blocks
that contains them. Then:
- for every odd prime p, Gamma mod p = G_p := {g in GL_3(F_p) : det g in <2>};
- for every odd squarefree q, Gamma mod q = {g in GL_3(Z/q) : det g in <2 mod q>}.

Since these images are finite, the semigroup generated by the blocks has the same images.

*Proof.*

*Step 1: two unipotents.* Take

    u1 = JOOZ^(-1) JOZZ = [[1, 0, 0], [-5/4, 1, -2], [-1/4, 0, 1]],
    u2 = JOOOZ^(-1) JOZZZ = [[1, 0, 0], [-27/8, 1, -6], [-3/8, 0, 1]].

Their images in U/[U, U] are (-1/4, -2) and (-3/8, -6), with determinant 3/4. Their commutator is
I - (3/4) E_(E,H), a transvection.

*Step 2: the Heisenberg group mod p (p >= 5).* Here 3/4 is a unit in F_p, so the images of u1, u2 span
U(F_p)/[U, U] = F_p^2. In the p-group U(F_p), p odd, [U, U] is the Frattini subgroup. By the Burnside basis theorem,
u1 and u2 generate U(F_p).

*Step 3: an opposite flag.* J has J_(H,E) = -1 and, since det J = 1, (J^(-1))_(H,E) = 1. Hence
J<e_E> ⊄ <e_lam, e_E> and <e_E> ⊄ J<e_lam, e_E>. So JF is opposite to F in every characteristic. By the Bruhat
decomposition, JF = b w0 F for some b in B_F(F_p). Then J U J^(-1) = b U^- b^(-1), where U^- is the opposite
unipotent group, and

    <U, J U J^(-1)> = b <U, U^-> b^(-1) = SL_3(F_p),

because elementary matrices generate SL_3 over a field. So Gamma mod p contains SL_3(F_p).

*Step 4: determinants.* det JOZ = 4 and det JOOZ = 8 generate <2>, and every block has determinant a power of 2.
Hence Gamma mod p = G_p.

*Step 5: p = 3.* The exhaustive search of wave 23 (scripts/modp_image.py) finds 11,232 = |GL_3(F_3)| elements, and
<2> = F_3^*. The Heisenberg route genuinely fails at 3: over all 155 same-length pairs with |w| <= 6, the odd part
of the gcd of the 2 x 2 determinants is exactly 3.

*Step 6: squarefree q.* SL_3(F_p) is perfect for p >= 3, and G_p/SL_3(F_p) is abelian. So the commutator subgroup
of Gamma mod q projects onto every factor SL_3(F_(p_i)). For distinct p these factors have no isomorphic nontrivial
quotient: each such quotient has PSL_3(p) as a composition factor, and |PSL_3(p)| determines p. By Goursat's lemma,
applied inductively, the commutator subgroup is all of prod_i SL_3(F_(p_i)) = SL_3(Z/q). Determinants give the
rest. ∎

## (D) The Haar local law

**Corollary D.**
- Let mu be a probability measure on invertible blocks whose support contains the seven blocks of (C).
- Let v be a state in Z^3 with v ≢ 0 mod p, for an odd prime p.

Then the law of M_(b_n) ... M_(b_1) v mod p, with the b_i i.i.d. mu, converges exponentially fast to the uniform
law on F_p^3 \ {0}. In particular the count H converges in law to

    nu_p(0) = (p^2 - 1)/(p^3 - 1),     nu_p(c) = p^2/(p^3 - 1)   (c ≠ 0).

The class 0 is slightly depleted.

*Proof.* The walk on the finite group G_p is irreducible by (C). Suppose supp mu lay in a coset gN of a proper normal
subgroup N. Then N would contain u1, u2 and the quotients J^(-1) JOZ (det 4) and J^(-1) JOOZ (det 8). The
G_p-conjugates of U(F_p) generate SL_3(F_p), so N would contain SL_3(F_p) and the determinants <4, 8> = <2>, so
N = G_p. Hence the walk is aperiodic and converges to Haar measure on G_p. Finally SL_3(F_p) is transitive on
nonzero vectors.

For p = 3, G_3 = GL_3(F_3) = SL_3(F_3) x {±I}, whose normal subgroups are 1, {±I}, SL_3(F_3) and the whole group.
J and JOZ are not ± each other mod 3 (their (1,1) entries are 1 and 0), so the blocks do not lie in a coset of
{±I}. Their determinants
1 (J), 4 ≡ 1 (JOZ) and 8 ≡ -1 (JOOZ) differ, so they do not lie in a coset of SL_3(F_3). The same conclusion
follows. ∎

**(E) Source law (immediate).** Invertible blocks are invertible mod p. So in the invertible root-change language
(a one-apex source followed by invertible blocks) a state is p-divisible after the first J iff the source state
was. Hence the limiting class law at p is

    pi_p delta_0 + (1 - pi_p) nu_p,

where pi_p is the share of words whose source is p-divisible. In particular no reset (THM-5120) and no
confinement can occur after the source.

The *prime singular series* of an orbit without p-divisible states is

    prod_{p odd} (1 - nu_p(0))/(1 - 1/p) = prod_{p odd} p^3/(p^3 - 1) = (7/8) zeta(3) = 1.05180.

The census (HYP-9682) measures 1.0526 on [2^23, 2^24).
