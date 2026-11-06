---
id: THM-4888
title: "Finite O/Z tail banks: (A) Q(m) = m^2 - r(m) (mod 3), r the number of runs of 1-bits, so H_s(m) = h + Rm - d m^2 + d r(m) - e(m mod 2) (mod 3) for every seed; (B) H_s(T(m)) - H_s(m) = (R - 2d)(T(m) - m) - e(par T(m) - par m) for Session 206's converse involution T, and the tail image of a physical seed is T-symmetric iff R = 2d and e = 0 iff the seed is the rooted edge; (C) KNOWN (Session 206 weighted-parameter-capacity): confinement to the classes h, h - e mod p | gcd(d, R) and the parity-split capacity, with the corollary that a J-born state over x with 3 | gcd(D, U) avoids 1 (mod 3) unless E = 2 (mod 3); (D) redundancy: JZ = OJ (KNOWN, Session 205 audit) nests W O J in W J; seeds agreeing in (D - E, U, d), resp. (D, U, d), share their odd-m, resp. even-m, images; R = 2d seeds are dilated root images (h - d) + d H_t(m) - e(m mod 2)"
status: >
  PROVED (new): (A), (B), (D)2, (D)3 (proofs in the packet PROOF.md, section 2). KNOWN / CITED, found independently
  in this session: (C) is contained in Session 206's audited weighted-parameter-capacity theorem (de303ea26d, pushed
  during this session; it gives the exact weights (g/p) c_b, g = gcd(R, d, p)); (D)1 JZ = OJ is Session 205's
  finite-oz-tail-bank audit identity (Session 199: ZJ = JO). Also KNOWN: the tail formula (Session 205), the capacity
  theorem, the root involution T with Q(T(m)) - Q(m) = 2(T(m) - m), the d = 1 collapse, root residue capacity and the
  example (5,6,2,3) (Session 206), and the complete source grammar and exact inverse (Session 198 bilinear-descent).
  Pre-registered as THM-4887 (d36f57c903), claims (A) and (B) only; renumbered after a reservation race
  (48ba576e1f). NOT CLAIMED: impossibility for every finite bank, an image-density limit, coverage; the quantitative
  laws are HYP-9480. Chan-Pak OPEN. AUDIT: PENDING.
session: kind-pasteur-2026-11-06-S13
source: 05-knowledge/results/finite_oz_bank_poisson_root_changes_20261106/PROOF.md
related:
  - HYP-9480
  - THM-4878
  - THM-4885
  - THM-4883
  - session205/finite-oz-tail-bank/ARITHMETIC_REDUCTION.md
  - session205/finite-oz-tail-bank/audit/REVIEW.md
  - session206/parameter-capacity/PARAMETER_CAPACITY.md
  - session206/weighted-parameter-capacity/WEIGHTED_CAPACITY.md
  - session206/root-oz-noncofiniteness/ROOT_OBSTRUCTION.md
  - session206/root-residue-capacity/ROOT_RESIDUE_CAPACITY.md
  - session206/exact-modular-image/DYADIC_SURJECTIVITY.md
  - session198/bilinear-descent/PROOF.md
  - session199/hidden-structure/dependencies/BALANCED_SOURCES.md
---

# THM-4888 - finite O/Z tail banks: mod-3 law, involution rigidity, redundancy

Setting: a physical unit-start seed s = (h, R, e, d) = (D, U, E, d) followed by an O/Z word with parameter m has
Hamiltonian-path count H_s(m) = h + Rm - dQ(m) - e(m mod 2) (Session 205), Q(0) = 0, Q(2m) = Q(m),
Q(2m+1) = Q(m) + m + (m mod 2). Moves: O(D,E,U,d) = (D,0,2U,d), Z = (D-E+U,d,2U-d,d), J = (D-E+U,D-E,D+U,D). The
rooted edge is t = (h, R, e, d) = (1, 2, 0, 1).

**(A) Mod-3 law (PROVED).** Q(m) = m^2 - r(m) (mod 3), where r(m) is the number of maximal runs of ones in the
binary expansion of m. Proof: with 2 = -1 (mod 3) the digit identity Q = sum b_j floor(m/2^(j+1)) + sum b_j b_(j+1)
becomes -(S^2 - s_2)/2 + a with S = sum (-1)^j b_j = m, s_2 the digit sum and a the number of adjacent 11 pairs;
-1/2 = 1 (mod 3) and s_2 - a = r.

**(B) Involution rigidity (PROVED).** H_s(T(m)) - H_s(m) = (R - 2d)(T(m) - m) - e((T(m) mod 2) - (m mod 2)). The
pairs 9 <-> 11 (same parity) and 4 <-> 7 (different parity) force R = 2d and e = 0 for a T-symmetric image.
Physically U = 2d and E = 0 give the positive coordinates (0, 0, 0, d). Under the exact inverse these have no
predecessor, and the only physical source without a predecessor is t itself (d = 1). So Session 206's halving, which
comes from the converse tournament, is unique to the root image (and its d = 1 subsets).

**(C) Residue confinement and parity-split capacity (KNOWN: Session 206 weighted-parameter-capacity).** If an odd
prime p divides d and R, then H_s(m) = h - e(m mod 2) (mod p). Even and odd parameters carry (X/d) ln(R/(R - d/2)) and
(X/d) ln((R - d/2)/(R - d)) (+ O_s(sqrt X)) respectively. Session 206 proves more: the exact weights (g/p) c_b on the
coset h - e b mod g, g = gcd(R, d, p). This session found the special case independently. Corollary (immediate): a
J-born state over x with 3 | gcd(D, U) has its whole tail image in the classes 0 and -E (mod 3), which is the
mechanism behind the hard class 1 (mod 3) of root-change and join families.

**(D) Redundancy.**
1. (KNOWN: Session 205 audit.) JZ = OJ, both giving (D + 2U, D, D + 2U, D), so the image of every seed W O J lies in
   the image of W J.
2. (PROVED.) The odd-m values (h - e) + Rm - dQ(m) depend only on (D - E, U, d), which also fixes the Z-child. The
   even-m values depend only on (D, U, d), which also fixes the O-child.
3. (PROVED.) If R = 2d, then H_s(m) = (h - d) + d H_t(m) - e(m mod 2): two dilated pieces of the root image.

These identities explain why finite banks are positively correlated (HYP-9480 (B)). In bank A, the 35 balanced seeds
are redundant, and the union misses six times the product of the members' miss rates.

Verification: packet `verify.py` (P5) checks (A) for all m < 2^26 and (B) on all 194 physical seeds with U <= 64 and
m < 2^16. Among all 94,010 physical states with U <= 4096, only the rooted edge has U = 2d and E = 0.
