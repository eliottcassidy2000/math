# Quotient Hamilton-cycle parity and the odd ambient correction

This is a separate structural corollary. It compresses an ambient correction test; it does not assert that a local separator exists for an arbitrary quotient. The main prime-five classifications do not depend on publishing this corollary first.

## 1. Every child-uniform physical cycle parity collapses

Let Q be any tournament on m nonempty child blocks A_i. Attach actual ports x->y with incidences constant on each whole child. Let G be the resulting physical extension and R the corresponding extension of the bare quotient Q, with the same two block incidence masks. Let Z(T) count directed Hamilton cycles up to cyclic rotation, in BOTH odd and even orders.

Then

    Z(G) = Z(R) (mod2).                                      (1)

Root each cycle at the unique physical port x. Its maximal child-run word is now unambiguously linearized and retains the closing arc back to x. If child i is used in k_i runs, its fillings are exactly one ordered k_i-path cover, contributing q_(k_i)(A_i). The two ports are single occurrences. Since q_k(A_i) is divisible by k! and H(A_i)=q_1(A_i) is odd, modulo two every term with a repeated child vanishes. The remaining words visit every child once and are precisely the Hamilton cycles of R rooted at x. Their filling factor is product_i H(A_i)=1 modulo two. This proves (1), including singleton children and arbitrary child interiors.

The unique external port is essential to this particular proof: it gives a distinguished cyclic starting position and prevents a cycle contained wholly inside one child. No cyclic-word orbit has been divided by its length.

## 2. A direct formula for the normalized two-cover functional

Write U=product_i H(A_i), an odd integer. Let Phi=H(G)/U and Psi=q_2(G)/U denote normalized residues modulo eight. Let epsilon be the parity of the total child order, equivalently the parity of |G|, and let Z=Z(R) modulo two.

The previously proved profile identity is

    q_2(G) = 2H(G)+2^|G|-4-4 epsilon Z(G) (mod8).

Here |G|>=3, so the power of two vanishes modulo eight. Using (1) and 4U^(-1)=4 modulo eight gives

    Psi = 2 Phi - 4 - 4 epsilon Z (mod8).                    (2)

The cycle count Z is computed only on the m+2-vertex quotient extension, regardless of child sizes. This formula retains the actual odd scalar U; it does not equate raw responses for differently scaled child profiles.

For the five profile types, epsilon=sum_i(s_i+c_i) modulo two. In the notation of coarse_mod8_affine, Phi modulo four is

    Phi = K0+2 sum_i(1-s_i)K_i (mod4).

Thus an equivalent version of (2) is

    Psi = 2K0+4 sum_i(1-s_i)K_i-4
          -4 Z sum_i(s_i+c_i) (mod8).                       (3)

The entire two-cover type dependence is therefore determined by the Hamilton base K0, the parities of the double-run coefficients K_i, and ONE quotient cycle parity Z. In particular all odd-cycle-bit coefficients coincide.

As a check on the independently established two-cover extraction coefficients B0,B_i,C_i, equations(2)-(3) imply

    B0 = 2K0-4-4(m mod2)Z (mod8),
    B_i = 2K_i+2Z (mod4),
    C_i = 2Z (mod4),                    for every i.        (4)

These identities concern normalized coefficient geometry, not arbitrary choices of artificial profiles.

## 3. Exact scope of the ambient criterion

Compare Q[A_i] with the converse quotient Q^op[B_i], on the same child name sets, where every corresponding complete child profile agrees. For a fixed child-uniform query let R and R' be the two bare quotient extensions. Their Hamilton-cycle parity difference equals that of the actual augmented cores by (1).

If the total child order is even, the known even-order modulo-eight transfer applies without a cycle test. If that order is odd and Z(R)=Z(R') modulo two, the known odd-order correction vanishes. In either case EVERY nonzero modulo-eight local Hamilton contrast survives every fixed ambient substitution as an odd multiple modulo eight. Both ports must be lifted with the original core module's exterior incidences.

In particular a local gap of four, together with either even total order or matching bare quotient cycle parities, has complete augmented profile gap (4,0,0,...) modulo eight and remains four in every fixed ambient. This is a sufficient test and does not manufacture the local nonzero gap.

If the bare quotient cycle parities differ and the total child order is odd, no ambient non-annihilation follows from this corollary. The known C3[TT2,TT2,TT3] query a=A1,c=A2 has a local gap of -68 but a cyclic-ambient gap of -600, losing its residue four. The verifier retains this negative control and its unequal bare quotient cycle parities. It is a modular loss, not an exact collision.

## 4. Controls and dependencies

The executable independently counts complete Hamilton cycles by subset dynamic programming, including the final closing edge. It checks the physical collapse on all 125 canonical cyclic child triples at three literal queries, and on larger arbitrary-quotient child controls. It checks every relation(4) against the exact small-model extraction data in both frozen affine packets and retains the odd-total modular-loss boundary.

Dependencies are the complete five-type normalization, the two frozen affine reductions, and mod8_transfer/PROOF.md. No existing frozen proof or certificate is modified.
