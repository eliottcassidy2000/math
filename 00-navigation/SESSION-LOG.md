## boxeph-2026-07-21-S217 -- class number = arithmetic entropy: hidden binary forms, and why 7 is rigid (HYP-8870)

**Owner:** look for even more hidden binary forms; think information theory.

**HIDDEN FORMS:** the disc -p anisotropic gate (S216) sits in a whole CLASS GROUP of binary forms (Gauss composition); the NON-PRINCIPAL classes are the hidden forms. Verified: h(-3)=h(-7)=h(-11)=1 (Paley/Heegner, principal only); h(-15)=2, h(-23)=3, h(-31)=3, h(-47)=5 (hidden forms enumerated).

**INFO THEORY (verified, hidden_binary_forms_class_number_is_arithmetic_entropy_boxeph_S217.py):** CLASS NUMBER h(D) = the ARITHMETIC ENTROPY = log2(h) bits BEYOND local Legendre (D/p) to decide which form represents p. disc -7 (h=1): p represented <=> (D/p)=1 (ONE Legendre bit determines all). disc -23 (h=3): (D/p)=1 primes SPLIT among 3 classes (principal 59,101,.. vs non-principal 13,29,31,..) -- which class = the Artin symbol in the Hilbert class field Cl(D)=Gal(H/K), invisible to any local test = log2(3)=1.58 hidden bits.

**WHY 7 IS RIGID:** LRC(14)=2*7 -> disc -7 -> h=1 -> ZERO arithmetic entropy. codex THM-2053's anisotropic gate residual is fully pinned by local S215 Legendre data -- NO hidden bits, no class-group slack. So (1) a counterexample has NOWHERE to hide; (2) the certificate must be the exact local (Euler/chi/Borsuk-Ulam, p=3mod4) one. Rigidity = why 7 is the first hard-but-tractable case (kps-S17). Heegner h=1 imag. quadratics = -3,-4,-7,-8,-11,-19,-43,-67,-163; -7 is LRC(14)'s.

**BONUS:** tournament score-distribution entropy separates the reify-ladder poles: transitive (0..n-1, entropy log2 n = MAX spread, the nullcone/rank-11 gate vertex) vs Paley/regular (all (n-1)/2, entropy 0, the symmetric disc -p pole).

**Honest:** class-group/Heegner/representation facts classical + verified; the contribution is the info-theoretic FRAMING (class number = arithmetic entropy; non-principal forms = hidden forms; S216 rigidity = zero hidden entropy for -7). Conceptual sharpening, not a proof step. Ties S215+S216+kps-S17. Artifacts: reflection class-number-is-arithmetic-entropy-...-boxeph-S217.md, HYP-8870, script (+.out).

## codex-2026-07-21-LRC-normal-fan-sail -- THM-2055 and MISTAKE-225

**PROVED REPAIR:** the THM-2053 determinant is
`D(d)=h_K(Rd)` for the centrally symmetric column polygon
`K=conv{+-c_i}`. It is a polyhedral norm. Only hull vertices can own gate
failure; their rational normal cones split the residual, and each owner cone
intersects one explicit tangent disk with local radius bound
`||d||<91||p||`. On the one-tail plane the 26 nominal signed-column disks
collapse to four hull vertices `+-(13,0),+-(12,1)`.

**CORRECTION:** HYP-8865's discriminant-`-7`/Heegner-classification and
isotropic-rank versus anisotropic-Euler equivalence are retracted in
MISTAKE-225. The left side is max-linear, the right quadratic form has
discriminant `-4`, and a sufficient gate's failure is only “unresolved.” The
Paley/Heegner data has no predicate-preserving map to an arbitrary star plane.

**NEXT CREATIVE TARGET (HYP-8870):** use ordinary Stern--Brocot/Klein sails in
each rational owner cone. Split only at hull-owner ties, positivity/collision
walls, pair-sum changes, Fejer resonances, or endpoint-owner exchanges. This
keeps a genuine two-dimensional continued-fraction address without importing
the repo's refuted multidimensional-CF slogans.
