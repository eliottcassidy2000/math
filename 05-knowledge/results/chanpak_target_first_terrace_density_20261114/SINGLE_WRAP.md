# One long reversal: the exact joint readout

This records a physical route beyond the first terrace, not a claim that
the wrap operation is new. Let (B) be any nonempty tournament. Put
singleton (u) before (B) and singleton (v) after it in an ordinal
sum, then flip (u\to v) to (v\to u). The result is strong: every
middle vertex lies on (u\to b\to v\to u). Define
[
 L_B(z)=sum_{S\subseteq V(B)}z^{|S|}H(B[S])H(B[V(B)\setminus S]),
 qquad H(\varnothing)=1.
]
Then the exact ordinary path count of the wrap is
[
 H(W(B))=H(B)+L_B(1).
]

**Proof.** In the odd-cycle packing expansion, packings that avoid the
flipped edge give (H(B)), since neither singleton can contribute.
A packing using it has exactly one cycle through it, determined by a
directed path on its middle vertex set (S); oddness of the cycle
requires (|S|) odd. Their contribution is (L_B(1)-L_B(-1)).
Applying the same expansion to both factors in (L_B) gives
[
 L_B(z) = sum_{\mathcal C} 2^{|\mathcal C|}
   (1+z)^{n-|V(\mathcal C)|}
   \prod_{C\in\mathcal C}(1+z^{|C|}),
]
over vertex-disjoint directed odd-cycle packings in (B). Each
cycle is assigned wholly to one side; each uncovered vertex is free.
At (z=-1) every summand vanishes, since (B) is nonempty and every
cycle is odd. This proves the formula.

For an ordinal sum (B=B_1\to\cdots\to B_k), both (H) and
(L(z)) multiply, since every induced path visits the blocks in
order. Write (h_i=H(B_i)), (ell_i=L_{B_i}(1)). The wrap therefore
turns independent blocks into the exact physical count
[
                          \prod_i h_i+\prod_i \ell_i .
]
The second product is not determined by the first: a directed
triangle and that triangle padded by two ordinal singleton blocks
both have count 3, but (L(1)=12) and (48), respectively.
A triangle and a transitive block of (r) vertices give 3 and
(12\cdot2^r) as the two products. As a separate family, (k)
triangles give (3^k+12^k), with one connecting reverse arc turning
the ordinal sum into a strong tournament.

For any (n\ge2), (4\mid L_B(1)). For even (n\ge4),
(8\mid L_B(1)). Indeed the summand in the packing expression
has exponent (n-|V(\mathcal C)|+2|\mathcal C|); check zero, one,
and at least two cycles. When (n\ge3) is odd, (L_B(1)/4)
modulo 2 is precisely the parity of the number of Hamiltonian
cycles of (B), since only a single spanning cycle contributes
exponent two. Thus the readout sees more than (H(B)), but it
still has residue constraints. This narrows the prospective
target-first problem to factoring (N) additively into two products
built from the **same physically attainable pairs**, rather than
selecting the products independently. No such general selection is
proved here.
