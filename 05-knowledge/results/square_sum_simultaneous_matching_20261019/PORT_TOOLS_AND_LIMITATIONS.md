# Elementary port tools and what they do not prove

## Two-cut path insertion

Let a square-sum Hamiltonian cycle contain two oriented edges a→b and c→d with four distinct endpoints. Let P=(p1,…,pr) be a square-sum path on new vertices disjoint from the host cycle. Suppose

    a+p1=u²,  pr+c=v²,  b+d=w².

Delete the two host edges and insert a–P–c and b–d. Starting at a, traverse P to c, then the old b-to-c segment backward, then b–d, and then the old d-to-a segment. This is one Hamiltonian cycle on the union of the vertex sets. Any previously protected edge remains if neither deleted edge was protected. This is elementary two-cut surgery, not a claimed new general extension theorem.

The conditions are sufficient and explicitly checkable. They do not guarantee that suitable host ports exist for an arbitrary P. To extend an integer interval by this operation alone, P must use exactly its new interval vertices.

## Affine descent criterion for one retained path segment

Let Q=q² be an odd square and h=(Q−1)/2. A path segment (p1,…,pr) can be represented as a lift of a path with first-position parity ε if and only if there is c∈[−h,h] such that

    p_j ≡ (−1)^(ε+j−1)c (mod Q)

and the quotients (p_j−(−1)^(ε+j−1)c)/Q are positive and pairwise distinct. If the original adjacent p sums are squares, the quotient adjacent sums are squares as well: Q divides each square sum, so q divides its square root. The centered residue chooses c uniquely once ε is fixed.

This is an algebraic representability criterion, not an existence theorem for a containing interval Hamiltonian cycle. The quotient path must still be present in a suitable input family, and the finite external-port glue states must exist.

## Stable prefix splitting

For a nice pair with a common fixed prefix, cuts within that prefix give pieces whose common-vertex membership agrees between the two input members. The endpoints of each piece have known original positional parities. A finite pair-glue certificate can therefore assign a separate parity shift to each stable piece and require that piece's shift to agree across the two output words. Shifts need not agree between different pieces of the same original block. Requiring one shared shift across all pieces, as our certified repairs do, is a stronger sufficient condition.

For arbitrary internal cuts whose orientations or induced vertex partitions vary between the input members, this stable-membership hypothesis must be proved separately. It cannot be inferred merely because both paths contain the cut edge.

## Current proved existence versus open extension

The exact seeded threshold theorem and the simultaneous matching corollary have complete existence proofs. The elementary insertion lemma and affine descent criterion only identify sufficient local conditions or representable pieces. They do not establish that every finite square-sum path eventually embeds, nor that every seeded Hamiltonian cycle can be extended. Such a theorem would need an additional port-availability or universal finite-glue existence argument.
