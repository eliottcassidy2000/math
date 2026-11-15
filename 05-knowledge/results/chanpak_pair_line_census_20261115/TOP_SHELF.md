# Highest source orders, without a whole high-order census

Write q=P₂+2H. At fixed order n, the pair coupling in the
positive odd-cycle expansion gives q/2^n≥1+τ/2. If the normalized
pair coordinate is strictly below 4 then τ≤5 and every nontrivial
SCC has size at most 7. The exact table for *all* such normalized
pairs at *all* n can therefore be built by taking products of
the strong profiles through seven and retaining their combined size.
Run `python top_shelf.py < profiles.tsv` with a complete profile
catalogue through seven in the same format as the order-10 census.
Uniform singleton padding produces each possible lifted order.

For target cap 59,999, this independently gives these complete
absolute sets at the top two orders of each reader:

| Reader | Source order | Hits at most 59,999 |
|---|---:|---|
| 6q−H | 13 | 49151 |
| 6q−H | 12 | 24575, 36861, 49147, 55287 |
| 3q−H | 14 | 49151 |
| 3q−H | 13 | 24575, 36861, 49147, 55287 |

The first coincidence is structural: doubling q by an ordinal
singleton and halving its coefficient gives exactly the same target.
More broadly (3q−H) on the padded source always duplicates
(6q−H) on its original source. This is one reason not to count
apparent coverage gains solely by the source orders on which the
same target reappears.

## A template-independent ceiling at this cap

A cyclic single-port reader means an ordinary tournament K with a
distinguished vertex x on a triangle, replacing x by an arbitrary
nonempty tournament A via uniform arcs. By physical path compression
its Hamiltonian count is `sum_j a_j P_j(A)`, with a_j nonnegative
integers and a₁=H(K). There is a walk `x,o,i,x` from a triangle
through x; in the blow-up with two copies of x along a transitive
edge, insert any other vertices outside the port one at a time in
this path by the standard tournament path insertion lemma. The
relative order x,o,i,x survives, so a₂≥1. Also a₁≥3 from the triangle
term in the positive Hamiltonian expansion. Therefore

    N ≥ 3H(A) + P₂(A) = H(A) + q(A).

Equality is attained by the triangle template (one compressed
two-visit pattern, no third visit). Thus if N≤59,999 then n≤15,
even with no size bound on the template and with arbitrary higher
visit coefficients. This ceiling is attained with n=15 and A
transitive. At n=15, q/2^n < 1.832, whence τ≤1. The only
source pairs are `(H,q)=(1,32768),(3,49152)`, for transitive A
or one C₃ with singleton components. In either case a₂=1 since
twice P₂ alone is too large. Also a_j=0 for j≥3. The visit support is an initial interval by the
[port matching theorem](../chanpak_reader_retention_20261115/MATCHING_AND_RETENTION.md),
so any such a_j>0 would force a₃>0, and the number
`P₃(A)≥3^15−3·2^15+3` of surjective assignments already exceeds
the cap. The reduction therefore forces the reader to be effectively
two-visit on these maximal-order sources, but it does not characterize
which a₁ are physically possible for a₂=1.
