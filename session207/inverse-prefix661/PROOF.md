# Session 207 continuation: complete inverse table through 661

Status: FINITE-EXACT inverse classification through count 661; PROVED construction bounds and structural reductions below. No new all-order completeness or strong-minimality theorem is claimed.

## Exact inverse classification on a full prefix

The independently reviewed `session207/size-spectrum/census.json` completely determines spectra S_1 through S_8. For each h in S_8, its first appearance is therefore its exact minimum m(h). The 320 counts in S_8 leave only nine allowed odd values at most661 unaccounted for:

    611,615,617,619,623,625,635,647,655.

`fill.cpp` finds nine-vertex witnesses for these nine counts by augmenting the previously published eight-vertex witnesses. `fill.out` records target count, nine-vertex adjacency code, eight-vertex base code, and new-vertex outgoing bitmask. The code convention is one bit per lexicographic pair i<j, with bit1 meaning i beats j. The exact codes are:

    611:45937614978, 615:45685416128, 617:45669179456,
    619:45115523264, 623:43805851840, 625:45650829504,
    635:45652394176, 647:45635608768, 655:45685956800.

`check.py` independently enumerates all9! permutations of each tournament, checks the claimed count and strong connectivity, and combines the exact previous spectra with these witnesses. Consequently each of these nine counts has m(h)=9, and `minimum_orders.json` gives the exact first order for every attainable h in [1,661]. There are329 such counts. All even counts and7,21 are impossible; the latter exclusions use the inherited global theorem, not this witness search. For each listed h, the complete inverse set is every integer n>=m(h), by universal-source padding.

This is not an exhaustive nine-vertex census. Exhaustiveness through eight provides the lower bound; nine explicit witnesses provide the matching upper bound.

## A genuinely infinite construction, and its minimality boundary

For k>=1 take a transitive chain A of k vertices and two vertices b,c with A->b->c->A. The resulting strong tournament has k+2 vertices and exactly2^k+1 Hamiltonian paths. There are three paths with A unsplit: A,b,c; b,c,A; c,A,b, where the chain ordering inside A is unique. Every other path is A_left,b,c,A_right, where A_left is any nonempty proper subset of A, both pieces in chain order. There are2^k-2 such subsets. Thus

    m(2^k+1)<=k+2.

The construction is strong because every chain vertex reaches every other via b,c whenever its direct chain direction is unavailable. The bound is an equality for k=1,2,3,4 by the exact spectra. It already fails at k=5: m(33)=6<7; at k=6, m(65)=7<8. Efficient construction is not a minimality theorem. This is the classical nearly-transitive tournament NT_(k+2): order the vertices as c, then the chain A, then b, and reverse just the first-last arc to obtain b->c. The formula H(NT_s)=2^(s-2)+1 is inherited from the Moon/Busch attribution in THM-4744; no novelty is claimed.

## Prime-power exponents as an additive cost problem

Fix an odd prime p. Any SCC factor of a tournament with H=p^k has count p^j. If c_p(j) is the smallest strong order realizing p^j, or infinity, then

    m(p^k)=min {sum c_p(j_i): all j_i>=1 and sum j_i=k}.

This follows from the proved strong-atom recurrence; singleton SCCs are omitted. It is a shortest additive composition problem in the exponent, not a claim that the prime p alone generates optimal witnesses.

For p=3 the exact table gives m(3^j)=c_3(j)=3,5,6,7,8 for j=1,2,3,4,5 respectively. Thus, writing k=5q+r with0<=r<5, dominance joins give the uniform infinite bound

    m(3^k)<=8q+d_r,    d_0=0,d_1=3,d_2=5,d_3=6,d_4=7,

for k>=1. This improves the naive3k triangle-stack bound substantially. It is only an upper bound once k>5; larger strong atoms may be more efficient.

More generally, for any attainable h>1, subadditivity m(h^(k+l))<=m(h^k)+m(h^l) implies that the limit of m(h^k)/k exists and equals inf_k m(h^k)/k, by Fekete's lemma. This gives a precise asymptotic efficiency invariant. Positivity, attainment of the infimum, or eventual periodic linear behavior is not established. The finite3-power witnesses show only that this limit is at most8/5 for h=3.

## Strong first occurrence: a precise unresolved question

Every attainable h<=661 has a strong witness at its first possible order: this is checked against the complete strong spectra through eight and the nine strong witnesses above. Therefore m(h)=a(h) throughout this finite range, where a(h) is the minimum strong order.

Question, not theorem: does every attainable h>1 have a strong witness on m(h) vertices? An affirmative answer would mean that dominance products never create a count earlier than strong tournaments do. The SCC product recurrence alone does not prove this. Composite atoms and the possibility of cheaper decomposed witnesses must remain in the exact algorithm until this stronger statement is established.

## Reproduction and scope

The package bundles an unchanged copy of the previous exact census as an explicit dependency (`census.json`), from the reviewed `session207/size-spectrum/` release at commit d3b91e52c6b756d32ab90fe7f1a406f23633f60d. `inputs.txt` lists its order-eight witness codes. Compile fill.cpp and run with inputs.txt to reproduce fill.out; then run check.py normally (assertions enabled) to reproduce check.out and minimum_orders.json. The checker is deliberately independent of the C++ forward subset DP: it uses direct path permutation enumeration for all new witnesses. No random search, sampling assumption or nine-vertex exhaustiveness assertion enters the result.
