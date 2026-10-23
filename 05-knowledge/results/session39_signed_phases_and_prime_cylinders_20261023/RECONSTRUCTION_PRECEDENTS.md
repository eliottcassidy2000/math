# Reconstruction precedents and the precise observer being separated

Status: literature reconciliation, with the elementary implications below stated separately. No claim of external priority for the Session39 construction is made.

## 1. Classical deck counterexamples already defeat unrooted counts

Stockmeyer's classical tournament counterexamples have equal full vertex-deleted decks despite nonisomorphism. Use Kocay's later constructive proof alongside the original paper. McKay records the families at orders 2^t+1 and 2^t+2, further extensions, and the self-complementary (self-converse for tournaments) odd-order family.

- P. K. Stockmeyer, 1977: https://onlinelibrary.wiley.com/doi/abs/10.1002/jgt.3190010108
- W. L. Kocay, 1985: https://onlinelibrary.wiley.com/doi/abs/10.1002/jgt.3190090406
- B. D. McKay, 2022, pp. 452–453: https://ajc.maths.uq.edu.au/pdf/83/ajc_v83_p448.pdf

The implication for unrooted motif counts is elementary. A particular induced k-vertex copy appears in exactly n−k vertex-deleted cards. Therefore the sum of its counts over all cards is (n−k) times its count in the whole tournament. Equal full decks imply equal unrooted k-censuses for every k<n.

Consequently, it would be misleading to advertise a linear lower bound for unrooted motif cutoffs as a new general reconstruction obstruction. Self-converseness alone does not supply such novelty either.

## 2. Ordinary refinement is an additional, substantive observer

Egrot–Hirsch identify the Stockmeyer tournament pairs in their discussion of D_(m,n) and D*_(m,n). Their Lemma 7.11 gives different tally spectra. See §§4 and 6–7 of:

https://arxiv.org/pdf/2008.01327

Tally refinement retains incoming and outgoing counts inside the vertex's own current equivalence class. Ordinary directed color refinement counts neighbors in every current class. By induction the latter refines the former: when a tally class is a union of color classes, its counts are sums of the individual color-class counts. Thus the different tally spectra for the D_(m,n), D*_(m,n) pairs treated there imply separation by ordinary unindividualized directed color refinement. This last sentence is the stated elementary inference, not a quotation of a separate theorem from the paper.

The classical deck counterexamples therefore do not automatically subsume the conjunction established here: pointwise equality of static singly-rooted motif profiles, blindness of every round of ordinary refinement initialized by those profiles, and repair by one actual distinguished vertex.

This focused comparison does not establish that no other published construction has the same conjunction. The contribution is presented as an explicit proved carrier and observer separation, not a priority assertion.

## 3. Four notions that must remain distinct

1. An unrooted k-deck is a multiset of induced k-vertex isomorphism types.
2. A static vertexwise rooted profile counts root-preserving induced types at each vertex, compared under one vertex correspondence.
3. Indexed k-hypomorphy requires T[X] and T'[X] to be isomorphic for every specified k-subset X under a common underlying vertex correspondence.
4. Individualization gives one actual vertex a distinguished color that is retained during subsequent refinement.

Equal counts do not provide the indexed correspondence in item 3. Giving every vertex its own rooted-count vector is not item 4.

Lopez's {1,…,6}-reconstructibility theorem concerns the indexed-hypomorphy notion. A modern primary research statement and the original references appear in Achour–Boudabbous–Boussaïri (2012), pp. 433–434 and references 9–11:

https://www.numdam.org/item/10.1016/j.crma.2012.04.022.pdf

The Session39 pairs are not asserted to be indexed ≤k-hypomorphic. Indeed, in the cycle-family notation, the same three vertices

    {(0,0),(0,1),(g−1,0)}

form a transitive triple in the +1-phase tournament and a cyclic triple in the −1-phase tournament. Their equal rooted counts arise from root-fixing bijections between collections of subsets; those bijections depend on which cells occur. They do not fix each indexed subset.

Also, “strongly self-dual” in some reconstruction literature means that every induced subtournament is self-dual. Here “strong and self-converse” means **strongly connected and isomorphic to its converse**. No hereditary self-duality is asserted.

## 4. Safe statement of the result

We give an explicit family at every order n≥11 of strongly connected, self-converse tournaments witnessing a linear-size obstruction for ordinary unindividualized color refinement seeded by vertexwise singly-rooted induced-motif profiles, together with a one-vertex individualization repair.

This is an observer-specific construction. It neither contradicts indexed reconstructibility nor claims a new general obstruction to tournament reconstruction from arbitrary information.
