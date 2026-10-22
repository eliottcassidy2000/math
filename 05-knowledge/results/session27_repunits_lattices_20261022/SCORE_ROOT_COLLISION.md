# A single E8 score root contains different Hamiltonian counts

Status: PROVED construction and FINITE-EXACT path counts, independently checked by full permutation enumeration. This is not a new classification of E8 or of regular tournaments.

Let Q be the Paley tournament on Z/7Z, with u->v when v-u belongs to {1,2,4}. The three vertices 0,1,3 form the directed triangle 0->1->3->0. Let Q' reverse precisely these three arcs. Both tournaments are regular: each triangle vertex loses one outgoing arc and gains one, and every other score is unchanged. Exact Hamiltonian counts are H(Q)=189 and H(Q')=171.

Adjoin a vertex labelled7 to both tournaments, with core vertex u pointing to7 exactly for u in {0,1,2}. The resulting T,T' have the identical labelled score vector

    (4,4,4,3,3,3,3,4),

but H(T)=657 and H(T')=643. Their respective outgoing-neighbor bitmasks are

    T:  [150,172,216,49,98,69,11,120],
    T': [156,165,216,50,98,69,11,120].

Bit j in mask i means i->j. This completely specifies both objects, not just their unlabelled score multisets.

## Exact lattice image and the obstruction

Subtracting (n-1)/2=7/2 from every score gives the same vector

    v=(1/2,1/2,1/2,-1/2,-1/2,-1/2,-1/2,1/2).

In the standard realization E8=D8 union ((1/2,...,1/2)+D8), this belongs to the half-integral coset: subtracting (1/2,...,1/2) gives an integer vector of even coordinate sum -4. Its squared norm is8/4=2, so it is an actual E8 root, not merely an arbitrary point in the lattice.

Therefore the marked score map from eight-vertex tournaments to E8 does not determine H, even on its minimum-norm score shell, even before any quotient by coordinate permutations, and even when the change of orientation is a single directed-triangle reversal. In particular no function of this score vector alone, nor of its norm or its E8 orbit, can reproduce all Hamiltonian counts.

The construction preserves the exact degree vector and its lattice image. It loses directed cycle/attachment information. Restoring a richer response observable is not optional if H is the target. This is a concrete obstruction to the specific score-lattice map; it does not rule out other maps from tournaments to lattices or modular forms.

## Reproduction and independent path

`score_root_hostile.py` computes the counts by subset dynamic programming. The Session27 response audit independently used every permutation of seven vertices for each core and every one of the40320 permutations for each eight-vertex tournament; it recovered the stated values and score vectors without using the DP.

The script's original import-free implementation is copied into the distribution version so it has no dependency on a private workspace path. Neither a local field of a Paley construction at composite order nor an unverified general Paley maximality assertion is used.

## The first information loss already occurs at order five

The score vector determines the number of cyclic triples:

    c3(T) = binom(n,3) - sum_v binom(outdegree(v),2).

Indeed every transitive triple has a unique source, and every choice of two outneighbors contributes exactly one transitive triple. For n<=4, the odd-cycle collection formula gives H=1+2c3, since no longer odd cycle or two disjoint odd cycles fit. Hence no two tournaments of order at most four with the same score multiset can have different H.

At n=5, the exact formula is H=1+2c3+2c5. Two explicit tournaments have outgoing masks

    [12,17,2,22,5] and [10,20,1,22,5].

Both have labelled scores (2,2,1,3,2) and c3=4. The first has c5=2 and H=13, while the second has c5=1 and H=11. They differ only by reversal of their core directed triangle. Thus order five is the exact first failure of the score-vector-to-H implication. The difference is identified: degree data retain the triangle count but lose the five-cycle contribution.

The companion `score_first_failure.py` checks both H values by all120 permutations and counts cycles directly, with rotations identified and orientations retained. The odd-cycle collection identity is inherited classical theory, not newly attributed here: Grinberg–Stanley, *The Redei–Berge symmetric function of a directed graph*, and Irving–Omar, *Revisiting the Redei–Berge Symmetric Functions via Matrix Algebra*, Corollary20, as pinned by the repository's THM-002-ocf.md.
