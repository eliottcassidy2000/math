# A single profile bit per even child controls coarse responses modulo sixteen

## 1. Complete profile reduction

Fix any named tournament quotient Q on m nonempty child sets. Every actual child A_i has even order at least four. Write H_i=q1(A_i), which is odd. Normalize its complete profile by H_i modulo16. The profile-domain theorem in `../even_cyclic_mod16/PROOF.md` gives

    q(A_i)/H_i=(0,1,b_i,c_i,8,d_i,0,...),
    b_i in {6,14}, c_i in {4,12}, d_i in {0,8}.

The actual transitive four-name tournament has residue profile

    e=(0,1,14,4,8,0,0,...).

Define epsilon_i=0 when b_i=14 and epsilon_i=1 when b_i=6. The other two profile bits will have no effect on the complete block-uniform Hamilton response modulo16.

For a physical block query (a,c), let F_4(Q;a,c) be the actual response of Q inflated by TT4 in every child, with the block masks expanded to those four-name sets. Let D_i(Q;a,c) count directed run words on the quotient labels and the two actual ports, in which i occurs exactly twice, every other old label exactly once, each port exactly once, and no two consecutive labels are equal. All transitions follow the specified quotient/port arcs, including x->y.

**Theorem.** For the complete actual child-uniform response,

    F_(Q[A_i])(a,c)
      = (product_i H_i)
          [F_4(Q;a,c)+8 sum_i epsilon_i D_i(Q;a,c)] (mod16).   (1)

This holds for every quotient and all block masks, without strong-connectivity or primeness assumptions. It is an identity of the entire coarse response bank, not a quadratic/cubic approximation.

## 2. Why the other two bits disappear

The exact run-word functional is multilinear in the complete child profiles, one profile factor per child. Expand each normalized profile around e. Each perturbation is 8 in coordinate2,3 or5. Products of two perturbations vanish modulo16. In one first-order perturbation term, only the parity of every other child factor matters. Only e1 is odd, so every other child must occur exactly once in a contributing run word.

Suppose the marked child occurs k>=3 times. Write the word as a prefix, its first occurrence of i, an intermediate segment, its second occurrence of i, another intermediate segment, and the remaining suffix. The two displayed return segments are nonempty: a Hamilton run word has no i->i transition. Their label supports are disjoint, because every other child and each port occurs at most once. Swap these two entire segments, including their joins between consecutive i occurrences. The directed joins remain legal because both return segments have the same i endpoints. Prefix, suffix, label multiplicities and all physical port uses are unchanged.

This is a fixed-point-free involution: the two nonempty segments have disjoint label supports and cannot be identical. It proves that the run coefficient is even. In particular the coordinate3 and coordinate5 perturbations disappear after multiplication by8. The coordinate2 perturbation is exactly D_i. Since +8 and -8 coincide modulo16, epsilon_i gives its required sign, proving (1).

This uses the previously established three-run segment-swap idea from `../parity_profile_transfer/PROOF.md`; the new implication here is its first-order modulo-sixteen profile reduction. The proof is for Hamilton run words. It must not be applied unmodified to ordered covers: a cut between identical child labels can make a return segment empty.

## 3. Small physical computation of the correction

Let F_1(Q;a,c) be the actual response with every quotient child a singleton. Let F_(2,i)(Q;a,c) instead inflate only i by a transitive pair, with its two port incidences uniform. All other children remain singletons. Their exact profiles show, over the integers,

    F_(2,i)(Q;a,c)=F_1(Q;a,c)+2D_i(Q;a,c).                  (2)

Indeed a singleton contributes only q1=1, while the pair contributes q1=1,q2=2 and no higher coordinates. Thus D_i can be obtained exactly from actual tournaments with m+2 or m+3 total names after adding ports. This supplies an independent physical check of the derivative coefficients.

Combining (1) for Q and Q^op gives the quotient-converse contrast

    Delta_actual
      = (product_i H_i)
          [Delta_4+8 sum_i epsilon_i(D_i(Q)-D_i(Q^op))] (mod16). (3)

A contrast that is already nonzero modulo8 in the TT4 baseline cannot be cancelled by any allowed profile bits. If the entire TT4 coarse bank is zero modulo8, divide its modulo16 contrast by8. Complete modulo-sixteen blindness then becomes the finite affine linear system over F2

    Delta_4(a,c)/8 + sum_i epsilon_i (D_i(Q;a,c)-D_i(Q^op;a,c))=0

for every block query. It has only m unknown bits, rather than3m. An inconsistent system proves coarse separation for every actual all-even-at-least-four inflation of that fixed quotient. A solution establishes only formal modular blindness; it is neither an exact collision nor proof of actual profile realizability.

## 4. Scope

The finite arithmetic domain has been reduced from 8^m to2^m for every quotient. No claim is made that the remaining linear system is inconsistent for all quotients. The observer is child-uniform; no equivalence of full named responses is asserted. No arbitrary ambient lift follows from this Hamilton-only reduction.

The normalization, the run-word functional and the older three-run involution are credited dependencies. Any ensuing quotient classification requires its own exact certificate and independent audit.
