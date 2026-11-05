# Independent review of the nonsymmetric rank-one boundary

**PASS, no correction.** This review binds `session197-prime-atoms/rank-one-boundary/PROOF.md` SHA256 `b7935c6f082ff3b6840b71e414a109a15509da3daa8fe91f336796cbac8ebc3b`, standalone verifier `642b46e03623d7243125cf8721591e518e5f677cd58e2ee8bab74717f1ec055f`, and frozen receipt `bb59884ed3dbbd3a24b9647c7aceec8f031462a35ed722622c7b919b408194f4`. No source file was changed.

## All-order arguments

For a fixed exterior ordering, the root-gap expansion correctly has separate initial and final factors and one two-choice factor per internal gap. At z=-1 the endpoint factors force an I-to-O word. Within a side the transition is unchanged. From I to O it becomes the negative reversed O-to-I incidence, while from O to I it remains B. Thus the entire signed Hamilton-word identity, including the a0 term, follows without an assumed exponential cancellation. The singleton case has the independently specified polynomial 1+z and zero boundary.

A nonempty bipartite graph of matching rank one is a star. If a chosen cover center is in O, a nonzero signed Hamilton word can enter O only at that center. It cannot later return to I and still finish in O without revisiting the center. Therefore it crosses exactly once and must use complete side Hamilton paths on either side of that crossing. Its negative sign and the two endpoint factors are exactly those in the proof. The dual argument for a center in I has the same sign. Both factors may vanish, and multiple valid centers give the same signed sum rather than conflicting formulas. Rank zero has no signed crossing and hence zero boundary.

The binary-clock interpretation is consistent with these definitions. A two-root word exists by starting from r,o,i,r for any crossing edge o->i and inserting every remaining vertex into that directed word; the ordinary tournament path-insertion argument also works with the two nonadjacent copies of r. Three root occurrences would require two disjoint O-to-I exterior fragments, impossible for a one-vertex cover. Hence a2 is positive and higher run coefficients vanish. The Stirling formula gives A=a2 and beta=a1-2a2, so H(D)=A+beta-K has the stated sign. This is a deletion identity, not an evaluation at a negative number of clones.

Both rank-two hostile examples have the claimed matching rank, complete run polynomial and opposite boundary signs. The result therefore does not imply a sign-definite rank-two boundary, a zero-boundary condition for unit ternary rays or any unproved unit-offset congruence.

## Independent exact controls

The independent `check.py` imports no source verifier. It enumerates every exterior permutation and every literal subset of root gaps, constructing and checking the resulting word rather than using a gap DP. A separate literal permutation sum evaluates the signed matrix directly. Side endpoint counts are also enumerated as literal side Hamilton words. Matching rank is computed by an independent exhaustive matching recursion.

- All 5,405 rooted labelled tournaments through order five are checked, plus 100 deterministic rooted controls at orders six and seven, for 5,505 roots in total.
- Every available star center is checked, yielding 4,975 exact endpoint factorizations.
- The controls contain matching ranks zero through three; singleton and empty-side boundaries are included.
- Both displayed rank-two run vectors and boundary values are reproduced exactly.

Independent normal and Python -O receipts agree. The supplied standalone C++ verifier was independently compiled in optimized and unoptimized assertion-enabled modes; both receipts match the frozen source output exactly, including its broader labelled/unlabelled and seeded universes. Its exact degree-class canonicalization is consistent with exhaustive one-vertex extension, and its stated signed 64-bit bound is valid at the tested orders.

The finite controls verify the implementation; the signed gap expansion and one-crossing decomposition establish the all-order conclusions. This review does not include any later generalization of the rank-two boundary.
