# Complete weight frontier for the provisional elimination template

**PROVED algebraic classification of a specified proof template.** This does not certify the candidate's circuit theorem, explicit disperser construction or advertised coefficient-four record. The challenge results were marked provisional when inspected on30 October2026.

The pinned [candidate, Section7](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/blob/4a61f024a4993592d97b49ca7d00979029b7fece/submission/elimination.md) uses a potential W=g s+i I+p q. Its worst-case row charges require g>0, i,p>=0 and

 i<=2g+p,
 t<=2g+2i-p,
 t<=3g+2p,
 t<=3g+i+p.

All other displayed rows are nonincreasing automatically for these nonnegative weights; the GP row with one child is its strongest member. If the underlying elimination invariant and coverage theorem are valid, telescoping these rows yields leading coefficient c=(t-i)/g and dimension penalty b=t/g. We classify the rowwise certificates for c>=0; we do not assume the unverified structural theorem in proving the following algebra.

## The exact feasible region

Normalize g=1. A pair(c,b) is feasible if and only if

 **0<=c<=4 and b>=max{c,2c-2,4c-7}.**

Every feasible pair has i=b-c,t=b and the complete allowed interval of p is

 **max{0,b-c-2,(b-3)/2,c-3} <= p <= b-2c+2.**

**Proof.** Substitute i=b-c,t=b into the four inequalities. Nonnegativity of i gives b>=c. Comparing each of the four lower bounds for p to its upper bound gives respectively b>=2c-2, c<=4, b>=4c-7 and b>=3c-5. The last is redundant: for c<=3 it is bounded by2c-2, and for c>=3 it is bounded by4c-7. These comparisons are necessary and sufficient for the displayed interval to be nonempty. Any p in the interval realizes the desired coefficients. This proves both directions and the complete weight classification.

The sharp lower envelope of the dimension penalty is therefore

 b_min(c)=c for0<=c<=2;
 b_min(c)=2c-2 for2<=c<=5/2;
 b_min(c)=4c-7 for5/2<=c<=4.

Explicit boundary weights are

 (i,p,t)=(0,0,c) on the first segment;
 (i,p,t)=(c-2,0,2c-2) on the second;
 (i,p,t)=(3c-7,2c-5,4c-7) on the third.

At c=4,b=9, the interval forces p=3 and i=5,t=9. Thus the quoted(1,5,3) potential is simultaneously optimal for coefficient and dimension penalty within this rowwise nonnegative template. Merely retuning these weights cannot exceed4 or reduce the coefficient-four penalty below9.

## What a genuine improvement must change

The ceiling follows already from the paired-cleanup and one-child-GP rows. To cross it one must improve a structural deletion estimate, enrich the state/potential, or prove a global incompatibility that prevents the individually worst rows from accumulating. The feasible-region calculation does not prohibit such improvements. It is not an upper bound on circuit complexity, a natural-proofs barrier, or a conclusion about all gate-elimination arguments.

The [source/model report](../SOURCE_MODELS.md) separates the current primary preprint baseline from provisional submissions. The restricted source inspection found no concrete counterexample to the candidate elimination lemma, but did not independently verify its disperser/family assembly or certify a new record.

`check_potential_frontier.py` uses exact rational arithmetic to verify the projection and all boundary witnesses. The displayed algebra proves the continuum statement; the grid is an independent regression test.
