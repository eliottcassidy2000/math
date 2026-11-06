# Independent audit of the one-sided multiplier classification

6 November 2026. **PASS in the stated one-sided context class.** This is not a spectrum-occupancy or threshold theorem.

The claim under review is that a fixed rooted core G can satisfy H(G*U)=a(H(U)+L(U))-start_root(U) for every rooted nonempty tournament U if and only if a=1 or a=2^m+1 for m>=1. Cross edges from the nonroot core to the nonroot input are all forward. Arbitrary mixed-mask contexts are not covered by this classification.

## Exact proof review

The inherited ordered-root-amalgamation formula gives the four coefficients of h,L,s,d exactly as stated, including the singleton input boundary. The four test profiles are correct and have nonzero determinant, so a universal scalar identity on physical rooted tournaments forces all four coefficients; there is no assumption that arbitrary coordinate tuples are realizable.

The forced core profile is (h_G,e_G,R_G,d_G)=(2a-1,a-1,2a,a). For a>1 it has zero merge slack. The inherited all-order zero-slack theorem therefore makes the root's outdegree exactly1: outdegree0 would give h_G=d_G. Let c be its sole out-neighbor and A=G-root. Removing the terminal root identifies e_root(G) with the paths of A not ending at c, so e_c(A)=1. The exact near-sink insertion formula then gives R_c(A)=a.

The second nonnegative-slack inequality is 0<=1-H(A-c). Its conclusion H(A-c)=1 is valid even when A-c is empty, by the stated convention. The slack is then zero. The inherited unit-deletion zero-slack spectrum applies to the actual rooted A, proving a=2^m+1. This is an all-order necessity argument, not a finite core search.

The positive core E(2^m+1), with its original marked near-sink root, has the required profile. The explicit binary word is1 followed by m zeroes, so its order is m+3. Its deletion count is a, its endpoint count is a-1, and zero slack gives R=2a. The edge handles a=1. Thus the classification has explicit physical representatives.

For a=3^k, the elementary equation3^k-1=2^m allows only k=1,2 when k>=1. The odd-k modulo8 argument is correct, and the even-k factorization gives two powers of two differing by2, necessarily2 and4. Hence the required a=27 core profile is excluded at every order. This is not exclusion of scalar H=53 and does not exclude mixed-mask alternatives or a different recursive state.

The formula on unit-start inputs gives the stated n to (2^m+1)n+2^m map, but does not establish a valid next source representation. The manuscript keeps this nonclosure explicit.

## Independent implementation

check.py imports no producer script. It uses endpoint subset dynamic programming and explicit subset convolution for the two-path responses. It independently constructs the six cores with a=1,3,5,9,17,33 and checks450 literal amalgamated output graphs against75 rooted labelled inputs per factor, covering all labelled tournaments through order4 with the root fixed at0. It also computes the four basis profiles and their determinant, and checks the inherited zero-slack characterization on1,099 rooted labelled inputs through order5 as an implementation control.

Normal and optimized runs are byte-identical. These finite controls are not the basis of the arbitrary-order exclusion. The load-bearing inherited structural statement is Section3 and Theorem3 of the frozen zero-slack proof; I read and checked the actual hypotheses before applying it.

The final packaged manuscript differs from the originally reviewed note only in inherited path normalization and replacing the submission note with the verifier scope and reproduction commands. Both versions are preserved in the audit. No mathematical repair was required. Status-only promotion is approved after preserving the reviewed original.
