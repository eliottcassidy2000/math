# Independent symbolic review of root's product law and class counts

Status: PASS analytically; exact finite mean/product controls and floating CDF diagnostics are separated. No literature-priority claim. This reviews root's later packet rather than claiming a new independent discovery.

## Product law

Use zero-based positions0..m-1. Set A=sum_j b_j 2^(m-1-j), B=sum_i(1-b_i)2^i. In AB/2^m the forward pairs j<i have exactly the ordinary inversion weights2^(i-j-1). Equal indices vanish. Reverse pairs contribute R with0<=R<m/2. Hence H=AB/2^m+1+#adjacent10-R and absolute error at mostm/2+1. The normalized product approximates H/2^m within(m/2+1)2^-m.

For uniform strong words, the first q=floor(m/2) binary digits of A/2^m and B/2^m use disjoint bit sets, with leading digits1. They are independent uniform dyadic grids in[1/2,1]. Couple each with an independent uniform variable in its dyadic cell; errors are at most2^-q, irrespective of the dependent lower source digits. The product error is O(2^-q). The product density is bounded by4log2, so the resulting CDF error is uniform O(2^-m/2), including the join at1/2 and endpoints.

The product CDF is4zlog(4z)-4z+1 on[1/4,1/2] and4z-4zlogz-3 on[1/2,1]. Its variance is(7/12)^2-(3/4)^4=55/2304. Exact mean of H is9*2^(m-4)+1/2 for m>=3: ordinary boundary-pair contributions are2^(m-2) and2^(m-2)-1; internal-pair mean is(2^(m-2)-m+1)/4; adjacent10 mean is(m+1)/4; add the constant1.

For X=2^k x,1<=x<2, root's exact length decomposition is correct: all words of lengths2..k contribute2^(k-1)-1, and only lengths k+1,k+2 can contribute further. The identity1+F(x/2)+2F(x/4)=2xlog2 holds exactly. Therefore the word count W(X)=Xlog2+O(sqrtX), uniformly in x.

Importantly, the distinct-H upper bound does not require the apex/isomorphism classification below. Reverse-complement word orbits already show at mostW(X)/2+O(sqrtX) values, since their fixed words total O(sqrtX) through the allowed length range. Thus missing odd values have count at least(1-log2)X/2-O(sqrtX). This strengthens the elementary dyadic bound in CUT_ADDENDUM.md. It does not establish a density limit for distinct values.

## Apex and class classification

Deletion-transitive apices form the intersection of all directed triangles; for a nontransitive tournament their number is at most3. If exactly3, there is one triangle. Every outside vertex must uniformly dominate or lose to that triangle, or it creates another triangle. The quotient must be transitive, so the only nontrivial strong core is C3.

If strong T has two deletable apices x,y, the remaining vertices form a unique transitive chain. Each apex inserts at a cut in that chain. Suppose x->y. Cyclic triangles are i->x->y->i, exactly at indices between y's cut and x's cut. Strongness forces y's cut0 and x's cut equal the chain length. Hence T=C3(1,1,n-2). Conversely this construction has exactly two such apices at n>=4.

An automorphism fixing an apex fixes the unique transitive deletion order pointwise; thus the action on apices is free. Tournament automorphism groups have odd order (an involution would reverse its swapped pair). For strong n>=4 the apex count at most2 forces rigidity. The exceptional C3 has automorphism order3. Therefore the rooted-word duplicate accounting and root's class formulas follow.

For strong n>=4, ordinary classes number2^(n-3)-1. Including source/sink padding gives2^(n-1)-binom(n-1,2)-1 ordinary classes at n>=2. Strong selfconverse classes are1 for even n>=4, and2^(r-1)+1 for n=2r+1>=5: the unique two-apex class is selfconverse but its two marked words are interchanged, while all other selfconverse classes are exactly reverse-complement-fixed words. Summing with equal source/sink padding gives n/2 for even n>=2 and2^r+r-1 for odd n=2r+1>=3. All arithmetic sums check.

## Controls and consequence

Run `node verify_asymptotic.js`. Exhaustive strong words at lengths3..18 verify the exact mean and pointwise dyadic error. The five-threshold CDF differences involve floating logarithms and are diagnostics, not exact proof gates. The all-order claims rest on the symbolic coupling and counting proofs above.

Partial summation also gives a meromorphic continuation of the WORD-weighted Dirichlet series to Re(s)>1/2 with sole pole at1, residue log2. The error integral is holomorphic there. This applies to words (and, using classification, ordinary classes), NOT to the unknown distinct-H-value counting series.
