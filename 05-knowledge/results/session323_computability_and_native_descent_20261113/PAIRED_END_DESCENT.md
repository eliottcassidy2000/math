# Exact target descent through arbitrary one-apex dyadic reservoirs

This extends the [finite factor-bank inverse](../session322_target_inverse_prime_channels_20261113/DYADIC_INVERSE.md) to an infinite physical grammar. It is complete for that grammar, not for the entire tournament spectrum.

Fix either of the two selective quotient charts and a literal nonempty source U. Let h=H(U) and A=−3h+6Q_U(2)+cQ_U(3), with c=2 or 6 according to the chart. Here Q_T(z) counts ordered covers by z paths with empty slots allowed. Allow the reservoir V to be an ordinal join of any number of strong one-apex-over-chain blocks, singletons, and, optionally, any specified finite bank of already certified dyadic blocks. Require Q_V(2) to be a power of two. The output is strong and its exact count is

    N=A·2^m−hD,  m=log₂ Q_V(2),  D=H(V).

The compiler below starts from N, constructs V if possible, and rejects exactly when there is no witness in this grammar.

## Paired-end arithmetic for a strong core

In a strong one-apex core the chain has length n≥2 and the apex's cut is a binary word x_1…x_n with x_1=1, x_n=0. These two conditions are necessary and sufficient: otherwise the first vertex is a source or the last a sink; with them the whole chain and apex lie on a directed Hamilton cycle. Let M_d count pairs i<j at distance d with x_i=1,x_j=0. Odd cycles all contain the apex, and selecting an even number of intervening vertices gives

    Q(2)=2^(n+1)+Σ_{d=1}^{n-1} M_d 2^(n−d)(3^(d−1)+1),
    H=1+Σ_{d=1}^{n-1} M_d(2^(d−1)+[d=1]).

These are the literal cycle-packing formulas used in the earlier [residue aperture](../session322_target_inverse_prime_channels_20261113/RESIDUE_APERTURE.md); here the maximal spans are inverted first.

There is exactly one maximal-distance pair. If n is even, it gives 2-adic valuation 2, while every other term has valuation at least 3. Hence Q(2) is never dyadic at even n. For n≡3 modulo 4 with n≥7, the top two terms are the only ones affecting the valuation below 5. If M_(n−2) is even the valuation is 3. If it is odd it equals 1, and the sum of the top terms is 2(5·3^(n−3)+3), of valuation 4 because n−3 is divisible by four. Neither case is a power of two, since Q(2)>2^(n+1). At n=3 the two cuts 100 and 110 both give (H,Q(2))=(5,32).

For the remaining lengths n=4J+1, reveal the cut from its two ends. Put a_i=x_(i+1), b_i=1−x_(n−i), and C_k=Σ_{i=0}^k a_i b_(k−i). Then M_(n−1−k)=C_k. Pair consecutive k. For the outer disjoint shells j=0,…,J−1, their contributions are

    8·4^j(α_j C_(2j)+β_j C_(2j+1)),
    α_j=(3^(n−2−2j)+1)/4,   β_j=(3^(n−3−2j)+1)/2.

Both coefficients are odd. Terms from later shells are divisible by the next power of four after factoring 8. Since a dyadic output must have valuation at least n+2, every revealed base-four digit must vanish. This is a carry decoder: if r is the normalized carry from the previous shells, accept exactly those new four bits with r+α_j C_(2j)+β_j C_(2j+1) divisible by four, then divide that sum by four to obtain the next carry. At j=0 this requires a_1+b_1=1, as a_0=b_0=1 and α_0+β_0≡0 modulo 4.

For j≥1, the four new bits are a_(2j),b_(2j),a_(2j+1),b_(2j+1). Their coefficients modulo 4 are α_j+β_j b_1, α_j+β_j a_1, β_j, β_j. Exactly three are odd; the even one is α_j+β_j, equal to 0 at even j and 2 at odd j. Indeed

    α_j+β_j = (5·3^(n−3−2j)+3)/4,

and the exponent is 2 modulo 4 at even j, 0 modulo 4 at odd j. Complement any of the three odd bits with coefficient −1 modulo 4; this shifts the target but turns their residue counts into (1,3,3,1) over the four values. If the even bit has coefficient 0 there are at most six choices; if it has coefficient 2 there are exactly four. Thus the carry tree has O(24^(n/8)) final candidates. These disjoint shells leave precisely the central bit undecided. Try both values and evaluate the full formulas, rejecting every candidate whose Q(2) is not a power of two or whose H is not required. This gives an exact necessary-and-sufficient word compiler, not just modular tests. No claim is made that large n≡1 mod 4 always produce dyadic cores.

## Actual target-decreasing recursion

As before, A>h and D≤2^(m−1). The target forces the unique m with A·2^(m−1)<N<A·2^m, and then the unique D=(A·2^m−N)/h. Reject if these constraints fail. Every SCC's Q(2) divides the dyadic product and hence is itself dyadic. An SCC with count d>1 can occur only if d divides D.

For each divisor d>1 of D use the carry compiler to discover all permitted one-apex factors of count d and dyadic cost at most m. There are at most two candidate chain lengths per d: the extreme pair gives H≥2^(n−2)+1 and the [cut-one binary count formula](../session317_exact_anchor_atlas_20261113/ONE_RETURN_STAR.md) gives H≤2^n−1 (write H=2B+1−K(B) with B<2^(n−1) and K≥0). Also Q(2)≥2(3^(n−2)+1) gives n≤2+(m−1)/log₂3. If a finite additional bank is desired, include its verified factors. For each d keep only the least cost c(d) and a literal word or graph realizing it; higher costs need not be kept because free singleton padding pays any unused power of two.

Now use the descending divisor recurrence

    cost(1)=0,
    cost(t)=min_{d>1, d|t} [c(d)+cost(t/d)],

with an empty minimum interpreted as infinity. Each successful step physically deletes one SCC of count d and slot response 2^c(d), leaving a smaller divisor and decreasing the remaining exponent. Accept exactly when cost(D)≤m, then assemble the chosen blocks and m−cost(D) singleton padding vertices. Conversely any witness gives a legal factorization of D, and replacing a factor by a cheaper factor of the same H only makes more room for padding. This proves completeness.

On a successful path, whenever the remaining reservoir is nonempty, the scalar transition is

    N → A·2^(m−c(d))−h·D/d.

It is literally realized by the same outer cut and the remaining ordinal reservoir. The shell decoder recovers the new two coordinates from the new absolute count. For a nontrivial SCC, d≥3 and 2^c(d)≥2d, so c(d)≥3. Since N≥(A−h/2)·2^m and the new count is less than A·2^(m−c(d)), the transition more than halves N. Deleting a singleton also gives N'=(N−hD)/2<N/2. Stop at the last literal factor or last singleton, before an empty reservoir. Thus this is an actual halving-speed target descent with legal physical steps, rather than a congruence itinerary that might lose its interface.

The carry search, given the divisor list, takes N^(log₂24/(8log₂3)+o(1)) arithmetic operations, with exponent about 0.362; the subpower factor accounts for divisor states and polynomial-in-log N work. Elementary trial factoring is a separate O(sqrt D) option. The replay checks the pruning against all cuts through chain length 14 and verifies actual absolute compiler outputs. This does not prove that a prescribed arbitrary target enters the grammar.
