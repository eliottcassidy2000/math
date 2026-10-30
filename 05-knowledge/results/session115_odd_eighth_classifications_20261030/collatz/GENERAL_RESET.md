# Every one-reset controller: exact maxima and complete mixed-radix entry classes

**PROVED, scoped.** For every a>1 and L>=0, the controller w=1^L a has an explicit maximal ordinary all-layer content mean over **all finite positive entries**, and all entries attaining that mean are classified. This includes all reset residue classes and all order stalls. No universal5/8 clock bound is claimed here; the separate a=2 and constant-letter bounds retain their stated scopes.

The proof extends [the a=2 classification](ONE_RESET_CLASSIFICATION.md). The new coordinate is the full cyclic-intersection modulus g_e. When the reset factor is a5-adic unit, g_e is the ordinary all-one period. Otherwise it changes the maximizing-phase radix, and dropping it gives false counts.

## 1. All parameters and the exact answer

Let

 P=3^(L+1), Q=2^(L+a), B=3^(L+1)−2^(L+1),
 h=v5(2^a−2), t=ord_5(P/Q),
 r=v5(P^t−Q^t) if t>1,
 r=min(v5(P−Q),v5(B)) if t=1,
 T_e=t·5^max(0,e−r).

All valuations in these parameter definitions are finite. Define the phase modulus g_e as follows:

 **t>1, e<=h:**  g_e=5^max(0,min(r,e)−1).

 **t>1, e>h:**   g_e=2·5^max(0,e−h−1,min(r,e)−1).

 **t=1:**        g_e=2·5^max(0,min(r,e)−1,e+r−h−1).

Then

 **max_u R(u,w)=R(empty,w)
   =(r+1/4)/t + sum_(e>=1) floor(L/g_e)/T_e.**         (1)

The correction sum is finite: g_e is nondecreasing under divisibility and eventually grows fivefold at each depth. Divide by L+a for the rate per original valuation clock. The exact ordinary-mean theorem and fixed-exit invariance are inherited from [Session114](../../session114_continuation_symmetry_20261030/collatz/PERIODIC_CONTENT_RATES.md); this result performs the previously unbounded entry optimization exactly for the whole family.

All maximizing entries, and only those entries, satisfy

 **1^j u=w^b,  0<=j<=L, b>=0,**                      (2)

for a phase j with

 **j mod g_e <= L mod g_e for every e>=1.**           (3)

Only depths with g_e<=L need checking. The case b=0 is precisely j=0,u empty. For positive b and j>0, u=1^(L−j)a w^(b−1); for j=0, u=w^b.

## 2. Complete mixed-radix description and phase count

Delete repeated values from the divisibility chain1,g_1,g_2,... . The successive ratios are the following explicit radix sequences:

- If h=0, or if t=1: **2,5,5,5,...**.
- If h>0,t>1,r<=h: **r−1 copies of5, then2, then5,5,...**.
- If h>0,t>1,r>h: **h−1 copies of5, then10, then5,5,...**.

Zero copies means the indicated initial block is absent. Let L_i and j_i be the ordinary positional digits in this same mixed radix. Then(3) is equivalent to **j_i<=L_i for every i**. Thus the exact number of maximizing phases is

 **product_i (L_i+1).**                              (4)

This is a finite product, since eventually all digits are zero. In the first case it becomes the earlier Lucas formula: L=2K+delta, j=2c+epsilon, epsilon<=delta and binom(K,c) nonzero modulo5; the count is(delta+1) product_i(K_i+1). The other two cases are exact mixed-radix no-borrow conditions, not automatically binomial-support identities in one prime characteristic.

## 3. Deriving the intersection modulus

Write A=F_1 and M=F_w. Modulo5^e, the order of A is h_e=2·5^(e−1), and that of M is T_e. Since<A> is cyclic, its intersection with<M> is exactly

 **{A^j : j=0 modulo g_e}**                          (5)

for a divisor g_e of h_e. We now derive the displayed formula rather than assuming the intersection is trivial.

For A^j let z_j=v5(3^j−2^j). For j>0,

 z_j=0 if j is odd, and z_j=1+v5(j) if j is even.

The affine commutator numerator is

 (P−Q)(3^j−2^j)−B(3^j−2^j)
   =−(2^a−2)2^L(3^j−2^j),

with valuation h+z_j.

If t>1, then t is2 or4 and the residue subgroup contains both1 and−1. The principal-unit test for A^j is

 1+v5(j)>=min(r,e),

because v5((3/2)^(tj)−1)=1+v5(j). The commutator test is h+z_j>=e. When e<=h the latter is automatic, giving the first g_e formula. When e>h it requires j even and1+v5(j)>=e−h, giving the second formula. The identity exponent j=0 satisfies every test.

If t=1, put s=v5(P−Q), b0=v5(B), so r=min(s,b0). Both the multiplier and translation depths of A^j equal z_j. The two pivot cases in the inherited exact cyclic-membership theorem consequently reduce to

 z_j>=min(r,e), and h+z_j>=e+r.

Their combined threshold J=max(min(r,e),e+r−h) is positive for these actual controllers. Indeed, if h=0, then r=0 and J=e. If h>0, necessarily a=1 modulo4; t=1 then forces L odd, giving s,b0>=1 and r>=1. Thus j must be even and1+v5(j)>=J, proving the third formula.

For completeness, r<=h whenever t=1: the difference(P−Q)−B has valuation h, so min(v5(P−Q),v5(B)) cannot exceed h. Reading the formulas now shows that g_e has exactly the distinct-value chains and radix ratios in section2. This also proves g_e divides g_(e+1), and eventual fivefold growth.

## 4. Why empty entry is globally optimal, and why the equality list is complete

The L+1 phases can be represented by A^j, j=0,...,L. Position0 replaces the completed-block phase M, since multiplying by M does not affect membership in<M>. For a fixed entry, successful A^j lie in the intersection of<A> with one left coset of<M>. If nonempty this is one coset of the intersection group(5), so successful positions form exactly one residue class modulo g_e.

The maximal number of such positions is1+floor(L/g_e), attained by residue0. Empty entry attains this at every layer. Each successful position has density1/T_e in controller time, so summing gives(1).

Equality requires equality at every layer because all layer deficits are nonnegative. An entry whose phase depths are all finite misses every sufficiently high layer and is strictly below(1). Thus an equality entry has an infinite-depth phase. Faithfulness of positive valuation words and the commuting-word lemma imply(2), because w has exactly one non-one letter and is primitive.

For an entry satisfying(2), its successful positions are j modulo g_e. The count is

 1+floor(j/g_e)+floor((L−j)/g_e).

It equals1+floor(L/g_e) exactly when(3) holds. This proves necessity and sufficiency. Finally, for any divisibility chain of place values, all prefix-residue inequalities are equivalent to digitwise comparison: an offending digit contributes at least its place value, greater than the total possible lower-digit deficit. Independent digit choices give(4).

## 5. The former obstruction becomes a retained coordinate

At a=5,L=2, w=(1,1,5), both<A> and<M> modulo5 equal

 {(1,0),(4,3)}.

All three positions0,1,2 are successful at the first layer. The unit-reset cap1+floor(2/2)=2 is false, but the general formula has g_1=1 and counts3 correctly. Its actual maximal block mean is69/40.

A visible radix switch occurs at a=5,L=8: h=1,t=2,r=3, so the radix sequence is10,5,5,... . Every j=0,...,8 is maximizing; the unit-reset Lucas condition would retain only five of those phases. The exact maximal block mean is45/8, with clock rate45/104. This is a mean per controller block, not a return probability.

## Verification and scope

`check_general_reset.py` tests the g_e formulas against literal cyclic-subgroup generation, the closed maximum against the inherited exact marked means, every small entry's bound and equality class, and the mixed-radix phase count. It includes high-h reset letters and both nontrivial-intersection examples. All gates remain active under Python optimization.

This is a complete maximal-entry classification for one-reset controllers, not a classification of every arbitrary controller or of all values attained by nonmaximizing entries. The general exact marked-rate theorem still evaluates those individual rates. No universal all-layer finite-word bound, general5/8 asymptotic bound, or Collatz convergence claim follows.
