# Complete entry classification for every one-reset controller

**PROVED, scoped.** For each L>=0 the controller w=1^L 2 has an exact all-layer maximal content mean over **every finite positive entry word**. We classify all entries attaining that mean, including every order-lifting stall, and prove a uniform strict clock bound below1/2. This does not prove the arbitrary-controller5/8 conjecture or Collatz convergence.

The exact marked periodic-mean theorem from [Session114](../../session114_continuation_symmetry_20261030/collatz/PERIODIC_CONTENT_RATES.md) is inherited. The new result removes entry optimization and every subgroup-membership test for this complete family. Its maximizing phases have an exact Lucas/no-carry classification.

## 1. Formula, bound and equality entries

For a positive valuation word x, let F_x(X)=(P_x X+B_x)/Q_x, and c(x)=v5(gcd(P_x−Q_x,B_x)). For fixed entry u and controller w, let R(u,w) be the ordinary mean per complete controller block, as defined and proved to exist in Session114. All layers are included. Fixed exits do not change this mean.

Set

 K=floor(L/2), delta=L mod2,
 r_L=v5(9^(L+1)+16·4^L),
 T_e=4·5^max(0,e−r_L),  h_e=2·5^(e−1).

Then r_L>=1 and

 **max_u R(u,1^L2)=R(empty,1^L2)
   = (r_L+1/4)/4 + sum_(e>=1) floor(L/h_e)/T_e.**       (1)

The sum is finite, since floor(L/h_e)=0 for h_e>L. The original clock per block is L+2, and uniformly

 **R(u,1^L2)/(L+2) <= (L/2+9/16)/(L+2) < 1/2.**       (2)

The first inequality is a convenient explicit bound, not asserted optimal. The finite tests find maximal clock value9/32 at L=0, but an all-L9/32 theorem is not claimed.

All equality entries in (1) are classified as follows. Choose a phase j in {0,...,L} satisfying

 **j mod(2·5^k) <= L mod(2·5^k) for every k>=0.**       (3)

Then choose an integer b>=0 such that the word equation

 **1^j u = (1^L2)^b**                                 (4)

holds. These, and only these, are the entries attaining (1). For b=0 this means j=0,u empty. For j>0,b>=1, it says u=1^(L−j)2(1^L2)^(b−1); for j=0 it says u=w^b.

Writing j=2a+epsilon, (3) is equivalent to

 epsilon<=delta, and every base5 digit of a is at most the corresponding digit of K.

Equivalently, epsilon<=delta and binom(K,a) is nonzero modulo5. Thus the number of maximizing phases is exactly

 **(delta+1) product_i (K_i+1),**                      (5)

where K_i are the base5 digits of K. Each phase gives an infinite family of entries through (4).

## 2. Why every entry has the same sharp layerwise cap

Let A=F_1 and M=F_w. The controller coefficients are

 P=3^(L+1), Q=2^(L+2), B=3^(L+1)−2^(L+1).

Modulo5 the multiplier P/Q is ±2, hence has order4. Only P²+Q² is divisible by5 among the factors of P⁴−Q⁴, so the exact affine order modulo5^e is T_e above.

The all-one map has coefficients (3^j,2^j,3^j−2^j) and order h_e modulo5^e. The commutator numerator between M and A^j is

 (P−Q)(3^j−2^j)−B(3^j−2^j)
       = −2^(L+1)(3^j−2^j).

Consequently, if A^j belongs to <M> modulo5^e, commutation forces3^j=2^j modulo5^e, and then A^j is identity. Thus

 **<A> intersect <M> is the trivial group modulo5^e.**  (6)

This holds at every depth, independently of order stalls.

The L+1 controller phases can be represented by A^j for j=0,...,L: phase0 replaces the completed-block phase M, which does not alter cyclic-subgroup membership. For a fixed entry, successful phases form the intersection of <A> with one left coset of <M>. By (6), this intersection contains at most one group element, hence at most one residue j modulo h_e. Among j=0,...,L, its multiplicity is at most1+floor(L/h_e).

Each successful phase contributes exactly1/T_e to the layer-e mean by the inherited marked-periodic theorem. Empty entry selects residue0 and attains the maximal multiplicity at **every layer simultaneously**. Therefore

 max_u R(u,w)=sum_(e>=1) [1+floor(L/h_e)]/T_e,

which is (1). No average over entries, phase independence, or finite truncation replaces the ordinary all-layer mean.

## 3. Complete maximizing-entry classification

Every layer cap has a nonnegative deficit. Equality in their sum requires equality at every layer.

If all phase depths of an entry are finite, then high layers contribute zero, while the empty-entry cap remains positive. Therefore an equality entry has an infinite-depth phase j. The inherited faithful positive-word representation and commuting-word theorem imply that 1^j u is a power of w: w is primitive, since it contains exactly one letter2. This proves the necessity of (4), including the phase0 interpretation.

For an entry satisfying (4), the successful phase positions at level e are exactly those congruent to j modulo h_e, by (6). Their number is

 1+floor(j/h_e)+floor((L−j)/h_e).

It equals1+floor(L/h_e) exactly when j mod h_e<=L mod h_e. This proves (3) and its sufficiency.

For h=2 the condition is epsilon<=delta. At h=2·5^k it reduces to a mod5^k<=K mod5^k. These inequalities for all k are equivalent to digitwise a_i<=K_i: if the first offending digit had a_i>K_i, its contribution5^i could not be cancelled by the total lower-digit deficit, whose magnitude is less than5^i. Conversely digitwise comparison gives every residue inequality. Finally, in characteristic5,

 (1+X)^K = product_i (1+X^(5^i))^(K_i),

so the coefficient of X^a is nonzero exactly for those digitwise choices. This proves the Lucas criterion and count (5) without importing a probabilistic interpretation.

The connection contract is exact: entry phases -> residue occupancies -> no-carry binomial support. It preserves the predicate of maximizing **every return layer**. The phase j and completion word (4) are retained; a bare count or binomial residue does not reconstruct the original entry.

## 4. A uniform bound surviving all stalls

The internal sum in (1) satisfies

 sum_(e>=1) floor(L/h_e)/T_e
 <= sum_(e>=1) (L/(2·5^(e−1)))/4 =5L/32.

Also

 5^r_L <=9^(L+1)+16·4^L <=25·9^L.

The elementary integer inequality9^8<5^11 gives r_L<=2+11L/8. Inserting this and the internal cap into (1) yields

 R <=r_L/4+1/16+5L/32 <=L/2+9/16.

This proves (2), with no assumption on the size or frequency of order stalls. The earlier hostile L=10 (which refutes a universal5/4 multiplier from first-layer mean to full mean) is covered.

## 5. Exact stall classification

Since9/4 lies in1+5 Z_5 with v5(9/4−1)=1, it generates that principal-unit group at every finite depth. Therefore for every e>=1 there is one residue ell_e modulo5^(e−1) such that

 r_L>=e iff L=ell_e modulo5^(e−1).

The residues are nested and are determined by9(9/4)^ell_e+16=0 modulo5^e. This classifies all stalls without a bounded-L scan; each next digit has a unique lift. It does not imply an untruncated natural-average formula in L without additional uniform tail control. Formula (1) and bound (2) need no such averaging.

## 6. Scope and verification

`check_one_reset.py` compares the formula against exact marked periodic means, explicitly enumerates the trivial subgroup intersections at finite depth, tests entry optimization and equality classification, and checks the Lucas equivalence. Normal and optimized gates remain active.

The all-one equality5/8 and the geometric tail1/8 are inherited. Here the genuine new synthesis is cyclic-coset occupancy with Lucas support, not an assertion that every appearance of1/8,3/8,5/8,7/8 measures the same object. General changing controllers and Collatz convergence remain open.
