# MRB parity, odd-cycle packings, and exact edge responses

## Scope

This is an exact transform of the tournament path-cover profile, inspired by alternating-series parity. It does **not** establish a complete spectrum of tournament path counts, nor imply that an analytic positivity statement gives a construction at every target count.

Let T be a finite tournament on n vertices. Let K range over packings of pairwise vertex-disjoint directed odd cycles, of length at least 3. Define

    d(K) = sum_C (|C|-1)/2,
    B_d(T) = sum_{K:d(K)=d} 2^{|K|},
    P_T(z) = sum_d B_d(T) z^d,
    F_T(r) = sum_K 2^{|K|} r^{n-2d(K)} = r^n P_T(r^-2).

The last display is a polynomial identity after clearing negative powers; the evaluation at r=0 is by the polynomial, not the quotient. B_0=1; for n>=1 the exponents n-2d(K) are positive. The path-cover theorem identifies F_T(1) with the number h(T) of directed Hamiltonian paths.

## MRB diagonal transform

Set a_m=m^{1/m}-1 and, for d>=0,

    lambda_d = sum_{m>=1} (-1)^m a_m / m^{2d}.

Each converges by the alternating-series test (after m=3, a_m strictly decreases to zero). The constant lambda_0 is the usual MRB series with this sign convention. For each fixed real z, the finite polynomial P permits termwise interchange, giving

    Lambda_T(z) := sum_{m>=1} (-1)^m a_m P_T(z/m^2)
                 = sum_d lambda_d B_d(T) z^d.                (1)

All lambda_d are strictly positive. For d>=1, group successive terms 2k and 2k+1: all terms after the first pair are positive because a_m decreases for m>=3; the first is positive since 9a_2>4a_3. For d=0, the sum of the first two pairs is 2a_2-a_3-a_5>0, and every later pair is positive. For example, bounds .4142<a_2<.4143, .4422<a_3<.4423, .3797<a_5<.3798 prove these inequalities by elementary powering. Thus (1) is an injective diagonal change of coordinates on the entire finite defect polynomial:

    B_d(T) = [z^d] Lambda_T(z) / lambda_d.

In particular Lambda_T(z)>=lambda_0 for z>=0, with equality at z>0 exactly when T has no directed odd cycle, equivalently is transitive. The forward implication uses the standard fact that every nontransitive tournament has a directed triangle. This produces a real analytic certificate of nontransitivity, but its coefficients are already combinatorial and the transform alone provides no new realizability theorem.

## Exact finite-difference cancellation

For K, put c(K)=n-2d(K). For k>=0, inclusion-exclusion and the Stirling numbers of the second kind give

    Delta^k F_T(0) = k! sum_K 2^{|K|} S(c(K),k).              (2)

Every summand is nonnegative; this is the number of surjective labelings of the c(K) abstract components by k labeled colors, summed with packing weight. The components are n-|union K| unmatched singleton vertices plus one abstract component for each selected cycle; no canonical graph structure on them is needed. Since every exponent has the parity of n,

    F_T(-r)=(-1)^n F_T(r).

Newton interpolation at -1 therefore yields the finite, exactly canceling identity

    sum_{k=0}^n (-1)^k Delta^k F_T(0)=(-1)^n h(T).           (3)

This parity cancellation is exact: alternating signs in (3) cancel nonnegative surjection counts, just as paired MRB terms cancel, but (3) is not obtained by substituting MRB values into a graph theorem. The finite differences determine F_T via F_T(r)=sum_k binom(r,k) Delta^k F_T(0), hence recover its exact defect spectrum. They do not assert which spectra are attainable.

## Signed odd-cycle edge response

Let T' be obtained from T by reversing e={u,v}. Consider oriented simple odd cyclic orders C using e and with every other arc directed correctly in T. If the e arc is correct in T', set epsilon(C)=+1; if it is correct in T, epsilon(C)=-1. Then

    F_T'(r)-F_T(r) = 2r sum_C epsilon(C) F_{T-V(C)}(r).     (4)

Proof: a packing unaffected by e cancels. In a packing that changes, its unique selected cycle containing e is precisely one of those C, and the other cycles form an arbitrary packing of T-V(C). A selected cycle contributes factor 2r while each untouched component contributes its usual factor. Since both endpoints of e are in C, deleting V(C) makes the residual tournament the same before and after reversal.

At r=1, each residual weight is the positive integer h(T-V(C)). Consequently, raw counts of created versus destroyed cycles do not determine the sign. Let c(T) be the number of directed odd cycles. Reducing the packing sum mod 4 gives

    h(T')-h(T) = 2(c(T')-c(T)) (mod 4).

More generally, if N_j(T) denotes the number of packings of exactly j cycles, then h(T)=sum_j 2^j N_j(T), and h(T) mod 2^k depends only on j<k.

### Explicit sign reversal

On vertices 0,...,6 take arcs

    0 -> 1,2,3,4;  5,6 -> 0;
    1 -> 2,3,5;    4,6 -> 1;
    2 -> 3,4,5;    6 -> 2;
    4 -> 3;        3 -> 5,6;
    4 -> 5,6;      6 -> 5.

Reversing 0->3 creates 8 directed odd cycles and destroys 6, but h changes 115->111. The corresponding sums of h(T-V(C)) are 10 for creation and 12 for destruction, hence (4) gives 2(10-12)=-4. The accompanying independently executable search.py enumerates oriented cyclic orders and disjoint triangle pairs; for seven vertices, a pair of disjoint directed odd cycles must be two triangles, and three cycles cannot be disjoint.

## An explicit parity-separated construction

Let T_n be transitive on 1,...,n except that 1->n is reversed, for n>=3. Every directed cycle must contain n->1, because all other arcs strictly increase the label. Its remaining arcs therefore traverse an increasing subsequence of the n-2 internal vertices. Odd cycles correspond exactly to odd subsequence sizes k. All these cycles meet at vertices 1 and n, so no packing contains two of them. Consequently

    F_{T_n}(r)=r^n + 2 sum_{k odd, 1<=k<=n-2} binom(n-2,k) r^{n-k-1},
    h(T_n)=1+2^{n-2}.

The second identity is the even/odd binomial decomposition, with the empty packing contributing 1. If tournaments T_1,...,T_s are arranged in an acyclic composition (all arcs go from each earlier block to each later block), no directed cycle can meet more than one block. Thus F=product_i F_{T_i}, and h=product_i h(T_i). In particular every finite product of 1+2^k for k>=1 is realized. This is an explicit subfamily of the spectrum, not a characterization.
