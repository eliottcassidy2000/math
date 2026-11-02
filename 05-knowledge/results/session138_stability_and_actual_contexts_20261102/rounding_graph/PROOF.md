# Independence-polynomial certificates for weighted tournament rounding

Status: independently reviewed and accepted within the stated domains. This extends the binary-core
argument in `../stability/MULTISTAR_STABILITY.md`; it makes no novelty claim.

Let A be a weighted tournament of order n>=3, H its weighted directed
Hamiltonian-path sum, u=det(I+A), v=(-1)^n det A, and Delta=H-u-|v|.
Write q_e=A_ij A_ji. Let G have precisely the nonbinary edges of A.
For a graph J on r>=1 vertices, write i_even(J),i_odd(J) for the numbers
of independent vertex subsets of even and odd size (including the empty
set in the even count), and define

    mu(J)=2^(1-r) min(i_even(J),i_odd(J)).

## Determinant lemma

For any weighted tournament B whose uncertainty graph is contained in J,

    det(I+B)-|det B| >= mu(J).

Indeed the augmented-skew identities express u+v as 2^(1-r) times the
sum of even principal Pfaffian squares on B, and u-v as the corresponding
sum of principal Pfaffian squares containing a binary observer. Every
independent subset of J has entirely binary skew entries. An even-size
binary tournament Pfaffian is an odd integer: its expansion has an odd
number (2l-1)!! of signed unit terms. Its square is at least one. For odd
independent subsets, adjoining the observer gives the same conclusion.
Retaining only these terms gives the two separate parity counts, and then
taking their minimum proves the assertion. The empty Pfaffian is one.

## An arbitrary edge-order certificate

Choose any deterministic order e_1,...,e_m of the edges of G. Set
G_t=G-{e_1,...,e_(t-1)} and J_t=G_t minus both endpoints of e_t. Then

    Delta(A) >= E Delta(T) + sum_t q_(e_t) mu(J_t),       (1)

where T is the product-law independent binary rounding of A.

For one edge e with probability p and pair-deletion remainder R,

    G_epsilon(A)=(1-p)G_epsilon(A^0)+p G_epsilon(A^1)
                 +p(1-p)(u_R-epsilon v_R),

with G_epsilon=H-u+epsilon v. Taking the minimum over epsilon gives
Delta(A)>=(1-p)Delta(A^0)+p Delta(A^1)+q_e(u_R-|v_R|).
At step t the uncertainty graph of R is J_t, independent of the previous
binary orientation outcomes. Its coefficient is at least mu(J_t) on
every branch. Iterating conditional expectations proves (1). The fixed
edge order may depend on the original probabilities, but not on future
rounding outcomes. n>=3 guarantees all J_t are nonempty vertex graphs.

Thus an edge order and exact independent-set counts give a finite,
orientation-independent stability certificate. The method does not imply
that the best such certificate is the optimal tournament inequality.

## Immediate elementary consequences

For n=3 each remainder has one vertex and mu=1. Hence Delta>=E Delta(T)+V.
Every binary triangle has Delta=0, and the direct determinant expansion
also gives Delta=V, recovering the exact identity.

For n=4, the remaining two vertices have mu=1 if their mutual edge is
binary and mu=1/2 otherwise. Partition the six possible edges into the
three complementary disjoint pairs. In each pair, round the smaller-q
edge before the larger-q edge (zero-q edges are already binary). The
contribution of a pair with variances a<=b is at least a/2+b, which is
at least 3(a+b)/4. Thus

    Delta(A) >= E Delta(T)+(3/4)V.

This proves the sharp unreserved order-four coefficient without a
multivariable polynomial certificate; two ordinal fair pairs attain it.
No new complete equality classification is asserted here.

For n=5, any pair-deleted remainder has three vertices. If it has s
uncertain edges, its parity counts give mu=1-s/4, for all s=0,1,2,3.
Order all uncertain edges by nondecreasing q. Each unordered disjoint
edge pair contributes the smaller variance exactly once to the loss.
Consequently

    Delta(A) >= E Delta(T)+V
                -(1/4) sum_{disjoint unordered edge pairs e,f} min(q_e,q_f).

This is a valid support-sensitive estimate, not the sharp universal
7/8 theorem. It can be weaker when uncertainty is dense.

## Matching uncertainty and the vertex-cover bound

For c disjoint uncertain pairs and t isolated vertices, the independence
polynomial is (1+2z)^c(1+z)^t. If t>=1 the even and odd subset counts are
both 3^c 2^(t-1), so mu=(3/4)^c. If t=0 and c>=1, the parity imbalance
has absolute value one, giving mu=(3/4)^c-4^(-c).

If J has a vertex cover of size c and a nonempty complement of size r-c,
its independent subsets include all subsets of that complement. Each
parity count is at least 2^(r-c-1), recovering mu>=2^(-c).
The case c=r uses the separate uniform determinant bound, not this
independent-subset count.

## Fixed-parameter evaluation

If a vertex cover C of J of size c is supplied, put I=V(J)\C. For every
independent S subset C, let d=|I\N(S)|. If d>0, add 2^(d-1) to each
parity count. If d=0, add one only to the parity of |S|. Every independent
set is counted once as S together with a subset of I\N(S). This computes
mu exactly in O(2^c poly(n)) operations (integer bit lengths are O(n)).
It gives a fixed-parameter certificate evaluator, not a polynomial-time
algorithm for general independent-set counting or Hamiltonian counting.

## Matching-support equality boundary

The matching-support constants are proved in the MULTISTAR appendix.
For n>=4, equality in the unreserved inequality Delta>=a_r V (with
isolates) or Delta>=b_r V (without isolates), and V>0, forces
E Delta(T)=0. The inherited binary equality theorem then forces every
supported completion to be transitive: order n>=4 excludes the exceptional
three-cycle. Thus A lies in the H=1 locus. Its scalar minimization implies
all r pairs are fair when r>=2. With r=1 the coefficient is one and any
fractional pair probability is allowed, provided the H=1 condition holds.
Order three must remain separate, since Delta=V holds identically even
when supported completions include a directed three-cycle. This argument
does not claim rigidity from equality in the stronger baseline-containing
inequality alone.
