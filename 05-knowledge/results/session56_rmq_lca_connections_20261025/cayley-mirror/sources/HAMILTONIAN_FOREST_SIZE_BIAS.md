# Exact Hamiltonian-size-biased forest law for random tournaments

Status: PROVED, with independent proof and finite-exact audits in mirror/HAMILTONIAN_FOREST_AUDIT.md. Asymptotics, merged-node masses and the full-polynomial hostile are independently audited in bridge-audit/HAMILTONIAN_ADDENDA_AUDIT.md. No novelty assertion against the literature. This is an ensemble identity, not reconstruction of a tournament or its metagraph from a forest statistic.

## Objects and target

Let T be uniformly sampled from the 2^{n choose 2} labelled tournaments on [n], n>=2. Write L_T=diag(outdegrees)-A_T and F_T(t)=det(I+tL_T). The directed matrix-forest theorem says its coefficient of t^j counts spanning directed forests with j arcs, every nonroot pointing toward its component root. Let H(T) count directed Hamiltonian paths, and κ(T)=2^{n-1}H(T)/n!. Then Eκ=1.

Set a=1+nt/2. We claim the exact polynomial identities

E F_T(t)=a^{n-1},
E[κ(T)F_T(t)]=product_{k=1}^{n-1}(a+i t cos(kπ/n)).

For t>=0 this gives the real positive form

E[κF]/EF=product_{k=1}^{floor((n-1)/2)}(1+(t/a)^2 cos²(kπ/n)).

It is strictly larger than1 if n>=3,t>0. Hence H and this rooted-forest partition function are positively correlated under the uniform tournament measure. This does not assert monotone ordering between individual tournaments.

## Proof by conditioning, not a false determinant expectation rule

Expand F by the directed matrix-forest theorem. Each forest uses each unordered edge at most once. Its arc indicators therefore are independent Bernoulli(1/2). Taking expectations replaces every directed edge weight by1/2 in THIS forest polynomial. The resulting Laplacian is (nI-J)/2, so EF=a^{n-1}. In general expectation does not commute with determinants; the forest expansion is the justification here.

Write H as the sum of the n! indicators of ordered directed Hamiltonian paths. Relabelling symmetry gives E[κF]=E[F | 1→2→...→n]. Conditioning fixes these n−1 edges and leaves every other orientation independent and fair. Forest expansion again lets us replace each edge by its conditional probability. The conditional mean weighted Laplacian is

(nI-J)/2+B,

where B has diagonal (1/2,0,...,0,−1/2), upper superdiagonal −1/2 and lower subdiagonal +1/2. B1=0. On the quotient by the constant vector, J vanishes and the mean Laplacian acts as nI/2+Bbar. Therefore

E[κF]=det(aI+tB)/a,

initially a≠0, extending polynomially to all t.

## Tridiagonal spectrum

Put Q=2B. If f_j(x) obeys f_0=1,f_1=x,f_j=x f_{j-1}+f_{j-2}, expansion of the tridiagonal determinant gives det(xI-Q)=x f_{n-1}(x). Indeed the leading interior minors are d_j=f_j−f_{j-1}; the final endpoint yields (x+1)d_{n-1}+d_{n-2}=x f_{n-1}. Since f_j(x)=i^j U_j(x/(2i)), its zeros are 2i cos(kπ/n),1<=k<=n−1. Thus B has the asserted spectrum, including the additional zero eigenvalue. No diagonalizability is required; for even n zero can have a nontrivial Jordan block.

Substitution yields the product, and pairing cosines of opposite signs proves positivity.

## The user's coefficient sequence is exactly present

Let a_m=(m+1)^{m-1}/m!, including a_0=1. Then a_1,a_2,... are 1,3/2,8/3,125/24,54/5,16807/720,... . At t=2,

E F_T(2)=(n+1)^{n-1}=n!a_n.

Thus the Cayley/forest coefficients do occur in a natural tournament ensemble, and the Hamiltonian-biased correction is explicitly known. The expectation forgets individual orientations, roots, path incidence and metagraph wiring. It cannot imply a Hamiltonian-count spectrum theorem.

## Asymptotic consequence

For fixed real t>0,

E[κF_T(t)]/EF_T(t)=1+1/n+O_t(n^{-2}).

Proof: sum_{k=1}^{n-1}cos²(kπ/n)=(n−2)/2; summing over one member of each nonzero pair gives (n−2)/4. Also (t/a)^2=4/n²+O_t(n^{-3}). Expand log of the finite product. The fourth-power contribution is O(n·n^{-4})=O(n^{-3}); hence log ratio=1/n+O_t(n^{-2}), and exponentiate. The t=0 boundary gives ratio1 exactly and is not uniform with fixed-t asymptotics as t approaches0.

## Hostile and positive controls

For n=3, a directed cycle has H=3,F(t)=1+3t+3t²; a transitive tournament has H=1,F(t)=1+3t+2t². At t=1, the cycle forest count is7 although H=3. This is an explicit failure of preserving the Hamiltonian gap7 under the forest map.

Required tests: exact enumeration of all labelled tournaments through n=5; direct forest-polynomial determinants, independent subset-DP Hamiltonian counts; exact n<=12 tridiagonal polynomial identities. The companion checker records the tested scope.

## Explicit recurrence and coefficientwise consequences

The product has an integer-binomial form:

E[κF_T(t)] = sum_{l=0}^{floor((n−1)/2)} binom(n−1−l,l) a^{n−1−2l}(t²/4)^l.

This is the continuant recurrence P_0=1,P_1=a,P_m=aP_{m−1}+(t²/4)P_{m−2}, evaluated at m=n−1 (a depends on the fixed n). The l=0 term is EF. Every remaining coefficient is nonnegative, and the l=1 term is (n−2)t² a^{n−3}/4. Therefore Cov(κ,[t^j]F_T)>0 for each 2<=j<=n−1; the constant and one-edge coefficients have zero covariance.

Let W(T)=sum_v(outdeg(v)−(n−1)/2)^2. Since tr L=n(n−1)/2 and tr L²=sum_v outdeg(v)^2 (no two opposite tournament arcs),

[t²]F_T = ((n(n−1)/2)^2−sum_v outdeg(v)^2)/2.

Thus the exact mixed-moment identity is

Cov(κ,W)=−(n−2)/2.

Equivalently, for signed score vector d=(A−A^T)1, Cov(κ,||d||²)=−2(n−2). This matches the score-energy mechanism of the prior Walsh-chaos work at the level of an independently derived exact covariance; it does not imply score energy determines H.

For fixed t>0, expansion of the product further gives

E[κF]/EF = 1+1/n−(3/2+4/t)/n²+(4/t+12/t²−10/3)/n³+O_t(n^{-4}).

Derivation: with u=(t/a)^2, use log ratio=u(n−2)/4−u²(3n−8)/32+O(nu³), using exact cosine second/fourth-power sums for n>=3, then expand u=4n^{-2}(1−4/(nt)+12/(n²t²)+O_t(n^{-3})). The implied error holds for fixed t, not uniformly at t=0.

## Primary sources for imported forest facts

Chebotarev–Agaev, *Forest matrices around the Laplacian matrix*, https://arxiv.org/abs/math/0508178 . Chebotarev–Shamis, *Matrix-Forest Theorems*, https://arxiv.org/abs/math/0602575 . The ensemble conditioning and tridiagonal calculation above are proved directly here. Literature-wide priority for their combination has not been established.

## Exact mass constraints for the merged tournament metagraph

The ensemble identity descends to a concrete constraint on the user's metagraph. For an unlabelled tournament class [T], let

h(T)=H(T)/|Aut(T)|.

This equals the number of fixed-Hamiltonian-path presentations (tilings) of the class: Aut(T) acts freely on the ordered Hamiltonian paths, and each orbit has a unique relabelling with path1→...→n. Consequently sum_[T] h(T)=2^d where d=(n−1)(n−2)/2, the number of free nonpath arcs.

To merge converse pairs, give a self-converse node weight w=h, and a distinct converse pair weight w=2h. Define its forest observable as Fbar=F for a self-converse node and Fbar=(F_T+F_T*)/2 for a distinct pair. This averaging is essential: out-rooted forest counts need not be assumed converse invariant.

The exact mass rule is

sum_merged_nodes w(node) Fbar_node(t)
 =2^d sum_{l=0}^{floor((n−1)/2)} binom(n−1−l,l)(1+nt/2)^(n−1−2l)(t²/4)^l.

Proof: the left is sum_classes H(T)F_T/|Aut(T)|=sum_labelled H(T)F_T/n!. Substitute the ensemble identity. Merging is just a regrouping, not a new probabilistic assumption.

In particular, with signed score energy W_s=||(A−A^T)1||², which IS converse invariant,

sum_nodes w(node) W_s(node) = 2^d [n(n−1)−2(n−2)] =2^d(n²−3n+4).

These are exact polynomial and moment constraints on a candidate metagraph inventory. They can reject a proposed inventory even when its unweighted node/edge totals agree. They are not sufficient realizability conditions: neither blue/black incidence nor which pair of tilings is joined is retained. A valid realization still requires those incidence coordinates. The statistics are not the signed linear-forest loop coefficients of Session30; the directed rooted-forest Laplacian here is a different object with branching allowed.

### Full-polynomial hostile: counts do not recover Hamiltonian multiplicity

At order5, use unordered pairs in lexicographic order, with bit0 meaning the smaller label points to the larger and bit1 reversing it. Masks40 and41 both have

Fout(t)=Fin(t)=1+10t+39t²+70t³+49t⁴,

but their Hamiltonian counts are respectively15 and13, and automorphism orders3 and1. Their fixed-path presentation multiplicities are therefore5 and13. The arc lists and all Hamiltonian paths are in `forest_collision.json`; a separate audit reconstructs the result without using these supplied path lists. Both whole directed forest polynomials, not merely their quadratic coefficients, therefore fail to reconstruct H or the tiling multiplicity. This is compatible with all ensemble identities above.
