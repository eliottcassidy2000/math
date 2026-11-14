# The full fixed-count inverse fiber and its neutral dimensions

The five-coordinate construction in [the Neményi comparison](BUSCH_NEMENYI_FERMAT.md) can be embedded in a complete rational generating function at each fixed prescribed h. This concerns the number of witnesses as order grows, not the existence of any witness.

Fix h>1. Let an ordered skeleton be a list S=(C_1,…,C_k) of isomorphism classes of strong tournaments of order at least three with product H(C_i)=h. Let s(S) be its total order. There are finitely many skeletons: every H(C_i)≥3, and [Busch](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) bounds |C_i| by 1+3 log_5 h. For h=1 use only the empty skeleton, allowing the empty tournament at order zero. Inserting g_i singleton components in its k+1 gaps is injective and surjective onto the tournaments with that skeleton. The ordinary class generating function is therefore

    A_h(z) = Σ_S z^s(S)/(1−z)^(k(S)+1).

In particular the class count is eventually a polynomial of degree exactly d=max_S k(S) whenever h is attained. Its leading coefficient is (# ordered skeletons of length d)/d!. This is the exact maximum number of independent singleton-gap coordinates minus one; the total-order constraint consumes one coordinate. If h is absent the series vanishes identically.

The maximal skeleton length equals the maximal number of factors >1 in a factorization of h into attainable counts. One direction reads off the SCCs. For the other, split each witness into its own nontrivial SCCs, then concatenate. Thus the top degree is a multiplicative invariant of the entire spectrum. In particular if every prime factor of h is attained, d=Ω(h), where Ω counts multiplicity. Otherwise at least one strong core necessarily carries more than one prime factor, and d<Ω(h).

## Exact converse folding

Converse sends the skeleton to (C_k*,…,C_1*) and reverses the gap vector. Self-converse witnesses come only from skeletons fixed by this involution. Each pair of mirrored gaps contributes a factor 1/(1−z²), and an unpaired central gap exists precisely for even k. Consequently the self-converse and converse-merged class series are

    B_h(z)=Σ_{S=S*} z^s(S)/[(1−z²)^floor((k+1)/2) (1−z)^[k even]],
    M_h(z)=(A_h(z)+B_h(z))/2.

The first formula proves the entire remaining periodicity has period dividing two after the transient, with no larger quasipolynomial periods coming from the padding. For squarefree h>1 any mirrored pair of nontrivial cores would contribute a square divisor; hence only k=1 can occur in B_h. Its unique gap parameter increments total order by two. For the product h=5·17·257·65537 the leading coefficient of A_h is the product of the numbers of strong realizers of those four counts. In fact the top-dimensional part consists precisely of their four independent core choices in any order, with five gaps. No such top-dimensional skeleton is self-converse.

This general formula tells which dimension belongs to the observation and which dimension belongs to the physical interface. For example, all prime-count fibers have eventual degree one, but a lower bound on the coefficient or a characterization of the available endpoint data requires a theorem about the strong cores themselves.
