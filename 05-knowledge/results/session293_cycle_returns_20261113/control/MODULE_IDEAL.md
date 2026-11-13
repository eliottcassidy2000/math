# The exact arithmetic ideal of a uniform tournament module

Status: PROVED in note; [independently reviewed](../audit/REVIEW_IDEAL.md). No claim that every positive member of the ideal is attainable.

## 1. Inherited cover kernel, and an extra excursion factorial

For nonempty C, let F_k(C)=k!pc_k(C) be the ordered spanning directed-path-cover counts. The general substitution kernel is THM-3121, and the full cover fingerprint is THM-1975/3134/3166. A **uniform context** means an arbitrary tournament exterior E with each exterior vertex either beating every vertex of C or losing to every vertex of C; all other exterior arcs are unrestricted. C is therefore a module. The empty exterior is allowed.

Every Hamilton path splits uniquely into k nonempty C-runs and exterior excursions. Thus

H(E[C])=sum_k a_k(E)F_k(C),   a_k(E) nonnegative integers.

But the coefficients have a stronger divisibility:

**Excursion factor:** (k−1)! divides a_k(E).

Indeed an exterior skeleton with k C slots has a fixed prefix, k−1 nonempty internal excursions, and a fixed suffix. Each internal excursion starts at an exterior vertex beaten by C and ends at an exterior vertex beating C. Its vertices are disjoint from the other excursions. Permuting these excursions in any order, leaving prefix and suffix fixed, gives another valid skeleton; no arc between different excursions is traversed. The symmetric group on the k−1 excursions acts freely because their actual vertex sets are nonempty and disjoint. This proves the factor. It also explains why a naive triangle probe with two transitive blocks does not have top coefficient1: unused within-block edges do not force the order of singleton excursions.

Set D_k(C)=(k−1)!F_k(C). Every uniform-context Hamilton count is a nonnegative integer combination of the D_k.

## 2. Exact probe basis, ideal and fingerprint

For m≥0, make an exterior total order

x1→y1→x2→y2→...→xm→ym,

and orient C→xi and yi→C for every i. Let P_m(C) be its full Hamilton count (P_0=H(C)). Then

P_m(C)=sum_{k=1}^{m+1} b_mk D_k(C),

with nonnegative integers b_mk and b_m,m+1=1.

Proof. Each internal excursion has at least two exterior vertices, so k≤m+1. Equality consumes all2m exterior vertices in m two-vertex excursions and leaves prefix/suffix empty. The only allowed pairs are xi→yj with i≤j. A perfect matching with i≤σ(i) for every i is identity, by summing the inequalities. Thus the excursions are exactly xi,yi, each used once, in any of m! orders. The coefficient of F_(m+1) is m!, hence that of D_(m+1) is1. Integrality of every other coefficient follows from §1. □

Consequently, for n=|C|,

**Exact uniform-context ideal:**

gcd_{all uniform E} H(E[C])
 = gcd(P_0(C),...,P_(n−1)(C))
 = gcd_{k=1}^n ((k−1)!k!pc_k(C)).

The last equality follows by integer unitriangular inversion; every D_k is an integer combination of the first k probes. The first equality follows because probes are actual uniform contexts and every context count is a combination of the D_k. This is an equality of generated ideals, NOT a positive semigroup coverage theorem.

The same triangular inverse proves **complete contextual observability**: for two modules of order at most n, equality of their first n probes is equivalent to equality of the entire ordered-cover profile and hence equality in every uniform context. The division by factorials is exact over the integers. This is a specific physical probe realization of the inherited complete-profile invariant, rather than a claim that scalar H suffices.

## 3. Symmetry is sufficient, not necessary, for the arithmetic obstruction

Every automorphism of C extends to every uniform context by fixing its exterior. It acts freely on Hamilton paths, so |Aut(C)| divides the ideal generator. The generator can be strictly larger, and can be nontrivial even when C is asymmetric.

Examples use n5 forward edge bits in lexicographic pair order(01,02,03,04,12,13,14,23,24,34); bit1 means i→j.

* mask8: H=9, |Aut|=1, ordered covers (9,66,186,240,120), ideal generator3.
* mask76: H=15, |Aut|=5, ordered covers (15,90,210,240,120), ideal generator15.

The excursion factor can strictly strengthen the ideal even on an asymmetric
core. For order eight, the outgoing bit rows
`(56,117,129,22,100,76,141,59)` give `H=513`, trivial automorphism group,
and ordered covers `(513,6102,32508,99432,187200,216000,141120,40320)`.
Their unscaled gcd is 9; after excursion scaling the exact context divisor
is 27. No minimal-order claim is made.

Hence every uniform context containing the asymmetric mask8 module has count divisible by3. Every uniform context containing the second module has count divisible by15. Removing visible automorphisms alone does not remove arithmetic constraints on future contexts.

If two uniform contexts use the same C, their count difference is divisible by this generator. Thus a uniform-module repair of size2 requires generator1. This is necessary, not sufficient; even a primitive ideal gives signed integer combinations, not a physical target selector.

## 4. Exact finite prime-power precision

To determine whether p^s divides the generator it suffices to check

v_p((k−1)!k!)+v_p(pc_k(C))≥s

only for k for which v_p((k−1)!k!)<s. All higher k vanish automatically. The cutoff is the first k with the displayed factorial valuation at least s. This is an excursion-specific strengthening of Session292's factorial precision principle, not an identification of their different coefficients. The k=1 condition keeps H(C) itself. Since H(C) is odd, the ideal generator is odd.

At prime precision p, the relevant k are1,...,p−1; their factorials are units, so p divides the ideal iff pc_k(C)=0 mod p for every k<p (with pc_k=0 beyond the order). At higher prime powers the extra(k−1)! can matter.

## 5. Consequence for return design

This gives a computable all-context obstruction before experimenting with a seam. A source with nontrivial module ideal cannot support a uniform-context ±2 repair, regardless of how creative or large the exterior is. A successful construction must choose a primitive source, break the uniform module by selective ports, change cores, or compensate simultaneously across modules. The sparse-guard grammar from Session292 gives a lawful way to consider the last two possibilities, but no unrestricted selector or Chan–Pak proof follows.
