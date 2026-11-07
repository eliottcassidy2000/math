# Session216: overnight tournament spectrum constructions

**Research ACTIVE.** The Chan–Pak target, every odd Hamiltonian-path count at least23, remains **OPEN**. This index records independently verified checkpoints while the remaining construction work continues. No unrestricted finite-range extension, positive-density theorem or eventual threshold is claimed.

## 1. Exact inverse for arbitrary apex-mask words

[Proof and implementation](inverse/README.md), [independent audit](audit/inverse/REVIEW.md), [release approval](audit/inverse/RELEASE_REVIEW.md).

For a fixed labeled k-apex tournament, arbitrary incidence words on a transitive chain admit an exact finite saturated-state inverse. It returns a literal witness or an all-length absence certificate in that class; a budget-limited run returns INCOMPLETE. Hidden saturated coordinates may later disappear and are retained correctly.

The two-apex implementation independently exhausts708 states for target7 and4626 for target21. These are class-relative controls, not new ambient exclusions. All4998 positive odd targets through10000 other than7 and21 have literal witnesses on at most twelve vertices, independently recounted. Bounded digit-grammar probes do not provide a closed induction. Transfer and saturation methods are inherited from Sessions215 and203; this is their nonuniform mask-language specialization.

Published [checkpoint1](https://github.com/eliottcassidy2000/math/commit/94e3ca9e045266b3f71351b799415f820a9bf268) contains this first packet.

## 2. Eleven adjacent odd counts at every dyadic scale

[Construction](interval/README.md), [independent audit](audit/arithmetic/REVIEW.md).

For each k>=1, explicit strong tournaments on at most k+7 vertices realize all eleven odd values from36*2^k-7 through36*2^k+13. More generally eleven fixed rooted strong contexts satisfy

    H(G_d[U])=36 Q_U(2)+d H(U),     d=-7,-5,...,13,

for every nonempty tournament U. Here Q_U(2) counts ordered two-bin spanning path covers, allowing empty bins. This is a uniform family of short intervals, not a continuous tail. Different values select different contexts; no single-edge staircase is asserted. All891 literal composites and every complete run-coefficient vector were independently checked.

## 3. Matching number classifies the available module runs

[Proof and physical compiler](star/README.md), [independent audit](audit/star/REVIEW.md).

For any marked tournament, the nonzero root-module run coefficients occur exactly at degrees0 through nu+1, where nu is the matching number of arcs from the root's outneighbors to its inneighbors. Matching controls support, not multiplicity. A matching supplies separated root excursions, and path insertion fills the remaining vertices without deleting their separators.

Exactly two possible nonempty runs correspond to a star-shaped cut. Removing its center gives an ordered L,root,R decomposition and a factored endpoint/path-cover formula. This explains the eleven contexts through several actual factorizations of36.

The same interface has a returning four-state append operation. Four explicit source states retain four adjacent odd choices after every common binary word and every transitive root-module size. The source partition is retained throughout; these are composable positive families, not a recurrence on arbitrary unmarked counts. Their width is preserved, not amplified. Independent checks include5405 matching cases and17055 literal composite counts.

## 4. A logarithmic improvement to the spectrum lower bound

[Audited derivation and reproduction](density/REPRODUCE.md).

If A(X) counts distinct tournament Hamiltonian-path counts at most X, then for every real X>=4096,

    A(X) >= X^(2/3)/(60 (ln X)^(2/3)).

The proof optimizes the splitting parameter in the inherited THM-4887 strip estimate. It improves the previous asymptotic logarithmic loss from power1 to power2/3; the power exponent2/3 is unchanged. The corresponding exact-energy bound is7.3*2^(4n/3)*n^(2/3) for n>=12. This does not prove linear energy, positive density, density one or completeness. An algebraic typo in the bound-only optimality discussion was repaired and preserved with its audit lineage; headline constants are unchanged.

## Continuing obligation

Further work tests source synchronization, richer two-apex constructions and interval-width amplification. The missing step is an all-target physical construction or proved occupancy theorem. Finite certificates, a fixed-width interval family and analytic support bounds are separate kinds of progress and do not supply that step by themselves.
