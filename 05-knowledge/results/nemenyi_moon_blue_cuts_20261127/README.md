# Neményi stress samples, a sharper Moon bound, and the blue cuts

2026-11-27. This extends the [five-stress atlas](../five_stress_connection_atlas_20261112/), in particular its determinant comparator. It does not identify stress trajectory curves with directed cycles, or imply Chan–Pak.

## 1. The geometry strengthens Moon inside the stress family

The classical isothermal result is attributed to [Neményi (1933)](https://doi.org/10.1002/zamm.19330130506); [Sáenz](https://link.springer.com/article/10.1007/s00707-008-0118-3) discusses the broader inverse problem, with its own assumptions, and states a two-to-five parameter result. In the local, body-force-free isothermal setting used in the [prior derivation](../five_stress_connection_atlas_20261112/STRESS.md), the deviator is q=(Z')²λ, where λ=A|z|²+2Re(Bz)+C is real. Its four coefficients and the invisible hydrostatic constant give five parameters. This counts stress fields, not tournaments. The sample tournament compares i→j iff Im(conj(q_i)q_j)>0, excluding zero and parallel deviators. On a connected region where λ never vanishes, the sign is constant, so the five parameters change none of these comparisons at fixed sites.

Sort the sampled nonzero rays of q by their angles. A subset is transitive exactly when its rays fit in an open semicircle. In one direction the angular order is transitive; in the other direction a transitive subset has a source and all other rays are in its forward open semicircle (one can include the source by a small shift). A subset with no such semicircle has all circular gaps strictly below π and its angular order traces a directed spanning cycle. This is the round-cycle case of [Brouwer's locally transitive model](https://ir.cwi.nl/pub/6805/6805D.pdf). For any strong sample tournament, the angular order is **intrinsic**: the immediate successor of a vertex is the first vertex in the transitive order on its out-neighborhood. Thus its angular spanning cycle is canonical, even if many other directed spanning cycles exist.

[Moon's 1966 theorem](https://doi.org/10.4153/CMB-1966-038-7) gives a directed k-cycle through every vertex of a strong tournament, for every k from 3 to n. Here one can start with **any specified directed triangle** and insert the remaining vertices **in any order**, obtaining a nested flag of actual cycles that never changes the cyclic order of the retained vertices. Each insertion just splits an angular gap of length less than π into two smaller gaps. For a chosen vertex, choose any triangle through it and apply this. The triangle property is essential for this argument; it does not say that an arbitrary prescribed arc of a general tournament extends.

Let d_i be vertex i's outdegree. The exact number S_k of strong induced k-subtournaments, equivalently of canonical angular k-cycles, is

    S_k = binom(n,k) − Σ_i binom(d_i,k−1),  3≤k≤n.

Indeed every non-strong subset is transitive and has a unique source i, and any k−1 vertices in the out-neighborhood of i give such a subset. The analogous rooted count at v subtracts binom(d_v,k−1)+Σ_(i→v)binom(d_i−1,k−2) from binom(n−1,k−1). These are not formulas for the total number of directed cycles. In particular a spanning subset contributes exactly one *canonical* cycle to S_n.

**Sharp strengthened Moon bound in this family.** For every strong sample tournament, S_k≥binom(n−2,k−2), simultaneously for all k. Indeed in any strong tournament, the sum of the j smallest degrees is at least binom(j,2)+1 for 1≤j<n: there must be an arc out of the corresponding proper nonempty vertex set. Reversing only the extreme arc of a total order attains all these lower prefixes (its degrees are the integers 0 through n−1 with 0 raised to 1 and n−1 lowered to n−2). Thus its degree vector majorizes all the others. The integer function f(d)=binom(d,k−1) is discretely convex by Pascal's identity. Its sum is maximized at the majorant: write f as an affine function plus a nonnegative combination of the hinges (d−t)_+, whose sum is the maximum over j of the sum of the j largest entries minus jt. Substitution gives the claimed lower bound, with equality at this tournament. In the extremizer, equivalently a cyclic triple with one transitive block of n−2 vertices and two singleton blocks, a subset is strong exactly when it contains both singleton blocks and k−2 vertices in the long block. Moon's general minimum n−k+1 counts a larger family than the locally transitive one considered here.

## 2. Exact blue canonical cuts

Mark the canonical cycle in its direction with indices i modulo n and retain the cyclic degree word d_i, not just the degree multiset. Vertex i beats exactly its first d_i successors, so the circular word recovers the tournament. Its automorphisms are the shifts preserving that word. Its anti-automorphisms (isomorphisms to the converse) are exactly reflections i↦r−i satisfying

    d_(r−i) = n−1−d_i   for every i.

Necessity is reversal of the intrinsic cycle and degree complementation. Sufficiency follows because reflection takes each consecutive block of d_i out-neighbors to exactly the d_i in-neighbors of the image. If a is the number of automorphisms, there are a such reflections in a self-converse class and none otherwise.

A Hamiltonian path cut from the canonical cycle starting at i and ending at i−1 gives a blue tiling (fixed by axis reflection in the fixed-path model) exactly when r=2i−1 is an allowed reflection. For odd n each allowed reflection contributes one cut. For even n, every allowed reflection has odd r: even r would fix two vertices and an anti-automorphism cannot fix both ends of an arc. Each allowed odd r contributes two cuts. Distinct path cuts give the same tiling exactly if related by a tournament automorphism; this action is free. Dividing the cut counts by a proves the exact rule:

| Canonical cycle cuts per class | not self-converse | self-converse |
| --- | ---: | ---: |
| blue cut tilings, n odd | 0 | 1 |
| blue cut tilings, n even | 0 | 2 |

These are **not** the total blue tilings; they are a canonically specified part of them. The bitwise tile flip reverses the closing arc while preserving the path edges, so it takes these tilings to ones whose marked path does not close. This distinguishes the even-order two-cut result from saying that the two cuts are a blue flip-pair. Arbitrary closing Hamiltonian paths need not cut the canonical cycle either.

## 3. Same net, same eight sites, same Moon profile; different symmetry and H

Take Z=z²/2 on the right half-plane, so h=z². Use eight fixed sites z_i=x_i(1+i·sqrt(−1)/100) for integer indices i=0,…,7, with x_i=2 at i=2,3; x_i=5 at i=1,4; and x_i=8 otherwise. Two admissible amplitudes on this one net are

    λ_A=|z−2|²−1/4,   λ_B=|z−5|²−1.

The signs at the sites are negative precisely at {2,3} and {1,4}, respectively; the argument of h increases inside a half-circle. Thus i→j for i<j exactly when the two sampled λ signs agree. All comparisons are nondegenerate; the continuously prescribed net is allowed to pass through isotropic loci between sites. This respects the fixed-site rigidity for everywhere nonisotropic fields in §1.

| Output | cyclic degree word | H | self-converse? | all blue tilings |
|---|---|---:|---|---:|
| A | 5,4,5,4,3,2,3,2 | 393 | yes | 23 (2 canonical cuts) |
| B | 5,5,4,4,3,2,2,3 | 347 | no | 0 |

The common sorted degrees are 2,2,3,3,4,4,5,5, so both have the same complete canonical Moon profile (16,40,44,26,8,1) at lengths 3 through 8. Yet forgetting the cyclic ordering loses the information needed to recover even self-converseness, let alone H. The source of the ambiguity is concrete: the two stress zero circles change different sign chambers without changing the prescribed principal axes or moving the sample sites.

Run **node replay.js** for exact integer stress determinants, explicit Hamiltonian path enumeration of the pair, and checks of the profile and blue cuts across 219 strong sign-word presentations through order eight. The universal assertions above follow from proofs, not the replay cutoff.
