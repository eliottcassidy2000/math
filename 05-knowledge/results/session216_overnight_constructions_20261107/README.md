# Session216: overnight tournament spectrum constructions

**Extended session complete; the Chan–Pak conjecture remains OPEN.** The target is every odd Hamiltonian-path count at least 23. This index separates all-order theorems, exact finite certificates, restricted obstructions and the remaining proof obligation. No unrestricted finite-range extension, positive-density theorem or eventual threshold is claimed.

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

## 5. Complete two-apex census and ternary rays

[Certificate, proof and replay](two-apex-advances/README.md), [independent audit](two-apex-advances/audit/AUDIT.md), [release review](two-apex-advances/RELEASE_REVIEW.md).

Every odd count through2,327,673 other than7 and21 is realized by a literal two-apex tournament on at most16 vertices. The complete mask-word census through14 chain vertices has first additional omission2,327,675; the explicit15-letter word120021002033010 realizes that boundary count. This establishes the minimum chain length for that target within this fixed class, not minimum tournament order globally. It does not improve the earlier unrestricted range.

The retained certificate contains1,163,835 explicit words. All were checked by an independent36-coordinate transfer;297 certificate entries and the boundary repair were independently checked by subset DP. The full357,913,941-node census was independently rerun. Lossless certificate parts reconstruct the original gzip with exact hashes.

Three exact all-length ternary families are60*3^m-1,120*3^m+1 and72*3^m+1 for m>=0. Bounded global-map and wrapper probes do not yield a digit induction. Their scopes are stated explicitly in the proof packet.

Published [checkpoint2](https://github.com/eliottcassidy2000/math/commit/a1ff23f5b22dcacf6bcc5b6d6b71a25683b1b7e2) contains Sections2–4.

## 6. Three fixed cores give a concrete arithmetic coverage target

[Construction and binary formula](three-core/README.md), [portable certificate replay](three-core/REPLAY.md), [independent audit](three-core/audit/REVIEW.md).

Three fixed strong cores, one with H=5 and two with H=9, provide only 80 core-cut states. Arbitrary cut-labeled transitive prefix and suffix padding followed by a fixed star closure realizes every odd count from 23 through 1,000,000 except 47 and 63, with at most 16 padding vertices. All 499,987 retained recipes and the complete 88,735,747-choice scan were independently checked, together with 3,920 small literal constructions. The largest retained tournament has 22 vertices. Separate explicit strong tournaments supply 47 and 63.

A positive-cover inequality proves that those two exceptions are absent at EVERY padding length within this particular grammar. The same bound gives a finite exact all-length membership algorithm for any specified target. These class-relative exceptions are not missing tournament counts.

The sufficient arithmetic conjecture that these three templates realize every odd N>=65 remains OPEN. Proving it, together with the already verified small cases, would settle the Chan–Pak target. Finite coverage and exact membership algorithms do not prove this infinite occupancy statement. This packet is a compact constructive grammar, not a new unrestricted coverage record.

## 7. Exact registers, persistent adjacent pairs and inverse numeration

[Register theorem and audit](register/README.md), [numeration proof, inverse and audit](numeration/README.md).

With the H=5 core fixed and a nonempty initial cut, selected/unselected sink appends realize exactly the register values alpha>=5 except 8. The joint endpoint-count support is also classified. The Hamilton count after closure still depends on a second coordinate, tau; register completeness is not count completeness.

Two small physical sources close to counts 49 and 51 and preserve their difference after every common append word. A stronger same-order pair closes to 4633 and 4635 and preserves adjacency under arbitrary common cut-labeled transitive padding on BOTH sides. These are returning interfaces, not scalar operations on arbitrary realized counts, and they preserve width rather than prove unbounded unit-gap amplification.

For the three-core formula, deterministic alpha stripping removes the free suffix search. A balanced-prefix lemma says that q^2*4^b>hN leaves at most one positive prefix candidate for a fixed suffix. Exact recipes repair the short-side scout omissions 4,295,971 and 9,948,407. Narrow balanced bands still have finite gaps, so no digit induction follows. Large retained graphs are certified by the exact independently evaluated construction formula and explicit Hamilton cycles; only the stated smaller controls use whole-graph subset DP. The optimized C++ band search requires its documented resource ceiling and signed64-safe intermediates.

Published [checkpoint3](https://github.com/eliottcassidy2000/math/commit/656b8a9248b1534d3a868026356edcdb6d4eb4a1) contains Section5.

## 8. Rooted reciprocity narrows the universal paired-source route

[Carrier theorem and rooted reciprocity](carriers/README.md), [reconstruction instructions](carriers/START_HERE.md), [complete order-ten audit](order10/README.md), [order-ten reconstruction](order10/START_HERE.md).

The rooted cover-polynomial identity F(1-k)=(-1)^(n-1)G(k) is proved in a squarefree ring. Its integer coefficients rule out every SAME-order-parity universal O/Z/J pair whose start counts differ by 2. A universal width-eight pair must instead have opposite order parity, even starts in one of (4,2), (10,8), (16,14), (22,20), and all eight shared nonstart coordinates divisible by 3. The later two-adic refinement gives C and b divisible by 24 and L,T congruent to 6 modulo 12.

These are necessary conditions for this particular universal source criterion, not restrictions on arbitrary tournament counts. They do not prove that every universal pair is impossible. The complete independently audited census proves only that no such width-eight pair has BOTH orders at most 10. Canonical augmentation verifies all 191,536 order-nine deletion classes; all 1,412 surviving order-ten profiles were independently recounted. The independent root compiler omits the producer's minimum-start pruning and still recovers exactly the same inventory.

Separate finite asynchronous-word and right-amalgamation probes have their own declared banks and depth limits. They are not all-order exclusions. The common right-factor theorem explains why a nontrivial common deletion factor cannot repair an already synchronized pair; it does not exclude asymmetric substitutions or the left-reset mechanism below.

## 9. Unequal-deletion reset pairs and a conserved defect

[Positive left-reset mechanism and structural controls](paired-reset/README.md), [reconstruction and local dependencies](paired-reset/START_HERE.md).

Adding z with root->z->all old nonroot vertices gives a physical reset operation Q. Explicit sources with unequal deletion counts 1 and 5 retain paired adjusted bases under every common repetition of Q. Their compiled counts are

    216*2^k-7 and 216*2^k-5,     k>=0.

This is a new composable source mechanism, not a claim that the displayed values were previously unattained. A five-coordinate interface gives a broader sufficient criterion for synchronization under every binary Q/Z word. No physical fully binary seed satisfying that criterion was found.

For nontrivial rooted sources in the stable state s=e=d, the nonnegative integer

    kappa=(h-a-rho+b-3d)/2

is conserved by Q and Z. Different deletion counts can in principle be compensated by different defects. The singleton is explicitly excluded: its value is -1/2. This boundary correction is preserved with the independent audit and original manuscript.

The sharper two-adic congruences form a closed class, so they do not force an infinite divisibility contradiction. Explicit tournaments of orders 17 and 18 share H=2025 and B(2)=5,971,968, disproving a tempting global-moment parity obstruction. Bounded macro and structured counter-rotation searches remain finite negative evidence only.

## 10. Final proof frontier and checkpoints

The strongest compact positive target now is the three-core formula in Section6: prove it represents every odd N>=65. Its 80 fixed cut states, explicit binary formula and exact all-length inverse make this a concrete sufficient arithmetic conjecture. Proving membership separately for many targets, even with certified bounds, does not prove universal occupancy.

The other open route is a physical binary paired-source family with adequate interval coverage. The all-order reciprocity restrictions eliminate same-parity universal seeds, while the Q/Z interface allows a different deletion/defect balance. A thin dyadic ray, a fixed-width adjacent menu, or normalized width with growing physical spacing does not supply that coverage.

Published checkpoints:

1. [Fixed-apex inverse and literal baseline](https://github.com/eliottcassidy2000/math/commit/94e3ca9e045266b3f71351b799415f820a9bf268).
2. [Star constructions and logarithmic density refinement](https://github.com/eliottcassidy2000/math/commit/a1ff23f5b22dcacf6bcc5b6d6b71a25683b1b7e2).
3. [Complete two-apex census and ternary rays](https://github.com/eliottcassidy2000/math/commit/656b8a9248b1534d3a868026356edcdb6d4eb4a1).
4. [Three-core arithmetic, registers and exact inverses](https://github.com/eliottcassidy2000/math/commit/79c2e3666ba1c78aecb2f0868933824a46ecadd2).
5. The current final packet adds the carrier/reciprocity, complete order-ten and left-reset results, together with this synthesis.

All promoted result packets carry independent audits, exact input freezes, replay instructions and checksums. Large retained certificates are losslessly split; reconstruct them as directed before replay. Historical frozen copies may retain their original pending-audit labels; the current release reviews govern promoted claims. No full Chan–Pak proof, eventual threshold, positive-density theorem or unrestricted finite-range record is claimed.
