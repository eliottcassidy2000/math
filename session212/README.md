# Session212: cycle counts are complete, while cycle sizes retain gaps

**Status: PROVED elementary all-order constructions and transfer lemmas; FINITE-EXACT censuses through order nine; explicitly INHERITED order-ten census. No literature-priority claim.**

Hamiltonian cycles are directed spanning cycles counted modulo cyclic rotation. H counts all directed Hamiltonian-path words. These are different observables: the all-order cycle-count support is the nonnegative integers, while the unrestricted Hamiltonian-path spectrum conjecture remains open and its omissions 7 and 21 remain valid.

## 1. The global answer and its structural strengthening

Every nonnegative integer is a tournament Hamiltonian-cycle count. Every positive integer has a strong realization. For every integer at least two, one can simultaneously require a sole-exit vertex, no degree-defined gate, and feedback vertex number exactly two.

The [binary construction](cycle-bridge/PROOF.md) uses a transitive chain and two marked extra vertices. Its two orientations of one apex arc give

    HC(B-) = q + tau,       HC(B+) = q + tau + 1.

At every order N>=5, these families realize the entire interval 0 through 2^(N-4)+1. The nongated positive interval supplied at that order is 1 through 2^(N-4); the extra top value in the forward family has a gate. For k>=2 a strong witness exists at every order

    N >= max(5, 4 + ceil(log2(k-1))).

The nongated sole-exit variant has the slightly weaker guaranteed threshold max(5,4+ceil(log2(k))). Neither threshold is claimed to be the exact minimum order.

This is a constructive all-order theorem: a marked insertion preserves every old arc and implements HC->2HC or HC->2HC+1. It does not infer universality from a numerical initial interval. Deleting the two apices leaves a transitive tournament. Conversely, deleting one feedback vertex leaves a unique possible Hamiltonian-path remainder, so feedback vertex number at most one permits at most one Hamiltonian cycle. This proves the sharp feedback threshold.

The [endpoint companion](cycle-bridge/ENDPOINT_BRIDGE.md) also realizes every positive integer as a specified-start Hamiltonian-path count of a strong tournament, and as a specified-start-and-end count of a strong tournament. Those marked counts do not equal total H in general.

## 2. Independent routes and inherited ingredients

The [selected-endpoint construction](cycle-spectrum/PROOF.md) independently reaches the same complete interval and logarithmic bound. A tournament with one reversed extremal arc supplies a binary endpoint vector; a forced-entry attachment selects any desired subset sum. It additionally supplies a strong bridge. This and the binary-comb construction are complementary proofs of one global count-support theorem, not separate novelty claims.

The interfaces are inherited and explicitly pinned:

- Session210 already proves that a degree gate a->b, d+(a)=d-(b)=1, gives the exact cycle/path correspondence HC(Gamma(T))=H(T), with size shift two.
- THM-4862 already gives the sole-exit attachment converting one endpoint-path count into HC.
- The selected binary endpoint vector was already recorded in the November 4 Chan-Pak increment audit.
- Session211 already supplied ungated examples with HC=7 and HC=21, and showed that scores do not recognize the gated class.

The present packets recover these mechanisms, prove unbounded target coverage with explicit grading and marked recursion, and distinguish total counts from endpoint restrictions. The recognizable-gate cycle spectrum still coincides with the open Hamiltonian-path spectrum. This session does not solve that restricted problem.

## 3. Fixed-order spectra are not generally intervals

The [complete spectra packet](cycle-spectrum/README.md) and the independent [cycle-level census](cycle-levels/CYCLE_LEVELS.md) agree through order nine:

- Orders 3 and 4: {0,1}.
- Order 5: 0 through 3; order 6: 0 through 5.
- Order 7: 0 through 15, together with 17 and 24. Thus 16 and 18 through 23 are missing at this order.
- Order 8: 0 through 54.
- Order 9: 0 through 249, except the 21 explicitly listed holes in both packets.
- Order 10: 0 through 961, together with 971. This is a checked extraction of the inherited complete THM-4889 census, not a fresh order-ten enumeration. The supplied maximum witness was independently recounted.

Fresh order-nine coverage uses every one-vertex extension of every order-eight isomorphism class: 6,880*256=1,761,280 presentations. Repeated presentations do not invalidate support coverage, but their frequencies are not isomorphism-class frequencies. Independent labeled enumeration through order seven and independent implementations further check the finite results.

## 4. Exact inverse size sets, not just eventual existence

The [inverse-size proof and witness inventory](cycle-levels/CYCLE_LEVELS.md) determine the complete all-order size set for every cycle count k from 1 through 593:

    {N: an N-tournament has HC=k} = {m_c(k), m_c(k)+1, ...}.

The minimum m_c(k) is given explicitly by the finite spectra, and every subsequent size is certified by a separate rooted endpoint witness plus count-preserving chain padding. The argument checks any intervening sizes rather than assuming a minimum-order cycle witness has a paddable vertex. Another 25 listed values above 593 have the same exact minimum-ten/full-tail certification.

This conclusion is not merely a bounded census. The finite evidence identifies an endpoint witness; the proved closure and padding then give an infinite tail. In contrast, adding a universal source or sink destroys Hamiltonian cycles. No global spectrum-nesting or universal minimum-order upper-ray theorem is claimed.

## 5. Replay and evidence

Full commands, universes, source pins, and controls are in each packet:

- [Cycle spectrum](cycle-spectrum/README.md): generate order-eight classes, run the C++ fixed-cycle and extension censuses, then run `python3 session212/cycle-spectrum/verify.py` and its `python3 -O` counterpart. The optional inherited extraction reads pinned git objects without modifying the checkout.
- [Unified bridge](cycle-bridge/README.md): `python3 session212/cycle-bridge/verify_binary_comb.py` checks 2,040 family cases, 248 literal growth steps, and 75 inherited-gate controls. `python3 session212/cycle-bridge/verify_endpoint.py` separately checks 511 endpoint instances, both apex orientations, 62 complete-permutation cycle instances, tiny boundaries, and strong fixed-endpoint realizers. Both normal and optimized runs agree.
- [Cycle levels](cycle-levels/CYCLE_LEVELS.md): `python3 session212/cycle-levels/verify.py` replays 1,166 saved witnesses and 1,779 physical closure/padding instances. Its companion C++ program exhausts every labeled tournament through order seven. Normal and optimized Python outputs agree.

The displayed Python commands are run from the repository root; full compilation and regeneration commands specify their packet-local working directories. Finite checks verify code and conventions. Infinite conclusions follow from the explicit bijections and padding proofs.

## 6. What remains open

1. Determine exact minimum orders for general cycle counts and whether every count's possible orders always form an upper ray from its minimum. This session proves that for the stated 593-count range and additional listed values, and proves eventual tails for every count.
2. Explain and classify fixed-order gaps beyond the verified spectra. Order-seven and order-nine gaps disappear at the next observed order, but this is not an all-order filling law.
3. Determine what cycle counts are possible under stronger connectivity or other restrictions that exclude the forced endpoint interfaces used here.
4. Solve the separate Hamiltonian-path spectrum completeness question. Arbitrary endpoint counts and complete cycle support do not supply an unrestricted path-count constructor.
