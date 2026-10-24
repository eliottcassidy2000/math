# Bounded atomic decomposition of compatible paired-charge returns

**PROVED (conditional native theorem).** For a supplied jointly compatible binary source inventory with paired-charge group G=(Z/2)^4 plus Z/128, every two accepted assignments of that same inventory can be joined by atomic balanced packets changing at most 264 slots each. Every differing slot changes exactly once. The source-group bound is independent of the number of supplied slots. It does not assert that arbitrary native supplies satisfying the hypotheses exist, nor that the same group or constant persists at larger tournament orders.

## Exact interface

Each labelled slot has authenticated U/V alternatives. Sources lie outside U,V, and the two targets have equal paired-charge label. Every alternative retains the same entire declared shared-coordinate record. Every source and kernel lift is explicit. All before/after normalized resource orbits in the supply are mutually disjoint, and disjoint from the unchanged remainder. Electrical acceptance is exactly equality of target margins and the original physical-row paired charges; any additional constraints must also be preserved by every subpacket, not silently omitted.

For the original order-eight atlas the complete paired group is G=(Z/2)^4 plus Z/128. This is the audited result in `session45_injective_multisource_moves_20261023/dyadic-constraints/THEOREM.md`, sections 2 and 4. It encodes inclusion of every original dyadic harmonic bank, not equality of critical groups and not odd-primary guards. In particular the full exponent-256 filtration, not unscaled evaluation of one integer X at arbitrary higher moduli, supplies the all-dyadic interpretation.

Given accepted assignments a,b, remove unchanged slots. Let A be slots changing U to V and B those changing V to U. Equal margins imply |A|=|B|=k; equal physical-row charges imply sum ell(A)=sum ell(B). Pair A with B arbitrarily, retaining the actual slot indices. For pair i put d_i=ell(A_i)-ell(B_i). Then sum d_i=0 in G.

## The finite-group bound, with proof

For a finite abelian p-group G=direct-sum C_(p^a_j), put

    D=1+sum_j(p^a_j-1).

Every D-term sequence in G has a nonempty zero-sum subsequence. This is the classical p-group Davenport theorem (Olson); a primary accessible statement is Zhi-Wei Sun, *Zero-sum problems for abelian p-groups and covers of the integers by residue classes*, arXiv:math/0305369v8, Theorem 1.2: https://arxiv.org/html/math/0305369 . No priority claim is made here.

For completeness, work in F_p[G]. Write t_j=[e_j]-1. The group algebra is the truncated polynomial algebra F_p[t_j]/(t_j^(p^a_j)), so its augmentation ideal J has J^D=0: its largest surviving monomial degree is D-1. For D elements g_i, the product of (1-[g_i]) therefore vanishes. If no nonempty subset sums to zero, the coefficient of [0] in that product is exactly 1, contributed only by the empty subset, a contradiction. The standard sequence having p^a_j-1 copies of each e_j is zero-sum-free, proving that D is sharp. Adding the negative of its total gives a minimal zero-sum sequence of length D.

For our group D=1+4+127=132. Whenever more than 132 paired differences remain, find a zero-sum subsequence among the first 132; when at most 132 remain, the entire remaining zero-sum sequence is available. Execute exactly those pairs as one atomic packet and remove them. Each event changes at most 2D=264 slots and at least two. Every slot is used once, there are at most k integer events, and the full return has no extra source-slot work beyond its endpoint changes. Actual owner-clock work is the sum of the supplied per-slot path-change costs.

Every intermediate assignment retains the original target margins and charge. Jointly fixed shared records and the globally disjoint resource supply ensure that each mixed old/new assignment is an actual native packet. Thus the additive argument compiles into legal atomic events under the exact interface hypotheses. If resource alternatives overlap, or an equality changes with the choices, this final implication is not licensed.

## Sharpness of the charge-only budget, and its native boundary

Let g_1,...,g_132 be the minimal zero-sum sequence consisting of the four binary basis elements, 127 copies of the cyclic generator, and their negative total. Supply 132 positive slots with these labels and 132 negative slots all with label zero. Any nonempty balanced charge-preserving subpacket must choose a zero-sum subsequence of the positive labels, hence all of them, and then all negative slots. The accepted charge relation has exactly the empty and full changes. Thus 264 is sharp for the abstract interface allowing repeated labels at distinct slots.

This construction does not authenticate 264 native sources, distinct resources, common shared pins or a compatible remainder in the actual order-eight atlas. Distinct native sources can have equal charge labels, but not every abstract multiplicity has been supplied. The theorem is a conditional all-size decomposition rule, not an unconditional native existence theorem. The presently authenticated finite supply is much smaller.

The result strengthens the earlier compatible-supply pigeonhole criterion in a different direction: that criterion finds *some* equal-sum alternatives; here any specified pair of accepted assignments is connected. It does not weaken the known 16-source or new 18-source primitive examples: those require large simultaneous events inside their own smaller complete inventories, while our uniform upper budget is 264.

## Connection to the earlier ternary compiler

Session44's `ternary-newton/MULTIVARIABLE_FILTRATION.md` proves the exact augmentation-ideal index 2(k+r)-1 for (Z/3^k)[F3^r]. At k=1 its index 2r+1 is precisely D(F3^r). The common algebra is literal finite translations and their differences. Here vanishing of a product forces a subset-sum collision by inspecting the identity coefficient; there it controls the depth of repeated finite differences and retains the top norm/mass datum.

The source and target of this connection are different representations of the same group algebra. Preserved data are the additive group and its faithful coefficient arithmetic. Lost data are native ports, source ownership and charge-to-row maps; those require the explicit interface above. This explains why nilpotency is useful for bounded packet extraction but cannot by itself manufacture a native route.

## Computation

`zero_sum.py` uses a source-indexed subset-sum DP to split 24 seeded zero-sum returns with up to 1000 pairs, plus every collision among all subsets of a ten-source slice of the authenticated twenty-source menu. A separately implemented counting DP verifies that the 132-term extremal sequence has exactly two zero-sum subsets including the empty subset. No DP complexity or small experiment is used in place of the algebraic proof. `NATIVE18.md` gives the new complete native circuit beyond sixteen.
