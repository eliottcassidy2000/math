# Complementary physical and scope review of the general prime-five theorems

**PASS. No mathematical corrections.** This review covers `prime_five_all_profiles` and `prime_five_all_ambient`. It complements, rather than duplicates, the separate independent reconstruction of their complete local and ordered-cover kernel banks. Final source pins are recorded in the paired receipts.

## Analytic scope and normalization

The five normalized child types cover all nonempty tournaments. The singleton type denotes a physical one-name child, never an absent child. The C type's normalized Hamilton coordinate is one, while the actual canonical C3 has Hamilton count three. The same odd product scalar multiplies both quotient orientations and every cover coordinate; arbitrary profile-preserving target replacements retain that scalar. The physical controls restore it explicitly.

The local Hamilton-only truncation is valid. With one three-run child, its two nonempty return segments can be exchanged freely. With two repeated children, one of multiplicity three, each return segment contains at least two distinct other symbols: a one-symbol return would require a directed two-cycle. If exchanging the two segments fixed the word, at least two other symbols would have to repeat, whereas only one other child may repeat. Thus the exchange is free and the even coefficient times two even child weights vanishes modulo eight. This argument is used only within one Hamilton path, never across cover cuts.

The first-degree-per-port bounds survive identifying physical bits within each child: identifying variables can reduce degree but cannot increase it. Consequently the 256 child-query evaluations with each mask of weight at most two determine the complete 1,024-query coarse Boolean response, retaining every cubic and quartic coefficient. The separately retained full banks and their independent interpolation audit certify the stated finite classification; this complementary audit does not claim to reconstruct those banks again.

The three prime representative orbits are complete. This audit independently enumerates every labelled five-name tournament, tests every proper nontrivial subset for modularity and reconstructs all physical relabelings of the three representatives. Their disjoint orbit sizes are 120,120,24 and their union is exactly the 264 prime tournaments. Transporting quotient names, child sets and masks together is essential and is respected.

## Ambient cases and profile tails

The singleton branch is justified modulo four. Other nontrivial children have the same normalized q2 residue, and their three-run terms vanish modulo four. Thus a fixed singleton-position pattern has the same baseline congruence regardless of the other nontrivial child types. The 31 nonempty position patterns for each representative are complete. Same-order q2/q3 residues depend only on order, and higher factorials vanish modulo four, so a Hamilton contrast two has full profile contrast `(2,0,0,...)` and survives any fixed ambient.

The Q41 nontrivial branch has its separate modulo-four anchor. For Q10/Q76 with even total order, only the inherited even-order modulo-eight transfer is used. For regular Q76 with all-E children, the special modulo-sixteen full-cover certificate is required and explicitly invoked; no general modulo-sixteen lift is inferred. Its q4/q5 tails depend only on the common order and all h>=6 vanish modulo sixteen.

For odd total order, the 512 assignments are exactly the odd-parity half of the four nontrivial types on five children. The 51 retained run shapes contain all multiplicities at most three with at most two repeated children. Every omitted child-factor product is zero modulo eight. Cuts permit arbitrary transitions, including the same child across a cut, and every path slot is nonempty. The three-coordinate kernels therefore compute the complete first-three cover contrasts. All h>=4 coordinates vanish separately by h! divisibility. An actual contrast `(4,0,0,...)` is preserved by the odd first occupied-run coefficient in every fixed exterior.

These branches exhaust every nonempty child assignment. The new packets keep local separation and arbitrary-ambient separation logically distinct.

## Physical lift and family quantifiers

Both ports use the original core's incidence to every exterior name, while retaining edge x->y. The augmented core is a common module; no exterior name or arc is conversed. The four source ambient controls use a cyclic outer quotient with two singleton exterior names. Independently reconstructed masks give the exact integer identity

    H(C3[augmented core,1,1]) = 3 H(augmented core)+q2(augmented core).

This directly verifies the producer's residue formula, including all odd normalization factors.

The family conclusion is source-relative. Operations are restricted to original modules inside the fixed core, and the whole prime core is selected oddly. Proper original modules lie within its children, and hereditary monotonicity yields coordinatewise larger target child profiles, even for crossing descendant families. Whole-augmented-profile duality produces opposite baseline contrasts when only core bits are complemented; exterior port bits remain fixed. Choosing the negative source-minus-baseline member makes all target profile increases preserve strict separation. The chosen member may depend on the source and fixed ambient, but not on the particular permitted target family. No pairwise injectivity among arbitrary outputs, arbitrary exterior changes or unrestricted reconstruction is asserted.

## Independent executable controls

The auditor imports no producer implementation. It reconstructs every graph from the quotient rows and actual canonical child tournaments, checks all physical masks, and recomputes the five complete canonical profiles by literal words and cuts.

- All 18 local physical controls reproduce exactly modulo sixteen, with full integer q1/q2/q3 values retained here.
- All four ambient physical controls reproduce exactly, including their actual core q1/q2/q3 residues and the exact strong-ambient identity above.
- All 93 singleton-position certificates are independently instantiated using literal singleton/pair children. Each source-target difference has complete first-three profile `(2,0,0)` modulo four.
- All 264 labelled prime quotients and their complete orbit coverage are independently recovered.

The C++ engine uses an endpoint/predecessor recurrence, with an independent cut term for preceding path slots. It computes full integers rather than the producer's byte residues. Every graph has at most 17 names, and 17!*binom(16,2)<2^64 bounds the first three cover coordinates. Both Python modes use explicit checks and produce identical receipts.

Frozen audit allowlist: `REVIEW.md`, `README.md`, `audit.py`, `count.cpp`, `normal.json`, `optimized.json`, `SHA256SUMS`.
