# Unsaturated four-choice modules and portable unit switches

**Status: PROVED construction and propagation statements; FINITE-EXACT physical witnesses. No unrestricted Chan–Pak coverage claim.**

## Inheritance and live mechanisms

The closest inherited mechanism is Session228's H19/Q256 equal-alpha pair; Session229 `overlap/PROOF.md` gives nearest-neighbor defects, coprime digit rigidity, and fixed-template entropy. The hostile is its binary menu versus radix 256/19. The corrected near miss is counting distinct masks rather than distinct transfers. The least-used sidecar here is the full endpoint signature (A,G,S,E), inherited from Session216 `three-core/PROOF.md` and `numeration/PROOF.md` under `05-knowledge/results/session216_overnight_constructions_20261107/`. Those packets already allow arbitrary binary prefix/suffix words. We do not present that grammar or its inverse search as new. During this session the context-stable switch construction independently supplied a saturated H15/Q120 unit switch; the present packet removes saturation with different physical blocks.

The concept board is accidental alpha equalities, rigid tournaments, full endpoint signatures, central transfer defects, and singleton context propagation. Anchor: enlarge genuinely different unsaturated choices. Niche: isolate context-independent gap two. Wildcard: explain why preserving symmetry is incompatible with unsaturation.

## Definitions

Outgoing-bitmask notation has bit j of B[i] equal to one exactly when i→j. A cut mask selects vertices beaten by an added root. H counts Hamilton paths and H(empty)=1. Q counts ordered pairs (P,R) of directed paths covering the block, allowing either path empty. A counts such covers with R nonempty starting in the cut; G counts covers with P nonempty ending in the cut; S and E count full Hamilton paths starting and ending in the cut. T counts covers with both P ending and R starting in the cut. Coordinates below are always (H,Q,A,G,S,E,T), not the different ordering used in some inherited files.

For ordered composition B⇒C, with coordinates (h,q,A,G,S,E,T) and (k,r,B,J,V,F,U), the inherited exact composition rule is

    (hk, qr, Ar+hB, Jq+kG, Sk, hF, AJ+Tk+hU+EV).       (1)

This follows by separating the empty/nonempty block portions of each cover. In particular, prefix/suffix singleton blocks have h=1, q=2 and coordinates (1,2,b,b,b,b,0), b in {0,1}.

## 1. A genuinely four-choice unsaturated module

Take

    B29 = [26,44,57,16,34,9].

Then H=29, Q=320, gcd(H,Q)=1, and the six A weights are

    [50,50,86,31,31,43].

The COMPLETE A=81 fiber has masks 9,10,17,18. Their (G,S,E,T) are respectively

    9:  (129,7,15,22)
    10: (117,7,13,16)
    17: (93,7,9,12)
    18: (81,7,7,14).

Thus these are four distinct physical transfers, not four isomorphic masks. With strong singleton flanks p⇒B29⇒q and root→p, q→root, the four literal nine-vertex Hamilton counts are 673,691,719,729. For these flanks the general scalar formula is

    D = 2Q+H+2A-G-S-T.                                 (2)

It follows either by the six-state transfer or by enumerating which paths contain the flanks. The checker independently counts the literal graphs by subset-path DP.

**All-order amplification.** The ordered composition of k copies of B29 has H=29^k, Q=320^k, hence remains unsaturated. Choosing independently one of the four masks at every copy gives a common A and at least 4^k distinct full transfers. Indeed

    A = 81 sum_(i=0)^(k-1) 29^i 320^(k-1-i),
    G = sum_(i=0)^(k-1) G_i 320^i 29^(k-1-i).

If two G words agree, reduction modulo 320 forces their first G difference to be zero: 29 is invertible modulo 320 and each difference has absolute value at most 48. Cancel 320 and induct. All G words are therefore different. This provides arbitrarily rich, rigid, unsaturated modules with a common A. It does NOT make the menu rich enough for an interval theorem: 4^k is still smaller than the effective radix (320/29)^k, and Session229's entropy barrier still applies to repetitions with fixed strong flanks. Furthermore all S values agree locally, so this four-choice family has no adjacent quadratic defects of its own.

## 2. Unsaturated full-interface unit switches

Here are two independently checked physical witnesses.

    B95 = [64,121,35,37,13,81,28]
    masks 12 and 20
    (H,Q,A,G,S,E) = (95,1276,312,270,24,18)
    T = 38 and 36.

    B389 = [20,85,72,195,204,159,33,71]
    masks 47 and 181
    (H,Q,A,G,S,E) = (389,5600,3369,3456,250,271)
    T = 1858 and 1860.

Both gcd(95,1276)=1 and gcd(389,5600)=1. These are stronger than same-A pairs: the entire endpoint interface agrees and only T changes by two. For the star closure L→v→root with root→cut, the inherited formula D=2A+H−T gives counts 681/683 for B95 and 5269/5267 for B389. Both cores are strong, so any nonempty cut gives a strong star closure.

**Context theorem.** If two cut blocks share (H,Q,A,G,S,E) and have T difference δ, then in any fixed ordered left and right context they continue to share those six coordinates, while their T difference is

    δ H(left) H(right).                                  (3)

Empty context contributes H=1. Apply (1) once on either side and iterate. Therefore arbitrary transitive prefix/suffix blocks preserve the EXACT difference δ, with arbitrary independently chosen cut bits. In particular the unit switches above yield consecutive odd star counts in every such context whenever that context is strong. Their strongness legality is identical: the deletion graph is unchanged and the cores are strong with nonempty cuts, so the inherited rule requires the first prefix vertex selected when a prefix is present, otherwise the core cut nonempty. This is a pair-completion theorem, not a theorem that every odd target belongs to one of those pairs.

**Central defect interpretation.** In the inherited six-state convention let w=e3−e5 and rho(z)=z0+z2+z4. The transfer difference of two such full-interface blocks is

    F = −δ w rho.                                        (4)

The only changed cover categories are (1,1):−δ, (2,1):δ, (2,2):−δ, (1,2):δ. Substitution in the update rule gives (4). Every nonempty physical context M satisfies Mw=H(M)w and rho M=H(M)rho. Hence MF=FM=H(M)F, and products of two central defects vanish because rho(w)=0. This explains why context propagation has no quadratic interaction. With k independent identical unit switches, the scalar variations are exactly 2mH^(k−1), m=0,...,k; an exponential menu must NOT be inferred from the 2^k words. Unsaturation removes a compulsory divisor of every scalar, but does not remove this spacing of switch differences.

## 3. Symmetry obstruction and bilateral rigidity

**Automorphism lemma.** For any nonempty tournament B, |Aut(B)| divides both H(B) and Q(B). An automorphism fixing a Hamilton sequence fixes every vertex. An automorphism fixing an ordered pair of covering path sequences likewise fixes every vertex. Both actions are free; orbit counting proves the divisibilities. Consequently gcd(H,Q)=1 forces B rigid. The converse is false in general; no converse is used. Preserving nontrivial exact symmetry cannot produce the desired unsaturated modules. The equal weights of B29 are accidental arithmetic equalities in a rigid block, not an unaccounted automorphism.

**Bilateral carry rigidity.** In a repeated block with gcd(H,Q)=1, equality of the aggregate A fixes the whole local A word, and equality of the aggregate G fixes the local G word. The latter follows from the same modular argument in reverse order; 0<=G_i<=Q−H<Q supplies the strict bound. Thus B29's four-choice menu cannot be combined across positions to cancel its differing G coordinates while preserving the full interface. This explains why the separate full-interface B95/B389 mechanism is useful rather than merely cosmetic.

## 4. Exploration and exact verification

These witnesses came from EXPLORATORY finite searches, not an exhaustive classification of all sizes or an optimized record. `search.py` has seed 2291209 and 1080 generated tournaments: n=7,...,12, 180 trials each, first 120 from perturbed odd circulants (or an even deletion), the remaining 60 uniformly labeled tournaments. `sparse_search.py` has seed 2291210 and 3600 near-transitive samples: n=7,...,12, 600 each, reversing 1+(trial mod6) sampled edges. Its n7 trial411 contains B29 followed by a sink. `interface_search.py` has seed 2291211 and 4000 uniformly labeled samples cycling n=7,8,9. B95 occurs at trial891 and B389 at trial3046. `exhaust6.py` separately scans all 32768 labeled six-vertex tournaments as a bounded exact search; no all-size conclusion follows. Logs and hit files are discovery provenance, not dependencies of the proof checker.

Run `python3 verify.py > CHECK.json` and `python3 -O verify.py > CHECK.optimized.json`, then compare byte-for-byte. The checker uses only the standard library, literal permutation enumeration of path covers (independent of the discovery subset DP), literal subset Hamilton DP for closure counts, all 4^k G words through k=7, all 720 permutations to verify B29 rigidity, and all length-four singleton cut contexts at every split for both unit switches. No Python assertions are used. The all-order claims follow from the proofs above, not finite tests.

The concrete gain is a four-way unsaturated transfer menu and two unsaturated full-interface gap-two switches, including a central-defect propagation law. The missing global step remains an overlap/coverage theorem across varying cores or templates; this packet does not solve it.

### Scoped extension follow-on

`extend_triples.py` exhausts every one-vertex extension of the context-stable triple cores [48,49,3,7,12,28] and [40,37,17,6,11,28]: 64 orientations for each core and all 128 cut masks. It filters gcd(H,Q)=1 and groups by all six interface coordinates across both core topologies. Exactly 48 of the 128 extensions pass that filter, producing 1772 interfaces. None has a pair of T values differing by two (maximum step-two run length one). `EXTEND_TRIPLES.json` records this FINITE-EXACT negative result. It excludes only this explicit extension universe, not seven-vertex cores generally: B95 already supplies a seven-vertex unsaturated pair outside it. Thus merely adjoining one vertex to those two fixed H45 templates does not transport their universal triple to the unsaturated setting. This bounded follow-on is a stopping reason for that operation, not a global obstruction.
