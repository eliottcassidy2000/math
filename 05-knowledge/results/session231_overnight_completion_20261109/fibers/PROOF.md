# Movable-core carries and endpoint-erasing repair fibers

Status: proof candidate with exact certificates; unrestricted Chan–Pak completeness remains OPEN. No external literature-priority claim. Work is scoped to physical ordinal deletion blocks and star closure.

## Inheritance and concept board

Closest mechanism: Session229 `dynamic/PROOF.md` exact seven-coordinate composition and pure-T context law, `rich/PROOF.md` H29/H95/H389 witnesses, and `overlap/PROOF.md` bilateral carry rigidity and scalar-entropy barrier. The canonical hostile is multiplication of switch gaps by the Hamilton counts of other nontrivial modules. The corrected near miss is treating many choices or gcd-one as interval occupancy. The least-used sidecar is the starting/ending pair (S,E), which singleton padding can erase rather than preserve. Live concepts: fixed-position decoding, movable-core carries, erased endpoints, non-ordinal extension, and physical scalar occupancy.

Coordinate convention throughout: (h,q,A,G,T,S,E), as in Session229 `dynamic`, not the permutation used in `rich`. For nonempty ordinal blocks,

    compose((h,q,A,G,T,S,E),(j,k,B,J,U,V,W))
      = (hj,qk,Ak+hB,Jq+jG,AJ+Tj+hU+EV,Sj,hW).

The selected/unselected singleton has state s_b=(1,2,b,b,0,b,b). The formula is for NONEMPTY blocks. Omit an empty context instead of pretending (1,1,0,0,0,0,0) is a two-sided identity: its S,E coordinates do not obey the nonempty endpoint law. An exploratory padding implementation initially made that error; it was corrected before final results, and the checker includes the singleton control. No output from the invalid version is retained.

Star closure adds v,r with L→v→r and r beating the cut in L. Its count is D=2A+h−T. Selecting the first singleton of L ensures strongness: r reaches that vertex and then all later blocks, while every L vertex reaches r via v.

## 1. Fixed-position decoding (PROOF)

Let one physical core have h odd, q even. Let a prefix of a singleton vertices have binary cut number P, and suffix of b singleton vertices have binary cut number R, each read in forward order, most significant bit first. Then

    A_total=2^b(qP+A_core)+hR.

Since 0≤A_core≤q−h and 0≤R<2^b, A_total mod 2^b determines R uniquely (h is odd). Dividing by 2^b gives qP+A_core, whose quotient and remainder by q give P,A_core uniquely. Thus for a fixed split and fixed h,q, equal A_total fixes BOTH entire padding cut words and the core A. Equality of G_total then fixes core G. This needs no coprimality of h and q.

Consequently, fixed-position padding cannot combine different core A/G fibers. It can, however, erase S/E and enlarge a full six-coordinate fiber by combining masks which agree in A,G but disagree in S or E.

## 2. Finite automaton for moving a core (PROOF)

Comparing C followed by singleton y with singleton x followed by C' preserves A,G exactly when

    A'=2A+hy−qx,
    G'=(G+qy−hx)/2.                                  (M)

Here h is odd and q is even, so integrality forces x=G mod 2. There are at most two outgoing choices y∈{0,1}. For a move across d singleton positions, intermediate pairs need not be physical core interfaces. This is crucial: restricting them to the physical core bank would miss genuine long moves.

Every path starting and ending in physical core A/G values stays within

    −h < A_i < q,   −h < G_i < q.

Proof: after i moves, G_i=(G_0+q R_i−h P_i)/2^i, with binary integers 0≤P_i,R_i<2^i. The strict bounds follow from 0≤G_0≤q−h. Apply the analogous backward formula to A_i from the physical final A_d to obtain the A bound. Thus the finite directed graph on integer pairs in this open square is complete for all padding lengths. Conversely every path whose endpoints occur in the physical bank supplies actual ordinal constructions C_0 Y and X C_d, where X,Y are its recorded x/y words, with equal A,G. Intermediates are bookkeeping states, never claimed physical tournaments.

For arbitrary comparisons with different core positions, equal A_total first recovers their shared suffix bits modulo 2^min(b,b'). It also recovers their common prefix: for the earlier core position, the residual after q2^b P is in [0,q2^b−h]. Thus the later prefix starts with the same a bits. Cancelling common singleton contexts reduces the comparison to exactly one path (M). Hence this graph is an all-length recognition mechanism, not just a bounded search.

### Certified fixed-bank consequences

`verify.py` derives every core state twice, by subset DP and literal permutation-and-split enumeration. It completely explores (M) from every distinct physical A/G pair, retaining every edge returning to the physical bank, including edges whose destination was previously visited.

- Original B29=[26,44,57,16,34,9]: 64 physical A/G pairs, 236 total source-counted reachable states, no nontrivial return.
- B389=[20,85,72,195,204,159,33,71]: 240 physical A/G pairs, 10,460 source-counted reachable states, no nontrivial return.
- B95=[64,121,35,37,13,81,28]: 96 physical A/G pairs, 100,488 source-counted reachable states; 68 return ingress edges from 20 sources. The shortest nontrivial return has length10.

The only ignored returns for the first two banks are the all-zero cut with all-zero crossing bits and the full cut with all-one crossing bits. These preserve the entire seven-coordinate state, because for a full cut A=G=q−h and S=E=h. Therefore original B29 has at most one T value in every full-interface fiber under arbitrary singleton padding and arbitrary core position. B389 has at most two: its same-A/G core masks already have identical S,E and at most two T values. These are exact finite-certificate/all-length conclusions for the named single-topology banks, not for all h29/h389 cores.

## 3. A real H95 four-value movable fiber (FINITE-EXACT + context propagation)

Let X,Y be the two columns of this length20 path:

    (1,1),(0,1),(0,1),(1,0),(0,0),(0,0),(0,1),(0,1),(0,0),(1,1),
    (0,0),(1,1),(0,1),(0,1),(1,1),(1,1),(1,0),(0,1),(0,1),(0,0).

Use the four deletion graphs

    s1 ⇒ (B95,mask14 or22) ⇒ Y ⇒ s1,
    s1 ⇒ X ⇒ (B95,mask46 or54) ⇒ s1.

They have 29 vertices and common (h,q,A,G,S,E)

    (95,5351931904,4184607987,3843289869,95,95).

Their T values are {3005019510,3005019518,3005019550,3005019558}. Strong star closure gives four actual counts on31 vertices:

    5364196511, 5364196519, 5364196551, 5364196559.

The offsets are 0,8,40,48. This is an actual cross-position full-interface enlargement, not a mask-count claim. Arbitrary common singleton context preserves these differences. It does not supply consecutive odd values.

The shortest-return-path experiment in `grow_windows.py` also checks every physical intermediate core position on each returned path, with both possible final singleton bits and initial singleton bit1. It finds maximum fiber size4 and maximum consecutive step-two run2 for H95. It is NOT an exhaustive scalar analysis of all walks: cycles might change weighted behavior even when reachability has already been exhausted. Tested physical self-return loops have zero drift when iterated with their matching interior boundary bits; this is a hostile control, not an all-loop theorem.

## 4. Endpoint erasure supplies a different unsaturated H29 unit switch

Take the genuinely nonisomorphic six-vertex cores

    C=[28,41,58,48,34,1], cut19,
    D=[20,41,58,49,34,1], cut19.

Their states are respectively

    (29,320,129,129,34,11,13),
    (29,320,129,129,34,13,11).

Thus they have the same A,G,T but exchanged S,E. In s1⇒core⇒s0 the endpoints become S=29,E=0 and T becomes G+T+S. The resulting full interface is

    (h,q,A,G,S,E)=(29,1280,898,287,29,0),
    T=174 or176.

Strong star closure gives counts1651 and1649 on10 vertices. The switch is unsaturated, since gcd(29,1280)=1. It remains a gap-two pair in every common singleton context. Neither core is isomorphic to the original B29; checking all720 label permutations certifies all three pairwise nonisomorphisms. This explicitly escapes the single-topology bilateral-rigidity conclusion, without contradicting it.

## 5. Saturated five-value H45 fiber and its exact hole

For C45=[48,49,3,7,12,28] and D45=[40,37,17,6,11,28], the padded choices s1⇒core⇒s0 with (core,cut)=(C45,9),(C45,3),(D45,3),(C45,6),(D45,10) share

    (h,q,A,G,S,E)=(45,1728,1122,303,45,0),
    T in {164,166,168,172,174}.

Their strong star counts on10 vertices are {2115,2117,2121,2123,2125}. The missing2119 would require T170. `erase_endpoints.cpp` exhausts ALL32768 labeled six-vertex tournaments and all64 cuts, for each final singleton bit0/1. At h45,q432,A=G129, the local T+S set is EXACTLY {35,37,39,43,45}, so T170 cannot be supplied by another six-vertex core with that interface. The entire padded census has maximal consecutive run3, and its unsaturated subuniverse has maximal run2.

`CORE_BANK.jsonl` retains all286 distinct seven-coordinate states at (h,q)=(29,320) or(45,432), with literal graph/cut witnesses. Running the movable-core experiment across the COMPLETE h29 bank gives maximal run2, while the complete h45 bank still gives run3. These are shortest-return-path scalar experiments, not all-walk scalar classifications.

## 6. Non-ordinal two-vertex extension control

`extend2.cpp` exhausts all32768 orientations of two new vertices added to B95 (all14 old-new arcs and the new-new arc), and every cut on each9-vertex graph. Of these19586 extensions satisfy gcd(h,q)=1. Exactly22 such topologies have a repeated full interface within that topology. No within-topology fiber has a consecutive triple. Unit-gap witnesses occur only for codes0,127,16383,16384,32640,32767: the six ordinary ordinal placements of two singleton vertices. The other16 duplicate topologies have no unit pair. This extends, but does not repeat or contradict, Session229's one-vertex extension census. It does NOT group interfaces across different extension topologies.

## Reproduction and next obligation

Run `python3 -B verify.py > CHECK.json`, repeat with `-O`, and compare; literal core counting and subset DP must agree. Run `python3 padding.py`, `python3 grow_windows.py`, `python3 grow_all.py` for exploratory scalar outputs. Compile the C++ census files with `g++ -O3 -std=c++17`; their exact encoding loops are in source. All claims distinguish graph reachability from weighted scalar coverage.

The positive next question is an endpoint-erasing bank with a genuinely four-consecutive T fiber, or a weighted movable-core cycle whose compatible iteration has nonzero unit drift. Reachability alone cannot establish either. In the H45 bank the target hole is now explicit rather than a vague demand for more cardinality. The H29 pair demonstrates that changing core topology while erasing endpoints is a real escape from single-template rigidity.

## 7. Smaller movable carriers: H5 pair and H9 triple

The moving mechanism does not require the large H95 core. Let C5=[8,1,3,6], with h5,q32. Compare the two ordered deletion objects

    s1 ⇒ (C5,cut4) ⇒ s1 ⇒ s1 ⇒ s0 ⇒ s0,
    s1 ⇒ s0 ⇒ s1 ⇒ s1 ⇒ (C5,cut2) ⇒ s0.

Their full interface is (5,1024,716,209,5,0), and T is142 or140. Thus strong star closure gives1295 and1297 on11 vertices. This is an unsaturated pure-T switch assembled from one familiar four-vertex core, not an arithmetically new seven/eight-vertex core.

Let C9=[16,1,3,7,14], with h9,q84. The three deletion objects

    s1 ⇒ (C9,cut13) ⇒ s0 ⇒ s1 ⇒ s0,
    s1 ⇒ s1 ⇒ (C9,cut4) ⇒ s1 ⇒ s0,
    s1 ⇒ s1 ⇒ s0 ⇒ (C9,cut22) ⇒ s0

have full interface (9,1344,1090,427,9,0) and T334,338,336. Their strong star counts are1855,1851,1853 on11 vertices. This is saturated (gcd(9,1344)=3), but it is a moving SINGLE-core realization of a universal triple. `verify_repairs.py` checks each actual closure by independent forward Hamilton subset DP and verifies strong connectivity, rather than merely substituting the claimed scalar formula.

## 8. Complete seven/eight-vertex endpoint-erasure census

The inherited Session230 `primes/classgen.c` is copied unchanged and regenerated. It returns456 classes at order7 and6880 at order8. The complete order7 result is independently cross-checked against ALL2^21 labeled graphs, using separate graph encodings. Every cut is included. Both report104432 erased interfaces, maximal consecutive step-two run3, and maximal unsaturated run2. Order8 reports3385144 erased interfaces and the same two maxima. Consequently no original full-interface fiber through order8 contains four consecutive step-two T values either: adding the same selected prefix and unselected suffix preserves such a fiber and translates T by the common G+S.

The completed order9 bucketed census covers191536 classes, every512 cuts and both suffix bits. It reports193541720 erased interfaces, maximal consecutive run3, and unsaturated maximum2. No inference is made beyond order9. The earlier order7/8 interface counts accidentally included one empty map entry created by a missing-target diagnostic. The diagnostic now uses `find()`, and reruns give the corrected104432/3385144 counts; extrema and witnesses never changed.

Completeness is by inductive vertex extension and canonicalization: every n-vertex tournament deletes to a represented (n−1)-vertex class; every orientation of its added vertex is generated; canonization maps isomorphic copies to one representative. The canonicalizer retains actual permuted adjacency codes and exhaustive individualization/refinement branches, so equality of final codes is an actual isomorphism, not a hash approximation. The parent reports an independent degree-canonical audit of the complete enumeration. `ORBIT_MASS.json` also verifies sum(n!/aut)=2^(n choose2) and expected counts; orbit mass is a secondary control, not a substitute for canonization validity. The independently exhaustive labeled-order7 scan agrees exactly with the class-based scan.

### Residue-support hostile

The conjecture that an erased interface always misses an even residue mod8 is already false at order8. For the classgen graph/cut witnesses

    (38838976,24),(9478848,20),(38838976,20),(38838976,18),(9478848,18)

one has common h477,q6192,A=G1473 and local T+S values321,345,351,381,387. After s1⇒core⇒s0 the common full interface is

    (477,24768,15330,3423,477,0),
    T={1794,1818,1824,1854,1860}.

All four even residues mod8 occur, but the actual gaps are24,6,30,6. The missing interval coordinate is therefore not explained by this simple modular support test. The complete order6/7 erased universes have no such residue-support counterexample; order8 is the first in this census.
