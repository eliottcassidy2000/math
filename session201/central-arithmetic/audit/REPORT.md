# Independent audit: bounded central-word arithmetic

**PASS for the physical compiler, exact finite inverses, finite-menu normalization, graph-order bound, and conditional implication.** The assertion that the length-eight family contains every odd target above 2647 remains **OPEN**. Exact finite coverage, even at the larger replayed bound, is not an all-target proof.

## Frozen sources and scope

`INPUTS.json` records the original eleven-file core plus the two supplemental reviewed dependencies. The core stayed unchanged throughout review. Principal hashes:

- `PROOF.md`: `797ff54fae332b917456a38330d659db2fce178b0732b49799c3e89d6a74dc0c`.
- `verify.py`: `a7cda7cc9dabdd8adb17f2101753157886e2d4b92bab0df2964e75db0a25ce5f`.
- `inventory.cpp`: `a2a5b627d22864f37de28ec32fde825d5b7169d6900b4b2ad4eacad45711798c`.
- Generated `central_forms.h`: `105eb8f7f929560097d9c96b267615f004ef744f668b68e5afbce5f8ca703844`.

The physical-adjoint dependency supplies the actual O/Z/J grammar and C1 characterization. The physical-existence dependency supplies the separately reviewed tighter source-pair inverse, not occupancy. The older forty-template packet explicitly keeps its own different occupancy lemma open. The monotone-C1 dependency and actual adjacency lists supply provenance for the small repair targets. The independent audit recounts those graphs directly.

The author's separately assembled exploratory archive is outside this source approval. Any later status promotion or additional scan should preserve its own exact inputs and limits.

## 1. Actual unit sources and physical compiler

For each integer u>=2, inverse cover recursion selects a unique binary source word: even u ends in O with predecessor u/2, and odd u>=3 ends in Z with predecessor (u+1)/2. This agrees with the author's u-1 binary-digit prescription. Every predecessor remains at least two and decreases until the seed. Thus tau(u) is an actual source, with g and parity fixed by u; there are no independently adjustable profile coordinates.

The word has length ceil(log2 u)-1 and source order 1+ceil(log2 u), including powers of two and the seed u=2. The excluded singleton u=1 is not silently used.

With chronological action matrices, the proposed graph word is exactly `B(v) W dagger(B(u)) J`. Adjunction turns its precursor F readout into `tau(u)^T Q M_W tau(v)`. The closing J is the inherited C1 construction. Its exact graph order is `1+ceil(log2 u)+ceil(log2 v)+|W|`. These are equalities of actual constructions, not cone or scalar-surjectivity assertions.

Conversely, strip the maximal initial and final O/Z blocks of any source word. The middle is empty, J, or begins and ends in J. Absorb the initial block into v and the dagger of the final block into u. This gives the stated unrestricted-central-length parameterization. Nothing in this argument bounds the remaining middle by eight.

## 2. Matrix certificates and complete inverse domain

Independent generation recovers 1095 central words through length eight and exactly 592 matrices. Matrix representatives are stored at their minimum central length. The full family is closed under dagger, and `A_W^T=A_Wdagger` is exact. Thus swapping u and v is licensed for the union of this complete transpose-closed family.

The cone coefficient calculation is exact. Each u0-coordinate is one or two, so multiplying the ratio by four makes every coefficient integral. The entrywise minimum proves `Psi_W(u,v)>=(b_W/4)uv`. Moreover the physical generators preserve the nonnegative cone and do not decrease U on that cone, so the inherited half-product comparison gives b_W>=2 for all words. The finite checker independently certifies every displayed matrix.

Consequently every representation of n in the complete family can be ordered u<=v, possibly replacing W by its dagger. Its indices then satisfy the stated isqrt and floor bounds. This is a necessary complete search domain, not a guarantee that it contains a solution. Completeness is target-dependent and all-height for the fixed finite menu; it does not depend on an empirical index cap.

For fixed W,u, the row coefficients alpha and beta are strictly positive for every u>=2. The nonnegative-column certificates in P^T A_W, with positive e and p rows, prove this because an actual unit source has e+p=1. The arithmetic progression filter is exactly the congruence obtained after writing v=2t+q and eliminating alpha*g(v) modulo alpha. Gcd divisibility, the reduced modular inverse, parity-dependent interval endpoints and modulus one are handled correctly. Rechecking the full integer equality is indispensable and is present. A congruence hit alone is never accepted.

## 3. The normal menu is an exact directed-index reduction

The local identity JZ=OJ holds on all four raw coordinates. Its finite same-length rewrite components terminate because there are only finitely many words at that length. The audit independently reconstructs every component and checks the complete raw matrix of every member.

If an equivalent word X has binary boundary blocks, write it chronologically as `L core R`, where L and R contain only O/Z. Then the exact index transformation is

    tau(u')=M_dagger(R) tau(u),
    tau(v')=M_L tau(v),
    Psi_X(u,v)=Psi_core(u',v').

Both transformed endpoints remain actual unit sources, with their new index equal to their actual U coordinate. The middle is strictly shorter and remains inside the original central menu. Iterating this shortening cannot continue indefinitely. Therefore discarding these components preserves the union of value sets. Keeping one exact representative from each remaining component also preserves that union. No graph-isomorphism identification or unproved global rewrite completeness is needed.

The independent component census recovers 235 retained words and 376 discarded components, with retained counts `1,1,1,1,3,8,21,55,144` by length. It additionally checks 1880 explicit directed-index absorption maps.

The reduced representatives are not presumed transpose-closed. The audited reduced search therefore allows both index orders. A diagnostic reduced search imposing u<=v happened to lose no target through 5001; this finite outcome does not itself justify that shortcut for every target. No finite counterexample to that shortcut is claimed here.

## 4. Graph-order estimate and length-eight boundary

The independent calculation of minima over all raw words, not only chosen matrix representatives, yields exactly `2,4,4,8,12,20,36,60,100`. Combining each minimum with the product bound and the exact graph order verifies all nine cases of `order<=floor(log2 n)+7` using integer powers of two, without floating logarithms. This is an upper bound for a represented target's compiler witness, not a lower bound on optimal tournament order.

For 1877 and 17383, the complete direct inverse domain and the independently solved progression inverse agree exactly. Their minimum central length is eight. Alternative words in the independent receipt are exact matrix-equivalent representatives; they do not indicate different arithmetic target sets. The named author witnesses JJJJJOJJ at (3,5) and (3,55) were separately replayed and counted literally. Their graphs have 14 and 17 vertices, the stated Hamilton-path counts, and exactly one Hamilton cycle. Thus 17383 genuinely defeats the proposed uniform length-seven family, without obstructing other physical source constructions.

## 5. Conditional coverage, not closure

The open statement quantifies over explicit integer indices and one finite central-word list. If it were proved for every odd n>2647, the physical compiler would supply C1 witnesses for all those targets, with the stated size bound. Below the threshold, the exact family inventory plus 27 literal repair graphs supplies the remaining positive odd existence claims other than 7 and 21; target one is separate and transitive.

This is a sufficient implication. The bounded-central family has not been proved equal to unrestricted C1 for all targets, and universal tournament coverage has not been proved equivalent to its occupancy lemma. The global nonrepresentability of 7 and 21 is also a separate inherited question/theorem, not established by this packet or by its finite C1 comparison.

## 6. Independent complete finite checks

Author normal and optimized Python runs at cap 5001 both match both frozen receipts exactly. The standalone `independent.py` imports no author code. It uses independent raw matrix generation, recursive binary source words, exact rational cone coefficients, explicit boundary-index absorption, increasing-height source buckets, an extended-Euclid progression inverse, and direct DFS enumeration of directed paths instead of subset dynamic programming.

Its normal and optimized receipts agree byte-for-byte. They independently recover:

- 592 exact generated header matrices, each with its correct minimum central length and product bound.
- 305,283 full-family arithmetic pair tests through 5001.
- 335,523 directed-index tests for the reduced menu, giving exactly the same values and minimum central lengths as the full family.
- 126,910 complete physical sources through U=5000, whose F image is exactly the same 2471 targets through 5001.
- The stated 29 gaps, all at most 2647. Their nonmembership in full C1 at these fixed values follows from the complete physical source comparison, not merely the central-length restriction.
- 180 progression membership controls, with full direct/progression candidate equality checked on the two critical targets and a separate small target range.
- All 27 repair adjacency lists, with their exact target Hamilton-path counts and two or three Hamilton cycles. Twenty-six have two cycles and the target 67 graph has three.

The direct source cutoff is complete because F>=U+1. Actual sources are reconstructed from the grammar; a relaxed cone is not the enumeration universe. Graph cycles are counted once by fixing their initial vertex. The critical 17-vertex graph's DFS examines 3,803,129 directed partial paths and verifies the final consequence directly.

## 7. Replayed C++ inventory and arithmetic safety

The supplied C++ source was independently compiled with g++ -O3 -std=c++17 and run at 10,000,001. Its output matches the frozen inventory receipt exactly: 1,663,896,007 tested index pairs and the same 29 length-eight gaps, with 31 gaps at length seven and 35 at length six. This is an independently replayed complete finite arithmetic inventory, not an all-order occupancy theorem.

Every header entry was reconstructed independently from the symbolic generators. The maximum absolute raw matrix entry is 64. In the documented cap range, g(u)<=u, each raw unit vector has coordinate sum at most 3u, and the product bound gives uv<=2N. Therefore the absolute sum of terms in any evaluated bilinear expression is at most `18*64*N`, at most 115,200,001,152 for N<=100,000,001. This is far below signed 64-bit range. The g-table recurrence stays within the 32-bit declared index cap, the index-product loop guards use 64-bit arithmetic, and array indices stay inside their allocations. Thus the claimed exact arithmetic has an explicit safety bound rather than an unchecked absence of observed overflow.

## Reproduction and limits

From this audit directory:

    python3 independent.py > check-normal.json
    python3 -O independent.py > check-optimized.json
    cmp check-normal.json check-optimized.json
    cmp check-normal.json INDEPENDENT_NORMAL.json
    sha256sum -c MANIFEST.sha256

The independent checker reads only the frozen generated header and repair adjacency lists. It locates them in the packet root when installed under its audit directory, or in the original neighboring source packet during this review; all source inputs are identified in `INPUTS.json`. The author's checker and C++ scan have the reproduction commands in the source README.

The compiler, finite inverse and normalization have been proved and checked. Their complete finite support is useful evidence. The universal occupancy assertion remains OPEN, and no extrapolation is included in this approval.
