# Prime five-name quotients with arbitrary even children are visible modulo sixteen

## 1. Statement

Let Q be any prime tournament on five quotient names, and let A_0,...,A_4 be arbitrary nonempty named tournaments, each of even order at least four. Compare Q[A_i] with Q^op[A_i], keeping all physical child names and interiors fixed. A block-uniform physical two-port query distinguishes their complete Hamilton responses modulo sixteen. Profile-preserving target interior replacements are also allowed.

At most two fixed block queries suffice once an isomorphism to a representative below is fixed. The representative labeling is an arithmetic coordinate system: every chosen quotient bit is expanded onto the original physical child set. No physical names are identified, removed or replaced in the query.

This packet proves the LOCAL quotient-converse theorem. It is not unrestricted reconstruction or pairwise orbit injectivity. An ambient extension is treated separately; it is not inferred merely from a nonzero modulo-sixteen gap.

## 2. Complete actual profile domain and its one-bit reduction

For every even-order child of order at least four, H=q1 is odd. The profile-domain theorem in `../even_cyclic_mod16/PROOF.md` gives the complete normalized upper domain

    q/H=(0,1,b,c,8,d,0,...),
    b in {6,14}, c in {4,12}, d in {0,8} modulo16.

For completeness, q2/H and q3/H follow from their modulo-eight formulas. Dividing ordered covers by k! and reducing modulo two gives S(N,k), because every induced Hamilton count is odd. The Stirling recurrence gives S(N,4) odd for even N>=4, so q4=8 modulo16; q5 is zero or eight. Every k>=6 has sixteen dividing k!. Multiplication by an odd inverse retains zero/eight. The domain need not be entirely realizable; it contains every actual profile.

The actual TT4 profile reduces to e=(0,1,14,4,8,0). Write epsilon_i=0 for b_i=14 and epsilon_i=1 for b_i=6. Expand the exact child-run functional around e. Each perturbation has coefficient eight, so products of two perturbations vanish modulo sixteen. In a first-order term, every other child must appear exactly once because only e1 is odd.

If the perturbed child occurs three or five times, exchange its first two intermediate return segments. They are nonempty, have disjoint label supports, and have identical quotient-child endpoints. Every other child and each port occurs only once. The exchange preserves all directed joins and is fixed-point-free, proving an even coefficient. These q3/q5 perturbations vanish after multiplication by eight. This is the inherited return-segment involution from `../parity_profile_transfer/PROOF.md`, applied as in the contemporaneous `../even_profile_linearization/PROOF.md`; the argument is reproduced here and is Hamilton-specific, not an ordered-cover claim.

Let D_i(Q;a,c) count run words with child i used twice, every other child once, both ports once, and all consecutive symbols following the quotient/port arcs. Then the whole block-uniform response satisfies

    F_(Q[A_i])(a,c)
      = Lambda [F_(Q[TT4])(a,c)+8 sum_i epsilon_i D_i(Q;a,c)] mod16,
    Lambda=product_i H(A_i), an odd unit.                 (1)

The same Lambda applies to the opposite quotient. No value of H modulo sixteen is lost. Define C_i=D_i(Q)-D_i(Q^op) modulo two and B=F_(Q[TT4])-F_(Q^op[TT4]) modulo sixteen. The complete contrast is therefore

    Delta_actual=Lambda[B+8 sum_i epsilon_i C_i] mod16.   (2)

Both full responses, not just their contrasts, are retained by this formula.

## 3. All prime quotient types and explicit certificates

Encode a labeled five-name tournament by the ten bits in increasing unordered-pair order, bit one meaning the smaller name points to the larger. Exhaustive module testing on all 1,024 codes gives 264 prime labeled tournaments. Explicit permutation orbits of the following representatives are disjoint and cover them, with sizes120,120,24:

    Q10=(20,1,2,7,14),
    Q41=(18,8,3,5,14),
    Q76=(24,17,3,6,12).

Rows are outgoing bit masks. Thus no prime quotient type is omitted. This finite completeness check is reproduced directly and is not an assumed census.

The complete exact run-word bank gives these affine contrast certificates before multiplication by Lambda:

- Q10 at (a,c)=(0,2): 12+8(epsilon0+epsilon2+epsilon4) modulo16.
- Q41 at (a,c)=(0,1): 6+8(epsilon0+epsilon1+epsilon2) modulo16.
- Q76 at (a,c)=(0,1): 8(epsilon2+epsilon3) modulo16.
- Q76 at (a,c)=(2,18): 8+8(epsilon2+epsilon3) modulo16.

The first two never vanish. In the last pair, exactly one value is eight. Multiplication by the common odd unit preserves nonvanishing, and preserves the value eight. This proves the theorem for every actual even-child profile assignment. It also proves visibility for all 8^5 formal upper-domain assignments, without asserting they are all realized.

## 4. Complete banks and executable proof

The run-word dynamic program includes both port orders, all empty prefix/suffix possibilities, and every child multiplicity contributing modulo sixteen. A run symbol cannot follow itself without a port or another child between them. TT4 run weights have valuations 0,1,2,3 at multiplicities1,2,3,4; multiplicity5 vanishes in the baseline. Thus only total valuation at most three is retained, exactly rather than approximately. The first-order profiles are separately checked at multiplicities2,3,5 with every other child used once.

The receipts retain all 1,024 block queries for each quotient and each orientation: baseline values modulo sixteen and all fifteen first-order coefficient parities, including the identically zero q3/q5 columns. They also retain the complete 256-coordinate Boolean coefficient bank of bidegree at most(2,2), for each orientation and every affine channel. Every full query is independently reconstructed from that bank. All cubic and quartic coarse interactions and zeros are kept; nothing is replaced by an endpoint-only or quadratic approximation.

All 32 remaining q2-bit states are checked against the full bank and the explicit battery. Their visibility covers all8^5 profiles by Section2. Eight independent literal physical counts on22 names verify the displayed TT4 baseline queries, using a predecessor-removal subset dynamic program rather than quotient runs. These counts are exact residues modulo sixteen, not claimed full integers.

`INTERPOLATION.md` describes a second small-physical route for independent replay of every baseline and correction channel. The retained complete banks are finite exact certificates, not a claim of exhaustive reconstruction of all tournaments of order twenty or higher.
