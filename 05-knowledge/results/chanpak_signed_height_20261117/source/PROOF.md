# Changing the one-arc source, changing the reader, and exact source filters

Write (H,P₂) for the ordinary and optional two-slot Hamilton counts. The single-return reader with rooted finished values (a,b) reads aH+bP₂. The arc in its source is the only differing edge; all reader contacts remain fixed. This note gives two new source phases on the previous reader, a wholly different reader height with a physical source, and bounds excluding natural attempts to make the source indefinitely taller.

## Two new phases for the branching reader

Start from the previously verified order-eight paired source with outgoing masks

    (174,124,120,240,65,80,1,118),

flipping {5,6} for its mate. Append two vertices independently with fixed outmasks to the old eight vertices, making the first new vertex beat the second. Among all 65,536 ordered contact-mask pairs, exactly two retain the oriented difference (Delta H,Delta P₂)=(10,−8). The lower source profiles are

| appended masks | H,P₂ | resulting lower count at reader coordinate h |
|---|---|---|
| (12,31) | (7999,157792) | 2462868h−655165 |
| (43,75) | (8021,157952) | 2465532h−655871 |

The inherited reader has (a,b)=(12h−3,15h−4), so these really give +2 at every legal h of its branching grammar. The new initial pairs (h=1) are 1807703/1807705 and 1809661/1809663. The same flipped old arc is the only difference at every reader height. For exact membership in these phases, solve the indicated linear equation for positive odd h, trying N and N−2 for either track, and apply the original ABC inverse to h. Merely checking its residue is not acceptance.

## An alternate m=6 stationary reader with a physical one-arc source

Use the line 63h−R=−1, initial transitive TT_7 rooted at its last vertex, (h,R,S)=(1,64,0). Padded migration of width 18 and address 4160 preserves S=0 and sends h→266272h+4161 on the line. This follows directly from F(4160)=4129 and rev_18(4160)=2080. After any number of repeats finish with width 2, address 2. The rooted read has (a,b)=(192h+3,255h+4), and the final rooted graph is strong because the final root reaches the first vertex of its two-vertex spine and every old vertex reaches the old ending port. The old mark has one full end by the root migration invariant.

Take lower order-ten source

    (102,260,88,851,131,286,50,367,597,247).

It has (H,P₂)=(5781,119040). Flipping {4,5} gives (5611,119056). Append three ordinary ordinal sinks to each; these multiply P₂ by 8 and leave H fixed. The resulting order-thirteen oriented difference is (−170,128), so the physical count changes by

    −170(192h+3)+128(255h+4)=2.

The lower count is 243951552h+3826623, beginning at 247778175, with a +2 companion. Exact inverse: compute h from the initial affine read, then repeatedly invert h→266272h+4161 by positive exact integer division until 1; failure of any step rejects. Heights from this ray are different from the original three-letter grammar. This is a changed profile and reader, not arbitrary source-order coverage.

More generally for even m the proposed direction (−(8·2^m−2)/3,2^(m+1)) requires source-order parity. The odd-cycle-pack law gives P₂≡(−1)^n H mod3, since a pack with d atoms has q-weight q^d times its fixed cycle multiplicity, and d≡n mod2 because each atom has odd size. Thus m≡4 mod 6 is completely impossible for **any** pair of the same order, m≡2 mod 6 needs even source order, and m≡0 mod 6 needs odd order. In particular m=4's (−42,32) cannot work even without the one-arc requirement.

## Structural source constraints and exact inverse for a fixed reader

For canonical (10,−8) the source order is even by the same parity identity. The changed arc stays inside a common SCC in the two source versions: a reversal between SCCs preserves every old internal odd-cycle pack and only adds new packs, so H and P₂ cannot move in opposite directions. Apply this also in reverse if an SCC would split. Every unaffected ordinal component therefore factors out of both counts. Its P₂ factor divides 8; the internal difference is divisible by 4 because every P₂ on at least two vertices is divisible by 4 (pair complementary subsets, whose counts are odd by Rédei), so the unaffected product P₂ is at most 2. As P₂≥2^n, there is at most one outside vertex. Thus a canonical unit source is strong in both versions, or differs from such a core by just one ordinal singleton.

For a fixed original finished reader with 5a−4b=1, every possible +2 source vector is

    (Delta H,Delta P₂)=(10+4t b, −8−4t a),   t integral.

The Bezout general solution reduces to multiples of 4 because Hamilton differences are even, optional-two-cover differences divisible by four, and a is odd. If t≠0, positivity of the second source profile implies the lower N≥4ab−9a (use H≥4b−9 for t negative or P₂≥4a+8 for t positive). Hence in the high-reader regime 4ab−9a>N only the canonical direction is possible. To decode a physical source profile at any fixed height, H≡5N mod b and P₂=(N−aH)/b. A supplied source-order bound n gives H≤N/(a+(n+1)b) from the injective full-path cuts into two labeled slots. If that interval is shorter than b it has at most one candidate; checking the integer profile or its feasible contact graph cannot be replaced by residue selection.

## No uniform one-port growth of the original source

For a tournament B replacing any nonswitched original vertex j uniformly, count path covers by the number r of nonempty B runs. For one or two labeled slots the skeleton has at most respectively 1+min(indegree,outdegree) or 2+min(indegree,outdegree) such runs: each additional run within the same slot consumes a distinct outgoing and incoming exterior vertex. Thus each read is a finite linear combination of Q_r(B), the ordered nonempty r-slot counts. These coefficients are recovered exactly by substituting TT_t for t=1,...,the corresponding bound: Q_r(TT_t)=r! S(t,r), triangular. The Delta H rows for j=0,1,2,3,4,7 are

    (10,−207,−170), (10,4), (10,4), (10,6), (10,118,54), (10,6).

For |B|≥2, Q₂(B)≥H(B) by splitting a full path at an internal cut. Thus the j=0 delta is at most −197H(B); every other delta exceeds 10. None is even ±10. In particular there is no hidden nontransitive substitute at this single port with canonical gap.

More generally, a fixed source quotient with one growing uniform TT_t and the switched arc elsewhere has H/P₂ differences that are finite exponential polynomials Σ c_j j^t, since the number of visits is bounded and r! S(t,r) has this form. If both differences are (10,−8) infinitely often, even on one parity of t, the leading nonzero exponential would dominate, so both functions must be constant at all positive t. But the order parity alternates, contradicting the necessary even order. This excludes arbitrary *independent* growth of one transitive module in a fixed source architecture.

Finally suppose two arbitrary uniform TT_t blocks are tethered at the same height, with all arcs between them in one direction, and original source unchanged. Across each projected path slot, between visits to fixed old vertices, their quotient order is monotone; all visit counts are bounded by the number of old vertices and slots. Thus the differences are finite exponential polynomials in t again (bases ij for the two modules). A constant infinite tail would force canonical (10,−8) for t=1 and t=2. All t=1 candidates are the two explicit order-ten masks above, but t=2 gives respectively (210,7104) and (−2710,−40216). Therefore this entire equal-width tandem enlargement fails, even though its parity is compatible. This does not rule out higher-order nonuniform stages or different linked growth directions.

Replay `python verify.py` checks literal source profiles, triangular coefficients, tandem controls, and stationary invariance. Compile `find.cpp` for the exhaustive finite two-new-vertex phase test.
