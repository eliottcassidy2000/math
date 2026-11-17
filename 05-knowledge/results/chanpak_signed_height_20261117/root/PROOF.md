# Change the signed invariant while preserving the very same arc repair

**All-order construction.** The previous bank fixed 3h−R=1. That is not necessary: a power-of-two address can choose between infinitely many signed lines, and a word of physical migrations can move between the lines. The source is unchanged. This extends the geometry of the reader but does not solve general absolute placement.

The rooted state and migration are from [bank fusion](../chanpak_bank_fusion_20261117/migration/PROOF.md). Write (h,R,S) for Hamilton paths, marked terminal covers with an optional second path, and full paths starting at the root; there is one full path ending at the root. On S=0 the transfer of a width-k address u has upper (h,R) matrix

    M(k,u) = [[F(u), u+1], [D−rev_k(u), D]],     D=2^k,

where F(0)=1, F(2v+1)=F(v), and F(2v)=F(v)+v+[v odd] for positive even inputs. A padded leading zero keeps S=0.

## 1. The Fermat-like signed lines and their exact connecting graph

For q=2^i, i≥2, set ell_q=(q+1,−1). The seed is the transitive tournament on i+1 vertices rooted at its sink: (h,R,S)=(1,q,0), so ell_q(h,R)=1.

Consider the condition

    ell_q M(k,u) = ell_Q,       Q also a power of two ≥4.

The second coordinate gives (q+1)(u+1)=D−1. The multiplicative order of 2 modulo q+1 is exactly 2i: the powers below q are less than q+1, then 2^i≡−1, and multiplying that −1 by powers below q cannot give 1 before exponent 2i. Thus k=2it with t≥1, D=q^(2t), and u=(D−1)/(q+1)−1; this indeed fits with leading zero. Write v=u+1. The base-q digits of u, padded to 2t places, alternate 0 and q−1 except that the final q−1 is lowered to q−2. Reversing the binary word gives exactly

    rev_k(u) = qv−D/2,       c=D−rev_k(u)=((q+3)v+3)/2.

As t goes up by one, u_new=q²(u+1)+q−2. Peel the trailing zero, the preceding block of i−1 ones, then the next i zero bits in F. Since u/2 is odd, this gives

    ΔF = (q²+2q−3)v/2 + (q+1)/2.

For t=1, F(q−2)=q/2+1. Hence the first-coordinate defect E_t=(q+1)(F−1)−c is −q/2 at t=1 and increases by exactly 2 at each step (use v_new=q²v+q−1 in c). In particular

    ell_q M(k,u) = (q+1 + 2(t−1)−q/2, −1).

**Exact edge rule:** from an input on ell_Q=1 to an output on ell_q=1, a single positive-width S=0 migration of this form exists if and only if q≤2Q. It is unique among the leading-zero single steps and has t=1+Q/2−q/4. Thus one can climb one exponent at a time but jump arbitrarily far down. The stationary letter on line q uses t=1+q/4. Multiple such edges compose literally; changing line is not a relabeling of the count.

The height update on this edge is h_new=m h−v with m=F+(Q+1)v. The arithmetic of the register is unusually uniform: modulo q², v≡q−1, and

    2(q+1)m = [2(q+1)(Q+1)+q+3]v + 2Q+5 ≡ 2q mod q².

Thus the valuation of m is exactly i, irrespective of the incoming Q.

## 2. A different strong finisher for every line

Use the old strong one-arc source (174,124,120,240,65,80,1,118) with (H,P₂)=(187,3736). Its chosen reversal has (197,3728). The inherited uniform single-return reader observes H h_f+P₂ R_f; the change is 10h_f−8R_f.

For every even K≡2 mod4 let w=(4·2^K−6)/5. It has leading bit 1, so this last migration makes the incoming graph strong: the new root reaches the spine and then all old vertices, which have a path to the unique old terminal root and its new return arc. It satisfies

    (10,−8) M(K,w) = (K+4,−2).

The second coordinate is immediate from w. For the first, w at width K+4 is 16(w+1)+2. Peeling its four trailing bits in F gives ΔF=21w/2+13; reversing the same append changes c=2^K−rev by 21·2^K/2. Thus the first coordinate increases by 4; its base at K=2 is 6. Taking K=2q−2 yields (2(q+1),−2)=2ell_q. Consequently any path along the signed-line graph, followed by the appropriate strong finisher, makes the original arc reversal add exactly 2.

Seed-and-finish counts at i=2,3,4 are 1169603, 553857347, and 69674591564099, with +2 partners. Both the number of signed lines and the number of paths on them are unbounded. The large widths are a real cost, not an asymptotic improvement disguised by the address notation.

## 3. Decode a prescribed height or certify a miss in this grammar

For a prescribed lower endpoint N, the final line has K=2q−2. Since R≥q and w+1≥2^(K−1), positivity gives N≥2^(K−1). Enumerate powers of two q under that bound. For each q the reader count is an affine strictly increasing function of h, so solve for a unique candidate positive integral h.

Now work backwards on the switching graph. h=1 is a seed. If an edge from Q led to the present (q,h), its multiplier and deposit give h≥m−v≥Qv. Also v≥q^(2(t−1)), so i(Q−q/2)≤floor(log₂ h). This gives a finite bound on Q, in addition to Q≥q/2. For each candidate power of two Q, require exact integral positive division (h+v)/m; it is strictly smaller than h. Recurse with memoization. This certifies exact membership or nonmembership within the described physical grammar; it does not guarantee a unique predecessor, uniform time complexity, or every target. A single stationary ray has the simpler exact test: for recurrence h_new=m h−v, the shifted register (m−1)h−v must be its initial value times a power of m.

## 4. The phase obstruction to changing *only* the invariant

For this fixed source, every integral two-coordinate reader with exact signed difference +2 necessarily obeys 5h_f−4R_f=1. Thus h_f=4s+1, R_f=5s+1, and the lower endpoint is identically

    N=3923+19428s.

This is independent of our migration assumptions. Several heights in the new grammar lie outside the previous narrower phase 42779 mod 116568, but no change of signed line alone can escape 3923 mod 19428 while retaining this source and a two-coordinate +2 read. The right conclusion is not “grow the line graph further”; it is “change the source, change the switch vector, or let a reader use a different observable.”

## 5. Even the complete changing-line grammar is sparse

Allow every starting q and every legal path, with no uniqueness assumption. At an edge put G=m−v. For h≥1 the new height is at least Gh; moreover G≥Qv. For a climb q=2Q, v=2Q−1. For every other edge q≤Q, t−1≥Q/4, hence v≥q^(2(t−1))≥2^Q. These estimates give a uniform contraction
     
    sum_(q≤2Q) G(Q,q)^(−1/6) < 31/32

where q runs through powers of 2 at least 4. At Q=4 the exact growth floors are 33 and 241, with reciprocals of sixth roots less than 9/16 and 13/32 respectively. At Q=8 the climb floor is at least 120, and the two other floors at least 104856 and 233016, giving a sum below 1/2+1/6+1/6. For Q=2^j≥16 the climb contributes less than 1/2 and the remaining j−1 choices contribute at most (j−1)/(Q2^Q)^(1/6) < 1/3 (the bound decreases with j and holds at j=4).

Starting at arbitrarily large Q does not break summability. For nonempty paths the sum over all first edges is finite: the climb contributions are O(Q^(−1/3)) over dyadic Q, and the nonclimb totals are O((log Q)2^(−Q/6)). Each subsequent layer is contracted by the uniform factor. Empty paths are summable after finishing since N≥2^(2Q−3). All readers have N≥h, so the sum of N^(−1/6) over all finite histories is finite. Hence there are **O(X^(1/6))** realized lower endpoints below X in this grammar, including all starting lines and all detours. This is an upper bound on this particular repair bank only, not on any reader satisfying 5h_f−4R_f=1.

## Sequence archaeology that actually informed the construction

* The repo's [Stern–Brocot/Berggren ordinal](../../01-canon/theorems/THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge.md) separates the primitive rational address from its scale. Here the signed line is the address, while h is the scale. Undo the edge on the address graph before solving an integer division.
* The [Zeckendorf guard register](../zeckendorf_guard_automaton_20261003.md) explicitly separates the integer and its shifted register. Here the backward test keeps the physical h and a separate bit-reversal read. Dropping the reversal would lose the linear defect 2(t−1)−q/2.
* The [toothpick four-comb](../../01-canon/theorems/THM-1127-four-comb-toothpick-self-similarity.md) makes varying scales sample one fixed object instead of treating each tail as independent. Here the same digit-block update governs all heights and produces the whole directed graph q≤2Q, rather than unrelated hand-picked widths.
* The [aliquot affine compiler](../aliquot_dyadic_frontier_20261014e.md) separates the reversed affine candidate from integrality and native admissibility. Here division, positivity, and the power-of-two line condition are all actual decoder guards. The [taxicab boundary work](../golden_taxicab_boundary_20261011.md) provides the caution: choosing which boundary coordinate is exposed does not make every algebraic solution a physical realization.

Run `python replay.py`: exact drift, bit reversal and reader-row controls; literal cycle-free seed and physical rooted DP; independent target inverses on 63 stationary samples and 115 arbitrary changing-line words. The formulas are proved for unbounded width above, not extrapolated from those tests.
