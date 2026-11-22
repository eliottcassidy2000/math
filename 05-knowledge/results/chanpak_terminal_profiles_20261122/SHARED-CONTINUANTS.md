# Shared tournament parts give exact compatible continuants

This extends the [transitive-gap backchain](../chanpak_absolute_descent_20261122/BACKCHAIN.md) to arbitrary tournament parts and to all ordered path-cover counts. It proves new exact scalar heights, an upper bound for the order of any compatible presentation, and exact separation of two different auxiliary counts. It does not prove that every sufficiently large odd height has a compatible state.

## 1. Join arbitrary parts along shared boundary vertices

Take ordered spine vertices \(b_0,\ldots,b_k\). Between \(b_{i-1}\) and \(b_i\), insert an arbitrary tournament \(C_i\), possibly empty. Direct all edges between lattice positions forward except the \(k\) fixed arcs \(b_i\to b_{i-1}\), and retain each \(C_i\)'s internal orientation. A boundary-reaching directed cycle must use all consecutive reversed arcs between two extrema \(b_a,b_b\), and then one forward path from \(b_a\) to \(b_b\) through exactly the gap blocks in that interval, allowing empty visitation in any block. Indeed each intervening cut has just one backwards arc, so the cycle uses it once. After the reversal from \(b_b\) through \(b_a\), all the intervening boundary vertices have already been used; every remaining block transition is forward. A cycle meeting no boundary is fully inside one gap.

Let \(Z_T(q)\) count ordered covers by \(q\) paths, with empty slots allowed. Thus \(Z_T(1)=H(T)\), and \(Z_T(0)\) is 1 for empty \(T\), else zero. From the [positive odd-cycle formula](https://arxiv.org/abs/2307.05569), \(Z_T(q)\) can be counted by choosing a disjoint packing of odd cycles, placing each odd cycle in one slot with weight \(2q\) and each uncovered vertex in one slot with weight \(q\). This follows by expanding the odd-cycle theorem separately in each of the \(q\) ordered slots; once the packing is fixed, all its cycles and all remaining singleton vertices are independent color atoms.

For any gap \(C\), the two reads for the portion selected by one spine cycle are exactly
\[
 \sum_{S\subseteq C} H(C[S]) Z_{C-S}(q)=Z_C(q+1),\qquad
 \sum_{S\subseteq C}(-1)^{|S|} H(C[S]) Z_{C-S}(q)=Z_C(q-1).
\]
For the unsigned equality add a distinguished slot. For the signed equality the same distinguished slot contributes \(-1\) per atom placed there, since every cycle and every uncovered singleton has odd order. Its \(q+1\) color choices thus give signed total \(q-1\) per atom. This also proves the equalities at zero-size gaps and the cancellation at \(q=1\) for a nonempty gap.

Write \(a_i=Z_{C_i}(q+1),\,b_i=Z_{C_i}(q-1),\,c_i=Z_{C_i}(q)\). Uncrossed gaps contribute \(c_i\). A maximal run of crossed gaps from \(a+1\) through \(b\) corresponds to one cycle at the level of spine intervals; summing its internal choices and parity has weight
\[
 q\left[\prod_{i=a+1}^{b}a_i+\prod_{i=a+1}^{b}(-b_i)\right].
\]
Indeed the cycle contains \(b-a+1\) spine vertices; its slot has \(q\) choices, and its factor 2 cancels the usual half in the odd parity projector. A single empty gap gives weight zero, consistently with the lack of a forward return across a single reversed boundary edge. A singleton spine vertex outside the cycles instead contributes \(q\).

## 2. Positive three-register formula for every \(q\)

Let \(F\) be the current total on a prefix of gaps, \(G\) the part where its last boundary is a singleton, and \(A,B\) the terminal interval contributions using the positive and negative signed products, respectively. Thus \(F=G+A+B\), initialized at \((G,A,B)=(q,0,0)\). For a new gap with data \((a,b,c)\) the exact update is
\[
 G'=qcF,\qquad A'=a(A+G),\qquad B'=-b(B+G).
\]
To extend a terminal interval, multiply by that same mode; to start one, use the previous last boundary as a singleton (already carrying one \(q\)). Leaving the new boundary singleton instead introduces \(qc\). Since distinct intervals cannot share a boundary, this gives all choices exactly once.

There is a nonnegative version: put \(D=A+B,\ U=-B,\ V=G+B\), and \(F=D+U+V\). Start at \((0,0,q)\) and use
\[
 \begin{pmatrix}D'\\U'\\V'\end{pmatrix}=
 \begin{pmatrix}a&2a&a-b\\0&0&b\\qc&qc&qc-b\end{pmatrix}
 \begin{pmatrix}D\\U\\V\end{pmatrix}.
\]
All coefficients are nonnegative since \(a\ge b\) and \(qc\ge b\). In particular this is one literal construction producing both \(Z(1)\) and \(Z(2)\) from the same ordered modules, not an arbitrary linear machine assigned values after the fact.

For \(q\ge2\) and a specified gap, the three-register map is invertible over the rationals, with determinant \(qabc>0\). In the \((G,A,B)\) coordinates recover \(x=G'/(qc),y=A'/a,z=-B'/b\), then \(G=y+z-x, A=x-z, B=x-y\). Integrality and membership of the recovered state remain necessary. For \(q=1\), any nonempty gap has \(b=0\) and loses that register; this is a precise difference between one-slot and two-slot compatibility.

## 3. Scalar continuants and exact partial inverse

At \(q=1\), write \(h_i=H(C_i),Q_i=Z_{C_i}(2)\) for nonempty gaps. They reset \(B=0\), so retain the pair \((F,G)\) at those boundaries. The successive positive step, directly or after one empty gap, is exactly
\[
 (F',G')=((h+Q)F,hF),\qquad
 (F',G')=(2(h+Q)F-hG,\ h(2F-G)),
\]
respectively. For the second formula the empty gap sends the old \((G,A,B)\) to \((F,F,-G)\). Ignoring the preceding second count would change the exact target.

Put \(\epsilon_i=1\) for a direct step, 2 for a step with one intervening empty gap. Then \(F_i-G_i=\epsilon_iQ_iF_{i-1}\), and eliminating \(G\) gives the nonnegative nearest-neighbor continuant
\[
 F_i=(h_i+\epsilon_i Q_i)F_{i-1}
    +(\epsilon_i-1)h_i\epsilon_{i-1}Q_{i-1}F_{i-2}.
\]
Start with \(F_0=G_0=1,F_{-1}=0\). When every gap is nonempty and direct, this is \(\prod_i(h_i+Q_i)\) on \(\sum_i|C_i|+k+1\) vertices, saving \(k-1\) vertices relative to gluing the corresponding cyclic wrappers by ordinal product. The full host is strong with all gaps nonempty, and also with isolated empty gaps strictly between nonempty ones: a zero gap is crossed forwards by jumping two boundary positions then returning on the adjacent reversed arc. All boundaries can move backwards along the spine, and every vertex inside a gap is bracketed by its boundaries.

In an isolated-zero word, after the first positive boundary \(0<G/F\le1/3\), since \(Q\ge2h\) and the bridge ratio is at most the direct ratio. For a gap with \(t=(h+Q)/h\), its last-step ratio \(G'/F'\) is \(1/t\) if direct and in
\[
 [5/(6t-1),\,1/t)
\]
if bridged from a previous non-initial positive boundary. This follows from the strictly decreasing map \(r\mapsto(2-r)/(2t-r)\) on \(r=G/F\in(0,1/3]\). The first step can always be chosen direct in this normalized language. In particular all transitive choices \(t=1+2^m\) and the nontransitive four-vertex count-five gap with \(t=37/5\) have pairwise disjoint intervals. For the inserted \(t=37/5\), \(1/9<25/217\) and \(5/37<5/29\) separate the neighboring transitive choices. A proposed paired state hence determines its last module type and whether that application was direct or bridged in this enlarged alphabet. A bridge recovers \(F=(F'-G')/(2Q)\), \(G=2F-G'/h\) exactly. A direct step recovers \(F=G'/h\) but forgets \(G\), leaving a genuine scalar membership problem below it.

**Order bound even with arbitrary zero runs.** Each nonempty gap multiplies \(H\) by at least \(h+Q\ge3\) and at least \(2^{|C|}\), since each assignment to the two slots admits paths in both induced tournaments. In nonnegative coordinates every two consecutive empty gaps take \(F=D+U+V\) to \(3D+5U+3V\ge3F\). Pair disjoint zero steps within each run, leaving at most one unpaired per run. The number of positive steps plus zero pairs is at most \(\log_3 H\), hence the number of gaps is at most \(2\log_3 H+1\), including no more than one extra unpaired zero. Thus any successful presentation in the entire generalized grammar has order at most \(2+\log_2 H+2\log_3 H\). This is conditional on compatibility, not an existence theorem for arbitrary heights.

## 4. Verified new heights and actual separation of second counts

A single gap with rows \((2,8,3,5)\) has \((h,Q)=(5,32)\), giving a strong six-vertex height \(37\) and pair \((37,5)\). Preceding this gap by a positive transitive width \(m\) and exactly one empty gap gives
\[
 (F',G')=\bigl(74\cdot2^m+69,\,10\cdot2^m+5\bigr)
\]
on \(m+8\) vertices. In particular \(217,365,661,1253,2437,19013,37957,75845\) from \(m=1,2,3,4,5,8,9,10\) (and \(37\) itself) are absent from the complete *old transitive-gap grammar* below \(100001\), verified by the exact original positive transfer. This is a strict architectural enlargement, not an asymptotic improvement to the proved old sparsity exponent or a claim of new unrestricted first births.

The modules \((2,4,65,23,39,15,59)\) and \((2,8,11,97,15,23,55)\), both strong, have \((h,Q)=(45,732)\) and \((37,740)\). They give the same height \(777\) but different compatible second registers \(45,37\); appending the same isolated-zero and positive transitive block separates their new heights by eight. For any \(|C|\ge2\), \(Z_C(2)\equiv0\bmod4\): the two empty allocations contribute \(2H\equiv2\bmod4\), and the remaining unordered bipartitions number \(2^{|C|-1}-1\), each with twice an odd product. Thus equal \(h+Q\) forces the linked second-register difference to be divisible by four (the one-vertex case has no second homometric module); this particular exchange cannot cause a two-unit gap through the stated coupling.

The three-register two-slot observer has strictly more information than even the compatible scalar pair. The following two strong eight-vertex gap modules have the same \((h,Q)=(345,5640)\) but different \(Z(3)\):
\[
 (80,13,25,97,74,23,166,63):45225,\qquad
 (128,129,27,51,99,71,15,124):44505.
\]
They are interchangeable in every scalar/compatible-prefix \(q=1\) chain. At \(q=2\) their one-gap wrappers already differ by \(1440\), since \(Z(2)=2[Z_C(3)+2Z_C(2)-H(C)]\). More generally, at the changed gap \(a\) is the only differing coefficient; its nonzero change affects \(D'\) by a strictly nonzero multiple \(D+2U+V\). Later nonnegative matrices with positive diagonal in \(D\) preserve the strict difference. Thus the genuine ordered two-slot count is not determined by a chain's scalar pair, and the adjacent levels of the slot hierarchy are actually distinguished by physical modules.

## 5. A prime-height restriction on compatible chains

A common divisor of the scalar registers cannot disappear under any future integer transfer. At a nonempty-gap boundary \(B=0\) and \((G,A,B)=(G,F-G,0)\), so this common divisor is exactly \(\gcd(F,G)\). An adjacent nonempty gap (no intervening zero) multiplies both new registers by the old \(F\). If the old \(F>1\), the new height is already strictly larger than this divisor, and heights never decrease thereafter. Thus a prime final height cannot use such an adjacent positive step. For an arbitrary intervening zero run, any new nonempty part with \(\gcd(h,Q)>1\) equally introduces a persistent divisor, since \(F'=Q(A+G)+hF,\ G'=hF\) with the old full state before the nonempty part. Thus prime-height presentations in this architecture require coprime gap profiles and at least one zero between each two nonempty gaps. This removes reset steps from the search for prime targets; it does not imply that all prime targets have such a search path.

For precisely one intervening zero and a primitive incoming pair, a gap whose \(Q\) is a power of two has the even sharper identity
\[
 \gcd(F',G')=\gcd(2QF,h(2F-G))=\gcd(F,h).
\]
The first equality subtracts the two outputs, and the second uses oddness and \(\gcd(F,2F-G)=1\). Hence transitive bridges preserve primitiveness exactly, whereas a nontransitive dyadic-slot bridge preserves it precisely when its one-slot count is coprime to the current height. For example, repeating the count-five gap \((h,Q)=(5,32)\) with isolated zeros after its first use never creates a common divisor, since its first height is 37, and \(F'\equiv4F\bmod5\) at each such step. This is a literal family of primitive continuant states, not a primality assertion.
