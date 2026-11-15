# Composition defect and the missing third slot

## What was literal in the old tower work

The repo's [Paley/Fano construction](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_mod6_20260922_paley_fano_octonion_design.md) recovers octonion signs from the order-seven Paley tournament **together with** one of its two Steiner triple systems. The tournament or the number seven alone does not specify the multiplication. The older [two-towers reflection](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/07-reflections/two-towers-cayley-dickson-and-hyperoperation-meet-at-apex7-but-the-octonion-structure-is-only-numerical.md) records the contrasting failed experiment: the Lonely Runner three-way emptiness cumulant did not single out the Fano-line triples. That is a negative experiment, not a general nonexistence theorem for all connections.

Here is the operative Cayley–Dickson phenomenon. For real, complex, quaternion, or octonion source, use the standard convention
\[
(a,b)(c,d)=(ac-\bar d b,\quad da+b\bar c),\qquad N(a,b)=N(a)+N(b).
\]
The four squared terms in the product norm give exactly the product of the two input norms. The only issue is the mixed-term difference
\[
N((a,b)(c,d))-N(a,b)N(c,d)
=2\Re((da)(c\bar b))-2\Re((ac)(\bar b d)).
\]
Associativity and cyclic invariance of the real part make these terms agree through quaternion sources, giving composition through octonion outputs. In an octonion source, the real part is insensitive to rebracketing three factors, and cyclic on three factors, so the displayed defect is exactly
\[
-2\Re(d[a,c,\bar b]),\qquad [a,c,e]=(ac)e-a(ce).
\]
Thus an associator invisible at the scalar level with three inputs becomes observable when paired with a fourth. For a concrete check at the next rung, with lexicographic doubling basis \(e_0=1\), \(e_i=(e_i,0)\) and \(e_{m+i}=(0,e_i)\) at doubling size \(m\),
\[
(e_1+e_{10})(e_4-e_{15})=0;
\]
the four basis products in reading order are \(e_5,-e_{14},-e_{14},e_5\), while the two input norms are both 2. The norm identity is lost despite the literal doubling formula continuing indefinitely. See [Baez, The Octonions, §2.2](https://math.ucr.edu/home/baez/octonions/node5.html) for the classical tower.

## Exact tower in the tournament reader: from two coordinates to three

Write \(q_k(A)\) for the number of ordered covers of a nonempty tournament \(A\) by \(k\) directed path slots, allowing empty slots. Let \(P_i(A)\) count the ordered covers by \(i\) **nonempty** paths: \(q_k(A)=\sum_i\binom{k}{i}P_i(A)\), \(q_0=0\), and \(q_1=H\). Let \(C(A,B)\) be the tournament on \(A\to B\to c\to A\) with \(c\) a singleton. Then for all \(k\ge1\),
\[
\frac{q_k(C(A,B))}{k}=
(q_{k+1}-q_k)(A)(q_{k+1}-q_k)(B)+q_k(A)q_k(B)
+(q_k-q_{k-1})(A)(q_k-q_{k-1})(B).
\]
**Proof.** Mark the slot containing \(c\), in \(k\) ways. The segment immediately before \(c\) may be empty, B-only, or A-then-B, encoded \(1+v+uv\). The segment after \(c\) may be empty, A-only, or A-then-B, encoded \(1+u+uv\). Each other slot has all four choices, as its cross-edge is always \(A\to B\). Thus the result divided by \(k\) is the evaluation of
\((1+v+uv)(1+u+uv)((1+u)(1+v))^{k-1}\) by \(u^i v^j\mapsto P_i(A)P_j(B)\). Set \(s=1+u,t=1+v\). Its polynomial equals
\[
(s^{k+1}-s^k)(t^{k+1}-t^k)+s^kt^k+(s^k-s^{k-1})(t^k-t^{k-1}),
\]
as expansion checks. The weights now evaluate to the stated optional-slot counts. This is an all-order bijective derivation.

This is exactly a **rank-three positive Gram form**, for every \(k\ge2\). It cannot be expressed as a bilinear form on two scalar features of each block: test \(TT_1,TT_2,TT_3\), for which \(q_j(TT_n)=j^n\). At fixed \(k\ge2\) the three geometric bases \(k-1,k,k+1\) are positive and distinct; the three columns at \(n=1,2,3\) are independent by the Vandermonde determinant. The displayed differences are an invertible change of coordinates. At \(k=1\), \(q_0=0\) makes two coordinates equal, collapsing the rank to two and recovering the earlier \(x_Ax_B+2y_Ay_B\) with \(x=q_2-q_1,y=q_1\).

### A mod-eight pattern survives the rank change

Every induced subtournament, including the empty one with convention \(H=1\), has odd Hamilton path count. Hence \(q_j(A)\equiv j^{|A|}\pmod2\) by writing the count as a sum of products over \(j^{|A|}\) ordered assignments of vertices to slots. Also \(q_2(A)\equiv0\pmod4\) for \(|A|\ge2\): complementary subsets pair equal odd products, with an even number of pairs. More generally **every even \(j\) has \(q_j(A)\equiv0\pmod4\)** for \(|A|\ge2\). In the falling factorial expansion terms from \(P_i(A)/i!\) for \(i\ge3\) are divisible by 4 when \(j\) is even (the factors \(j,j-2\) occur). The first two terms are
\(j(H+(j-1)P_2/2)\), divisible by 4 since \(P_2/2\equiv H\pmod2\).

Consequently on the diagonal \(A=B\), the normalized count \(q_k(C(A,A))/k\) is **3 mod 8 for odd k, 2 mod 8 for even k** when \(|A|\ge2\). For odd \(k\) all three Gram coordinates are odd; for even \(k\) the middle one is divisible by 4 and the other two are odd. For singleton \(A\), the normalized count is \(k^2+2\) instead (giving 6 mod 8 at \(k\equiv2\pmod4\)). The same odd-grade mod-eight phase survives although the actual positive form has changed rank. One must not infer that rank two or the old exact prime-norm characterization survives.

## What a second wrap can read that the first cannot

For \(R(A)=C(A,TT_1)\), the Gram formula specializes to the inherited singleton transfer
\[
q_k(R(A))=k(q_{k+1}(A)+kq_k(A)-q_{k-1}(A)).
\]
See [the prior reader-retention session](https://github.com/eliottcassidy2000/math/blob/4ff6bb2a9e6def238d5b0d87264b5f44a23fa33b/05-knowledge/results/chanpak_reader_retention_20261115/README.md) for its role in observability. Here it gives an explicit physical loss of closure of the norm pair even among equal-size sources. For lexicographic pair coding (bit 1 points from smaller to larger label), take the ordinal sums
\[
A=(5,76)\Rightarrow(6,88),\qquad B=(5,12)\Rightarrow(6,408).
\]
The four block triples \((q_1,q_2,q_3)\) are \((15,120,525),(33,360,2025),(11,100,465),(45,432,2277)\). Optional-slot counts multiply across an ordinal sum, since each slot is its A part followed by its B part. Thus both sources have \(q_1=495,q_2=43200\), but distinct \(q_3=1063125,1058805\). They are indistinguishable to **every unmarked first cyclic Hamilton-count reader** by the rank-two theorem. A *second* identical singleton wrap separates them: the two counts are 2341755 and 2333115. The difference is 8640, precisely twice the third-slot difference. Even the binary norm profile and source size together are not a continuation interface.

Iteration yields \(H(R^t(A))=t!q_{t+1}(A)\) plus a linear combination of \(q_1,\dots,q_t\), directly by multiplying the leading coefficients at successive grades. Thus all future scalar H readouts under one fixed literal wrap recover the entire optional-slot profile by triangular inversion. The linear table with rows \(t\) and columns \(A=TT_n\) has infinite rank: any dependence among its first m rows would give a nonzero linear combination of \(1^n,2^n,\ldots,m^n\) vanishing for all \(n\ge1\), contrary to Vandermonde. This rules out a universally sufficient **fixed-dimensional linear readout state** for unbounded wraps. It does not rule out nonlinear encodings, special bounded families, or preserving a patch only for a finitely prescribed set of future readers.

## Exact criterion for when the analogy works

The common mechanism is **defect observability, not the number eight**. Cayley doubling pairs a previously unobserved associator with a new input; the tournament wrapper makes a previously unobserved third-slot count visible. Compensation on the visible scalar alone can be valid for one layer and fail at the next, with the failure pinned to a specific operation and coordinate. By contrast, the odd-grade congruence survives the rank change because it is proved directly from the parity of *each newly visible coordinate*. That is the test to apply to proposed Chan–Pak closure: prove that the repair is neutral on the continuation coordinates the next allowed wrapper reads, or build a second controllable repair for their defect. This theorem does not contribute a universal exact target selector.
