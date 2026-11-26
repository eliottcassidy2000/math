# What the small Hamiltonian fibers change in the two-port and shifted-square constructions

This follows the [small-fiber classification](../chan_pak_count_fiber_grammar_and_hidden_cycle_trades_20261126/README.md) through the earlier [23-path seed](../chanpak_h23_golay_chains_20261124/README.md), the [unequal-chain inverse](../chanpak_h23_golay_chains_20261124/UNEQUAL-CHAINS.md), and the distinction in the [cycle seam](../../../session212/cycle-spectrum/PORTS-AND-GOLAY.md) between closing a path and selecting which contacts it reads.

## 1. Open ports explain exactly what is lost at 7, 9, 15, 17, 21, and 23 in this host

Take the six-vertex tournament with outgoing bit rows \((10,4,9,18,39,15)\), and replace vertices 1 and 5 by modules \(A,B\) with uniform external arcs. For a nonempty module \(M\), write \(D=H(M)\) and \(Q=P_M(2)\). A word can use at most two runs of each module, since the first module has just one possible external successor and the second just one possible external predecessor. The run weights are \((0,D,Q-2D)\), indexed by \(0,1,2\); crucially, the weight for an empty module is instead \((1,0,0)\). By enumerating the allowable words in the four distinct fixed vertices and \(i,j\) occurrences of the ports, the whole run tensor is
\[
S=\begin{pmatrix}5&9&4\\13&23&9\\6&11&4\end{pmatrix}.
\]
The old nonempty two-port formula is recovered from the bottom-right block. The boundary values are
\[
H(\emptyset,\emptyset)=5,\qquad H(\emptyset,B)=4Q_B+D_B,\qquad
H(A,\emptyset)=6Q_A+D_A.
\]

This makes the shifted-square boundary precise. Setting the two chain lengths to zero in the *nonempty* formula gives the formal value \(7=3^2-2\); actually deleting both ports gives **5**, because the support of the run tensor changes. Single-port deletion, with the other port a chain, gives the rays \(4\cdot2^v+1\) and \(6\cdot2^u+1\). In particular **9 and 17** occur on the first ray, at widths 1 and 2; the next count there is 33. The other ray starts at 13 and 25. With both ports occupied the minimum is **23**. Indeed \(D_M\ge1,\ Q_M\ge2^{|M|},\ Q_M\ge2D_M\), and \(4Q_AQ_B+3Q_AD_B+D_AQ_B-D_AD_B\) is strictly increasing in each coordinate when compared with the singleton pair. For a fixed cardinality of 1 or 2 the module is transitive. Thus none of **7, 15, 21** occurs anywhere in this host for arbitrary fillers, including empty ports. The count 15 exists globally by the previous classification: the additional hole at 15 is a control against interpreting a host-specific gap as a global obstruction.

The triangular-number connection here is exact but limited: if \(s^2=H+2\) then \(H+2=T_{s-1}+T_s\), two adjacent triangular layers. The smallest formal odd-root case lies in the wrong occupancy chamber, whereas the next odd-root case \(s=5\) gives the actual seed at 23.

## 2. Equal heights propagate through some interfaces and split at others

For any fixed host with **one** substituted homogeneous module \(M\), a skeleton of \(j\) runs reads exactly
\[
A_j(M)=\sum_{\ell=0}^j(-1)^{j-\ell}\binom j\ell P_M(\ell).
\]
This follows by assigning the vertices of \(M\) to the \(j\) nonempty ordered blocks, each carrying a Hamiltonian path. Thus even all possible homogeneous hosts of this kind cannot distinguish two modules with the same whole \(P_M(z)\). A root-marked operation within the module need not respect this equivalence.

Fill port \(A\) of the old seed with the two modules from each small fiber, and port \(B\) with a chain of width \(v\ge1\). The height is \((4Q+D)2^v+3Q-D\). The following table uses the concrete rows in the replay, including the source-padded six-vertex module at 15. The last column compares the sets of heights obtained by replacing **one chosen vertex of the module** with a two-vertex chain, over all possible choices.

| Original \(D\) | \(Q\) for the two modules | Heights in the old host | Rooted clone readouts |
|---:|---:|---|---|
| 9 | 84, 84 | Both \(345\cdot2^v+243\) | \(\{17,23\}\) vs. \(\{15,19,29\}\) |
| 15 | 216, 216 | Both \(879\cdot2^v+633\) | \(\{15,33,37\}\) vs. \(\{25,27,33,37,59\}\) |
| 17 | 224, 224 | Both \(913\cdot2^v+655\) | \(\{33,47\}\) vs. \(\{29,35,51\}\) |
| 23 | 296, 272 | \(1207\cdot2^v+865\) vs. \(1111\cdot2^v+793\) | \(\{47,51,137\}\) vs. \(\{41,45,53,59,65,81\}\) |

The height-17 trade is invisible to **every homogeneous substitution host**, but the single rooted clone separates the two modules for every root choice, even if the roots are chosen independently. This distinguishes a global color count from local geometric placement, just as the old Golay construction distinguishes the binary code from where its contacts sit. At height 9 the two modules differ by one arc. That arc therefore gives a neutral reversal at every height \(345\cdot2^v+243\). These are genuine different strongly connected isomorphism classes: the seed is strong and has no nontrivial proper module (checked on its 56 nontrivial proper subsets). In any substitution into this seed, a module meeting two or more component parts must meet all parts, since its support in the quotient is a module. If it misses a vertex in one of the parts, that vertex would have a uniform relation to every other part, making the corresponding quotient vertex a source or sink. Hence the parts are the maximal proper modules and any isomorphism must carry them to one another. Only \(A\) is nontransitive, so it must map to \(A\), which is impossible for the two nonisomorphic height-9 modules.

## 3. A new exclusion band for the unequal-chain shifted squares

For nonempty chains \(u,v\), the earlier formula gives
\[
H+2=4\cdot2^{u+v}+3\cdot2^u+2^v+1.
\]
The diagonal always works. Apart from the known \((u,v)=(7,3)\), **same-parity off-diagonal widths** can give a square only when
\[
v\ge3u+14 \quad (v>u),\qquad u\ge3v+6 \quad (u>v).
\]
This is an exact exclusion region, not a classification of the remaining far-off-diagonal region.

To see it, take \(X=2^{\min(u,v)}\), \(k=|u-v|/2\ge1\), and write the hypothetical odd root as \(2^{k+1}X+e\) with \(e\ge1\) odd. Then
\[
(C-2^{k+2}e)X=e^2-1,
\]
where \(C=2^{2k}+3\) when \(v>u\), and \(C=3\cdot2^{2k}+1\) otherwise. For \(k=1\), the first case is strictly between offsets 0 and 1, while the second is strictly between offsets 1 and 2. For \(k=2\) in the first case the only possible odd offset is 1 and fails; in the second case the possible offsets are 1 and 3, and only 3 works, giving \(X=8\). For \(k\ge3\), write \(C=2^{k+2}q+c\), with \((q,c)=(2^{k-2},3)\) or \((3\cdot2^{k-2},1)\). The equality forces \(e\le q\); since \(q\) is even, \(e\le q-1\). Hence \(X<q^2/2^{k+2}\), giving the stated two bounds by integrality.

For opposite parities, \(u\) even and \(v\) odd is ruled out modulo 3. If \(u\) is odd and \(v\) even, the four possible pairs modulo 4 leave only \(u\equiv1,\ v\equiv0\pmod4\), by testing squares modulo 5. Modulo 8 then rules out \(u=1\).

An exact interpretation of the earlier factorization \(4H+7=(2^{u+2}+1)(2^{v+2}+3)\) is a factor exchange: a square \(H+2=s^2\) requires repartitioning the prime powers into the two coprime bins \(2s-1,2s+1\), which differ by 2. The diagonal needs no exchange. The exceptional root 67 has original factors \(513=27\cdot19\) and \(35=5\cdot7\); the square witness exchanges their partners to give \(133=7\cdot19\) and \(135=27\cdot5\). That is a concrete connection to the repo's 19/67 theme, with no implication about class numbers. The independent integer-square replay finds no other off-diagonal pairs with both widths at most 256.

## 4. Collatz comparison: an exact comparison of readings, not an explanation by numerology

For the odd-step map \(x\mapsto(3x-1)/2^{v_2(3x-1)}\), the small heights separate into three backgrounds: \(9\to13\to19\to7\) enters the \(7,5\) cycle; \(15\to11\to1\) enters the fixed point; and \(21\to31\to23\to17\) enters the seven-node cycle starting at 17. Thus 7 and 21 do not share a Collatz orbit. The actual tournament relation between them is the multiplicative core obstruction \(21=3\cdot7\) together with the absence of a strong count-21 core. On the constructive side the two-port ray exactly exposes a fixed odd offset: \(h_{v+1}=2h_v-b\) for \(h_v=a2^v+b\). Like an inverse branch with a fixed residue condition, it is exactly invertible on that ray, but there is no reason it should cover all integers. In place of identifying the two dynamics, the rooted readouts give an actionable analogue: the apparent state \(H\) does not determine what a marked move will do next.
