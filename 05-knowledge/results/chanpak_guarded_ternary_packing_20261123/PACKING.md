# Five ternary seed residues in the eleven-state golden packing reader

This supplies an exact finite packing controller for the centered ternary arithmetic modulo 5, with a distinguished role for the sixth seed \(1\). It also proves a separate physical no-go for the one-reversed-bridge tournament language underlying the prime support \(105=3\cdot5\cdot7\). Neither is a claim that the tournament-spectrum conjecture is solved.

## 1. The five initial states come from normalized heights, not raw residues

At the deleted ternary branch \(7\to21\), the nontrivial seeds in planar order are \(19,61,63,65,23\), respectively at relative depths \(1,2,2,2,1\). Normalize a height by \(u=3^{-r}H\) at depth \(r\), then reduce modulo 5. The five seed values become
\[
 (19/3,61/9,63/9,65/9,23/3)\equiv(3,4,2,0,1)\pmod5,
\]
a bijection with \(\mathbb Z/5\), whereas raw heights have repeated residues. The arithmetic mirror is \(u\mapsto14-u\), so on these residues it is \(i\mapsto4-i\), the pentagonal reflection fixing \(i=2\) (height 63).

Use the eleven states of the inherited [cyclic guarded reader](../session91_packing_golden_carriers_20261029/operators/LUCAS_FIVE_GOLDEN_OPERATOR.md), the independent sets of \(C_5\). Place the five normalized residues on its five singleton states \(\{i\}\). One legal nonempty double-toggle moves
\[
 \{i\}\longrightarrow\{i,i\pm2\}\longrightarrow\{i\pm2\}.
\]
Each insertion obeys the guard that no adjacent site is occupied. A shift by \(1\) or \(-1\) is two such double-toggles in the opposite sign, since \(\mp4\equiv\pm1\bmod5\). These are the shortest paths avoiding the empty state: parity makes a path between singletons even, and any two-toggle path avoiding empty has the displayed independent pair as its middle.

For digit \(d\in\{-1,0,1\}\), the absolute arithmetic would update \(H'=3H+2d\). The normalized residue updates by
\[
 u'=u+2\,3^{-(r+1)}d\pmod5.
\]
Hence the elementary nonempty toggles above enact every centered ternary word modulo 5, with a four-phase clock \(r\bmod4\). Each digit can use exactly four elementary toggles, padding by a legal add-delete loop when needed. A shift by \(\pm1\) occurs and needs four, so that common budget is sharp under these constraints. If the five singleton coordinates at each clock phase are required bijective and the neutral digit has no net movement, a periodic clock of length \(L\) requires \(3^L=1\bmod5\): composing the neutral condition \(\phi_{r+1}(3H)=\phi_r(H)\) once around the clock proves this by injectivity. Therefore four phases are minimal in this model.

The five-state tape does **not** remember an absolute value. In the original integer controller, a word determines infinitely more than a finite field residue. For example at depth zero different positive heights already collide modulo 5; after a common suffix their actual difference is multiplied by a power of 3.

## 2. The literal 11 and the distinction between five and six

The inherited golden observation \(\rho(I)=\sum_{i\in I}4^i\bmod11\) sends the singleton states to the nonzero squares \(QR_{11}\). Independent pairs are precisely \(\{i,i+2\}\) as unordered pairs, and map to \(6\cdot 4^i\), exactly the nonsquares; the empty state maps to 0. Thus the controller's step on singleton golden labels is literal multiplication by \(4^{\,2\,3^{-(r+1)}d}\), where the exponent is taken modulo 5. Each nonempty double-toggle passes through the nonsquare sector, one legal insertion and one legal deletion; it is not claiming that every such transition is an oriented Paley edge.

Among the five residues \(4^i\), Paley orientation (difference is a square) is \(i\to i+1,i+2\): divide the differences by \(4^i\) and check \(4-1=3,\ 4^2-1=4\). The pentagonal reflection is the projective involution \(x\mapsto3/x\) on \(QR_{11}\), an anti-automorphism of this induced regular five-tournament. It fixes precisely \(x=5\) in the square subgroup, the residue of the normalized 63 state.

The external seed \(1\) can be represented by the empty packing state 0, consistently distinguishing it from the five nonempty initial states; using its normalized residue instead would cause a collision with 23. In the induced six-vertex Paley host \(0\cup QR_{11}\), vertex 0 is a universal source. The integer \(H\) of the five- and six-vertex hosts is 15, not one of the five seed values. However, the five-state mirror cannot extend as a converse isomorphism of the six-vertex host, which has a source but no sink. In projective coordinates the obstruction is literal \(0\leftrightarrow\infty\): if one adjoins a formal seventh vertex \(\infty\) as universal sink, the same involution extends to a converse isomorphism. This separates the role of a distinguished vacuum port from the cyclic state symmetry.

## 3. The 105 contraction cannot be the packing update

In the [seven-letter golden bridge language](../chanpak_55_105_golden_20261110/golden/NOTE.md), \(105\) is exactly the squarefree prime support of uniform divisibility of the three-coordinate observer. In particular the literal word \(W_5=q_2 c c q_2 s s s s q_2\) has matrix \(5U_5\), with \(U_5\) unimodular. Its raw mod-5 action is constant zero, whereas every packing digit update, including its clock, is a permutation. More formally, if a state-observation map \(\pi\) into the packing controller intertwines this raw contraction with any word of its digit permutations, then \(\pi(0)=\sigma(\pi(v))\) for every \(v\), with \(\sigma\) bijective. Hence \(\pi\) is constant on the set of inputs on which the intertwining is required. This particular collapse is avoided if the normalized output \(vU_5\), rather than raw mod-5 zero, is retained. It does not show that the resulting normalized coordinates transport the packing metric or produce tournament ternary transitions.

There is in fact a stronger all-order absolute obstruction to this physical mechanism. Take a marked source \(T,a\) with \(H=H(T)>1\), let \(E\) count Hamilton paths ending at \(a\), and let \(C\) count rooted ordered two-slot covers in which a designated slot starts at \(a\) (the other slot may be empty). Then \(C\ge H\): split any Hamilton path just before \(a\), recovering it by concatenation in the original order. Append any nonempty tournament \(B\), making all cross arcs point old to new except \(b\to a\) for a fixed \(b\in B\). In the inherited [single-bridge formula](../golden_hamilton_transfer_20261103/RECURRENCE.md), put \(h=H(B), S=\#\) full paths of \(B\) starting at \(b\), and \(D=\#\) pointed two-slot covers with designated slot ending at \(b\). Then
\[
 H(T\star B)=hH-ES+CD,
\]
where \(D\ge h\) by the same split argument, \(S\le h\), \(E\le H\).

- If \(h\ge5\), the correction \(CD-ES\) is nonnegative, so the result is at least \(5H\).
- If \(h=3\), \(B\) has a unique nontrivial SCC, a directed triangle, with singleton SCCs before or after it. If \(b\) is the initial singleton, \(S=D=3\), and the result is \(3(H-E+C)\), a multiple of three. If \(b\) is in the initial triangle then \(S=1\), \(D=3\): the triangle has one full path ending at \(b\) and two two-slot covers with designated path ending there, and all later singletons belong to the other slot. The result is at least \(5H\). Otherwise \(S=0\) and the result is at least \(6H\).
- If \(h=1\), \(B\) is transitive. Let \(k\) be the number of vertices before \(b\). Then \(S=[k=0]\), \(D=2^k\): each earlier vertex independently joins the designated slot or the other slot, and no later vertex can join the first. For \(k\ge2\) the result is at least \(5H\). The only two low possibilities are \(H-E+C\) and \(H+2C\) at \(k=0,1\), respectively.

Since \(3H+2<5H\), an exact negative digit \(3H-2\) must use the first low possibility and forces \(C=2H+E-2\). Then \(H+2C\ge5H-4>3H+2\) for \(H\ge5\); the same first low possibility cannot give both signs. Thus **no one fixed marked source of count at least 5 can have both signed exits even when arbitrary blocks \(B\) are allowed**, as long as there is exactly one reversed cross contact from the entire appended block to the old source. This covers the seven-letter language regardless of word length, by taking all fresh vertices together as \(B\): later one-bridge edits internal to the appended part do not change its single reversed contact with the old source. It does not rule out reinstalling the old port, extra reversed contacts, or a different renewable ready class using more general physical operations.

## 4. Golden eigenvalues and exact two-step mixing

Let \(R\) cyclically shift the singleton labels, and choose the three digits independently and uniformly. The unsigned magnitudes of the permitted residue shifts alternate between 1 and 2 across clock phases, so the stochastic matrices alternate
\[
 K_1=(I+R+R^{-1})/3,\qquad K_2=(I+R^2+R^{-2})/3.
\]
The two nontrivial eigenvalues are \(\varphi/3\) and \(-\varphi^{-1}/3\), each twice, with \(\varphi=(1+\sqrt5)/2\); the roles on each fixed nontrivial Fourier mode are exchanged between \(K_1\) and \(K_2\). Their product on every such mode is \(-1/9\). Equivalently the nine two-letter words produce signed sums ranging from \(-3\) to \(3\): zero occurs once and each of the other four residues occurs twice after reduction modulo 5 (the two sums \(\pm1\) already each occur twice). Thus
\[
 K_1K_2=(2J-I)/9
\]
where \(J\) is the all-ones matrix. Starting at any singleton, after \(2k\) random digits the law is exactly uniform plus \((-1/9)^k\) times its initial deviation from uniform. The golden eigenvalues are therefore forced by the exact arithmetic/packing transfer, and their product eliminates the pentagonal distinction at two steps except for an explicitly decaying, alternating remainder. This is mixing of the finite guard controller, not a distribution result for attainable Hamilton heights.

## A golden mixing law, not merely a golden coordinate

Let \(R\) rotate the singleton states by one, and give the three
centered digits equal weight. At a phase whose absolute normalized step
is \(\pm j\), the transition is \(K_j=(I+R^j+R^{-j})/3\).
The four-phase clock alternates the two step classes \(j=1,2\).
On the four nonconstant Fourier modes the eigenvalues are
\(\varphi/3\) and \(-\varphi^{-1}/3\), swapping which modes carry them
between phases. Indeed \(1+2\cos(2\pi/5)=\varphi\) and
\(1+2\cos(4\pi/5)=-\varphi^{-1}\).
If \(J\) is the all-ones five-by-five matrix, then exactly
\(K_1 K_2=(2J-I)/9\): among the nine two-digit words, one stays
and two reach each other state. Thus after \(2k\) steps the
probability differs from uniform by \((-1/9)^k\) times its
initial centered law. This is an exact mixing law for finite
packing states, not a bound on the distribution of absolute counts.
