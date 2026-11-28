# Exact density of square shells without a one-apex +2 flip

2026-11-28. The [branching Mills theorem](../chanpak_branching_mills_shells_20261128) gives two distinct strong one-apex counts inside every odd square shell for \(h\ge5\). The same two counts cannot always be required to differ by a single arc reversal and by two in height. The limiting fraction of odd parents without any such **one-apex** pair is
\[
\frac38-\frac{1}{2\sqrt2}=0.0214466094\ldots .
\]
The restriction is to tournaments each becoming transitive upon deletion of one vertex. Their deletable vertices are **not** assumed to be the same; the lemma below proves that every one-arc pair in the class has a common choice. No obstruction is asserted for a pair in which either endpoint is outside this class.

## 1. A common apex always exists

**Lemma.** Any two one-arc neighbors in the one-apex class have a common deletable vertex. When the apex is initially marked differently at the two ends, one can choose a common mark.

Choose a deletable \(a\) for \(T\) and a deletable \(b\) for \(T'\); they are distinct. The changed arc \(ij\) cannot touch either \(a\) or \(b\), because then that vertex would be deletable at both ends. The two restrictions to \(C=V\setminus\{a,b\}\) are transitive orders differing only at \(ij\). Thus \(i,j\) are consecutive there, in order \(i,j\) in \(T\) and \(j,i\) in \(T'\).

In the transitive order \(T-a\), vertex \(b\) must sit between \(i,j\); if it did not, reversing their arc would retain transitivity after deleting \(a\), contrary to the hypothesis. In \(T'-b\), vertex \(a\) must similarly sit between \(j,i\). All four vertices therefore have uniform earlier and later transitive flanks, and the core contains the unchanged cycle \(i\to b\to j\to a\to i\). If \(a\to b\), deleting \(j\) from either tournament leaves the transitive core \(a\to i\to b\) with \(a\to b\). If \(b\to a\), deleting \(i\) leaves \(b\to j\to a\) with \(b\to a\). In both cases the changed arc disappears and the flanks remain transitive, so there is a common deletable vertex after all.

## 2. Complete large-height normal form

Use a common apex and bit words on a transitive spine. Let \(F(x)\) denote the count of the ordinary binary word of \(x\), with \(F(0)=1\). For an apex-contact reversal, write \(u0v,u1v\); \(x\) is the value of \(u\), and \(y\) the value of the reversed complement of \(v\). The [middle-contact flip formula](../chanpak_branching_mills_shells_20261128) gives the signed difference \(y+[y\text{ odd}]-x-[x\text{ odd}]\). Reverse-complementing exchanges the roles of \(x,y\) and negates the difference; take it to be +2. If \(x=0\), then \(y=1\) or 2 and the counts are bounded; \(y=0\) is impossible with the chosen sign. Otherwise uniquely write \(x=2p-a,\ y=2p+2-b\) with \(p\ge1,\ a,b\in\{0,1\}\). Leading zero or trailing one padding can be deleted, as it supplies only transitive flanks. The possible lower counts are exactly
\[
L_{ab}(p)=2(2p-a)(2p+2-b)+2p+F(2p-a)+F(2p+2-b)-1.
\]
Indeed appending the central zero adds \(x+[x\text{ odd}]=2p\); its doubled even prefix then contributes \(2xy\) to the suffix count. Every displayed count is produced by the literal apex edit with ordinary, unpadded binary words.

If a spine arc changes, the two transitive spine orders can differ only by swapping consecutive sites. If their contact bits agree, there is no height change. If not, compare \(u10v\) and \(u01v\), with the same definitions of \(x,y\). The prefix recurrence gives a difference \(x+2-[x\text{ odd}]\) from the two central positions, plus \(y-[y\text{ odd}]\) from the suffix, for a total \(x+y+2-[x\text{ odd}]-[y\text{ odd}]\). This is 2 only if \(x,y\le1\). After removing flanks, these have at most four sites on the spine and so give only bounded heights. This completes the classification apart from bounded small cases that cannot affect the density or the example below.

## 3. Four parity choices make one cluster

Let \(D(n)=F(n+1)-F(n)\). From \(F(2m)=F(m)+m+[m\text{ odd}]\) and \(F(2m+1)=F(m)\) one gets
\[
D(2m)=-m-[m\text{ odd}],\qquad
D(2m+1)=D(m)+m+1+[m+1\text{ odd}],
\]
hence \(-n/2-1\le D(n)\le n+O(\log(n+2))\) by induction. In the quartet at \(p\), the initial step from \(L_{11}\) in either coordinate is respectively \(4p+2+A\) and \(4p-2+C\), where
\[
A=p+[p\text{ odd}]+D(p-1),\qquad
C=p+1+[p+1\text{ odd}]+D(p).
\]
Both \(A,C\) are between \(p/2-O(1)\) and \(2p+O(\log p)\); one has upper bound \(p/2+O(1)\) because one of \(p-1,p\) is even. The step gains an additional 2 if the other coordinate has already increased. The minimum and maximum are \(L_{11}\) and \(B_p=L_{00}\); for all \(p\ge1\) the two initial steps are positive by the same lower bound for \(D\). The first internal gap is the smaller step, the middle is their absolute difference, and the last is the smaller step plus 2. Thus all internal gaps are \(<5p\) for large \(p\).

In contrast the inter-cluster gap is exact:
\[
B_p=8p^2+12p+1+F(p)+F(p+1),\qquad
L_{11}(p+1)-B_p=6p+6.
\]
The formulas follow by direct substitution of the binary recurrence. In square-root coordinates the inter-cluster gap has limiting width \(3/(2\sqrt2)>1\), while every internal gap has width at most \(5/(4\sqrt2)+o(1)<1\). Only the gaps between clusters can swallow an entire large square shell.

## 4. Phase averaging and the density

The sequence \(\sqrt{B_p}\) is equidistributed modulo 2. To see this, freeze \(p\bmod Q=v\) with \(Q=2^k\) fixed. Write \(Y(v)=Q-1-\operatorname{rev}_k(v)\), where the reversal uses exactly \(k\) bits. Applying the prefix rule to the last \(k\) bits of \(p\) gives
\[
F(p)/p=Y(v)/Q+O(1/Q+Q/p).
\]
Indeed, if \(p=xQ+v\), the cross term is \(xY(v)\); the two word counts contribute \(O(x+Q)\) and the boundary correction is at most one. The same lemma applies to \(p+1\), using its residue \(v+1\bmod Q\). Expanding \(\sqrt{B_p}\) now gives \(2\sqrt2\,p\) plus a constant depending only on \(v\), up to \(O(1/Q+Q/p)\) with absolute constants. On each residue class the leading phase modulo 2 is an irrational rotation and therefore uniform: the averages of each nonconstant Fourier mode tend to zero by the geometric-series formula. Continuous periodic tests can be approximated by trigonometric polynomials. Their uniform continuity allows first taking \(p\) to infinity, then \(Q\) to infinity, and approximating interval indicators by continuous tests. This proves the required phase law without an independence assumption on the binary digits.

For an odd parent \(h\), the strict shell admits a lower count \(L\) exactly when \(h^2<L<(h+1)^2-3\): the upper count is \(L+2\) and its upper floor endpoint is \(L+3\). For large \(h\), absence means the open interval is contained in an inter-cluster gap, since in square-root coordinates it has width \(1+o(1)\). Each gap captures exactly the odd integers in
\[
[\sqrt{B_p},\sqrt{L_{11}(p+1)+3}-1].
\]
Its width tends to \(\delta=3/(2\sqrt2)-1\); at most one odd integer lies there eventually. Endpoint perturbations of size \(o(1)\) do not affect the limit by the phase law, so the fraction of gaps that capture a parent is \(\delta/2\). There are \(X/(2\sqrt2)+O(1)\) gaps relevant up to parent \(X\), versus \(X/2+O(1)\) odd parents. The desired fraction is \(\delta/(2\sqrt2)=3/8-1/(2\sqrt2)\).

## 5. Exact example and replay

The first cubic pair in the earlier construction is \(B_{177}=252{,}965\) and \(252{,}967\). The next quartet starts at \(L_{11}(178)=254{,}033\). For odd parent 503 the square shell has endpoints \(253{,}009\) and \(254{,}016\) and lies entirely inside the gap. It still has the two distinct strong witnesses from the quadratic construction, but no one-arc +2 pair for **any** two one-apex hosts.

Run node replay.js. It independently verifies 1,000 literal apex-contact pairs using the inversion sum, checks the exact inter-cluster identity through \(p=5{,}000\), gives 1,077 misses in the complete quartet scan through parent 100,000 and 10,752 inter-gap misses through parent 1,000,000, and exhaustively checks the apex-change lemma on labeled tournaments through order six. The arithmetic scan concerns the large-count quartet catalogue; bounded exceptions are handled in the proof.
