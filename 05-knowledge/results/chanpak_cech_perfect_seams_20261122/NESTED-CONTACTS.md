# Exact classification of nested-contact unit seams

**PROVED.** Use the positive odd-cycle packing identity from [the overlap ledger](OVERLAP.md).

Let \(u\to v\) be the old arc, and let \(R=T-\{u,v\}\). For \(i\in R\), write \(a_i=[u\to i]\), \(b_i=[v\to i]\). Assume the contacts are nested and choose the names and direction of comparison so that \(a_i\ge b_i\) for all \(i\). The old arc \(u\to v\) is then the lower count, as follows. The two compared hosts can always be ordered this way, regardless of which orientation was supplied first.

**Theorem.** Write \(D=\{i:a_i>b_i\}\). Reversing \(u\to v\) to \(v\to u\) has normalized defect at least \(|D|\), and a positive unit defect occurs exactly when:
1. \(D=\{t\}\), and \(S=R-t\) is transitive;
2. on the transitive order \(i<j\) of \(S\), putting \(c_i=[t\to i]\), one has
\[
 c_i\le b_j,\qquad b_i\le c_j\quad(i<j).
\tag{6}
\]
After removing common initial and terminal singleton components, every such unit seam has order at most five. Its lower ordinary count is one of \(1,3,9,13\), and each occurs.

### Exact positive transfer

The only odd cycles affected by the flip are the old and new cycles containing the changed arc. Each such packing contains only one affected cycle, since they all use \(u,v\). For a directed path \(P\) **on an odd positive number of outside vertices** in \(R\), with initial vertex \(i\) and terminal vertex \(j\) (possibly the same singleton), there is an old odd cycle using \(P\) if \(b_i(1-a_j)=1\), and a new one if \(a_i(1-b_j)=1\). Its possible disjoint companions contribute exactly \(H(R-V(P))\). Hence
\[
 \frac{H(T')-H(T)}2
 =\sum_{P}
 \left(a_i(1-b_j)-b_i(1-a_j)\right)H(R-V(P)).
\tag{7}
\]
The bracket equals
\[
 (a_i-b_i)(1-b_j)+b_i(a_j-b_j)\ge0.
\]
Equivalently each lost cycle is injected into a gained one with the same vertex support, by exchanging the roles of \(u,v\) along the outside path. The positive weights are identical on their companions. Each singleton \(i\in D\) contributes \(H(R-i)\ge1\). If \(D\) is empty all terms vanish; otherwise the jump is positive.

A unit therefore forces one mismatch \(t\) and \(H(R-t)=1\). A tournament with \(H=1\) has no directed triangle by the packing formula, so it is transitive. The only other possible positive terms in (7) are paths on at least three outside vertices starting at \(t\) and ending at a vertex with \(b_j=0\), or starting where \(b_i=1\) and ending at \(t\). In the transitive \(S\), such a path can be shortened to three vertices: in the first case retain its initial outside neighbor and its final vertex; in the second retain its first vertex and the predecessor of \(t\). So these terms vanish exactly when no \(t\to i\to j\) with \(b_j=0\) and no \(i\to j\to t\) with \(b_i=1\) exists for \(i<j\) in \(S\). This is precisely (6). Conversely (6) leaves the singleton \(t\) as the only extra term, of weight 1, proving sufficiency.

### The forced two-tape shape

Condition (6) says the ordered pairs \((b_i,c_i)\) have the form
\[
 00^*\;\bigl(\text{at most one }10\text{ and at most one }01,
                  \text{ in either order}\bigr)\;11^*.
\]
A \(00\) cannot follow any non-\(00\), and a \(11\) cannot precede any non-\(11\), by the cross inequalities. Repeating \(10\) would violate \(b_i\le c_j\); repeating \(01\) would violate \(c_i\le b_j\). Conversely every displayed tape obeys (6). All \(00\) vertices beat \(u,v,t\) and their later outside vertices, giving common initial singleton components. All \(11\) vertices give common terminal singleton components. The remaining core is \(u,v,t\) and at most two mixed vertices. In the upper graph the three vertices \(v\to u\to t\to v\) form a directed triangle, and a \(10\) vertex joins its strong component via \(v\to i\to t\), while a \(01\) vertex joins via \(t\to i\to u\). Thus the reduced upper graph is strong and there is no further common cut.

To read the small core, use \(u\to t\to v\), old \(u\to v\), and the printed tape contacts and outside order. The old cycle inventories and the resulting counts are:

| mixed tape | old triangles | old directed 5-cycles | old and new \(H\) |
| --- | ---: | ---: | --- |
| empty | 0 | 0 | \(1,3\) |
| \(10\) | 1 | 0 | \(3,5\) |
| \(01\) | 1 | 0 | \(3,5\) |
| \(10,01\) | 4 | 2 | \(13,15\) |
| \(01,10\) | 3 | 1 | \(9,11\) |

For an explicit check of the last two rows, label \(u,v,t\) by \(0,1,2\) and the two mixed vertices by \(3,4\) in outside order. In the \(10,01\) case the triangles in cyclic order are \(024,034,132,134\), and the 5-cycles are \(01324,02134\). In the \(01,10\) case the triangles are \(023,142,234\) and the 5-cycle is \(01423\). No two odd cycles are disjoint at these orders, so the old count is 1 plus twice the total; the new count is exactly two more by (7). Thus any literal \(+2\) seam other than these four lower heights necessarily has **incomparable** flipped-endpoint out-neighborhoods within the outside tournament. It must use genuine signed cancellation rather than a one-sided nested-cycle injection.

## 3. Any bounded positive nested jump has a logarithmically bounded primitive core

The unit classification can be made robust beyond a jump of two. Suppose the normalized positive jump in (7) is \(m\ge1\), and the pair has no common proper ordinal cut. Then
\[
 |V(T)|\le 8+5\log_2(m+1)+\log_3 m+3\log_5 m.
\tag{8}
\]
Thus at each bounded jump there are only finitely many primitive nested-contact pairs, even without a bound on the starting scalar count.

Choose any mismatch \(t\). Its singleton term in (7) proves \(H(R-t)\le m\). Nontrivial strong components of \(R-t\) have path counts \(q_i\ge3\) whose product is \(H(R-t)\), and [Busch's inequality](https://doi.org/10.37236/1141) bounds their sizes by \(1+3\log_5 q_i\). Their total number of vertices is at most \(\log_3m+3\log_5m\). All remaining vertices of \(R-t\) are singleton components, transitively ordered. For these use \(c_i=[t\to i]\) and divide them into five classes: all-zero or all-one triples \((a_i,b_i,c_i)\); the two common-contact types \(a_i=b_i\) with \((b_i,c_i)=(1,0),(0,1)\); and the mismatches \(a_i>b_i\).

If \(R-t\) is nonempty, its first strong component contains some vertex \(i\) for which \((a_i,b_i,c_i)\ne(0,0,0)\): otherwise that component would beat the rest of \(R-t\) and also \(u,v,t\), forming a common cut. Every all-zero singleton lies in a later component. Fix the witness \(i\). Any odd nonempty subset of these all-zero singletons, in its increasing order, can appear in a positively contributing outside path. If \(c_i=1\), use \(t,i,\text{subset}\); if \(c_i=0\), then \(a_i=1\), so use \(i,\text{subset},t\). The arcs follow the component order and the all-zero contact with \(t\); the bracket in either case is 1 by (7). Thus \(r\) such singletons give at least \(2^{r-1}\) units of the normalized defect when \(r>0\). The symmetric last-component argument applies to all-one singletons: fix \(i\) in the last strong component not of type \(111\); if \(c_i=0\) use \(\text{subset},i,t\), and otherwise \(b_i=0\) so use \(t,\text{subset},i\).

For a common-contact \(10\) type, each nonempty even subset of that type in increasing order followed by \(t\) is a positively contributing path, giving \(2^{r-1}-1\) terms. For common-contact \(01\) place \(t\) first instead. Finally, every odd nonempty subset of mismatch singletons is itself a positively contributing outside path, because its initial vertex has \(a_i=1\) and its final vertex \(b_j=0\). All these bounds hold separately, so each of the five singleton classes has size at most \(1+\log_2(m+1)\). The case \(R-t=\varnothing\) needs none of these arguments. Adding the ports and \(t\) proves (8).

If a nested pair is not primitive, first split off common ordinal components. Each removed component multiplies the whole defect by its positive path count, so the middle defect is no larger than the original one. The statement therefore bounds the genuinely active nested interval, not the arbitrary common padding outside it. It does not assert such a bound for incomparable contacts; equation (7) would then be a cancellation of signed, not positive, contributions.
