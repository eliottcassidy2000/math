# Exact endpoint decoding and four more physically supplied scales

## Every positive integer is an exact marked endpoint count

Let \(N=\sum_{i=1}^r b_i2^{i-1}>0\) be its ordinary binary expansion. Take an increasing chain \(c_0,\ldots,c_{r+1}\), a vertex \(y\) beating this whole chain, and an apex \(x\). Set \(x\to y,c_0\), \(c_{r+1}\to x\), and \(c_i\to x\) iff \(b_i=1\) for \(1\le i\le r\). Deleting \(x\) gives a transitive order. The whole tournament is strong by the cycle \(x,c_0,\dots,c_{r+1},x\) and the access through \(x,y,c_0\).

A Hamiltonian path **starting \(c_0\) and ending \(c_{r+1}\)** must use \(x\to y\) consecutively, since \(y\)'s only in-neighbor is \(x\) and it cannot start. The last chain vertex before \(x\) is some \(c_i\) with \(b_i=1\). Every \(c_j\) with \(1\le j<i\) can choose its side of \(x,y\) independently, later vertices must be on the far side, and both runs are increasing. This gives exactly \(2^{i-1}\) paths for that bit. Summing over the bits gives \(N\) exactly on \(r+4\) vertices. This marked construction is simpler than the separate three-run/perfect-matching encoder in [absolute marked targets](ABSOLUTE_MARKED_TARGETS.md). It still does not count all paths in the tournament.

## A coefficient that really becomes the whole Hamilton count

Let \(U\) have a marked exit \(z\), and add a vertex \(v\) whose only out-neighbor is \(z\). Let
* \(H=H(U)\), \(h=H(U-z)\);
* \(E\) be the number of Hamilton paths on \(U\) ending at \(z\);
* \(L\) be the number of ordered two-slot covers on \(U\) with the second slot nonempty and starting at \(z\), the first allowed empty (the earlier rooted relaxed response on \(U\)).

In a compressed Hamilton path, a marked occurrence of \(v\) must either be last or be followed by \(z\), so no third occurrence is possible. The compressed one-visit count is \(H-E+L\): a terminal mark has \(H-E\) valid preceding paths, and otherwise there is a split with suffix starting at \(z\). The two-visit count is \(L-h\), because the second mark must be last and the path starting at \(z\) on the other side may not be just the singleton \(z\). Replace \(v\) by a length-\(s\) transitive block. The count of allocations to one or two nonempty slots is \(1\) or \(2^s-2\), respectively. Thus for every \(s\ge1\) the **entire** count is
\[
C+D2^s,\qquad D=L-h,\quad C=H+2h-E-L.
\]
The same formula with entry and exit reversed applies by converse to a sole-entry mark.

These literal strong frames all give \(C=1\). Rows refer to the full frame (including its marked vertex) and the displayed four readings refer to its exterior \(U\) rooted at the unique exit.

| \(D\) | Full frame rows | mark and exit | exterior \((H,h,E,L)\) |
|---:|---|---|---|
|135| \(258,276,9,291,301,327,415,319,4\) | \(8\to2\) | \((115,41,20,176)\) |
|405| \(262,28,424,337,293,331,279,379,2\) | \(8\to1\) | \((379,83,56,488)\) |
|1215| \(760,1013,769,518,44,652,956,540,697,16\) | \(9\to4\) | \((1095,309,188,1524)\) |

Thus \(1+D2^s\) for each of these three \(D\)'s and every \(s\ge1\) is an **unmarked strong count**, continuing the earlier offset-one inverse chain of slopes \(15,45\). When \(D\) is divisible by three and \(s\ge2\), the \(k=1\) odd-chart inverse of \(1+D2^s\) is \(1+3D2^{s-1}\); its greedy guard is valid, since \((H+1)/2\equiv1\bmod3\) and the predecessor menu contains 1. At the old slice's coordinate \(a=9\) this gives the actual chain \(961,1441,2161,3241,4861\). This is finite progress on a recursive lane, not a proof that all future slopes are supplied.

The eight-vertex strong frame \((2,8,3,21,103,79,143,63)\) has a sole **entry** at mark 7. Its one- and two-visit counts are 69 and 34, respectively, so its literal transitive inflation gives \(1+34\cdot2^s\). The coefficient 34 is consequently a physical shift, not just the label \((69-1)/2\). Deleting the two vertices 3 and 6 from this frame (or any such inflation) leaves a transitive spine. Their orientation word, in transitive order, is the clone block \(01\), followed by \(10,00,11,11,01\). The only possible common three-run open interval would start after the first \(10\) and end before the following \(00\): they are consecutive. Thus the whole rail has a two-hub realization without ternary sites, even though other 69 realizations use them.

Literal replays verify the listed first- and second-visit readings, several full inflations and arbitrary small endpoint decoding. The identities give the unbounded conclusions.
