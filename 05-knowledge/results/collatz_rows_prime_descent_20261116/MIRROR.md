# Two physical edits can cancel the large responses exactly

This closes the two-seam test in the [main note](README.md), not the
arbitrary target problem. Let \(U\) be a binary spine word beginning in 1,
\(V=\operatorname{RC}(U)\) its reverse-complement, and take the two
central words
\[
A=010010=01\,00\,10,\qquad B=100001=10\,00\,01.
\]
In the same natural vertex labeling, flip the old arcs inside the two
disjoint adjacent pairs, then read the transitive spine in its new order.
This is exactly the move \(UAV\to UBV\): the apex contacts have not changed
as arcs, only their reading order. Both tournaments are strong.

Use the [concatenation rule](../one_apex_all_residues_20261116/README.md).
The central profiles are \((X,Z,H)(A)=(18,45,15)\) and
\((X,Z,H)(B)=(33,30,17)\). Their lengths and \(X+Z=63\) agree.
Since \(X(U)=Z(V)=x\), the exterior contribution from the central
position values is \(63x\) in either word, and the exterior
cross term is \(64x^2\). The two boundary adjacency contributions
sum to \(p\), where \(p\) is the last bit of \(U\), in both words.
Thus the exact lower and upper counts are
\[
64x^2+63x+2H(U)+13+p,\qquad
64x^2+63x+2H(U)+15+p.
\]
This is a physical \(+2\) move at arbitrarily large counts,
despite the \(\Omega(\sqrt H)\) lower bound for a single adjacent
swap. More generally, with unrelated outside words the difference
is \(15(Z(V)-X(U))+3-p-q\), where \(q\) is the first bit of \(V\).
Matching the absolute exterior values, not merely their residues,
is what cancels the large response. The same reasoning localizes
changes between any two equal-length palindromic central words on
a mirrored exterior: \(X+Z=2^{|A|}-1\) and the boundary totals agree.

At \(U=1100\), \(x=12\), the graphs on fifteen vertices have twin
prime counts 10007 and 10009. The replay verifies both by full
Hamiltonian-path subset DP, reconstructs the two literal old arc
reversals, and checks primality by trial division. This does not
prove infinitely many twin primes. Even varying every \(U\), this
particular six-bit switch gives only \(O(\sqrt B)\) windows below
\(B\), since \(x\) is their ordinary positive binary address.
The local adjustment is exact, but positioning its baseline at an
arbitrary prescribed prime remains unsolved. In fact even the full
mirror-palindrome template \(U\,W\,\mathrm{RC}(U)\), with arbitrary
nonempty palindromic \(W\), has only \(O(\sqrt B\log B)\) outputs up to
\(B\). For length \(m\), the count is at least \(2^m x^2\), with
\(x\ge1\), and there are at most \(2^{\lceil m/2\rceil}\) palindromes.
Multiplying by the at most \(\sqrt{B/2^m}\) choices of \(x\) gives
\(O(\sqrt B)\) per \(m\), and there are at most \(\lfloor\log_2 B\rfloor\)
lengths. It is the equality of reflected exterior values that makes
the correction exact and also restricts where it can live.

Run python3 mirror.py from this directory.
