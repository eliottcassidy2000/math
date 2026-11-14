# An exact unmarked smooth rail at every exponent pair

**Theorem.** Every target
\[
N=1+2^A3^B\qquad(A\ge3,\ B\ge0)
\]
is the Hamiltonian-path count of a **strong** tournament. For \(B\ge3\) it has a construction with just two hubs and a transitive spine. Thus a complete unbounded two-parameter family of absolute targets, not only their residues or marked contributions, is obtained. This does not yet assert closure under adding several different smooth terms.

## The two-variable cancellation that gives the rail

On an increasing spine put two hubs \(x\to y\). At each spine position use the digit \(2[x\to v]+[y\to v]\) in \(\{0,1,2,3\}\). The spine word is
\[
11\,0^k\,2^\ell\,311,\qquad k,\ell\ge1.
\]
The exact evaluation is
\[
H(k,\ell)=1+36\,2^k3^\ell=1+2^{k+2}3^{\ell+2}.
\]

Here is a finite certificate *for the identity at every size*, rather than a numerical pattern fit. In any ordering of the hubs a homogeneous spine block can enter any chosen nonempty subset of the three path runs. With its occupied-run set fixed, all outside vertices are earlier or later than the whole block, and every vertex in the block has the same hub adjacencies. Consequently either all assignments using exactly that occupied set work or none do. For a set of size \(r=1,2,3\), the number of assignments of a length-\(m\) block is the number of surjections onto \(r\) labels. For \(m>0\) these are \(1,\ 2^m-2,\ 3^m-3\cdot2^m+3\), respectively. Thus independently in \(k\) and \(\ell\), the exact output is a linear combination of \(1,2^m,3^m\). No higher coordinate or polynomial prefactor is possible. The values at lengths 1, 2, 3 give an invertible \(3\times3\) evaluation matrix in each variable.

The replay counts the literal tournaments by subset Hamiltonian DP and checks the entire \(3\times3\) grid:
\[
\begin{pmatrix}
217&649&1945\\
433&1297&3889\\
865&2593&7777
\end{pmatrix}.
\]
Each entry equals \(1+36\cdot2^k3^\ell\). By the already proved degree restriction this gives the identity at all lengths, exactly. The base \(k=\ell=1\) is strong (rows \(254,252,504,368,96,192,128,280,115\)); every growing word is obtained by transitive inflation of its two marked positions. Substituting nonempty modules into a nontrivial strong tournament preserves strongness. Choosing \(k=A-2,\ell=B-2\) proves the claimed sector.

The cancellation here is real: one block has only its dyadic dependence left, the other only the ternary dependence, and the constant is **one**. It is not the assertion that individual Hamiltonian frame monomials can be turned on without companions. In particular the ternary-depth lower bound does not conflict with the factor 36 in this particular board.

## The three low ternary exponents

* \(B=0\): the near-transitive block \(R_A\) has \(1+2^A\) paths.
* \(B=1\): inflate vertex 0 of the strong count-13 frame \((2,12,9,17,7)\) ([previous use](../session332_hidden_digits_and_arithmetic_diamonds_20261114/ODD_COORDINATE_PREDECESSORS.md)). Its fixed marked module can be visited once or twice, with literal visit coefficients \((13,6)\). At length \(s=A-1\) this gives \(13+6(2^s-2)=1+3\cdot2^A\), even for \(A=2\).
* \(B=2\): inflate vertex 4 of the strong frame \((2,4,41,19,7,27)\). The compressed visit coefficients are \((19,9)\), and none beyond two. Thus a length-\(A\) transitive inflation gives \(19+9(2^A-2)=1+9\cdot2^A\), even for \(A=1\). The replay verifies the compressed visits by enumerating exterior orders and marker positions.

These lower cases also stay strong. The theorem is an exact target compiler: a target of the stated form supplies its two exponents and thereby the full word. The size for the two-hub sector is \(A+B+3\). The set is still sparse because it is an offset-one copy of a restricted smooth monoid. The missing inverse is a way to combine distinct terms and reestablish the single baseline, rather than multiply their individual \(1+\text{term}\) counts.
