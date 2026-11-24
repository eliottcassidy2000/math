# A direct weighted Hamilton-cycle to tournament-path projection

**PROVED.** The generic directed Hamilton-cycle polynomial on \(m\ge3\) labels is a \(0,1,\) variable projection of the independently weighted Hamilton-path polynomial of one fixed **strong** tournament on \(2m\) vertices. The identity holds over every commutative ring. This transfers any arithmetic-circuit or border lower bound for the former polynomial to the latter, at a factor of two in vertex size. It gives no Boolean or unweighted Chan–Pak lower bound.

The attachment *VP ≠ VNP over Every Field* claims such lower bounds in Theorems 1.1 and 1.4, and in §2.5 explains the Boolean/formal distinction. We use the supplied claim **conditionally**; neither this reduction nor its replay audits the claimed theorem or its Lean development.

## Directed split vertices, on a fixed strong orientation

For each label \(i=0,\ldots,m-1\), use vertices \(I_i,O_i\). Fix arcs \(I_i\to O_i\) and \(O_i\to I_j\) for \(i\ne j\), and orient arcs within each of the I and O sets transitively in any order. This is a tournament. It is strong: between any distinct index pair the four vertices form the directed cycle \(I_i,O_i,I_j,O_j\); these strongly connected pairs overlap.

In the full weighted Hamilton-path polynomial, specialize:
\[
 I_i\to O_i:\ 1 \ (i\ne0),\quad I_0\to O_0:\ 0,\quad
 O_i\to I_j:\ x_{ij}\ (i\ne j),
\]
and give all within-type arcs weight zero. A term not identically zero must start at \(O_0\) (no support arc comes in) and end at \(I_0\) (no support arc goes out). At each other I vertex its only supported exit is to the matching O. Hence the word has form
\[
 O_0,I_{i_1},O_{i_1},\ldots,I_{i_{m-1}},O_{i_{m-1}},I_0
\]
for a permutation of \(1,\ldots,m-1\). This is exactly the directed Hamilton cycle cut at label 0; it has the identical monomial, and each cycle appears once, not \(m\) times. The proof is a termwise bijection, not an identity inferred only on 0/1 inputs.

A formal-arithmetic circuit over a field specializes by replacing input gates, with no extra arithmetic gates. The same argument specializes the coefficientwise approximation over the Laurent-series field. If the attachment's bounds hold for the Hamilton-cycle polynomial, then they hold at size \(2m\) for this fixed tournament's independently weighted Hamilton-path polynomial. Setting some arc weights zero is not the same operation as reversing arcs in the ordinary unweighted puzzle.

## Unconditional obstruction to a fixed two-read separation

The [SL₂ transfer](SL2-TRANSFERS.md) is exact only when the states really separate the chosen sides of the interface. A clean algebraic test is the rank across a partition of the variables (the matrix of coefficients, with monomials on each side as its indices). Consider an ordinal join of \(r\) directed triangles. In triangle \(i\), weight its three cycle arcs \(x_i,y_i,1\) and set all join arcs to 1. Then the formal weighted polynomial is
\[
 F(X,Y)=\prod_{i=1}^r(x_i+y_i+x_i y_i).
\]
Across the \(X|Y\) partition, the matrix on squarefree monomials is the \(r\)-fold Kronecker power of \(\left(\begin{smallmatrix}0&1\\1&1\end{smallmatrix}\right)\), so has rank \(2^r\) over every field (the factor has determinant \(-1\)). Thus any exact expression \(\sum_{j=1}^d f_j(X)g_j(Y)\) for this partition needs \(d\ge 2^r\); a two-dimensional weighted interface does not cover all placements. This is a separation-rank lower bound, **not** an arithmetic-circuit lower bound: the very same polynomial has a short product circuit. It isolates which feature of the attached question the transfer does and does not address.

Replay: run sl2-transfers-check.py with Python 3.
