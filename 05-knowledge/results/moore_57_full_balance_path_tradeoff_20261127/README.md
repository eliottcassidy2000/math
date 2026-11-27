# Exact full balance and Hamiltonian counts in the rooted Moore encoding

2026-11-27. This is about the dense tournament encoding of a
hypothetical degree-57 Moore-perfect matching table, not about
existence of the underlying graph. The
[off-branch bound](../moore_57_slotwise_balance_20261127/README.md)
and the full maximum over distinct slots measure different things.
The full optimum is exactly **2912 for every quotient orientation**.
With near-regular interiors it is exactly **2914**. In the Paley
benchmark, attaining the off-branch target 160 would force more
than \(10^{1260}\) Hamiltonian paths, even before counting paths
that return to a branch they have left.

## 1. Same-branch entries fix the full objective

Use 57 branches of size 56, internal skew matrices \(C_i\),
and cross blocks \(S_{ij}(J-2P_{ij})\) for a quotient tournament
\(S\) and the Moore matchings. Let \(A\) be the full signed
matrix. Each cross block contributes \(-52J-4I\) to its diagonal
block of \(A^2\), since \(P_{ij}\) is a permutation matrix.
Thus for distinct slots \(a,b\) in the same branch,
\[
 (A^2)_{ia,ib}=(C_i^2)_{ab}-2912.
\]
Some \((C_i^2)_{ab}\le0\): its ordered off-diagonal mean is
\(1-\|C_i\mathbf1\|^2/(56\cdot55)<1\), and all entries are
even. This proves the 2912 lower bound independently of \(S\).

For the witness, take the bordered Paley skew conference
matrix \(P\) of order 28 from
\(\mathbb F_{27}=\mathbb F_3[x]/(x^3-x-1)\). The polynomial has
no root in \(\mathbb F_3\). With quadratic-residue arrows,
\(P^2=-27I\); the replay verifies this directly. Double it:
\[
 C=\begin{pmatrix}P&I+P\\-(I+P)^\top&-P\end{pmatrix}.
\]
It is a skew tournament matrix with \(C^2=-55I\), by multiplying
the commuting blocks. In the replay's vertex convention,
switching it by bitmask \(\mathtt{1cd175b1704278}\) (hexadecimal)
bounds every signed row sum by 11, with the same square.
Use this switched matrix in every branch.

For any quotient, \(|(S^2)_{ij}|\le55\). The prior
[slot formula](../moore_57_slotwise_balance_20261127/README.md)
bounds every off-branch entry in absolute value by
\(52|(S^2)_{ij}|+4+2R+4\) when all interior signed row sums
have magnitude at most \(R\). Here the first 4 uses the unique
Moore route to a slot. With \(R=11\) this is at most 2890,
strictly below 2912. This proves exact attainment for every
quotient; for the Paley-59 deletion it caps the off-branch part
at 186. In the Paley benchmark all-transitive interiors instead
have full maximum exactly 2966: consecutive slots have
\((C^2)_{ab}=-54\), all internal off-diagonal squares are at
least \(-54\), and the cross cap with \(R=55\) is 274.
A single transitive branch already forces the lower bound 2966.

## 2. The exact cost of requiring near-regular interiors

Here near-regular means \(R\le3\). We prove that any such
56-vertex skew tournament matrix \(C\) has an off-diagonal
entry of \(C^2\) at most \(-2\). Suppose otherwise and put
\(M=C^2+55I\). This is symmetric and entrywise nonnegative,
with zero diagonal. Every eigenvalue of \(M\) has even
multiplicity by real skew-symmetry of \(C\) and its
nonsingularity. Indeed any even principal tournament skew
matrix is invertible modulo two, where it is \(I+J\).
The support of \(M\) must be disconnected by Perron
simplicity for a connected nonnegative symmetric matrix.

For internal outdegrees \(d_a\),
\((C^2)_{ab}\equiv -56+2(d_a+d_b)\pmod4\).
This follows by counting the disagreements between the two
sign rows off \(a,b\); their parity is \(d_a+d_b-1\).
If degree parity were mixed, all opposite-parity entries
of \(M\) would be positive and the support would be
connected. So all degrees have equal parity. The signed
row sums, which have one common residue modulo four, now
take either \(-3,1\) or \(-1,3\). There are 14 of magnitude
three and 42 of magnitude one, because the sum is zero.
Thus \(\|C\mathbf1\|^2=168\), giving all-one Rayleigh quotient
52 for \(M\), and its spectral maximum is strictly below
55 by nonsingularity of \(C\).

Let \(I\) be a connected support component attaining that
maximum. Since \(C\) commutes with \(M\), \(C_{II}\) commutes
with \(M_{II}\). The latter's Perron line is one dimensional,
so \(C_{II}\) kills it by skew-symmetry. Hence \(I\) has odd
size and \(C_{II}\) has a one-dimensional rational kernel
by its even principal minors. The Perron eigenvalue is
therefore rational. Because \(M/4\) has integral entries
(all degrees have equal parity), it is a multiple of four.
It lies in \([52,55)\), so equals 52. Equality of the
all-one Rayleigh quotient gives \(M\mathbf1=52\mathbf1\).

Now all support components attain 52 and have odd size.
On their constant subspace \(C\) acts by \(K\) with
\(K^2=-3I\), since \(C^2=M-55I\). For distinct components
\(i,j\) its entries \(k_{ij}\) are signed row sums into
the odd-sized component \(j\), hence nonzero odd integers.
They obey \(|I_i|k_{ij}=-|I_j|k_{ji}\), and the diagonal
entries of \(K\) vanish. Thus
\(\sum_{j\ne i}|k_{ij}k_{ji}|=3\). The number of components
is even and is either two or four. With two, the sizes
are 14 and 42, from the product three and the size ratio.
With four, each product is one, so the sizes are all 14.
Either case contradicts odd sizes. This proves the negative
off-diagonal entry, at most \(-2\) by parity.

For attainment take the induced 56-vertex submatrix of
Paley 59 deleting vertices 56,57,58. Its signed sums have
magnitude at most three, and its square entries off the
diagonal are one plus three sign products, hence at least
\(-2\). The replay checks the full census. Its cross cap is
2874 for any quotient, or 170 for the Paley benchmark.
So the exact full optimum under \(R\le3\) is 2914 for
every quotient. In particular a Paley benchmark attaining
off-branch maximum 160 would require \(R\le3\) by the prior
score screening rule, and cannot attain full optimum 2912.
The true off-branch optimum is not determined here.

## 3. Hamiltonian counts when off-branch balance is 160

The [positive odd-cycle packing formula](https://arxiv.org/abs/2307.05569)
counts Hamiltonian paths as a sum of weights \(2^t\) over
disjoint packs of \(t\) directed odd cycles. Use only triangles.
In an \(n\)-vertex tournament with maximum absolute signed
row sum \(R\), each induced \(s\)-vertex tournament has row
sums bounded by \(b=\min(s-1,R+n-s)\). They sum to zero.
Their squares sum to at most \(s b^2\) for even \(s\) and
\((s-1)b^2\) for odd \(s\); the odd case follows by maximizing
on the zero-sum box, whose extreme points have all but one
coordinate at \(\pm b\) and the remaining coordinate zero.
The exact score formula for triangles then gives
\[
 L_s=\max\!\left(0,\left\lceil
 \frac{s(s^2-1)-3(s-(s\bmod2))b^2}{24}
 \right\rceil\right).
\]
The number of unordered packs of \(t\) disjoint triangles
is at least \(\left\lceil\prod_{j=0}^{t-1}L_{n-3j}/t!\right\rceil\),
by counting ordered choices. Summing with weights \(2^t\)
gives \(L=10089842422501174206273\) for an interior
\((n,R)=(56,3)\); for the Paley quotient \((57,2)\) it gives
\(L_Q=60601093599556238992453\).

**At most half the Hamiltonian paths can start at a specified
vertex with an in-neighbor.** Delete that initial vertex \(u\).
In the resulting path insert \(u\) after the last vertex that
beats it. Such a vertex exists, and any following vertex loses
to \(u\). The new path does not start at \(u\). Deleting \(u\)
recovers the original tail, so this is an injection from paths
starting at \(u\) into paths not starting there.

Along any quotient Hamiltonian path, visit the branches as
contiguous blocks. The previous endpoint rules out just one
initial slot in the next block, the vertex across the reversed
matching. By the injection there are always at least
\(\lceil L/2\rceil\) internal choices for that next branch.
All interior vertices have in-neighbors by near-regularity.
Distinct quotient paths give disjoint families of such paths,
since the branch order is read directly from the full path.
Thus reaching off-branch 160 in the Paley benchmark implies
\[
 H(A)\ge L_Q\,L\,\lceil L/2\rceil^{56}>10^{1260}.
\]
This connects the signed balance tradeoff to the same
disjoint-triangle mechanism that crossed 21, without claiming
that every large odd count is reachable or that the encoding
attains off-branch 160.

Run `node replay.js` for the skew conference and Paley deletion
checks and exact path-count arithmetic. The triangle-packing
lower bound is also checked against direct Hamiltonian counts
on all labelled tournaments of orders 3 through 6.
