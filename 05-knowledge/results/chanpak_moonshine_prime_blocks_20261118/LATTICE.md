# A literal Gaussian lattice for tournament switching

This addendum turns the lattice analogy into an exact construction. It recovers the switching class and its projective symmetry group, but **not** the Hamiltonian count without a frame. This distinction is directly relevant to choosing a physical block at a prescribed prime.

Write $A_{uv}=1$ for the arc $u\to v$, $A_{vu}=-1$, and $A_{uu}=0$. On $\mathbf Z[i]^n$ use the Hermitian Gram matrix
$$
K_A=(2n-1)I+iA.
$$
Indeed $z^*K_Az\ge n\|z\|^2$: the sum of the cross-term absolute bounds is at most $(n-1)\|z\|^2$. The shortest nonzero vectors are exactly unit multiples of the coordinate vectors, because they have norm $2n-1$, while any vector with squared coordinate norm at least 2 has lattice norm $\ge2n$. Any **Gaussian-linear** isometry is therefore monomial with Gaussian-unit entries. Comparing a single off-diagonal entry forces the ratio of the unit phases at its two ends to be real, hence $\pm1$. Up to a common Gaussian phase the isometry is exactly a signed permutation. The Gaussian-linear isometry class is thus equivalent to the switching class up to relabeling. Choosing a shortest coordinate frame with purely imaginary off-diagonal entries, and retaining the signs, recovers the tournament signs. Isometry classes here do not refer to all isometries of the underlying real lattice.

Let $G$ be the switching stabilizer, i.e. the permutations taking the signs into their switching class. Its lift to real signed permutations has central kernel $\{\pm I\}$; the Gaussian-linear isometry group modulo common Gaussian phases is $G$. For a nonidentity involution $g\in G$, a swapped pair forces the corresponding two signs to be opposite, since the permutation reverses the ordered pair in a skew matrix. Its lift squares to $-I$ (the square is in the kernel, so testing the swapped pair determines it). In particular an involution cannot fix a vertex. Two distinct commuting involutions therefore lift to anticommuting matrices, because their product also lifts to something whose square is $-I$. Three independent commuting involutions are impossible: the first lift commutes with the product of the other two lifts, contradicting the relation for the first involution and the nontrivial product. This short proof of the rank-two bound is a special case of the stronger [Babai–Cameron classification](https://doi.org/10.37236%2F1516).

The bound is attained by the four-vertex $3+1$ frame, a directed triangle dominating a singleton. Its switching stabilizer is $A_4$. In mask order $(10,12,9,0)$ and zero-based labels, the double transpositions $(01)(23),(02)(13),(03)(12)$ lift using cut masks 6, 3, 5 respectively. The three-cycle on the triangle needs no switch. Checking the 24 permutations gives exactly the 12 even ones. The inverse image of the Klein four is $Q_8$, and the rotation acts nontrivially with order three on it. The real signed lift is the binary tetrahedral group $Q_8\rtimes C_3$, or $\mathrm{SL}_2(3)$. This is an actual central double cover on tournament data, not an identification with the Monster's $2A$ centralizer $2.B$.

## What is conserved, and why it does not choose a prime

Fix labels and sum over the $2^{n-1}$ switch frames. Put $C_k$ for the ordered covers by $k$ **nonempty** paths. For a fixed permutation of the vertices and placement of $k-1$ cuts, the required arcs form a forest with $k$ components. Each component has one free switch sign, because along its edges specifying one sign determines all the others. Modding out the common sign leaves $2^{k-1}$ compatible frames. Consequently
$$
\sum_{\text{frames}} C_k = 2^{k-1}n!\binom{n-1}{k-1},
\quad
\sum_{\text{frames}}h=n!,
\quad
\sum_{\text{frames}}B=2n\,n!.
$$
These formulae count labeled frames, not relabeling orbits; the last equality also holds at $n=1$ with $C_2=0$. For several prescribed covers in the same frame, use the union of their forests. An edge prescribes the product of the two endpoint signs. If prescribed products disagree on an edge, or a cycle has product $-1$, the joint contribution is zero. Otherwise it is $2^{c-1}$, where $c$ counts all components including isolates. This supplies exact moment constraints on a profile bank; knowing an average does not suffice to select a frame.

For example, the transitive masks $(14,12,8,0)$ and the cut at vertex 1 giving $(12,1,10,2)$ have isometric Gaussian lattices. Their $(h,B)$ profiles are $(1,16)$ and $(5,32)$, respectively. Over all eight labeled frames, four give each profile, and the sums are 24 and 192 as predicted. The information forgotten by the lattice is already enough to change the answer from 1 to a prime. For Chan–Pak, a useful lattice argument must therefore either produce the marked shortest frame from the target or prove its existence via correlations, not just select the unmarked lattice class.

The Gaussian Gram construction and double-cover computation are independent of the finite source catalogue used in the main note. They were recovered in the background research pass on this question. The following short literal replay verifies the order-four computations (run in this directory alongside replay.py):

~~~python
from itertools import permutations
from replay import profile

def cut(R, s):
    return tuple(r ^ ((15 ^ s) if s >> i & 1 else s)
                 for i, r in enumerate(R))

def relabel(R, p):
    return tuple(sum((R[p[i]] >> p[j] & 1) << j for j in range(4))
                 for i in range(4))

T, Q = (14,12,8,0), (10,12,9,0)
assert cut(T,2) == (12,1,10,2)
assert sorted(profile(cut(T,s)) for s in range(8)) == [(1,16)]*4+[(5,32)]*4
for p in permutations(range(4)):
    choices = [s for s in range(8) if cut(Q,s)==relabel(Q,p)]
    inversions = sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))
    assert len(choices)==1-(inversions%2)
assert [cut(Q,s)==relabel(Q,p) for s,p in [(6,(1,0,3,2)),(3,(2,3,0,1)),(5,(3,2,1,0))]] == [True]*3
~~~
