# Seven recovered mechanisms and four honest bridges to hard problems

Session 89, 2026-10-28. Repository source pin: `03408f84ac43fcdcfb0887e5d6bc1fc8792415cf`. The new contribution in this lane is an exact repair of a historical codec and a typed bridge/obstruction ledger, not a claim to solve isomorphism, reconstruction, or polynomial identity testing. “Inherited” means the scoped statement in the named repository source; exploratory output is not silently promoted to theorem status.

## 1. Consecutive triangles really tile a square — after fixing the old map

**Old sources:** March 24 `04-computation/tiling_square_s325.py`, its output, and `04-computation/new_sequences_s327.py`. **Status:** historical coordinate proposal; new elementary correction proved in [REPAIRED_SQUARE_CODEC.md](REPAIRED_SQUARE_CODEC.md).

- **Source / target / map:** fixed-path tilings of sizes `m+2,m+1` map to `m×m` matrices by `a_ij -> M[j-2,i]` and `b_ij -> M[i,j-1]`.
- **Preserved:** every bit, Hamming distance, coordinate flips, uniform measure. Via THM-830 these tilings are free labelled *gap* tournaments of orders `m+1,m`.
- **Lost:** nothing at the labelled level; natural matrix symmetry and natural tournament relabeling have not been identified.
- **Correction:** the old upper map uses `j-2`, discards `m-1` bits, and has image dimension `m²-m+1`. Its own side-three output shows only seven filled cells. Its rank census is therefore not a census of all binary square matrices.

This supplies the closest genuine predecessor to the user's observation, rather than claiming the triangular-number identity was just discovered.

## 2. Observer coordinates carry an affine relabeling cocycle

**Old source:** `05-knowledge/results/tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md`, §1. **Status:** inherited elementary proof.

- **Source / target / map:** a labelled tournament `T` with fixed linear observer `pi` maps to the simple graph `G_pi(T)` of edges reversed relative to `pi`.
- **Preserved:** every labelled arc; converse becomes graph complement, arc reversal becomes XOR, and switching becomes addition of a cut.
- **Relabeling law:** `G(gT)=gG(T) xor I(reference,g reference)`. Observer changes satisfy `I(pi,sigma)=I(pi,rho) xor I(rho,sigma)`.
- **Lost by forgetting the observer:** the order-dependent affine term. Ordinary graph isomorphism is not tournament isomorphism under this codec. The empty graph is fixed by all permutations, but no tournament is fixed by a transposition.

For paired-tournament encodings this is the exact warning: conjugating an action through a bit codec gives a valid action, but not necessarily the visually natural matrix action. The action, not only the array, is part of the instance.

## 3. Tree normalization exposes both cycle syndrome and gauge multiplicity

**Old sources:** the same trijection packet §§2–3; THM-1415 and THM-1467 are its credited predecessors. **Status:** inherited all-size proof.

- **Source / target / map:** switching classes `x+B` of tournament bits map to representatives with prescribed bits `a` on a constraint graph `F`.
- **Preserved:** the switching class, provided the normalization succeeds. If `F` has `c` components and cycle rank `beta`, success is equivalent to zero cycle syndrome of `x+a`; every nonempty fibre has `2^(c-1)` representatives. Coverage is a fraction `2^(-beta)` of the switching classes.
- **Lost:** forgetting the original switch loses `n-1` bits; forgetting the cycle syndrome incorrectly asserts existence; forgetting component phase incorrectly asserts uniqueness.
- **Exact trijection:** a spanning tree is precisely the case `c=1,beta=0`, giving a transversal. Different tree charts have invertible linear transition maps obeying a cocycle law.

This is a reusable certificate mechanism: a small signed cycle can certify inconsistency. It is not a general canonical-label algorithm, because switching and isomorphism are different equivalence relations.

## 4. Paired faces must retain common overlap and deletion seams

**Old sources:** THM-796, §§1–3; THM-830, §§1–4. **Status:** inherited all-size coordinate and gluing proofs; their finite non-lumpability tests retain their finite status.

- **Source / target / map:** a fixed-path tiling `X_n` maps to its two endpoint-deleted faces and apex bit:
  `X_n ≅ (X_(n-1) ×_(X_(n-2)) X_(n-1)) × F_2`.
- **Preserved:** every arc through the actual common interior restriction and apex. The gap map `t_(a,b) -> G_(a-1,b)` gives a lossless labelled tournament of order `n-1`.
- **Lost:** replacing face objects by separately sorted isomorphism classes loses their common labeling. Passing to complement lines loses a coherent `C_2` endpoint phase for `n>=5`. Internal vertex deletion also needs the newly adjacent shortcut bit; bare gap contraction discards that seam.
- **Minimal hostile:** in `0→1→2→0`, deleting vertex 1 does not leave the putative forward observer path `0→2`. Deletion of a marked-path tournament cannot simply forget a row without checking or repairing the seam.

This is the right precedent for reconstruction decks: retain coupled cards and overlap maps, not just independently canonicalized marginals. It does not establish unlabelled graph-deck reconstruction.

## 5. Signed triangle curls recover a tournament only with the cut coordinate

**Old source:** `05-knowledge/results/metagraph_reconstruction_20261023.md`, §1 and linked port-rigidity proof. **Status:** inherited all-order statement, with finite wiring experiments distinguished there.

- **Source / target / map:** signed adjacency `a_ij∈{±1}` maps to signed curls `c(i,j,k)=a_ij+a_jk+a_ki` and one distinguished vertex's outdegree.
- **Preserved:** all arcs are reconstructible when the score coordinate is retained. Without it, equal curl tensors identify cyclic rotations of the ordered strong components. A fixed Hamiltonian path also supplies the missing cut.
- **Lost:** triangle counts or unlabeled multisets forget labels and signs. Even the complete labelled curl tensor alone forgets the cut of the component necklace.
- **Mechanism:** two equal-curl arrays differ by a gradient. Since differences are in `{−2,0,2}`, any nonconstant potential has two levels and reverses a uniformly directed cut. Those cuts occur between strong components.

This creates a promising *specific* target for the pair-matrix corpus: determine exactly which row/column observations retain the corresponding cut and coupling. It does not say a handful of counts is a complete invariant.

## 6. Transfer-matrix symmetry is a constrained polynomial identity

**Old source:** March 6 THM-030, `THM-030-transfer-matrix-symmetry.md`. **Status:** inherited induction proof, with symbolic finite checks stated separately.

- **Source / target / map:** weighted tournament arcs `t_ij=r+s_ij`, `s_ij=−s_ji`, map to the signed subset sum of Hamiltonian-path endpoint products defining `M_W[a,b]`.
- **Preserved:** the constant pair-sum relation and path combinatorics imply `M_W[a,b]=M_W[b,a]`; the proof shows only even powers of `r` occur.
- **Lost:** the result does not assert injectivity of this transfer observable, nor does it apply to arbitrary directed adjacency arrays or arbitrary square-bit encodings.

Thus “the matrix is symmetric” can mean (i) raw transpose symmetry, (ii) symmetry of an induced transfer matrix, or (iii) a permutation automorphism. These are not interchangeable. A useful PIT extension would specify the exact polynomial and quotient relations before importing a symmetry identity.

## 7. Numerical fingerprints and cycle indices are filters, not canonical labels

**Old sources:** `07-reflections/tournament-matrix-classic-results-atlas-codex-s210.md`, items 51–52; `07-reflections/the-quarter-tiling-model-and-the-reconstruction-wall-at-n7-kps.md`; `05-knowledge/results/session80_tournament_symmetry_20261027/SYMMETRY_INFORMATION.md`. **Status:** atlas is a research dictionary; the n=7 reflection is explicitly a finite sampled obstruction; the later cycle-index packet distinguishes proved information-loss laws from finite certified examples.

- **Source / target / map:** adjacency maps to characteristic/OCF polynomials, determinants, or a permutation-group cycle index.
- **Preserved:** the chosen scalar/spectral/orbit-count observations. Matrix conjugacy by permutation is the actual isomorphism relation, not equality of those observations.
- **Lost:** the sampled n=7 note records collisions of `(OCF,determinant)` and skew spectrum. The cycle-index packet gives actual prime tournaments at order 81 with equal cycle indices but different binary-decoration stabilizer counts: 80 versus 26 orbits of stabilizer order 27.
- **Lesson for pairs:** separately canonicalizing two pieces loses the relative alignment; retaining a cycle index loses subgroup-level simultaneous fixed data. A coupled stabilizer intersection requires actual actions, not only group orders or spectra.

No all-size lower bound on every possible invariant is inferred from a finite collision. The proper claim is that these named summaries are insufficient in the exhibited scope.

# Four bridges with explicit complexity boundaries

## A. A transported pair action can retain general graph isomorphism — not remove it

Let `F` be any explicit polynomial-time square-bit-matrix/pair codec. Define pair equivalence by

```
F(M) ~ F(N)  iff  N=P M P^T for a permutation matrix P.
```

**Source / target / map:** binary square matrices under simultaneous row/column relabeling map to paired tournaments carrying the **transported** action `F(PMP^T)`. **Preserved:** isomorphism and stabilizers exactly. **Lost:** the target action need not be natural tournament relabeling.

This problem is GI-complete under polynomial-time reductions. GI reduces to it by taking symmetric zero-diagonal adjacency matrices. Conversely encode a matrix as a four-colored simple graph with vertices `u_i,r_i,c_i`, and `w_ij` for each one entry; add edges `u_i-r_i`, `u_i-c_i`, `r_i-w_ij`, `c_j-w_ij`. Colors are `U,R,C,W`. Any isomorphism uses the same permutation on the three indexed sorts, and preserves exactly the one entries. Colors can be removed by attaching one distinguished hub for each color, with different sufficiently large numbers of pendant leaves; the original graph has no isolated vertices, so leaves and original vertices remain distinguishable. The construction has polynomial size.

This is an elementary reduction/certification of the *transported* problem, not a new TI hardness theorem. For the natural diagonal action on a tournament of order `n` and a rooted tournament of order `n+1`, stabilizers are subgroups of the first tournament's odd-order automorphism group. The zero matrix, by contrast, has stabilizer `S_n`. For `n>=2`, no stabilizer-preserving equivariant bijection between these natural actions exists. The parent packet's n=2 orbit comparison (matrix 10, natural rooted pair 8) is another finite obstruction; the checker here independently verifies the matrix count 10.

**Current primary-source status:** Babai proves general GI quasipolynomial. Grohe–Neuen's February 2026 paper still reports the unrestricted TI bound `n^{O(log n)}`, while proving `k^{O(log k)}n^{O(1)}` for twin width `k`. It also supplies bounded-twin-width tournaments requiring linear-order WL dimension, so fixed-dimensional refinement cannot be promoted to a universal TI proof. [Babai](https://arxiv.org/abs/1512.03547), [Grohe–Neuen](https://theoretics.episciences.org/17560/pdf). A newer September 2026 primary preprint proves FPT isomorphism for `F`-free tournaments parameterized by the forbidden tournament size, hence polynomial-time isomorphism on every fixed nontrivial hereditary tournament class; this is still not an unrestricted TI algorithm. [Neuen](https://arxiv.org/abs/2609.24804).

## B. Intertwiner determinant is a GL-similarity certificate, not a permutation certificate

For same-sized matrix tuples `A=(A_i), B=(B_i)` over a field `K`, compute

```
L(A,B)={X : X A_i=B_i X for every i}.
```

This is a linear space. If `E_1,...,E_d` is a basis, its generic determinant is
`D(z)=det(sum_j z_j E_j)`.

**Source / target / map:** a tuple-comparison instance maps by linear equations to an intertwiner space and determinant polynomial. **Preserved:** an invertible intertwiner is exactly simultaneous similarity; if `K` is infinite, its existence is equivalent to `D` not being the zero polynomial. Over finite `K`, actual `K`-rational nonvanishing must be distinguished from formal nonzeroness for arbitrary matrix spaces; passing to a sufficiently large extension is the usual PIT setting. **Lost:** the coordinate basis/permutation condition.

The tiny binary hostile is

```
A=[[1,0],[0,0]], B=[[1,1],[0,0]], X=[[1,-1],[0,1]].
```

Here `XA=BX` and `det X=1`, but no permutation conjugates the diagonal `A` to the nondiagonal `B`. This already defeats replacing permutation intertwiners by unrestricted invertible matrices. For synchronized same-sized pairs, require the same unknown `X` for both components; testing each component separately loses even more information. For consecutive-size rooted tournament pairs, first specify a common-space padding and the root constraints rather than silently treating two different matrix sizes as one tuple.

Chistov–Ivanyos–Karpinski give deterministic polynomial-time algorithms for this structured simultaneous-similarity/module-isomorphism problem over finite and algebraic number fields. Their paper explicitly separates it from finding a nonsingular matrix in an **arbitrary** linear matrix space. Quiver representation classification can be wild without this particular isomorphism decision problem thereby being NP-hard or open. [Primary paper, introduction and Theorem 2](https://theory.cs.uni-bonn.de/marek/publications/85169-cs.pdf).

**Concrete research target, not theorem:** identify a sharply defined paired-tournament family for which an invariantly specified extra condition forces every surviving intertwiner space to contain a permutation, or exhibit systematic hostiles to that forcing. Solving only GL similarity does not meet that target.

## C. Dyadic monomial separation gives a certificate cost, not general PIT derandomization

A coordinate codec does not lower the symbolic degree or preserve cancellations unless that is proved. For a multilinear polynomial in `d` variables, the substitution

```
x_i=t^(2^i), i=0,...,d-1
```

is injective on monomials because binary subset sums are distinct. Its resulting degree can be `2^d−1`. Without the multilinearity promise it is not injective at all: `x_1−x_0²` becomes zero. For individual degree at most `D`, base `D+1` instead gives unique exponent digits, with potentially exponential univariate degree. Fixed numeric evaluation can also cancel, e.g. `x_1−2x_0` at `(2,4)`.

**Source / target / map:** bounded-individual-degree polynomial coefficients map to univariate exponent addresses. **Preserved:** formal coefficients only under the stated digit bound. **Lost:** polynomial-size degree/explicit-evaluation cost; without the bound, even formal identity information. A brute univariate degree-plus-one evaluation test is not polynomial in `d` here. This is an obstruction to that proposed route, not a lower bound against every algorithm.

For a symbolic rank question, use the relevant minors or generic determinant and keep field, degree, circuit size, and evaluation model explicit. A numerical rank sample certifies a lower bound on generic rank, not exact generic rank by itself. Andrews–Kush–Tell's 2025 result gives polynomial-time PIT **under specified circuit lower-bound assumptions**, not an unconditional solution. Kabanets–Impagliazzo prove that deterministic polynomial-time PIT over integer arithmetic circuits would imply a major circuit lower-bound disjunction. [Andrews–Kush–Tell](https://arxiv.org/abs/2504.06044), [Kabanets–Impagliazzo](https://www2.cs.sfu.ca/~kabanets/papers/poly_derand.pdf).

The bounded archive searches did **not** recover an old repository theorem asserting a general dyadic PIT algorithm. None is attributed to the corpus here.

## D. Counting and reconstruction transfer only with their actual observables

**Counting source / target / map:** define `PerPair(F(M))=per(M)` by decoding. This is a polynomial-time bijective re-encoding of the zero-one permanent problem, so its exact integer counting hardness is inherited. It is not a reduction to counting tournament Hamiltonian paths, and says nothing about such a reduction without a construction preserving the counted witnesses. This decoded permanent is a function of the labelled codec, not automatically of natural pair isomorphism classes: at n=2 the zero matrix and its natural paired transposition image [[0,1],[1,0]] represent isomorphic natural pairs but have permanents0 and1. Determinant counts signed permutation terms; permanent counts unsigned ones. For the all-ones 2×2 matrix these values are 0 and 2; modulo two they coincide, losing the integer count. Valiant's theorem establishes #P-completeness for zero-one permanent under its counting reductions. [Valiant's original paper](https://sites.math.washington.edu/~billey/colombia/references/valiant.permanent.1979pdf.pdf).

**Deck source / target / map:** a labelled matrix deletion operator can always be transported as `F_(n-1) d_v F_n^{-1}`. This preserves the *defined* labelled deletion operation. It does not make that operation equal to deleting one naturally named tournament vertex in each component; the tiling seam and overlap examples above explain why that equality requires proof. Nor does it identify unlabelled deck equivalence under two different group actions.

This distinction is indispensable: ordinary tournament reconstruction from vertex-deleted cards is false in infinitely many orders by Stockmeyer; ordinary simple-graph reconstruction remains a separate conjecture. Recent primary work proves semi-random graph results, not the unrestricted conjecture. Thus paired tournament decks are a prospective information-rich research carrier, not a shortcut through an already false tournament reconstruction principle. [Stockmeyer](https://onlinelibrary.wiley.com/doi/abs/10.1002/jgt.3190010108), [Asilis–Chen–Hansen–Teng, ITCS 2026](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ITCS.2026.12).

# Verification and publication scope

`verify_history_bridges.py` passes 79,312 always-active checks in ordinary and optimized Python, with identical receipts. It independently checks the repaired square map and the historical fibre defect, Hamming cube edges, the explicit matrix-to-graph reduction through order three, the GL-versus-permutation hostile, the matrix n=2 orbit count, and elementary dyadic/counting hostiles. It does not reprove every historical theorem or claim to exhaust all relevant history.

The live codec repairs do not validate or regenerate the old rank/determinant census, and no unrelated historical claim is promoted by this packet.

## Late recovered connection and concurrent arithmetic instance

[SELF_CONVERSE_PRIOR_STATUS.md](SELF_CONVERSE_PRIOR_STATUS.md) pins the exact S22 boundary between a count equality and a structural lift, and gives an independently reproved odd-order self-converse TI reduction with no novelty claim. [concurrent/BANK_MATRIX_BOUNDARY.md](concurrent/BANK_MATRIX_BOUNDARY.md) independently replays the supplied nine-prime/F64 bank and identifies its square-sigma subsets with a five-bit `3×3` rectangle code. A natural pair relabeling sends the empty subset to `{59,7}` and changes square sigma to nonsquare sigma, giving a literal arithmetic/action boundary.
