# Ternary composition: retain boundary types, ownership, and channels

**Session 77 — mathematical synthesis.** Three-way branching recurs in tournament substitution, arithmetic trees, and SU(3) color tensors. The useful connection is not the number three alone: it is the information each construction must retain so that later composition remains exact.

Two companion notes make this principle concrete:

- [Typed color contractions, alignment memory, and Petersen](COLOR_FUSION.md): finite-dimensional SU(3) representations, exact six-port tensors, and a signed Gram realization of Petersen.
- [Ternary composition remembers more than fifteen Hamiltonian paths](TERNARY_MEMORY.md): the audited path-cover substitution law and an explicit scalar collision exposed by a fixed attachment.

The representation theory and composition theorem have established sources. The finite examples supply explicit checks and comparisons, not a new QCD dynamics, confinement, or mass-gap theorem. No identification of tournament composition with the color-tensor category is asserted.

## 1. Which ternary object?

Several distinct constructions must remain separate.

| Construction | What branches or composes | Data that must survive |
|---|---|---|
| Berggren tree | Three transformations of primitive Pythagorean root pairs | Signed/oriented chart, slope guards, ancestry |
| 3-adic residue tree | Three refinements of each residue class | Compatible residues, carries, integer clock |
| Tournament substitution | Three factors in a directed-cycle quotient | Quotient orientation and factor path-cover profiles |
| Color tensor network | Typed fundamental and dual ports | Representation types, port order, channels, signs |

The [Sep25 arithmetic bridge](https://github.com/eliottcassidy2000/math/blob/02d1218a9a48986564693cae1ae65a050fb9f0e6/05-knowledge/results/ternary_bridge_20260925.md) gives a particularly sharp warning. An inverse-Collatz ray samples a Berggren branch at heights

\[
 H(j)=2^{k_0-1}\frac{4^j-1}{3}.
\]

For a fixed admissible integer k₀≥1, this sampling is sparse among ordinary nonnegative heights, but
\(v_3(H(j)-H(i))=v_3(j-i)\): it induces compatible bijections at every ternary residue depth. Residue completeness does not establish ordinary orbit coverage or Collatz convergence.

There is nevertheless an exact tournament bridge. In [THM-4754](https://github.com/eliottcassidy2000/math/blob/02d1218a9a48986564693cae1ae65a050fb9f0e6/01-canon/theorems/THM-4754-three-adic-sierpinski-tournament-collatz-arc-flip-steiner-midpoint.md), orient \(x\to y\) on \(\mathbb Z/3^k\) when the first nonzero ternary digit of \(y-x\) is 1. The resulting tournament is the iterated substitution \(C_3^{[k]}\). Maps \(x\mapsto3x+b\) give its three child blocks. Units congruent to 1 modulo 3 preserve arcs; those congruent to 2 reverse them. This is a specified map preserving a specified relation—not a shared-number analogy.

## 2. Fifteen paths do not determine the response to attachment

Let \(\operatorname{pc}_T(r)\) count spanning covers of a tournament by \(r\) unordered, vertex-disjoint directed paths, including singletons. Its first entry is the Hamiltonian-path count \(H(T)\).

The inherited, audited [THM-3121](https://github.com/eliottcassidy2000/math/blob/02d1218a9a48986564693cae1ae65a050fb9f0e6/01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md) proves a closed substitution law for the **complete profile**. For quotient adjacency matrix \(A_Q\), put

\[
 W_Q=\mathbf1^TX(I-A_QX)^{-1}\mathbf1.
\]

Then

\[
 \operatorname{pc}_{Q[S_i]}(d)
 =\sum_c\frac{\prod_i c_i!}{d!}[x^c]W_Q^d
          \prod_i\operatorname{pc}_{S_i}(c_i).
\]

The factorials remember ownership: occurrences of a factor label receive distinct actual factor paths. Only after that decoration is division by \(d!\) justified. Repeated bare quotient words do not themselves carry a free permutation action.

The companion exhibits strong five-vertex tournaments \(A,B\) with profiles

\[
 (15,39,33,10,1),\qquad(15,45,35,10,1).
\]

Thus \(H(A)=H(B)=15\), but the same cyclic attachment of two named singleton factors detects their different two-run counts:

\[
 H(C_3[T,K_1,K_1])=3H(T)+2\operatorname{pc}_T(2),
\]

so the outputs are **123 and 135**. This is a compositional obstruction, not a coincidence involving fifteen. The factors have automorphism orders 3 and 5, so it is not a collision of \(H/|\operatorname{Aut}|\). The profile is sufficient for these counting operations; neither global minimality nor efficient arbitrary profile computation is claimed.

## 3. The precise QCD-inspired repair: graded words, not Lie closure

Earlier exploratory physics analogies are not dependencies of these proofs. The [Oct23 Hall/Bott note, §6](https://github.com/eliottcassidy2000/math/blob/02d1218a9a48986564693cae1ae65a050fb9f0e6/05-knowledge/results/metagraph_hall_bott_20261023/README.md) already extracts a useful exact mechanism and identifies its loss.

An arc \(i\to j\) maps to \(e_i-e_j\). Summing all such roots retains score deviations but loses cycle placement. Keeping matrix units is richer:

\[
 [E_{ij},E_{jk}]=E_{ik}\qquad(i,j,k\text{ distinct}).
\]

Yet **every strongly connected tournament on \(n\) vertices generates \(\mathfrak{sl}_n(\mathbb C)\)**. A simple directed path produces its endpoint matrix unit by iterated commutators. Strong connectivity supplies both directions; their brackets produce all traceless diagonal matrices. Consequently, the ungraded generated algebra cannot distinguish strong tournaments of the same size.

The repair is to retain grades, distinct-vertex restrictions, actual endpoints, and shared-port ownership. These are exactly the kinds of information that quotient-word decoration protects in the path-cover theorem. The relevant algebra here is complex \(\mathfrak{sl}_n\), not an assertion that directed matrix units lie in compact \(\mathfrak{su}(n)\).

## 4. A genuine A2 connection—and a necessary type check

SU(3) has **rank two**, fundamental dimension three, and adjoint dimension eight. Its A2 weight/root quotient has order three. Prior [THM-2774](https://github.com/eliottcassidy2000/math/blob/02d1218a9a48986564693cae1ae65a050fb9f0e6/01-canon/theorems/THM-2774-tree-path-smith-index-ladder-and-binary-ternary-lattice-defects.md) supplies a marked integral version: the frame

\[
 \begin{pmatrix}1&-1&0\\0&1&-1\\1&1&1\end{pmatrix}
\]

has Smith form \(\operatorname{diag}(1,1,3)\), and, among integral Euclidean isometries of the ambient coordinate lattice, the stabilizer of the marked diagonal is \(S_3=W(A_2)\). This is an integral index and symmetry connection, not merely a real two-dimensional picture.

The color calculation supplies another exact bridge. In
\(V^{\otimes3}=1\oplus8\oplus8\oplus10\), the two octet channels have a centered A2 multiplicity plane. Cyclic permutation acts there, and \((C-C^{-1})^2=-3\), recovering the triangular/Eisenstein operator. Re-pairing mixes the two channels: retaining the label “8” alone loses information.

But a rooted three-child color operation \(V^{\otimes3}\to V\) is impossible: the center acts trivially on the source and nontrivially on the target. Epsilon contracts three fundamental ports to a **scalar**, not another fundamental port. A three-child rooted node has four incident ports; an epsilon vertex has three. Typed webs, rather than untyped ternary trees, provide the composition-closed framework.

## 5. Six ports, five dimensions, and Petersen

For six named fundamental ports, take the ten tensors obtained by splitting the ports into two unordered epsilon triples. Their Gram matrix \(G\) has diagonal 36, off-diagonal entries \(\pm12\), and

\[
 G^2=72G,\qquad \operatorname{rank}G=5.
\]

They span the full five-dimensional invariant space. Thus ten descriptions are neither ten independent channels nor ten interchangeable labels. A signed four-term relation is pointwise exact; dropping signs destroys it. The norm of a linear combination is \(a^*Ga\), not an independent sum of coefficient norms.

After an explicit sign gauge, the associated Seidel matrix yields the Petersen graph, equivalently \(KG(5,2)\): vertices correspond to edges of \(K_5\), adjacent when disjoint. This is a verified Gram/switching construction, **not** a Petersen-to-\(K_5\) minor contraction or a physical graph prediction. Its frame-theoretic setting is classical; conventions and exact equations are given in [COLOR_FUSION.md](COLOR_FUSION.md).

The common lesson is now operational: **a summary is adequate only when every permitted later composition can be computed from it**. Tournament counts need run ownership; color contractions need boundary types, multiplicity channels, and signed pairings. The explicit 123/135 and six-port examples show where smaller summaries fail, while leaving the two mathematical composition laws distinct.
