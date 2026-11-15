# The exact tape of a uniform tournament module

**PROVED (all orders).** The complete one-hole continuation state for Hamiltonian-path counting is the path-cover vector weighted by two factorials. Fixed-modulus arithmetic admits a finite state machine, with an explicit smallest reader for any first differing grade. Absolute arithmetic admits neither a fixed finite-state readout nor a fixed-dimensional linear readout, even for tapes whose scalar path count is always one.

## 1. Two independent permutation factors

Let \(D\) be a nonempty tournament and \(U_j(D)\) count *unordered* covers of its vertices by \(j\) nonempty directed paths (each path carries its direction and vertex order). Put
\[
F_j=j!(j-1)!,\qquad Z_j(D)=F_jU_j(D).
\]
A reader \(R\) is a tournament with one marked vertex \(v\); write \(R[D]\) for substitution, where every outside vertex's orientation towards all of \(D\) is inherited from its orientation towards \(v\).

Every Hamiltonian path in \(R[D]\) splits into \(j\) maximal visits to \(D\). For each unordered internal cover, the \(j!\) orderings have identical validity once the outside pieces have been placed. The outside pieces are a possibly empty prefix and suffix and \(j-1\) nonempty interior pieces. Each interior piece is a directed path starting at a vertex beaten by \(D\) and ending at a vertex beating \(D\). Permuting the interior pieces freely is another independent symmetry: no pair of adjacent outside pieces is joined directly. The action has size \((j-1)!\) because the pieces are nonempty and vertex-disjoint. Thus
\[
H(R[D])=\sum_{j\ge1}a_{R,j} Z_j(D),\qquad a_{R,j}\in\mathbb Z_{\ge0}.
\]
The coefficient depends only on \(R\). Moreover every interior piece needs at least two vertices: a single outside vertex cannot both beat and be beaten by the module. Consequently a reader with \(m\) outside vertices has support only on \(j\le1+\lfloor m/2\rfloor\). The second factorial is a reader symmetry, distinct from merely making the internal paths unordered.

## 2. An optimal physical probe for each tape coordinate

For grade \(k\), build a test \(B_k\) with \(k-1\) outsiders of each type \(b_i,c_i\), the relations \(D\to b_i\), \(c_i\to D\), and
\[
b_i\to c_\ell\quad\Longleftrightarrow\quad i\le\ell.
\]
Within each type choose a transitive order of increasing indices. At \(k\) visits, the outside budget is exhausted by \(k-1\) two-vertex interior pieces, so the prefix and suffix are empty. A perfect matching of the pairs must be the diagonal: its permutation \(\pi\) satisfies \(\pi(i)\ge i\) for every \(i\), and their sums are identical. The pieces may be permuted freely, but there is exactly one orbit. Thus \(a_{B_k,k}=1\), and
\[
H(B_k[D])=Z_k(D)+\sum_{j<k}a_{B_k,j}Z_j(D).
\]
This is an integer unitriangular decoder. It works without division over every \(\mathbb Z/M\mathbb Z\). Calibrated rows at \(k=1,2,3,4\) are \((1),(3,1),(9,10,1),(31,83,25,1)\).

**Exact contextual equivalence.** Two nonempty modules \(D,E\) have identical path counts in every one-hole reader if and only if \(Z_j(D)=Z_j(E)\) for every \(j\). Over a fixed modulus \(M\) the same equivalence holds with each equality replaced by congruence. If their first differing grade is \(k\) (exactly or modulo \(M\), respectively), the least number of outside vertices needed to distinguish them is exactly \(2k-2\). Both statements follow from the support bound and the attaining triangular probe.

This is not an isomorphism fingerprint, and does not govern a context with several exposed ports from the same module or with nonuniform contacts.

## 3. What the finite machine really knows

Let \(K(M)=\min\{j\ge1:M\mid F_j\}\) for \(M\ge2\). Because \(F_j\mid F_{j+1}\), the first \(K(M)-1\) coordinates suffice for **every** one-hole reader, with no order bound. The possible residue vectors number at most \(\prod_{j<K(M)}M/\gcd(M,F_j)\); not all vectors need be realized by physical modules. For \(M=p^a\), Legendre's formula applied to the two factorials gives
\[
K(p^a)=\tfrac12(p-1)a+O_p(\log(a+1)).
\]
For example, only grades 1 and 2 matter modulo \(3\), while the first grade with \(256\mid F_j\) is \(7\). For any finite repertoire of one-hole instructions \(D\mapsto R[D]\), the modular vector is a deterministic finite state. Indeed, executing a probe after an instruction is again a one-hole reader: expand the two outer graphs first and leave the innermost port unfilled. Its state transitions are explicitly linear modulo \(M\), by expanding these composed readers and applying the unitriangular decoder. No approximation or unproved continuation assumption is present.

In contrast, the absolute table with columns \(B_k\) and rows \(TT_n\) has unbounded rational rank. For these tapes \(U_j(TT_n)=S(n,j)\), because a set partition into \(j\) blocks determines the unique directed path in each block. For \(1\le n,k\le r\), the evaluation table therefore factors as a unitriangular Stirling matrix, a diagonal with entries \(F_j\), and a unitriangular reader matrix (with the appropriate transpose). Its determinant is \(\prod_{j=1}^r F_j\ne0\). Every \(TT_n\) has \(H=1\), but \(B_2[TT_n]\) has \(H=2^n+1\), so these physical tapes also require infinitely many absolute equivalence states even without a linearity assumption.

One can replace \(TT_n\) by \(C_3\Rightarrow TT_n\) if a nontransitive tape is desired: the common immediate count is now \(3\), whereas \(B_2\) reports different powers of two with nonzero coefficient. Neither this nor the rank argument asserts a bounded-state obstruction for an entirely different selective construction.

## 4. Bridge from the old Rule 30 and Jacobian threads

[THM-3456](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary.md) separates an invertible, triangular free-boundary trace compiler from selecting the actual boundary of the distinguished single seed. The reader family above gives an exact tournament counterpart for *decoding*: its unique diagonal matching plays the role of the permutive corner. It does **not** invert the physical source problem. Coordinates of a module are simultaneous positive path-cover counts of a common tournament, not freely programmable cellular-automaton inputs.

[THM-3463](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-3463-rule30-mealy-section-suffix-parity-current-and-complexity-boundary.md) has an exact finite Mealy codec on a projected observer, alongside a separate growing-input complexity bound. The modular versus absolute split above gives a corresponding distinction with a literal tournament machine and an independent determinant test. More unexpectedly, the factorial truncation in [the old planar-Jacobian/Gaussian-moment Hankel note THM-2810](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go.md) is structurally parallel: full rank in characteristic zero, sharp observer loss in positive characteristic. Here the two factorials come from two independent path orders, not from Gaussian moments.

## 5. Consequence for a Chan–Pak compiler

For a proposed one-port repair, there is now a finite, universal modular test for continuation preservation. At the same time, fixed-modulus equality cannot certify any absolute lift. In particular, the strategy “alter the module but preserve its continuation against every reader exactly, then read off a changed count” is inconsistent: the replacing modules have identical \(Z\), hence there is no changed count in this class. A productive compiler must preserve only the reads that its next branch actually uses, or change the exposed boundary/reader itself. A sharp intermediate target is to construct physically different modules sharing several exact low grades but separating at a specified later grade; the minimal-size probe above tells exactly when such a delayed difference can first be read. No general physical synthesis of these profiles or Chan–Pak proof is claimed.

Reproduction: **python ca_bridge_state_check.py** calibrates the four matrix rows using transitive inputs, verifies each diagonal is one, then independently compares their predictions with subset-DP Hamiltonian counts for all 75 labeled inputs of orders 1–4 (300 comparisons).

## 6. An exact FRACTRAN-style guard which cannot decrement at one uniform port

Put \(g(D)=\gcd_{j\ge1} Z_j(D)\), a positive integer because \(Z_1=H(D)>0\). Every uniform readout is divisible by \(g(D)\). Conversely, if a number divides every physical readout, it divides each coordinate by the integral decoder. Therefore
\[
g(D)=\gcd_R H(R[D]).
\]
In particular \(g(D)\mid g(R[D])\): any reader of the modified tournament is a composite reader of the original. This invariant is useful for a FRACTRAN analogy because it identifies a genuine, not just desired, divisibility guard carried by the physical tape. An instruction repeatedly nesting uniform readers cannot decrement any prime valuation of \(g\). That says nothing about the individual valuation of the current displayed \(H\) above the guaranteed guard.

**The robust factor of three.** Let \(Q_D(t)\) be the optional ordered slot polynomial; its parity \(Q_D(-t)=(-1)^{|D|}Q_D(t)\) follows from the Rédei–Berge odd-cycle packing expansion (see [Frozen grades](https://github.com/eliottcassidy2000/math/blob/a2e72288d7e083245263444b8ff62458e517ecfe/05-knowledge/results/chanpak_two_three_archive_20261115/FROZEN_GRADES.md)). Since \(Q_D(2)=2H(D)+Z_2(D)\) and \(2=-1\bmod 3\), modulo 3 the state is \((Z_1,Z_2)=(H,0)\) for odd order and \((H,-H)\) for even order. All later factorials are zero. The one- or two-vertex transitive module has \(H=1\) and the corresponding parity. Thus for every reader
\[
H(R[D])\equiv H(D)H(R[TT_\epsilon])\pmod 3,\qquad
\epsilon=\begin{cases}1& |D|\text{ odd},\\2& |D|\text{ even}.\end{cases}
\]
In particular, \(3\mid H(D)\) if and only if \(3\mid g(D)\); the zero state is absorbing even under adaptive choices from *all* uniform readers. No automorphism hypothesis is required. From a divisible-by-3 core one needs a genuinely nonuniform contact to reach a target not divisible by 3. This narrows Chan–Pak strategies that attempt to move from a selected \(3\)-divisible seed back to units by switching only uniform continuations.

Higher powers do not simply inherit this property from \(H\). On five vertices, with pair order \((0,1),(0,2),\ldots,(3,4)\) and a 1 bit meaning lower-to-higher, the module of code 8 has \(H=9,\ Z_2=66\), so the cyclic 3-vertex reader \(B_2\) returns \(3\cdot9+66=93\) and \(g=3\). The one-port operation lowered the displayed 3-adic valuation, but did not lower the robust valuation of its continuation guard.
