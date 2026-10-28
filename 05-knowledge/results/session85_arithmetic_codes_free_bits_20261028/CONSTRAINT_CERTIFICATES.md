# A bit without another measurement: constraint certificates, codes, and legality

Date:2026-10-28. **PROVED elementary finite-domain statements and explicit arithmetic example.** Boolean coordinate rings, linear parity inference and factor graphs are standard tools; no priority is claimed. The purpose is to move from asking what must be stored to constructing certificates that a stored or missing coordinate can be recovered.

## 1. What 'free' can mean

Fix a legal state set V, an observation O, and a requested bit b. The bit is recoverable without another measurement iff

\[
O(x)=O(y)\Longrightarrow b(x)=b(y)\qquad(x,y\in V).
\]

This does not create information. A constraint defining V or a previously available observation eliminates the alternative. The cost of obtaining that constraint, proving its applicability and computing the reconstruction must remain explicit.

For a linear state constraint Hx=s over F2 and a linear observation Ax, a linear target ell*x is recoverable iff

\[
\ell\in\operatorname{rowspan}\begin{pmatrix}H\\A\end{pmatrix}.
\]

A row-span expression is a reconstruction certificate, including its known affine constant. If no such expression exists, a vector v in the joint kernel with ell*v=1 gives an explicit indistinguishable pair x,x+v for every solution x. This statement assumes the affine constraint set is nonempty.

For several requested linear bits Lx, the minimum number of **additional arbitrary linear measurements** is

\[
\operatorname{rank}\begin{pmatrix}H\\A\\L\end{pmatrix}
-\operatorname{rank}\begin{pmatrix}H\\A\end{pmatrix}.
\]

Choose a basis of the missing row-space directions for the upper bound; each new scalar linear measurement increases rank by at most1 for the lower bound. This is not a theorem about arbitrary measurement costs or nonlinear legality constraints.

## 2. Codes plus a conflict complex

For pairwise coprime arithmetic atoms a_i, the parity vectors of prime exponents in sigma(a_i) form columns of a binary matrix H. A selected product has square sigma exactly when Hx=0. With composite atoms, coprimality can fail. Retain a conflict edge ij whenever gcd(a_i,a_j)>1; legality imposes x_i*x_j=0.

Thus a useful exact carrier for this finite selection problem is

\[
V=\{x\in\{0,1\}^m:Hx=0,\ x_ix_j=0\text{ for each conflict edge }ij\}.
\]

This is the intersection of a linear code with an independence complex. It is usually not a vector space: symmetric difference can leave the legal complex. The conflict graph has a different meaning from a graph whose cycle matroid is supposed to represent all arithmetic dependencies. A Fano obstruction to the latter does not forbid this richer representation.

This construction is exact for the stated **legal coprime selection language**. It does not assert that forbidden products fail to exist as integers, or that their actual sigma can still be computed by the illegal product of atomic sigmas.

## 3. An explicit one-bit gain from legality

Use the inherited seven-atom bank, in this order:

\[
(a_0,\ldots,a_6)=(4,2,5,3,7,13,11).
\]

In the square-class basis(2,3,7), the column masks are

\[
(4,2,3,0,1,5,2).
\]

There is exactly one conflict edge, between atoms4 and2, so x_0x_1=0. There are96 legal selections,16 formal zero-syndrome selections, and12 legal zero-syndrome selections. These counts and the illegal symmetric-difference example are inherited from [the corrected incidence packet](../incidence_attractors_20261028/THEORY.md).

Erase the four coordinates E={0,1,2,5}, corresponding to atoms{4,2,5,13}, while observing the others. The four columns have rank3. Without the legality constraint, every compatible observation has two completions, differing by1111 on E: exactly one bit is missing.

Take the observed pattern x_3=1,x_4=x_6=0, so atom3 is selected and atoms7,11 are not. The three parity equations on the erased coordinates become

\[
x_2+x_5=0,\quad x_1+x_2=0,\quad x_0+x_5=0.
\]

Hence all four erased bits are equal. Coprimality forbids x_0=x_1=1, so all four must be zero. The actual source is3, with sigma4 and aliquot edge3→1. The alternative formal selection is1560=4*2*5*3*13; its product of atomic sigmas is7056=84², but its actual sigma is5040. The alternative is removed by a known guard, not by an extra measurement.

The same recovery works with atom3 absent, yielding the unit source1. Including3 above makes clear that the example does not depend on allowing the unit as a positive square-sum edge.

### A hostile context where the bit is not free

Instead observe x_3=0,x_4=1,x_6=0: select7 and omit3,11. There are **two legal completions**:

- atoms{2,5,7}, giving70 with sigma144;
- atoms{4,7,13}, giving364 with sigma784.

Both obey all stated parity and coprimality constraints. Any method that claims to choose between them from these observations alone is incorrect. If an additional known predicate says F=U+S, it selects364; that predicate is extra information, not something already implied by square sigma.

More generally, the erased conflict coordinates obey x_0+x_1=x_4+x_6. If the observed right side is0, legality forces both to0 and determines the erased bits. If it is1, both opposite assignments are legal. This is an exact conditional-recovery rule, not a uniform four-erasure correction theorem.

Under the uniform distribution on the12 legal code selections, four observation fibers are singletons and four have size2. Therefore recovery without another measurement occurs with probability1/3 and the conditional entropy of the erased block is2/3 bit. For the unfiltered16-state linear code it is1 bit. The distributions and legality filter are part of these statements.

## 4. A checkable nonlinear certificate

Work in the Boolean coordinate ring

\[
R=\mathbb F_2[x_0,\ldots,x_{m-1}]/(x_i^2-x_i).
\]

Let I be generated by parity equations, conflict monomials and the observed-coordinate equations. In the recovering context above, the three parity equations imply x_0+x_1 belongs to I. Since x_0x_1 also belongs to I,

\[
x_0=x_0(x_0+x_1)+x_0x_1\pmod{x_0^2-x_0}
\]

is an explicit certificate that x_0 is forced to zero. The other erased coordinates follow from the linear relations. This combines a linear code with a nonlinear legality condition; a linear rank calculation alone cannot see the gain.

### Finite Boolean completeness theorem (standard algebra)

Evaluation gives an isomorphism

\[
R\cong\prod_{x\in\mathbb F_2^m}\mathbb F_2.
\]

For any finite list of polynomial constraints f_j=0, the quotient by their ideal is the ring of all F2-valued functions on their common legal zero set V. Indeed in the product ring the ideal has full component at every excluded point (some f_j is1) and zero component at every legal point.

Consequently, if O observes a subset of coordinates, a target coordinate x_k is recoverable on V iff there is a polynomial g in the observed variables with x_k+g in I. Sufficiency follows by evaluation. For necessity define g on each attained observation by its unique target value, extend arbitrarily to unattained observations, and use the Boolean truth-table polynomial representation.

This supplies both a language for higher-order arithmetic constraints and exact certificates of recovery. It promises neither small certificates nor an efficient algorithm: truth tables and elimination may grow exponentially. It does not decide unbounded arithmetic truth or solve open number-theory problems by encoding them.

## 5. What ordinary graphs do and do not exclude

A cycle-space representation restricts every graph vertex to an even-incidence law. Fano/Hamming dependencies do not fit that graphic-matroid restriction. A Tanner/factor graph instead carries named parity checks and equality rules, so its graph comes with different local relations. Such a graph can represent the same binary code exactly without being a cycle-space realization.

The useful extension is therefore **relations plus their composition laws**, not merely a larger unlabeled graph. Signed integer incidence adds lattice integrality and torsion; an LCA/tree arrangement adds separator messages; the code-plus-conflict example adds nonlinear legality. None of these resources can be silently discarded when claiming a bit was recovered for free.

## Reproduction

Run `python3 check_constraint_certificates.py` and with `-O`. The standard-library checker enumerates all128 selections, checks actual divisor sums on all96 legal selections, the exact erased-coordinate fibers, both named contexts, and the Boolean certificate at every Boolean point satisfying the context. It also exhaustively verifies the linear recoverability/rank criterion on small binary systems. All checks remain active under optimization.
