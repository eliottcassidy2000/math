# THM-441 — Convolution, correlation, and the adjoint with the carrier retained

**Status:** PROVED elementary identities below. **SUPERSEDED:** the historical
self-converse=self-adjoint, automatic-extremality, cross-modulus-adjoint, and
unqualified energy/unit-distance identifications. Targeted correction2026-10-11.
**Original source:** opus-2026-06-07-S706, from the user's convolution/correlation
seed. Historical artifacts and companion reflection are retained for provenance.

## D1–D3. Exact identities on one finite abelian group

Fix a finite abelian group G, inner product <a,b>=sum_x conjugate(a(x))b(x),
and the unnormalized Fourier convention
hat f(chi)=sum_x f(x) conjugate(chi(x)). Define

\[
 (h*g)(x)=\sum_y h(y)g(x-y),\qquad
 (h\star g)(x)=\sum_y\overline{h(y)}g(x+y).
\]

Then C_h^*=R_h, where C_h g=h*g and R_h g=h star g.
Indeed, change variables z=x-y in the finite inner-product sum. Also,
with (sigma h)(x)=h(-x),

\[
 h\star g=(\sigma\overline h)*g,\qquad
 \widehat{h*g}=\widehat h\,\widehat g,\qquad
 \widehat{h\star g}=\overline{\widehat h}\,\widehat g.
\]

These prove the adjoint, reflection, and conjugate-symbol laws directly.
C_h is self-adjoint precisely when h(x)=conjugate(h(-x)).

## M1. Sum and difference multiplicities retain the group

For an indicator1_S on this same G, 1_S*1_S counts ordered sums and
1_S star1_S counts ordered differences. Reflection relates them exactly.
The repository's shell modulo2n-1 and clock modulo n use **different
finite groups**. Their convolution identities are two instances of this
construction; they are not literally adjoints on one space without a
specified map between the spaces. Coprimality does not supply that map.

## M2. Converse is transpose; self-converse is conjugacy, not equality

For a circulant tournament, its real adjacency A is a convolution
operator and its converse has matrix A^T=A^*. A self-converse tournament
has a permutation matrix P with PAP^{-1}=A^T. A multiplier carrying its
connection set H to -H is one sufficient explicit such isomorphism.
This says **unitary equivalence to its adjoint**, not A=A^*.

The directed triangle is the minimal hostile:

\[
 A=\begin{pmatrix}0&1&0\\0&0&1\\1&0&0\end{pmatrix}.
\]

Swapping vertices1 and2 conjugates A to A^T, but A!=A^T. In fact the
adjacency of any tournament on at least two vertices is nonsymmetric,
because exactly one of each pair of off-diagonal entries is1. Its
skew-adjacency A-A^T is skew-Hermitian, a separate exact assertion.
Self-converse symmetry therefore does not force a real spectrum or
extremality of an unrelated objective.

## M3. Energy and distance counts are different functionals of correlation

For S subset G, finite Parseval with the declared normalization gives

\[
 E(S)=\#\{(a,b,c,d)\in S^4:a+b=c+d\}
     =\|1_S\star1_S\|_2^2
     =\frac1{|G|}\sum_\chi|\widehat{1_S}(\chi)|^4.
\]

The factor1/|G| is required under D1–D3's convention. Autocorrelation has
nonnegative Fourier transform, so its convolution operator is positive
semidefinite. This holds for every S; it does not by itself identify a
maximizer of energy or a tight LRC configuration.

For a finite Euclidean point set P, define the finite-support displacement
count r_P(z)=#{(a,b) in P²:b-a=z}. The **unordered** unit-distance count is

\[
                         U(P)=\tfrac12\sum_{\|z\|=1}r_P(z).
\]

Additive energy instead sums r_P(z)² over every displacement. Both use the
same correlation, but they are different functionals, with no equivalence
of their optimization problems proved here. A finite-group model also
needs an explicit embedding and a retained distance predicate.

## Surviving research use and correction lineage

The reusable operation is to retain the object and inner product before
passing to its adjoint, and to retain a chosen marking before identifying
isomorphic observations. The repaired
[THM-406, covering-depth master object](THM-406-covering-depth-master-object-factorial-moments-and-spectral-identity.md)
likewise distinguishes the cyclic subspace from the ambient space and
full finite moments from a proper truncation. Neither theorem excludes
all moment or measure arguments for existence.

The original document's slogans identifying self-converse with Hermitian,
self-adjoint with tight, and all repository objects with convolutions are
retracted as unproved or false. The directed-triangle control is in
[the Vitali recovery](../../05-knowledge/results/vitali_atom_recovery_20261011.md)
and its exact companion script. No new LRC, colouring, or Collatz bound is
claimed. The underlying classical D1–D3 identities survive unchanged
apart from making Fourier normalization explicit.

Historical artifacts: `04-computation/convolution_correlation_adjoint_s706.py`
and its output; reflection
`07-reflections/convolution-correlation-adjoint-unifies-clock-shell-converse-s706.md`;
HYP-2311. Their original numerical summaries are not a fresh audit here.
