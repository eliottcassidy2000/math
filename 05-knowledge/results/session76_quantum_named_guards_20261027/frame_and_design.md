# One integer defect controls a tournament frame and a symmetric design

The positive construction is classical:
[Renes, *Equiangular Tight Frames from Paley Tournaments*](https://arxiv.org/abs/math/0408287).
The point of the derivation below is the exact positivity test for an
arbitrary regular input and a short hostile certificate. No priority claim
is made for these elementary equivalences.

## Positivity is stronger than regularity

Let `K` be the real skew adjacency matrix of a regular tournament of odd
order `n≥3`: `K_ij=1` means `i→j`, and `K1=0`. Let `J` be the all-ones matrix.
Define

    P = ½(I+J/n−iK/√n),       D = K²+nI−J.

The following are equivalent:

1. `P` is positive semidefinite;
2. `P` is an orthogonal projection;
3. `D=0`;
4. the tournament is doubly regular, with `(n−3)/4` common out-neighbours
   for every distinct pair.

To prove the spectral part, `P` has eigenvalue one on the all-ones line.
On its orthogonal complement, pair the eigenvalues of the real skew matrix
as `±iσ_j`. The corresponding eigenvalues of `P` are
`(1±σ_j/√n)/2`. Since `tr(−K²)=n(n−1)`, the average of the `(n−1)/2` numbers
`σ_j²` is `n`. Positivity requires each at most `n`; hence they are all
exactly `n`. This is precisely `K²=J−nI`, and then `P` is a projection of
rank `d=(n+1)/2`.

For the combinatorial equivalence, on an arc `i→j`, partition the other
vertices into common out-neighbours `a`, common in-neighbours `b`, paths
`i→k→j` counted by `c`, and reverse paths counted by `e`. Regularity gives
`a=b`, `e=c+1`, and therefore

    D_ij = n−3−4a,    D_ii=0.

In particular the nontrivial positive case requires `n≡3 (mod 4)`.

## The classical frame and the exact parameter relation

When `D=0`, set

    G = (2n/(n+1))P = (nI+J−i√n K)/(n+1).

This is a positive semidefinite rank-`d` matrix with diagonal one and
nonzero eigenvalue `n/d`. Thus it is the Gram matrix of `n` unit vectors
forming an equiangular tight frame in complex dimension `d`. Its squared
off-diagonal overlaps are `1/(n+1)`.

Paley11, using quadratic-residue steps `{1,3,4,5,9}`, gives an `(11,6)`
frame with squared overlap `1/12`. No two distinct frame vectors are
orthogonal, so it is not a union of mutually unbiased orthonormal bases.
The relation to the score-distance number is only the exact parameter
identity

    n=2d−1 ⇒ (n²−1)/8 = d(d−1)/2 = dim so(d).

For `n=11`, this produces `d=6` and `15`. It does not identify fifteen
native moves with fifteen independent Lie generators.

## The same defect is design nonuniformity

Let `A=(J−I+K)/2` be the zero-one adjacency matrix. Regularity implies

    AAᵀ=AᵀA = ((n+1)I+(n−3)J−D)/4,
    P²−P = −D/(4n).

Writing `H0=((n+1)I+(n−3)J)/4`, these combine into the literal matrix identity

    P²−P = (AAᵀ−H0)/n.

Hence, with ordered distinct pairs in the final sum,

    ||P²−P||_F² = tr(D²)/(16n²)
                 = n^(−2) Σ_(i≠j) (c_ij−(n−3)/4)².

This is a precise Frobenius-defect bridge: failure of projection and
failure of pair-intersection uniformity are the same matrix error at a
specified scale.

At order eleven, `D=0` gives `AAᵀ=3I+2J`. The out-neighbour blocks form the
classical symmetric `2-(11,5,2)` biplane, of design order `5−2=3`. Its
complement has parameters `2-(11,6,3)` and Gram `3I+3J`. Thus `2,3,6,11`
occur as linked structural parameters. The construction from the residue
set is recorded in
[O'Reilly Regueiro, §1.2.1](https://www.matem.unam.mx/~eugenia/papers/paper1.pdf).

## An integer negative-vector certificate

`D` is symmetric with zero diagonal and zero row sums. If it is nonzero,
some entry `D_ij=d>0` exists. Put `x=e_i−e_j`, `y=Kx`, and
`z=√n x+iy`. Then

    ||y||²=2n+2d,      z*Pz=(n||x||²−||y||²)/2=−d<0.

This avoids numerical eigenvalue thresholds. The regular cyclic order-eleven
tournament with steps `{1,2,3,4,5}` has `tr(D²)=3520`. Its squared projection
defect is `20/11`; at `(i,j)=(0,4)` the above witness gives `z*Pz=−4`.
It has the same fifteen-move transitive-class distance as Paley11, but
fails this frame recipe.

Regularity is essential: the transitive tournament of order three is an
irregular input for which `P` is positive definite despite `D≠0`. Failure
of this particular recipe does not exclude other encodings or frames.
`check_frame_defect.py` reproduces the integer identities and the named
Paley/cyclic examples, not the full spectral/design literature.

## External quantum claims are not premises

The supplied MUB manuscript asserts `N(6)=3` under an explicit binary64
arithmetic/compiler contract and a reported completed computation. Its
[official formalization scope](https://github.com/openai/math/blob/main/lean/docs/266.md)
explicitly does not certify exclusion of four arbitrary bases; it describes
a weaker at-most-five bound and supporting cancellation statements. No
full native MUB or Lean replay is claimed here.

The [26 October QCD announcement](https://openai.com/index/building-scientific-understanding-with-academic-communities/)
reports existence results, but says the physics artifacts are not yet
public. An existence announcement is not an inspectable proof of a mass
gap or confinement. The finite-dimensional SU(3) identities in the companion
note import none of those continuum assertions.

The other supplied manuscripts concern local smoothing and ETH. Their
[selected formal scopes](https://github.com/openai/math/blob/main/lean/docs/079.md)
and [complexity scope](https://github.com/openai/math/blob/main/lean/docs/378.md)
must be distinguished from sharper manuscript statements. Littlewood–Paley
terminology is not a Paley tournament map. A naive quadratic-variable SAT
encoding gives no automatic speedup over exponential-in-vertices tournament
subset dynamic programming.
