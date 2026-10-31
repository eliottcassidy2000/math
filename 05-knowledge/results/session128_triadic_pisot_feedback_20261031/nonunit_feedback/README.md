# Nonunit feedback: the missing finite-place geometry

**Status: PROVED scoped obstruction; FINITE-EXACT checks; CITED representation-space context.**
This separate addendum leaves the earlier Pisot/shear packet unchanged.

## Exact finding

The integral rational-resolvent cases recovered in the sibling triadic
packet include operators `2C^4` and `−103C^17`, where C represents
multiplication by the Tribonacci unit beta. The positive operators

    N_(2,4)=2M^4,          N_(103,17)=103M^17

are primitive substitution incidence matrices with positive Pisot expansion
factors `2 beta^4` and `103 beta^17`. The second is a **sign flip** of the
raw resolvent operator, not a positive-incidence realization of the raw
negative operator. The direct positive Pisot/substitution assertion applies
to `+103C^17`, not `−103C^17`.

Their graph-directed internal systems in the complex plane exist by
contraction. But **if their edge pieces are essentially disjoint, all
component areas must be zero**. Finite nonzero area and essential
disjointness cannot coexist. The exact obstruction is their nonunit norm:
`8` or `103^3=1092727`, rather than1.

## 1. Operator, substitution, and contraction checks

Use

    f(t)=t^3−t^2−t−1,
    beta>1 real, alpha its upper-half-plane conjugate,
    |alpha|^2=beta^(-1),
    M=[[1,1,1],[1,0,0],[0,1,0]],
    sigma=(12,13,1).

For positive integers q,n define the ordered substitution

    tau_(q,n)(j) = sigma^n(j) repeated q times by concatenation.

Its incidence matrix is `N=qM^n`. The original positive length vector
`ell=(1,beta−1,beta^(-1))` satisfies

    ell N=(q beta^n) ell.

Thus the expansion is `lambda=q beta^n`. On the complex internal space the
linear part is multiplication by `lambda*=q alpha^n`. Prefix digits come
from the actual repeated word, not merely from the incidence matrix. For
copy index `r=0,...,q−1`, an original prefix digit t in `sigma^n(j)` becomes
`t+r alpha^n ell_j(alpha)`. This specifies the ordered digit law without
materializing millions of letters.

Both displayed N are strictly positive. Since beta is a cubic unit and
beta>1, no positive power beta^n is rational: otherwise it would be a
rational algebraic-integer unit, hence±1. As the field has prime degree3,
`q beta^n` still has degree3. Its other two conjugates are `q alpha^n` and
its complex conjugate.

The earlier exact bound beta>9/5 yields

    (9/5)^4−2^2 = 4061/625 >0,
    (9/5)^17−103^2 = 8583157041463444/762939453125 >0.

Therefore `q^2 beta^(−n)<1` in each case. Each positive lambda is a Pisot
number; the negative raw eigenvalue `−103 beta^17` is not a Pisot number
under the standard positive-root definition, although its two remaining
conjugates contract. Contraction plus the finite prefix graph gives a
unique compact graph-directed tuple in C in either positive case. There
is no assertion here that it is a positive-area planar tile.

The exact matrices are

    2M^4 = [[14,12,8],[8,6,4],[4,4,2]],

    103M^17 = [[2009839,1686831,1092727],
               [1092727,917112,594104],
               [594104,498623,323008]].

Their characteristic polynomials are respectively

    x^3−22x^2−20x−8,
    x^3−3249959x^2+2875039x−1092727.

Their field norms and matrix determinants are q^3. They are nonunits, so
these changes cannot be conjugacies of the original unit operator: even
the determinant differs. This is distinct from the norm-preserving shear
family in the earlier packet.

## 2. General finite-measure obstruction

Let a finite graph have nonnegative incidence matrix N with a strictly
positive left eigenvector ell and eigenvalue lambda>0. Suppose measurable
sets R_i of finite measure satisfy a graph-directed equation with a common
invertible linear part B and translations indexed by every graph edge,
retaining multiplicity. Suppose the pieces **within each target equation**
are disjoint up to measure-zero intersections. Put `m_i=area(R_i)` and
`delta=|det_R B|`. Additivity and change of variables imply

    m=delta N m.

Multiplication by ell gives

    ell m = delta lambda ell m.

If m is nonnegative and nonzero, strict positivity of ell gives `ell m>0`,
so necessarily `delta lambda=1`.

For the present complex-plane systems,

    delta=|q alpha^n|^2=q^2 beta^(−n),
    lambda=q beta^n,
    delta lambda=q^3.

Hence q=2 and q=103 violate the necessary condition. All compact sets have
finite area, so any essentially disjoint compact solution has m=0. If at
least one component has positive area, positive-measure overlap of some
edge pieces is unavoidable. The argument does not decide which alternative
occurs for these planar attractors. It concerns edge pieces in one target
component, not an unsupported global claim about all component intersections.

## 3. Why finite places repair precisely this volume deficit

**CITED setup:** Minervino–Thuswaldner,
[*The geometry of non-unit Pisot substitutions*](https://arxiv.org/abs/1402.2002),
§3 and Lemma3.1, constructs the internal representation space using the
non-Perron Archimedean completions together with finite completions at prime
ideals dividing the inflation factor. Their Haar-scaling identity is the
reciprocal of the Perron expansion. Published version:
[Annales de l'Institut Fourier 64 (2014), DOI10.5802/aif.2884](https://numdam.org/articles/10.5802/aif.2884/).

Here beta is a unit also in the maximal order O_K, so

    (q beta^n)=(q) in O_K.

Thus the relevant finite places are precisely those over2 in the first
case and those over103 in the second; this note does not need a prime-ideal
factorization or an identification `O_K=Z[beta]`. In normalized local Haar
moduli their total scaling factor is

    product_(p | q) |q beta^n|_p = q^(−3).

Indeed beta is a unit at every finite place, and the norm of the ideal(q)
is q^3. Adjoining those finite factors changes the internal scaling from
`delta` to

    delta_full=delta/q^3=1/(q beta^n)=lambda^(-1).

Consequently the necessary measure-balance equation becomes consistent:
`delta_full lambda=1`. This is the exact missing geometric sidecar exposed
by the feedback transformation. Removing this obstruction does **not**
prove a tiling, coincidence, positive measure, or pure dynamical spectrum
for our two substitutions. No such conclusion is imported here.

## Reproducibility and boundaries

    python3 verify.py > normal.json
    python3 -O verify.py > optimized.json
    cmp normal.json optimized.json
    sha256sum -c SHA256SUMS

The dependency-free verifier checks matrix powers, positive entries,
characteristic polynomials, determinants, positive-length eigenvector
identities, exact contraction inequalities, and the two scaling balances.
It also independently verifies

    (I−gamma C)R=C

for `(gamma,R)=(0,C),(1/2,2C^4),(56/103,−103C^17)`. It does not re-prove the
sibling packet's global classification of rational integral resolvents.
There is no giant-word generation or numerical area estimate.
