# Independent audit: coordinate symmetry and cyclic-parent ownership

**Session 114, 30 October 2026. Verdict: the scoped statements audited below
are valid.** This is an independent algebraic review and a separately written
finite checker, not formal proof-assistant certification or a literature
priority assessment. No major open conjecture is declared solved.

## 1. One-coordinate Keller theorem

Reviewed `jacobian/README.md`, Theorem 1, Corollary 2, Theorem 3 and
Corollary 4. The field has characteristic zero, the mixed weights are positive
coprime integers `(p,-q)`, and the determinant is a nonzero constant.

The decisive argument is valid without any symmetry or degree restriction
on the other coordinate. If `Q=y B(x^q y^p)`, the bracket with the
integer-weight-r part of P has weight `r-p`. Only the weight-p component
`x A(x^q y^p)` can supply the nonzero constant. Its bracket is

    AB + q t A'B + p t AB'.

The top coefficient is `(1+q deg A+p deg B) lead(A) lead(B)`, which is
nonzero. Thus A and B are nonzero constants. The original equation then
forces `P=a x+f(y), Q=b y`. Projection is used on a polynomial identity;
there is no unsupported claim that coefficient averaging preserves Keller
maps. The symmetric assertion follows by the same decomposition.

The finite sufficient cutoffs are correctly strict:

    Q: m > max(q(d-1), pd+q),
    P: m > max(qe+p, p(e-1)).

These prevent nonzero multiples of m in the supported weight-difference
intervals. They are not asserted optimal in every parameter. At balanced
weights the shear `Q=y+x^d` shows why weakening the strict inequality can
invalidate the specified triangular direction.

At the sharp full-map boundary `m=wD+1`, the classification requires
`D>=2`. For `w>1`, Q satisfies the one-coordinate no-alias criterion. For
`w=1`, the two off-zero determinant channels force
`alpha B'=beta A'=0`; simultaneous nonzero aliases then leave an uncancelled
`-alpha beta D^2 t^(D-1)` term. This proves triangularity. The excluded
degree-one case really has scalar action and allows every invertible linear
map. At balanced degree at most m, the extra degree adds no further alias
monomial, so the same proof applies for `m>=3`.

The audit checker independently enumerates all supported monomials for its
declared finite ranges, differentiates directly in x,y, extracts actual
weight-zero and off-zero channels, and verifies mixed-weight leading
coefficients. It does not import the producer's implementation.

## 2. Cyclic-parent count and collective phase

Reviewed `charts/PARENT_MULTIPLICITY.md`, sections 1--4. In
`W=C_m wr C_r`, an element `(t,z)` is a regular `rm`-cycle exactly when
`gcd(t,r)=gcd(sum z_i,m)=1`. There are
`phi(r)phi(m)m^(r-1)` such elements, and each regular cyclic subgroup has
`phi(rm)` generators. Therefore the parent count is

    K = phi(r)phi(m)m^(r-1)/phi(rm)
      = m^(r-1) phi(g)/g,  g=gcd(r,m).

The generator-independent slope is `rho=(sum z_i)/t mod g`, and the
number of parents at each unit slope is `m^(r-1)/g`. Conjugating the
inherited two-parent generators by `(i,j)->(i,nu j)` transports the kernel
`sum z_i=t mod g` to `sum z_i=rho t mod g`. Thus the phase-group assertions
correctly inherit the earlier unused-fibre and pair-closure proofs.

**Arithmetic correction before publication:** when `m=3,r=3^k`, the exact
count is `2*3^(r-2)`, not `3^(r-2)`. The final producer document contains
the corrected factor two. The smallest case has six parents, not three.

For a fixed legal path, all parent charts use the same invertible
complement mask. Intersection therefore reduces to common-source
intersection, and the two inherited generators already force the full
W-invariant bank. Consistent pins leave exactly `2^(d-q(P))` sources.
The resulting ownership is exactly K **within the specified roster**;
additional parents outside W are not excluded.

This is a count of distinct named cyclic subgroups, not a claim that all
chart supports are distinct. A singleton intersection of the entire roster
does not imply pairwise singleton intersections. The source-cover lower
bound follows by counting legal-source sets and does not bound a fixed
output family when different paths use different masks. These distinctions
are explicit in the audited statement.

The independent checker enumerates literal wreath-product permutations,
tests their cycles by traversal, forms cyclic subgroups from powers,
verifies generator-independent slopes, and generates the stated phase
groups by permutation composition. It directly tests every product-bank
tournament against every counted parent in four finite cases. Expected
parent counts are 6,25,81,49 for `(r,m)=(3,3),(3,5),(5,3),(3,7)`; natural
common-source counts are 2,4,4,8, and all selector intersections are one.

## 3. Reproduction and limits

Run from this directory:

    python check_symmetry_independent.py
    python -O check_symmetry_independent.py

Both use explicit exception gates. `symmetry_independent.json` contains
the deterministic receipt. The checker requires Python 3 and SymPy. It
checks 3,130 exact Jacobian gates, plus the explicitly reported finite
parent universes and all invariance tests in those universes. The finite
tests are controls for the all-parameter proofs, not exhaustive tests of
all Keller maps or tournaments.

The work does not manufacture coordinate symmetry for an arbitrary Keller
map, establish a universal chart-cover lower bound, or resolve planar JC
or unrestricted tournament enumeration.
