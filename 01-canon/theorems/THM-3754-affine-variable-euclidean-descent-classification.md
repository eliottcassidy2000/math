---
id: THM-3754
title: "Affine-variable smoothness and Jacobian-mate Euclidean descent"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED INHERITED
  SYNTHESIS/REFINEMENT.  Over an algebraically closed characteristic-zero
  field, the nonsingular polynomials
  Q=f(X)+Tg(X) are classified exactly by the constant boundaries or, for
  nonconstant g, by the condition that every root of g is multiple and no
  root of g is a root of f'.  A polynomial Jacobian mate exists exactly when
  g is a nonzero constant, or when g=0 and f is affine nonconstant.  For
  every nonconstant g a top-T logarithmic-derivative equation and subtraction
  of powers of Q descend any hypothetical mate to the impossible equation
  g p'=c.  Euclidean division makes f mod g the retained source-shear
  sidecar; X+X^2+X^3T is the first nonradial pure-power example.  This
  concisely repackages THM-2063, THM-3326, and THM-3406; it is not a newly
  closed JC stratum.
source: root + jc_sparse_direct_search / 2026-08-23
audit: >
  PASS.  An independent hostile audit checked every rootwise smoothness case,
  both constant boundaries, the radical/gcd reformulation, the top-T sign,
  constants of k(X), strict power-subtraction descent, terminal equation,
  Euclidean source gauge, and the first nonlinear residue example.  Normal,
  optimized, and frozen output agree; script/output/semantic hashes and
  CHECKS=6622 match.  The audit also checked the inheritance boundary against
  THM-2063, THM-3326, and THM-3406; the concise field-uniform proof and exact
  controls refine their presentation rather than claim a new stratum.
depends_on:
  - THM-2063-one-fiber-linear-planar-keller-pairs
  - THM-3326-linear-in-z-unit-response-trichotomy-and-jet-torsion
  - THM-3406-affine-modification-power-jets-and-principal-part-transgression
related:
  - THM-3716-monomial-broughton-hamiltonian-obstruction-family
  - THM-3741-radial-two-charge-keller-component-classification
script: 04-computation/jc2_affine_variable_euclidean_descent_thm3754.py
output: 05-knowledge/results/jc2_affine_variable_euclidean_descent_thm3754.out
script_sha256: d17141a3db956f541758a53690be8dfb5a7ceb581f238643a3107cffdc315339
output_sha256: fd4de5bd112fc2132f6b55c39c0b9ef863e34c2286de186c3572a9b6df40d1d5
semantic_sha256: d4d1e92e034e1eeda583d1effc0b833c1c24ebd71dfc0d704c8cfe836e9223d8
hash_basis: raw LF bytes
---

# THM-3754 -- affine-variable components descend by their own powers

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED INHERITED
SYNTHESIS/REFINEMENT.**  THM-2063 already closes full Keller pairs with an
output affine in a source direction; THM-3326 gives the same rootwise
smoothness and unit-response boundary; and THM-3406 records the corresponding
affine-modification filtration and residue sidecars.  The value here is a
short field-uniform synthesis: it puts the smoothness classification, direct
power-subtraction descent, Euclidean residue gauge, and exhaustive exact
controls on one surface.  It does **not** close a new JC stratum.

After THM-3741 closes every radial two-charge component, this inherited class
is still a useful boundary: keeping one source variable affine while letting
its coefficient and constant term vary independently remains completely
rigid.  Powers of the component itself form the exact kernel needed to
eliminate the highest affine-variable degree of a prospective mate.

Let `k` be an algebraically closed field of characteristic zero, take
`f,g in k[X]`, and put

```text
Q(X,T)=f(X)+Tg(X).                                      (1)
```

## 1. Complete smoothness and mate classifications

The polynomial `Q` has no critical point in `k^2` exactly in the following
cases:

```text
(C*) g is a nonzero constant; f is arbitrary;
(C0) g=0 and f=aX+b with a!=0;
(N)  deg g>=1, every root alpha of g is multiple,
     and f'(alpha)!=0 at every such root.               (2)
```

Writing `rad(g)=g/gcd(g,g')` up to a nonzero scalar, the nonconstant
criterion can equivalently be recorded without listing roots as

```text
rad(g) divides g',                    gcd(g,f')=1.       (3)
```

Among the nonsingular list, a polynomial `P` with `J(P,Q) in k*` exists if
and only if `Q` lies in `(C*)` or `(C0)`.  Consequently every member of `(N)`
is a nonsingular noncoordinate polynomial, and the complete nonconstant-`g`
family is empty as a planar Jacobian-counterexample search space.

## 2. Smoothness is a rootwise vertical-fibre test

This is the concise algebraically closed form of THM-3326's root-support
trichotomy.

The gradient is

```text
Q_X=f'(X)+Tg'(X),                  Q_T=g(X).            (4)
```

If `g` is a nonzero constant, the second entry never vanishes.  If `g=0`,
the first entry has no zero over the algebraically closed field exactly when
it is a nonzero constant, proving `(C*)` and `(C0)`.

Suppose now that `g` is nonconstant and let `alpha` be one of its roots.  If
`g'(alpha)!=0`, the point

```text
(X,T)=(alpha,-f'(alpha)/g'(alpha))                     (5)
```

is critical.  If `g'(alpha)=0`, then the whole vertical line `X=alpha` is
safe exactly when `f'(alpha)!=0`; otherwise every point of that line is
critical.  These alternatives prove `(2)` root by root.  Over characteristic
zero, `g'(alpha)=0` at every root exactly when every root is multiple, which
is the first divisibility in `(3)`; the second condition in `(3)` is the
remaining avoidance statement.

## 3. The highest T-row is an exact logarithmic derivative

This is THM-2063's coefficient descent written directly for the chosen
component and over an arbitrary algebraically closed characteristic-zero
field; the proof is retained here to make the synthesis self-contained.

The no-mate result is stronger than the smoothness theorem: it holds for
every nonconstant `g`, whether or not `(1)` is nonsingular.  Suppose toward a
contradiction that

```text
J(P,Q)=c in k*,                                        (6)
```

and expand

```text
P=sum_(j=0)^N p_j(X)T^j,                p_N!=0.         (7)
```

If `N>=1`, the coefficient of `T^N` in `(6)` is

```text
g p_N'-N g'p_N=0.                                     (8)
```

In the rational function field `k(X)`, this is precisely

```text
(p_N/g^N)'=0.                                         (9)
```

The constants of this derivation are `k`, so

```text
p_N=lambda g^N,                    lambda in k*.       (10)
```

But `Q^N` has top `T`-coefficient `g^N` and
`J(Q^N,Q)=0`.  Replacing

```text
P by P-lambda Q^N                                     (11)
```

therefore preserves `(6)` and strictly decreases `deg_T P`.  Repeating
`(8)--(11)` finitely often leaves `P=p(X)`.  Equation `(6)` has then become

```text
g(X)p'(X)=c.                                          (12)
```

This is impossible when `g` is nonconstant.  That proves the universal
obstruction.

On the two surviving boundaries the mates are explicit:

```text
g=d in k*:          P=X/d;
g=0, f=aX+b:        P=-T/a.                            (13)
```

Both have `J(P,Q)=1`, completing the mate classification.  Notice that the
descent discovers its own correction basis: it does not guess an ansatz for
the lower rows of `P`; it quotients them successively by powers of the very
component whose Hamiltonian derivation is being inverted.

## 4. Euclidean residue is the missing sidecar

This elementary gauge is the source-coordinate shadow of THM-3406's fuller
affine-modification and principal-parts sidecar.

Divide

```text
f=qg+r,                              deg r<deg g.       (14)
```

The triangular source automorphism `T_tilde=T+q(X)` rewrites `(1)` as

```text
Q=r(X)+T_tilde g(X).                                  (15)
```

Thus the residue class `f mod g`, not the full polynomial `f`, is the datum
retained after source-shear gauge reduction.  This quotient does not lose
the smoothness test.  At a multiple root `alpha` of `g`, both
`g(alpha)` and `g'(alpha)` vanish, while differentiating `(14)` gives

```text
f'(alpha)=r'(alpha).                                  (16)
```

For the pure-power coefficient `g=X^m`, smoothness is simply `f'(0)!=0`,
and `(15)` reduces `f` to degree below `m`.  When `m=2`, that residue is only
affine and recovers the monomial Broughton boundary of THM-3716.  The first
genuinely nonlinear retained residue occurs at `m=3`:

```text
Q=X+X^2+X^3T=X(1+X+X^2T).                             (17)
```

It is a quartic nonsingular noncoordinate with a reducible zero fibre.  Its
extra `X^2` charge is not removable by the triangular source gauge, yet the
Euclidean descent still proves that it has no mate.  This example separates
two lessons that a search must retain simultaneously: leaving the radial
profile creates real moduli, but remaining affine in one variable keeps an
exact descending obstruction.

## 5. Exact controls and the next boundary

Reproduce the exact audit surface with

```bash
python3 -B 04-computation/jc2_affine_variable_euclidean_descent_thm3754.py
python3 -B -O 04-computation/jc2_affine_variable_euclidean_descent_thm3754.py
```

The assertion-free companion verifies the universal top-row equation and
power subtraction through seven symbolic depths; exhausts all 6,561 pairs
of coefficient vectors of degree at most three over `{-1,0,1}` by direct
two-variable Groebner ideals; checks Euclidean residues, smooth and singular
boundaries, the quartic `(17)`, both positive mate boundaries, and 24 bounded
hostile mate systems.  These are exact controls; the finite descent
`(8)--(12)` proves the arbitrary-degree mate theorem.  **QED.**

A candidate component must now be nonlinear in both variables after every
triangular source reduction.  Equivalently, no coordinate direction may
leave a finite filtration whose leading coefficient is generated by powers
of `Q`.  Quadratic dependence in the second variable is the first unresolved
filtration boundary; multi-charge substitutions such as
`Q=X+F(X^mT)` provide another controlled way to cross it.
