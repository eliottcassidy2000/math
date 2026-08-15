---
id: THM-3418
title: "One-monomial nonlinear-fiber Keller classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY DERIVATION/TYPE/REPLAY-AUDITED.
  Let K be a characteristic-zero field, d>=2,
  P=f(x)+g(x)z^d, and Jac(P,Q)=kappa in K*.  Then necessarily
  f=ax+b and g=c are constant/affine in x, with a!=0, and every mate is
  Q=(kappa/a)z+H(P).  Conversely all such pairs have Jacobian kappa,
  an explicit polynomial inverse, and a tame factorization.  Before the
  recurrence is used, gradient unimodularity has the exact boundary
  f'=a in K* and g|(g')^2 (for nonzero g), equivalently every geometric
  root of g is repeated.  For nonconstant g the residue-one coefficient
  tower has degree k(deg(g)-1) and a never-zero leading multiplier, so a
  polynomial mate is impossible.  This uniformly closes the sparse
  top-plus-constant fiber stratum for every d>=2; it does not close any
  fiber polynomial with an intermediate z-coefficient and does not prove
  JC(2).
source: root-2608-jc-nonlinear-monomial-fiber-2026-08-15
audit: exact coefficient derivation; algebraic-closure gradient gate; Fraction polynomial recurrence/telescoping/sector/inverse replay; normal and optimized outputs byte-identical; independent gradient/recurrence/nontermination/kernel/inverse/tameness/sector-colimit/module-type/hash/routing/scope audit clean
depends_on: []
related:
  - THM-2063-one-fiber-linear-planar-keller-pairs
  - THM-2071-quadratic-fiber-square-parity-gate
  - THM-2118-all-degree-cubic-faber-boundary-flux-coprimality
  - THM-3412-hamiltonian-principal-part-differential-and-prufer-torsion-arms
script: 04-computation/jc_one_monomial_nonlinear_fiber_thm3418.py
output: 05-knowledge/results/jc_one_monomial_nonlinear_fiber_thm3418.out
script_sha256: 8ae459feb969e00cbb92acd718cdeb49b0d20d86323cdb98e6c31a8fd6a4fe2f
output_sha256: 4bce9383f966c05468d4c61de5e7bdd3ce77f4a0ce9a8e997e0fed4cfc5dfe05
semantic_sha256: 90d08019ad87638fecc375d76258e282bc6a4de17f53e01ca87802c8f05c5c92
hash_basis: LF-normalized bytes
---

# THM-3418 -- one-monomial nonlinear-fiber Keller classification

**PROVED + VERIFIED-EXACT.**

## 1. Inheritance and connection contract

[THM-2063](THM-2063-one-fiber-linear-planar-keller-pairs.md) uses
coefficient descent to classify every planar Keller pair with a component
linear along one source fiber.  THM-2071 and THM-2118 close the much larger
quadratic- and cubic-fiber strata.  The present sparse family keeps only the
constant and top fiber coefficients but allows the top degree to be
arbitrarily large:

```text
P=f(x)+g(x)z^d,                         d>=2.             (1)
```

The closest corrected near miss is the linear-`z` response tower of
[THM-3412](THM-3412-hamiltonian-principal-part-differential-and-prufer-torsion-arms.md).
There `z=0` is not a critical slice because `P_z=g`; here the factor
`z^(d-1)` changes the gradient gate completely.  The least-used sidecar is
the exact `Z/dZ` fiber grading of the Hamiltonian operator.  Its character
one chain, not a total-degree estimate, supplies the obstruction.

| item | exact content |
|---|---|
| source | a planar Keller equation with `P` as in `(1)` |
| target | a complete normal form for `Q`, or a typed obstruction |
| map | the Hamiltonian differential `D_P=Jac(P,-)` |
| preserved | fiber exponent modulo `d`, polynomial degree in `x`, and leading coefficient |
| destroyed by total fiber degree | the forced residue-one chain and its terminal defect |
| required sidecar | the coefficient index modulo `d` |
| cheapest hostile | `P=x+x^2z^2`: its gradient is unimodular, but its forced chain never terminates |

This is a genuine additional planar-Jacobian stratum, but it is a narrow
one.  A polynomial with any intermediate `z`-coefficient is outside the
statement.

## 2. Classification theorem

Let `K` be a field of characteristic zero, let `d>=2`, and take

```text
f,g in K[x],             P=f(x)+g(x)z^d,              Q in K[x,z].   (2)
```

The polynomial `g` is allowed to be zero.  For a fixed `kappa in K*`, the
following are equivalent:

1. `Jac(P,Q)=kappa`;
2. there are `a in K*`, `b,c in K`, and `H in K[T]` such that

```text
f=ax+b,                  g=c,
Q=(kappa/a)z+H(P).                                            (3)
```

Every pair in `(3)` is a tame polynomial automorphism.  In particular, a
hypothetical counterexample to `JC(2)` cannot have, in any affine
source/output-pencil chart, a component of the two-term form `(2)` with
fiber degree at least two.

## 3. The exact gradient boundary

If `Jac(P,Q)=kappa`, then the gradient row of `P` is unimodular, because

```text
P_x Q_z-P_z Q_x=kappa.                                      (4)
```

For `(2)`,

```text
P_x=f'(x)+g'(x)z^d,              P_z=d g(x)z^(d-1).         (5)
```

Extend scalars to an algebraic closure.  This does not change whether the
ideal `(P_x,P_z)` is the unit ideal.  On the slice `z=0`, the second entry
vanishes, so `f'` has no zero.  Hence

```text
f'=a in K*.                                                (6)
```

Now let `alpha` be a root of a nonzero `g`.  The second entry of `(5)`
vanishes for every `z`, while

```text
P_x(alpha,z)=a+g'(alpha)z^d.                               (7)
```

If `g'(alpha)!=0`, algebraic closure supplies a `d`th root of
`-a/g'(alpha)`, producing a common zero.  If `g'(alpha)=0`, `(7)` is the
nonzero constant `a`.  Conversely, every zero of `P_z` has either `z=0` or
`g(x)=0`, and the preceding two checks then exclude a common zero.  Thus,
for nonzero `g`,

```text
(P_x,P_z)=K[x,z]
iff f'=a in K* and every geometric root of g is repeated
iff f'=a in K* and g divides (g')^2.                       (8)
```

The last equivalence is the factor-valuation test in characteristic zero.
For `g=0`, only `f'=a in K*` remains.  Notice the sharp contrast with the
linear-`z` gate: the new slice `z=0` forces the whole base polynomial `f` to
be affine.

## 4. The unavoidable residue-one tower

Assume `(4)`.  By `(6)`, write `f=ax+b`.  Expand

```text
Q=sum_(n=0)^N q_n(x)z^n.                                   (9)
```

With `D=D_P=Jac(P,-)`, direct differentiation gives

```text
D=(a+g'z^d) partial_z-dg z^(d-1) partial_x.                (10)
```

The coefficient of `z^m` in `D(Q)=kappa` is

```text
a(m+1)q_(m+1)
 + 1_(m>=d-1) ((m-d+1)g' q_(m-d+1)-dg q'_(m-d+1))
 = kappa 1_(m=0).                                         (11)
```

In particular `q_1=kappa/a`.  The equations with exponents divisible by
`d` involve only input exponents congruent to one modulo `d`.  Put

```text
A_k=q_(1+kd),                  n_k=1+kd.                  (12)
```

Then `(11)` forces

```text
A_0=kappa/a,
a n_(k+1) A_(k+1)=d g A'_k-n_k g' A_k.                   (13)
```

Suppose `g` is nonconstant, with

```text
r=deg(g)>=1,                     gamma=lc(g).             (14)
```

Induction in `(13)` gives

```text
deg(A_k)=k(r-1).                                           (15)
```

Indeed, if `ell_k=lc(A_k)`, the candidate top coefficient in the numerator
of `(13)` is

```text
gamma ell_k (d k(r-1)-(1+kd)r)
 =-gamma ell_k(dk+r).                                     (16)
```

It is never zero in characteristic zero.  More precisely,

```text
ell_(k+1)
 =-[gamma(dk+r)/(a(1+(k+1)d))] ell_k,       ell_0=kappa/a. (17)
```

Thus every `A_k` is nonzero.  This contradicts the finite `z`-degree in
`(9)`, proving that a Keller mate forces `g` to be constant.

There is a useful exact failure identity.  The truncation of the mandatory
character-one chain,

```text
Q_K^(1)=sum_(k=0)^K A_k z^(1+kd),                         (18)
```

satisfies

```text
D(Q_K^(1))
 =kappa-a(1+(K+1)d)A_(K+1)z^((K+1)d).                    (19)
```

Every internal coefficient cancels; the obstruction moves outward but never
dies.  Formula `(19)` explains the failure mechanism more sharply than the
bare degree contradiction.

For the smallest gradient-unimodular hostile `P=x+x^2z^2`, the chain begins

```text
A_0=1,              A_1=-(2/3)x,              A_2=(8/15)x^2, ...    (20)
```

The nonsplit rational control `g=(x^2+1)^2` has the same never-terminating
mechanism; splitting roots is not used in `(13)--(17)`.

## 5. Constant `g`: every mate and the inverse

It remains to take `g=c in K`, including `c=0`.  Then

```text
P=ax+b+cz^d,                  K[x,z]=K[P,z],             (21)
```

because `x=(P-b-cz^d)/a`.  In the coordinates `(P,z)`,

```text
D(P)=0,                       D(z)=a,
D=a partial_z.                                             (22)
```

Characteristic zero makes the kernel of `partial_z` on `K[P,z]` exactly
`K[P]`.  Therefore `D(Q)=kappa` is equivalent to the unique form `(3)`.
Conversely `(22)` verifies its Jacobian immediately.

For target coordinates `(u,v)=(P,Q)`, the inverse is

```text
z=(a/kappa)(v-H(u)),
x=(u-b-cz^d)/a.                                            (23)
```

It is polynomial.  The tame factorization is also explicit: first apply the
affine scale/translation in `x` and the elementary shear
`x |-> ax+b+cz^d`, then scale `z` and apply the target shear
`z |-> (kappa/a)z+H(P)`.  This proves the equivalence and tameness.  QED.

## 6. Exact Hamiltonian-sector presentation

The proof above uses only one sector of a larger exact response object.  This
section records the full elementary presentation without claiming a closed
form for the cokernel.

Continue with `f=ax+b`, `a!=0`, but allow arbitrary `g`.  Give `K[x,z]` its
`Z/dZ` grading by the exponent of `z`.  Then `P` has degree zero and `D` in
`(10)` has degree `-1`.  Hence

```text
C_P=K[x,z]/D(K[x,z])=direct_sum_(j mod d) C_j,             (24)
```

as a `K[P]`-module, where the target sector `C_j` is the quotient of the
`j`th fiber-exponent sector by the image of sector `j+1`.

For `n>=1`, define the first-order operator

```text
L_n(q)=(d/(an))gq'-(1/a)g'q.                              (25)
```

For `1<=s<d`, the relation obtained by applying `D` to
`qz^(s+kd)` is exactly

```text
[q]_k=[L_(s+kd)(q)]_(k+1).                                (26)
```

Consequently

```text
C_(s-1)
 ~= colim (K[x] --L_s--> K[x] --L_(s+d)--> K[x] --> ...). (27)
```

Here and below the displayed colimit is first taken in `K`-vector spaces.
Its `K[P]`-module structure is induced rather than stagewise: writing
`u=ax+b`, multiplication by `P` is

```text
P [q]_k=[u q]_k+[g q]_(k+1).                            (27a)
```

This is compatible with every transition.  Indeed, for `n=s+kd`,

```text
L_n(uq)+gq-uL_n(q)=((n+d)/n)gq,
L_(n+d)(((n+d)/n)gq)=gL_n(q).                           (27b)
```

Thus `(27a)` gives exactly the action inherited from multiplication by
`P=u+gz^d` on the cokernel, and `(27)` is an isomorphism of `K[P]`-modules
with this action.

The wrap sector has one additional boundary relation.  Since

```text
D(q(x))=-dgq' z^(d-1)                                    (28)
```

and differentiation is onto `K[x]`, its first stage is `K[x]/(g)`.  The
exact presentation is

```text
C_(d-1)
 ~= colim (K[x]/(g) --L_d--> K[x]/(g^2)
             --L_(2d)--> K[x]/(g^3) --> ...).             (29)
```

The maps are well-defined because

```text
L_((k+1)d)(g^(k+1)h)=g^(k+2)h'/(a(k+1)).                 (30)
```

Since differentiation is onto `K[x]`, `(30)` also shows inductively that
the relation ideal at the next stage is exactly `(g^(k+2))`, not merely a
subideal.  The action `(27a)` applies to `(29)` as well.  It is independent
of representatives because changing `q` at stage `k` by `g^(k+1)h` changes
its two terms by multiples of `g^(k+1)` and `g^(k+2)`, respectively; its
compatibility with the transition maps is again `(27b)`.

Equations `(24)--(30)` are a Hamiltonian-response presentation, not a new
Keller argument: the classification already follows from the single chain
`s=1`.  They are also structurally different from THM-3412's linear-`z`
affine-modification tower.

## 7. Boundary and scope audit

- `d=1` is excluded.  In that case `P_z=g`, so `z=0` does not force `f'` to
  be constant; THM-2063 supplies the correct larger classification.
- `d=2` and `d=3` are already contained in the broader THM-2071 and
  THM-2118 closures.  The new point is the one-line mechanism uniform in
  every `d`.
- The gradient gate alone is not the conclusion.  `P=x+x^2z^2` passes it
  but has no mate; the first failed implication is termination of `(13)`.
- Total fiber degree does not retain the proof.  The decisive coordinate is
  the input exponent modulo `d`.
- No statement is made for
  `f(x)+g_1(x)z+...+g_d(x)z^d` with any nonzero intermediate coefficient,
  for characteristic `p`, or for `JC(2)` as a whole.

## 8. Exact computational referee

The standard-library companion uses exact `Fraction` polynomial arithmetic
and no `assert` truth gate.  It checks the gradient boundary on the affine,
simple-root, repeated-root, nonsplit repeated, constant, and zero controls.
It then verifies `(13)--(19)` independently by direct bivariate Jacobian
calculation for `35` nonconstant `(d,g)` cells, checking `245` nonzero
degree/leading/telescoping coefficients.  Another `525` exact checks replay
the sector relations and wrap-quotient identity `(25)--(30)`.  Finally it
checks `21` constant-`g` normal forms and both coordinates of `(23)`.
Ordinary and optimized runs are byte-identical.

Reproduce with

```text
python3 04-computation/jc_one_monomial_nonlinear_fiber_thm3418.py
python3 -O 04-computation/jc_one_monomial_nonlinear_fiber_thm3418.py
```

Artifact and semantic hashes are pinned in the frontmatter.
