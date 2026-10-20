# Independent audit: actual-source conductor and collision multiplicity

20 October 2026. **PASS within the explicit monic, primitive, actual complex
Keller-source hypotheses.** No mathematical correction is required for
Theorems A–D as stated. Two proof-clarity improvements requested during
review were incorporated and checked in the final frozen version.
This verdict does not prove existence of a noninvertible Keller map,
exclude every such map, or solve the complex Jacobian conjecture.

## Reviewed sources and reproducibility

The primary reviewed proof is
`session13-trace-lattice/ACTUAL_SOURCE_CONDUCTOR_AND_COLLISION_MULTIPLICITY.md`.
Its frozen SHA-256 is:

`45494620a95261528932376918d078be3d0e8bf7c3ab9a31f9a5643f8f43bd84`

The initial reviewed proof hash was
`5e9611a16252d60aa673e0fb0b112bf777b5a9f24bfb080b9b86bb8f1de4e25f`.
The final revision adds characteristic-polynomial specialization, explicit
scalar bookkeeping, and publication-relative predecessor links. The
mathematical statements and checker are unchanged.

The reviewed checker is `session13-trace-lattice/verify_conductor.py`:

`6904185214088b4973c85081dd4aa194397bfc9a5661dea43902f267f6b2e69c`

The producer's output `session13-trace-lattice/results.json`:

`76ba0bbdb8f38296ce17d247804d6ae3254f9570eebaab16a107ff94baa7ea09`

The original checker was copied, without edits, into `replay/`. Both

    python3 replay/verify_conductor.py
    python3 -O replay/verify_conductor.py

passed **293 checks over 30 parameterized source pairs**. Their JSON outputs
are byte-identical to each other and to the reviewed producer output.

The independent checker is `independent_audit.py`. Both

    python3 independent_audit.py
    python3 -O independent_audit.py

passed **125 checks**, including **14 curve families**, and produced
byte-identical JSON. It does not import the producer checker. In particular:

- Collision lengths are computed from standard monomials of a Groebner basis,
  not by reading the degree of the resultant.
- Norms and Fitting determinants are computed with multiplication matrices.
- For seven coprime-leading-degree curves, delta is computed by counting
  missing leading degrees of a free power basis, not by taking the degree of
  the inclusion determinant.
- A sharp failure of image/conductor specialization is tested explicitly.

These are finite exact controls, not exhaustive evidence for hypothetical
noninvertible Keller maps. The universal claims are justified by the proofs
audited below.

## Inheritance and adversarial posture

Closest proved mechanism: the actual finite graph normalization and the
collision-gap unit in the predecessor source-line/collision packets.

Canonical hostile: the immersed nodal/tacnodal cylinder with nonconstant
Jacobian. Corrected near miss: singular support is not a divisor with its
contact multiplicities, and a unit gap does not imply properness.

Least-used sidecar retained here: the multiplication on the actual
`A`-orders, through the inclusion matrix, trace dual, and collision
presentation. A bare trace form or determinant does not establish
realizability.

Concept board: finite normalization; characteristic-polynomial
specialization; trace-dual orders; collision contact length; conductor
equalizers. The independent negative probes target the information lost in
each passage.

## 1. Finite normalization and the two free orders

Put `A=C[x,u]`, with `u=P(x,y)`, and `B=C[x,y]`.

1. Monicity of `P` of degree `n>=1` in `y` makes `B` a free `A`-module with
   basis `1,y,...,y^(n-1)`. There is no hidden generic-freeness assertion.
2. `S=A[Q]` is integral over `A`. The assumed equality
   `C(x,P,Q)=C(x,y)` gives the same fraction field as `B`.
3. The minimal polynomial of `Q` has coefficients in `A`: `A` is integrally
   closed and `Q` is integral. Its degree over `Frac(A)` is exactly `n`.
   Thus its monic polynomial is the primitive hypersurface equation `R`
   with the stated variable ordering, and `S` is free over `A` on
   `1,Q,...,Q^(n-1)`.
4. Since `B` is integral over `S`, normal, and has its fraction field, it is
   exactly the normalization of `S`. An element integral over `S` is also
   integral over `B`, so no larger normalization has been overlooked.

All these arguments use the actual source inclusion, not an abstract
polynomial that merely resembles an inverse equation.

## 2. Every source slice remains primitive

For each `a`, monicity keeps `P_a` nonconstant of degree `n` and makes the
slice parametrization finite. Keller's condition makes
`(P_a'(y),Q_a'(y))` nowhere zero.

The normalization of its image is `A^1`: a finite dominant map from `A^1`
extends to projective completions; the target completion has genus zero,
and the complement can have only one point because all its preimages lie
at infinity in the source. The induced self-map of `A^1` is polynomial.
If it had degree greater than one, its derivative would have a zero, and
the chain rule would contradict immersion. Hence the original
parametrization is birational.

Consequently `C(P_a,Q_a)=C(y)` and
`[C(P_a,Q_a):C(P_a)]=n` for every `a`, not just generically. Each individual
normalization branch is smooth. A unibranch singularity would violate the
analytic inverse-function theorem applied to a component with nonzero
derivative.

### Incorporated algebraic hardening of the reduced-section step

Let `T_Q` denote multiplication by `Q` on the free `A`-module `B`. Global
primitivity implies

    R(u,v,x)=det(v I-T_Q).

Specializing the matrix at `x=a` commutes with the determinant. The slice
primitivity just proved says that `Q_a` has minimal polynomial of degree
`n`; its characteristic polynomial has the same degree and is therefore
that minimal polynomial. Thus `R(u,v,a)` is monic and reduced irreducible
directly. The final reviewed note uses this argument. It replaces the
valid analytic local-graph argument cited in the original draft, rather
than repairing a false conclusion.

In particular, the subsequent conductor calculation for a slice does not
assume general conductor base change.

## 3. Trace-dual conductor identity

For a monogenic order `A[theta]` with monic separable equation `f`,
the pairing

    (g,h) -> Tr(g(theta)h(theta)/f'(theta))

is perfect over `A`. Lagrange interpolation identifies the trace with the
coefficient of degree `n-1` in the remainder modulo `f`. On the power
basis, the matrix has zeros below its anti-diagonal and ones on the
anti-diagonal, so its determinant is a unit. This proves the integral
lattice formula, not just its generic-vector-space version.

The differentiated graph identities give exactly

    R_v(P,Q,x)=c P_y,
    R_u(P,Q,x)=-c Q_y,
    c=R_Z(P,Q,x)/J.

The sign is correct for the convention
`J=P_x Q_y-P_y Q_x`. Since `J` is a nonzero constant, `c` is a polynomial.
Separability and `P_y!=0` show `c!=0`.

Using the perfect dualities,

    bB subset S
    iff b S^dagger subset B^dagger
    iff (b/c) S subset B
    iff b/c belongs to B.

The final implication uses `1 in S` in one direction and `S subset B` in
the other. Thus the conductor is exactly `cB`, including as an ideal of
`S`. Repeating this calculation over `C[u]` for each slice proves the
specialized claim independently and shows `c(a,y)` is never the zero
polynomial. There is no invalid flatness or conductor-base-change step.

## 4. D-resultant, scalar normalization, and Fitting ideal

The generic roots `t_i` of `P(x,t)-P(x,y)` other than `y` give

    R_v(P,Q,x)/P_y
      = product_i (Q(x,t_i)-Q(x,y))/(t_i-y).

Both numerator and denominator change sign in each factor if the roots
are written in the opposite order, so no residual sign remains.
Monicity of `U` identifies this product with `Res_t(U,V)` exactly.
The identity extends polynomially through all nongeneric points.

`E=B[t]/(U)` is free of rank `n-1`. Multiplication by `V` gives a square
presentation of its cokernel, exactly the collision module. Its
determinant is the resultant, so

    Fitt_0^B(M_col)=(c).

No radical, support-only replacement, or annihilator/Fitting
identification is used. The `n=1` case has `U=1`, zero module, empty
determinant `1`, and `c=1` in this normalization.

The preceding explicit gap identity makes `t-y` invertible in the full
collision algebra, including its nilpotents. It is therefore legitimate
to interpret its scheme length as ordered distinct collisions with
contact multiplicity under the Keller hypotheses.

### Incorporated scalar clarification

All displayed exact equalities are valid as written because `R` is
monic in `v`. If instead `R'=lambda R`, then

    c'=lambda c,
    Res(U,V)=c=lambda^(-1)c',
    Norm(c')=lambda^n I^2.

The conductor and Fitting ideals do not change. It would be incorrect to
retain the literal equation `Norm(c')=I^2` without the monic
normalization. The theorems retain that normalization, and the final
reviewed preface now states this rescaling rule explicitly.

## 5. Norm square and exact delta length

The inclusion matrix satisfies
`Gram_S=M^t Gram_B M`. The discriminant/norm signs for the two
rank-`n` power bases cancel, yielding exactly

    Norm_(B/A)(c)=det(M)^2.

The field trace, the two ranks, and the derivative norms all refer to
the same extension. There is no unnoticed norm between different
fields.

After specializing, `M(a,u)` is still a full-rank inclusion matrix by
slice primitivity. Smith normal form over `C[u]` gives

    dim_C(B_a/S_a)=deg_u det(M(a,u)).

For any nonzero `h(y)`, each linear factor `y-r` has norm
`(-1)^n(P_a(r)-u)`. Therefore the norm has `u`-degree `deg_y h`,
including constants. Combining this fact with the norm square yields
the first equality in (9).

Independently, the specialized multiplication-by-`V` presentation is
square over `C[y]` with nonzero determinant `c(a,y)`. Smith normal form
gives the cokernel length as `deg_y c(a,y)`. This proves the second
equality:

    collision scheme length = deg_y c(a,y) = 2 delta_a.

For smooth local branches of a reduced plane curve, the conductor
component on branch `i` is the restriction of the product of the other
branch equations. Its order is the sum of the pairwise intersection
multiplicities. This proves formula (11) and includes tangencies, not
just nodes.

## 6. The specialization gates and degree-two classification

The external one-line injectivity theorem is used only to require a
nonempty collision on every slice of a noninvertible Keller map. The
strict lower bound `deg c_a>=2`, the nonvanishing leading-coefficient
conditions, and the even degree-drop conditions follow correctly.
Without the noninvertibility hypothesis, a degree-zero specialization
instead forces an automorphism by that theorem.

In the generic degree-two sector, every slice has degree two and delta
one. Smoothness of the individual branches rules out a cusp. The local
conductor formula allows exactly two branches with intersection
multiplicity one at a single singularity. Thus the two conductor roots
are distinct on every slice.

Since neither the leading coefficient nor the discriminant of the
quadratic can vanish anywhere on `A^1_C`, both are nonzero constants.
Completing the square and constant rescaling are therefore polynomial
coordinate changes preserving `x`, and make `(c)=(y^2-1)`.

The two roots map to the same unique node on every slice, giving
polynomial identities `P(x,1)=P(x,-1)` and `Q(x,1)=Q(x,-1)`.
Consequently `S` lies in the equalizer. Conversely, the conductor gives
`(y^2-1)B subset S` and `C[x] subset S`. Polynomial division proves the
equalizer is precisely

    C[x]+(y^2-1)B.

Finally, `w=y^2-1` and `z=y(y^2-1)` generate that ring over `C[x]`,
and satisfy `z^2=w^2(w+1)`. The relation is irreducible and is the
entire kernel. The asserted nodal-cylinder classification follows.

This is a reduction, not an exclusion. Membership of `P,Q` in the
equalizer does not by itself imply the full generator equality.
Neither condition alone constructs a constant-Jacobian pair.

## 7. Independent hostiles and what they prove

### Explicit requested node

For

    P=y^2-1, Q=x+y(y^2-1),
    R=(v-Z)^2-u^2(u+1),

the conductor is `y^2-1`, the collision length is two, and the two branch
tangent vectors have nonzero determinant. But `J=-2y` and
`R_Z(P,Q,x)=-2y(y^2-1)`. This confirms the nodal ring while rejecting
the attempted identification of the raw derivative with the conductor
outside the Keller setting.

### Higher contact and ramification

For `P=y^2`, `Q=x+y h(y^2)`, the order is

    C[u] + y h(u) C[u] subset C[u] + y C[u].

Direct membership and multiplication by `y` show that its conductor is
`h(y^2)B`. Directly `delta=deg h`, and eliminating `t` with `t+y=0`
makes the collision algebra `C[y]/(h(y^2))`. Its length is
`2 deg h`. Repeated roots away from zero give higher contact at
distinct parameters; roots at zero introduce diagonal/ramification
contributions. Thus identical numerical parity can represent different
geometry unless immersion is retained.

### Genuine base-change failure without the slice hypotheses

For `P=y^2`, `Q=xy`, the global order is
`C[x,y^2,xy]`, with conductor `xB` and monic equation
`v^2-x^2u`. At `x=0`:

- The abstract tensor-product ring is `C[u,v]/(v^2)`.
- Its image in `C[y]` is `C[y^2]`.
- The parametrization is no longer birational.
- The conductor polynomial specializes to zero.
- The collision algebra is not finite-dimensional.

This is exactly the obstruction avoided by proving slice primitivity
before using specialized trace-dual and characteristic-polynomial
formulas.

## 8. Primary sources and final scope

The primary text of Gwoździewicz,
[Injectivity on one line, Theorem 1.1](https://arxiv.org/pdf/alg-geom/9305008),
was inspected directly and matches the claimed use over `C`.

The primary text of van den Essen–Yu,
[The D-resultant, singularities and the degree of unfaithfulness](https://hub.hku.hk/bitstream/10722/156137/1/content.pdf),
was inspected directly. Section 1 defines the same divided-difference
resultant; Proposition 1.1 gives the derivative identities with the
raw-resultant sign, removed here by choosing the monic minimal
polynomial. The note's attribution to the classical D-resultant and to
Abhyankar's Taylor resultant is appropriate.

The Campbell–van den Essen DOI did not yield accessible full text
during this audit. No uninspected assertion from that paper is needed
for the proof verdict.

**Final verdict:** the conductor/Fitting/resultant identities, exact
norm square, all-slice `2 delta` length formula, and degree-two
nodal-cylinder classification pass adversarial review in the stated
scope. The remaining existence/exclusion problems are correctly left
open.
