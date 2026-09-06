# Independent audit of vertical torsion and the response connection

**Status: INDEPENDENT ANALYTIC/SOURCE AUDIT PASS + EXACT FROZEN REPLAY.**
This audit accepts the universal arguments in
[the primary note](planar_jc48_sep06_torsion.md) under its explicit
hypotheses: an algebraically closed characteristic-zero field, a
polynomial `P` with unimodular gradient, and rational constant field
`ker D=k(P)`. It does not infer a Keller inverse or any chart entry.

## 1. The complete principal-part classification

The automatic vertical-pole argument is correct. A polynomial derivation
preserves each localization `R_(g)`. If `g` is a polar prime and `g` does
not divide `Dg`, the numerator of the differentiated highest pole is a
unit modulo `g`, since the characteristic-zero integer pole order is
nonzero. Thus `Df` cannot be polynomial. Once `g` is invariant, the
induced derivation on its curve function field is nonzero: otherwise
both components of the unimodular gradient vanish modulo `g`. Its
constant field is the algebraically closed ground field. One way to
check this last step is that a transcendental constant would make the
curve function field a finite separable extension of a field killed by
the derivation, forcing the whole derivation to vanish. Consequently
`P` is constant on that prime, so every pole is vertical.

Smoothness makes each fiber reduced, with disjoint irreducible components.
At each component, `u=P-c` is a uniformizer. The derivative of `u^m f`
has zero residue at a pole of order `m`, so its residue is a scalar.
Subtracting that scalar times `u^-m` and iterating gives the claimed
unique principal part. This procedure uses the specified function `u`;
it needs no choice of a coefficient-field splitting or local coordinate.

The class map is well-defined on the polynomial response quotient.
Changing a polynomial representative adds a polynomial primitive;
changing a rational primitive adds a rational function of `P` by the
stated constant-field hypothesis. The latter changes every component
over a fixed value by the same principal part. Conversely, a diagonal
tuple can be removed by a finite sum of rational functions of `P`.
The remaining primitive has no affine polar prime and is polynomial by
normality of `k[x,y]`. This proves the exact kernel, including the
diagonal quotient rather than a choice of one distinguished component.

The CRT surjectivity argument also passes. The ideals `(f_i^n)` are
pairwise comaximal, and the polynomial `u^n L_i(u)` is a legitimate
residue prescription. Dividing the resulting numerator by `u^n` gives
the requested principal part at every component because the other
factors of `u` are units in its local DVR. Each ideal is `D`-stable,
so differentiating the congruences makes `Dh` divisible by their
product `u^n`. The response is polynomial globally. Finite sums handle
different target values without introducing additional poles.

It follows that the full torsion is exactly the direct sum of
`(k^{I_c}/k*1) tensor u^-1 k[u^-1]`. The number of arms is `r_c-1`.
Multiplication by `P-c` is invertible on every other primary summand,
so the formula `dim ker(P-c)^n=n(r_c-1)` holds in the entire response
module, not merely inside a selected cyclic class. No implicit bound
on the set of reducible fibers or on a pole order is used.

## 2. Canonical connection, signs, and exact heights

With `D=(-P_y,P_x)`, the volume convention gives
`i_D(dx wedge dy)=-dP`. Applying a Bézout vector field `V(P)=1`
therefore gives `[V,D]=-(div V)D`, with the sign in the primary note.
Thus `(V+div V)Dh=D(Vh)` and the operator descends to `C_P`.

Two Bézout fields differ by `hD`, because the gradient row is unimodular.
Since `div D=0`, their two response operators differ by
`hDa+(Dh)a=D(ha)`. This verifies choice independence on the whole
quotient. The commutator with multiplication by `P` is the identity.
It is an algebraic connection, not a `k[P]`-linear endomorphism.

For a rational primitive `f`, the primitive of the connected response is
`Vf`. The vector field preserves each local DVR even though it need not
preserve its maximal ideal. Hence it creates no pole from the regular
remainder. Since `V(u)=1` and the principal coefficients are scalars,
the componentwise principal part is exactly the ordinary derivative.
This proves the canonical connection formula without choosing global
CRT selectors compatible with differentiation.

A nonzero highest coefficient in the component quotient stays nonzero
after multiplication by its negative integer pole order. Each derivative
therefore raises the exact primary order by one, with unchanged primary
support. This verifies both the annihilator exponent formula and the
linear independence of the derivative sequence. It also excludes every
nonzero finite-dimensional or finitely generated torsion submodule that
is stable under the connection.

For `theta=[1]` and `mu=[div V]`, the identity `mu=nabla theta` is direct.
The converse implication in `theta=0 iff mu=0` is independently checked:
if `div V=Dh`, then

```
omega=A dy-B dx+h dP,
d omega=(div V-Dh) dx wedge dy=0,
omega(D)=A P_x+B P_y=1.
```

The polynomial Poincaré lemma gives `omega=dQ` and `DQ=1`.
This is existence of a mate for the fixed `P`, not its polynomial
invertibility. The proof correctly avoids inferring that `theta` is
torsion merely because its derivative is torsion.

## 3. Controls and classical overlap

The displayed one-root family has the stated Bézout field and rational
primitive `Q0=x^(1-r)/(lambda*(r-1))`, with `DQ0=1`. Its top coefficient
relative to `u=P` is nonzero on the component `x=0` and zero on the
other component, giving exact order `r-1`. The equal-value example
`P=x^2+(x^2-1)^2 y` is smooth and has exactly three disjoint components
at value one, hence two arms. The other fibers are irreducible by
primitivity of the linear polynomial in `y`. In both examples
`k(x,y)=k(P)(x)` and the induced derivation is a nonzero rational
multiple of `partial_x`, verifying the rational constant-field hypothesis.

The singular hostile `P=xy` has **nonzero** `theta`: the derivation acts
on `x^i y^j` with eigenvalue `j-i`, so no constant is in its image.
The draft's opposite word was flagged before this audit and was repaired
by the author; the current statement was checked. This example is excluded by the smoothness
hypothesis; its intersecting component ideals show why the CRT step
cannot be transported to a singular fiber.

The literature comparison was checked against the primary text of
[Bonnet, Relative exactness modulo a polynomial map](https://arxiv.org/html/math/0602223).
Proposition 1.2 makes topological relative-exactness classes torsion;
Theorem 1.5 gives vanishing under its quasi-fibered hypotheses.
The introduction explicitly credits Gavrilov's connected reduced plane
case and a comprehensive plane torsion analysis by Bonnet--Dimca,
*Relative differential forms and complex polynomials*, Bulletin des
Sciences Mathématiques 124 (2000), 557–571. No priority claim for the
general torsion phenomenon or connected-fiber consequence is justified.
The earlier paper has not been compared theorem by theorem in this audit.

The relation to that literature is mathematically typed: under the unit
gradient, contraction with `D` identifies
`Omega^1/(R dP+dR)` with `R/D(R)`. Surjectivity uses a polynomial
Bézout form evaluating to one on `D`, and its kernel consists of
multiples of `dP`. This does not replace the primary note's direct
component-jet and connection proof.

## 4. Source acceptance

The full [standalone source](../../04-computation/planar_jc48_sep06_torsion.py)
was independently read. It imports SymPy, not another research producer.
The four one-arm families have `r=2,3,4,5`, with derivative orders zero
through three. Their rational primitives, polynomial responses, exact
nonzero highest principal coefficients, and polynomial annihilating
witnesses are checked separately. Three polynomial controls per family
check the Hamiltonian commutator, connection commutator, and Bézout gauge
identity including its divergence term.

For the three-component fiber, the source constructs both line-component
selectors through pole order five, retaining zero principal parts on the
third component. The confluent congruences verify the entire principal
jet, not just its leading coefficient. Ordinary plane-ring differentiation
checks that their responses are polynomial. It verifies connection
transitions on orders one through four modulo an explicit polynomial,
and module transitions on orders two through four. The determinant with
the diagonal column certifies that the two component directions remain
independent over the coincident target value. These tests supplement the
universal CRT argument; they are not an extrapolation to other fibers.

The coordinate controls have genuine polynomial primitives in the
`(P,y)` monomial basis. The singular `xy` monomial calculation and
noncomaximal-ideal check are correctly kept outside the smooth theorem.
All gates remain active under optimization.

Independent replay commands:

```bash
python3 -B 04-computation/planar_jc48_sep06_torsion.py > /tmp/planar_jc48_torsion_audit_normal.out
python3 -B -O 04-computation/planar_jc48_sep06_torsion.py > /tmp/planar_jc48_torsion_audit_optimized.out
cmp /tmp/planar_jc48_torsion_audit_normal.out /tmp/planar_jc48_torsion_audit_optimized.out
cmp /tmp/planar_jc48_torsion_audit_normal.out 05-knowledge/results/planar_jc48_sep06_torsion.out
```

Both runs pass **278 exact gates**, byte-identical to each other and
to [the author's frozen output](planar_jc48_sep06_torsion.out).

```
source SHA256
3c01beaeca3abf4567de0a31d2939ae7e1858b0b85746986ccf674f6b336332c
frozen and both independent replay outputs SHA256
4798029319d16097ce37fe1db3e1379255f7fe13e74f414b85bf348a2abfdda3
```

The final frozen comparison and source hash recheck passed. No
mathematical correction remains. This audit is frozen and did not edit
the primary source or output.
