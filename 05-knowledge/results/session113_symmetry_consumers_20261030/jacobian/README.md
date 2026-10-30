# Session 113: finite symmetry certificates and the Jacobian cost of acquiring symmetry

**Status:** proved elementary special-class statements and a sharp obstruction. The planar Jacobian conjecture remains open. No NC2/GMC implication is used. No claim of literature priority is made.

## 1. The operation and its scope

A symmetry can simplify a Keller map only if it is an **actual symmetry of that map**, not merely a symmetry imposed on its coefficients afterward. Here we give a finite, degree-aware certificate for an actual symmetry, and compute exactly why the tempting imposition operation fails.

Work over C (or an algebraically closed field of characteristic zero). A Keller map is F=(P,Q) in C[x,y]^2 with det JF=c in C*, and its degree is max(deg P,deg Q). Fix integers w>=1 and m>=1, a primitive m-th root zeta, and

    g(x,y)=(zeta*x,zeta^(-w)*y).

**Equivariance always means F∘g=g∘F in the displayed coordinates**, with the same source and target action. Neither an unknown nonlinear conjugacy nor equivariance only on a quotient is supplied by our results.

## 2. Finite-to-continuous rigidity theorem

**Theorem.** If F is a g-equivariant Keller map of degree at most D>=1 and

    m > wD+1,

then F=(a x,b y), with ab=c. The strict numerical inequality is uniformly sharp: for every w>=1 and D>=2, at m=wD+1 the nonlinear automorphism

    (x+y^D,y)

is g-equivariant and Keller.

**Proof.** Equivariance says that each monomial x^i y^j of P has weight i-wj congruent to 1 modulo m, and each such monomial of Q has weight congruent to -w. All supported weights lie in [-wD,D]. Their differences from the required weights have absolute value at most

    max(wD+1,D+w)=wD+1,

where (wD+1)-(D+w)=(w-1)(D-1)>=0. Thus there are no nonzero multiples of m among those differences: congruence is equality. The map is equivariant under the entire multiplicative group, not just the finite cyclic group.

The nonnegative integer solutions of i-wj=1 and i-wj=-w give

    P=x A(t), Q=y B(t), t=x^w y,

for A,B in C[t]. Direct differentiation gives

    det JF = AB + w t A' B + t A B'.                 (1)

Neither A nor B is zero, since c is nonzero. If their degrees are a,b with leading coefficients alpha,beta, then the leading coefficient of (1) is

    (1+w a+b) alpha beta

in degree a+b. Characteristic zero makes it nonzero. Constancy therefore forces a=b=0. Finally, when m=wD+1, y^D has weight -wD congruent to 1, proving the sharp example. ∎

For w=1 this reads **D<m-1 implies diagonal linear**, with first possible nonlinear degree exactly m-1. This is an exact cutoff, not an asymptotic estimate.

**Sharper support certificate.** The sufficient degree bound can be replaced by the exact finite test that no supported monomial has the correct weight modulo m but a different integer weight. It may succeed when the coarse degree bound fails. We do not claim that the coarse bound is optimal for each individual pair (m,w).

**Multiple clocks.** Simultaneous equivariance for cyclic subgroups of orders m_1,...,m_k, in the same displayed coordinates and weights, implies equivariance for their generated subgroup of order L=lcm(m_1,...,m_k). Hence L>wD+1 is sufficient. Orders are combined by lcm, not by their product unless coprime. This makes the certificate finite: a bounded-degree map does not need an infinite family of symmetry tests.

**Relation to corrected canon.** The exact-torus endpoint is consistent with `01-canon/theorems/THM-4727-mixed-weight-torus-reduction-of-keller-maps-ramification-exponent-twisted-planar-jacobian-equation-and-the-dimension-profile.md`. Our elementary leading-coefficient calculation avoids importing its higher-dimensional claims. The extra operation here is the finite-order/degree lifting certificate and its matching aliasing boundary. This is not a reduction of arbitrary Keller maps to that class.

## 3. The missing conservation law for Reynolds projection

Write a polynomial H as the sum of its cyclic weight components H_r, r in Z/mZ: retain precisely monomials x^i y^j with i-wj≡r. Put {P,Q}=P_x Q_y-P_y Q_x. Differentiation subtracts source coordinate weights, so

    {P_r,Q_s} has cyclic weight r+s+w-1.

Consequently the Keller equation is equivalent to the complete set of identities

    sum_{r+s ≡ k+1-w} {P_r,Q_s} = c if k=0, else 0.    (2)

These are polynomial identities, not scalar counts. Reynolds projection onto equivariant maps gives

    R(F)=(P_1,Q_{-w}).

It is equivariant, but its determinant retains only one summand of the k=0 identity. Exactly,

    det J R(F) = c - sum_{r+s ≡ 1-w, (r,s)≠(1,-w)} {P_r,Q_s}.  (3)

Thus projection preserves the original Keller constant **if and only if the omitted balanced-character brackets sum to zero**. For preservation of some nonzero constant, the omitted sum must instead be constant and different from c. Merely retaining the right weight, normalizing the derivative at the origin, or checking support cannot replace this obligation.

### Hostile control: an automorphism whose projection is not Keller

Take w=1 and

    F=(x+y², y+(x+y²)²).

It is a composition of two triangular automorphisms, so det JF=1. For every m>=6 its equivariant Reynolds projection is exactly

    R(F)=(x, y+2xy²),
    det J R(F)=1+4xy.

The projected map is even equivariant under the full torus; it still fails to be Keller. The omitted balanced pair is

    {y²,x²}=-4xy,

which precisely compensates the retained determinant defect. The other nontrivial cancellation is {x,y⁴}+{y²,2xy²}=4y³-4y³=0, in a nonzero character channel. The original map and its projection both have derivative identity at the origin: first-order normalization misses the obstruction. At (1,1), the projected determinant is 5, while at the origin it is 1.

This is stronger than just saying averaging is nonlinear: it gives the exact missing algebraic observable, including its sign and character channel. The discarded data contain essential volume cancellation even when the retained map has perfect continuous symmetry.

## 4. What this contributes and what remains missing

The useful interface is:

    actual finite equivariance + bounded degree
      -> no character aliasing -> continuous equivariance
      -> one-variable top-degree obstruction -> linear Keller map.

The failed interface is:

    arbitrary Keller map -> average coefficients -> equivariant Keller map.

Equation (3), not equivariance alone, is the gate for the latter. A route toward the general conjecture would need a Keller-preserving way to manufacture adequate symmetry, or a new theorem controlling the omitted bracket channels. Neither is supplied here, and the two-shear example rules out unconditional Reynolds averaging. The result does not settle JC2, establish a new general degree bound for counterexamples, or transfer a special-image NC2/GMC theorem to arbitrary maps.

## 5. Reproducible verification

Run `python check_symmetry.py` and `python -O check_symmetry.py` (Python 3 + SymPy). Explicit exception gates keep every check active in both modes; publication reruns give byte-identical output. `check_results.json` records the exact run. Controls include:

- every balanced modulus m=3,...,80: no alias through degree m-2, actual alias and nonlinear Keller example at degree m-1;
- 204 weighted no-alias tests and 192 sharp weighted boundary maps;
- direct symbolic formula (1) at weights 1,...,7;
- exact componentwise identities (2), discarded term (3), and failed averaging for all m=6,...,30;
- the origin-vs-(1,1) Jacobian contrast.

The proof establishes all parameter values; these finite exact tests are regression and hostile controls, not a replacement for the proof. No Git operations were performed by this lane.
