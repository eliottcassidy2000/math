# Four different meanings of fourth order, and the continuum boundary

Status: elementary proved distinctions and cited problem scope. These are supporting calculations, not new continuum Yang–Mills results. They do not depend on the attached-paper audits.

## 1. Four notions that must not be identified

1. **Walsh interaction degree:** a product of four independent sign coordinates. Session33's Fano count has precisely degree-four nonconstant Walsh terms.
2. **Field-polynomial degree:** the classical Yang–Mills action has a quartic term in the gauge potential because F_A=dA+g A∧A and the action is quadratic in F_A. The full expansion also has quadratic and cubic pieces. The isolated quartic piece is not the whole locally gauge-invariant action.
3. **Differential operator order:** squaring a second-order kinetic operator produces a fourth-order operator. An inequality involving this square is not an identity with a quartic potential or a fourth Walsh moment.
4. **Spacetime dimension:** four real coordinate dimensions in the Clay problem. A four-cycle, four-face tetrahedral term, or four-variable monomial does not supply a four-dimensional continuum field theory.

There are actual connections among some of these objects, but each needs a map. The Fano rainbow Wilson-loop identity in the companion note supplies one such map between specific finite edge observables and degree-four line characters. It does not identify all four notions.

## 2. Why dimension four is structurally special for the classical action

On R^d set A_lambda(x)=lambda A(lambda x). The derivative and quadratic terms scale together, giving F_(A_lambda)(x)=lambda² F_A(lambda x). Therefore

    integral_R^d |F_(A_lambda)|² dx
       = lambda^(4-d) integral_R^d |F_A|² dx.

Thus d=4 is the scale-invariant classical dimension. This is a dimensional calculation, not a construction of a renormalized quantum measure or proof of a quantum mass scale. In particular the finite-code parameter J/h is not identified with a continuum lattice spacing or a renormalized coupling.

## 3. Quartic flatness does not decide a quantum gap

Consider the one-dimensional operator H_g=−d²/dx²+g x^4 for g>0. Under the unitary coordinate scaling x=g^(−1/6)y, it is equivalent to g^(1/3)H_1. The confining potential tends to infinity, so the usual compact-resolvent argument gives discrete spectrum; the one-dimensional ground eigenvalue is simple and the next eigenvalue is strictly larger. Hence

    gap(H_g)=g^(1/3) gap(H_1)>0.

Yet the potential Hessian at x=0 vanishes. A missing quadratic Taylor term therefore does not imply a missing quantum gap. This is the standard anharmonic-oscillator mechanism, not a new theorem about Yang–Mills.

Conversely on R² take K_g=−partial_x²−partial_y²+g(x−y)^4. Orthogonal coordinates u=(x−y)/sqrt2, v=(x+y)/sqrt2 give

    K_g=(−partial_u²+4g u^4) + (−partial_v²).

The free v-direction makes the spectral threshold continuous: there is no positive spectral interval above the bottom excluded from the spectrum, and no normalizable full ground state. The potential is still nonnegative quartic with vanishing Hessian at the origin. Quartic degree and positivity alone distinguish neither example. Confinement, flat directions and the actual Hilbert-space realization matter.

The finite code model supplies a third distinction: even when every finite system has a unique ground and a positive gap, a parameter limit can make that gap vanish, and an invariant-sector gap can hide that behavior.

## 4. A tournament does not specify a unique simplicial field theory

A tournament gives an orientation on every pair of vertices. Completing every clique to a simplex gives the full simplex, which is contractible; its triangle faces and tetrahedral boundaries obey ordinary simplicial cancellation. Selecting only transitive directed faces instead gives a different complex. For the three-cycle, the full completion fills the triangular disk while the transitive-face complex leaves the boundary circle. These are not two descriptions of the same topology.

A gauge field additionally needs a structure group, edge transport law, face holonomies, an action/measure and a declared gauge action. A vertex algebra additionally needs fields and locality/associativity identities. Neither structure follows from drawing arrows on the one-skeleton alone. The companion winner-algebra example proves exactly how associativity can erase tournament information, while the nontrivial Z2 associator example shows actual arithmetic on simplex labels satisfying higher coherence.

## 5. What a finite simplicial result must still supply for the prize target

The [official Clay description](https://www.claymath.org/millennium/yang-mills-the-maths-gap/) concerns a nontrivial quantum Yang–Mills theory in four dimensions with a positive mass gap, and still lists the problem as unsolved. An exact finite partition function or a finite matrix gap does not establish this.

A serious route must identify the physical gauge-invariant states/observables, construct appropriate continuum and infinite-volume limits with the requisite axioms, and control the spectrum above the actual vacuum in that limiting theory. Uniform estimates in the relevant physical scaling are essential; a fixed finite mesh or a fixed perturbative order is not that uniform statement.

Simplicial discretizations themselves are established mathematics: Christiansen–Halvorsen, [A simplicial gauge theory](https://arxiv.org/abs/1006.2059), proves a discretely gauge-invariant action and consistency for its specified construction. That is a useful primary model for the discretization layer, not a theorem resolving the continuum mass gap.

The productive next question is therefore not whether a repeated four proves a gap. It is which precisely defined representation sector, gauge constraint and limiting estimate our fourth-order construction controls.
