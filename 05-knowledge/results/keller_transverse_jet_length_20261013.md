# Transverse jets retain the entire Keller cokernel length

2026-10-13. **PROVED:** the elementary jet identity and its scoped Keller
consequences. **FINITE-EXACT:** nonconstant-connection and residue controls.
Self-audited; no priority, peer audit or planar-Jacobian solution claimed.

## 1. Inheritance and the proposed length comparison

[THM-4620](../../01-canon/theorems/THM-4620-keller-dual-field-module-kernel-ladder-and-de-rham-cokernel.md)
identifies `coker D_f` with relative de Rham cohomology.
[THM-4626](../../01-canon/theorems/THM-4626-surjective-dual-fields-are-locally-nilpotent-and-hyp-9255-is-the-planar-jacobian-conjecture.md)
already proves the stronger all-dimension criterion: generic rank zero
forces a dual field to be locally nilpotent and its cokernel to vanish.
We do not reprove a weaker finite-jet version of that theorem.

The supplied Selmer manuscript's Proposition 3.7 uses a first-order lifting
obstruction to bound a specified Artin-ring cohomology length independently
of depth; its other construction seeks a length growing with depth.
The natural new question is whether ordinary thickening of a Keller
component can supply a similar upper bound. It cannot, for an exact reason.

Closest mechanism: the new transverse Weyl action on the cokernel.
Hostile: `dz/z` on a punctured line. Corrected near miss: all formal point
germs are integrable although a global residue survives. Least-used sidecar:
the [Pruefer versus finite-Jordan distinction, THM-3477](../../01-canon/theorems/THM-3477-jacobian-effectivity-moving-numerator-jordan-prufer-response-normal-form.md).

## 2. An elementary canonical jet lemma

Let `M` be a vector space over a characteristic-zero field, with operators
`t` and `nabla` satisfying `[nabla,t]=1`. There is no finite-generation
assumption. For every integer `b>=1`,

\[
 J_b:M/t^bM\longrightarrow (M/tM)\otimes_K K[t]/(t^b),\qquad
 [m]\longmapsto\sum_{j=0}^{b-1}\frac{[\nabla^j m]}{j!}t^j             \tag{1}
\]

is a canonical `K[t]/(t^b)`-module isomorphism.

**Proof.** Differentiating `t^b u` fewer than `b` times leaves a multiple
of `t`, so every coefficient is well defined. The commutator gives
`[nabla^j(tm)]/j! = [nabla^(j-1)m]/(j-1)!`, so (1) intertwines
multiplication by `t` with the truncated shift. For injectivity, a zero
constant coefficient writes `m=tu`; the remaining zero coefficients reduce
to the depth-`b-1` statement for `u`. For surjectivity, lift the desired
constant coefficient arbitrarily, subtract its entire jet, and lift the
remaining shifted depth-`b-1` jet inductively. The case `b=1` is the identity.

In particular, whenever `d=dim_K(M/tM)` is finite,

\[
                       \dim_K M/t^bM=bd.                \tag{2}
\]

The operator `nabla` generally **does not descend** to one fixed quotient
`M/t^bM`; equation (1) does. Thus this is not a finite-dimensional
representation of the Weyl relation, which would contradict its trace.

## 3. Apply the lemma to a marked planar Keller pair

Normalize `Jac(f,g)=1`, put `A=K[x,y]`, `D_f={f,.}`, and

\[
 M=A/D_fA,\qquad t=f-a,\qquad \nabla=-D_g.
\]

Here `M` is an additive quotient, not a quotient ring. Multiplication by
`f` and `-D_g` preserve `D_fA`; the dual fields commute, and
`[-D_g,f]=1`. Hence (1) applies. Since `D_f(f)=0`, ordinary right exactness
also identifies

\[
 M/t^bM=\operatorname{coker}\bigl(D_f:A/t^bA\longrightarrow A/t^bA\bigr).
                                                                  \tag{3}
\]

At an ordinary base value, where the Gauss–Manin module is a finite-rank
connection of rank `mu`, equations (2)-(3) give length **exactly `b mu`**.
If the pair is noninvertible, THM-4620/4626 gives `mu>0`. Merely increasing
transverse precision therefore never gives a depth-independent upper
bound, or even a sublinear one, on this same module.

There is a sharper fixed-component consequence. Every other polynomial
mate has the form `g'=g+h(f)` by the inherited planar centralizer theorem.
On `M`,

\[
 -D_{g'}=-D_g-h'(f)D_f=-D_g.                            \tag{4}
\]

Changing the mate by **any** polynomial target shear leaves not only the
length but the transverse connection and every map (1) unchanged.
This precisely identifies a saturated changing-generator route. Changing
the component `f` or the admissible global complex is a different operation
and needs its own comparison preserving the Keller/source predicate.

## 4. Why local success and unchanged rank do not pay the obstruction

On the explicitly different source `K[t,z,z^-1]`, take `D=d/dz`.
The cokernel is `K[t]` times the residue class of `z^-1` and
`nabla=d/dt`. Its transverse quotients have lengths `b`. Around any
point `z=c!=0`, however, the completed local coordinate gives a formal
primitive `log(1+(z-c)/c)` for `1/z`, so the local completed cokernel is
zero. Passing to point germs destroys the global puncture sidecar.
This is a control on a punctured source, not a counterexample on `A2`.

Conversely a nonzero Pruefer module
`K[t,t^-1]/K[t]` satisfies `M/t^bM=0` for every `b`: multiplication by
`t` is surjective. THM-3477 already records this phenomenon and its actual
transition arrows. Therefore a zero test at one bad support point cannot
replace the **ordinary generic** rank-zero criterion.

The [Smith-jet theorem, THM-2960](../../01-canon/theorems/THM-2960-local-smith-jet-fitting-barcode-and-negative-depth-chamber-atlas.md)
already supplies the finite-presentation comparison:
`length(M/t^bM)=rho*b+sum min(b,alpha_i)` for free rank `rho` and finite
Smith bars `alpha_i`. Its eventual slope is the persistent free rank.
The present Weyl identity explains why the ordinary Keller carrier has
pure linear growth, rather than a bounded exceptional bar of the sort
the Selmer obstruction seeks.

## 5. Rerouted target and checks

**Later continuation:** [the inverse transverse clock](keller_inverse_clock_20261015f.md)
gives an explicit source-polynomial inverse to the transverse operator on
the cokernel, a factorial representative formula and its exact return
cost. It also identifies special-fibre torsion heads and their Pruefer
arms; the inverse clock can leave that torsion. Neither analytic decay nor
the absence of an entirely torsion cokernel implies cokernel vanishing.

A useful next target must change the admissible global carrier while
retaining the exact comparison maps. One can ask for a source-preserving
map between two component cokernels with a certified kernel/cokernel
defect and an independently decreasing invariant. An invertible basis
change, a polynomial mate shear, or the same module at higher jet depth
cannot supply that decrease. The previously excluded nonlinear cotangent
carriers remain excluded; this argument does not reopen them.

The concurrent [dual-field inverse compiler](jc_bilinear_reconstruction_20261013a.md)
provides a concrete positive alternative: an exactly closed finite-dimensional
space of **actual source polynomials**, containing the source coordinates,
yields an explicit inverse. This is not finite length of `M/t^bM`.
It is finite closure under all dual derivations before taking a quotient.
The source inclusion is the distinguishing coordinate. Its hostile control
also closes a tempting detour: characteristic-`p` dual fields always have
`D^p=0`, so increasing the prime and the tested depth together gives no
characteristic-zero closure certificate.

There is a second useful comparison with the concurrent
[Selmer growth formula](bsd_selmer_growth_20261013a.md): after the exact
torsion correction it reads `u_b=b*s+sum min(b,e_i)`, the same elementary
free-rank/finite-bar shape as the inherited Smith formula. Adjacent *complete*
Selmer levels can bound `s`; the ordinary Weyl carrier here has no finite
bars and slope `mu` from its very first layer. The shared invariant is the
growth of fully specified modules, not an identification of the modules.

[Code](../../04-computation/experiments/marked_switches_20261013.py) and
[output](marked_switches_20261013.json) check the exact jet map for the
nonconstant rank-two connection
`d/dt + [[t,1],[1,t^2]]` at depths one through eight, its multiplication
shift and quotient well-definedness; Laurent residue primitives; mixed
Smith bars; and explicit preimages on a divisible Pruefer arm.

Connection contract: a Keller pair maps to its marked global cokernel
and then to (1). The map preserves every transverse jet and its length;
it discards divisible support and does not preserve the global source
under completion at an individual point. The required sidecars are the
ordinary base locus, punctures/residues, and the actual transition maps.
