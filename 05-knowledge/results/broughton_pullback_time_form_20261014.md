# Broughton pullback and the toy clutch's logarithmic time form

**Date:** 2026-10-14. **Status:** Exact synthesis; classical methods credited. This supplements the frozen toy-clutch checkpoint without changing it. No priority claim is made.

## Exact relationship, rather than analogy

Let `F(U,Y)=U+U^2Y`, the classical Broughton-type polynomial. For the critical-free toy clutch, put `t=2z+1`. Then

`P_*(x,(t-1)/2)-8 = F(Phi(x,t))`,

`Phi(x,t)=(U,Y)=(xt-8,(t-1)^2/16)`,

`Jac_(x,t)(Phi)=t(t-1)/8`.

Thus the toy polynomial is an explicit polynomial pullback of the classical example. The map `Phi` is generically two-to-one, and is not a polynomial coordinate change. On the hyperbola `xt=8`, it maps to the line `U=0` by the degree-two polynomial `Y=(t-1)^2/16`; its domain coordinate is Laurent (`t != 0`).

For a smooth planar fibre, write the Hamiltonian derivation as `D_F=F_U partial_Y-F_Y partial_U`, and let its time form satisfy `omega_F(D_F)=1`. Under a polynomial map `Phi` with determinant `J`, the chain rule gives

`D(F composed Phi) (H composed Phi)=J (D_F H) composed Phi`.

Consequently, on smooth fibres wherever `J != 0`,

`omega_(F composed Phi)=Phi^*omega_F/J`.

The identity extends meromorphically to the completed fibre. This is a pullback identity divided by the Jacobian determinant, not an ordinary pullback identity for time forms.

On the base component `U=0`, the Broughton time form is simply `dY`, which is exact. On `xt=8`, the transformed time form becomes

`((t-1)dt/8)/(t(t-1)/8)=dt/t`.

The zero at `t=1` cancels, leaving a nonzero logarithmic residue at the puncture `t=0`. In the original `(x,z)` volume coordinates the time form is `dt/(2t)`, because `t=2z+1`. This explains why the toy clutch retains a period obstruction even after its finite critical point has escaped.

## An independent generic-fibre check

For the whole half-slope family `B=b+x/2`, put

`U=xt+2b`, `r=U(t-1)`, `q=P+2b`.

Then `q=U+r^2/16`, and exact differentiation gives

`Jac_(x,t)(P,r)=U+r=q+r-r^2/16`.

For generic fixed `q`, the fibre is rationally parametrized by `r`:

`U=q-r^2/16`, `t=(U+r)/U`, `x=U(U-2b)/(U+r)`.

These formulas are used only off their finitely many excluded parameter values. The time form in the original coordinates is

`omega=dr/[2(q+r-r^2/16)]`.

When `q != -4`, its denominator has distinct roots `8 +/- 4 sqrt(q+4)`, with nonzero residues `-/+ 1/sqrt(q+4)`. Thus the generic-fibre rational time form has no rational primitive. This is an alternative all-half-slope no-mate proof, including the critical-free exception. It does not rely only on the topology of a special fibre.

## Literature boundary and primary sources

Broughton's original source is S. A. Broughton, *On the topology of polynomial hypersurfaces*, Proceedings of Symposia in Pure Mathematics 40, Part 1 (1983), 167–178, DOI [10.1090/pspum/040.1/713056](https://doi.org/10.1090/pspum/040.1/713056). His later article is *Milnor numbers and the topology of polynomial hypersurfaces*, Inventiones Mathematicae 92 (1988), 217–241; [institutional publication record](https://scholar.rose-hulman.edu/math_fac/7/).

The Hamiltonian time-form method and its period/residue obstruction are established tools, not a new general principle claimed here. A directly relevant primary treatment is A. Bustinduy, L. Giraldo, and J. Muciño-Raymundo, *Jacobian mates for non-singular polynomial maps in C^n with one-dimensional fibers*, Journal of Singularities 9 (2014), 27–42, DOI [10.5427/jsing.2014.9b](https://journalofsingularities.com/volume9/article2.html). Its [full paper](https://journalofsingularities.com/volume9/bustiunduy-giraldo-mucino-raymundo.pdf) develops time forms in Section 2, states the period/residue primitive criterion in Remark 4, and treats the sign-equivalent Broughton polynomial `z1-z1^2 z2` in Section 3.5, Example 1. That example also warns that fibrewise zero periods alone need not produce a globally polynomial mate: reducible-fibre extension can fail.

The contribution of this supplement is the exact `Phi` above and the explicit residue calculation for this clutch architecture. Residue cancellation at ramification and newly exposed punctures give a concrete test for further polynomial pullbacks.
