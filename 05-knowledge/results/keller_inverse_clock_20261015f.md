# A source-preserving inverse clock on the planar Keller cokernel

2026-10-15. **PROVED:** the explicit inverse, polynomial transport identities,
and compact-set bound below. **FINITE-EXACT:** the declared symbolic controls.
**OPEN:** vanishing of the actual Keller cokernel and JC(2). Self-audited;
no historical novelty claim is made for polynomial de Rham exactness.

## 1. Inheritance: change the operation, not the excluded carrier

The closest proved mechanism is the transverse Weyl action in
[transverse jet length](keller_transverse_jet_length_20261013.md), together
with [THM-4620, Keller dual-field modules](../../01-canon/theorems/THM-4620-keller-dual-field-module-kernel-ladder-and-de-rham-cokernel.md).
Those jets have length exactly `b mu` at an ordinary point. Increasing their
depth, or replacing a mate by `g+h(f)`, does not reduce that obstruction.
The positive alternative already available is the
[finite source-polynomial closure compiler](jc_bilinear_reconstruction_20261013a.md).

The canonical hostile here is a regular-singular Laurent module with
nontrivial Kummer character, not a proposed Keller counterexample. The
corrected near miss is to regard a small representative after repeated
integration as a small representative of the *same* class. The least-used
coordinate is the inverse transverse clock, including the operator needed
to return to the original class.

We also recovered the stronger old exclusion
[THM-3981, centered cusp quadrature](../../01-canon/theorems/THM-3981-centered-cusp-quadrature-genus-two-transcendence.md):
its nonexact time differential remains nonexact on every finite algebraic
cover. Thus the finite-cover idea is not a new escape from that carrier.
[The genus-two transfer-word classification](genus_two_transfer_words_20261011.md)
remains closed in its stated scope. The present construction starts instead
with an actual arbitrary planar Keller pair.

## 2. The transverse operator is bijective on the cokernel

Let `K` have characteristic zero, `A=K[x,y]`, and `Jac(f,g)=1`. Write

\[
 D=f_x\partial_y-f_y\partial_x,\qquad
 E=g_y\partial_x-g_x\partial_y,\qquad M=A/D(A).
\]

This is an additive quotient, not a quotient ring. The actual polynomial
derivations commute, so `E` induces an operator on `M`. Also

\[
 dh=(Eh)\,df+(Dh)\,dg,\qquad df\wedge dg=dx\wedge dy.       \tag{1}
\]

**Proposition J1 (PROVED).** The operator `E:M -> M` is bijective. An
explicit inverse on classes is

\[
 E^{-1}[r]=[R_f r],\qquad
 R_f r=f_x(x,y)\int_0^x r(s,y)\,ds.                       \tag{2}
\]

All integrals in this note are polynomial antiderivatives. In particular
there is no analytic continuation or local replacement of the source.

**Proof.** Set `P=integral_0^x r(s,y) ds`, `a=-g_x P`, and `b=f_x P`.
The polynomial form `P dy=a df+b dg` has differential `r dx wedge dy`.
Consequently

\[
                         E b-D a=r.                     \tag{3}
\]

This proves surjectivity on classes. If `E b=D a`, the form
`a df+b dg` is closed. Every closed polynomial one-form on the affine
plane is a polynomial differential: integrate its `dx` coefficient and
then the remaining function of `y`. Write it as `dh`; (1) gives `b=D h`.
Thus the induced operator is injective. Its inverse is therefore exactly
(2), and (2) is well defined modulo `D(A)`. This also proves, without
assuming a quotient projection commutes with integration, that
`R_f(D(A)) subset D(A)`. QED.

The alternative expression `f_y integral_0^y r(x,s) ds` gives the same
inverse class. Changing the lower `x` anchor changes (2) by `f_x a(y)`,
which equals `D(integral a(y)dy)`. The inverse class is anchor independent,
although a concrete polynomial representative retains its chosen anchor.
Like the inherited transverse connection, it is unchanged by a polynomial
mate shear. Its formula uses the actual component and Cartesian source.

## 3. An exact factorial formula, with its transport cost retained

For `m>=1`, repeated use of (2) gives

\[
 R_f^m r(x,y)=\frac{f_x(x,y)}{(m-1)!}
 \int_0^x \bigl(f(x,y)-f(s,y)\bigr)^{m-1}r(s,y)\,ds.       \tag{4}
\]

Induct and interchange polynomial integrals; the inner integral uses
`d f(s,y)=f_x(s,y) ds`. No inverse of `f` or of `f_x` is used, so critical
zeros of `f_x` cause no exception to this identity.

Over the complex numbers, on a fixed closed polydisc `|x|<=a, |y|<=b`,
let `L=sup |f_x|` and `H=sup |r|`. The iterated Volterra expression over
the ordered integration simplex gives the sharper bound

\[
                 \sup |R_f^m r|\le H\frac{(aL)^m}{m!}.   \tag{5}
\]

The straight radial paths stay in the same polydisc. This bound is on the
displayed representatives; it does not define a separated norm on `M`.

The retained clock says

\[
               [R_f^m r]=E^{-m}[r],\qquad
               [E^m R_f^m r]=[r].                       \tag{6}
\]

There is even a literal polynomial invoice for the latter equality. Put
`A_g h=-g_x integral_0^x h(s,y) ds`. Equation (3) is `ER_f-I=D A_g`.
Telescoping, using `[E,D]=0`, gives

\[
 E^mR_f^m-I
   =D\sum_{j=0}^{m-1}E^j A_g R_f^j.                     \tag{7}
\]

Thus factorial smoothing has not annihilated an obstruction. It has moved
it along an invertible clock. Returning it to the original class requires
the `m` powers of `E` in (6), and those powers cannot be bounded away by
the smallness of (5). Cartesian differentiation normally has factorial
Cauchy costs, exactly the scale that could consume this apparent gain.

## 4. A sharp hostile with regular singularity and finite monodromy

On the explicitly different module `K[z,z^-1]`, take

\[
 t\cdot h=zh,\qquad \nabla=\partial_z+\frac{1}{2z}.
\]

It satisfies `[nabla,t]=1` and

\[
 \nabla z^n=(n+\tfrac12)z^{n-1}\quad(n\in\mathbb Z).
\]

So `nabla` is bijective, yet the module is nonzero. Its inverse orbit of 1 is

\[
 \nabla^{-m}1=\frac{z^m}{\prod_{j=1}^m(j+1/2)}.          \tag{8}
\]

Every term is nonzero, but it tends to zero uniformly on each fixed compact
disc after interpreting the displayed polynomial. The connection has only
regular singularities at 0 and infinity; a horizontal solution `z^-1/2`
has order-two monodromy. Hence replacing an irregular exponential hostile
by a geometric finite-character example does not remove the obstruction.
This is an abstract module control, **not** a realization as `A/D_f(A)`
for a noninvertible Keller map.

There is a useful **conditional** corollary: if `M` is torsion-free over
`K[f]`, then `K(f) tensor M` has no nonzero rational horizontal section.
After clearing denominators such a section `h` satisfies `s(f)h in M`;
successive derivatives give `s^(j)(f)h in M`. These polynomials generate
the unit ideal, so `h in M`, contradicting injectivity of `E`. The
torsion-free condition is additional: THM-4626 excludes a nonzero cokernel
that is *entirely* torsion; it does not exclude torsion inside a module of
positive generic rank. The distinction leads to the following exact repair.

## 5. A special fibre supplies a real Pruefer arm, if it disconnects

Now let `K` be algebraically closed. Put `t=f-a`, and let the smooth fibre
`f=a` have `r_a` connected components. The inherited planar centralizer
gives `ker D=K[f]`. The kernel/cokernel sequence for multiplication by `t`
on `A --D--> A` therefore gives

\[
 M[t]\simeq \ker(D:A/(t)\longrightarrow A/(t))/K
       \simeq K^{r_a-1}.                                \tag{9}
\]

Each smooth affine curve component has only constant functions killed
by its nonzero derivation; `Dg=1` ensures it is nonzero there. This proves
the second identification. In particular, torsion at this point measures
actual disconnectedness, rather than an arbitrary presentation defect.

For any Weyl module with `[E,t]=1`, its `t`-primary torsion has the direct
vector-space decomposition

\[
 T_a=\bigoplus_{j\ge0} E^j M[t],\qquad
 t E^j v=-jE^{j-1}v.                                   \tag{10}
\]

Induct on the least power of `t` killing an element to prove spanning;
applying the highest power of `t` proves independence. Thus it is a sum
of `r_a-1` Pruefer arms as a `K[t]`-module. Multiplication by `t` is
surjective on this torsion, so `T_a/t^b T_a=0` at **every** depth `b`.
This realizes the previously inherited jet-invisible mechanism on the
actual cokernel, conditional on a disconnected special fibre.

J1 makes the inverse clock more informative still. For nonzero `v in M[t]`,
put `w=E^-1 v`. Then

\[
          E^{-m}v=\frac{t^{m-1}}{(m-1)!}w\quad(m\ge1).   \tag{11}
\]

Indeed `[t,E^-m]=m E^{-(m+1)}`. Every term is nonzero since `E` is
bijective. Therefore `w` is not `t`-primary torsion; in fact it is not
`K[f]`-torsion at all. Torsion decomposes into its primary supports, `E`
preserves them, and the head `v` is not in `E(T_a)` by (10). A torsion
preimage is impossible. Equation (11) explains exactly how a torsion
head transfers to a free direction under the inverse clock. The factorial
denominator is real and still does not annihilate that direction.

An exact abstract control is `K[z,z^-1]` with `E` multiplication by `z`
and `t=-partial_z`. Here `E` is bijective, `[E,t]=1`, the `t`-primary
torsion is `K[z]`, its head is `v=1`, and `w=z^-1` generates the free
direction modulo that torsion. Formula (11) becomes `z^-m=z^-m`.
This example shows directly that bijectivity on the whole module does
not give bijectivity on its torsion submodule.

No disconnected fibre of an actual noninvertible Keller pair is supplied
here. Equations (9)--(11) are implications for any such pair and retain
the unresolved existence question. They also explain why the conditional
horizontal-section corollary cannot silently discard torsion.

## 6. What this changes in the proof search

The new usable object is an explicit inverse-clock representative and its
exact image correction (7). A proposed contraction argument can now be
tested against that clock instead of an unspecified primitive. A successful
next step would control the *recovery cost* in (6) using an invariant of
the actual polynomial source, or produce the existing finite dual-field
closure certificate. Mere analytic decay, finite-cover decoding, and
unchanged transverse jets do not supply either condition.

Connection contract: `(f,g,r,m)` maps to `R_f^m r` together with the clock
`m`, return operator `E^m`, and polynomial correction in (7). It preserves
the source class **after return**, not before return. Forgetting the return
operator destroys that predicate. Compare the native signed terminal
transport and discounted debt in
[the reset continuation](collatz_reset_debt_20261015f.md): the similarity
is the need to price transport, not an identification of these operators.

[Exact controls](../../04-computation/experiments/marked_debt_controls_20261015f.py)
check (3), preservation of `D(A)`, and 60 instances of (4) on three actual
Keller maps, 45 literal return invoices (7), 20 inverse-clock Kummer
iterates, and 29 Laurent Weyl/torsion controls. Normal and `-O`
runs produce [the same output](marked_debt_controls_20261015f.json).
