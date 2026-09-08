# Independent audit of base-only quartic pole closure and root gluing

**Status: INDEPENDENT ANALYTIC + EXACT SOURCE + REPLAY AUDIT PASS.**
The [primary note](planar_jc48_sep08_quartic_boundary.md) proves its
base-only Jacobian extension and the resulting global square prefix on
the specified DG surface. Neither is a closure of the general quartic
Keller stratum or JC(2). The proposed mate has unrestricted degree.

## 1. The changed hypothesis is actually paid

I read the full source and the load-bearing canonical passages:
[THM-2129 / quartic-faber-three-coefficient-boundary-classification](../../01-canon/theorems/THM-2129-quartic-faber-three-coefficient-boundary-classification.md),
[THM-2180 / quartic-first-pole-stratum-is-empty](../../01-canon/theorems/THM-2180-quartic-first-pole-stratum-is-empty.md),
[THM-2189 / nonsplit-quartic-deck-forces-the-remaining-pole-congruence](../../01-canon/theorems/THM-2189-nonsplit-quartic-deck-forces-the-remaining-pole-congruence.md),
and [THM-2202 / uniform-all-degree-quartic-pole-closure](../../01-canon/theorems/THM-2202-uniform-all-degree-quartic-pole-closure.md).
The targeted correction check retains MISTAKE-248: the last divisibility
condition is one congruence, not two independent conditions. The later
uniform theorem governs over the older notes' historical open-pole text.

For J(P,Q)=kappa(x)!=0, the coefficient of fibre degree n+3 is unchanged.
Actual constant target shears remove each mate degree divisible by four;
the leading coefficient identity supplies exactly the scalar required.
A degree-zero residual mate would give -P_z Q'(x), which cannot be both
nonzero and fibre-independent. The final positive degree is therefore
not divisible by four, so valuations of the leading coefficient make the
quartic leading coefficient a square. This is an unbounded descent in
the proposed mate degree, not a finite bank assumption.

On the quadratic deck the Jacobian becomes kappa/U. The triangular
Faber expansion has constant coefficients because any nonconstant top
coefficient contributes -4a_d' Z^(d+3). At each successive subtraction
the Jacobian has degree at most two, so the same argument remains valid,
including the last degree-zero coefficient. The constant field on a
chosen deck component is C. Hence the first two flux combinations have
derivative zero. Only the final primitive equation has changed.

I checked the input use in all pole cases. THM-2180 uses only the regular
original fibre boundary and its uniquely deepest nonzero binomial term.
At the second possible pole, the negative-scale boundary/flux triple
uses constancy of the two fluxes; each scale exponent is negative, so
exact vanishing from nonsplit parity is unnecessary. THM-2129 excludes
odd reduced degrees. The cases of zero linear face coefficient are
covered separately by the cubic recurrence or the nonzero even binomial
coefficient; no active face is lost.

For n=2, the three quantities in THM-2189 §5 force C=A and a nonzero
lower linear seed, then valuation(Lambda)=A. The unique term of order
-3A in the second flux cannot cancel with the two terms of order -2A.
Neither the final primitive nor invertibility of kappa at that place is
used. For every higher twice-odd degree, I checked both chambers and
their equality boundaries in THM-2202. Below and on the cusp, the top
first-flux order is negative (use C<2A and rL<=(2r-1)C), so a unique
matching odd seed is indeed required; its second flux is L orders
shallower. Above the cusp, the uniquely required odd boundary seed has
a uniquely negative first flux. The case Lambda=0 and every lower even
or odd seed are included. Again the final primitive equation is absent.
Thus zeros of kappa at the tested place create no missing case.

## 2. Full-chart gluing and its exact consequence

The literal change x=1/r, t=-r^2-r^4b has Jacobian r^2. The second-chart
pair therefore satisfies the newly justified hypothesis at lambda r^2,
including r=0. Its fibre degree is still four. Polynomiality of the
leading coefficient r^16 V(1/r)^2 forces deg V<=8, so the chosen leading
square root r^8 V(1/r) is polynomial. Matching the three highest square
coefficients uniquely identifies its canonical root with the transported
source root. This establishes regularity on the entire chart, not just
at the generic boundary point.

Both resulting auxiliaries are global and have the declared degrees:
F=H^2+L, H in L_2 and L in L_1. If L were constant, the source identity
J(F,G)=2H J(H,G) would make a nonconstant H divide a nonzero constant.
Thus L is nonconstant. No step proves that H or L is an output-pencil
member, or that either inherits a Keller mate. The already proved L_3
exclusion cannot be applied to an auxiliary without that missing input.

The hostile F=t^4+2x^4t^2 is global but its root t^2+x^4 has exact
boundary pole order four. Its critical source line excludes a mate,
so it tests precisely the missing hypothesis in an automatic-gluing
claim. The separate [boundary-linear theorem](planar_jc48_sep08_boundary_linear.md)
provides the sharper F=b^4-v control: even a global root and an
affine-source critical-free coordinate with a rational mate do not pay
polynomial regularity. That result is a connection, not a dependency of
the present proof.

## 3. Frozen computational evidence

I independently ran normal and optimized Python, checked their equality
with the frozen output, and read the actual source. All 126 always-active
gates pass in each mode. The symbolic affine-root identity and literal
chart multiplier test the transport itself. The finite Faber bank with
independent coefficient derivatives checks both fluxes, while the
unbounded argument is the analytic dependency transfer above. The six
base-only examples explicitly allow the right side to vanish where
the leading coefficient vanishes. The fibre-dependent-Jacobian hostile
fails the theorem's actual hypothesis and has the asserted pole.

    python3 -B 04-computation/planar_jc48_sep08_quartic_boundary.py
    python3 -B -O 04-computation/planar_jc48_sep08_quartic_boundary.py

Source SHA256: `de2e1125bd6dc2446ac38b46b19e2f3e59ac71410b7965ff4b1716b49a86a439`.

Output SHA256: `e3c2a76c515933f92a13b4305a05561910491cf6c1c5a887a1c21c6e7ac8a939`.

The primary may be promoted to PROVED with this scope. No external
priority claim or full conjecture consequence is inferred.
