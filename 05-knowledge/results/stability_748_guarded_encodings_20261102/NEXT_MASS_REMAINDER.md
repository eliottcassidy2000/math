# The next retained coordinate: variation inside the light groups

**Status: PROVED local inequality; its use to close any order below314 is OPEN.**
This keeps the exact squares discarded in the current proof. It is a
concrete continuation target rather than another independent worst-case bound.

Use the notation of STABILITY_314.md, but allow any positive t,d,lambda
with d=8t and lambda=t/8. There are L>=1 light groups with0<w_i<=t,
total selected light mass U and discarded light mass D_l. Let D_h=D-D_l.
The heavy matching has mass M_h and size at most m; all selected weights
are at most1/4. Define

    beta=max(0,M_h-(m-1)/4-t),
    r=U-t floor(U/t),                  0<=r<t.

Then the literal disjoint group wedge sum satisfies

    W >= tD-lambda(dL-4U)+(t/4)U+J,

    J >= (D_l+U-2tL)^2/(4L)+r(t-r)/4+beta D_h.         (1)

**Proof.** For a light group the exact quadratic identity leaves

    (S_i/2-t+w_i/2)^2+w_i(t-w_i)/4

after the tw_i/4 credit used by the314 proof. Cauchy--Schwarz sums the
squares to the first term of(1). With fixed total U and0<=w_i<=t,
the concave sum sum w_i(t-w_i) is minimized at an extreme point of that
box slice: all but at most one weight are0 or t. The remaining weight
is r, giving r(t-r). Permitting zero weights only weakens this lower bound.
Every heavy selected weight is at least M_h-(m-1)/4 as the others are
at most1/4. Its extra selected-edge cost above tS_i is therefore at
least beta S_i. The groups' counted wedge pairs remain disjoint. QED.

Even D_l may be eliminated WITHOUT forgetting its compatibility with U.
The minimum of the right side over0<=D_l<=D is obtained at

    D_l=clamp(2tL-U+2beta L,0,D).

This follows by differentiating its convex quadratic. Thus(1) has a
compact computable lower bound in the actual aggregate data(D,U,L,M_h,m).
It loses individual group orientations and correlations, so feasibility
of its scalar constraints would not be an actual tournament witness.

Keeping the baseline-BR as in the314 proof, a sufficient replacement
for its tail estimate is

    J >= lambda [4V(3/4)^m-(c-d)L-2U]_+.               (2)

With adaptive t=4/R,lambda=1/(2R), the final baseline contradiction is
automatic once(2) holds: kappa_n<4/R and lambda R=1/2. At R=156,
the present slot count only forces m>=8 and the old scalar tail test
can fail there. Equation(1) is the next precise test: a small heavy
matching near its mass capacity forces heavier owners, while the light
mass and discarded mass cannot both make their quadratic remainder zero.

This does not prove(2). Next work should optimize or give an exact dual
certificate over the compatible(V,B,D,U,L,M_h) domain, with explicit
hostiles. If a feasible abstract point survives, inspect whether actual
edge assignments and directed triangle incidence can realize it before
treating it as an obstruction to the mathematical theorem.
