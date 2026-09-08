# Independent audit: all-finite quartic boundary type (4,3,1)

**Status: PROVED / INDEPENDENT ANALYTIC AND SOURCE AUDIT PASS.**
Root accepts the full all-finite theorem in [four_three_one](planar_jc48_sep08_four_three_one.md): for the declared complete global square-prefix class F=H²+L, a rational Jacobian mate forces L constant; no polynomial mate exists, with no bound on its degree. The proof concerns the fixed graph-complement surface and does not assert the general planar Jacobian conjecture. Fourfold infinity is outside this theorem.

I read the full proof, its 107-gate source, and the exact inherited local statements it uses. Both independent replays reproduce the frozen output. I also reconstructed the final primitive equations directly in the original (u,t) coordinates, independently of the producer's (u,q) computation. Both systems and the same-partition rational control agree exactly. A second agent independently audited the local divisor, generic integrality and Riemann--Roch mechanism before the source freeze.

## Complete entry, coordinates and first exactness rows

The complete global L2 and L1 coefficient spaces, not an arbitrary truncated polynomial family, are explicitly retained. Their second-chart restrictions are polynomial. If the t-leading coefficient of L vanishes then globality forces L constant. This is essential: arbitrary polynomials of degree zero in t are not all global on the surface.

The leading-field exactness theorem and the complete degree-eight classification impose the weighted midpoint relation on the three finite roots. The actual source scaling x=dX, t=d^-2 T normalizes their nonzero separation and has volume multiplier 1/d. Constant rescalings normalize the leading scalar. The remaining u=x-p translation is a coefficient change only; no automorphism of the graph-complement surface is inferred.

The inverse-series coefficient T2=-M/(4N) is correct. It follows after all preceding coefficients of the square prefix have been included. The corresponding coefficient of a rational mate has derivative T2/2. Field trace descends the rational differential M dx/N to C(x), paying its ordinary rational residues. Formal exactness is used only as a necessary condition here.

The first-jet supplier is genuinely a rational-mate obstruction: its finite tangent branches have nonzero logarithmic residues. It is not merely a source critical-point test. I re-read [shared_roots, Section 3](planar_jc48_sep08_shared_roots.md) and checked this direction explicitly. At a finite singular shared point, a mate therefore requires m>=3, normal order j>=2 and lower-section order n>=2. Normal-unit and simple-point cases remain separate and regular.

The fourfold point must be active. Without any active point there are no primitive poles on any compact component. If the triple alone is active, the exact rescaling w=tau², s=tau³Z gives four simple nonzero determinations paired into two normalized branches. Each differential pole has order two, so their total possible primitive pole degree is two. There are no other poles. The restriction F|D has degree four; a generic transverse D point gives a differential zero of order two and hence primitive local degree three. Comparing on its own component already contradicts the total upper bound two. No irreducibility is needed at this stage.

The fourfold active jets make M divisible by u². The simple-point residue requires M(-1)=0. The remaining residue at zero eliminates exactly its cubic free coefficient, yielding M=lambda u²(u²-1), lambda nonzero. Its displayed rational primitive differentiates correctly. In particular M(3)=72lambda, so the triple is now an M-unit point. The complete induced L row retains the term -2p lambda u; omitting it would change the moving residue.

## Moving residues, zeros and the anchored split

At the fourfold blow-up q=1+u²t, the actual volume is u^-2 du wedge dq. The normalized branches are original generic roots of f(q)=zeta. Differentiating their first displacement gives the displayed residue (g'f'-gf'')/(f')³. Since exactness must hold for every generic fibre and every one of those roots, its numerator is identically zero. Thus g=Cf' with a scalar C. Comparing the full polynomial coefficients gives C=2p, B=-108p and Q1=2pb.

The proof correctly treats p=0 before division. There A is free; it is not the specialization of the p!=0 coefficient formula. Both original source derivatives vanish on u=0 in that anchored family, which already rules out polynomial mates. The later rational exclusion retains the whole anchored family as well.

At the active fourfold point the exact leading quartic has nonzero constant, and a repeated nonzero root satisfies an equation independent of the fibre value. That equation is nonzero because its constant is -4a². Consequently its four generic roots are simple and nonzero; they exhaust the actual local degree. Each gives a differential pole of order two. Thus a hypothetical primitive has total pole degree at most four.

If P(3) is nonzero, the M-unit triple has normal order zero and its cancellation pair forms one branch with u-3=tau². The differential order is four, not two: the Puiseux factor and the derivative of the normalized parameter both contribute. Its primitive local degree would be five, greater than the entire possible pole degree four. This is a zero-of-the-form obstruction, stronger here than residue vanishing. It forces P(3)=0.

For p!=0 this yields b=-36p-54 and P'(3)=864p. For p=0 it yields A=30-b/9 and P'(3)=-6(b+54). The b=-54 case has normal order at least two and is not silently treated as balanced order one. Both remain in the proof.

## Geometric integrality and the complete differential divisor

The noncomposition proof is complete in the t-degree. A nontrivial outer polynomial must have degree two or four. In degree two, completing its square produces (K1-H)(K1+H)=L-constant. Choosing the top sign makes K1+H have degree two in t, while the right side has degree at most one. Therefore K1=H and L is constant, a contradiction. In outer degree four, the inner polynomial has t-degree one, so N² must be a fourth power up to scalar. The odd multiplicities of N make that impossible.

The nonsquare sidecar in the latter argument is necessary. The inherited example F=h⁴+h has nonconstant L and is composite. The proof records this failed shortcut and its repaired hypothesis rather than importing the false implication that nonconstant L alone proves noncomposition.

The classical [Arzhantsev--Petravchuk, Theorem 1 and Lemma 3](https://arxiv.org/pdf/math/0608157v2), independently read in the primary earlier in this research session, connect noncomposition, closed polynomials and relative algebraic closure of C(F). In characteristic zero this gives a regular function-field extension and geometrically integral generic fibre. The local canon route is also identified by both ID and path: [THM-3827, generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases](../../01-canon/theorems/THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases.md). No assertion about all special fibres is needed or made.

At the M-unit triple of normal order at least two, s has order two in u-3 and the leading cubic has three simple roots, giving a unit differential. In normal order one, a simple leading root also gives a unit. The cubic cannot have a triple root by its square-coefficient identity. At a double root its moving critical value has order at most two generically, because the fibre derivative has a nonzero quadratic leading term. Order one normalizes to a unit differential; order two gives forbidden nonzero logarithms. Thus every surviving triple branch has differential order zero under the hypothetical mate.

At the shared simple point use the full section numerator calN=N+sP+s²Q as local coordinate, not the one-variable boundary polynomial N. The two normalized branches have calN of order two in s and the relative form has order zero. At D there are four generic transverse points, each of differential order two. The point S intersect Dbar is absent from the generic curve since N8=1. At all other points of W generic smoothness and the nonvanishing source volume make the form a unit. Therefore its complete divisor has four double poles and four double zeros, with nothing omitted. Its canonical degree is zero on the geometrically connected compact normalization; the genus is one.

## Full primitive basis and descent

Let E be the sum of the four points above u=0. Any primitive has poles at most simple there and nowhere else, so it belongs to L(E). Its degree is four and genus is one, hence Riemann--Roch gives dimension exactly four. The displayed functions 1, 1/u, J and (u-3)H/u all lie in this space.

Their pole checks cover the whole curve. At the fourfold point q and H are finite and the denominator u allows only a simple pole. At the triple, t has order -2 in u-3 and H has order -1; the respective squared and single numerator factors cancel them, including on a ramified Morse branch. At the shared simple point the numerator of J vanishes at least to the necessary order and H is finite. At D, J and H are global and the displayed ratios have finite limits. The affine line u=0 has constant F and is absent from the generic fibre; no extra affine denominator remains.

Independence over the geometric generic constant field follows from t-degrees two, one, then zero, since a nonzero polynomial of degree below four cannot vanish modulo the geometrically integral generic equation. The same independence proves uniqueness of the expansion. Since every basis function is defined over C(zeta), the coefficients of a source rational mate descend to C(zeta). No algebraic constant-field enlargement is hidden in the final coefficient comparison.

Consequently the primitive is A(zeta)/u+B(zeta)J+C(zeta)(u-3)H/u+D(zeta). This is a complete pole-space consequence, not a finite mate-degree ansatz.

## Both coefficient contradictions independently recovered

I reconstructed the original full H,L in (u,t), used the ordinary source Jacobian directly, and reduced by F-zeta as a polynomial in t. I then imposed the two necessary leading coefficient relations from the primary. The cubic and quadratic remainders vanish exactly in both the nonanchored and anchored cases. The next original-coordinate remainder, divided by the actual u² transformation factor, gives precisely the primary's final obstruction.

For p!=0 its selected coefficient is

    48(h+192p²)zeta/lambda -(32p+48)(h+192p²)-lambda.

It is a nonzero rational function of the independent fibre parameter: if its slope vanishes, the remaining constant is -lambda. For p=0 the corresponding factor is

    108(b²+108c)zeta+2b lambda(b²+108c)-243lambda².

Again zero slope leaves a nonzero constant. Hence B=0 in either case, then A=C=0. A constant on the fibre cannot have Jacobian one. This proves the rational conclusion. When L is constant, the nonconstant polynomial factor 2H excludes polynomial mates. These are different final steps, and the proof retains their different scopes.

The independent original-coordinate calculation is reproduced below. It imports no producer implementation:

```python
import sympy as S
u,t,p,h,f,lam,b,c=S.symbols('u t p h f lam b c')
N=u**4*(u-3)**3*(u+1)
J=u*(u-3)**2*(u+1)*t+u*u-(2*p+5)*u

def jac(F,G):return S.cancel(S.diff(F,u)*S.diff(G,t)-S.diff(F,t)*S.diff(G,u))
def need(name,v):
 if S.cancel(v)!=0:raise RuntimeError(name)
 print(name,'PASS',flush=True)
bs=-36*p-54
cs=h-(2*p-1)*bs+108*p*p-108*p+27
P=u*u*(2*u**4+(-4*p-16)*u**3+(16*p-18-bs)*u*u-108*p*u+bs)
Q=u**4+(-4*p-8)*u**3+(4*p*p+16*p-36-bs)*u*u+2*p*bs*u+cs
H=N*t*t+P*t+Q
F=H*H+lam*(u*u*(u*u-1)*t+u*u-2*p*u)
G=3*(2*p+3)/u+J-4*(h+192*p*p)*(u-3)*H/(lam*u)
rem=S.Poly(S.cancel(S.rem(jac(F,G),F-f,t)),t)
need('Original t-coordinate cubic after elimination',rem.coeff_monomial(t**3))
need('Original t-coordinate quadratic after elimination',rem.coeff_monomial(t**2))
need('Original t-coordinate nonanchored obstruction',S.expand(S.cancel(rem.coeff_monomial(t)/u**2)).coeff(u,4)-(48*(h+192*p*p)*f/lam-(32*p+48)*(h+192*p*p)-lam))
P=u*u*(2*u**4-16*u**3+(30-b/9)*u*u+b)
Q=u**4-8*u**3+(12-b/9)*u*u+c
H=N*t*t+P*t+Q
F=H*H+lam*(u*u*(u*u-1)*t+u*u)
G=-b/(6*u)+J.subs(p,0)-(b*b+108*c)*(u-3)*H/(27*lam*u)
rem=S.Poly(S.cancel(S.rem(jac(F,G),F-f,t)),t)
need('Original t-coordinate anchored cubic after elimination',rem.coeff_monomial(t**3))
need('Original t-coordinate anchored quadratic after elimination',rem.coeff_monomial(t**2))
need('Original t-coordinate anchored linear obstruction',rem.coeff_monomial(t)/u**2-(u-3)**3*(u+1)*(108*(b*b+108*c)*f+2*b*lam*(b*b+108*c)-243*lam*lam)/(243*lam))
Z=u*u*(u-3)*t+u-2*p-3
H=(u-3)*(u+1)*Z*Z+c
GH=(u-2)/(12*u*(u-3)*Z)
need('Independent literal same-partition rational mate',jac(H,GH)-1)
print('Independent 431 original-coordinate audit PASS')
```

All seven independently derived identities pass. The same-partition rational family is not merely formal: H=(u-3)(u+1)[u²(u-3)t+u-2p-3]²+kappa is a global L2 section for every finite p, and its displayed rational mate has literal Jacobian one. Squaring H gives the constant-L rational control. This verifies that excluding rational mates also at constant L would be false.

## Source audit, pins and final scope

I read every source gate and checked both execution modes. The 107 gates include complete symbolic coefficient spaces and chart expressions, inverse coefficients, the actual residue equations, all parameter reductions, both full primitive systems and the rational controls. Local numeric order identities are controls of the written valuation proof; they do not replace analytic branch exhaustion. There is no search over a finite set of possible mate degrees.

Reproduce from the worktree root:

    python3 04-computation/planar_jc48_sep08_four_three_one.py
    python3 -O 04-computation/planar_jc48_sep08_four_three_one.py

Both reproduce the same 426 output bytes with 107 always-active gates.

| Artifact | Bytes | SHA256 |
|---|---:|---|
| Source | 14646 | `2d4051d677c1c1ad72d603d87a8da881cd8bf30a10c38aa3438952e3a0c3ef8d` |
| Output and both replays | 426 | `ca712f951b47665fb2314accfa8a231bfc9d03b99cb895823529d7c9f0bffc34` |
| Final primary before promotion | 17080 | `d94cf9a2cf9e6c874514e7275ed0681bad6525731940beeb1acb2467c0442b6f` |

Before freeze, root requested explicit full-section numerator notation and verified the theorem-ID path. These were clarifications, not changes to the mathematical conclusion. The separate anchored family and rational hostile were also fully incorporated before this final frozen audit. No frozen source or output was changed.

The accepted theorem covers every normalized finite p and every lower coefficient allowed by the full global L2/L1 spaces. A triple or simple point at infinity is already excluded by the odd-degree leading classification, but the fourfold-at-infinity case is not included in this primary's proof. This audit does not promote that remaining location or any wider quartic/Jacobian claim.
