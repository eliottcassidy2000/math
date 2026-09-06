# Polynomial charts, Student response, and delayed modular depth recognition

**Status: PROVED + FINITE-EXACT; [independent audit accepted](continuing8_20260906_connection_boundaries_audit.md).** This is a bounded incoming-work synthesis. It proves a natural scalar-intertwiner obstruction and an all-prime-power family with an exact positive-characteristic recognition clock. It supplies no characteristic-zero Keller obstruction, no new valuation-11/source-memory classification, and no general Mahler-orbit consequence.

## 1. Inheritance and the connection board

The incoming baseline is `45995e6d7`. The complete current board and the two relevant audited proofs were read:

- `05-knowledge/results/planar_jc_long_20260906_board.md`;
- `05-knowledge/results/planar_jc_long_20260906_smatrix.md`, especially the actual unrestricted Student operator, integral saturation and the separate generic mixed-minor tester;
- `05-knowledge/results/planar_jc_long_20260906_depth.md`, especially the full integral diagonal image, fixed-support recognition cutoff and characteristic-three hostile.

The inherited source normalization is **THM-4308 / `01-canon/theorems/THM-4308-source-normal-bracket-hasse-truncation-through-row-eight.md`**. The diagonal basis is inherited from **THM-4369 / `01-canon/theorems/THM-4369-source-packet-pascal-circuit-kernel-and-boundary-basis.md`**, completed and used in the incoming depth report. The source-memory correction remains **THM-4426 / `01-canon/theorems/THM-4426-source-normal-row-fourteen-weight-eighteen-memory-repair.md`**: a frozen or restricted response is not the complete source operator.

Our audited `continuing8_20260906_coin_rational_handoff_audit.md` proves that nonconstant polynomial substitution over a residue field preserves a nonzero denominator polynomial. Our audited `continuing8_20260906_mahler_reset_cost.md` instead concerns ordinary integer launches, actual binary termination and a fixed reduced denominator `3^j`. Neither statement identifies a JC source map automatically.

| Live concept | Exact operation and predicate | Missing coordinate or decisive control |
|---|---|---|
| Coin denominator content | Compose rational integral germs by charts; preserve nonzero reduced denominators | A chart becoming constant modulo the denominator prime; the nonunit affine hostile |
| Actual Student response | Apply `(x^2+6) derivative -2mx` in its declared integral lattice | Scalar chart intertwining, quadratic parameter, torsion; test the source functions 1 and x |
| Complete depth image | Signed diagonal extraction and division by `(1-z)^rho` | Fixed support and prescribed zero tail; actual-source Frobenius family below |
| Modular nonvanishing versus observation depth | Substitute `z -> z^(p^a)` | Vanishing order multiplies although the polynomial remains nonzero |
| Mahler native termination | Convert fixed carry prefix to actual integer height | Positivity, ordinary address, and reduced phase denominator; these have no supplied JC map |

No new literature claim is made about the S-matrix paper. Its results are not dependencies of either theorem below.

## 2. Exact obstruction to the natural Student chart transfer

Let K have characteristic zero, let `m` be nonzero in K, and put

`L_m f=(x^2+6) f' -2mx f`.

**Theorem 1.** Suppose a nonconstant polynomial `phi in K[x]` and a rational function `h in K(x)` satisfy

`L_m(theta(phi(x)))=h(x) (L_m theta)(phi(x))`                 (1)

for every polynomial theta. Then either `phi=x,h=1` or `phi=-x,h=-1`. Both possibilities do satisfy (1).

This is an identity for the unrestricted differential expression. For the actual Student indices `m>=2`, the two test functions used in the proof already belong to the allowed source space `deg theta<=m-1`. A nonlinear chart can also violate that degree cap; the theorem finds a failure even before imposing the cap. No claim excludes a more general matrix/gauge transport or a changed operator.

**Proof.** Apply (1) to theta=1. Since `2m!=0`, it forces `h=x/phi`. Applying it to theta=x and cancelling the identical zeroth-order terms yields

`(x^2+6) phi phi'=x(phi^2+6)`.                            (2)

If phi has degree d and leading coefficient a, the leading coefficients in (2) are `d a^2` and `a^2`, both at degree `2d+1`. Characteristic zero forces d=1. Write `phi=ax+b`, with a nonzero. The left side minus the right side in (2) is

`-ab x^2 +(6a^2-b^2-6)x+6ab`.

It vanishes only if b=0 and `a^2=1`. Conversely direct substitution proves the two permitted identities. QED.

The assumption on m is sharp. At `m=0`, every nonconstant polynomial phi satisfies (1) with

`h=(x^2+6)phi'/(phi^2+6)`.

The denominator is a nonzero polynomial. Thus a proof that cancels m without recording it would be false.

For the coin-safe chart `phi=x^2`, the scalar forced at nonzero m is `h=1/x`. The remaining theta=x residual is exactly

`x^3+12x-6/x`,

which is nonzero. Preservation of the coin's denominator predicate therefore does not preserve the Student response equation under this natural map.

There is an exact surviving covariance if the quadratic parameter is retained. With `L_(m,c)=(x^2+c) derivative -2mx`, every nonzero a satisfies

`L_(m,c)(theta(ax))=(1/a)(L_(m,a^2 c)theta)(ax)`.

Changing the source chart changes the parameter and potentially the integral coefficient lattice. Setting c=6 and suppressing that change recovers the obstruction of Theorem 1.

## 3. An exact delayed-recognition family in every prime characteristic

Use a separate symbol pi for the source generator, reserving p for a prime. The variables x and t are algebraically independent indeterminates. Over the integers define

`u=x^2 t`, `pi=t(1+u)`, `y=xt pi`.

At depth zero the actual monomial source module is the algebra `Z[pi,y]`. Over a field k use its corresponding monomial span `P_0(k)=k[pi,y]`. Projection `pr_N` retains all t-rows through N. The polynomial being tested is fixed, including all its zero rows beyond its stated support.

Let `p` be any prime, `a>=1`, and write

`P=p^a`, `L=floor(P/2)`, `epsilon=P-2L in {0,1}`, `T=3L+epsilon`.

**Theorem 2.** Over any field of characteristic p, the fixed polynomial `t^T` first fails projected membership in `P_0` at exactly

`N_p=T+P`.

Thus it passes through `N_p-1` and fails at `N_p`, although it never belongs to `P_0`. The corresponding characteristic-zero universal sharp cutoff for support T and depth zero is

`N_0=floor(4T/3)+1=T+L+1`.

The delay is

`N_p-N_0=ceil(P/2)-1`,

which is unbounded as a grows for every fixed prime. At P=2 the delay is zero; the theorem does not claim strict delay in that smallest case.

**Actual lift.** The exact integer source identity

`pi^3-y^2=t^3(1+u)^2`

gives

`F=pi^epsilon (pi^3-y^2)^L=t^T(1+x^2t)^P`.

In characteristic p, Frobenius reduces this to

`F=t^T+x^(2P)t^(T+P)`.                                  (3)

Consequently an actual source polynomial, rather than an ambient formal compensator, matches `t^T` at every row before `T+P`.

**No alternative source can postpone failure.** The signed diagonal of intercept `ell=2T` uses the coefficients of `x^(2n-2T)t^n`, multiplied by `(-1)^(n-T)`, and records them at `z^(n-T)`. The target `t^T` has diagonal polynomial 1. Every source monomial lies on one full diagonal, so other intercepts cannot affect this test.

All depth-zero source monomials on this intercept are precisely

`pi^(T-3j)y^(2j)`, `0<=j<=L`.

Their traces are

`(-1)^j z^j(1-z)^(T-j)`

`=(1-z)^P [(-1)^j z^j(1-z)^(L-j)]`.

The bracketed polynomials form an integral triangular basis of all polynomials of degree at most L, with diagonal entries `(-1)^j`. Reduction modulo p keeps this basis invertible. Hence the **full** diagonal source image, not merely the trace of (3), is exactly

`(1-z)^P k[z]_(<=L)`.                                   (4)

In characteristic p, `(1-z)^P=1-z^P` and

`(1-z)^(-P)=1+z^P+z^(2P)+...`.

Matching the target 1 through relative order R means that a polynomial q of degree at most L agrees with this inverse through R. Because `L<P`, q=1 works through `P-1`, but agreement through P would require `q_P=1`, impossible. This proves the exact first failure at absolute row `T+P` and excludes every alternative truncated source. Nonmembership at all depths also follows directly because 1 is not divisible by the nonconstant factor `(1-z)^P`. QED.

This strengthens the incoming isolated characteristic-three control to an explicit all-prime-power family. It is a theorem about the declared reduction of the actual source module. It does not change the incoming characteristic-zero recognition theorem, bound the support of a Keller pair, or assert that a modular source lift lifts to a characteristic-zero bracket solution.

## 4. The usable transfer and the saturation boundary

The same substitution `z -> z^P` is injective on a polynomial ring over a field and multiplies the order of every zero at the origin by P. The coin argument uses only the first property. A finite depth observer needs the second coordinate too. This is the first failed implication in a proposed transfer from “nonzero modulo p” to “detected at the characteristic-zero observation row.” Equation (3) is the exact decisive control.

There is a sound one-sided modular consumer for the incoming full depth modules. For an integral target prefix, if its reduction modulo p fails the complete projected depth test, then the original prefix fails over Q. The reason is not generic modular rank: on each diagonal the full source basis has a leading unit triangular block. Every finite projection therefore has a saturated integral image (either the whole ambient diagonal or the graph of an integral map over the first independent coordinates). Rational membership of an integral prefix consequently supplies an integral lift, whose reduction would be a modular lift. This is a corollary of the incoming full integral image. Modular acceptance, by contrast, is not a characteristic-zero certificate, and Theorem 2 makes the delay arbitrarily large.

That reasoning cannot be transferred to the Student bank without its torsion sidecar. For the actual index-two integer operator

`A_2=((0,6),(-4,0),(0,-3))`,

the integral target `2x=(0,2,0)^T` has the unique rational source `theta=-1/2`. It passes the equation modulo 2 with the zero source, but no source solves it modulo 4: the middle response coordinate is always divisible by 4. Thus a modulo-4 rejection here coexists with characteristic-zero compatibility. The incoming Smith group records exactly this missing integral membership coordinate. A primitive denominator, a saturated image and a nonzero integer determinant are different predicates.

The Mahler actual-termination theorem supplies a related methodological restriction, not a ring map: its reset lag bound needs a fixed reduced denominator `3^j`, the ordinary address and its actual binary height. The Frobenius exponent P here is a characteristic-dependent source multiplier, with no supplied relation to that reduced phase denominator or to positive ordinary orbit states. Neither its growing delay nor the Mahler logarithmic reset bound therefore transfers to the other system. A claimed future bridge must exhibit a source map carrying those coordinates explicitly.

The strongest next object is a finite measurement packet retaining the actual operator, the declared support and zero tail, and its integral saturation/valuation data. For depth modules the unit basis supplies saturation and allows sound modular rejection. For the unrestricted Student response the incoming Smith factors supply the arithmetic sidecar. For coin charts one must retain both reduction degree and row capacity. These are concrete separate consumers, not a new scalar invariant shared by all three problems.

## 5. Reproducible exact controls

The adjacent `continuing8_20260906_connection_boundaries.py` imports no research producer. It verifies the complete linear-chart coefficient equation, both surviving charts on monomials through degree eight, five nonlinear/affine hostiles, the sharp m=0 exception, and the parameter-covariant scaling identity. It verifies the actual source identity over Z before reduction.

Its declared Frobenius universe is p in `{2,3,5,7,11,13}`, a in `{1,2,3}`, with `p^a<=125` (14 families). Every family is checked against the full triangular source basis, every Frobenius coefficient, the inverse-series first obstruction and both cutoffs; all families with `P<=13` additionally use literal bivariate expansion of the actual source expression. The script also checks the actual Student `2x` target modulo 2 and all source residues modulo 4. No field-rank fit, prime-bank inference or generic source compatibility is used.

Run with Python and with Python `-O`; both modes pass **210 always-active exact gates** and produce identical LF output. The universal claims follow from the proofs above, not extrapolation from the bounded controls. The parent owns filing and status promotion; no Git or maintained document was edited by this scout.
