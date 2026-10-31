# Session 128: audit of the attached Amanda manuscript

## Source and scope

Audited all 23 pages of Ryan E. S. Charmley and Steve Konrad, *The Konrad backbone in an Amanda loop*, version 1.0, December 22, 2025, DOI [10.6084/m9.figshare.30934796](https://doi.org/10.6084/m9.figshare.30934796). Attachment filename: `2025-Konrads_Backbone_Amanda_Loop-Charmley-v1.0.pdf`; 1,213,371 bytes; SHA-256 `de4c94232203e09156ccd233ff0a374540ad1daeca0e1e31ac3f88dea21193a3`.

Text was extracted and actual rendered pixels were inspected for pp. 7, 15, 17, 18 and 23, resolving lost absolute-value bars and the genuine factor-two discrepancy. Page 23 is artwork, not a blank page or an additional mathematical theorem. The PDF is an input, not a republication deliverable. This note contains an original audit and derived mathematics, not an endorsement of the manuscript's ontology. The manuscript itself identifies its repository publication and lack of assured external peer review on p. 22.

## Immediate decisions

1. Retain the exact Fibonacci identities, cubic polynomial/trigonometric identities, Fourier inverse, and positively seeded Tribonacci ratio theorem with explicit hypotheses.
2. Repair Eq. (60) to match Eqs. (69) and the projection identity in Proposition 3.2.
3. Do **not** import the claimed automatic long-run symmetry: the native loop has a rigorously certified persistent nonzero parity example below.
4. Treat categorical equivalence, Pisot substitution attractors, and exact linear feedback absorption as **new constructions with new proofs**, not results supplied by the attachment. The native Step F is a nonlinear threshold loop.

## Page-by-page mathematical audit

### Pages 1-5: algebraic backbone versus conserved invariant

The ontological symbols in Eq. (1) are stipulated metaphors/axioms, not division operations in a defined field. No physical or foundational conclusion follows from the subsequently correct elementary identities.

On p. 3, Eqs. (3)-(6) correctly reduce

`K(a,b,c) = (a^3-a*b^2-c^3)/(b*c)`

to `2b+3c` on `a=b+c`, provided `bc != 0`. The Fibonacci specialization needs `n>=2`, since at `n=1` the stated quotient has `c=F_0=0`. One can separately extend the linear expression to `n=1`; this is not an evaluation of the original quotient there.

The word “invariant” is not a conservation law under the shift `T(a,b,c)=(a+b,a,b)`. For example, `K(2,1,1)=5`, whereas `K(3,2,1)=7`. What is invariant is the additive constraint plane. The sequence of functional values is Fibonacci-recursive, and the ratio limit on the standard Fibonacci orbit is correct (p. 4, Eqs. (7)-(8)). Assertions that arbitrary recursive feedback is optimal or converges to the golden ratio are not established.

The golden-angle resonance numbers on p. 5 are numerically consistent: `E=1.649063135827318...`, relative discrepancy from `sqrt(e)` is `0.000207351681127...`. This is a numerical proximity, not an identity, dynamical invariant, or explanation of the later dynamics. The classical Fourier identity in Eq. (11) is correct.

### Page 6: triadic sine closure

Equations (12)-(14) and Proposition 2.11 are correct. With `q(theta)=sin(theta)` and `m(theta)=3theta mod 2pi`, one has the exact semiconjugacy

`q o m = T3 o q`, where `T3(s)=3s-4s^3`.

The zero grids of `sin(3^k theta)` are nested and their union is dense. This does not assert density of every forward orbit, contraction, or convergence to a preferred attractor. Counterexamples are the fixed point `s=0` and the two-cycle `1 -> -1 -> 1`. The map `q` is not injective, so it does not itself define an equivalence/conjugacy. It becomes a conjugacy only after the precisely defined quotient `theta ~ pi-theta (mod 2pi)`; tripling respects this equivalence relation. The sine polynomial is the **negative** of the cosine Chebyshev polynomial `C3(s)=4s^3-3s`; their shared degree does not identify the dynamics.

### Pages 7-8: cubic and quadratic seals

Equations (15)-(23), their stated Fibonacci domains, Cassini signs, the penta-seal, and square-lift recurrence are correct. Rendered p. 7 confirms that Corollary 2.18 includes the absolute value needed for its final value `3`; that bar was lost in extracted text.

A cubic-seal audit is not a complete converse to additive closure over arbitrary real triples. Factorization gives

`a^3-b^3-c^3-3abc = (a-b-c) * [(a+(b+c)/2)^2 + 3(b-c)^2/4]`.

Thus `(-1,1,1)` passes this seal while violating `a=b+c`. Positivity of all three coordinates repairs the converse. The square-scale recurrence `L_{t+1}=3L_t-L_{t-1}` actually holds for every exact additive/Fibonacci-type backbone, not only the standard integer Fibonacci seed. Conversely it does not certify that the seed was the standard Fibonacci seed.

### Pages 8-10: analytic closures

The Dirichlet sawtooth Eq. (24) is correct on `0<x<2pi`; at multiples of `2pi` the sine series equals `0`. Equation (26) is a legitimate explicitly chosen `n=0` extension at `x=1`, not the ordinary value of an undefined quotient without that convention.

The Borwein plateau and first-deviation fraction are correct. The exact fraction in Eq. (29) was independently checked by a rational simplex-tail calculation in `verify_persistent_bias.py`. If independent `U_k` are uniform on `[-a_k,a_k]` and `a_0 >= sum_{k>=1} a_k`, then convolution at zero is `1/(2a_0)`, giving the stated integral. For the seven odd reciprocals, put `delta=sum_{k=1}^7 1/(2k+1)-1`. Since `0<delta<2/15`, the two tails are simplices and the normalized integral is exactly

`1 - delta^7 / [2^6 * 7! * product_{k=1}^7 (1/(2k+1))]`.

It equals twice the rational coefficient of `pi` printed in Eq. (29). Merely failing a sufficient condition would not prove deviation; this tail calculation does. The arctangent integrals and alternating series are correct; the endpoint termwise-integration claim should be justified by a limiting argument rather than uniform convergence of the geometric series on `[0,1]`.

### Pages 10-13: product normalization, Viete, and Tribonacci

The manuscript correctly rejects the malformed unnormalized finite-limit fragment in Eq. (32). The normalized product in Eqs. (35)-(37) is correct. Binet's formula yields

`R_k/F_{k+1} -> phi^(-1) * product_{n>=1} (1+phi^(2-4n))/(1-phi^(-4n))`.

Absolute summability of the corrections proves convergence. Its reported value `1.09294258499686746022...` agrees with independent evaluation. The accumulator update prose on pp. 11 and 19 is ambiguous/off by one if “odd steps” refers to the old index `n`: precisely `Q_{n+1}=Q_n*F_{n+1}^{(-1)^(n+2)}`, hence multiply when the **new index** `n+1` is odd, divide when it is even.

Viete's trig formula in Eq. (39), under `4p^3+27q^2<0`, is correct. The broad “real radicals cannot express” prose needs an irreducibility/base-field hypothesis. For example, `t^3-t=0` lies in the three-distinct-real-root regime and all its roots are rational. The discriminant condition alone is not casus irreducibilis over a stated coefficient field.

The Tribonacci matrix, characteristic polynomial, constant, and radical formula in Eqs. (40)-(47) are correct. A precise repair to Theorem 2.43 is: for a real initial state `x_0=(u_2,u_1,u_0)`, require nonzero Perron projection

`u_2+(alpha-1)u_1+alpha^(-1)u_0 != 0`.

Then denominators are eventually nonzero and `u_{n+1}/u_n -> alpha`. Positive seeds suffice. In fact every nonzero rational seed suffices by irreducibility of `x^3-x^2-x-1` and `alpha^(-1)=alpha^2-alpha-1`. “Nondegenerate” must mean this spectral condition, not merely nonzero real seeds. For instance `u_n=Re(beta^n)` for a nonreal root `beta` is a real nonzero-seed solution with no Perron component and cannot have ratio limit `alpha>1`.

The closure depth in Eqs. (50)-(51) is a first-passage statistic; it is not a uniform iteration bound or a permanent-error certificate. A positive example can hit zero residual once and immediately leave: let `u_2=1` and `u_0=u_1=(alpha-1)/2`. Then `u_3/u_2=alpha`, but `u_4/u_3=3/2+1/(2alpha) != alpha`. No substitution, star embedding, Rauzy tile, Pisot-attractor, or pure-discrete-spectrum theorem appears here.

### Pages 14-20: exact native Amanda loop

The additive backbone shift preserves its plane. The positive-seed Tribonacci adapter has the stated limiting ratio. The convex gain blend is well-defined once its denominators are nonzero; calling its scalar a “contraction” does not establish contraction of the whole coupled state evolution.

Equation (57) specifies amplitudes `A_j=-kappa*rho_t^* w_{j-1}p_j^2 sigma_j`. The firing schedule is an arbitrary finite subset; the paper gives no law selecting it. Irrational rotation equidistribution concerns the full sequence (or specified prefixes), not every arbitrarily selected sparse subsequence.

**Actual error in Eq. (60), visually verified on p. 15:** its exponential coefficients `1/2,-1/4` give `cos(3theta)-(1/2)cos(6theta)`, twice the printed cosine expression. To preserve Eq. (69) and Proposition 3.2, choose

`K(theta)=(1/2)cos(3theta)-(1/4)cos(6theta)`

and change the exponential coefficients to `1/4,-1/8`. Doubling the cosine and subsequent field instead is another consistent normalization, but one must not mix them. The persistent-parity counterexample is unaffected by this positive factor of two.

With normalized convolution `dtheta/(2pi)`, the source-correct Poisson operator is `Delta_r=r^(-2)d^2/dtheta^2`. On mean-zero fields its inverse has Fourier multiplier `-r^2/n^2` for `n!=0`. Equations (63)-(69) are correct with the repaired normalization. In particular,

`phi(theta)=sum_j A_j[-r^2/18 cos(3(theta-theta_j))+r^2/144 cos(6(theta-theta_j))]`.

The carrier `E=span_R{cos3theta,sin3theta,cos6theta,sin6theta}` is exact, not an approximation. Step C's neutralization is mathematically redundant for this kernel but useful as a numerical audit.

**Basel/gauge correction:** the piecewise quadratic `G` in Eq. (67) has zero mean, despite its nonzero constant term in that pointwise formula. For every source `f`, `Delta_r(G*f)=f-mean(f)`, and `G*1=0`. Thus the statement after Eq. (67) that the Basel constant becomes “audible” when neutrality is perturbed is not true of this mean-zero inverse. A nonneutral right-hand side has no periodic Poisson solution. Adding a kernel constant is a gauge change, not a repair of this obstruction.

Audits are measurements, not automatic enforcement mechanisms. J1 vanishes identically in exact arithmetic for the prescribed source. Rendered p. 17 confirms J3 is `|1-(W+tail)|`; its advertised tail equals the actual tail for the natural choice `w_k=phi^(-k-2)`. Merely normalizing arbitrary weights does not supply that specific tail. Requiring J3 to vanish for every cap forces precisely these geometric weights. J4/J4b are constant numerical discrepancies, not changing feedback diagnostics. J5's projection identity is correct with the repaired kernel but its normalization is undefined if all event amplitudes vanish (including an empty firing set). J6-J9 are exact identities under their premises; they do not prove the premises or global dynamics. J10 and J11 are asymptotic diagnostics with the qualifications above.

The paper calls J1-J3 required structural invariances, but defines no reject/retry/stop rule, acceptance tolerance, or correction of the firing schedule. J2 only measures parity, and J5 is explicitly optional. The counterexample below concerns the explicit update equations, not a hypothetically added hard acceptance filter. Such a filter could reject the witness rather than correcting it; the manuscript supplies no proof that accepted states are eventually reached.

The source's actual Step F (p. 20, Eqs. (82)-(83)) is

`b_{t+1}=(1-lambda)b_t+lambda*Pi_t`, `b_0=0`, `0<lambda<=1`,

`sigma_j(t+1)=sgn_+(sin(theta_j)-b_{t+1})`,

where `Pi_t=(1/(2pi))*integral sgn(phi_t')`. The optional seed is `sin(3theta_j)`. The memory filter is linear in its two inputs, but the coupled field-to-sign-to-field loop is nonlinear and discontinuous. It cannot silently be replaced by a linear resolvent equation.

### Pages 21-23: interpretation boundary

The discussion repeats the true component identities, but its claim of enforced long-run symmetry is false without further hypotheses or a different controller, as the certified example below proves. “Context lattice” is an interpretation; no partial order/lattice operations or categorical equivalence functors are defined. Page 22 supplies publication/attribution metadata. Page 23 contains artwork and no further mathematical result.

## Repaired native theorem and certified counterexample

The following are new deductions from the exact native equations, not claims proved by the attachment.

### A. Exact carrier, bounded parity, scalar-gain invariance

For any finite firing set and repaired kernel, the unique zero-mean Poisson solution lies in `E`. Write `g=phi'`. Its modes are `+-3,+-6`, so

`g(theta)+g(theta+2pi/9)+g(theta+4pi/9)=0`.

Unless `g` is identically zero, it has finitely many zeros. Almost every triple therefore includes at least one positive and one negative value; its three signs sum to `+1` or `-1`. Integrating gives `|Pi|<=1/3`. If `g=0`, use the usual `sgn(0)=0`, so `Pi=0`. Convexity then gives `b_t in [-1/3,1/3]` for all `t`. This parity argument was independently supplied by the exact-feedback collaborator.

For a fixed schedule, geometry, tier weights, and sign state, `phi` and `g` scale by the common positive factor `kappa*r^2*rho_t^*`. Hence their sign statistic is exactly independent of these positive scalar gains. Blending the Fibonacci and Tribonacci gain alone cannot tune native parity. A negative or zero gain, different geometry, changing event schedule, or changing relative weights is outside this assertion.

### B. Two-harmonic zero-count lemma

Let

`g(x)=R sin(x-delta)+S sin(2x-epsilon)`, with `R>0`, `S>=0`, `5S^2<R^2`.

Every zero lies in the two arcs where `|sin(x-delta)|<=S/R`. On the arc around `0`, the derivative is strictly positive, since its fundamental contribution is at least `sqrt(R^2-S^2)>2S`; on the arc around `pi`, it is strictly negative. Values at the endpoints have the signs required by the dominant fundamental, so continuity gives one zero in each arc. There are exactly two simple zeros per period. The argument includes `S=0` directly; strict inequality excludes degeneracy.

### C. Native counterexample valid for every allowed memory parameter

Start the standard Fibonacci backbone at `(a_0,b_0,c_0)=(3,2,1)` (here the middle backbone coordinate is distinct from the scalar axis memory); then its gain ratio is strictly between zero and one at every step. Choose the unchanging firing set `S_t={1,2}` and canonical normalized weights `w_k=phi^(-k-2)`. Keep the paper's golden-angle locations and use initial signs `(sigma_1,sigma_2)=(+1,-1)`, consistent with its zero-memory seed. The certified seed bounds are

- `h_1=sin(theta_G) > 0.67549 > 1/3`;
- `h_2=sin(2theta_G) < -0.99617 < -1/3`.

Part A forces `|b_t|<=1/3`, so neither sign can ever switch. Dropping common positive scalar factors and using `x=3theta`, the derivative has the same sign as

`g(x)=sum_{j=1}^2 A_j[4sin(x-3jtheta_G)-sin(2x-6jtheta_G)]`,

where `A_1=-4phi^(-2)` and `A_2=9phi^(-3)`. If its fundamental coefficients are `u,v` and second-harmonic coefficients are `w,z`, then the certificate establishes

`R^2=u^2+v^2 in [46.3714429355552027073977,46.3714429355552027073978]`,

`5S^2=5(w^2+z^2) in [42.6687811555978962251999,42.6687811555978962252000]`.

Part B proves there are exactly two zeros. Rigorous endpoint signs locate them at

`r_1 in [3.0450368193027045,3.0450368193027046]`,

`r_2 in [5.5300956371054050,5.5300956371054052]`.

Since `g(0)<0`, the positive set is exactly the interval between those two roots, modulo the period. Therefore

`Pi=(r_2-r_1)/pi-1`

lies in

`[-0.2089812105451970362491,-0.2089812105451969407560]`.

(The substitution `x=3theta` repeats the same period three times and leaves this mean unchanged.) In particular parity never vanishes. For **every** `0<lambda<=1`,

`b_t=Pi*(1-(1-lambda)^t)`.

This is a global, forever-valid counterexample to convergence under the explicit update equations, not a finite simulation or a claim that all native trajectories fail. If parity is instead imposed as an external hard acceptance criterion, the example is rejected, not repaired; no such acceptance procedure is specified in Steps A-F. It works with every positive scalar gain trajectory allowed by the positive Fibonacci/Tribonacci adapter, and with either consistent normalization of the factor-two kernel.

## Verification and transfer notes

`verify_persistent_bias.py` uses only the Python standard library. All pass/fail gates use exact rational interval arithmetic: outward rounding at 450 decimal places; Machin's formula with alternating rational arctangent bounds; an integer square-root bracket; degree 210/211 Taylor polynomials with Lagrange remainders; and the analytic zero-count lemma. It does not infer root completeness from samples. Normal and `python -O` outputs agree byte-for-byte. The exact-feedback collaborator independently reviewed the interval logic, coefficient signs, and analytic zero-count argument; an explicit in-period root-order gate was added on that review. The rational Borwein check is separate and exact.

Keep the own-work checker, its JSON certificate, and this audit eligible for incorporation into Session 128. Do not include extracted source text, rendered source pages, helper scripts, or the original attachment in a repository publication. No repository modifications or external mutations were made by this audit worker.
