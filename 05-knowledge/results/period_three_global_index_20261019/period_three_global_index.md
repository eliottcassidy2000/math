# The period-three residue from a global fixed-point budget

19 October 2026. **PROVED, scoped:** the residue-budget application, its resultant/discriminant explanation, and the two-cycle coalescence calculation below. **INHERITED:** the period-three cubic, discriminant 49, cyclotomic realization, and the already established iterative residue 51/49. **FINITE-EXACT:** the companion controls. This note gives a new proof and structural explanation within this research sequence, not a historical-priority claim.

## 1. Inheritance and the question

The [quadratic return note](../quadratic_clock_return_20261012a.md), Section 2, establishes the parabolic three-cycle of `f_c(x)=x²+c` at `c*=−7/4`, with

\[
 P(x)=x^3+\tfrac12x^2-\tfrac94x-\tfrac18,
 \qquad f_{c*}^{\circ3}(x)-x=(f_{c*}(x)-x)P(x)^2. \tag{1}
\]

Its three distinct real roots form one primitive cycle, and `disc(P)=49`. The [seventh-root bridge](../arithmetic_braids_20260917_geometry.md), Section 2, already identifies this cycle with the reversed doubling action on seventh-root traces. The [difference-visibility note](../difference_visibility_norms_20261004.md), Section 4, already records the real cubic cyclotomic field, its discriminant, and angle tripling. These arithmetic facts are not new here.

The [cotangent/parabolic companion](../parabolic_cotangent_integrability_20261019/parabolic_cotangent_invariant_barrier.md), Section 5, computes the iterative residue `51/49` by reducing the quadratic and cubic return coefficients modulo `P`. The present question is why that local number is forced by the other fixed points, and whether its denominator has a proved connection to the existing discriminant.

The closest mechanism is the classical holomorphic fixed-point index formula: see [Milnor, *Dynamics in One Complex Variable*, 1991 revised preprint, Section 9, Theorem 9.2 and Problem 9-1](https://abel.math.harvard.edu/archive/118r_spring_05/docs/milnor.pdf#page=63). We give the needed polynomial argument directly. The concept board is: local index / global residue / cycle weight / collision resultant / cyclotomic action. The hostile is confusing algebraic multiplicity with the number of distinct orbit points, or treating a shared denominator as an arithmetic theorem without a map connecting the quantities.

## 2. A precise polynomial index budget

For a polynomial `F` of degree at least two and a finite fixed point `p`, define

\[
 I_F(p)=\operatorname{Res}_{z=p}\frac{dz}{z-F(z)}.
\]

**Lemma 1 (classical polynomial index formula).** Summing over distinct finite fixed points, including multiple ones,

\[
 \sum_{F(p)=p} I_F(p)=0. \tag{2}
\]

If `F′(p)≠1`, then `I_F(p)=1/(1−F′(p))`. If
`F(p+u)−p=u+a u²+b u³+O(u⁴)`, with `a≠0`, then

\[
 I_F(p)=b/a^2. \tag{3}
\]

**Proof.** The rational one-form `dz/(z−F(z))` has no finite poles except fixed points. Its coefficient is `O(z^(−deg F))` at infinity, so its residue at infinity is zero. The residue theorem proves (2). The two local formulas follow by Laurent expansion. In particular, in (3) the expansion begins `−1/(a u²)+(b/a²)/u`. ∎

The sum is zero for finite points of a polynomial. The holomorphic fixed-point index at the point infinity on the sphere is one; that is a different quantity from the residue at infinity of the affine one-form just used.

The index is invariant under a local holomorphic coordinate change `h`. Indeed, the difference

\[
 \frac{h'(z)}{h(z)-h(F(z))}-\frac1{z-F(z)}
\]

is holomorphic near a fixed point: Taylor expansion in `F(z)−z` cancels its apparent pole. Thus its residue is zero.

**Corollary 2 (cycle-weighted budget).** Let `f` be a polynomial of degree at least two, and put `F=f^{∘n}`. Partition the finite fixed points of `F` into primitive `f`-cycles `C` of lengths `ℓ_C` dividing `n`, with primitive multipliers `μ_C`. Then

\[
 \sum_{\mu_C^{n/\ell_C}\ne1}
       \frac{\ell_C}{1-\mu_C^{n/\ell_C}}
 +\sum_{\mu_C^{n/\ell_C}=1}\ell_C I_C=0. \tag{4}
\]

Here `I_C` is the index of `F` at any one point of `C`. For the second sum, `μ_C≠0`, so `f` locally conjugates the return germs at consecutive cycle points and the index is constant along the cycle. For the first sum, the return multiplier is the same at every cycle point, including when `μ_C=0`. This proves (4) from (2).

There is no algebraic-multiplicity factor in (4). For a fixed point of multiplicity `k`, write
`F(p+u)−p−u=u^k H(u)`, with `H(0)≠0`; its contribution is the single residue

\[
 I_F(p)=-[u^{k-1}]\frac1{H(u)}. \tag{5}
\]

Formula (3) applies only to `k=2`. For example, `F(z)=z+z³+z⁵` has index `1` at its triple fixed point zero and indices `−1/2` at each of `±i`; their sum is zero.

For a simple parabolic germ as in (3), its formal logarithm starts
`v(u)=a u²+(b−a²)u³+O(u⁴)`. Its time-form residue is therefore

\[
 \rho=\operatorname{Res}_0\frac{du}{v(u)}=1-I_F(p). \tag{6}
\]

Consequently, if all multiple fixed points in (4) are simple parabolic, their contributions can be written `ℓ_C(1−ρ_C)`. In particular, if there is exactly one multiple cycle, and it is simple parabolic of length `ℓ`, then
`ρ=1+(sum of the simple-point indices)/ℓ`. This statement concerns formal time-one return units; it does not assert analytic flow embedding.

## 3. The global calculation at −7/4

Let `q_c(x)=f_c(x)−x=x²−x+c`, and let its two roots be `α₁,α₂`. The corresponding multipliers of `f_c` are `λ_i=2α_i`, hence

\[
 \lambda_1+\lambda_2=2,\qquad \lambda_1\lambda_2=4c.
\]

Whenever both are simple fixed points of `F=f_c^{∘3}`, elementary symmetric algebra gives

\[
 \begin{aligned}
 S(c)&=\frac1{1-\lambda_1^3}+\frac1{1-\lambda_2^3}\\
 &=\frac{6(4c-1)}{(4c-1)(16c^2+4c+7)}
  =\frac6{D(c)},\qquad D(c)=16c^2+4c+7. \tag{7}
 \end{aligned}
\]

The two-simple-point derivation excludes `c=1/4` and the two complex roots of `D`. The apparent singularity at `c=1/4` is removable for the cluster sum: it has value `2/3`, which is also the index of the single coalesced fixed point of `f_{1/4}^{∘3}`. At a root of `D`, a primitive cube-root multiplier causes a collision with period-three points; the separated formula has a genuine pole and must not be used there as a finite isolated budget.

At `c*=−7/4`, equation (1) gives precisely two simple fixed points and three double fixed points. Each root of `P` is double because `P` is squarefree and coprime to `q_{c*}`. The three parabolic indices equal one number `I`, by local conjugacy along the cycle. Since `D(c*)=49`, the global budget gives

\[
 0=\frac6{49}+3I,
 \qquad I=-\frac2{49},
 \qquad \boxed{\rho=1-I=\frac{51}{49}}. \tag{8}
\]

No cubic Taylor expansion of the degree-eight return map is required. The degree-eight multiplicity count is `2+3·2=8`; the index count is `2+3` distinct contributions. Multiplying the parabolic index by six would be an error.

## 4. Why the denominator and the discriminant both equal 49

Put `Φ₃(x,c)=(f_c^{∘3}(x)−x)/q_c(x)`. It is a monic polynomial of degree six in `x`. In the quotient by `q_c`,

\[
 \Phi_3(x,c)\equiv1+2x+4x^2\equiv1-4c+6x.
\]

Thus its resultant with the fixed-point polynomial is

\[
 \operatorname{Res}_x(q_c,\Phi_3)
 =\prod_{q_c(\alpha)=0}(1+2\alpha+4\alpha^2)
 =16c^2+4c+7=D(c). \tag{9}
\]

This identity is polynomial in `c`, including exceptional parameters. At `c*`, the factorization `Φ₃=P²` gives

\[
 D(c*)=\operatorname{Res}(q_{c*},P)^2. \tag{10}
\]

Now order the roots so that `f_{c*}(r_i)=r_{i+1}`, cyclically. Since both polynomials are monic and their degrees have even product,

\[
 \operatorname{Res}(q_{c*},P)
 =\operatorname{Res}(P,q_{c*})
 =\prod_{i=1}^3 q_{c*}(r_i)
 =\prod_{i=1}^3(r_{i+1}-r_i).
\]

Squaring this cyclic gap product gives every unordered pair gap exactly once, squared. Therefore

\[
 \boxed{D(c*)=\operatorname{Res}(q_{c*},P)^2
        =\operatorname{disc}(P)=49.} \tag{11}
\]

One can check the unsquared sign and value without radicals: `P mod q_{c*}=x+5/2`, so `Res(q_{c*},P)=q_{c*}(−5/2)=7`.

Equations (7)–(11) are the precise connection: the remaining fixed-point indices supply a denominator that is a collision resultant; the square factorization turns it into a squared gap norm; for a three-cycle this is exactly its discriminant. The number `51/49` is thus forced by a global residue budget together with the period-three factorization. No general rule equating an iterative-residue denominator with a field discriminant is asserted.

A sharp scaling hostile is `y=2x`. The monic orbit polynomial becomes `8P(y/2)=y³+y²−9y−1`, with discriminant `64·49=3136`, while the conjugacy-invariant index and iterative residue remain `−2/49` and `51/49`. The underlying number field is unchanged and still has field discriminant 49; this hostile concerns the coordinate-dependent polynomial discriminant. Equation (11) is stated in the specified centered monic quadratic coordinate; the local index and iterative residue are conjugacy invariant.

## 5. The inherited cyclotomic action, with its Galois interpretation

Set `t=−1/2−r` and `M(t)=t³+t²−2t−1`. Directly,

\[
 P(-\tfrac12-t)=-M(t),\qquad
 t\longmapsto1-t-t^2\pmod M. \tag{12}
\]

If `ζ` is a primitive seventh root of unity and `t=ζ+ζ^(−1)`, then
`ζ³ M(ζ+ζ^(−1))=1+ζ+⋯+ζ⁶=0`. Also

\[
 1-t-t^2=t^3-3t=\zeta^3+\zeta^{-3}\pmod M. \tag{13}
\]

Thus the action is specifically the Galois automorphism `ζ↦ζ³` restricted to the real cubic field. The exponent classes modulo sign cycle as
`[1]→[3]→[2]→[1]`. Equivalently it is the inverse of exponent doubling, exactly as in the inherited seventh-root bridge. The polynomial `M` is irreducible by its rational-root test, and its other roots `t²−2` and `1−t−t²` lie in `Q(t)`, giving the cyclic cubic splitting field directly.

The cyclepoint index belongs to this field because it is a rational expression in the local return coefficients. Galois conjugation permutes these coefficients by permuting the three roots, while local conjugacy makes their indices equal. Hence this index descends to `Q`. This explains rationality; equation (8) determines the value. The separate [rational neutral-cycle lemma](../quadratic_cycle_phase_20261015b.md), Section 2, Lemma R1, already proves that a rational quadratic has no rational neutral cycle of distinct period at least two, by a complete valuation argument; its only neutral rational cycles are the fixed points `(x,c)=(1/2,1/4)` and `(−1/2,−3/4)`. The explicit action is a Galois reading of the inherited cyclotomic model, not a new classification of quadratic cycles or an analytic conjugacy on a neighborhood of the whole map.

## 6. A second global check: coalescing cycle indices

Use the trace parameter `v` with `c=−(7+v²)/4`. For `s=(−1+v)/2`, the inherited cubic-cycle factors are

\[
 P_s(x)=x^3-sx^2+(c-s-1)x-cs-c-1,
 \qquad \Phi_3=P_sP_{-1-s}.
\]

Writing `A(v)=v²+v+7`, the same cyclic-gap identity, or direct polynomial division, gives

\[
 \operatorname{Res}(q_c,P_s)=A(v),\qquad
 \operatorname{disc}(P_s)=A(v)^2,\qquad
 D(c)=A(v)A(-v).
\]

These are polynomial identities in `v`; at degenerate parameters their values may be zero. The multipliers are eight times their root products, namely

\[
 \mu_\pm=1-v^2\pm v(7+v^2).
\]

The branch with multiplier `1−v²−v(7+v²)` obeys `1−μ=vA(v)`. For sufficiently small nonzero `v`, the two cycles are distinct and their fixed points under `f_c^{∘3}` are simple. Their two per-cyclepoint indices satisfy

\[
 \frac1{1-\mu_+}+\frac1{1-\mu_-}
 =-\frac2{v^4+13v^2+49}
 =-\frac2{D(c)}. \tag{14}
\]

As `v→0`, one root from each cubic approaches each root of `P`. A small fixed contour about any one limiting root contains exactly those two nearby simple fixed points. Its residue integral varies continuously, so the sum in (14) tends to the parabolic index `−2/49`. This also checks the sign and the weight in (8). It does not make either individual divergent summand finite.

## 7. The same 49 is the local unfolding scale

Fix a root `r` of `P` and keep it fixed when differentiating the parameter. Write

\[
 a_r=\tfrac12\partial_x^2 f_{c*}^{\circ3}(r),\qquad
 \alpha_r=\left.\partial_c f_c^{\circ3}(r)\right|_{c=c*}.
\]

Exact reduction modulo `P` gives

\[
 a_r=14r^2-28r-\tfrac{21}{2},\qquad
 \alpha_r=2r^2-\tfrac{11}{2},\qquad
 \boxed{a_r\alpha_r=49}. \tag{15}
\]

With `ε=c−c*` and the fixed coordinate `z=a_r(x−r)`, the return map is

\[
 G_\varepsilon(z)=z+z^2+49\varepsilon
       +O(z^3,\varepsilon z,\varepsilon^2). \tag{16}
\]

This scale is also consistent with the exact multiplier splitting in Section 6: `v²=−4ε` and `μ_±=1±7v+O(v²)`, while the local equation predicts splitting `1±2√(−a_r α_r ε)+O(ε)`.

**Proposition 3 (local delay only).** Choose a sufficiently small fixed real `η>0`, start `z₀=−η`, and let `N_ε` be the first iterate of `G_ε` reaching `z≥η`. Then, as `ε↓0`,

\[
 \sqrt\varepsilon\,N_\varepsilon\longrightarrow\frac\pi7. \tag{17}
\]

The corresponding number of original quadratic steps for these complete returns is `3N_ε`, with leading coefficient `3π/7`.

**Proof.** Put `Q_ε(z)=z²+49ε`. On a sufficiently small fixed neighborhood, (16) implies
`G_ε(z)−z=Q_ε(z)(1+O(|z|+ε))`, so increments are positive and comparable to `Q_ε`; shrinking the neighborhood also ensures `G_ε′>0`. Fix `0<δ<η`. The orbit portions outside `[-δ,δ]` have `O_δ(1)` steps uniformly as `ε↓0`, since their increments are bounded below by a positive multiple of `δ²`.

On the central portion, use

\[
 T_\varepsilon(z)=\frac1{7\sqrt\varepsilon}
             \arctan\frac{z}{7\sqrt\varepsilon}.
\]

An increment of size `O(Q_ε(z))` changes `Q_ε` by relative `O(δ)` for `|z|≤2δ` and sufficiently small `ε`. Since `T_ε′=1/Q_ε`, each central step changes `T_ε` by `1+O(δ)`, with constants independent of `ε,δ`. The entry and exit points of this portion are `−δ+O(δ²+ε)` and `δ+O(δ²+ε)`. For fixed small `δ` they stay on their respective sides of zero as `ε↓0`; hence their total clock difference, multiplied by `√ε`, tends to `π/7`. Summation bounds the liminf and limsup of `√ε N_ε` within relative `O(δ)` of `π/7`. Finally let `δ↓0`. ∎

This is a specified local bottleneck, not global escape of the critical orbit. For these real parameters `c∈[−2,1/4]`, put `β=(1+√(1−4c))/2≤2`. The interval `[-β,β]` contains zero and is forward invariant because `c≥−β` and `β²+c=β`; the critical orbit therefore does not escape past 2. No full-period-window itinerary or global delay is claimed.

## 8. Scope, controls, and stopping point

The new transfer is: global return polynomial → local parabolic index via the finite residue budget. It retains the return-time unit and primitive cycle length; it forgets full local higher coefficients and analytic conjugacy data. The resultant calculation adds a genuine arithmetic explanation for this particular denominator. The cyclotomic model supplies an explicit Galois action and rationality, not a claim of analytic flow embedding, global Hamiltonian integrability, a Keller map, or a major-conjecture consequence.

Run `python3 verify_period_three_index.py` and `python3 -O verify_period_three_index.py`. The exact rational/symbolic controls cover the stated factorization, discriminant/resultant chain, fixed-point symmetric identities, the trace family, quotient-ring Galois action, one direct local-derivative comparison, multiplicity and rational-map hostiles, the unfolding product `a_r α_r=49`, and the removable boundary `c=1/4`. Explicit exceptions remain active with optimization. Finite checks support the formulas; the residue theorem and the proofs above establish the unrestricted statements.
