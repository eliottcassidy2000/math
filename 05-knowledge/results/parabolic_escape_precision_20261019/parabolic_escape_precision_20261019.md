# An exact clock for the quadratic cusp: bounded error, phase, and band widths

2026-10-19. **PROVED HERE, independently and elementarily:** the exact clock identity, an explicit bounded-error estimate, the threshold-phase limit, adjacent-band asymptotics, and the Cauchy occupation law. **CITED / KNOWN:** the leading π-law and its bounded-error refinement. **OPEN in this note:** a proposed next-order phase coefficient. These results are not claimed to be new in the literature.

## 1. Scope, inheritance, and attribution

For ε>0, let

\[
x_0=0,\qquad x_{k+1}=x_k^2+\tfrac14+\epsilon,
\qquad N_R(\epsilon)=\min\{k\ge0:x_k>R\},\quad R>\tfrac12.
\]

The inequality is **strict**. This convention matters at the band boundaries. The familiar escape radius is R=2, but the following real one-dimensional arguments need only R>1/2.

Aaron Klebanoff, *π in the Mandelbrot Set*, Fractals 9(4) (2001), 393–402, proves √ε N₂(ε)→π. His paper describes Dave Boll's 1991 experiments and earlier real-axis heuristics. Primary PDF: <https://www.doc.ic.ac.uk/~jb/teaching/jmc/pi-in-mandelbrot.pdf>.

Thies Brockmoeller, Oscar Scherz, and Nedim Srkalovic, *Pi in the Mandelbrot set everywhere*, arXiv:2505.07138v1 (2025), Theorem 6, already gives an O(1) additive remainder at primitive parabolic parameters. The authors credit a sharper 2023 result of Paul Siewert. Their Sections 2–3 relate escape to perturbed Fatou-coordinate phases. The 2023 thesis itself was not obtained in this pass. Primary paper: <https://arxiv.org/html/2505.07138v1>. Its broader complex-parameter statements are not needed below.

The repository's period-three graph note distinguishes the multiplier-one orbit at c=−7/4 from the rational orbit at c=−29/16; see `05-knowledge/results/period3_graph_square_bridge_20261019.md`. That is a genuine parabolic connection through a local return map, but neither the critical orbit, the first-return scaling, nor its parameter unfolding can be transferred just by matching the word “parabolic.” The current note treats only the explicit cusp family at 1/4.

Concept board: exact iteration clock / summable defect / attracting-repelling phase / integer staircase / parameter resolution / base digits. The key hostile is replacing a real asymptotic by an exact floor or a congruence.

## 2. The exact one-step clock

Put δ=√ε, y=x−1/2, a=−1/2, b=R−1/2, and

\[
f_\delta(y)=y+y^2+\delta^2,\qquad q_\delta(y)=y^2+\delta^2.
\]

Define

\[
\Psi_\delta(y)=\frac1\delta\arctan\frac y\delta
                  +\frac12\log(y^2+\delta^2).
\tag{1}
\]

For every δ>0 and y≥−1/2,

\[
\Psi_\delta(f_\delta(y))-\Psi_\delta(y)=1+D_\delta(y),
\tag{2}
\]

where

\[
D_\delta(y)=\frac1\delta\arctan\frac\delta{1+y}
 +\frac12\log((1+y)^2+\delta^2)-1.
\tag{3}
\]

Moreover,

\[
0\le D_\delta(y)\le 2(y^2+\delta^2).
\tag{4}
\]

**Proof of (2).** Directly,

\[
q_\delta(f_\delta(y))=q_\delta(y)((1+y)^2+\delta^2).
\]

The tangent subtraction formula gives the arctangent increment in (3). There is no hidden branch jump: 1+y>0, and the angle increment is in (0,π/2).

**Proof of (4).** Write w=1+y≥1/2 and call the right side of (3) dδ(w). Then

\[
\frac{d}{dw}d_\delta(w)=\frac{w-1}{w^2+\delta^2}.
\]

Thus dδ has its minimum at w=1. Set

\[
J(\delta)=d_\delta(1)=\frac{\arctan\delta}{\delta}
                 +\frac12\log(1+\delta^2)-1.
\]

Since J(0)=0 and

\[
J'(\delta)=\frac{\delta-\arctan\delta}{\delta^2},
\qquad 0\le J'(\delta)\le\delta/3,
\]

we have 0≤J(δ)≤δ²/6. On the interval between 1 and w, the denominator t²+δ² is at least 1/4 when w≤1, and at least 1 when w≥1. Therefore

\[
0\le d_\delta(w)-J(\delta)\le 2(w-1)^2.
\]

This proves (4). For y≥0 one has the better bound Dδ(y)≤y²/2+δ²/6. QED.

**The decisive telescoping.** Along an orbit from y₀ to yₙ,

\[
n=\Psi_\delta(y_n)-\Psi_\delta(y_0)-E_n,
\quad E_n=\sum_{k=0}^{n-1}D_\delta(y_k),
\quad 0\le E_n\le2(y_n-y_0).
\tag{5}
\]

The error is bounded by **distance travelled**, not time spent. This is why even a very long bottleneck does not accumulate a logarithmic full-passage error.

## 3. Uniform full-passage error

For each fixed R>1/2,

\[
N_R(\epsilon)=\frac\pi{\sqrt\epsilon}+O_R(1),\qquad\epsilon\downarrow0.
\tag{6}
\]

Indeed yN lies in (b,b+b²+δ²], so the endpoint is bounded away from 0 and bounded above. Formula (5), and

\[
\Psi_\delta(v)-\Psi_\delta(a)
=\frac\pi\delta-\frac{\arctan(\delta/v)+\arctan(2\delta)}\delta
 +\frac12\log\frac{v^2+\delta^2}{1/4+\delta^2},\quad v>0,
\tag{7}
\]

give (6) immediately.

For R=2 and 0<ε≤1/4, a deliberately simple explicit bound is

\[
-11<N_2(\epsilon)-\pi/\sqrt\epsilon<3.
\tag{8}
\]

To check it, v=yN is in (3/2,4], E≤9, the angular subtraction in (7) is between −8/3 and 0, and the logarithmic term is between (log 5)/2 and log 8. These constants are convenient, not sharp.

## 4. Exact threshold phase and adjacent-band width

For each n≥1 there is a unique εn>0 such that xₙ(εn)=R. Existence and uniqueness follow because xₙ is continuous and strictly increasing in ε≥0, with xₙ(0)<1/2<R and xₙ(ε)→∞. The εn decrease strictly to 0. Strict decrease follows from xₙ₊₁(εn)>xₙ(εn); convergence to 0 follows because every positive ε has a strictly increasing unbounded orbit.

Let δn=√εn. At εn, equality is **not** escape, so

\[
N_R(\epsilon_n)=n+1,
\qquad N_R(\epsilon)=n\ \Longleftrightarrow\
\epsilon_n<\epsilon\le\epsilon_{n-1}\quad(n\ge2).
\tag{9}
\]

Let f(y)=y+y² and let

\[
g(y)=\frac{\sqrt{1+4y}-1}{2}
     =\frac{2y}{1+\sqrt{1+4y}}
\]

be its nonnegative inverse. Define

\[
h(y)=-\frac1y+\log|y|,
\quad D_0(y)=\frac1{1+y}+\log(1+y)-1.
\]

The convergent, unperturbed endpoint sums define

\[
C_R=h(b)-h(a)
 -\sum_{k=0}^{\infty}D_0(f^k(a))
 -\sum_{j=1}^{\infty}D_0(g^j(b)).
\tag{10}
\]

**Theorem.**

\[
n-\frac\pi{\sqrt{\epsilon_n}}\longrightarrow C_R,
\tag{11}
\]

and consequently

\[
\epsilon_n=\frac{\pi^2}{(n-C_R+o(1))^2},
\qquad
\epsilon_n-\epsilon_{n+1}\sim\frac{2\pi^2}{n^3}.
\tag{12}
\]

**Proof.** At δ=0, (4) gives summability on both tails: for uₖ=fᵏ(a), Σuₖ²=−a; for vⱼ=gʲ(b), Σⱼ≥1vⱼ²=b. These identities just telescope f(y)−y=y².

At εn the orbit starts exactly at a and ends exactly at b. For every fixed m, its first m points tend to fᵏ(a), k=0,…,m−1, as n→∞. Its last m points tend to gʲ(b), j=1,…,m. This uses continuity of the positive inverse branches; for fixed m they remain positive for sufficiently small δn.

The middle block has defect at most

\[
2\bigl(y_{n-m}-y_m\bigr)
\longrightarrow2\bigl(g^m(b)-f^m(a)\bigr).
\]

The latter tends to 0 as m→∞. Since every defect is nonnegative, first taking n→∞ with m fixed and then m→∞ proves convergence of the total defect to the two sums in (10). For fixed nonzero endpoints, (7) gives

\[
\Psi_\delta(b)-\Psi_\delta(a)-\pi/\delta
\longrightarrow h(b)-h(a).
\]

Now (5) proves (11).

Write tn=π/δn=n−C_R+rn with rn→0. Then tn₊₁−tn=1+o(1), so

\[
\epsilon_n-\epsilon_{n+1}
=\pi^2\frac{(t_{n+1}-t_n)(t_{n+1}+t_n)}{t_n^2t_{n+1}^2}
\sim2\pi^2/n^3.
\]

This is a valid adjacent-difference argument because rn=o(1). Merely knowing εn∼π²/n², or even εn=π²/n²+O(n⁻³), would not justify it. QED.

**Fatou-coordinate interpretation.** The normalized attracting and repelling coordinates are

\[
\Phi_-(a)=h(a)+\sum_{k\ge0}D_0(f^k(a)),\quad
\Phi_+(b)=h(b)-\sum_{j\ge1}D_0(g^j(b)).
\]

Their Abel equations are Φ±(f(y))=Φ±(y)+1 where defined, and Φ±(y)−h(y)→0 at the relevant side of 0. Thus C_R=Φ₊(b)−Φ₋(a). In particular changing b to f(b) increases C_R by exactly 1.

## 5. The staircase cannot have a single additive limit

The cluster set of

\[
N_R(\epsilon)-\pi/\sqrt\epsilon\quad(\epsilon\downarrow0)
\]

is **exactly** [C_R,C_R+1]. In particular its liminf and limsup differ by one.

**Proof.** On the nth band, t=π/√ε ranges over [tₙ₋₁,tₙ), so n−t ranges over (n−tₙ,n−tₙ₋₁]. Its endpoints tend to C_R and C_R+1 by (11). Every interior value occurs on every sufficiently late band; the two endpoints are limits of band endpoints. QED.

Therefore a notation such as N=π/√ε+C+o(1) is false for a real parameter tending continuously to zero. A formula with a floor/ceiling needs a quantified phase error and careful handling of equality; the threshold formulation (9)–(11) is safer.

## 6. A certified value of the R=2 phase

For m≥1 let u=fᵐ(−1/2)<0 and v=gᵐ(3/2)>0, and put

\[
C^{(m)}=h(v)-h(u)+2m.
\]

Using D₀(y)=∫₀ʸ t/(1+t)² dt gives the rigorous tail enclosure

\[
C^{(m)}+\frac{u}{2(1+u)^2}-\frac v2
\ \le C_2\le\
C^{(m)}+\frac u2-\frac{v}{2(1+v)^2}.
\tag{13}
\]

This follows by bounding the attracting tail between (−u)/2 and (−u)/(2(1+u)²), and the repelling tail between v/(2(1+v)²) and v/2. The resulting interval has width O(m⁻²); it does not require extrapolating convergence.

Directed interval evaluation with m=10,000 and 50 decimal digits yields the deliberately outward-rounded enclosure

\[
\boxed{-2.463836664<C_2<-2.463836643.}
\tag{14}
\]

The companion script checks these endpoints using `mpmath.iv`. Floating-point samples and the interval certificate are labeled separately.

An immediate consequence is that for **all sufficiently small** ε,

\[
N_2(\epsilon)\in
\{\lfloor\pi/\sqrt\epsilon\rfloor-2,
  \lfloor\pi/\sqrt\epsilon\rfloor-1\}.
\tag{15}
\]

This follows from the staircase limits and the strict enclosure −3<C₂<−2 and −2<C₂+1<−1. The switch between the two alternatives is controlled asymptotically by the fractional part of π/√ε crossing −C₂−2≈0.46383665. No unqualified eventual floor identity is asserted at phases approaching that boundary.

## 7. Half passages, moving thresholds, and a Cauchy law

For a fixed real λ, let Nλ(ε) be the first n with

\[
x_n>\tfrac12+\lambda\sqrt\epsilon.
\]

For all sufficiently small ε this is a threshold above the starting point. Then

\[
N_\lambda(\epsilon)
=\frac{\pi/2+\arctan\lambda}{\sqrt\epsilon}
 +\frac12\log\epsilon+O_\lambda(1).
\tag{16}
\]

Indeed the endpoint is y=λδ+Oλ(δ²), so its Ψδ value is

\[
\frac{\arctan\lambda}{\delta}+\log\delta+O_\lambda(1).
\]

The starting value is −π/(2δ)+h(a)+O(δ²), while (5) bounds the summed defect uniformly. This proves (16).

At λ=0, the time spent approaching the fixed-point location is π/(2√ε)+(log ε)/2+O(1). Starting exactly at y=0 and passing to fixed b>0 instead takes π/(2√ε)−(log ε)/2+O(1). The full passage cancels these opposite logarithms. This diagnoses why a coarse one-sided comparison can suggest an unnecessary O(log(1/ε)) total error.

For the full orbit k=0,…,N_R−1 define its empirical scaled-position measure

\[
\mu_\epsilon=\frac1{N_R(\epsilon)}
 \sum_{k=0}^{N_R(\epsilon)-1}\delta_{(x_k-1/2)/\sqrt\epsilon}.
\]

By monotonicity, its distribution function at λ is Nλ/N_R, up to at most one index convention. From (6) and (16),

\[
\mu_\epsilon((−\infty,\lambda])\longrightarrow
\frac12+\frac1\pi\arctan\lambda.
\tag{17}
\]

Thus με converges weakly to the standard Cauchy distribution, with density 1/[π(1+u²)]. Equivalently the angle arctan((x−1/2)/√ε) becomes uniformly occupied on (−π/2,π/2). This is a distribution of time spent along one deterministic orbit, not a claim that its iterates are independent random variables.

## 8. Every base works, with a digit firewall

For any integer base B≥2 and m≥1, choose ε=B⁻²ᵐ. Then

\[
N_R(B^{-2m})=\pi B^m+O_R(1),
\qquad B^{-m}N_R(B^{-2m})=\pi+O_R(B^{-m}).
\tag{18}
\]

Decimal notation plays no dynamical role. More generally ε=A B⁻²ᵐ, A>0 fixed, gives N_R/Bᵐ→π/√A. Odd decimal perturbation exponents give alternating scale factors because of the square root.

Consequences and limits:

1. Every fixed finite base-B prefix of π is eventually recovered by the normalized approximation. This uses irrationality of π to keep that fixed prefix away from a terminating base-B boundary.
2. A uniform O(B⁻ᵐ) error means **numerical accuracy**, not guaranteed agreement of the first m digit symbols. Carries near a base-B boundary can alter earlier displayed digits.
3. Given a certified error |A−π|≤η, k fractional base-B digits are certified exactly when the whole interval [A−η,A+η] lies in one interval [jB⁻ᵏ,(j+1)B⁻ᵏ). For R=2, (8) supplies η=11B⁻ᵐ; a sharper clock enclosure improves it.
4. Even bounded additive differences between integer sequences do not preserve residues, primality, or twin-prime status. Equation (15) explicitly exhibits the phase-dependent one-or-two-unit shift here.
5. The direct iteration takes Θ(Bᵐ) steps to produce m-scale accuracy. This is an explanatory π phenomenon, not a competitive π algorithm.

## 9. Parameter precision is cubic in the desired escape count

At escape count n, the band width is asymptotic to 2π²/n³, while its distance to the cusp is asymptotic to π²/n². Thus resolving an individual band requires absolute parameter resolution on the n⁻³ scale, or relative resolution on the n⁻¹ scale. A uniform guarantee is impossible at the exact threshold: an arbitrarily small perturbation can change the escape integer there.

For ε=B⁻²ᵐ, the relevant scale is B⁻³ᵐ. Storing c=1/4+ε without losing ε and resolving which escape band it occupies are different demands. Working directly with yₖ₊₁=yₖ+yₖ²+ε avoids the initial subtraction/cancellation of adding a tiny ε to 1/4. Interval propagation, increasing precision when an interval straddles the threshold, is the honest way to certify an escape count.

## 10. A next-order question, not yet a theorem here

The exact complex fixed points are ±iδ, with multipliers 1±2iδ. The local Koenigs time period at the upper fixed point has expansion

\[
\frac{2\pi i}{\log(1+2i\delta)}
=\frac\pi\delta+\pi i+\frac\pi3\delta+O(\delta^2).
\]

This motivates the threshold-phase hypothesis

\[
n=\frac\pi{\delta_n}+C_R+\frac\pi3\delta_n+o(\delta_n).
\tag{H}
\]

The coefficient is supported by small numerical threshold computations, but the local complex period does not by itself prove the real endpoint matching error is o(δ). A valid proof must control the perturbation of both normalized Fatou coordinates, or sharpen the elementary split-orbit comparison. No stronger spacing-error rate is claimed from (H).

The discrete-versus-continuous obstruction is useful: the model clock (1) has the correct divergent period and logarithmic residue, while the discrete multiplier contains the next correction. A formal vector-field embedding or a local normal form must not be promoted to a global embedding of the quadratic polynomial.

## 11. Reproduction and decisive tests

Run the companion `verify_parabolic_escape_precision_20261019.py`. It performs a directed-interval phase enclosure using (13), high-precision threshold bisections, numerical adjacent-spacing controls, and base-2/base-10 escape samples. The exact proofs above are independent of those finite checks. The controls distinguish threshold equality from strict escape, use y=x−1/2 directly, and label uncertified sample values as numerical.

The strongest next useful test is not a larger digit table. It is a proof or refutation of (H), with an explicit endpoint-matching error, followed by whether the expansion gives certified band membership from endpoint data without Θ(n) iterates.
