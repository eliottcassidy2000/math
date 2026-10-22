# MRB residue data versus the actual defect

Status: PROVED analytic obstructions and identities; VERIFIED interval certificates; OPEN modular explanation of the actual MRB defect.

## Result in one paragraph

The eighth residue has genuine honeycomb and CM identities, and the MRB constant has a genuine Abel–Plana/odd-Laplace decomposition. These facts do not yet connect the actual defect to a modular form. There is a sharp obstruction: even all negative-axis branch jumps, any prescribed finite initial segment of the MRB summands, and the complete positive-axis inverse-power asymptotics fail to determine the alternating constant in an Abel–Plana-admissible class. Explicit entire perturbations preserve those data while shifting the constant arbitrarily. The numerical defect and the next greedy index are certified below by an Euler-transform interval computation with an analytic remainder below 2.654×10^-84.

## 1. Objects, inheritance, and exact bridge inventory

Let

- f(z)=exp(Log(z)/z)−1, using the principal logarithm;
- C=Σ_{n≥1}(−1)^n f(n);
- A(m)=(2/π)m^(−1/m)sin(π/m), m≥4;
- D=A(8)−C;
- k=√2−1 and q=exp(−π√2).

The closest proved mechanisms are alternating Abel–Plana summation, the principal branch jump, the AGM transformation, and the honeycomb connective-constant theorem. The hostile examples are entire functions invisible to branch-jump data. The corrected near miss is D≈q^4, not D=q^4. The useful sidecar is the full right-half-plane function, rather than its negative-axis discontinuity or one residue.

No MRB entry was recovered in the local repository HEAD by targeted filename searches and the constants/frontier routes. The supplied `documentMRBconnections.pdf` discusses the gain-law interpretation and paired convergence; it is not the contour paper. The supplied `document24.pdf` is the square-block paper; its sign issue is handled in the separate square-block audit. The mathematical assertions here are derived independently below. Burns’s separately retrieved September 2026 manuscript is primary-source context for the already-known odd-Fourier and residue constructions, not evidence that D is modular.

Live concept board:

1. Alternating sum: full positive-integer profile; observable C.
2. Branch jump: negative-axis boundary difference; observable A(m).
3. CM modulus: elliptic-period ratio; observable k and q.
4. Honeycomb SAW: exponential combinatorial growth; observable μ.
5. Odd Laplace modes: smooth continuous density g(t); algebraic mode tail.
6. Entire perturbations: invisible to branch jumps; nontrivial alternating functional.

The arrows A(8)→k and A(8)→μ preserve exact algebraic values after normalizing π. They discard the remaining values f(n), the full contour profile, and the alternating functional. No map transporting the MRB defect to a modular coefficient has been supplied.

## 2. The branch jump is exact, but not a complete reconstruction datum

For x>0, the boundary values are

f(−x±i0)=x^(−1/x) exp(∓iπ/x)−1.

Consequently

f(−x−i0)−f(−x+i0)=2i x^(−1/x)sin(π/x)=iπ A(x).  (2.1)

Thus the residue magnitudes have a genuine analytic origin, not merely a decimal fit. But the boundary values individually grow extremely rapidly as x→0+: their modulus before subtracting one is x^(−1/x). A naive keyhole-contour argument that suppresses the small circle solely because its radius tends to zero is invalid. Equation (2.1) also does not determine f: one may add an entire function without altering any branch jump.

This last issue admits an exact, much stronger obstruction.

### Theorem: branch jumps plus any finite prefix do not determine the sum

For every integer K≥1, put

h_K(z)=e^(−z)∏_{j=1}^K(z−j),

and f_{K,c}(z)=f(z)+c h_K(z), for any real c. Then:

1. f_{K,c}(j)=f(j) for j=1,…,K.
2. Every negative-axis branch jump of f_{K,c} equals the corresponding jump of f.
3. h_K(x) is smaller than every inverse power of x as x→+∞. Therefore the two functions have identical complete asymptotic expansions in powers of 1/x and log x on the positive axis.
4. The perturbation is admissible for alternating Abel–Plana: it is analytic in the right half-plane, decays exponentially on the positive real axis, and has only polynomial growth along a vertical line.
5. Nevertheless,

Σ_{n≥1}(−1)^n f_{K,c}(n) = C+c(−1)^(K+1)K!/(e+1)^(K+1).  (2.2)

Proof of (2.2): let x=−e^(−1). The first K perturbation terms vanish, and

Σ_{n≥K+1}x^n∏_{j=1}^K(n−j)
=K! Σ_{n≥K+1}binom(n−1,K)x^n
=K! x^(K+1)/(1−x)^(K+1).

The final expression equals the coefficient in (2.2), which is nonzero. Since c is unrestricted, the alternating constant can be shifted by any prescribed real amount while preserving all the listed data. □

For K=1 this is especially transparent: h_1(z)=(z−1)e^(−z), and the shift is c/(e+1)^2. Its contour-side identity is

Im h_1(1+it)=e^(−1)t cos t,

∫₀∞ t cos t/sinh(πt) dt = (1/4)sech²(1/2),

so the Abel–Plana integral changes by exactly c/(e+1)^2 as well.

Scope: this theorem refutes determination from branch jumps/residue data and finite positive samples, not the possibility of an exact identity special to the full function exp(Log z/z)−1. It does not assert that the perturbed functions remain nth-root profiles, nor that they satisfy f(z)=O(log|z|/|z|) uniformly throughout the right half-plane. The missing global growth/profile sidecar is exactly what a successful uniqueness or contour reconstruction theorem would need to exploit.

## 3. A convergent contour representation and its actual modes

Alternating Abel–Plana, shifted to the regular point z=1, yields

C=∫₀∞ g(t)/sinh(πt) dt,

g(t)=Im[exp(Log(1+it)/(1+it))−1].  (3.1)

The endpoint behavior is g(t)=t+O(t³), and at infinity g(t)=O(log t/t), so (3.1) is absolutely convergent. This avoids crossing the singular point z=0.

The exact kernel identity

1/sinh(πt)=2Σ_{j≥0}exp(−π(2j+1)t)

gives

C=Σ_{j≥0}F_{2j+1},

F_ell=2∫₀∞ exp(−π ell t)g(t)dt, for odd ell≥1.  (3.2)

Absolute interchange follows from integrating |g(t)|/sinh(πt). A first-quadrant rectangle, followed by integration by parts, also gives

F_ell=2∫₁∞f(x)cos(π ell x)dx
=−2/(π ell)∫₁∞f′(x)sin(π ell x)dx.  (3.3)

These are two representations of the same modes, not independent modular data.

Watson’s elementary Laplace scaling applied to g(t)=t−t³/2+O(t⁵) gives

F_ell = 2/(π² ell²)−6/(π⁴ ell⁴)+O(ell^−6).  (3.4)

Therefore, for J retained odd modes,

C−Σ_{j=0}^{J−1}F_{2j+1} = 1/(2π²J)+O(J^−2).  (3.5)

This is a precise obstruction to treating q^4 as the leading term omitted by the natural odd-mode truncation: that truncation has an algebraic tail, not a geometric nome tail. It does not exclude a different, as-yet-undiscovered modular transformation.

Numerical controls (100-decimal working precision):

- F_1 = 0.141552078623057607079056043660564…
- F_3 = 0.021620631821963423061242124317148…
- F_5 = 0.007997816172196490016772048250549…
- F_201 = 0.000005015738287589290917333072853449…
- F_201·π²·201²/2 = 0.99999247576999306086706235373046…

A second hostile control tests the most literal exponential-scale explanation: truncate the integral at T=4√2, so exp(−πT)=q^4. The discarded tail satisfies numerically

∫_T∞g(t)/sinh(πt)dt / q^4 = −0.212046476966226782476635657038…,

not −D/q^4=−0.988450743895774519977687688583…. Thus A(8) is not the integral truncated at the height whose kernel exponential is q^4.

For a simple analytic upper bound on any such vertical tail, |f(1+it)|≤e^(π/2)−1, so

|∫_T∞g(t)/sinh(πt)dt|
≤2(e^(π/2)−1)e^(−πT)/[π(1−e^(−2πT))], T>0.  (3.6)

The exponent here is imposed by the arbitrary cutoff T. It is not a derivation of the distinguished value T=4√2.

### An exact positive-theta-tail obstruction

Using the conventional theta series,

[theta_3(0,q)+theta_4(0,q)−2]/4 = Σ_{r≥1}q^(4r²) > q^4, for 0<q<1.

The certified MRB defect instead has 0<D<q^4. Therefore the unweighted positive even-theta tail cannot equal D. This is a sign obstruction, not a numerical fit. Signed or normalized theta combinations are not excluded; a proposed normalization must be fixed before comparing with the MRB value.

## 4. The r=2 singular modulus is exactly correct

Use the modulus convention

K(k)=∫₀^(π/2)(1−k²sin²θ)^−1/2 dθ,
k′=√(1−k²), and K(k′)/K(k)=√r.

For k=√2−1 one has k′=√(2k) and 1+k=√2. Let M be the arithmetic–geometric mean. Its first step and homogeneity give

√2 M(1,k)
=√2 M((1+k)/2,√k)
=M(1,√(2k))
=M(1,k′).

Since K(k)=π/[2M(1,k′)] and K(k′)=π/[2M(1,k)], it follows exactly that K(k′)/K(k)=√2. Hence

q=exp(−π√2),
k=theta_2(0,q)^2/theta_3(0,q)^2,
λ(i√2)=k²=3−2√2.  (4.1)

Do not confuse the elliptic modulus k with the parameter m=k²: the latter is what mpmath.ellipk accepts. In particular, calling ellipk(√2−1) would test the wrong value.

The eighth residue obeys

k=2^(1/4)(πA(8))²,
A(8)=2^(−1/8)√k/π.  (4.2)

These exact statements have no C on either side. They fix the residue’s arithmetic interpretation but supply no estimate on its difference from C.

## 5. The exact honeycomb identity and its limit

Duminil-Copin and Smirnov proved that the honeycomb self-avoiding-walk connective constant is

μ_hex=√(2+√2)=2cos(π/8).

Since sin(π/8)cos(π/8)=√2/4,

A(8)=2^(1/8)/(π μ_hex).  (5.1)

This is exact. The map is a trigonometric/algebraic identification of one residue with the reciprocal of a combinatorial growth constant. No correspondence between MRB summands and honeycomb walks, no transported alternating observable, and no honeycomb remainder estimate for D has been established. The entire-perturbation theorem shows precisely why the shared A(8) value alone cannot supply one.

## 6. Rigorous numerical enclosure without trusting quadrature

The script `certify_mrb.py` uses mpmath interval arithmetic at 120 decimals for a finite Euler transform, plus the following analytic remainder.

For integers M≥4 and N≥1, set

S_M=Σ_{j≥0}(−1)^j f(M+j),
Δf(x)=f(x+1)−f(x),
L=log(3M/2)+π/6,
b=2L/M.

Repeated averaging gives

S_M=Σ_{j=0}^{N−1}(−1)^j Δ^j f(M)/2^(j+1)+R_(M,N),

|R_(M,N)|≤2N!e^b [L/M^(N+1)+(L/N+1/N²)/M^N].  (6.1)

Proof of the bound: on |z−x|=x/2, x≥M, one has |arg z|≤π/6, |z|≤3x/2, |z|≥x/2, and

|f(z)|≤(2e^b/x)(log(3x/2)+π/6).

Cauchy’s estimate then gives

|f^(N)(x)|≤2^(N+1)N!e^b[log(3x/2)+π/6]/x^(N+1).

Express Δ^N f(M+j) as the N-fold unit-cube integral of f^(N), use the decreasing majorant, and bound the resulting sum by its first term plus its integral. The Euler remainder has the prefactor 2^−N, giving (6.1). □

At M=256,N=100 the bound is less than 2.654×10^-84. Interval arithmetic encloses every finite radical, difference, π, logarithm, sine, and exponential. The resulting interval certifies, in particular,

1.891807127880073527054533938441983110587504636897727612419589120×10^-8
< D <
1.891807127880073527054533938441983110587504636897727612419589122×10^-8.

It also certifies

D/q^4 = 0.988450743895774519977687688583150115109467139287502977396141875111…,

with interval width below 2.8×10^-76. Thus D=q^4 is decisively false; the relative discrepancy is about 1.1549256104 percent.

The same interval determines the second nearest-residue index n=105719004. Both midpoint margins are strictly positive:

D−[A(n)+A(n+1)]/2 > 1.1531280076483094×10^-16,

[A(n−1)+A(n)]/2−D > 6.36339108938903×10^-17.

Together with strict monotonicity of A(m) for m≥4, these certify the unique nearest integer index. For completeness, A′(x)/A(x)=[log x−1−πcot(π/x)]/x². The inequality cot t>1/t−t for 0<t≤π/4 follows by differentiating t cos t−(1−t²)sin t. Thus πcot(π/x)>x−π²/x>log x−1 for x≥4: the final difference x+1−log x−π²/x is positive at 4 and strictly increasing. Hence A′(x)<0. The inverse root for the second step is numerically

m*=105719003.85560256733977518406685302945…,

and the interval certificate proves A(105719003)>D>A(105719004). Therefore 105719004 is the ceiling of the inverse root and the nearest-amplitude choice, not a floor-overestimate choice. By contrast, the initial continuous inverse is 8.00000099077209685…, so index 8 is a floor choice and overshoots C. A fixed floor rule or a fixed ceiling rule does not reproduce both selections. Nearest amplitude is the convention tested here.

The new residual magnitude is

D−A(105719004)=2.5839445781801954275830241249248822113390313828449…×10^-17.

As a cross-check only, 100-decimal adaptive quadrature of (3.1) and accelerated summation agree with the interval result. The script keeps these numerical checks explicitly separate from the certificate.

## 7. Cheapest next decisive test and stopping condition

A proposed modular explanation must now exhibit at least one new piece of data not preserved by h_K: for example, a full-profile modular transformation, a differential equation with a uniqueness condition, a contour identity with all boundary terms controlled, or a map of the entire MRB sequence into a modular family. It must then derive the defect coefficient and quantify a remainder without fitting the coefficient to D.

Merely identifying A(8), rediscovering its CM modulus, reproducing the honeycomb constant, or observing q^4’s numerical scale cannot cross this obstruction. This pass stops at the proved obstruction and certified numbers, rather than promoting an unsupported explanation. The arithmetic nature of C and a modular formula for D remain OPEN.

## Sources and reproduction

- Richard E. Crandall, *Unified algorithms for polylogarithm, L-series, and zeta variants* (2012), §7.5: https://www.marvinrayburns.com/UniversalTOC25.pdf . Primary context for the existing diagonal eta-derivative representation, not used in the interval proof here.
- Marvin Ray Burns, *A Convergent Greedy Residue with Geometry: Ramanujan and Michele Nardelli’s Version* (September 2026), §§3–4,9: https://www.academia.edu/175280204/A_Convergent_Greedy_Residue_with_Geometry_Ramanujan_and_Michele_Nardellis_Version . Primary context for the existing odd-mode/residue bridge and the explicit statement that modularity of the residual remains open. The author’s listing is https://indiana.academia.edu/MarvinBurns .
- NIST DLMF §19.8, particularly the AGM formula 19.8.5: https://dlmf.nist.gov/19.8#E5 .
- NIST DLMF §23.15, nome, modulus, and lambda conventions: https://dlmf.nist.gov/23.15 .
- H. Duminil-Copin and S. Smirnov, *The connective constant of the honeycomb lattice equals √(2+√2)*, Annals of Mathematics 175 (2012): https://arxiv.org/abs/1007.0575 ; author-hosted paper: https://www.unige.ch/~smirnov/papers/saw.pdf .

Reproduce from the workspace root:

    python session31-mrb/certify_mrb.py > session31-mrb/certify_mrb.out

All interval assertions and all independent hostile-control checks pass. No repository files or git state were modified by this worker.
