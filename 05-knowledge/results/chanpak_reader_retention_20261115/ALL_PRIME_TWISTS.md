# Finite-jet transversality for the alternating floor chains

Follow-up to the alternating AB/BA analysis for THM-5202 (Wave 28)
and THM-5212(A) (Wave 29). This is a theorem about depth-averaged
twists, not about mixing between independent words at a fixed size.

Use T=[[0,1],[8,9]], with initial column v=(9,89) or (5,53).
The first coordinate of T^j v is the Hamiltonian count of the
corresponding alternating chain at depth j+1.

**Theorem.** For every fixed odd prime p, the complete primitive
p^r twist is O_p(p^(−r/2)) for both chains, uniformly in the primitive
numerator. There is an effective constant depending only on a
two-by-two logarithmic matrix. For the discriminant prime p=113 there
is the stronger conclusion that each chain is exactly equidistributed
modulo every 113^r: each nontrivial primitive twist is zero.

This shows that square-root Fourier decay at all single-prime towers,
or exact equidistribution at a ramified prime, need not come from a
large family of Hamiltonian counts. The two chains have only O(log X)
distinct counts up to X. Earlier dyadic analysis has its own residual
mod-16 obstruction.

## Odd-prime proof: bound the critical jets

Since det T=−8, choose m with T^m=I mod p. Put L=log(T^m), an integral
matrix divisible by p. For each phase i mod m the count as a function
of t is F_i(t)=e₁ᵀ exp(tL) T^i v on Z_p. These are convergent power
series with integral coefficients, with all formal derivatives
valid. The same holds after Taylor shift by t_0∈Z_p.

The matrix B with rows e₁ᵀL and e₁ᵀL² is invertible. Indeed over a
splitting field T has two distinct eigenvalues α,β, satisfying
x²−9x−8=0; both eigenvectors have first coordinate nonzero. Its two
logarithmic eigenvalues log(α^m), log(β^m) are nonzero and unequal.
If either log were zero, or if they were equal, then α^m, β^m, or
(α/β)^m would be 1, because log is injective on 1+p times the ring
of integers of the splitting field for odd p. Each alternative
is excluded at the real embedding: |α|>1>|β|>0. In that eigenbasis
det B is a nonzero product of the two eigenvector coordinates, the
two logarithmic eigenvalues, and their difference.

Let d=v_p(det B), a computable nonnegative integer. For all i and t
the vector exp(tL)T^i v is primitive, since v is primitive and both
factors are invertible over Z_p. Consequently at least one of F_i'(t)
and F_i''(t) has valuation at most d. If both entries of B times that
vector were divisible by p^(d+1), multiplication by the integral
adjugate would force both entries of the vector to be divisible by p.

For primitive conductor p^r and s=floor(r/2), split any complete
t sum mod p^r as t=x+p^(r−s)y. Higher Taylor terms vanish in the
phase. The complete inner y sum cancels unless F_i'(x)=0 mod p^s.
If s>d, partition the possible x into residue classes modulo p^(d+1).
In any class with a candidate x₀, set e=v_p(F_i''(x₀))≤d. For x,y in
that class, expansion at y gives
v_p(F_i'(x)−F_i'(y))=v_p(x−y)+e:
the valuation of F_i'' stays e and the integral higher coefficients
make the nonlinear remainder strictly smaller p-adically. Therefore
all solutions in this class agree modulo p^(s−e). The total measure
of candidates is at most p^(2d+1−s). When s≤d the estimate is replaced
by 1. Averaging the equally weighted m phases proves the theorem,
with for instance C_p=p^(2d+3/2).

The constant does not require infinite-precision logs. Let U_m be the
AB Lucas sequence (U_0=0,U_1=1,U_(m+2)=9U_(m+1)+8U_m). Then
`d = v_p(det(T^m−I)) + v_p(U_m)`.
To see this, L=(T^m−I)J for an integral J≡I mod p, and L=aI+bT
since it commutes with T and the centralizer of this cyclic
two-by-two matrix is its polynomial algebra. Thus det B=b det L.
In the splitting field b=(log α^m−log β^m)/(α−β). The logarithm
is an isometry on 1+p times the ring of integers, making v_p(b)
equal to the valuation of (α^m−β^m)/(α−β)=U_m. This yields the
finite integer formula.

## Ramified prime: the nilpotent jump is useful

At p=113 the repeated characteristic root is λ=61 mod p, of order
56. Put M=T^56, B₀=M−I. Modulo p one has
B₀=[[57,25],[87,56]], B₀²=0, and B₀≠0.
Its image is the λ-eigenline of T and has nonzero first coordinate.
For the two starting columns, e₁ᵀB₀v is 26 and 28 mod p, respectively.
Thus it remains nonzero after any phase shift T^i.

By the binomial theorem M^(p^s)≡I+p^s B₀ mod p^(s+1) for every s≥0.
For s=1, the middle binomial terms have a factor p and B₀² has a
factor p; the final B₀^p is divisible by p² because B₀² is
divisible by p. Induction for s≥1 is then the ordinary binomial
lifting from I+p^s B₀. Also M^t B₀=B₀ mod p for every nonnegative t.

It follows that changing t by p^s k changes the count in phase
i+j·56, as a function of j=t, by
p^s k e₁ᵀ B₀ T^i v modulo p^(s+1). At s=0 this is true even when
k is not one, by B₀²=0 mod p; at higher s it follows from the lift.
This gives exactly one choice of depth digit to hit every prescribed
next output digit. Each of the 56 phases individually visits every
residue exactly once as t ranges mod p^r. Matrix order mod p^r is
56·p^r, so this is an exact periodic histogram and not only an
analytic-limit description. This is a rank-two analogue of a locally
permutive digit rule; nilpotent parts need not be information loss
when their coupling to the output is kept.

In fact the exact vanishing survives a mixed conductor 113^r q
for any q coprime to 113 whose count period has no factor 113. Fix j
modulo both 56 and the q period. The remaining depth index for the
113 part still runs uniformly mod 113^r, by CRT, and the q character
is constant. This applies to every q=2^a 3^b 5^c by the earlier
explicit period bounds.

## Randomized older-repo tangents

Seed 37329, update s←1664525s+1013904223 mod 2^32, sample index s mod
current pool size and remove the selected entry. Initial ordered pool:
Rule 30, Yang–Mills, square-sum, aliquot, transversality, Leech lattice,
MRB constant, Hadwiger Nelson, Selberg, Gauss sums, Jordan, primitive
triples, FRACTRAN, profinite, p-adic logarithm, Cayley Dickson.
The six draws were FRACTRAN, MRB constant, transversality, Hadwiger Nelson,
primitive triples, and Rule 30.

- [Four machine lenses](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/lenses_tm_ca_fractran_tournaments_20261101.md):
  preserve the observable, not only the update. Here the two rows of
  logarithmic jets are the needed observable; rank is not sufficient.
- [September Collatz transversality foundry](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/collatz_procgen_20260922_transversality_foundry.md):
  exact 2-adic repetition and isometry arguments suggested separating
  time congruence into residue phases and then lifting digits.
- [MRB draw → October tangent-Legendre certificate](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session54_reframed_numerology_20261025/tangent-cm/README.md):
  an exceptional arithmetic parameter merits a separate exact analysis,
  not extrapolation from generic cases. Here and there the equations
  differ; the comparison is methodological.
- [Hadwiger-Nelson draw → D19 parity-shift synthesis](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/d19_parity_shift_synthesis_20261014.md):
  a special locus can have a radically different phenotype from a larger
  family. No claimed correspondence between the geometric obstructions.
- Primitive-triples search returned no usable old packet on this pass.
- [Rule 30-inspired local gate](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/chanpak_target_carry_20261115/XOR_CARRY_GATE.md):
  fixed guard bits make a chosen next output bit permutive. At 113,
  fixed depth phase plays the guard, and the exact Jordan kick makes
  each next output digit permutive. This is a design analogy, not a
  simulation of Rule 30.

## A finite-field recipe for the exact critical mass

Let p≥5, let a=min_{i,j}v_p((T^m−I)_{ij}), and put
D=(T^m−I)/p^a mod p. If d=3a, then det D is nonzero and
the row pair (e₁ᵀD,e₁ᵀD²) is independent over F_p. This follows
either from the determinant formula or by reducing (p^(−a)L) mod p
and evaluating det B / p^(3a). The exceptional fraction is then
the simple finite number
`γ(v,p) = #{0≤i<m : e₁ᵀD T^i v=0 mod p}/m.`

Every regular phase is Haar on one coset mod p^a. An exceptional
phase is either Haar on a finer coset or, up to a translation,
p^(2a) times a unit multiple of a Haar-uniform square on Z_p.
Indeed on a critical phase divide the
analytic difference F(t)−F(0) by p^(2a) **if** its linear coefficient
is divisible by p^(2a). This always holds for a=1 but need not hold
for a>1. If the linear coefficient instead has valuation b<2a
(with b≥a+1), its leading term dominates and that phase is Haar on
a coset mod p^b; its high-conductor contribution vanishes. When
the division is allowed, the reduced quadratic coefficient is
e₁ᵀD²T^i v/2, a unit, and every term of degree at least three is
divisible by p. The exact square form follows by solving F'=0
via digit lifting after normalization, and taking the convergent
square root of the unit Taylor factor as in the earlier quinary
proof. Thus for r>2a the universal theorem improves to
`|M_v(h/p^r)| ≤ γ(v,p) p^{−(r−2a)/2}`.
For a=1 no extra derivative test is needed and there is exact
cancellation at every r≥2 if the phase orbit avoids the critical
line. No assertion is made that d=3a holds for every p.
