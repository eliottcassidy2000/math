# Log-scale sampling: exact finite poles, aliasing, and the infinite-limit boundary

## Conventions and classical input

Let b>1 be an integer with at least two distinct prime factors, such as10 or6. Then b^m is never a prime power for m>=1, eliminating endpoint half-jump corrections. Put

a_(m+1) = [psi(b^(m+1))-psi(b^m)]/[(b-1)b^m], m>=1.

The classical explicit formula is

psi(x)=x-sum_rho^* x^rho/rho-log(2pi)-(1/2)log(1-x^-2),

at non-prime-power x>1. The star is the symmetric height-truncation limit, with multiplicities. It is not an absolutely convergent unordered sum. Classical background: Ian Petrow's Introduction to Analytic Number Theory, explicit-formula chapter, https://metaphor.ethz.ch/x/2018/hs/401-3100-68L/IANTnotes2018.pdf ; the sharp psi formula is also stated in equation1.12 of the Manchester doctoral text https://pure.manchester.ac.uk/ws/portalfiles/portal/274923527/FULL_TEXT.PDF . These inputs are standard theorems, not new results of this note.

Direct subtraction gives the COEFFICIENTWISE identity

a_(m+1) = 1 + sum_rho^* c_(rho,b) lambda_(rho,b)^m + t_(m,b),

lambda_(rho,b)=b^(rho-1),
c_(rho,b)=-(b^rho-1)/[(b-1)rho],
t_(m,b)=-log[(1-b^(-2m-2))/(1-b^(-2m))]/[2(b-1)b^m].

Here t_(m,b)<0. The constant log(2pi) cancels. This identity is obtained by subtracting two convergent symmetric sums, with the same cutoff. It requires no interchange of a zero sum and a generating series.

## 1. Finite-zero generating functions are rational and their poles survive

Let S be any finite multiset of nontrivial zeros, including multiplicity. Its zero contribution has ordinary generating function

Z_(b,S)(z) = sum_(rho in S) c_(rho,b) lambda_(rho,b) z/[1-lambda_(rho,b) z].

This equality follows from finitely many geometric series near0, then rational continuation. Every pole has location

z_(rho,b)=b^(1-rho), so |z_(rho,b)|=b^(1-Re rho).

Aliased zeros have the same lambda. Necessarily their real parts coincide and their imaginary parts differ by integer multiples of2pi/log b. In one alias class, b^rho is constant. Its combined coefficient is

C=-(b^rho-1)/(b-1) * sum_(rho' in class)1/rho'.

Because0<Re rho'<1, every Re(1/rho') is strictly positive. Also |b^rho|>1, so b^rho-1 is nonzero. Therefore C is nonzero. Consequently EVERY distinct finite alias class contributes a genuine SIMPLE pole. Multiplicity changes its residue, not its order. This excludes finite pole cancellation for this particular arithmetic coefficient system; cancellation would be possible for arbitrary signed weights.

Conditionally on RH for the zeros in S, all these poles have modulus sqrt(b). Conversely, if every distinct pole of this finite model has that modulus, then every zero used in that model has real part1/2. This is only a finite-model equivalence.

## 2. Two bases recover a matched zero, but do not perform the matching

For a single zero, one observed pole in base b determines

Re rho = 1-log|z_b|/log b,
Im rho = [-arg z_b+2pi j]/log b, j in Z.

If b,c are multiplicatively independent positive integers, log b/log c is irrational. Suppose the SAME zero has been matched to z_b and z_c. If rho and rho' yield the same two poles, then

rho-rho'=2pi i j/log b = 2pi i k/log c.

Irrationality forces j=k=0. Thus the pair determines rho uniquely. Bases6 and10 have no prime-power endpoints and are multiplicatively independent by unique prime factorization.

Important boundaries: this does not identify which pole from one unlabeled spectrum matches which pole from the other. Several aliased zeros at one base may have merged residues. Approximate pole locations also do not provide a numerically stable recovery theorem without a height bound and separation estimates. Irrational rotations can make distinct large heights arbitrarily close in both phases; uniqueness with exact infinite precision is not uniform stability.

## 3. The elementary correction has a separate, fully controlled pole tower

Expanding the logarithm absolutely, for m>=1,

t_(m,b)=sum_(j>=1) d_(j,b) mu_(j,b)^m,

d_(j,b)=(b^(-2j)-1)/[2j(b-1)], mu_(j,b)=b^(-(2j+1)).

Therefore T_b(z)=sum_(m>=1)t_(m,b)z^m equals

sum_(j>=1) d_(j,b) mu_(j,b) z/[1-mu_(j,b)z].

The double series is absolutely convergent for |z|<b^3. The right side then extends meromorphically to C with simple poles at b^(2j+1), j>=1. Proof: on any compact set avoiding these poles, the tail is bounded by a constant times sum_j b^(-2j)/j. Every residue is nonzero and the pole locations are distinct. Thus trivial-zero corrections and nontrivial-zero finite modes are rigorously distinguishable by pole locations.

## 4. What follows for the ACTUAL infinite arithmetic generating function

Define A_b(z)=sum_(m>=1)a_(m+1)z^m. The prime number theorem gives a_(m+1)->1. Hence the ordinary series has radius exactly1 and

(1-r)A_b(r)->1 as r increases to1.

This is an Abelian leading behavior, not by itself a proof that A_b has a meromorphic simple pole at1. For the centered series E_b(z)=A_b(z)-z/(1-z), the elementary estimate a_(m+1)=O(m) already gives analyticity in |z|<1. Under the standard RH consequence psi(x)-x=O(sqrt(x)(log x)^2), its coefficients satisfy a_(m+1)-1=O(m² b^(-m/2)), so E_b is analytic in |z|<sqrt(b). No boundary-pole or natural-boundary conclusion follows from that lower bound on its radius.

### Unconditional local uniform identity inside the unit disk

We can close the interchange obligation here using a quantitative classical theorem, rather than coefficientwise convergence alone. Montgomery–Vaughan, *Multiplicative Number Theory I*, Theorem12.5, p400, gives for fixed c>1, x>=c and T>=2:

R(x,T) = O_c(log(x) min(1,x/[T dist(x,prime powers except x)]) + x log²(xT)/T),

where R is the remainder after the sharp explicit formula's sum over all nontrivial zeros with |Im rho|<=T. This is unconditional; the theorem statement requires no special choice of cutoff height. The proof may choose an auxiliary height internally, but that is not a hypothesis of the result. Primary text: https://www.ndl.ethernet.edu.et/bitstream/123456789/23715/1/Hugh%20L.%20Montgomery.pdf , Theorem12.5, equations12.3–12.4. The full PDF fetch timed out, but the indexed primary text exposed the theorem, assumptions, remainder and beginning of its proof. The similar display in Mazur–Stein's Appendix A is written under their simplifying RH convention, so it is NOT our source for unconditionality.

For x=b^m, its distance to every prime power is at least1. Since min(1,x/[T distance])<=x/T, we have

|R(b^m,T)|/b^m <= C_b [m log b+(m log b+log T)^2]/T.

Subtracting the two endpoints and dividing by (b-1)b^m gives, after absorbing the linear term,

|e_(m,T)| <= C'_b [((m+1)log b+log T)^2]/T,

e_(m,T)=a_(m+1)-1-t_(m,b)-sum_(|Im rho|<=T)c_(rho,b)lambda_(rho,b)^m.

This holds for every m>=1 simultaneously. No moving m cutoff is needed. Put L=log b and fix0<r<1. Elementary geometric sums yield

sum_(m>=1)r^m |e_(m,T)| <= (C'_b/T) [L² r(1+r)/(1-r)^3 + 2L(L+log T)r/(1-r)^2 + (L+log T)^2 r/(1-r)].

Thus this error is O_(b,r)((1+log² T)/T), tending to0. For each finite T, summing its finitely many zero modes as geometric series is legitimate in |z|<1 because |lambda|<1. Consequently

E_b(z)=lim_(T->infinity) Z_(b,S_T)(z)+T_b(z),

locally uniformly in the UNIT DISK, unconditionally. This is the desired analytic identity, now justified by an explicit summable bound rather than by a formal interchange. The new step is the direct geometric summation of the classical remainder at separated integer endpoints.

### What remains outside that disk

The quantitative estimate above diverges as r approaches1, and does not apply to r>1. It gives no meromorphic continuation across the unit circle, no pole persistence of the rational approximants, and no natural boundary. Under RH the independent coefficient estimate still proves analyticity of E_b in |z|<sqrt(b), but extending the ZERO-SUM convergence into that larger disk requires a further theorem.

Poles of approximants need not survive a limit: z^N/(1-z) has a pole at1 for every N but tends locally uniformly to0 in |z|<1, whose analytic continuation is entire. The finite alias noncancellation proof is therefore not an infinite-pole persistence theorem. Accumulation of pole locations, conditional summation, merged infinite alias classes and continuation across their closure remain distinct issues. No global RH equivalence, natural boundary, or recovery of all zeros from one sampled sequence is asserted.

## 5. Useful research target

The rigorous map now includes an unconditional locally uniform unit-disk identity. Zero real parts control finite-model radial pole locations, heights control phase, and logarithmic sampling causes aliasing. The remaining substantial target is continuation or pole persistence beyond the unit disk with all sums justified. A two-base method additionally needs a matching theorem and quantitative separation, not merely irrationality of logarithms. This note proves no such matching or persistence theorem.
