# Independent review: marked strata and exact denominator laws

Verdict: PASS for the final mathematical proof, conditional on the explicitly cited inherited tournament inputs. The review does not reprove the effective strong-order bound, rigid alphabet closure, or strong h=3 classification. The new denominator and marked-series arguments are sound. The strengthened non-D-finiteness statement is also sound.

## Main findings

1. The marked SCC decomposition is type-correct. Strong ordinary atoms, fixed strong centers, and the single empty center remain different objects. Ordinary half-words are not converse orbits. Every additive, converse-invariant marker must be squared along with z in the half substitution; squaring z alone is incorrect. Direct products of ordered SCC automorphism groups justify prime-valuation additivity. The h and b exponents in E_(x,d) are x²a and x beta, not the center values alone.
2. The finite-polynomial assertion uses actual finite alphabets. At x=1, finiteness comes from the positive-deficit strong-order bound, not division by x(x−1). At x>1, x(x−1)a≤d bounds the center h. The resulting finite union of alphabets indeed gives a finite set of possible Aut/H primes on each fixed-degree stratum, apart from the already included unbounded rigid prime 3. This does not compute that prime set or any unknown coefficients.
3. The degree-zero marked formula is correct: rigid halves have x=1, while the center is empty, K1, C3, or R5. The empty tournament has h=b=|Aut|=H=1 and contributes y t, not one, when h and b are marked. Its raw automorphism marker is t as well.
4. Both irreducibility arguments are valid. In particular, q=1−z²−z⁶ is NOT irreducible modulo 2: q mod 2=(1+z+z³)². Its Q-irreducibility follows from the cubic-extension argument and the contradiction 2b=a², b²−2ac=1, c²=1. Also q mod p=2z(1−z), giving coprimality.
5. For F_j, positivity at rho makes the pole order exactly r(j)+1, then irreducibility propagates this exponent to the whole p factor. The corresponding positive sum at sqrt(rho) proves the exact q exponent of every nonempty S_d and J_j. Empty strata correctly have denominator one. These arguments never assume every odd h is realized.
6. For B_(3^r), the p exponent is r+1 because the unique length-r factorization uses only R5. At q, the half bound is t≤floor(r/2) except r=1, where t=1 is allowed. The leading J numerator is z^(10m)(1+z+z³) for even r and z^(10m+5) for odd r. If q divided the proposed cancellation numerator over Q, primitive polynomial divisibility would imply divisibility after reduction mod 2. That is impossible: the leading numerator contains at most one copy of 1+z+z³, whereas q mod 2 contains two. This covers r=0 without relying on an unknown positive center coefficient. For r=1, the q² term from the x=3 rigid-center family survives directly.
7. The exact reduced denominator therefore gives minimal eventual homogeneous recurrence order 3(r+1)+6e(r), with e(1)=2 and e(r)=floor(r/2)+1 otherwise. Unknown strong atoms and centers affect numerator, initial conditions, and onset, not this order. General odd-degree q cancellation is correctly left open.
8. The raw-aut series is correctly interpreted in Q(z)[[t]]. Its only nonzero coefficient indices are 3^c, with coefficients z^(3c)/(1−z)^(c+1). Any D-finite relation gives a finite-width polynomial-coefficient recurrence. Its last nonzero coefficient polynomial is invertible at all sufficiently large integer indices; a sufficiently late zero block therefore forces an identically zero tail. This contradicts the infinitely late 3^c coefficients. The width-zero case is explicitly covered. Thus the series is not D-finite in t over Q(z), and in particular not rational. No unsupported assertion about differential finiteness in z at fixed t is made.

## Inheritance and novelty boundary

The prior fixed-degree packet already contains the SCC factorization and positive-pole argument. The symmetry-width packet, Section 3, already contains the h=1 triangle-marked rational series, binomial coefficient formula, and C3^c automorphism law. The earlier leaf and Session22 trivalent packets already establish exact recurrence orders nine and eighteen. These are now explicitly credited as controls, not new findings.

The supported additions are the finite multivariate extension, exact general F/S/J denominators expressed through actual atom support, the all-r power-of-three exact denominator family beyond the inherited first two cases, and the raw-order non-D-finite obstruction. This is a repository-relative review, not a literature-priority certification.

## Independent executable controls

The standalone checker in this directory shares no code with the author's checker and uses no assert statements. It uses an eleven-letter involutive formal alphabet. Only K1, C3, R5, and S4 have certified tournament signatures; the other letters deliberately test repeated signatures, distinct converse mates, and higher-weight unknown-atom behavior. They are not claimed to exist as tournaments.

The checker verifies:

- Literal word enumeration through vertex order 16 with no h-pruning, compared against independently computed weighted convolution formulas.
- SCC, nu_3(|Aut|), nu_5(|Aut|), h, and b marker evaluations at 2,3,5,7,11 respectively.
- F and J at h≤200; S at even d≤200; B at k≤100; the separate actual marked B0 boundary.
- Exact rational F denominators through h=150 and S denominators through even d=100.
- Exact B_(3^r) denominators for r=0,...,8, producing recurrence orders 9,18,21,24,33,36,45,48,57.
- Mod-2 noncancellation and the half-exponent feasibility bound for r=0,...,100.
- The raw-aut empty term and coefficient-degree growth through order 30.

These are controls for the exact algebra and conventions, not a census of large-order tournaments. The all-order proof rests on the mathematical arguments above, not finite test extrapolation.

No repository, Git, remote, or publication state was modified by this review. See normal.json, optimized.json, and receipt.json for the frozen proof hash and final check counts.
