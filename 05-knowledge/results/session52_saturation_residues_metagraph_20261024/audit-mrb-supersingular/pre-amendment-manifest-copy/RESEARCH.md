# MRB, Ogg primes, and the mod-15 filter: a genuine CM bridge and its boundary

24 October 2026. **PROVED** elementary identities, periodic-filter theorem, and finite residue obstructions. **CITED** Ogg's classification, CM reduction theory, elliptic-function conventions, and the honeycomb connective constant. **FINITE-EXACT** integer certificates in `results.json`. **VERIFIED-INTERVAL** MRB enclosures with an analytic remainder. **OPEN** a modular formula for the actual MRB defect. No Fortune, LRC, or tournament-counting conjecture is solved.

## 1. Inheritance and the precise new connection

The inherited mechanisms are:

- Session31, `mrb/MRB_CONTOUR_AUDIT.md` and `analytic-audit/AUDIT.md`: the exact branch jump, the CM/honeycomb identities of A(8), a certified MRB defect, and an entire-perturbation obstruction. These files survive in the archived Session31 publication bundle, although the later checkout's matching route is absent. They are credited rather than re-presented as discoveries.
- Session32 and the audited note `monster_levels_square_sum_7_21_20261022.md`: the true Monster/Baby Monster/Fi22 level filtration 1 | 2 | 6. MISTAKE-668 prevents confusing an order quotient with a subgroup index or automatic nesting with a new correspondence.
- `supersingular_divisor_shapes_qseries_20261022.md`, repaired by MISTAKE-666: total versus prime-field supersingular counts, and the importance of preserving the exact level and CM discriminant.
- Session51, `prime-offsets/THEOREM.md`: the two signs of a primorial offset have the same small-prime wheel, not independent screens; exact prime supports and first-success minimality require extra data.

The canonical hostile is a residue-surviving object without the target witness. The corrected near miss is D approximately q^4, not D=q^4. The least-used sidecar here is the actual Frobenius trace of a concrete CM curve, rather than its shared vocabulary with Ogg's primes.

The new synthesis supplies a concrete, nonnumerological route:

    A(8) -> k=sqrt(2)-1 -> lambda=3-2sqrt(2) -> j=8000
         -> E: y^2=x^3+4x^2+2x -> supersingular reduction modulo p.

It preserves a CM elliptic curve and its reduction predicate, not the MRB alternating functional. Its prime filter is modulo **8**, and its quantifier is **one specified curve**. Ogg's predicate concerns **all** supersingular j-invariants in a characteristic. Sections 2–3 prove these distinctions and exhibit counterexamples in both directions.

Live concept board: MRB alternating functional; a single branch-jump amplitude; CM discriminant -8; absolute Ogg rationality; mod-15 unit screen; two-sign primorial support. No cosmetic tournament is imposed on these objects.

## 2. The residue A(8) really reaches a CM curve

Put

    f(n)=n^(1/n)-1,
    C=sum_{n>=1} (-1)^n f(n),
    A(m)=(2/pi) m^(-1/m) sin(pi/m), m>=4 in the residue dictionary tested here,
    D=A(8)-C.

The definition of C is the convergent MRB series, not the oscillating unregularized sum of (-1)^n n^(1/n). It appears in [Crandall, Section 7.5](https://www.marvinrayburns.com/UniversalTOC25.pdf) and [OEIS A037077](https://oeis.org/A037077). Since f decreases to zero after n=3, convergence follows by the alternating test. It is not absolute: f(n) is asymptotic to log(n)/n.

**PROVED algebra; CITED special-function conventions.** Let k=sqrt(2)-1 and k'=sqrt(1-k^2)=sqrt(2k). The AGM step and homogeneity give

    sqrt(2) AGM(1,k)=AGM(1,k').

Using K(k)=pi/[2 AGM(1,k')], one obtains K(k')/K(k)=sqrt(2). Thus the elliptic nome is q=exp(-pi sqrt(2)) and lambda(i sqrt(2))=k^2=3-2sqrt(2). Here k is the elliptic modulus, whereas lambda=k^2 is the parameter. See [DLMF 19.8.5](https://dlmf.nist.gov/19.8#E5) and [DLMF 23.15](https://dlmf.nist.gov/23.15).

Elementary trigonometry gives

    k=2^(1/4) (pi A(8))^2,
    pi 8^(1/8) A(8)=2 cos(pi/8) k.

The identification 2 cos(pi/8)=mu_hex is the [Duminil-Copin–Smirnov honeycomb theorem](https://arxiv.org/abs/1007.0575), not a consequence about C.

Now lambda^2-6lambda+1=0. Therefore

    1-lambda+lambda^2=5lambda,
    (1-lambda)^2=4lambda,
    j(lambda)=256(1-lambda+lambda^2)^3/[lambda^2(1-lambda)^2]=8000.

The lattice Z+Z i sqrt(2) is preserved by multiplication by i sqrt(2), so its elliptic curve has CM by Z[sqrt(-2)], the maximal order of discriminant -8. The displayed rational curve E has

    c4=160, Delta=512, j=c4^3/Delta=8000.

Over the complex numbers equality of j identifies its geometric isomorphism class with that lattice curve. The exact lambda and Weierstrass calculations are independently checked using integer pairs in Q(sqrt(2)), with no floating-point equality test.

By Deuring's reduction theorem, at every prime p>=5 the good reduction of E is supersingular exactly when p is inert in Q(sqrt(-2)). Since Delta=2^9, every odd prime is good, and the quadratic residue law gives

    E mod p supersingular <=> (-2/p)=-1 <=> p=5 or 7 modulo 8, p>=5.

For the theorem and its hypotheses see [Mula–Murru–Pintore, Theorem 4.1](https://link.springer.com/article/10.1007/s10231-024-01528-x); [Elkies–Ono–Yang, Section 3](https://people.math.wisc.edu/~tonghaiyang/deuring.pdf) gives the reduction mechanism. The supplementary p=3 calculation is finite and ordinary, not an extension inferred from a citation stated for p>=5.

**FINITE-EXACT.** For all 167 odd primes below 1000, the script counts E(F_p) by a Legendre-symbol sum and verifies trace zero exactly in the stated residue classes. For odd p<=107 it independently counts square residues with multiplicity. In particular:

| p | j mod p | #E(F_p) | trace | E supersingular? | Ogg prime? |
|---:|---:|---:|---:|---|---|
|17|10|24|-6|no|yes|
|37|8|38|0|yes|no|
|47|10|48|0|yes|yes|

Thus neither Ogg membership implies supersingularity of E nor supersingularity of E implies Ogg membership. Also 17 and 47 have the same residue modulo 15 and opposite E-reduction types, even though both are Ogg primes.

## 3. What modulo 15 can and cannot detect

**CITED.** Ogg's set is

    O={2,3,5,7,11,13,17,19,23,29,31,41,47,59,71}.

For a prime p, membership means all supersingular j-invariants in characteristic p lie in F_p, equivalently genus(X_0^+(p))=0. It also gives exactly the prime divisors of the Monster's order. These are the relevant statements in [Aricheta, Introduction](https://arxiv.org/html/1809.07421v2). Counting fifteen members does not make 15 a congruence mechanism.

**PROVED finite obstruction.** The following covers every unit residue modulo 15. The second column retains the small prime 2 where it occurs; deleting primes <=5 still leaves every unit residue represented.

| residue | Ogg primes | non-Ogg prime | E ordinary prime | E supersingular prime |
|---:|---|---:|---:|---:|
|1|31|61|211|31|
|2|2,17,47|107|17|47|
|4|19|79|19|79|
|7|7|37|67|7|
|8|23|53|83|23|
|11|11,41,71|101|11|71|
|13|13|43|43|13|
|14|29,59|89|59|29|

All listed numbers are certified prime by trial division. The Ogg/non-Ogg decisions are independently checked by reduced-form class numbers, and through 107 also by counting one nonsingular elliptic curve for every j in F_p. Specifically, for p>3,

    s(p)=floor(p/12)+epsilon_p,
    epsilon_p=0,1,1,2 for p=1,5,7,11 modulo 12,
    r(p)=[h(-4p)+1_(p=3 mod 4) h(-p)]/2,
    g_plus(p)=[s(p)-r(p)]/2.

These counting formulas are cited rather than newly proved: see [Delfs–Galbraith, Section 1](https://arxiv.org/pdf/1310.7789) for the total count and the equivalent piecewise rational-count formula. Here h(D) counts primitive positive-definite reduced forms, also for the nonfundamental discriminant -4p when p=3 modulo 4; it is not silently the class number of the maximal order in that case. The earlier exact implementation is recorded in Session51. This pass checks all primes 5<=p<=211, and the independent all-j point-count route covers p<=107. It does not extrapolate a finite census into Ogg's classification.

Consequences:

1. Any mod-15 residue filter that retains every Ogg prime >5 must retain all primes >5. Thus there is no nontrivial sound residue rejection at this modulus.
2. There is no function b:Z/15Z->{0,1} with b(p mod 15)=1_(p in O) for all primes p. Every unit residue has an explicit positive and negative witness.
3. There is no mod-15 determination of E's CM reduction type. All eight classes have explicit ordinary and supersingular witnesses.

**PROVED CRT statement.** For each unit r modulo 15, its four odd lifts modulo 120 correspond bijectively to residues 1,3,5,7 modulo 8. Two are ordinary and two supersingular for E. The certificate records all 32 lifts and one prime in each, with the largest least-prime representative 409. This is a CRT decomposition plus a finite witness bank, not a prime-density claim or independence heuristic. Modulo 120 is the correct common carrier for these two residue observables; it is not a new characterization of Ogg's set.

There is also a genuine **level 15**, distinct from reduction modulo 15. Aricheta's Appendix A gives rationality primes {2,7} for 15+ (the full Atkin–Lehner quotient X_0(15)/<w_3,w_5>) and {2} for 15+15 (Fricke-only X_0(15)/<w_15>), always at primes not dividing the level. These curves and lists are not the projection of O onto Z/15Z. The distinction also avoids conflating the Fricke-only plus convention with the full quotient at composite level.

## 4. A real mod-15 MRB theorem: retain parity or lose convergence

This section derives a new synthesis rather than an Ogg/MRB identification.

### Theorem 1. Periodic means control root-profile divergence

Let a_n be a real periodic sequence of period L with mean mu. Then

    sum_{n<=N} a_n f(n) = (mu/2)(log N)^2 + O(1).

If mu=0, the partial sums converge. If mu!=0 they diverge with the displayed leading term.

**Proof.** Taylor's formula gives f(n)=log(n)/n+O((log n)^2/n^2); the error is summable. Hence sum_{n<=N} f(n)=(log N)^2/2+c+o(1). The sequence a_n-mu has bounded partial sums, and f(n) decreases to zero eventually. Summation by parts therefore makes sum(a_n-mu)f(n) convergent. Combining proves the statement. QED.

### Theorem 2. Odd-modulus filters have a necessary parity sidecar

If w is any real 15-periodic sequence, then

    C_w=sum_{n>=1} (-1)^n w(n) f(n)

converges. Its coefficient sequence has period 30 and zero mean, because translating by 15 reverses its sign. A smaller fundamental period is possible for special w; for the unit mask below the fundamental period is exactly 30.

There is a justified residue decomposition

    C_w=sum_{r=1}^{15} (-1)^r w(r)
          sum_{k>=0} (-1)^k f(15k+r)
       =sum_{r=1}^{15} (-1)^r w(r)
          sum_{j>=0} [f(30j+r)-f(30j+r+15)].

Each inner alternating sum converges. Each paired inner sum converges absolutely, because f'(x)=exp(log(x)/x)(1-log x)/x^2 and the pair differences are O(log(j)/j^2). Finite residue sums may therefore be taken after limits; this is not an arbitrary rearrangement of a conditional series.

By contrast, splitting into the 30 one-sign residue series separately gives divergent series. Pairing r with r+15 is essential; it does not authorize subtracting two separately assigned infinite values.

### Corollary 3. The unit mask has an exact inclusion-exclusion law

Let u(n)=1 if gcd(n,15)=1, and zero otherwise. For odd d define

    M(d)=sum_{k>=1} (-1)^k [(dk)^(1/(dk))-1].

Then, with ordinary convergent sums,

    C_u=M(1)-M(3)-M(5)+M(15),
    C_(1-u)=M(3)+M(5)-M(15).

**Proof.** The pointwise identity u(n)=1-1_(3|n)-1_(5|n)+1_(15|n), together with (-1)^(dk)=(-1)^k for odd d, gives the finite-N identity with floor(N/d) upper limits. Every component converges, so taking the limit is legitimate. QED.

The unit mask has mean 8/15; deleting parity gives the sharp divergence law

    sum_{n<=N, gcd(n,15)=1} f(n)=(4/15)(log N)^2+O(1).

The complementary unalternated mask has leading coefficient 7/30. These density constants do **not** multiply the conditional constant C. Recomputed analytic-plus-interval certificates show

    C_u=0.5375844239373953410247043689688958465874...,
    C_(1-u)=-0.3497247814753282207761864349146226165315...,
    C_u-(8/15)C > 0.437.

The positive unit contribution even exceeds C, because the complementary signed contribution is negative. Therefore a density-based split of the constant is directly refuted, not merely unsupported.

An exact generating function makes the carrier visible. For |z|<1,

    R(z)=sum_{n>=1} (-1)^n u(n) z^n
        =[-z+z^2+z^4-z^7+z^8-z^11-z^13+z^14]/(1+z^15)
        =-z/(1+z)+z^3/(1+z^3)+z^5/(1+z^5)-z^15/(1+z^15).

This is the rational generating function of the **mask**, not of the MRB coefficients f(n). Its coefficient identity is checked through degree 600. The all-degree proof is the geometric series on each residue class. The 15-point quotient forgets the sign at n and n+15; n=1 and 16 are the smallest explicit witness.

## 5. Primorial two-sign filters and the actual shared operation

For a primorial P divisible by 15, every prime divisor ell of P satisfies

    P+d=0 mod ell <=> d=0 mod ell <=> P-d=0 mod ell.

Thus the local supports for the two signs coincide. Modulo 15 the survivors are the eight units, with fraction 8/15, not (8/15)^2. For P=30030 the full wheel has 5760 survivors and fraction 192/1001. These are exact counts of residue classes, not probabilities that one or both endpoints are prime.

Inherited decisive hostiles remain:

    30030-19=30011 is prime; 30030+19=151*199,
    30030-23=37*811;       30030+23=41*733.

Both offsets pass the full wheel, but only one sign succeeds at 19 and neither succeeds at 23. Endpoint-equals-sieving-prime exceptions must be excluded when using this as a rejection argument; the displayed large endpoints have no such exception.

The legitimate common operation is a mask on integer indices: u(d) screens offsets, while u(n) projects MRB summands. The map does not preserve endpoint primality, Fortunate minimality, the MRB sum, or a native metagraph completion witness. As in the inherited conditional-hole analysis, a superset of feasible support is not an actual witness set, and OR-versus-AND changes the complement operation. This is a usable information-loss diagnosis, not a new prime or tournament theorem.

## 6. The MRB defect remains an open identification problem

The certified facts are

    D=A(8)-C=1.89180712788007352705453393844... * 10^-8,
    D/q^4=0.98845074389577451997768768858...,
    0<D<q^4.

Thus D=q^4 is false. The exact even-theta identity

    [theta_3(0,q)+theta_4(0,q)-2]/4
       =sum_{r>=1} q^(4r^2)>q^4

also rules out that unweighted positive tail as D. This says nothing against every possible signed or normalized theta identity; any proposed normalization must be specified without fitting it to D.

The certificate recomputes C, not merely its recorded digits. For M>=4 and N>=1 the finite Euler transform is

    sum_{j>=0} (-1)^j f(M+j)
      =sum_{j=0}^{N-1} (-1)^j Delta^j f(M)/2^(j+1)+R,

where Delta f(x)=f(x+1)-f(x). With L=log(3M/2)+pi/6 and b=2L/M,

    |R| <= 2 N! exp(b) [L/M^(N+1)+(L/N+1/N^2)/M^N].

To justify the bound, use Cauchy's estimate on |z-x|=x/2, where |arg z|<=pi/6 and |f(z)|<=2 exp(b)[log(3x/2)+pi/6]/x. Express the Nth finite difference as an integral of f^(N) over the N-dimensional unit cube, and bound the remaining positive decreasing sum by its first term plus its integral. The factor 2^-N from Euler's remainder gives the displayed result. At M=256,N=100 it is below 2.654e-84.

For g_d(x)=f(dx), replace L by log(3dM/2)+pi/6, b by 2L/(dM), and multiply the displayed bound by 1/d. This certifies all four M(d) used in Section 4. Finite summation uses exact binomial weights and 130-decimal directed interval arithmetic. These are analytic-plus-computational certificates, not a formal proof-kernel claim about the interval library.

Both nearest-amplitude indices 8 and 105719004 pass strict neighboring-midpoint tests. The latter is a ceiling of the continuous inverse, while the former is a floor. Consequently “greedy” here must mean nearest amplitude; a single fixed floor or ceiling convention does not produce both. Monotonicity for x>=4 follows from

    A'(x)/A(x)=[log x-1-pi cot(pi/x)]/x^2<0.

For completeness, cot t>1/t-t on (0,pi/4], obtained by differentiating t cos t-(1-t^2)sin t, gives pi cot(pi/x)>x-pi^2/x>log x-1. The last inequality holds at x=4 and its difference increases thereafter. So neighboring-midpoint tests are global, not merely local choices.

The inherited obstruction remains exact: for any K>=1,

    h_K(z)=e^(-z) product_{j=1}^K(z-j)

is entire, vanishes on the first K positive integers, and is smaller than every inverse power on the positive axis. Adding c h_K to f changes neither the branch jumps nor those finite samples nor the complete positive-axis inverse-power asymptotics, but changes its alternating constant by

    c (-1)^(K+1) K!/(e+1)^(K+1).

The formula follows by summing binomial(n-1,K)(-e^-1)^n. These perturbations are not nth-root profiles and need not obey stronger uniform half-plane decay. Hence the obstruction only rejects reconstruction from the stated incomplete data, not a special identity using the full MRB profile. The actual modular identification of D and the arithmetic nature of C remain OPEN.

## 7. Typed connection ledger and decisive tests

| source | target | map | preserved predicate | lost data / needed sidecar | decisive test |
|---|---|---|---|---|---|
|A(8)|CM elliptic curve|normalize to k, square to lambda, compute j|j=8000 and CM discriminant -8|C and every remaining summand are absent|exact lambda polynomial and j calculation PASS|
|CM curve E|supersingular reductions|reduce at a good prime, apply Deuring|one-curve supersingularity; chi_-8|all-j rationality is not preserved|p=17 and 37 refute both Ogg implications|
|Ogg primes|residues modulo 15|p maps to p mod 15|a residue label only|class numbers, Frobenius pairs, height|all eight classes have positive/negative witnesses|
|mod-15 and mod-8 labels|common residue carrier|CRT to mod 120|both residue observables|does not recover Ogg predicate or MRB sum|32 exact lifts, prime witnesses <=409|
|primorial offsets|unit mask|reduce d modulo 15|small-prime rejection, with endpoint exception excluded|large factors and first success|d=19,23 near 30030|
|MRB summands|unit-filtered functional|multiply by u(n), retain (-1)^n|exact inclusion-exclusion|forget parity and convergence is lost|period 30; 4/15 log^2 divergence; scalar-density claim fails|
|theta parity tail|actual MRB defect|no established map|only a common small exponential scale|coefficient derivation and remainder needed|positive theta tail>q^4>D|
|these arithmetic objects|tournament/metagraph counts|no target-preserving map supplied|none asserted|compatible full records, multiplicities, automorphisms|do not infer A000568 values or native completion|

## 8. Reproduction, scope, and remaining questions

From this directory run:

    python3 check.py > check.out
    cp results.json results.normal.json
    python3 -O check.py > check_optimized.out
    cmp results.json results.normal.json
    cmp check.out check_optimized.out
    python3 test_rejections.py > rejections.out
    python3 -O test_rejections.py > rejections_optimized.out

All verification predicates use explicit exceptions, not removable assertions. `test_rejections.py` mutates the Ogg list, CM sign rule, parity generator, and analytic error target and checks that normal and optimized executions reject each mutation. `MANIFEST.json` supplies SHA-256 hashes; `SOURCES.json` distinguishes retrieved primary sources from the archived proofs and current checkout. No publication is performed by this bundle.

Stopping boundary: the concrete CM bridge, exact congruence obstruction, periodic-filter theorem, and certificates answer this lane. To go beyond them, specify a full-profile identity transporting C or D to a modular object, or a native map transporting an arithmetic witness to a metagraph record. More coincidences of counts, densities, or first terms do not cross either boundary. No literature-wide novelty claim is made.
