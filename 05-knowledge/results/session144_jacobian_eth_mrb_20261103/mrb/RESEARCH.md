# MRB constants, hidden entire profiles, periodic resonance, and precision

**Status.** PROVED elementary analytic statements below; CITED classical gamma/zeta conventions and Crandall's prior MRB representation; FINITE-EXACT and VERIFIED-INTERVAL checks labelled in receipts; numerical cross-checks are not interval certificates. The modular identification and arithmetic nature of MRB remain OPEN. Nothing here supplies a Jacobian counterexample classification or a complexity-class collapse.

## 1. Inheritance before extension

Put

    f(z)=exp(Log(z)/z)-1,  Re z>0,
    C=sum_{n>=1}(-1)^n f(n),
    A(m)=(2/pi)m^(-1/m)sin(pi/m),
    D=A(8)-C,  q=exp(-pi sqrt(2)).

The sum defining C converges in its stated order. The unshifted terms (-1)^n n^(1/n) do not tend to zero, so their ordinary infinite series is not a definition of C.

Recovered source paths in snapshot `5e1cedf`:

1. `05-knowledge/results/session31_mirror_mrb_20261022/mrb/MRB_CONTOUR_AUDIT.md`: branch-jump formula; shifted Abel-Plana representation; algebraic odd-Laplace mode tail; entire-perturbation nonuniqueness; certified D; exact A(8)/CM/honeycomb identities.
2. `05-knowledge/results/session31_mirror_mrb_20261022/ENTIRE_MRB_DEFORMATION.md` and `analytic-audit/AUDIT.md`: entire parameter deformation, uniform remainder, independently certified nearby zero, and certified Euler computation.
3. `05-knowledge/results/session31_mirror_mrb_20261022/GREEDY_INVERSE_CLOCKS.md`: universal nearest-amplitude expansion and convergence, not an MRB-specific modular theorem.
4. `05-knowledge/results/session52_saturation_residues_metagraph_20261024/mrb-supersingular/RESEARCH.md`: the concrete A(8)->CM discriminant -8, j=8000 route; parity-preserving mod-15 filters and their fourteen omitted Fourier modes.

These already establish

    D=1.8918071278800735270545339384419831105875... * 10^-8,
    D/q^4=0.98845074389577451997768768858315...,
    0<D<q^4.

Thus D=q^4 is false, as is equality with the unweighted positive even-theta tail sum_{r>=1}q^(4r^2). Identifying one residue with CM and honeycomb data does not transport the full alternating functional. This packet does not rediscover those inherited results.

The prior obstruction used exp(-z)P_K(z), where P_K(z)=product_{j=1}^K(z-j). It correctly noted that this perturbation did not preserve a uniform O(log|z|/|z|) bound across the entire right half-plane. Section 2 closes precisely that leftover possibility. The old theorem remains true and did not assert sufficiency of growth alone. The stronger growth condition that it explicitly left outside its example is now preserved by the new example as well.

**Concept board.** Full MRB profile and alternating functional; branch jump; finite prefix; uniform half-plane asymptotics; periodic Fourier mean/resonance; cutoff versus Abel regularization; certified numerical precision. Each connection below states the map and the information lost rather than transferring a conclusion on vocabulary alone.

## 2. Stronger nonuniqueness: even uniform half-plane asymptotics are blind

### Theorem 1 (smooth-Laplace invisible perturbations)

Fix K>=1 and a nonzero nonnegative C-infinity real function phi supported in (1,2). Define

    P_K(z)=product_{j=1}^K(z-j),
    H_K(z)=P_K(z) integral_1^2 phi(t)exp(-tz) dt.

Then:

1. H_K is entire and H_K(j)=0 for j=1,...,K.
2. H_K has zero discontinuity across the negative real axis. Adding cH_K preserves every branch jump of f.
3. For every integer R>=0 there is B_{K,R}<infinity such that

       |H_K(z)| <= B_{K,R}(1+|z|)^(-R),  Re z>=0.

   In particular all algebraic asymptotic orders, uniformly in this half-plane, remain unchanged. Along the positive real axis the perturbation is exponentially small as well.
4. Its alternating functional is explicitly nonzero:

       L(H_K)=sum_{n>=1}(-1)^n H_K(n)
             =(-1)^(K+1)K! integral_1^2 phi(t)/(e^t+1)^(K+1) dt.

Consequently for each prescribed real delta there exists a real c such that f+cH_K has every datum in (1)-(3) unchanged and alternating constant C+delta.

**Proof.** Integrating an entire exponential over compact t-support is entire in z, by uniform convergence on compact z-sets. The finite zeros come from P_K, and entire functions have no branch jump.

For z!=0, integration by parts m times has no boundary terms because phi and all derivatives vanish near the endpoints. It yields

    integral phi(t)e^(-tz)dt = z^(-m) integral phi^(m)(t)e^(-tz)dt.

On Re z>=0 the right side has magnitude at most |z|^(-m)||phi^(m)||_1. For |z|>=1, |P_K(z)|<=product_{j=1}^K(1+j)|z|^K. Choose m>=K+R. The compact remaining half-disk is bounded by continuity. This proves (3). Positive-axis exponential decay follows directly from t>=1 and polynomial growth of P_K.

The summands decay exponentially on positive integers, so absolute summation/integration interchange is valid. For x=-e^(-t), |x|<1, the binomial generating identity gives

    sum_{n>=1} P_K(n)x^n
       = K! sum_{n>=K+1} binom(n-1,K)x^n
       = K! x^(K+1)/(1-x)^(K+1)
       = (-1)^(K+1)K!/(e^t+1)^(K+1).

Integration proves (4). Its integral is strictly positive. Choose c=delta/L(H_K). QED.

**Scope and failure boundary.** This is nonuniqueness inside an explicitly enlarged analytic class, not a claim that f+cH_K is still the nth-root profile exp(Log z/z)-1. Preserving the actual global defining equation would prohibit these perturbations, but would also supply far more data than branch jumps or asymptotics. The perturbations can change growth in the left half-plane; no matching left-half-plane growth condition is claimed. No contour closure is justified across zero merely from this result.

The theorem needs no Abel-Plana invocation: the alternating shift is absolutely convergent and computed directly. It nevertheless remains admissible on the shifted right-half-plane contour because it has stronger vertical decay than required.

### Explicit elementary hostile family

For each integer r>=1, the fully explicit function

    H_{K,r}(z)=e^(-z)P_K(z)[(1-e^(-z))/z]^(K+r),

with the removable quotient at zero equal to 1, is entire and satisfies all finite-prefix/branch/asymptotic conditions except that its uniform half-plane decay is O_K,r(|z|^(-r)), not all orders at once. Along positive x it is exponentially flat.

Since (1-e^(-z))/z=integral_0^1 e^(-tz)dt, its alternating shift equals

    (-1)^(K+1)K! integral_[0,1]^(K+r)
           [e^(1+t_1+...+t_(K+r))+1]^(-(K+1)) dt.

The integrand is strictly positive and bounded. In particular the absolute shift has the explicit bounds

    K!/(e^(K+r+1)+1)^(K+1)
       <= |L(H_{K,r})| <= K!/(e+1)^(K+1).

No high-dimensional quadrature is needed for the proof. The code independently compares direct summation with the cube integral reduced to one dimension by the Irwin-Hall density. The r=1 family already preserves the formerly suggested uniform O(log|z|/|z|) condition.

## 3. Complete classification of periodic filters

Let a=(a_n) be a complex L-periodic sequence and mu=L^(-1)sum_{r=1}^L a_r. Write b_n=a_n-mu.

### Theorem 2 (mean is the exact convergence obstruction)

There is a finite constant K_a such that

    sum_{n<=N}a_n f(n) = (mu/2)(log N)^2+K_a+o(1).

The ordinary series converges if and only if mu=0. Its zero-mean convergence space has complex dimension L-1. In particular this is a complete classification of convergence for the specified L-periodic filters, not a classification of their closed forms.

**Proof.** Put h_n=log(n)/n. Taylor's remainder gives

    0<=f(n)-h_n<=exp(H)h_n^2/2,
    H=log(3)/3,

and sum h_n^2<infinity. The Stieltjes convention is

    gamma_1=lim_N[sum_{n<=N}log(n)/n-(log N)^2/2].

Thus with R=sum_{n>=1}(f(n)-h_n) and kappa=gamma_1+R,

    sum_{n<=N}f(n)=(log N)^2/2+kappa+o(1).

The periodic zero-mean sequence b has bounded partial sums; f(n) decreases to zero for n>=3. Summation by parts gives convergence of sum b_n f(n). Therefore

    K_a=mu*kappa+sum_{n>=1}b_n f(n).

If mu!=0 the nonzero log-squared coefficient prevents a finite limit, also for complex mu. QED.

**All-modulus parity classification.** For an m-periodic weight w, consider a_n=(-1)^n w(n).

- If m is odd, a has period dividing 2m and zero mean for every w: all such filtered MRB series converge.
- If m is even, a has period dividing m and

      mu=[sum_{r even}w(r)-sum_{r odd}w(r)]/m.

  Convergence is equivalent to equal even/odd weight totals, one explicit codimension-one condition.
- For the unit mask w(n)=1_{gcd(n,m)=1}, every odd m converges, but every even m has mu=-phi(m)/m and hence divergence

      sum_{n<=N, gcd(n,m)=1}(-1)^n f(n)
         =-phi(m)/(2m)*(log N)^2+K_a+o(1).

The even-modulus hostile is already m=2: selecting odd integers destroys all alternating cancellation. Session52's mod-15 mechanism extends to all odd moduli but cannot extend unchanged to even ones.

### Fourier boundary form

For omega=exp(2pi i/L), use

    a_hat(j)=L^(-1)sum_{r=1}^L a_r omega^(-jr),
    G(z)=sum_{n>=1}f(n)z^n, |z|<1.

Then

    sum_{n>=1}a_n f(n)z^n=sum_{j=0}^{L-1}a_hat(j)G(omega^j z).

It is a coefficient identity of an absolutely convergent series in the disk. Every G(xi) for |xi|=1, xi!=1, is an ordinary convergent series by Dirichlet's test. The coefficient a_hat(0)=mu is therefore precisely the resonant mode at z=1.

More generally at boundary z=omega^(-j), divergence occurs exactly when a_hat(j)!=0; subtracting its log-squared radial singularity leaves a finite limit as described below. At all other circle points every mode is nonresonant and converges. This classifies the boundary convergence set of each periodic filter. It does not assert holomorphic continuation through those boundary points or a modular transformation for G.

## 4. Regularization is a specified operation, not a free subtraction

### Theorem 3 (cutoff/Abel finite-part conversion)

For real t>0 define

    A_a(t)=sum_{n>=1}a_n f(n)e^(-tn),
    ell=log(1/t),
    Q=gamma^2/2+pi^2/12.

Then

    A_a(t)=mu*(ell^2/2-gamma*ell)+K_a+mu*Q+o(1), t->0+.

Thus the cutoff finite part K_a and this specified Abel finite part differ by mu*Q. They agree precisely for zero-mean filters (or trivially when mu*Q=0). If mu!=0, subtracting only mu*ell^2/2 from the Abel series leaves the divergent term -mu*gamma*ell.

**Direct proof.** Let F_a(x)=sum_{n<=x}a_n f(n), extended as zero for 0<=x<1. Stieltjes summation gives

    A_a(t)=t integral_0^infinity F_a(x)e^(-tx)dx
          =integral_0^infinity F_a(u/t)e^(-u)du.

For x>=1, write F_a(x)=mu log^2(x)/2+K_a+r(x), where r(x)->0 and is bounded. Integrating over u>=t, the bounded remainder tends to zero by dominated convergence. The missing interval 0<u<t contributes o(1) to each polynomial in ell and log u because t*(1+ell^2)->0 and the logarithmic moments near zero are integrable. Finally

    integral_0^infinity e^(-u)du=1,
    integral_0^infinity e^(-u)log(u)du=-gamma,
    integral_0^infinity e^(-u)log^2(u)du=gamma^2+pi^2/6.

Expanding (ell+log u)^2 proves the formula. The two logarithmic moments follow by differentiating Euler's gamma integral at 1, or from the standard gamma expansion. QED.

A second independent derivation uses the classical polylog expansion [DLMF 25.12.12](https://dlmf.nist.gov/25.12#E12): differentiate at s=1 after cancelling the gamma/zeta poles. It gives

    sum_{n>=1}log(n)/n * e^(-tn)
      =ell^2/2-gamma*ell+gamma_1+Q+o(1).

The absolutely summable nonlinear remainder f(n)-h_n tends to R. This independently checks the sign of gamma and the finite-part conversion. The Stieltjes convention matters; [DLMF 25.2.4-5](https://dlmf.nist.gov/25.2#E4) fixes it.

The operation “discard the mean” gives the canonical convergent functional sum(a_n-mu)f(n). It is not the same as declaring every divergent subseries to have an intrinsic value. A finite-part prescription has to be named, and a change of cutoff scale or damping kernel can change its value.

## 5. Numerical precision is an input size, not an ETH theorem

The inherited Euler certificate uses M>=4,N>=1,

    sum_{j>=0}(-1)^j f(M+j)
       =sum_{j=0}^{N-1}(-1)^j Delta^j f(M)/2^(j+1)+R_M,N,

with Delta f(x)=f(x+1)-f(x). For

    L=log(3M/2)+pi/6, b=2L/M,

it proves

    |R_M,N| <= 2N!e^b [L/M^(N+1)+(L/N+1/N^2)/M^N].

The prior proof uses Cauchy circles |z-x|=x/2, then the unit-cube formula for finite differences and an integral comparison. It does not assume complete monotonicity of f. The old M=256,N=100 certificate and its exact claim are replayed; the independent `audit/` establishes an explicit precision rule and accounts for the bit model and stability.

A useful exact conditioning identity is already elementary. Flatten the finite Euler sum into

    sum_{k=0}^{N-1}(-1)^k v_k f(M+k),
    v_k=sum_{j=k}^{N-1}binom(j,k)/2^(j+1)>0.

Then

    sum_{k=0}^{N-1}v_k
      =sum_{j=0}^{N-1}2^j/2^(j+1)=N/2.

Thus if every required f-value has absolute error <=epsilon, the whole prefix-plus-tail linear form has error at most (M-1+N/2)*epsilon. There is no exponential amplification after combining the weights exactly. A naive floating-point difference triangle can conceal this stability and is not the claimed bit algorithm.

**Exact scope:** the fixed real number C can be approximated with certified absolute error 2^(-p) using polynomially many bit operations in the requested output precision p, with elementary log/exp implemented by polynomial-time rational algorithms as detailed in the audit. This is uniform in p, not a theorem about arbitrary filters or input formulas. The code's directed-interval implementation is a concrete certificate subject to the correctness of its arithmetic library, not a formal proof kernel. Distinguish this elementary computable-analysis fact from the inherited Crandall series and from the unrelated complexity of SAT.

Crandall's 2012 paper [Section 7.5](https://www.marvinrayburns.com/UniversalTOC25.pdf) already gives the diagonal eta-derivative MRB series and discusses rapid computation. Hence fast computation of C is not a new surprise or a literature-priority claim. The contribution here is an explicit audited precision/stability contract for this repository's Euler method.

## 6. Typed cross-lane synthesis and stopping boundary

| Source | Target / map | Preserved | Lost or required | Decisive test |
|---|---|---|---|---|
| Full f | branch jumps + finite samples + uniform asymptotics | the recorded data | entire rapidly decreasing component | Theorem 1 changes C arbitrarily |
| Periodic a | mean mu / zeroth Fourier coefficient | exact convergence predicate | finite value K_a | Theorem 2 classifies convergence, not value |
| Finite-cutoff sum | exponential damping | mean and main log-squared growth | lower-order counterterm and finite shift | Theorem 3 adds -mu gamma ell and mu Q |
| MRB analytic formula | p-bit rational interval | certified value at precision p | symbolic/arithmetic classification | exact tail + conditioning + elementary bit model |
| A(8) | CM j=8000 curve | algebraic CM identification | full f and L(f) | inherited perturbation, strengthened here |
| Greedy A(m) dictionary | convergent signed amplitude series | starting real value | special MRB arithmetic origin | inherited expansion works for every small positive starting value |

### A minimal general certificate contract

For an object set X, an observation map O:X->Y, and a target T:X->Z, a function T_bar on O(X) satisfying T=T_bar composed with O exists if and only if T is constant on every fiber of O. Proof: necessity follows by evaluation; for sufficiency define T_bar(y)=T(x) for any x with O(x)=y, which is well-defined by the fiber condition. This set-level factorization test does not assert an effective algorithm, a runtime bound, or reconstruction of x.

For Theorem 1, X is the explicit affine family f+R*H_K, O consists of all listed branch/prefix/asymptotic observations, and T is the alternating sum. O is constant on X but T is nonconstant. Thus these observations cannot classify even this target value. By contrast the mean map on L-periodic filters does classify the convergence predicate by Theorem 2, while still not determining the convergent value. This is an actual positive/negative pair for a common certificate contract.

The shared research lesson with a local-to-global Jacobian analysis is limited and exact: before calling data a complete classifier, test whether two objects in the declared universe share those data while their target predicates differ. Here the universe and deformation are explicit and the target is a scalar alternating functional. There is no map carrying these analytic profiles to polynomial Keller maps, and no reduction to or from SAT has been supplied.

This yields a stronger version of the intended “radical picture”: a constant can be rapidly and rigorously computable, represented by several exact transforms, and possess striking algebraic shadows while its arithmetic identity and the sufficiency of its descriptive data remain wholly different questions. Numerical precision, representational richness, and complete classification are separate achievements with separate certificates.

**Remaining questions and stopping condition.** A proposed MRB modular identity must retain some global data not invariant under Theorem 1, or explicitly use the exact defining profile. Further precision, more finite samples, more asymptotic coefficients, and the same right-half-plane decay bound cannot by themselves fix the missing constant. No irrationality, transcendence, modularity, or nontermination theorem for the MRB greedy expansion is claimed. This lane stops after the strengthened obstruction, complete periodic convergence criterion, finite-part conversion, and certified precision audit; no speculative universal classification is promoted.
