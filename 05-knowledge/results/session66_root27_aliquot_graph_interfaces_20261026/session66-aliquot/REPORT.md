# Session 66: inverse aliquot fibers, certified searches, and the density boundary

**Status:** PROVED elementary reductions; CITED prior algorithms and density bounds; FINITE-EXACT computations. **OPEN:** a second odd untouchable, uniqueness of 5 among odd untouchables, and existence/exact value of the natural density. No new solution or algorithmic priority claim.

## 1. The user's association is exactly an inverse-fiber construction

Put `s(n)=sigma(n)-n` for positive integers, and `F(m)={n>=1:s(n)=m}`. The supplied examples read:

|m|F(m)|
|--:|:--|
|0|{1}|
|1|all primes|
|2|empty|
|3|{4}|
|4|{9}|
|5|empty|
|6|{6,25}|
|7|{8}|
|8|{10,49}|
|9|{15}|
|10|{14}|
|11|{21}|
|12|{121}|
|13|{27,35}|

These are **complete fibers**, not just listed witnesses. The checker establishes completeness using the bounds below. In particular, the transcript's `6 and 25 with 6`, `49 and 10 with 8`, and `35 and 27 with 13` are consistent, not arithmetic errors. Fiber 1 is the exceptional infinite fiber. Untouchable means a **positive target** with empty fiber; 0 is not an untouchable.

## 2. Inheritance and concept board

Read-only repository pin: `63f52936f8dd8592d4f2fee869748baad2e1997e`, clone `session65-repo`. Startup references read: `AGENTS.md`, `00-navigation/START-HERE.md`, `00-navigation/RESEARCH-PROTOCOL.md`, `01-canon/ACTIVE-GUARDRAILS.md`, and `05-knowledge/reference/CORE-PAPERS.md`. No mutation or publication of that clone was performed.

Closest inherited mechanism: `05-knowledge/results/aliquot_dyadic_frontier_20261014e.md`, Section 2, already gives the guarded prime-extension identity `s(ap)=s(a)p+sigma(a)` and stresses that rational/integer candidates require primality and coprimality. The least-used sidecar here is the **prime-exponent-one versus squarefull** partition. The inherited parity computation in `aliquot_parity_squares_20261014.out` concerns cycles; its finite cycle census does not answer inverse-fiber completeness or target density.

Concept board: (1) inverse fiber versus forward orbit; (2) source cutoff versus target cutoff; (3) prime-extension guards; (4) squarefull complement; (5) density of occupied targets rather than multiplicity; (6) multiplicative monoid versus a noncongruent fiber quotient. Anchor: certify or eliminate odd-target candidates. Niche: exact even-target census. Wildcard: the precise part of the F1 analogy.

## 3. Elementary complete-search bounds

### Lemma A (sharp global bound)

If `m>1` and `s(n)=m`, then `n` is composite and

`n <= (m-1)^2`.

Let p be its least prime factor and q=n/p. Both 1 and q are distinct proper divisors, so `q<=m-1`; also `p<=q`, proving the bound. Equality forces p=q=m-1 and is attained by `n=p^2`, `m=p+1` for every prime p.

If n is not the square of a prime, 1,p,q are distinct proper divisors, so `p+q<=m-1` and

`4n=4pq <= (p+q)^2 <= (m-1)^2`.

This proves actual finiteness of every fiber except F(1), not a merely empirical search window. It proves F(2) is empty. For m=5 the non-prime-square bound gives n<=4, while the only composite there,4, is a prime square; the prime-square possibility would require prime p=4. Hence F(5) is empty.

### Lemma B (parity and linear even-source bound)

The factorization formula shows sigma(n) is odd iff n is a square or twice a square: the factor for 2 is odd; an odd-prime geometric sum is odd exactly when its exponent is even. For even n>=6, the distinct proper divisors 1,2,n/2 give `n<=2s(n)-6`. Consequently an even preimage of an odd m>=5 is a square or twice a square and is at most `2m-6`.

For an even target, an odd source must be an odd square `a^2`. If a is prime its image is a+1. If a is composite, let p be its least prime; `p<=sqrt(a)` and `a^2/p` is proper, giving `s(a^2)>=a^(3/2)` and `a<=m^(2/3)`.

**Priority:** this last prime/composite root split and its near-linear census application are explicitly in Pomerance–Yang, Section5, not new here. [Primary manuscript](https://math.dartmouth.edu/~carlp/uupaper3.pdf).

## 4. A simple complete all-parity inverse certificate

Every composite n either has a prime p occurring to exponent exactly one, or is squarefull (every nonzero prime exponent is >=2).

* In the first case write n=ap with gcd(a,p)=1. For m>1, a>=2, and the distinct proper divisors1,a give a<=m-1. The inherited identity forces
  `p=(m-sigma(a))/s(a)`.
  Scan `2<=a<=m-1`, require the quotient to be an integer prime and p not dividing a, and verify the source.
* Every squarefull n has a **unique** expression `n=u^2 v^3` with v squarefree. At each prime, an even exponent is placed entirely in u², while an odd exponent>=3 leaves exponent3 in v³. **There is no gcd(u,v)=1 condition.** For example32=2²·2³.
  Enumerate `v<=cuberoot((m-1)^2)` squarefree and `1<=u<=sqrt((m-1)^2/v^3)`, then test s(n)=m.

This is complete by the partition, with at most m-2 core candidates and at most

`(m-1) sum_{v>=1} v^(-3/2) <= 3(m-1)`

squarefull pairs. This is an **O(m)-sized candidate universe**, not an unexplained claim that all arithmetic operations take unit cost. The supplied SPF-backed implementation factors only u,v<=m in the squarefull path, and its elementary arithmetic costs can be bounded by polynomial logarithmic factors per candidate. It is a transparent certificate generator, not the fastest known even-fiber algorithm.

**Known stronger even-target algorithm:** Booker proves `O_epsilon(m^(1/2+epsilon))` time for complete even-target inversion (Theorem2.2), and also uses the following two-prime divisor trick. [Primary paper](https://arxiv.org/pdf/1610.07471).

### Fixed-core hyperbola trick (known, rederived and checked)

For a>=2, set S=s(a), T=sigma(a). For distinct primes p,q coprime to a,

`s(apq)=m  iff  (Sp+T)(Sq+T)=Sm+aT`.

Thus divisors d,e of the right side, with d*e=Sm+aT and both congruent to T modulo S, generate candidates p=(d-T)/S and q=(e-T)/S. Retain p<q prime and gcd(a,pq)=1. The identity follows by expanding `s(apq)=S pq+T(p+q+1)`. Booker's version divides through by gcd(S,T); this is the same mechanism, not a new identity.

For prime-power cores a=q^k, the one-prime identity becomes

`m-1=(1+q+...+q^(k-1))(p+q)`,  p prime, p!=q.

This produces additional explicit odd witnesses besides the k=1 Goldbach witnesses. Failing this family does not establish untouchability.

## 5. The density question: what is and is not known

Let U(X) count untouchables in [1,X]. Almost all odd targets are touchable; the distinct-prime Goldbach conjecture would imply that5 is the only odd untouchable. Ordinary Goldbach does **not directly** supply a distinct pair: `s(p^2)=p+1`, not2p+1. The implication is one-way; ruling out every Goldbach pair does not rule out all other source shapes. The published verification to4·10^18 therefore excludes another odd untouchable up to4·10^18+1. Booker's footnote2 makes the distinctness detail explicit: the large computation finds one prime<=9781, reducing potential equal-pair exceptions to2p for primes5<=p<=9781. We independently checked all1204 of those small prime cases, but did **not** rerun the large computation. [Booker, footnote2](https://arxiv.org/pdf/1610.07471).

Chen–Zhao Theorem1 gives explicit constants g_M with `liminf U(X)/X>=g_M`; their displayed choice gives `g_M>0.0602757`. This is a lower-density theorem, not existence of a density. [Primary paper](https://publi.math.unideb.hu/paper/1551/download/10_5486_PMD_2011_4820.pdf).

### Certified stronger instantiation of the Chen–Zhao bound

For `M=2^6 3^5 5^4 7^3 11^2 13·17·19·23·29·31·37·41`, the script evaluates their explicit sum by rounding each **nonnegative rational term down** to a multiple of10^(-24), entirely with integers. Adding one unit of that scale per positive term yields a strict upper bound (the generic implementation adds one unit if there are no positive terms). It obtains

```
0.060275780518518619919468 <= g_M
                              < 0.060275780518518620563502;
0.060452936134095005620815 <= g_(43M)
                              < 0.060452936134095006909185.
```

Consequently the cited theorem gives the rigorously certified conclusion

`liminf_(X->infinity) U(X)/X > 0.060452936134095`.

This is a new **session calculation/instantiation of an established theorem**, not a claim of the best bound in the literature and not a new asymptotic-density formula. The complete rational enclosure is recorded in `results/density_lower_bound.json`;1,290,240 divisors are used for43M.

The improvement has an algebraic reason. Write `h(d)=max(0,1/(2d)-1/s(2d))` and `w_M(d)=phi(M/d)/(M/d)`. If q is an odd prime not dividing M, then

`g_(qM)-g_M = sum_(d|M) w_M(d) [h(qd)-h(d)/q] >= 0`.

Indeed `s(2qd)=q s(2d)+sigma(2d)`, so each bracket is nonnegative. Whenever s(2d)>2d, its exact value is

`sigma(2d)/[q s(2d) (q s(2d)+sigma(2d))] > 0`.

More generally, let `H(d)=d h(d)=max(0,1/2-d/s(2d))`. The same constant has the exact representation

`g_M = E[H(gcd(R,M))]`, where R is uniform modulo M.

If M divides M', lifting R uniformly modulo M' gives `gcd(R,M)|gcd(R,M')`. The ratio sigma(2d)/(2d) is nondecreasing under divisibility of d, so H(d) is also nondecreasing. Hence **g_M<=g_M' whenever M|M'**. This coupling proof (root-agent synthesis) explains the entire divisibility-directed lower-bound search.

Since6 divides the chosen M and s(12)=16>12, the increase is strict for every new odd prime. This explains the successful43 test and supplies a principled next search move, rather than numerically guessing a density from a finite prefix.

Luca–Pomerance prove a positive lower proportion of **even** touchable targets. Combined with the density-zero odd exception set, this gives `limsup U(X)/X<1/2`; the 2012 Pomerance–Yang statement that this strict inequality was unknown is superseded. [Primary paper](https://math.dartmouth.edu/~carlp/LucaPomeranceAA.pdf).

Pollack–Pomerance propose the model constant

`Delta = lim_{Y->infinity} [1/log Y] sum_{a<=Y, 2|a} exp(-a/s(a))/a`.

They prove **this model limit exists**. It is conjectural that U(X)/X converges and equals Delta; their model/evidence points near0.17. Computing Delta accurately would not prove the actual density formula. The published errata concern other propositions, not this separation. [Paper, Section3](https://math.dartmouth.edu/~carlp/btran10.pdf); [errata](https://www.pollack-math.net/reversal-errata.pdf).

The large published census by Chum–Guy–Jacobson–Mosunov reaches2^40; our100000 census below is a small independent validation, not a record. [Primary report, Section3](https://arxiv.org/html/2110.14136).

### A useful density reduction and a finite-core warning

Let E(X) be the number of even m<=X with no **even-source** preimage. Lemma B and the prime/composite-root split imply

`0 <= E(X)-#{even untouchables<=X} <= pi(X-1)+floor(X^(2/3))`.

The right side is o(X). Since odd untouchables also have density0, `U(X)=E(X)+o(X)`. Existence/value of the actual density is therefore equivalent to the limiting missing-target proportion of the even-source map, whose sources need only extend linearly to2X. **This is not a proof that either limit exists.**

Every fixed finite bank of prime-extension cores a<=C produces only `O_C(pi(X))=o(X)` targets. The even targets reached from squarefull sources also have density0: even squarefull sources are at most2X and number O(sqrt(X)); odd sources are the already bounded odd squares. Thus a finite core bank plus the squarefull complement cannot determine the positive-density even-touchable image. The unbounded core coordinate and collisions of its affine prime families are genuine missing data, not a numerical precision issue. This is a rigorously identified obstruction to a tempting finite-template density argument.

## 6. F1: precise surviving content and a decisive hostile

The positive integers under multiplication are the free commutative monoid on the primes. Divisors of `n=product p^e_p` correspond to the product of exponent chains `[0,e_p]`. Their multivariate enumerator is `product_p(1+X_p+...+X_p^e_p)`; subtract the top monomial and evaluate X_p=p to obtain s(n). This supplies a precise monoid/divisor-combinatorics perspective, without asserting an actual one-element field.

But the relation `n~n' iff s(n)=s(n')` is **not a multiplicative congruence**: s(2)=s(3)=1, whereas s(2·2)=3 and s(3·2)=6. Therefore the fiber labels cannot simply inherit multiplication from the integers. An F1 analogy that forgets the prime weights or this failure does not preserve the inverse-image question. The set of fibers is a useful quotient of a set, not automatically a quotient monoid, field, or scheme.

## 7. Reproducible validation and adversarial controls

Run from this directory:

```
python3 scripts/verify.py
python3 -O scripts/verify.py
python3 scripts/density_lower_bound.py
python3 scripts/verify_density_monotonicity.py
```

`verify.py` uses explicit exceptions, not removable assertions. Artifacts:

* `results/small_fibers.json`: complete F(m),2<=m<=500, matched against an independent divisor sieve through the proved common bound249001.
* `results/verification.json`:458461 composite sources<=500000 checked for inverse bounds/parity;5771 hyperbola cases a<=30,m<=200 matched against complete fibers.
* `results/untouchables_100000.json`: complete prefix of13863 untouchables;13862 even plus5. Every other odd target has an explicit checked distinct-prime witness (plus the small exceptions1,3,7). Fraction0.13863 is only a finite fraction.
* A naive input scan n<=100000 misses witnesses for9166 even targets<=100000, demonstrating why finite source and target cutoffs cannot be confused.

The independent audit caught a generic-input boundary in the initial density evaluator: odd M lacked the factor sigma(2)=3. It is repaired, with the exact positive control g_9=1/126; the two even-M headline calculations were unaffected. A second audit control ensures the strict upper endpoint remains strict when the exact sum is0 (for example M=1).

Hostiles retained: omitting p∤a makes a=2,p=2 falsely claim s(4)=5; forcing gcd(u,v)=1 misses32∈F(31); treating squarefull as square misses27∈F(13); applying a linear bound to all sources misses121∈F(12); multiplying fiber representatives is not well defined.

The provenance/source ledger and integer-certified Chen–Zhao evaluation are under `sources/` and `results/`. The independent audit is separately authored by the audit agent (parent-owned integration may relocate it). No computational run in this packet proves a global absence of a second odd untouchable or a density limit.
