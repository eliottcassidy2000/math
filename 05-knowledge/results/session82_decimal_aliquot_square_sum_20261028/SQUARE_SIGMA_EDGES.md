# A classical exceptional triple: prime-power aliquot square edges

Date: 2026-10-28. Status: elementary reductions and constructive infinity **PROVED**; global prime-power classification **COROLLARY OF CITED CLASSICAL THEOREM**; small inverse fibers **FINITE-EXACT WITH PROVED CUTOFF**. No novelty claim.

## 1. Result first

For prime p and integer a>=1, an aliquot step p^a -> s(p^a) is also an edge of the positive-integer simple square-sum graph **exactly** in these three cases:

| Source n | s(n) | n+s(n)=sigma(n) |
|---:|---:|---:|
| 3 | 1 | 4=2² |
| 81=3⁴ | 40 | 121=11² |
| 343=7³ | 57 | 400=20² |

The user's B-orbit is exactly 1001->343->57->23->1->0. Its only positive-label square-sum edge is 343--57. Its successive sigma values are 1344,400,80,24,1. The final 1->0 has square total but is not an edge in the graph whose vertices are positive integers. s(0) is not needed or defined here.

## 2. Precisely what is imported

Ljunggren's square-repunit theorem states that integers x>1,y>1,m>2 satisfying

    (x^m−1)/(x−1)=y²

are exactly (x,m,y)=(3,5,11),(7,4,20).

The statement and unrestricted integer domain are explicitly confirmed in the introduction, printed p257, of [Bugeaud–Mignotte–Roy, *On the Diophantine equation (x^n−1)/(x−1)=y^q*, Pacific J. Math.193(2),2000](https://msp.org/pjm/2000/193-2/pjm-v193-n2-p.pdf). Its reference11 is W. Ljunggren, *Noen setninger om ubestemte likninger av formen (x^n−1)/(x−1)=y^q*, Norsk Mat. Tidsskr.25(1943),17–20. This is the solved q=2 case, **not** the conjectural classification for arbitrary q. The original 1943 scan link could not be retrieved; the modern primary research article's explicit statement was inspected. No claim to have independently checked the deep original proof.

Apply this with x=p,m=a+1 for a>=2. For a=1, p+1=y² implies p=(y−1)(y+1), so primality forces y=2,p=3. Direct substitution verifies all three rows. This proves the classification, not a finite search extrapolation.

The number 1 is excluded from “prime power” here (a>=1). It has sigma(1)=1 and s(1)=0, not a positive graph edge.

## 3. Elementary cubic reduction and its exact boundary

For prime p,

    sigma(p³)=1+p+p²+p³=(p+1)(p²+1).

For an arbitrary positive integer x, study the **polynomial** F(x)=(x+1)(x²+1). It need not equal sigma(x³) when x is composite; even x=1 gives F(1)=4 but sigma(1)=1.

Euclid's algorithm gives gcd(x+1,x²+1)=gcd(x+1,2).

* If x is positive and even, the factors are coprime. Their product being square would make x²+1 square, impossible because x²<x²+1<(x+1)².
* If x is odd, x²+1 has 2-adic valuation exactly1. After dividing each factor by2, the two quotients are coprime. Thus F(x) is square iff

      x+1=2u², x²+1=2v²

  for positive integers u,v. Equivalently v²=2u⁴−2u²+1, with x=2u²−1. The solutions u=1,2 yield x=1,7. This reduction does **not** by itself prove that there are no other u; the global uniqueness uses §2 with m=4.

Consequently positive polynomial solutions are x=1,7. Over all integers they are x=-1,0,1,7: x<-1 makes F(x)<0, while the boundary values give0,1. The y-values are0,±1,±2,±20 as appropriate.

## 4. The exact general edge predicate and multiplicative transport

For every n>=1, n+s(n)=sigma(n). Hence the arithmetic sum is square iff sigma(n) is square. To call it an edge of the positive simple graph one additionally needs s(n)>0 and s(n)!=n.

Let sf(k) denote the squarefree kernel (product of prime factors having odd exponent). For coprime positive integers a,b,

    sigma(ab)=sigma(a)sigma(b),
    sigma(ab) square <=> sf(sigma(a))=sf(sigma(b)).

Thus square-sigma factors may be combined when coprime; more generally two matching nontrivial parity signatures cancel. Example: sigma(5)=6, sigma(23)=24, neither square, but sigma(115)=144. Multiplication transports the edge-total predicate, not the aliquot value:

    s(ab)=a*s(b)+b*s(a)+s(a)*s(b)   (gcd(a,b)=1).

The coprimality guard is essential: sigma(3)=4 and sigma(81)=121 are square, yet sigma(243)=364 is not. Powers are not a closure rule either: sigma(343²)=137257 lies strictly between370² and371². A square-edge step does not persist along the orbit: 81->40 is square-sum, while40->50 has total90.

## 5. Infinitely many coprime building blocks: constructive proof

[Beukers–Luca–Oort, *Power Values of Divisor Sums*, American Mathematical Monthly119(5),2012,373–380](https://dspace.library.uu.nl/bitstream/handle/1874/272370/amer.math.monthly.119.05.373.pdf) proves infinitude in Theorem2.1 and explains the exponent-parity method. No new attribution is claimed. Here is a self-contained version sufficient for the coprimality needed here.

**Claim.** For each integer B>=2 there is a squarefree m>1 all of whose prime factors exceed B and whose sigma is square.

For j>=2 take the odd primes B<p<2^j. Every prime factor of p+1 is at most2^(j−1). Encode the exponents of p+1 modulo2 as a vector in F₂^r, where r=pi(2^(j−1)). There are pi(2^j)−pi(B) vectors. If this number exceeds r, a nonempty subset has vector sum zero. Its distinct-prime product m is squarefree, supported above B, and sigma(m)=product(p+1) is square.

Such a j exists: the dyadic increments pi(2^j)−pi(2^(j−1)) are unbounded. Otherwise a uniform bound C would give

    sum over primes p of 1/p <= finite initial part + sum_j C/2^(j−1) < infinity,

contrary to Euler's divergence of prime reciprocals. For completeness, that divergence follows from Euler products: if the prime reciprocal sum converged, product_p(1−1/p)^(-1) would be bounded, but the partial products dominate every initial harmonic sum once all its prime factors are included. The harmonic series diverges. This establishes the claim without a prime-distribution conjecture.

Repeat with B larger than all previously used primes to obtain pairwise coprime blocks m₁,m₂,... . Every nonempty finite product of distinct blocks has square sigma. This is an infinite family by genuinely new prime support, not repeated powers of one seed.

Taking B>=7 also gives infinitely many n=343m with sigma(n)=400 sigma(m) square. These n are odd >1 with sigma(n) even, so s(n)>0 and n!=s(n): equality would make the square sigma(n)=2n congruent2 modulo4. Thus all are genuine positive simple square-sum aliquot edges.

The checker constructs 12 disjoint blocks by exact binary Gaussian elimination. The first is m=11*17*23=4301, sigma(m)=5184=72², s(m)=883. Another is97*127=12319, sigma=12544=112², s=225. These examples do not assert a fixed number of factors in all blocks, prime targets, or an infinite family of prime cubes.

## 6. Observation fibers and the untouched-prime boundary

Define F_t={n>=1:s(n)=t} and F_t^square={n in F_t:n+t is square}. For t>1 every source is composite. If p is its least prime factor, n/p is a proper divisor at least sqrt(n), distinct from1. Therefore

    t>=1+n/p>=1+sqrt(n), hence n<=(t−1)².

This proves a complete finite search bound; it is not a guessed cutoff. Exact divisor-sum sieves give:

| t | complete F_t | F_t^square |
|---:|---|---|
| 3 | {4} | empty |
| 7 | {8} | empty |
| 11 | {21} | empty |
| 13 | {27,35} | empty |
| 40 | {44,74,81} | {81} |
| 57 | {99,159,343,559,703} | {343} |

Also F₁ consists exactly of primes, so F₁^square={3}, by the elementary a=1 argument. Thus these exceptional sources are unique square-edge sources **within their own complete fibers**, not unique aliquot sources generally.

The observation t alone loses the square-edge predicate: s(159)=s(343)=57, but their sigma values are216 and400. A sidecar such as sf(sigma(n)) retains the missing information.

A prime target can have a square-edge source: s(115)=29 and115+29=144. Conversely the touchable prime3 has no square-edge source, because its only source4 has total7. Thus square-edge touchability is strictly stronger than ordinary touchability, and its classification cannot prove the conjecture that2 and5 are the only untouchable primes. Untouchability excludes every source shape; §1 restricts only prime-power sources, and §5 does not control targets.

A further exact bridge exists in the distinct-semiprime branch for odd target t. Set h=(t−1)/2 and primes p=h−r<q=h+r, with integer r>=1. Then n=pq=h²−r² has s(n)=p+q+1=t, and

    n+t=z² <=> r²+z²=(h+1)².

This is a Pythagorean equation with a separate primality guard on h±r. It is not Goldbach, not a claim that every Pythagorean triple produces primes, and it omits repeated-prime and other factorization branches. For example t=29,h=14,r=9 gives p=5,q=23 and9²+12²=15². The exact map and all its guards, rather than a visual resemblance, supply the connection.

### Requested inverse fibers: a structural exclusion, not just a filter

The independently audited generated `results.json` contains complete inverse fibers of targets1001,1000001,1001000, with respective sizes29,5408,2. Generate it with the adjacent `check_aliquot_fibers.py` before running this checker. The checker validates target counts and exception lists, regenerates every full fiber using the producer function, and records the generated input hash. Filtering every listed n by exact integer-square-root arithmetic gives **empty F_t^square for all three targets**. Completeness is inherited from that lane's proved inverse algorithm and audit, not re-proved by this filter.

For the two odd targets there is a structural explanation of the entire distinct-semiprime branch. If c has no prime divisor congruent1 modulo4, then c cannot be the hypotenuse of a nondegenerate integer right triangle. Indeed divide the three coordinates by their common gcd to obtain a primitive triangle. Its hypotenuse c₀>1 is odd; for any prime q dividing c₀, neither leg is divisible by q, so their squared ratio is−1 modulo q. Hence q=1 modulo4, contradiction. (The usual parity argument shows a primitive hypotenuse is odd.)

Here c=(t+1)/2 is501=3*167 or500001=3*166667. All displayed factors are primes congruent3 modulo4, with primality certified by trial division. The Pythagorean equation in §6 therefore excludes every distinct-semiprime square-edge source. The remaining exceptions are checked separately: one source11191 at t1001; six sources2660799,13011919,16390159,30001459,77519461,190124401 at t1000001. The even target1001000 has only sources718720 and1001998998001; both fail the same exact square test. The certificate records bounding consecutive square roots for each exception.

The hypotenuse condition alone is not sufficient after primality guards. For c=5 (t=9,h=4), the positive leg choices r=3,4 give centered pairs(1,7),(0,8), neither a prime pair. For c=13 (t=25,h=12), r=5 gives primes7,17: n=119,s(n)=25,n+s(n)=144=12². Thus the geometric obstruction and the prime constraints are separate and both necessary.

### A guarded 3–4–5 subfamily

Scale the right triangle by k: c=5k,r=3k,z=4k. Then h=5k−1, p=2k−1 and q=8k−1. Whenever p,q are prime (k>=1, hence distinct),

    n=(2k−1)(8k−1), s(n)=10k−1, sigma(n)=16k².

At k=3 this gives115->29; at k=6 it gives517->59. Both target values happen to be prime. Requiring p,q **and**10k−1 all prime is a three-linear-form primality condition: infinitude is not proved here. In particular this guarded identity must not be presented as an unconditional infinite prime-target family. The small fibers3,7,11,13 above already have empty square subfibers despite being touchable primes, unlike the proposed ordinary prime-untouchable list{2,5}.

## 7. Inheritance and verification

The recovered repository note `05-knowledge/results/aliquot_table_untouchable_primes_chains_seven_twentyone_20261026.md` was recovered and inspected. It already identifies s(p³)=1+p+p², warns against using composite bases in that prime-power formula, and preserves the distinct-prime-Goldbach condition for the global untouched-prime statement. Its numeral coincidences with tournament forbidden values are explicitly not equivalences. This report does not alter those boundaries. The translated-square-sum-triple companion supplies the distinct-label positive graph convention and translated-triple route, not a proof of this number-theoretic classification.

Reproduce from this directory:

    python3 check_aliquot_fibers.py
    python3 independent_prefix_audit.py
    python3 check_square_sigma.py > certificates.json
    python3 -O check_square_sigma.py > certificates.optimized.json
    cmp certificates.json certificates.optimized.json

The standalone small tests use only the standard library; the inverse-fiber producer requires SymPy. Keep the three scripts adjacent, or pass `--data-dir DIRECTORY` to this checker to locate the producer and generated results. Missing dependencies produce a clear error rather than accepting a partial result. Explicit checks remain active under -O. Tests include polynomial x from−100 to10^6; all prime-power sources n<=10^6; independent divisor enumeration through5000; gcd/Pell reductions through10000; 12 coprime squarefree blocks; all B-orbit edges; hostile multiplicativity examples; and complete fibers3,7,11,13,40,57 with proved bounds. The finite checks do not certify Ljunggren's global theorem: that dependency remains visible and cited.

## Independent audit

The cited Ljunggren theorem domain and its square-case solutions were checked directly in the primary research PDF. A separate audit verified the finite fibers, all5439 requested sources, exact primality and Pythagorean obstructions, and the self-contained coprime-family proof. The final generated-input checker was rerun normally and under optimization with identical output. The infinite classification depends on the explicitly cited theorem, not the finite tests.
