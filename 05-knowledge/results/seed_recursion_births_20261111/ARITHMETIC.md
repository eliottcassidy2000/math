# Birth fibres, diagonal forests, and square graphs

**PROVED IN NOTE for the elementary deductions below; CITED recovery is
identified. No prime-infinitude or tournament-coverage claim follows.**

## A1. What happens when integer n is born

Keep parent pairs, rather than only arrows from individual parents. For
addition, parents are positive integers a,b with a+b=n. For causal
multiplication, use proper nonunit parents 2<=a,b<n with ab=n; the pair
{1,n} cannot generate n from earlier nodes. Weak permits equal parents;
strict requires distinct parents. The swap (a,b)->(b,a) gives the counts

| At new node n | Weak unordered parents | Strict unordered parents | Removed diagonal |
|---|---|---|---|
| Addition | floor(n/2) | floor((n-1)/2) | {n/2,n/2} if n even |
| Proper multiplication, n>=2 | (tau(n)-2+1_square(n))/2 | (tau(n)-2-1_square(n))/2 | {sqrt(n),sqrt(n)} if n square |

These fibre formulas are CITED from
[THM-2422 / operation-fibres-summand-closure-and-twin-center-ancestry](../../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md).
At n=1 the proper multiplicative fibre is empty; do not apply its n>=2
formula there. The diagonal is the fixed-point set of the same swap
symmetry, with targets 2a in addition and a^2 in multiplication.

Let D(N)=sum_(n<=N) tau(n)=sum_(a<=N) floor(N/a). Summing the multiplicative
formulas over 2<=n<=N gives the exact cumulative birth counts

    weak:   (D(N)-2N+floor(sqrt(N)))/2,
    strict: (D(N)-2N-floor(sqrt(N))+2)/2.                    (A1)

The difference is floor(sqrt(N))-1, one for each nonunit square.
The additive cumulative strict count is floor((N-1)^2/4); the weak count
is floor(N^2/4), differing by floor(N/2). These are finite-sum identities,
not asymptotic heuristics. A divisor sieve supplies every multiplicative
birth count through N in O(N log N) elementary increments; the divisor
summatory function can also be grouped by its O(sqrt(N)) distinct quotients.

## A2. Two diagonal forests, with different roots

Recovering [arithmetic_braids_20260917_summand.md, Sections 1–2](../arithmetic_braids_20260917_summand.md):
the weak additive loopless shadow is every x->z with x<z. The strict
shadow removes exactly x->2x. The removed forest has rays rooted at each
odd u, namely u,2u,4u,... . The weak multiplicative loopless shadow is
proper divisibility; its strict version removes x->x^2 for x>=2.
Those rays are rooted at nonsquares b>=2, namely b,b^2,b^4,... .

The factorization coordinates explain the distinction. If
n=product p^(e_p), write k=v_2(gcd{e_p:e_p>0}). Then n=b^(2^k) with
b nonsquare, uniquely. This root b is **not** the squarefree kernel:
12 is nonsquare, while its squarefree kernel is 3. The two quotients
have different preserved information and must not be conflated.

The shadow's edge x->z certifies only one parent. Closure needs the
other parent too. On labelled integer graphs one can reconstruct it
(z-x or z/x), but ordinary reachability still imposes an OR condition
where generation requires both parents. Tournament seed construction
has a similar logical boundary: a scalar arrow cannot ensure that its
next required marked source exists. This is a typed analogy, not an
isomorphism or a reduction proving coverage.

## A3. Minimal generators and the strict {2,3} comparison

For causal weak multiplication on integers >=2, the indispensable
generators are exactly the primes. For strict multiplication they are
exactly primes together with prime squares. Indeed these and only these
have no distinct proper factor pair: if n is composite and not p^2,
its smallest prime factor p satisfies n/p>p. Both parents are smaller
than n, so induction proves sufficiency of the stated generators.

Thus the weak/strict diagonal decision changes the minimal generator
family, not merely a local edge count. With starting set {2,3}, strict
multiplicative closure is exactly

    {2,3} union {2^a 3^b : a,b>=1}.                        (A2)

Necessity: no new prime support can appear; no pure 2-power or 3-power
beyond the initial prime has a first legal distinct-parent birth.
Sufficiency: start with 6, then multiply repeatedly by 2 and 3, each
distinct from the mixed parent. In contrast strict additive closure
from {2,3} is N\{1,4,6}, as proved in THM-2422. The infinite holes in
(A2) and the three additive holes have different mechanisms.

The same proof gives, for arbitrary distinct 2<=a<b,

    strict product closure({a,b})={a,b} union {a^i b^j:i,j>=1}.

No multiplicative independence is needed; different exponent pairs may
represent the same value. Concurrent
[Session259, Section 2](../session259_arithmetic_graph_growth_20261111/SYNTHESIS.md)
proves the corresponding strict sum formula
{a,b} union {ia+jb:i,j>=1}. The exact shared mechanism is that every
nonseed expression contains both seed types; adding or multiplying by
either smaller seed then constructs every allowed positive coefficient
pair. Numerical collisions and birth times differ. Session259 also
independently identifies prime-square atoms, so that overlap is credited.

This explains how a missing diagonal can create new required roots.
It motivates recording marked seed deletion roots, but supplies no
claim that their sequence is the prime or prime-square sequence.
The tournament's first roots even include c=1, so a simple identification
of its roots with arithmetic primes is already false.

## A4. Square-sum and square-product graphs have exact growth laws

These are undirected graphs on {1,...,N}, with distinct vertices joined
when their sum or product is a square. There is no intrinsic orientation
between a pair, so no cosmetic tournament is imposed.

For the sum graph Q_N, the new vertex n has older neighbours s^2-n,
where n<s^2<2n. Its birth degree is therefore

    floor(sqrt(2n-1))-floor(sqrt(n)).                      (A3)

The prior [three-chain recovery](../session104_descent_exterior_20261029/summands/THREE_CHAINS_RECOVERED.md)
proves that Q_12 consists of

    10--6--3--1--8;  2--7--9;  11--5--4--12.

At 13 the neighbours 3,12 merge the first and third chains; at 14 the
neighbours 2,11 connect all vertices. For every n>=6, the next square
after n is <2n (since 2sqrt(n)+1<n). Hence each later vertex has an
older neighbour, so Q_N stays connected for all N>=14. Founders are
{1,2,4}, not the additive-closure holes {1,4,6}. Excluded loops occur
at n=2j^2 because 2n is a square; they do not characterize leaves.

For the product graph R_N write n=d a^2 with d squarefree, uniquely.
Then

    xy is square iff sf(x)=sf(y),                         (A4)

where sf is the squarefree kernel (prime-exponent parity vector). Thus

    R_N = disjoint union_(d<=N squarefree) K_floor(sqrt(N/d)).

This completely classifies its components at every N. At birth n=d a^2,
the vertex joins exactly a-1 older vertices of its clique. A new component
is born iff n is squarefree. Consequently

    components(R_N)=sum_(d<=N) mu(d)^2,
    edges(R_N)=sum_(d<=N squarefree) binom(floor(sqrt(N/d)),2).

Unlike Q_N, R_N acquires infinitely many components: each prime begins
one. Modulo square classes, multiplication adds exponent-parity vectors,
doubling toggles the prime-2 coordinate, and squaring sends every vector
to zero. To recover size and birth degree one must also retain a in
n=d a^2. This is a direct, rigorous use of a finite bit pattern for each
integer, not a finite-state description of all prime supports at once.

`arithmetic_births.py` checks all formulas and both {2,3} closures through
4096, with 512 independent direct pair-enumeration controls through 128.
At N=4096, R_N has 2491 components and 8614 edges; vertex 4096 joins
63 older vertices. The JSON records the first 512 full birth profiles
and the 4096 endpoint. These finite checks accompany the proofs above.

## A5. What transfers to the tournament enumeration program

| Source and map | Preserved predicate | Lost information and needed sidecar | Cheap decisive test |
|---|---|---|---|
| Parent pairs -> one-parent shadow | an eligible partner exists | simultaneous availability of both parents; keep full pair fibre | strict {2,3} closures |
| Multiplicative fibre -> diagonal forest | equal parents | all other factors; retain divisor fibre | 12 versus 16 |
| Integer -> squarefree class | product is square | size and birth time; retain square multiplier a | 3 versus 12 |
| Seed (B,O) -> odd coefficient c | its eventual tail address | valuation, port response, graph identity; retain all three | c=1 singleton versus five-vertex deletion root |
| Seed -> marked isomorphism class | every graph construction predicate | labels only; retain a relabelling witness | actual automorphism action |
| Fixed core -> full inflation tensor | every permitted block-size count | exterior alternatives; retain template identity | 23->45 but 59->129 |

For A000568, the practical improvement here is enumerating strong marked
cores and genuine automorphism orbits once, then applying the proved
padding convolution. This is an improvement for the refined seed
problem; it does not supply a faster general enumeration of all ordinary
tournaments or a new A000568 term.

For [A135244](https://oeis.org/A135244), the largest inverse of the proper
divisor sum, factor fibres should retain their **sum of divisors**, not
only tau(n), squarefree class, or the existence of a factor pair. The
existing [complete inverse work](../uniform_return_metagraph_20261111/arithmetic/PROOF.md)
already separates prime-square, distinct-semiprime and higher-factor
sources. Its exact bridge to [A002496](https://oeis.org/A002496), primes
r^2+1, is

    r^2+1 prime iff A135244(r^2+2)=(r^2+1)^2.

That result is inherited, not newly proved here. The new birth bookkeeping
explains why the square branch needs its own record; it does not replace
the existing inverse algorithm or its primality certificates. There is
no claim here of new terms, a faster measured runtime for those sequences,
or infinitude of primes r^2+1.

## A6. Integration with the incoming arithmetic growth work

During this session, [Session259's synthesis](../session259_arithmetic_graph_growth_20261111/SYNTHESIS.md)
arrived on the shared branch. It separates parent arrival, target
arrival, synchronous generation, and topological lifetime. This improves
the present dictionary: a seed deletion chain is a generation clock,
its vertex count is a size clock, and the value H is neither of those.
Likewise the graphs Q_N and R_N use a maximum-vertex clock; their edges
are not births of the sum/product targets themselves.

Its sum-product transport keeps the pair via (s,p)=(a+b,ab); the
discriminant s^2-4p=(a-b)^2 recovers the unordered parents exactly.
For 2<=a<b, both projections are causal in numerical size. A future
tournament counterpart must likewise retain a reconstructible interface;
keeping H alone fails at R5's explicit example. No such universal
tournament counterpart is supplied here.

The incoming divisor topology uses rad(n), the product of all distinct
prime factors. Our square-product graph uses sf(n), the product of
prime factors with **odd** exponent. They already differ at 12:
rad(12)=6, sf(12)=3. Neither is the root of the squaring forest, which
is 12 itself. This three-way separation prevents transferring a theorem
between incompatible quotients merely because all use squares or primes.

The next incoming checkpoint,
[Session260's modular clone note](../session260_quasiperfect_chanpak_20261111/bridges/RESEARCH.md),
gives a further precisely scoped bridge: for m>=2, replacing any rooted
tournament vertex by TT_m preserves H modulo m for every source exactly
when m is prime or a Carmichael number. This is CITED from that local
note, not used as a dependency here or independently audited in this
packet. It uses the same Stirling coefficients and provides a congruence
filter, whereas R1 requires exact defect cancellation. Carmichael numbers
in this statement are distinct from Carmichael's totient-multiplicity
conjecture; neither that conjecture nor exact Hamilton coverage is settled
by the congruence criterion.

The final rebase brought the
[finite-seed continuation](../summand_multiplicand_growth_20261111/README.md)
at 8cc23df026. Its general positive-weight monoid theorem extends A3's
two-seed closure to arbitrary finite seed sets; its prime-power map
a->p^a carries additive doubles to multiplicative squares while preserving
synchronous generation. These are CITED incoming extensions, not inputs
to R1–R11. They sharpen the next tournament obligation: an arithmetic
closure description needs a graph operation preserving the marked seed
predicate before its finite generator proof can transfer.
