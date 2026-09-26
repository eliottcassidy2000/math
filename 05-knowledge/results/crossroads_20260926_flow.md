# Source-to-hub flow: affine offsets, trimmed capacity, and the sharp descent price

**Status:** PROVED elementary cylinder, energy, affine-offset and capacity lemmas below; FINITE-EXACT experiments and independent inverse-tree controls. The sharp-price consequence uses the explicit near-critical tube population construction in section 6. No Collatz convergence claim follows. This lane does not reserve a theorem identifier or alter a hypothesis status.

## Inheritance, portfolio, and concept board

Anchor: the actual source-to-hub joint law in Collatz descent certificates. Niche: recover the repository's word-stratified Hall and transportation machinery. Wildcard: affine cohomology, where a multiplicative cocycle leaves a small additive offset.

Closest proved mechanism: [THM-4475, price of provable descent](../../01-canon/theorems/THM-4475-price-of-provable-descent-tends-to-zero.md) and [THM-4477, Cauchy-Schwarz price](../../01-canon/theorems/THM-4477-cauchy-schwarz-price-bound-and-am-fair-criticality.md). The canonical hostiles are the odd cycles `5 -> 7 -> 5` for `3n-1` and `13 -> 33 -> 83 -> 13` for `5n+1`. Corrected near miss: conditioning an energy on bad sources improves finite constants but still leaves a tilted diagonal obstruction. Least-used sidecar: the normalized ordered affine offset, retained together with time and odd-step count.

Recovered older mechanisms:

* [THM-2545, word-stratified Hall arrival](../../01-canon/theorems/THM-2545-word-stratified-hall-arrival-criterion-and-owner-word-transportation-hostile.md): separate margins do not determine target incidence.
* [THM-2549, future-pullback neutrality](../../01-canon/theorems/THM-2549-future-pullback-target-neutrality-and-cemetery-hall-boundary.md): chronology alone does not retain the target label.
* [THM-3044, pointed correspondence and Hall dual](../../01-canon/theorems/THM-3044-pointed-correspondence-determinant-hall-dual-and-cycle-boundary.md): close an actual correspondence before using a loop invariant.
* [THM-2177, unsplittable cost flow](../../01-canon/theorems/THM-2177-planar-counterexample-to-goemans-unsplittable-cost-flow-conjecture.md): fractional marginal feasibility does not imply a compatible unsplit routing.
* [Bernoulli boundary note](bernoulli_boundary_20260925.md): exact endpoints can contain the entire surviving integer even when bulk density tends to zero.

No theorem from LRC is asserted to imply Collatz. The transferred mechanism is to retain the actual incidence table and its capacity cuts. No tournament is constructed.

| Concept | Preserved information | Hostile / next test |
|---|---|---|
| Hub moments | source mass and aggregate service | heavy ancestors can be irrelevant to most sources |
| Bad-ancestor capacity | actual source-to-hub incidence | tilted diagonal remains |
| Bounded slope tube | actual incidence with weights bounded | population might lose its entropy exponent |
| Affine offset | exact integer source interval at fixed slope | audit repeated ancestors and source positivity |
| Hall-style cut | every source must hit a legal modified pair | does not itself assert unmodified descent |
| Signed and multiplier controls | cycle/escape distinction | plus and minus remain indistinguishable asymptotically |

## 1. Exact source cells and their joint incidence

Use the single-halving map

    T(n)=n/2              if n is even,
    T(n)=(3n+1)/2         if n is odd.

For a length-L parity word, let e_k count its odd symbols in its first k positions and put

    w_k=3^e_k/2^k,   w_0=1.

Call the word slope-undecided if `w_k>1` for every `1<=k<=L`. Each parity word corresponds to one source class `a mod 2^L`; this can also be obtained by exact iteration of the representatives `0<=a<2^L`. Its actual prefixes have the affine form

    T^k(n)=w_k n+c_k=w_k(n+h_k),   h_k=c_k/w_k.           (1)

A source class `a+2^L Z` maps at time k to exactly the endpoint class

    T^k(a) + (3^e_k 2^(L-k)) Z.                          (2)

Thus endpoint phase is explicit, not a free choice. A point in (2) has exactly one source in that particular source class. For sufficiently large positive endpoints that source is positive.

For any selected set S of slope-undecided source classes, define its exact bad-ancestor capacity

    Q_S(v)=sum_(a in S) sum_(k<L)
             w_k(a) 1_(v in endpoint class (a,k)).        (3)

It is periodic with period dividing `2^L 3^(L-1)`. If `rho=|S|/2^L`, its spatial mean is exactly

    E Q_S=L rho,                                        (4)

because each summand in (3) has mean `w_k/(3^e_k 2^(L-k))=2^-L`.

Two cells `(k,e,r)` and `(ell,f,s)` intersect if and only if

    r=s mod 2^(L-max(k,ell)) 3^min(e,f).

When they intersect, their contribution to `E Q_S^2` is exactly

    3^min(e,f) / 2^(L+max(k,ell)).                        (5)

This provides an exact joint-law calculation without enumerating the full period.

**Connection contract.** Source = an actual finite parity itinerary and its integer source class; target = its endpoint congruence with an affine weight. Preserved predicate = a modified map must hit one of these endpoint pairs to alter the source's first L steps. Destroyed information if only the law of Q is kept = which hubs lie on one source's actual path. Sidecar = the source residue and time. Decisive test = compare (5) to direct pointwise summation over the entire CRT period.

## 2. Conditioned energy improves constants, but not the exponent by itself

For all slope-undecided classes, replacing the unrestricted THM-4477 hub weight W by Q gives

    upper flip density >= rho_L^2 / E Q^2.               (6)

The same first-hit/logarithmic-weight argument proves (6); Q only charges ancestors that really belong to the selected bad population.

Exact examples:

| L | rho_L | E Q^2 | lower bound from (6) |
|---|---|---|---|
| 8 | 19/256 | 164793/32768 | 361/329586, about 0.001095 |
| 10 | 1/16 | 1005725/131072 | 512/1005725, about 0.000509 |
| 14 | 367/8192 | 945280621/67108864 | 134689/945280621, about 0.0001425 |

At L=8 the inherited unrestricted Cauchy-Schwarz bound is about 0.0000449 and the inherited supremum bound about 0.000747. The improvement is a concrete consequence of retaining actual incidence. It is not a proof of a better asymptotic exponent.

**PROVED hostile to the tempting hypothesis `E Q^2 <= poly(L) rho_L`.** Retain only the diagonal terms of (5) at time `k=L-1`. Every admissible prefix of length L-1 extends by an odd step to a slope-undecided length-L word. Therefore

    E Q^2 >= (1/2) P_(p=3/4)(w_j>1 for 1<=j<=L-1).     (7)

The right side is uniformly at least 1/8. Here is an elementary bound. Under independent parity symbols with odd probability 3/4, `R_j=2^j/3^e_j` is a nonnegative martingale because its expected multiplier is

    (3/4)(2/3)+(1/4)2=1.

The first odd step has probability 3/4 and leaves R=2/3. The maximal inequality bounds the conditional probability of ever reaching R>=1 by 2/3. Hence survival has probability at least `(3/4)(1/3)=1/4`. The same finite-time bound is enough; no limiting stochastic result is required.

Since rho_L decays exponentially, `poly(L)rho_L` tends to zero. The proposed energy hypothesis is false. Exact tilted survival at L=256 is about 0.380051, giving the much larger diagonal lower bound 0.190025. The mechanism is that squaring tilts rare high-slope ancestors toward odd frequency 3/4, above the critical frequency log_3(2).

A stronger retained-incidence dual remains valid:

    H_S = E_(source)[1_S / max_(k<L) Q_S(T^k n)].         (8)

For each endpoint v the total dual load is at most one, because every ancestor's denominator is at least Q_S(v). Thus flip density is at least H_S. Also `H_S >= rho^2/E Q_S^2`: use Cauchy-Schwarz and the exact change-of-variables identity

    E_(source)[1_S sum_(k<L) Q_S(T^k n)]=E Q_S^2.

The corrected exact enumeration runs through L=6. Its endpoint capacity period is `P=2^L 3^(L-1)`, but the joint source observable generally needs the larger period `2^(L-1)P`. A first experiment incorrectly reused P as the source period; those harmonic values were withdrawn before checkpoint. The repaired code also checks the displayed change-of-variables identity over the full source period. The exact values in the linked output are finite results, not asymptotic assertions.

## 3. The affine-offset multiplicity lemma

The key additional coordinate is h in (1). Its update is exact:

    h_(j+1)=h_j                 at an even step,
    h_(j+1)=h_j+1/(3w_j)        at an odd step.           (9)

For a slope-undecided prefix all w_j are at least one. Therefore

    0<=h_k<=k/3.                                        (10)

**PROVED.** Fix a positive integer endpoint v, a time k, and an odd-step count e. Among actual integer ancestors n with `T^k(n)=v`, exactly e odd steps, and every prefix slope at least one, there are at most

    floor(k/3)+1                                        (11)

distinct sources.

*Proof.* Their common slope is `w=3^e/2^k`; (1) and (10) put every source in the same real interval

    v/w-k/3 <= n <= v/w.

An interval of length k/3 contains at most floor(k/3)+1 integers. One source cannot carry two different length-k trajectories under a deterministic map. No independence or marginal estimate is used.

The statement includes k=0. Repetition of a source at different times is permitted; later bounds sum over times, so repetition can only overcount service. The signed map `(3n-1)/2` has `-k/3<=h_k<=0` and the same count. For an odd multiplier q and sign sigma, the interval length is k/q and the count is floor(k/q)+1.

This is a polynomial multiplicity bound for the slope-undecided part of the inverse graph. The entire inverse graph need not have polynomial size: contracting prefixes can make h arbitrarily larger than the bound (10).

## 4. A bounded slope tube gives a uniform hub capacity

Fix K>=1 and retain only words satisfying

    1<w_k<=K for every 1<=k<=L.                          (12)

For a fixed time k there are at most

    floor(log_3 K)+1

possible integer values e, because `1<=3^e/2^k<=K` puts e in an interval of length log_3 K. Consequently every integer hub is visited by at most

    M_L(K)=(floor(log_3 K)+1)
           sum_(0<=k<L)(floor(k/3)+1)                    (13)

distinct selected sources, counting a source again at each possible time. The corresponding weighted capacity satisfies `Q_S(v)<=K M_L(K)`.

Equation (13) controls all off-diagonal collisions as well as the diagonal. In particular,

    E Q_S^2 <= K M_L(K) L rho.                           (14)

The parent coordinator's crucial synthesis was to trim the high-slope tail instead of conditioning on badness alone. The elementary interval (10) then closes the remaining collision bound.

The fixed tube K=2 is a hostile: it contains no length-5 word, since the only integer odd-count allowing w_5>1 gives w_5>=81/32>2. A tube must widen with L to retain the full undecided entropy exponent. Experiments with K=4,8,16,32 confirm (14), including all affine offsets rather than only the moment totals.

## 5. A direct counting cut avoids moments entirely

Let F be the set of flipped pair indices in the THM-4470/4475 pairing family, with pair i equal to `{2i-1,2i}`. Suppose the modified map makes every sufficiently large source descend below itself within L steps.

Every selected source n satisfying (12) must meet a flipped pair among its first L original Collatz vertices: otherwise the modified trajectory agrees with T there, and positive affine carry gives `T^k(n)>n` throughout.

For `n<=X`, every such vertex satisfies

    T^k(n)=w_k(n+h_k)<=K(X+L/3).                         (15)

By (13) each vertex can serve at most M_L(K) selected sources; one flipped pair can serve at most twice that many. Therefore the exact finite inequality is

    # {selected n<=X}
       <= 2 M_L(K) # {i in F: i <= [K(X+L/3)+1]/2}.     (16)

Integer floors may be applied to the upper endpoint. Divide by X and let X grow. Since the selected source set is periodic with density rho, (16) proves

    lower asymptotic flip density >= rho/[K M_L(K)].     (17)

This is stronger in quantifier than merely a bound on upper density. It needs neither logarithmic averaging, Haar measure, second moments, nor a fractional matching. It is a capacity cut on actual source-to-pair incidence. Certificate reuse is allowed; the finite capacity comes from deterministic affine geometry, not an artificial rule that a certificate may be used once.

## 6. The sharp-price consequence and its exact scope

Write

    alpha=log_3 2,  h=H_2(alpha),  eta=1-h.

The coordinator supplies the following elementary population construction, repeated here to make the consequence explicit. Choose a block length b with `b=o(L)` and exactly `e=ceil(alpha b)` odd symbols per block. Then the block's final slope is strictly between 1 and 3. Rotating immediately after a minimum partial log-slope produces a word with every prefix slope greater than one. Each resulting word has at most b preimages under this rotation selection, so there are at least `binomial(b,e)/b` good blocks.

Concatenate `t=floor(L/b)` such blocks and fill the remainder with odd symbols. Every prefix slope is at least one and is bounded by

    K_L=3^(t+1) (3/2)^b=exp(O(L/b+b)).                   (18)

The number of selected words is at least `(binomial(b,e)/b)^t`, hence elementary Stirling bounds give

    rho_L^tube >= 2^(-eta L-O(b+(L/b)log b)).             (19)

Strictness at nonempty prefixes follows from irrationality of log_3 2: no nonempty prefix has slope exactly one. Concatenation is injective because block boundaries are fixed.

Choosing `b` of order `sqrt(L log L)` balances the remainder and block-count costs. Equations (17)--(19) give

    lower flip density >= 2^(-eta L-O(sqrt(L log L))).    (20)

Together with the inherited THM-4475 upper bound `delta_L<=2 rho_L<=2^(1-eta L)`, the consequence is

    -(1/L) log_2 delta_L -> eta=1-H_2(log_3 2).           (21)

This closes the stated sharp exponential-rate target for the price of enforcing uniform bounded-time descent in the pairing family, subject to the coordinator's independent audit and status update. It does **not** prove that the unmodified Collatz map descends. Indeed no fixed L can work for that map because of arbitrarily long initial rising runs.

The result is consistent with THM-4477's distribution-only square-loss theorem: (17) uses actual source intervals at fixed ordered affine slope, information absent from the unrestricted hub-weight law. It is also consistent with the minus sheet: its negative normalized carry has the same interval width, and for any fixed finite family all sufficiently large sources with slopes above one still fail to descend. Thus the price exponent can agree even though the minus map has nontrivial positive cycles.

### Arbitrary edit maps have the same exponential price

Let G be any deterministic map from positive integers to positive integers and let `D={n:G(n)!=T(n)}`. Require every n>1 to descend below itself within L G-steps. Every selected tube source must meet D in its first L T-vertices. The same finite count, now with individual edited vertices instead of pairs, yields

    lower density(D)>=rho/[K M_L(K)].                    (22)

For the upper bound define the actual finite-horizon bad set

    B_L={n>1:T^j(n)>=n for all 1<=j<=L}

and set G(n)=1 on B_L and G(n)=T(n) elsewhere. A source in B_L descends immediately. Any other source follows its original descending path until it either descends or encounters B_L and is sent to 1. Thus G satisfies the required L-step property.

The density of B_L is exactly rho_L. Indeed, each length-L parity class is affine at every prefix. A class with a prefix slope below one descends there above a finite threshold; a class with all prefix slopes above one never descends under the positive carry. No nonempty prefix has slope exactly one. Therefore B_L differs from the periodic slope-bad set by only finitely many integers.

Consequently, if epsilon_L is the infimum upper edit density over these arbitrary G, then

    2^(-eta L-O(sqrt(L log L))) <= epsilon_L
                                  <=rho_L<=2^(-eta L).   (23)

The rate is therefore the same for arbitrary edits and for the more constrained pairing family. The shortcut T itself has no uniform finite L; a zero edit set cannot satisfy this property. Neither the existence of arbitrarily sparse successful edits nor their exact rate proves convergence of T.

## 7. Parameter form: affine tubes beyond the Collatz alphabet

The mechanism is more general than its entropy application. Let `T_(q,sigma)` use n/2 on even n and `(qn+sigma)/2` on odd n, where q>=3 and sigma are odd integers. Consider any selected length-L actual source paths whose slopes satisfy

    0<A<=w_j=q^e_j/2^j<=K,    0<=j<=L,    A<=1<=K.

Their offsets have one sign and obey

    |h_k|<=k |sigma|/(q A).

At fixed time k the exponent e has at most `floor(log_q(K/A))+1` possible values, and each fixed slope has at most `floor(k |sigma|/(q A))+1` integer ancestors at any hub. Thus replace (13) by

    M_(q,sigma,A,K,L)
      =(floor(log_q(K/A))+1)
        sum_(k<L)(floor(k |sigma|/(q A))+1).              (24)

Any edit set required to hit every selected path satisfies the same density cut `density(D)>=rho/(K M)`, provided the selected source set has density rho; finitely many small sources that leave the positive domain may be omitted. The endpoint range is now bounded above by `K(X+L |sigma|/(q A))`.

For A<1 this is only a **hitting-set theorem**: membership in the tube no longer guarantees lack of numerical descent. The reason the edit set must hit these paths must be supplied separately. No drift conclusion follows merely from (24). The q=5 control retains its positive cycles, and the theorem has no automatic population estimate for a new multiplier.

The lower slope bound is load-bearing. For q=3 and sigma=1 the three sources `4,5,6` all reach hub 2 in exactly five shortcut steps with exactly two odd steps. Their prefix slopes stay above 1/16, but the A=1 multiplicity bound would give only `floor(5/3)+1=2`. The repaired bound (24) permits these three. The failure is caused by accumulated normalized offset during the contracting prefixes, exactly the coordinate A controls.

An even broader finite-alphabet version needs no prime powers. Suppose a deterministic integer map has d positive slopes `a_1,...,a_d` and finitely many affine branches `a_i n+b_i`, with `|b_i/a_i|<=B`. Along a path with every product slope at least A, the normalized offset lies in `[-kB/A,kB/A]`. At time k there are at most `binomial(k+d-1,d-1)` possible products of slopes. Hence a valid uniform capacity is

    sum_(k<L) binomial(k+d-1,d-1)(floor(2kB/A)+1).       (25)

If all intercepts have the same sign, replace 2kB/A by kB/A. This is polynomial in L for a fixed branch alphabet and fixed A,B. A slope cap K again converts it to an edit-density bound. Allowed branch words, actual integer guards, source population density, and the need to hit a path remain explicit inputs; an arbitrary affine alphabet does not inherit Collatz's parity-word count.

## 8. Fixed-slope injectivity is REFUTED; the carry interval survives

**REFUTED candidate, never used in the price proof:** at fixed time k and odd count e, the slope-undecided integer source-to-hub map need not be injective. The temptation was to replace `floor(k/3)+1` in (13) by one. The counterexample below retains two actual positive sources, every slope prefix, the endpoint, and both ordered carries. Minimality is not claimed.

The candidate first survived two finite probes: all sources through 1000000 and times through 40 (2388567 surviving prefixes), and all symbolic slope-undecided words through length 24 (654279 words, covering all integer heights at those depths). A structured carry decoder then found a collision at

    k=233,  e=153,  source difference=4.

### Exact positive source construction

Let u be the least positive inverse of `3^153` modulo `2^80`, and set

    u=186937257649965781719819,
    N=2^153 u-1,
    M=N-4,
    v=(3^153 u-1)/2^80.

Both sources are positive odd integers. Explicitly,

    M=2134446157293545680116155684904813624866473205075599455753325089128443,
    N=2134446157293545680116155684904813624866473205075599455753325089128447,
    v=1544714368800273392784179563875656116569530364984522231576309181404945461.

Their actual shortcut parity words, where 1 means an odd step, are

    N: 1^153 0^80,
    M: (110)^51 t,
    t=01110111100111110111011110110011111111110010011110000100101111010011001101010101.

The tail t has length 80 and contains 51 ones, so both complete words have length 233 and contain 153 ones. Every nonempty prefix in both words satisfies `3^e_j>2^j`. The actual integer trajectories both finish at v; this is explicitly replayed, rather than inferred solely from an affine equation.

Put `A=3^153`. The ordered carries are

    C_N=A-2^153,
    C_M=5A-2^153,
    C_M-C_N=4A.

Thus `A(N-4)+C_M=A N+C_N`, exactly explaining the merger. The normalized offsets differ by the integer 4. The safe interval capacity remains valid; the failed implication was that two admissible offsets could never differ by an integer source separation.

The construction also recovers earlier negative-center structure: `N=-1 mod2^153` initially follows the all-odd center -1, whereas `M=-5 mod2^153` follows the negative shortcut cycle `-5 -> -7 -> -10 -> -5`, with word 110 for 51 copies. The common parity-prefix lengths with these respective negative trajectories are exactly 153: each positive word next takes an even step where its negative reference takes an odd step. Their normalized final offsets are exactly `1-(2/3)^153` and `5-(2/3)^153`. This is a precise shared-prefix and carry statement, not a claim that either positive source lies on a negative orbit.

### How the hostile was found

For odd positions `0=p_1<p_2<...<p_e`, the carry is

    C=sum_(i=1)^e 3^(e-i) 2^p_i.

Given C and e, a candidate position list is decoded greedily: take the two-adic valuation of the remaining carry as p_i, subtract `3^(e-i)2^p_i`, and continue. Retain only increasing positions with `2^p_i<3^(i-1)` for i>=2, zero final remainder, and final slope above one. The bounded family `C=(d+1)3^e-2^e`, `3<=e<=200`, and `d=4,8,...<=e/3` compares against the dense-odd carry `3^e-2^e`. It produced the displayed e=153,d=4 witness without a huge source enumeration.

### Strongest surviving injectivity statement: all heights through time 31

Every slope-undecided source for k>=2 has its first two symbols odd, so its source class is `3 mod4`. At fixed k,e the minimum normalized carry is

    h_min=1-(2/3)^e,

attained by putting all odd symbols first. The largest legal positions are

    p_1^max=0,
    p_i^max=min(ceil((i-1)log_2 3)-1, k-e+i-1),
    h_max=sum_(i=1)^e 2^(p_i^max)/3^i.

These bounds are simultaneous: the position ceilings increase strictly, satisfy every prefix guard, and leave room for the remaining odd symbols. Hence they give the exact offset width for fixed k,e.

For all 199 possible `(k,e)` pairs with `1<=k<=31` and `3^e>2^k`, the width is strictly below 4. The maximum is

    13805179460/3486784401 = 3.9592868...,

at k=31,e=20. Two distinct sources in the same `3 mod4` class differ by at least 4, so no collision is possible through time 31 at any height. The case k=1 has only one word. At k=32,e=21 the width becomes `43561973452/10460353203>4`; this is an opening of the interval bound, not a collision assertion.

A separate obstruction explains why one-position experiments miss the phenomenon. Moving just the i-th odd position by d changes C by `3^(e-i)2^p(2^d-1)`. Divisibility by `3^e` requires d even and `d>=2*3^(i-1)`. The slope guard instead gives `d<(i-1)(log_2 3-1)`, an impossibility. Thus a collision needs coordinated changes of several odd positions, as in the explicit witness.

Reproduce with [crossroads_20260926_flow_collision.py](../../04-computation/experiments/crossroads_20260926_flow_collision.py); exact [normal output](crossroads_20260926_flow_collision.out) and [optimized output](crossroads_20260926_flow_collision_optimized.out). Both runs pass and are byte-identical. Validation uses explicit checks that remain enabled under `python3 -O`, including both actual parity words, every strict prefix slope, final endpoint, odd counts, carry difference, and the modular-inverse source construction. Independently, the root coordinator reconstructed the least source class as `(-C_N * (3^153)^(-1)) mod2^233`, obtained the same N, and directly replayed both positive shortcut trajectories with all prefix guards and both carries. The conjectured injectivity is retired; the polynomial interval capacity and sharp-price theorem are unchanged.

## 9. Controls, scope, and reproduction

The script [crossroads_20260926_flow.py](../../04-computation/experiments/crossroads_20260926_flow.py) uses Python integers and Fraction only. Output: [crossroads_20260926_flow.out](crossroads_20260926_flow.out).

    python3 04-computation/experiments/crossroads_20260926_flow.py

Universes:

* full exact source-cell energy through L=14 for `3n+1`, through L=10 for `3n-1` and `5n+1`;
* pointwise verification of cylinder-intersection energy on every CRT residue through L=4 for all three maps;
* exact full-period source/hub harmonic dual through L=6, including the independent change-of-variables check;
* an independent inverse-tree implementation for hubs 1..128 and depths 1..18 on all three maps, checking the affine-offset interval and fixed-(k,e) multiplicity;
* a separate forward implementation on sources 1..8192, depths through 18, and three signed/multiplier parameter settings, checking 209730 source groups with lower slope A<1 and retaining the explicit `4,5,6` hostile above;
* exact Bernoulli(3/4) tilted survival through L=256;
* slope caps 2,4,8,16,32 and L=6,8,10,12,14, retaining every exact source offset.

The plus and minus conditioned energy tables agree through the tested range, as sign conjugacy predicts. The q=5 normalized unrestricted conditioned energy `E Q^2/(L^2 rho)` rises from 1 at L=1 to about 8.42 at L=10; q=3 reaches about 1.60 at L=14. The exact odd cycles `5<->7` and `13->33->83->13` are checked separately. Their affine fixed points are retained; residue self-loops alone are not interpreted as positive integer cycles.

The initial untrimmed-energy hypothesis failed and its repaired form (14) retains an explicit tube parameter. The remaining Collatz question is still pointwise termination of the unmodified map; no averaging, price limit, or small-density modification identifies that truth value.
