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

The coordinator supplies the following elementary population construction, repeated here to make the consequence explicit. Choose block length `b=floor(sqrt L)` and exactly `e=ceil(alpha b)` odd symbols per block. Then the block's final slope is strictly between 1 and 3. Rotating immediately after a minimum partial log-slope produces a word with every prefix slope greater than one. Each resulting word has at most b preimages under this rotation selection, so there are at least `binomial(b,e)/b` good blocks.

Concatenate `t=floor(L/b)` such blocks and fill the remainder with odd symbols. Every prefix slope is at least one and is bounded by

    K_L=3^(t+1) (3/2)^b=exp(O(sqrt L)).                  (18)

The number of selected words is at least `(binomial(b,e)/b)^t`, hence elementary Stirling bounds give

    rho_L^tube >= 2^(-eta L-O(sqrt L log L)).             (19)

Strictness at nonempty prefixes follows from irrationality of log_3 2: no nonempty prefix has slope exactly one. Concatenation is injective because block boundaries are fixed.

Since `K_L M_L(K_L)=exp(O(sqrt L)) poly(L)`, (17) and (19) give

    lower flip density >= 2^(-eta L-O(sqrt L log L)).     (20)

Together with the inherited THM-4475 upper bound `delta_L<=2 rho_L<=2^(1-eta L)`, the consequence is

    -(1/L) log_2 delta_L -> eta=1-H_2(log_3 2).           (21)

This closes the stated sharp exponential-rate target for the price of enforcing uniform bounded-time descent in the pairing family, subject to the coordinator's independent audit and status update. It does **not** prove that the unmodified Collatz map descends. Indeed no fixed L can work for that map because of arbitrarily long initial rising runs.

The result is consistent with THM-4477's distribution-only square-loss theorem: (17) uses actual source intervals at fixed ordered affine slope, information absent from the unrestricted hub-weight law. It is also consistent with the minus sheet: its negative normalized carry has the same interval width, and for any fixed finite family all sufficiently large sources with slopes above one still fail to descend. Thus the price exponent can agree even though the minus map has nontrivial positive cycles.

## 7. Controls, scope, and reproduction

The script [crossroads_20260926_flow.py](../../04-computation/experiments/crossroads_20260926_flow.py) uses Python integers and Fraction only. Output: [crossroads_20260926_flow.out](crossroads_20260926_flow.out).

    python3 04-computation/experiments/crossroads_20260926_flow.py

Universes:

* full exact source-cell energy through L=14 for `3n+1`, through L=10 for `3n-1` and `5n+1`;
* pointwise verification of cylinder-intersection energy on every CRT residue through L=4 for all three maps;
* exact full-period source/hub harmonic dual through L=6, including the independent change-of-variables check;
* an independent inverse-tree implementation for hubs 1..128 and depths 1..18 on all three maps, checking the affine-offset interval and fixed-(k,e) multiplicity;
* exact Bernoulli(3/4) tilted survival through L=256;
* slope caps 2,4,8,16,32 and L=6,8,10,12,14, retaining every exact source offset.

The plus and minus conditioned energy tables agree through the tested range, as sign conjugacy predicts. The q=5 normalized unrestricted conditioned energy `E Q^2/(L^2 rho)` rises from 1 at L=1 to about 8.42 at L=10; q=3 reaches about 1.60 at L=14. The exact odd cycles `5<->7` and `13->33->83->13` are checked separately. Their affine fixed points are retained; residue self-loops alone are not interpreted as positive integer cycles.

The initial untrimmed-energy hypothesis failed and its repaired form (14) retains an explicit tube parameter. The remaining Collatz question is still pointwise termination of the unmodified map; no averaging, price limit, or small-density modification identifies that truth value.
