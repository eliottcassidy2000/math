# Degree-local tournament count channels

Status: PROVED channel extension and FINITE-EXACT executable checks. This packet newly restores reproducibility of the Session21 channel statements. It neither recovers lost source bytes nor claims old timings or a new self-converse formula. No performance improvement beyond the inherited degree-local kernel is claimed or benchmarked here.

## Definition, carrier and weight

For an odd partition lambda of n, let k_l be its multiplicities, K=sum k_l,
z_lambda=product l^k_l k_l!, and

    E(lambda)=sum_l k_l(l-1)/2 + sum_l l binomial(k_l,2)
              +sum_(l<j) gcd(l,j) k_l k_j.
    F_n(a,b)=sum_(lambda odd, |lambda|=n) 2^(a E(lambda)+b K)/z_lambda.

The API accepts literal integers a>=1, b>=0. The raw unrestricted channel counts simultaneous-isomorphism classes of a ordered tournaments on the same n vertices, with each vertex carrying one of 2^b colors. A permutation fixing this object must have odd cycle lengths because it fixes the first tournament. Its E unordered-pair orbits independently carry a orientation bits, and each of its K vertex cycles carries b color bits. This proves the formula by Burnside. In particular F_n(1,0)=A_n, the usual unlabelled tournament count. The raw channel is an actual integer orbit count; no claim is made that its summands are integers.

The optional no_fixed flag restricts the sum to partitions without 1-parts; it is not a count of tournaments having no fixed vertices or trivial automorphism group. The optional core flag substitutes x_l -> x_l/(1+x_l) in this cycle-index expression before x_l=z^l. Its coefficients may be signed. Only the inherited a=1,b=0 interpretation has been established as pair-free-core enumeration elsewhere; arbitrary channel core coefficients are presented solely as exact substituted series, not newly identified object counts.

K counts permutation cycles, not tournament directed cycles. Neither b nor E measures a tournament's cycle excess, Hamilton paths, strong components or metagraph degree. Aggregation loses individual partitions and individual objects; it preserves the specified weighted Burnside totals and future response action.

## Degree-local transfer proof

Freeze the published incoming kernel as incoming_local_scale.py. At descending odd l, a state has used degree d, rational weight w, future gcd-response c at l, and integer W=O(d)w, where O(d) is the odd part of d!. Add k actual cycles of length l and optional substitution excess e, so q=d+l(k+e). The new exponent contribution is

    a [k(l-1)/2 + l binomial(k,2) + k c] + b k.

Thus the transfer uses the same integer

    T=O(q)/(O(d) l^k O(k))

and the same excess multiplier (-1)^e binomial(k+e-1,e), with e=0 required when k=0. The left-shift exponent is

    a [k(l-1)/2 + l binomial(k,2) + k c] - v2(k!) + b k.

It is nonnegative because a>=1, b>=0 and v2(k!)<=binomial(k,2). Integrality of T follows by taking odd parts of the integer

    q!/(d! l^k k!)
    = binomial(q,lk) [(lk)!/(l^k k!)] (q-lk)!/d!.

All denominators and response updates are unchanged. Only the shift is modified: the response charge is multiplied by a, not b, and only actual cycles k receive b, not the substitution excess. The packed response remains an unweighted gcd sum; scaling it by a would double-count the charge. The incoming remaining-degree mask and carry-safe packing width are unchanged. Merging weights at a common degree remains valid, and terminal division is by O(d).

## Verified self-converse reduction, with source correction

The source is 01-canon/theorems/THM-283-sc-burnside-formula.md, whose statement sums over cycle lengths 2 mod 4, plus exactly one fixed point for odd n. Its numeric prefix n=2..19 is tested here. The source's explanatory sentence claiming that allowed cycles have only even-length pair orbits is incorrect: a transposition fixes its unordered antipodal pair. The actual criterion is the parity of pair-orbit length together with endpoint swap, as also described in HYP-2078.

For an anti-automorphism g, a pair orbit of length L is consistent precisely when the endpoint-swap parity of g^L equals L mod 2. Odd vertex cycles of length >1 fail on their internal pairs. Two fixed vertices fail on their mutual pair. A vertex cycle of length divisible by 4 fails on its antipodal pair: L is even but the endpoints swap. A cycle of length 2l with l odd passes: nonantipodal internal pair orbits have even length and no endpoint swap, while its antipodal pair orbit has odd length l and swaps endpoints. Between allowed even cycles the orbit lengths are even without swaps; edges from the optional fixed point to an even cycle also have even orbit length. This proves necessity and sufficiency.

Write the allowed cycle type as 2lambda plus epsilon fixed points, epsilon in {0,1}, with lambda an odd partition of m. Its pair-orbit count is

    Q=sum_l l k_l + sum_l 2l binomial(k_l,2)
      +sum_(l<j) 2 gcd(l,j) k_l k_j + epsilon K
     =2E(lambda)+(1+epsilon)K.

Its permutation centralizer has size 2^K z_lambda, cancelling the extra 2^K in 2^Q. Twisted Burnside gives

    C_(2m)=F_m(2,0),    C_(2m+1)=F_m(2,1).

Twisted Burnside counts isomorphism classes fixed by converse: within a self-converse class each labelled representative has an anti-isomorphism coset of size |Aut|, giving total n! across the class; a non-self-converse class contributes zero. This also justifies C_0=C_1=1 without exceptions.

## Merged and strong channels

Converse acts as an involution on ordinary classes, so M_n=(A_n+C_n)/2. Let S_n count nonempty strong classes and U_n count strong self-converse classes, with the singleton counted strong and S_0=U_0=0. Every tournament has a unique linearly ordered sequence of strong components, hence as ordinary generating functions

    A(z)=1/(1-S(z)).

Converse reverses that component sequence and replaces each component by its converse. A fixed sequence consists of an arbitrary left sequence, its forced converse right sequence, and either no center or one strong self-converse center. Therefore

    C(z)=A(z^2)(1+U(z)).

These are ordinary generating functions of isomorphism classes, not labelled exponential generating functions. The API extracts S recursively from A*S=A-1 and U recursively from C=A(z^2)(1+U). Strong converse-merged counts are (S_n+U_n)/2. No cancellation assumptions or recognition of individual strong tournaments is needed for these identities.

## API and reproduction

- channel_dp(bound,a=1,b=0,core=False,no_fixed=False) returns the exact Fraction coefficient vector and inherited row/state/transition/transfer telemetry.
- counts(bound) returns integer vectors ordinary, self_converse, merged, strong, strong_self_converse and strong_merged.
- python check.py writes receipt.json; python -O check.py --receipt optimized-receipt.json repeats assertion-independent controls.

The controls compare five channels (1,0),(1,2),(2,0),(2,1),(3,2), all four core/no-fixed combinations, and every bound 0..20 against a separately unmerged rational odd-partition enumerator. The (1,0) results and full telemetry must equal the frozen incoming kernel. All permutations acting on unordered edges independently check ordinary and anti-Burnside counts through n=7 without the odd-partition reduction. A literal labelled-tournament orbit census through n=6 verifies every output channel, with strong connectivity computed by reachability and self-converseness by complement membership in the relabelling orbit. The entire source prefix C_0..C_19 is checked, invalid input rejection is checked, and all six count vectors through n=40 are included. Both normal and python -O runs execute the same explicit exception-based checks.

Hashes pin the inherited theorem, the related monodromy explanation, the recovered Session21 statement, and the exact incoming kernel. Those sources are provenance, not substitutes for the proofs and fresh tests above. The raw multi-channel carrier interpretation, the corrected monodromy proof, and all code are supplied in this packet; no claim of novelty for the classical counts is made.
