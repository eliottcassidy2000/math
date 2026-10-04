# Binary debt and ternary sibling addresses in one Collatz selector

2026-10-04. **PROVED:** the scoped identities, ternary address isometry,
natural-density bounds and product law below. **FINITE-EXACT:** the bounded
word search and executable controls. **OPEN:** universal grounded coverage.
No literature-priority claim or new global convergence result.

This integrates the concurrent commits `9b320b7701` and `7c15ea9559` with
the [inverse shields and family lifts](collatz_join_shields_and_lifts_20261004.md).
The useful synthesis is two independent coordinates with different jobs:
binary digits select exact forward words; ternary digits select a smaller
common-future source. Neither coordinate replaces the integer or its
certificate. Their guards combine by CRT even though the underlying
affine Collatz operations do not commute.

## 1. Inheritance and the three lanes

Closest proved mechanism: the incoming
[debt carry decoder](debt_word_join_search_20261004.md), which retains
the alternative word's length and exact carry. Canonical hostiles: 7, 27
and703 still fail both selected banks. Corrected near miss: topological
address coverage is not integer coverage, as the incoming
[joint-address construction](joint_home_address_20261004.md) already shows.
Least-used sidecar: a real size bound on the possible ternary addresses.

Anchor: extend guarded smaller-child rewrites to debt eight. Niche: compute
overlap between the binary and ternary banks. Wildcard: recognize the
sibling index itself as a ternary coordinate. The live board is
**original source / ordered carry / binary guard / ternary address /
size bound / supplied child certificate**. This gives a typed connection
to the repository's address constructions, with the missing height bound
made explicit. No new general meta-pattern is promoted.

## 2. Eight additional debt families

Use the incoming notation

    Y=2^s*3^e*M+1,   s in{1,2}, e>=1, M positive odd.

For a positive valuation word u of length l, cost S and carry B,
the unique candidate word v with the same affine endpoint on M must have

    length(v)=l+e, cost(v)=S-s, carry(v)=(3^l+B)/2^s.

The decoder and exact-cylinder proof are inherited, not new here. The new
search extends the declared debt range from four to eight. Besides the
eight incoming rows, it finds these eight. The last exponent a can be
any positive integer.

| s | e | u on Y | v on M | Least source cost at a=1 |
|---:|---:|---|---|---:|
| 2 | 5 | (5,10,a) | (1,1,1,1,1,1,5,a+2) | 16 |
| 1 | 5 | (1,1,5,11,a) | (2,1,1,2,1,1,2,1,2,a+4) | 19 |
| 2 | 6 | (3,12,3,a) | (1,2,1,1,1,1,2,1,2,a+4) | 19 |
| 2 | 7 | (13,5,a) | (1,1,1,1,2,1,2,1,2,a+4) | 19 |
| 1 | 6 | (1,1,17,a) | (1,4,1,1,1,1,2,1,2,a+4) | 20 |
| 1 | 7 | (1,5,12,3,a) | (1,3,1,1,1,1,2,1,2,1,2,a+4) | 22 |
| 2 | 8 | (5,17,a) | (3,2,1,1,1,2,1,2,1,2,a+4) | 23 |
| 1 | 8 | (1,11,11,a) | (2,3,1,2,1,2,1,1,2,1,2,a+4) | 24 |

**PROVED all-height transport.** For every row and r>=1, the source word
and child word are

    (1^r, 2^e, u),   (1^(r-1), 2e+s, v),

where 2^e in the word means e copies of the letter two. On its exact
source cylinder the child is (n-1)/2 and the endpoints coincide. This
follows by the incoming affine comparison and exact final-oddness guard.
The two lengths are equal; source cost exceeds child cost by one.
Increasing both final exponents by a-1 preserves the identity.

The displayed least costs are **FINITE-EXACT within this word model**.
Every positive composition of every cost below24 is examined; cost24
stops at the first witness for the last missing pair. Guards discard
impossible first letters, and the length/cost inequality discards impossible
decodings. These are necessary filters, not an assumed valuation cutoff.
This does not minimize all possible certificate representations and does
not prove that a row exists for every debt e.

The source needs a supplied child certificate. These are affine common-
future transports; no additional assertion about preservation of first-hit
certificate length is needed here. An implementation can use the earlier
rooted suffix if either route encounters the root before its displayed end.

## 3. Sibling height is an exact ternary address

Let S(n)=4n+1 and let ell(k) be the least positive ell for which

    3^ell > 2^(ell+2k),   k>=0.

The previous note proves that the class

    n = R(k) mod3^ell(k),
    R(k)=-(4^k+2)/(3*4^k),                            (1)

gives a positive odd smaller source

    m=2^ell(k)*(S^k(n)+1)/3^ell(k)-1 < n.

It takes ell(k) exponent-one steps to S^k(n), which has the same next
odd iterate as n. In (1), the rational denominator after cancelling3 is
a unit modulo every power of3, so R(k) is a well-defined ternary integer.

**PROVED isometry and complete finite addresses.** For h>k,

    R(h)-R(k) = (2/3)*4^(-h)*(4^(h-k)-1),
    v3(R(h)-R(k)) = v3(h-k).                         (2)

For completeness, write h-k=3^b*d with 3 not dividing d.
The binomial expansion gives v3(4^d-1)=1. If x=1 mod3,
then x^2+x+1 has valuation exactly1, so cubing x adds exactly one
to v3(x-1). Induction gives v3(4^(h-k)-1)=b+1 and proves (2).
Thus k mod3^a maps bijectively to R(k) mod3^a for every a. Compatible
residues extend R to a distance-preserving bijection of the ternary integers.

Equivalently, a supplied n has a unique ternary address K(n) satisfying

    4^K(n) = -2/(3n+1),

where the power is defined by those compatible finite residues. The
certificate guard is

    k = K(n) mod3^ell(k).                            (3)

No real or complex logarithm is intended. This is a reparametrization of
the exact divisibility guard, not a new termination argument.

**Exact overlap.** For j<k, the k guard is contained in the j guard iff

    k = j mod3^ell(j).

Otherwise they are disjoint. This follows from (2) and the increasing
depths ell(k). Remove contained guards before adding their densities.
For k<=20 the retained indices are

    0,1,2,4,5,7,8,10,11,13,14,16,17,19,20.

This is a proved initial list, not a claim that excluding multiples of3
is the entire infinite overlap rule. For example later k=82 is contained
in k=1 because82=1 mod81.

## 4. The size sidecar turns ternary measure into natural density

The union of these guards is open and dense in the ternary integers.
Indeed, in any ball modulo3^a choose its unique address k in[0,3^a).
If ell(k)>=a, its guard lies inside that ball; if ell(k)<a, its guard
contains the ball. This statement alone gives no density conclusion on
ordinary integers.

The missing ingredient is positivity. For an applicable positive odd n,

    t=(S^k(n)+1)/3^ell(k)

is a positive even integer. Hence

    2^(ell(k)+1)-1 <= m < n,
    ell(k)+1 < log2(n+1),                            (4)

and ell(k)>2k. There are only O(log X) eligible k for n<=X. This is
an actual bound on candidate addresses, not a heuristic sparsity claim.

**PROVED natural density.** Let P be the indices not contained in an
earlier guard. Among positive odd integers the union has natural density

    delta = sum_(k in P) 3^(-ell(k)).                 (5)

For fixed K this is a finite disjoint union of residue classes. Above K,
there are only O(log X) relevant classes by (4); their count through X
is at most `(X/2)*sum_(k>K)3^(-ell(k))+O(log X)`.
The depths satisfy ell(k+1)>=ell(k)+3: minimality at k implies
`(3/2)^(ell(k)+2)<= (27/8)*4^k < 4^(k+1)`.
Thus the remaining series has a geometric tail, proving (5) by squeezing.

In particular delta<=9/26<1. Retaining the primitive guards through k=20
and bounding the rest gives an interval of width

    1/21694014376608093878068290546077358.

Exact rational endpoints are stored in the JSON. Numerically delta is
approximately0.34613647. The successful guard union is therefore both
ternary-dense and of natural density strictly below one. This is a precise
instance of the address/coverage distinction, with the real size sidecar
making the integer statement possible.

## 5. Binary and ternary exclusions multiply

If there were a positive counterexample, let n be the least one. Immediate
descent excludes first valuation at least2. The inherited unbounded reset
switch excludes a run of ones ending at a reset at least3. Their remaining
necessary domain is first-reset-two, of odd-relative density1/4.

The sixteen selected debt families are disjoint within that domain.
For a given row let

    c=2e+sum(u without its last letter)+1.

At fixed r, union over every terminal a is one exact binary prefix class
of modulus2^(r+c); summing r>=1 gives odd-relative density2^(1-c).
Distinct r have different initial runs of ones, distinct e different
following runs of twos, and s is distinguished by the next valuation.
Thus the remaining binary necessary domain has exact density

    D2 = 1/4 - sum_(sixteen rows) 2^(1-c).             (6)

**PROVED product law for these specified tests.** The sibling guards
depend only on powers of3. CRT makes every finite binary/ternary pair
independent among odd integers. The binary run tails tend to zero, and
the ternary tails have the uniform size bound from section4. Finite
truncation followed by these tail bounds therefore gives

    density of the combined necessary domain = D2*(1-delta).       (7)

The exact value is `D2=134156760173/549755813888`; (7) is approximately
`0.15956213718786713`. Exact rational interval endpoints are saved in the JSON.
This is a necessary domain for the **least** counterexample under these
specified rules. It is not the density of counterexamples, the unresolved
set under every repository rule, or a bound on stopping times. The eight
new rows make only a small additional change in density; their main value
is extending the checked mechanism to further debt levels.

In particular, the combined single-stage tests leave positive density.
They cannot exhaust all positive inputs by direct membership alone.
Iterating through actual forward frontiers is a different proposal and
still needs the original-source comparison; a smaller child of a later,
larger frontier need not be smaller than the original source.

One finite control connects this to the frozen239-source benchmark:
the incoming eight rows select exactly6783, with child3391 and join6113.
The sibling family at144615 supplies a separate route while the incoming
eight-row selector does not match it. These show distinct applicability;
they do not claim a new online learner's final seed count.

## 6. What changes in the proof strategy

Store an integer together with its exact binary debt state, the finite
ternary address needed for a candidate k, its original-source size bound,
and all available child certificates. The source projection remains the
same integer. The extra coordinates expose both valid rewrites and why a
particular attempted rewrite cannot apply.

The order justified here is: reuse a grounded shared suffix; test the
binary carry-decoded rules; enumerate only the finitely feasible ternary
addresses from (4); apply the inverse shield before further reverse search;
then advance the actual forward frontier while retaining all observations.
These are alternative proof routes and must not discard one another.

The precise open problem is still a well-founded rule for the frontier
left after these tests. A dense address bank or a larger finite debt table
does not supply that rule. The useful new coordinate is the **address
together with the height it costs**: either coordinate alone loses the
applicability condition. Hostiles7,27 and703 remain explicit next inputs.

Connection contract: source = checked sibling/debt diagrams; target =
mixed guarded certificate states; map = exact source identity plus the
two residue coordinates; preserved predicate = common future with a
supplied smaller-child proof; lost by quotienting = carry, height and
child identity; required sidecars = (4), ordered words and grounded suffix;
cheapest decisive tests = literal replay, exact CRT counts and the retained
hostiles. The operations themselves have not been shown commutative.

## Reproduction

Run `python3 -B 04-computation/experiments/collatz_binary_ternary_guard_fusion_20261004.py`
and repeat with `python3 -O -B`. The
[script](../../04-computation/experiments/collatz_binary_ternary_guard_fusion_20261004.py),
[JSON](collatz_binary_ternary_guard_fusion_20261004.json) and
[stdout](collatz_binary_ternary_guard_fusion_20261004.out) retain the full
search universe, sixteen words, guards, density fractions and controls.
No assertion statements disappear under optimization. The script performs8,765,636 explicit checks, examines12,580,862 words,
and independently checks585 finite CRT combinations. The new parameterized
selector is tested on every displayed family instance and compared with
the incoming selector on all10000 odd sources through19999. The lower-cost
search is exhaustive; the isometry and infinite-density statements follow
from the proofs above, not extrapolation from a finite census.
