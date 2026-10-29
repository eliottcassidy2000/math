# Pair-color palettes: partial comparisons change the symmetry-pole grammar

Status: PROVED below; finite exact exhaustive controls are separate.
Scope: unlabeled complete pair-colorings with endpoint-reversal involution.
No originality or continuous probabilistic-palette claim. This is a bounded
extension of `ALL_POLES_CORE_EXPANSION.md`, not a claim that all partial
metagraph observers are determined by these aggregate counts.

## 1. Objects and exact enumerator

Let C be a finite set of q>=2 colors with an involution τ. Write
`t=|Fix(τ)|`; necessarily 0<=t<=q and q-t is even. A structure on n labeled
vertices assigns a color c_ij to each i<j, with the reverse interpretation
c_ji=τ(c_ij). Every pair is colored; a neutral/undecided color, if present,
is an ordinary τ-fixed color. There are q^N labeled structures, N=n(n-1)/2.
Relabeling applies τ exactly when the order of a pair's endpoints reverses.
Let a_n(q,t) count S_n orbits, and define

    b_n(q,t)=n!a_n(q,t)/q^N,   B_(q,t)(z)=Σ b_n(q,t)z^n.

This depends only on (q,t), because any two involutions with those parameters
are conjugate as color permutations. As before, b_n=E|Aut| for the uniform
LABELED distribution, not the uniform-unlabeled distribution.

For a permutation cycle type λ=(l_1,...,l_c), put

    u(λ)=Σ_i floor((l_i-1)/2)+Σ_(i<j)gcd(l_i,l_j),
    v(λ)=#{i:l_i is even},
    z_λ=product_d d^(multiplicity_d)(multiplicity_d)!.

Then its exact number of fixed colorings is

    F_λ(q,t)=q^u(λ) t^v(λ),
    a_n(q,t)=Σ_(λ partitions n) q^u(λ)t^v(λ)/z_λ.         (1)

Use 0^0=1, so t=0 retains precisely all-odd types.

Proof: follow an unordered-pair orbit under the permutation. If returning
to its initial unordered pair preserves its ordered endpoints, it has q
color choices. If it reverses the endpoints, it has t choices. Within one
l-cycle there are floor((l-1)/2) nonreversing orbits and one reversing
antipodal orbit exactly when l is even. Between two different vertex cycles
there are gcd(l_i,l_j) orbits, none reversing, because a vertex cannot leave
its own cycle. Multiplying choices and applying Burnside proves (1).

## 2. Exact rational core expansion

A core μ has all parts >=2, omitting every fixed vertex. Let m=|μ|,
c=length(μ), a=m-c and

    C_μ(q,t)=(m!/z_μ) F_μ(q,t)/q^binomial(m,2).

The empty core has C=1. Adding f fixed vertices multiplies F by
q^[cf+binomial(f,2)] and the centralizer by f!. Therefore

    B_(q,t)(z)=Σ_μ C_μ(q,t) z^m/(1-q^(-a)z)^(m+1).       (2)

The fixed-point tail calculation is the inherited THM-514 mechanism, now
with signed pair-orbit counts. Coefficient extraction is exact and finite:

    b_n(q,t)=Σ_(μ:m<=n) C_μ(q,t)binomial(n,m)q^[-a(n-m)]. (3)

## 3. Normal convergence: the possible reversing fixed pairs cause no failure

For a core permutation on m vertices there are no fixed vertices. An
unordered pair fixed setwise must be exactly one of its two-cycles; hence
there are r<=m/2 such pairs. All other unordered-pair orbits have size >=2.
Consequently the TOTAL number of pair orbits, reversing or not, is at most

    r+[binomial(m,2)-r]/2 <= m²/4.

Since 0<=t<=q, each orbit has at most q choices, and thus

    F_μ(q,t)/q^binomial(m,2) <= q^[-m(m-2)/4].            (4)

Also c<=m/2, so a>=m/2. On |z|<=R, for large m,
|1-q^(-a)z|^(-m-1) <= exp[2(m+1)R q^(-m/2)],
which tends to 1 and is bounded. Summing m!/z_μ over size-m cores counts a
subset of permutations of m vertices, and is at most m!. The absolute
size-m contribution to (2) is therefore bounded by a fixed constant times

    m! R^m q^[-m(m-2)/4].                               (5)

The ratio of consecutive majorants is
(m+1)R q^[-(2m-1)/4], tending to zero. The tail converges normally on compact
sets. Finite initial rational terms give a single-valued meromorphic
function on C, with possible poles at q^a. Near z=0 the normal convergence
allows coefficient extraction and agrees with (3), establishing the
continuation of the normalized generating function itself.

The looser quadratic bound (4) is deliberate: a two-cycle fixes its
unordered antipodal pair, so the stronger tournament no-fixed-pair bound
would be FALSE for this general palette.

## 4. Exact pole classification

### At least one self-reversing color: t>0

All core weights are strictly positive. For a fixed defect a=k, each
μ_i=d_i+1 defines an ordinary partition (d_i) of k. There are p(k) types,
m=k+c<=2k, and the unique core of maximum size is μ=2^k. Every term in
this finite bank has pole at q^k, of order m+1. The unique maximal one
cannot cancel, and all other banks are holomorphic locally by normal
convergence. Thus the poles are EXACTLY

    z=q^k, k>=0, with exact order 2k+1.                  (6)

For μ=2^k, the k antipodal pairs are reversing, while the 2*binomial(k,2)
cross-cycle orbits are nonreversing. Hence

    F_(2^k)=t^k q^[k(k-1)],
    C_(2^k)=(2k)!/(2^k k!) t^k q^(-k²),
    lim_(z->q^k)(1-z/q^k)^(2k+1)B_(q,t)(z)
      =(2k)!/(2^k k!) t^k q^(k²)>0.                   (7)

### No self-reversing color: t=0 (necessarily q even)

Every even cycle contributes a reversing antipodal orbit, so all cores
with even parts have zero weight. The remaining odd cores have defect
2k and maximal size 3k, uniquely at μ=3^k. The tournament proof applies
verbatim with 2 replaced by q. The exact poles are

    z=q^(2k), k>=0, with exact order 3k+1,               (8)
    lim_(z->q^(2k))(1-z/q^(2k))^(3k+1)B_(q,0)(z)
      =(3k)!/(3^k k!) q^(3k²+k)>0.                    (9)

In BOTH cases the radius of the normalized Taylor series is exactly 1,
and its coefficients are not P-recursive, since the meromorphic function
has infinitely many genuine poles. The original unnormalized OGF and EGF
have radius zero and the original sequence is non-P-recursive by
q^N/n!<=a_n<=q^N and the growth argument in `LIMITS_AND_NORMALIZATION.md`.

## 5. What disappears when self-reversing colors are removed?

This is an exact comparison among admissible FINITE palettes, not a
continuous probability model. Fix an even integer q. Its allowed t are
0,2,4,...,q. Relative to any allowed positive t, setting t=0 removes each
core term containing an even part, since its factor t^v becomes zero.
Consequently:

* Every pole q^(2h+1) disappears altogether: an odd defect requires an odd
  number of even cycle lengths, hence at least one even cycle.
* At the surviving location q^(2h), the exact order drops from 4h+1 to
  3h+1, by h, for h>=1.
* The pole at 1 remains simple with amplitude 1.

For odd q there is no t=0 palette; its admissible t are 1,3,...,q. The
polynomial core weights can be formally specialized, but that formal
operation does not manufacture a finite odd-q palette without fixed colors.

## 6. Three concrete structure layers

| Structures | q | t | Poles | Orders |
|---|---:|---:|---|---|
| Tournaments, one of two opposite directions | 2 | 0 | 4^k | 3k+1 |
| Simple graphs, edge/nonedge both unchanged by reversal | 2 | 2 | 2^k | 2k+1 |
| Partial comparisons, forward/backward/undecided | 3 | 1 | 3^k | 2k+1 |

Thus allowing an undecided pair creates even-cycle fixed structures and
changes the symmetry-defect banks, not merely the count of colors. This is
an exact bridge from partial-relation states to analytic singularity type.
It does NOT assert that completion counts or vertex-growth transitions
close on (q,t), on defect alone, or on pair-count histograms. Those observers
can discard the attachment correlations required by the metagraph grammar.

## 7. Validation and provenance

Run `python3 session98_growth/limits/verify_palette_poles.py`.
The checker uses explicit non-removable checks. For every admissible
palette with 2<=q<=4, it exhaustively checks permutation fixed-color counts
and canonical unlabeled orbit counts through n=4, including the neutral
three-state partial relation. It compares full partition enumeration to
rational-core coefficients through n=12 and verifies pole-leading cores
through k=6. These are FINITE-EXACT controls; normal convergence and all-order
pole classification are supplied by the proof above.

The provenance is classical Burnside/Pólya pair-action enumeration,
THM-514's unit-tail core kernel, and this session's independently audited
normalized all-poles argument. The signed-orbit expression is derived in
full here; no claim of discovering the graph or partial-digraph enumerator
or proving literature novelty is made. Only the task's specified palette
extension is covered; no other metagraph processes are inferred.
