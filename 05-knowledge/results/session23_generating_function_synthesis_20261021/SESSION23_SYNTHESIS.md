# Session23: generating functions, symmetry and recoverable structure

Session started 2026-10-21 19:23 UTC. This report records the mathematical conclusions and their scope. Publication verification is recorded separately in the enclosing checkpoint; the report itself is not an execution-status notice.

## What changed

The session followed the user's request to prove additional generating functions and recover older repo concepts. It did not attempt to turn numerical coincidences into a universal tournament classification. The strongest advances come from three explicit carriers: ordered strong-component words, canonical maximal singleton/triangle-word modules, and actual Hamilton-path orbit coordinates.

### 1. Exact recurrence orders across an infinite black-degree family

The [refined metagraph packet](../refined_metagraph_generating_functions_20261021/README.md) proves multivariate rationality retaining SCC count, prime valuations of automorphism order, presentation count h and blue count b. The square substitution for a mirrored half applies to every additive marker, not just the vertex variable.

For p=1-z-z³, q=p(z²), the degree-3^r counting series has exact reduced denominator p^(r+1)q^e, where e=floor(r/2)+1 except e(1)=2. Its minimal eventual recurrence orders start9,18,21,24,33,36,45,48. The first two exact orders were inherited; the all-r classification is the extension. Unknown atom/center coefficients and recurrence onset remain separate.

A counter-boundary matters: marking the raw automorphism integer as an exponent already gives a series that is not D-finite in that marker over Q(z), already on the rigid word family. Marking its prime valuations instead is rational. Multiplicative information has to be put in additive coordinates.

### 2. Geometry and sampling without a complete atom inventory

The [renewal packet](../renewal_tournament_limit_laws_20261021/README.md) proves that, at fixed realizable h, the maximal nonrigid skeleton dominates. Its neutral-gap proportions approach a uniform simplex. For h=3^r, the exceptional R5 components approach uniform order-statistic positions.

The neutral-triangle count has mean mu*n+O(1), variance sigma²*n+O(1), and a Gaussian limit, where rho+rho³=1, mu=rho²/(1+3rho²), sigma²=rho²(1+rho²)/(1+3rho²)³. On each supported vertex parity, self-converse fixed-deficit slices have twice the variance: the two halves are mirrored, not independent. The raw Hamilton-path count grows exponentially in probability at fixed path-orbit count h, because automorphisms themselves grow exponentially.

An exact sampler on K1/C3/R5 words gives an asymptotically uniform sampler of all fundamental degree-3^r nodes, with total-variation error O_r(1/n). This is not unrestricted tournament sampling; no explicit finite-n error constant is asserted without the remaining finite inventories.

### 3. Marked almost-Hamiltonian paths expose rational count channels

The [marked-path packet](../marked_path_tournament_series_20261021/README.md) counts actual isomorphism classes with an ordered directed k-vertex path, leaving m vertices unmarked. Its exact count is2^binom(k-1,2) F_m(1,k), with the earlier odd-cycle Burnside channel F. Removing the free-arc prefactor gives an exact rational generating function whose simple poles are indexed by possible odd-partition cycle counts of m. Converse symmetry yields a separate parity channel with bases4^K.

The action on marked paths is free for omitted-vertex deficit at most2; it already fails at deficit3. The first correction is a genuine3-cycle term, not an assumed probability of asymmetry. Every odd self-converse tournament has one converse-fixed automorphism orbit of roots, not necessarily one distinguished vertex.

The accompanying narrow correction repairs THM-260's mixed-cycle centralizer denominator: equal-length cycle multiplicities, not the total number of cycles, supply the factorials. Its existing(5,3)table entry and rooted prefix remain unchanged.

### 4. An older open generating-function question is closed

The [fixed-symmetry-slack packet](../fixed_symmetry_slack_generating_functions_20261021/README.md) proves rationality for every fixed d=3f(T)-s(Aut(T)), where f is minimum feedback and s counts vertices moved by at least one automorphism. Older work had a noninjective bounded-exceptional grammar and explicit formulas only at low slack.

The repair is canonical: maximal modules in W=SEQ(K1,C3) form a unique partition. The quotient is bounded; short words and the actual quotient stabilizer are retained. Every remaining long-word position is fixed individually, so the grammar becomes an exact finite class compiler. Its ordinary and self-converse series are rational; no efficient inventory of all finite configurations is claimed.

This is a different finiteness mechanism from fixed h. Strong tournaments C3[K1,W,K1] can have fixed slack3 and unbounded order. The diagonal-C3 synchronization example is preserved by the canonical quotient C3[11,11,11], rather than incorrectly splitting its group into independent rotations.

### 5. A common limit law does not mean the same tournament family

The [separately reviewed slack-renewal transfer](../fixed_symmetry_slack_generating_functions_20261021/session23-root/SLACK_RENEWAL_TRANSFER.md) shows that every nonempty fixed-slack stratum has the same feedback-density constant mu and Gaussian variance sigma². The moved-support fraction tends to3mu, approximately0.5827620121. The self-converse slices again double the variance on their supported parity classes.

The explicit map is to finitely decorated word positions over vertex lengths1 and3. It preserves the marked counting statistic, not the original tournaments, their wiring, or any Collatz trajectory.

### 6. Sharp recurrence orders for every nonempty slack layer

The [sharp-pole follow-on](../sharp_symmetry_slack_recurrences_20261021/README.md) combines an actual-orbit separator upper bound with an injective nested-wrapper lower family. Write d=3a+4b with a≥0 and b∈{0,1,2}; this represents every nonempty slack layer, namely all d≥0 except1,2,5. Set m=2a+b+1 and k=ceil(m/2). The exact reduced denominators are:

- Ordinary classes: p^m, minimal eventual recurrence order3m.
- Self-converse classes: q^k, order6k.
- Converse-merged classes: p^m q^k, order3m+6k.

For d≥1, put H=4d²+d. Safe numerator-degree bounds are H+3m, H+6k, and H+3m+6k, respectively. Thus the corresponding recurrence is guaranteed for coefficient index greater than that bound. These are sufficient bounds, not the actual numerator degrees or earliest onsets. They do not compute the initial coefficients. They give linear-in-d recurrence orders and a quadratic-in-d safe initial range.

The lower-family operation is explicit: Wrap(U)=C3[K1,U,K1] preserves the actual automorphism group, adds one feedback edge and three units of slack, and respects converse, provided |U|≥2. Word gaps around nested wrappers attain the upper pole bounds. This connects an earlier single-backedge construction to an all-slack enumeration theorem.

### 7. Two square constructions and an explicit Fibonacci-product bijection

The [split-word spectrum packet](../split_word_metagraph_spectra_20261021/README.md) retains the actual Hamilton-path orbits of the slack3 strong core C3[K1,W,K1]. With a singleton letters and t triangle letters, its path-orbit set is one direct state plus a product of2-choice and4-choice split states. Therefore H=3^t(1+2^(a+2t)), |Aut|=3^t, and h=1+2^(a+2t). The vertex and spectral-exponent coordinates differ: the two letters contribute(1,1) and(3,2), respectively.

At fixed spectral exponent e≥2 the strong-type polynomial is

P_e(z)=sum_(0≤t≤floor(e/2)) binom(e-t,t) z^(e+t+2).

The binomial is binom(e-t,t), not binom(e,t). Ordinary neutral padding gives P_e/p². The explicit converse polynomials then yield the complete fundamental-degree support WITHIN symmetry slack3: powers2^(e-1) for e≥2, triangular numbers binom(2^r,2) for r≥1, and odd numbers1+2^e for e≥2. The corresponding nonzero self-converse degree slices have exact recurrence order6; the non-self-converse slices have order12. This is not a replacement for unknown contributions from other slack layers.

For a reflected even word, the split-state carrier is literally X×X and reflection swaps its coordinates. Blue states are the diagonal plus the direct state; black local incidences correspond to two-element subsets of X. This is a path-level bijection, not merely a count identity.

A different CLASS-level square comes from cutting a word at the midpoint of its spectral exponent. Words of half exponent r are Fibonacci-counted. A total-exponent2r word either has this cut, or has a unique triangle letter straddling it. These two cases yield two Cartesian squares, with converse swapping their factors. The separate empty-tournament unit plus the balanced-cut cores up to half-exponent budget r are explicitly bijective to the older pair-of-tilings family. Their cumulative counts are1,2,6,15,40,104,273. The classical Fibonacci identity and earlier sequence observation are inherited; the new result is the specified tournament map, its inverse, and the cut information it retains.

The distinction from complete wiring is demonstrated, not just cautioned. The order-10 cores W=1331 and W=3113 share (n,H,|Aut|,h,b)=(10,585,9,65,9) and28 local black incidences each. Their actual merged destination sets are disjoint,28 destinations each. The ordered source word and literal path/port remain indispensable. This does not claim their abstract unlabelled neighbourhoods are nonisomorphic; it proves these scalar invariants do not determine the actual destinations.

## What this session does not prove

- Complete metagraph wiring, black-core connectivity, or all bridges.
- Complete strong-atom or center coefficient inventories.
- Uniform153N coverage, a uniform long-horizon R' bound, or universal Collatz entry. No credit ledger was altered in this generating-function session.
- A general Lonely Runner or Hadwiger–Nelson result.
- A fastest general A000568 algorithm. Rational slice recurrences and marked channels are enumeration tools, not a global complexity theorem.

## Review and reproducibility

Every promoted packet distinguishes actual tournaments from formal/synthetic controls, retains source pins, and has an independent proof review with normal and optimized replays. The checks support the all-order proofs; finite numerical agreement is not used as proof of a limit or a general classification. Full commands, exact universes, receipts and remaining boundaries are in the linked packets.

The recovered source trail includes the earlier SCC/center laws, the continuation-width symmetry grammar, THM-284’s proved single-backedge formula, and the older metagraph-spectrum/torus-session Fibonacci discussion. No global literature-priority claim is made. Future work should retain these source paths rather than infer novelty from a search miss.

## Next proof obligations

1. Compute the finite configuration numerators in useful fixed-slack layers and exploit the proved recurrences, without treating a denominator as the initial data.
2. Retain the literal split-word/port carrier in a wiring classifier; the order-10 hostile specifies information that scalar summaries discard.
3. Test which fixed-h and fixed-slack intersections admit similarly explicit slice polynomials beyond slack3.
4. Keep the separate 153N source/clock ledger unchanged until a genuine all-horizon, non-recycled payment theorem is proved. Shared word grammars alone do not discharge it.
