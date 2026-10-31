# Independent audit: ternary schedules, rays and certified cycle multiplicity

## Verdict

**PASS. No mathematical correction is needed** in the reviewed `TERNARY_SCHEDULE_THEOREM.md`, including its exact-frontier appendix. The conclusions follow from the independently audited A/B/C universal catalogs and the common weak B invariant. They concern a specified coded family, not arbitrary prescribed paths, all square-sum rays, or an optimal count of all Hamiltonian cycles.

The operators are

- A: Q25, c12, H=(1,3);
- B: Q25, c8, H=(1,3,6,19), with independent shared head colors;
- C: Q49, c4, H=(1,3,6,19).

## New independent controls

`audit.py` imports no producer code. It computes prefixes using the accumulated affine coefficients rather than the producer's recursive prefix routine. It independently expands mixed Hamiltonian pairs and streams cycles, checking every vertex, edge, closing edge and complete prescribed prefix.

Normal and optimized receipts match byte-for-byte with **31,877,061 explicit checks**:

-30 mixed cycle controls: the producer's27 plus three extra odd-input controls;
- all27 common cycle hashes agree exactly;
- closed-form positivity, simplicity, square sums and exact inverse decoding for all3,279 nonempty ternary words through depth7;
- independent explicit-frontier traversal agrees with a bottom-up dynamic count on distinct integer products;
- all published frontier examples are reproduced, including35 at N16,601,562 and145 at N1,000,000,000; the independent additional N10^12 control gives1481;
- exact full-ternary frontier weight conservation is checked by the polynomial identity obtained on substituting v=1-2u, not by floating-point approximation to alpha.

The new mixed controls use certified small input pairs and the already audited universal tables. They do not replace those universal audits or claim to enumerate every mixed schedule. Consumed certificate hashes are recorded in the receipts.

## Uniform typed composition

All three paired catalogs consume and produce the same weak B type: common prefix1,3, required even/odd terminals8/10, complete consecutive intervals and common-label positional parity. The head-color relaxation in B preserves exactly this type, so it can feed A or C without any provenance restriction.

For a finite word w=a v, write Q=Q_a and h=(Q-1)/2. At every n at least the required word threshold, quotient decomposition gives

`m=floor((n-h)/Q)`, `d=n-Qm in[h,3h]`.

Inductively construct the v-prefix pair at m and apply the a paired table. The outer terminal table instead gives the cycle. The empty suffix can use the audited A-only all-size weak B supply. Thus the logical mixed theorem can use the minimal base supply42 through1061 plus the three new catalogs. The reviewed implementation's additional C fallback for small out-of-domain requests may use bases through2081; this is valid but not required for proving the prefix theorem on its promised domain.

The threshold transforms by `n -> Qn+(Q-1)/2`. Its shifted form is `n+1/2 -> Q(n+1/2)`, so for any operator order

`S_w=(85*R_w-1)/2`, `R_w=product Q_a`.

Order affects the vertices, not this threshold. For counts a,b,c of A,B,C, respectively, the prefix has2+2a+4b+4c vertices and threshold `(85*25^(a+b)*49^c-1)/2`.

## Affine formulas and simplicity

With outer-to-inner operators, define `R_t=product_(i<=t) Q_i` and `C_t=sum_(i<=t) c_i*product_(k<i)Q_k`. The block supplied by the next operator is `R_t*H_j+(-1)^j*C_t`. Append `(R_l-C_l,3R_l+C_l)` for a finite word. Every fixed prefix has even length, so local and global signs align.

Within an A block the root multiplier is2; within B/C it is2,3,5. Boundary multipliers are4,6,8 for A,B,C, because `last(H)+Q-c` equals16,36,64. Multiplication by R_t preserves squareness, since R_t is itself a square.

Each operator maps simple positive inputs to simple positive outputs. Opposite-sign lifted collisions would require Q(x-y)=±2c, impossible because0<2c<Q; same-sign collisions require the original labels equal. Its prefix literals avoid the two image residue classes, so no new prefix/lift collision occurs. Composition and nesting then prove all finite prefixes simple and produce simple infinite rays.

## Exact decoding and incompatible prefixes

A can be identified by third vertex13. Otherwise third vertex6 identifies B/C, with fifth vertex17 or45 distinguishing them. Removing the fixed prefix and applying the inverse affine map continues the decoding. This is a prefix code in the operator alphabet, despite the variable numbers of vertices added.

Consequently two operator words in an antichain of the word-prefix order yield vertex-prefixes that are mutually incomparable. An infinite schedule can also be decoded uniquely. The family of rays has cardinality continuum because infinite ternary schedules do, and the map is injective. Each finite initial segment lies in some P_w and hence extends to a Hamiltonian cycle at every sufficiently large interval order.

This is not simultaneous embedding of all branch choices. For example B and C prefixes would force vertex19 to have neighbors6,17,45 in one cycle. More generally, an undirected cycle containing edge1-3 has exactly one traversal beginning1,3, so it cannot realize two incomparable prescribed prefixes. This proves that antichain schedules force distinct **undirected** cycles, without orientation overcounting.

No infinite Hamiltonian-ray conclusion is drawn. These coded rays occupy very few integer labels relative to their exponentially growing block scales.

## Fixed-depth cycle counts

The coding and all-size theorem give2^l distinct cycles for N>=(85*25^l-1)/2 by using A/B words, and3^l for N>=(85*49^l-1)/2 by using all words. At fixed depth l with exactly c copies of C, the number is `binomial(l,c)*2^(l-c)` at threshold `(85*25^(l-c)*49^c-1)/2`. These follow from mutually incomparable prefixes, not an assumption that different construction programs must always yield different cycles.

## Weighted frontier

There is a unique positive alpha satisfying `2*25^(-alpha)+49^(-alpha)=1`: the left side is continuous, strictly decreasing from3 to0. Let X=(2N+1)/85 and N>=2082. Expand a node of product R into all three children exactly when49R<=X. The resulting finite leaf frontier is nonempty-word, prefix-free, and has X/49<R<=X.

Each expansion preserves the sum of weights R^(-alpha). It starts at1 and therefore ends at1. Since each leaf weight is at most `(X/49)^(-alpha)`, the number of leaves is at least `(X/49)^alpha`. Every leaf is affordable and gives a distinct cycle. The use of an implicitly defined real exponent is exact; numerical evaluation of alpha is unnecessary.

The independent polynomial check uses u=25^(-alpha), v=49^(-alpha), and2u+v=1. The sum over a full-ternary frontier is identically1 after replacing v by1-2u. This is an exact integer-coefficient control of the conserved weight law.

## Exact maximum coded frontier

The full affordable-word tree consists of words with R<=X. Every internal node has two children A/B and optionally a third C. Its leaves form an antichain. Any other antichain injects into those leaves by choosing a descendant leaf for each node: antichain subtrees are disjoint. Therefore the leaves maximize cardinality **within this coded tree and cost bound**. This is not an optimality statement for all graph cycles or all possible prescribed paths.

For1<=X<25 the leaf count f(X)=1. For25<=X<49 it is2f(X/25), and for X>=49 it is2f(X/25)+f(X/49). When N>=1062, the root is not a leaf, so every selected word is nonempty and f((2N+1)/85) is a valid cycle lower bound. The full affordable frontier refines the complete-ternary frontier, so its cardinality is no smaller.

A computation memoized by products25^a*49^b has O((log N)^2) arithmetic states. This does not bound the bit complexity by the same expression or imply that all cycles can be written in polylogarithmic time. The publication correctly separates state count from output materialization cost.

## Portable constructor review

The reviewed `radix25_publication/core.py` and `bases.py` implement the intended semantics correctly:

- the empty suffix uses A only and hence does not load Q49 certificates for A/B-only schedules;
- Q25 recursion begins at1062; C recursion begins at2082 when serving its broader weak-pair interface;
- promised prefix domains recurse deeply enough and use the correct suffix word;
- all consumed certificates and base chunks are hash-pinned;
- root-byte offsets, decoded byte-count formula and B invariants are consistent;
- direct endpoints are8 for A/B, or8/15 for even/odd C outputs.

No semantic error was identified. This is a mathematical/constructor review, not a claim of a general security audit.

## Files

`audit.py`, `normal.json`, `optimized.json`, `producer_comparison.json`.

The independent cycle checker takes explicit directories for A/B/C paired and terminal certificates and an explicit small-base JSON file. It does not import any package loader.

## Final coded-family upper bound

The later Kraft upper bound is also valid. For any finite operator-word antichain, extend each selected word to a common larger operator depth in the complete ternary tree. Its descendant weights sum to its own weight because2*25^(-alpha)+49^(-alpha)=1. The descendant sets are disjoint, so `sum R_w^(-alpha)<=1`.

For affordable words R_w<=X, each weight is at least X^(-alpha). Thus every affordable antichain has at most X^alpha words. Together with the complete-frontier lower bound,

`(X/49)^alpha <= f(X) <= X^alpha` for X>=49.

Therefore f(X)=Theta(X^alpha) is an exact asymptotic-order classification of the **maximum certificate-coded prefix antichain**, with alpha defined above. It is emphatically not an upper bound on the total number of Hamiltonian cycles of the graph. No additional certificate assumption is needed.
