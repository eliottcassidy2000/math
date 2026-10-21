# Independent audit of the frozen zero-heavy and palette packet

Audit date: 2026-10-21. Frozen producer manifest SHA-256:
`a5490b17a42062db222676a44a1eb2ef1c3f0b55499b2389cc70691886f9112c`.
The audit concerns exactly `zero_heavy.md` and `palettes.md` and the files
and dependencies authenticated by that manifest. It excludes all later
corridor exploration or evolving addenda. No repository or remote is changed.

## Mathematical verdict

**PASS in the producer's stated scopes.** The all-order results follow from
the proof arguments below. The separate deterministic receipts are finite
controls, not computational proofs of unrestricted statements.

### 1. Deficit classification and exact weighted classes

The inherited nonnegative deficit identity holds for *every* good weighted
order. For deficit below `2k+1`, its first summand forces `r=k`. Since the
`k` positive integer backward-edge costs sum to at most `k`, every cost is
one and every incident vertex has weight one. Thus `W=t`, and the identity
becomes `delta=(k+1)(2k-t)+kh+U`. In particular the correction from shared
endpoints at `k` to first shared endpoints at `k+1` is necessary and valid.

For `k>=3` and `0<=delta<=2k-3`, the only feasible pairs `(d,h)` are
`(0,0)`, `(0,1)`, and `(1,0)`. With no isolated incident vertices, `d=0`
means a matching. At `d=1`, the sum of positive degree excesses is one,
so the graph is exactly one P3 and `k-2` independent edges. The endpoint
count, actual occupied gap count, and full weighted tournament are retained.
No clone contraction is used to rescue an invalid record.

The three terms of formula (8) are correct:

- Full matching: `m` occupied gaps, capped deficits `0,...,k-1`. Since
  `delta<2k`, inclusion-exclusion stops after the single violating part,
  giving `C(m+delta-1,delta)-m C(m+delta-k-1,delta-k)`.
- One missing gap: `mE-J` admissible marked omissions; residual deficit
  `delta-k` goes into `m-1` occupied gaps. It is at most `k-3`, so the
  cap is inactive and the binomial is `C(m+delta-k-2,delta-k)`.
- One P3: there are `m-1` occupied gaps; residual deficit `delta-k-1`
  is at most `k-4`. The binomial is `C(m+delta-k-3,delta-k-1)`.

The negative-lower-index convention is essential and consistently applied.
The two boundary specializations (2), (3), (6), and (7) follow by elementary
simplification. In particular the forbidden matching incidences are
`J=F_k` generally, and `J=(2k-3)I_(k-1)` strongly. Deleting an adjacent
edge of a cut-crossing matching leaves a cut-crossing matching; the edge
cannot occupy either exterior end. Insertion in one of `2k-3` interior
positions is the inverse, so the strong correction has no missing factor.

The P3 recurrence is also correct. There is exactly one irreducible cut
block containing the P3; the prefix and suffix are arbitrary matchings.
This gives total `(k-1)F_k = sum P_j F_a F_b` over `j+a+b=k`. The two
matching factors are genuinely unrestricted, not irreducible matching
counts. The first values `3,24,246,3012` are correct.

### 2. Pair freeness, quotient uniqueness, and strongness

In a full-gap record, every endpoint has a backward partner and every two
endpoints are separated by an actual gap vertex. These distinguish all
possible pairs. For one omitted gap, a homogeneous pair exists exactly
when the two newly adjacent endpoints are matched to each other. It is
not claimed that this is the *only* homogeneous pair in a rejected record;
neighboring gap vertices can participate in a larger transitive clone chain.

Every weight-above-one vertex avoids all backward edges in every good
order. Their induced tournament is transitive, so their total order is
intrinsic and every weighted isomorphism preserves each position. Their
incidences partition the unit vertices into intrinsically ordered cells.
The full-matching corridor has at most one unit gap because
`delta<2(k-1)`. Apart from its exceptional cell, all cells are singletons.
At an exterior unit gap, the external matching partner fixes the two roles.
At an interior unit gap, either fixed external partners distinguish the
three roles, or the local tournament is C3. A C3 rotation changes a local
presentation but does not create a different matching/weight record.

Missing-gap records have every remaining gap weight at least three;
P3 records have every gap weight at least four. Hence their decoding has
no extra unit-gap ambiguity. A retained missing interior gap has a forward
edge between its two endpoints, which fixes their roles. The three families
have different actual vertex counts and therefore cannot be isomorphic.
These facts establish class counting, not just presentation counting.

After the clone rejection, every consecutive displayed arc is forward.
The resulting Hamilton path makes SCCs contiguous. A component cut is
one-way exactly when the backward graph does not cross its corresponding
proper endpoint cut. Exterior gaps are omitted in strong records. The
strong filter is therefore precisely cut-crossing, not connectedness of
the sparse backward graph.

The explicit scopes are respected: `delta=k` is proved by this argument
for `k>=3`, `delta=k+1` for `k>=4`, and formula (8) for the stated corridor.
The `k=2` values `42/11` and `81/15`, and `k=3,delta=4`, are finite-exact
controls outside that one-unit-gap proof. The pinned complete `k=2`
compiler receipt agrees. Multiple-unit-gap all-order claims, unrestricted
A000568, and feedback-vertex counting are not implied.

### 3. Exact source gate and physical-port repair

The independent implementation reconstructs the first seven coordinates
from the pinned source's literal coordinate table, then reconstructs all
nine ports using its reflection triples. It does not import the producer's
field arithmetic or accept its edge list as a definition. Direct exact
squared distances recover 16 distinct points and all 31 unit edges,
including the actual boundary contacts. A separate fixed-variable-order
colouring recursion obtains the source vector `(24,24,24,0,0,0)` and the
port-7 repair vector `(60,84,84,8,8,8)`.

The note's literal 16-entry root-3 colouring and all six receipt witnesses
are parsed and checked on that induced graph, with exact port names and
colours. A repeated colour is never treated as a geometric identification.

For a fixed root colour `c`, each copy excluding `c` needs at least one
changed physical port, so the lower bound is `N-m_c`. The stated one-port
repair attains it independently in each excluding copy. Root-only incidence
allows the extensions to combine. Thus the global radius is
`N-max_c m_c`. The explicit repair at fixed `c` has exactly
`24^(m_c) 8^(N-m_c)` root-`c` extensions. This counts internal extensions
of *that one repaired boundary*, not distinct repairs and not total
extensions over every possible root. If `c` does not maximize its load,
“nearest” is understood with that root prescribed, as the preceding
sentence of the source states.

Since the six loads sum to `3N`, the largest is at least `ceil(N/2)`;
complementary palettes attain equality. Therefore `floor(N/2)` is the
optimal radius in this architecture. Dividing by the actual `9N` ports
gives at most `1/18`, attained for even `N`. The case `N=1` correctly has
zero radius. The counts `15N+1` vertices and `31N` edges follow from the
actual one-root union, not from an abstract palette quotient.

### 4. Actual all-finite geometric realization and forcing limitation

For nonzero old and new points `a,b`, each prohibited equation
`|a-rho_t b|^2=e`, with `e=0` or `1`, becomes a nonzero polynomial of degree
at most two after clearing `1+3t^2`. It cannot vanish identically because
`a conjugate(b)` is nonzero and its real projection varies under rotation.
There are at most four excluded positive integer parameters per nonroot
cross-pair. At stage `j`, `(15j)(15)` such pairs give the valid bound
`900j+1`. Rotations fix the root, preserve all old points and named ports,
and add no cross-copy unit contact or coincidence. Induction realizes
every finite N in the stated field.

The independent control parses all 91 literal coordinates, reconstructs
every rotation, verifies the *first* surviving parameter, and checks all
4,095 unordered distances. The resulting 186 edges are precisely the
union of six source edge sets; all 54 boundary ports are distinct. Every
literal nearest-repair colouring and the four-colouring is checked on that
actual graph. Small whole-graph root fibres are counted directly without
using the product formula in the enumeration.

Root-only union remains four-colourable by aligning one root colour.
Consequently unbounded absolute conditional robustness cannot certify
unconditional six-colour forcing. The fresh-sixth-colour obstruction to
forcing an equality in a five-colourable context is valid. The proposed
five-colour circle-saturation target still requires the universal
neighbour-colour predicate for an actual source graph and remains OPEN.
There is no hidden six-chromatic or universal synchronizer conclusion.

## Independence and finite scope

`independent_audit.py` imports no producer or inherited checker. Its
canonicalizer is generic directed weighted colour refinement followed by
exhaustive permutations of all unresolved invariant cells. It does not
assume the large-gap decoder whose uniqueness is being tested. Literal
clone checks inspect every vertex pair. Actual strongness is computed by
Warshall reachability and compared to an independent cut-coverage array.
Weighted feedback uses an independent append-first subset dynamic program.

The finite universe is deliberately explicit:

- all matching and one-P3 graphs through `k=6`;
- formula (8) versus bounded polynomial coefficients throughout every
  corridor for `k=3,...,32`;
- all weighted records throughout the entire `k=3,4` corridor, plus the
  `k=2` two boundary layers and `k=3,delta=4` exception;
- the pinned complete `k=2` compiler's exact mass coefficients;
- the source gate, full six-copy geometry, six literal nearest repairs,
  small direct full-graph fibre counts, and every palette multiset through
  `N=4` compared with independent repaired-subset optimization.

All runtime gates use explicit exceptions, not removable assertions.
Normal and optimized execution and their exact totals are recorded in
`receipt.normal.json`, `receipt.optimized.json`, and `verification.json`.
The audit manifest pins these outputs and all audit artifacts. Frozen
input hashes are checked both before and after each full independent run.
