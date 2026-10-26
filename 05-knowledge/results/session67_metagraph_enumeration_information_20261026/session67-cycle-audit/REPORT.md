# Session 67 independent core/cycle audit

**Outcome: PASS within the stated proof and computation scopes.** The proposed
cycle-index inversion is correct, but is already proved and implemented in the
repository. It must not be presented as a new all-order enumeration theorem.
The useful remaining distinction is between counting all classes and classifying
or transporting their actual metagraph families.

Source snapshot: `b4234f728c0516b5fc521c001d69e26f3d2508d5`.
Read-only source checkout: `../session67-repo` (sparse checkout; omitted files
were read with `git show HEAD:<path>`). This audit changes no repository files.

## 1. Inherited result map and priority

- `01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md`,
  section 2: homogeneous-pair graph components are the canonical maximal
  transitive modules.
- `05-knowledge/results/weighted_core_transport_20261020/core_count.md`,
  sections 1–2: unique weighted-core class decomposition and actual
  `Aut(Q)` weight orbits; each fixed core's rational response.
- `.../weighted_core_transport_20261020/inheritance.md`, sections 2–4 and 7–8:
  extension dichotomy, alternating core generation, exact labelled transition
  kernel, full cycle-index substitution and odd-cycle hierarchy.
- `.../defect_boundary_20261020/boundary.md`, section 3: explicit signed
  ordered-list coefficient inversion, including the phase powers of cycle length.
- `.../tournament_future_boundary_20261020/README.md`, sections 2–4:
  horizon-dependent response DP; ordinary tournament counts through 200,
  ordinary core counts and fixed-point-free core mass through 100.

There is also direct classical terminology and priority. Boudabbous–Pouzet,
[*The morphology of infinite tournaments. Application to the growth of their
profile*](https://arxiv.org/pdf/0801.4069), Theorem 1, gives the unique
lexicographic decomposition into acyclic blocks indexed by an acyclically
indecomposable tournament. For finite tournaments, “acyclically indecomposable”
is exactly “no homogeneous pair”: a nontrivial transitive module has consecutive
homogeneous vertices, and a homogeneous pair itself is a transitive module.
The paper credits the terminology to Culus–Jouve (2005). Do not claim historical
priority for the canonical decomposition. The present audit independently
checks the finite argument, rather than importing infinite-case machinery.

## 2. Species and coefficient verdict

Let `T` be all tournaments and `C` the homogeneous-pair-free tournaments,
including the empty object in each. Let `L+` be nonempty linear orders. Then

    T = C o L+,
    Z_T(p) = Z_C(p_i/(1-p_i)),
    Z_C(p) = Z_T(p_i/(1+p_i)).

These identities hold in rational formal series completed by weighted degree
`deg(p_i)=i`. No analytic convergence is claimed. The proof is a natural
bijection on labelled sets: canonical transitive blocks determine the quotient,
and a quotient isomorphism lifts exactly when it respects the block orders.
For an automorphism i-cycle of quotient vertices, equal chain length k creates
k vertex i-cycles; this gives `sum_(k>=1) p_i^k`. Since the two substitutions are
formal inverses, the inverse identity follows coefficientwise.

For odd partition `mu` with multiplicities `m_d`, let `F(mu)=2^e(mu)` denote
Davis's fixed-tournament count and let `C(mu)` denote fixed pair-free tournaments.
Writing

    L(m,k) = m!/k! * binom(m-1,k-1),  1<=k<=m,
    L(0,0)=1,

with every other case zero, the exact formula is

    C(mu) = sum_(k_d<=m_d) F((d^k_d)_d)
               product_d (-d)^(m_d-k_d) L(m_d,k_d).

Replacing `-d` by `d` gives the inverse transform. In particular
`C(mu)=z_mu [p_mu] Z_C` is a nonnegative integer because it is literally a
fixed-point count, not because arbitrary signed series are positive. The
cycle-index coefficient itself need not be integral: type `(3,3)` has
`F=32`, `C=20`, and coefficients `16/9`, `10/9` respectively.
Any even permutation cycle gives zero for both classes.

Boundary checks:

- Empty object: `T_0=C_0=1`; both substitutions preserve the constant.
- One vertex: exactly one core.
- Two vertices: no pair-free core; every tournament is one chain.
- One odd cycle, or any type with no repeated cycle length: no nontrivial
  clone correction; `C(mu)=F(mu)`.
- Ordinary generating functions cannot be substituted as though these were
  labelled scalar compositions; the full permutation cycle action is essential.
- `Z_C(0,x^2,x^3,...)` counts automorphism mass with no fixed vertices, not
  the cardinality of a subclass of cores, and does not determine the full
  multivariate boundary by itself.

The existing future-response DP updates its response by the number `k` of
actual cycles, not the signed-substitution excess `s`; that distinction is
correct and load-bearing. Its necessity claim is wordwise distinguishability
for permitted continuations, not minimality among every symbolic algorithm.

## 3. Exact limitations and useful immediate corollaries

These are elementary deductions from the inherited mechanism, not a new claim
that a compact arithmetic counting algorithm is impossible.

### 3.1 One added vertex suffices for pair-free universality

Every tournament `T` of order `n>=2` is an induced subtournament of a pair-free
tournament of order `n+1`.

Proof: if T has b canonical chains, alternating the new incidence bits on
each chain gives, after excluding clone insertions,

    2^b-2b  if n=b,
    2^b-1   if n=b+1,
    2^b     if n>=b+2.

The first case necessarily has `b>=3` when `n>=2`, since no pair-free tournament
has order two. All three counts are positive. The exceptional parent of order
one has no pair-free one-vertex extension, so the size qualification matters.

Consequently the hereditary closure of pair-free cores is the entire class of
tournaments. Pair-free does not mean the binary relation has become simple.
If `a_n` and `c_n` count ordinary classes, forgetting a distinguished vertex
in these extensions gives

    a_n <= (n+1)c_(n+1),  n>=2.

The factor is a bound on vertex-orbits, not an assertion that different deleted
vertices always give different parents.

For comparison only, the stronger modular notion has a classical two-added-
vertex embedding bound: Brignall–Ruškuc–Vatter,
[*Simple extensions of combinatorial structures*](https://research-repository.st-andrews.ac.uk/handle/10023/1997),
Theorem 2.1, recalling Erdős–Fried–Hajnal–Milner. The pair-free corollary above
requires neither that theorem nor its exceptional-case classification.

### 3.2 A one-pair contraction counting bound

Every non-pair-free n-class arises by replacing one vertex of some
(n−1)-class by a transitive pair. The replacement is unique up to isomorphism
once the marked parent vertex is chosen. Thus

    0 <= a_n-c_n <= (n-1)a_(n-1),  n>=1.

The surjection is not a bijection: different marked parents may inflate to
the same child. This is another reason naive family-transition multiplicities
cannot simply be summed as class counts.

### 3.3 No fixed finite chain-inflation bank can cover all orders

If a fixed finite bank has largest core order B, its total contribution at
order n is at most a constant times `binom(n-1,B-1)` (or a sum over smaller
orders), hence `O(n^(B-1))`. But

    a_n >= 2^binom(n,2)/n!.

The proportion represented by that fixed bank therefore tends to zero. This
rules out a finite list of fixed cores with only chain-weight variation as a
complete A000568 classification. It does not rule out a recursive grammar
that creates new cores, nor the existing response DP.

### 3.4 Pair-free cores retain quadratic information scale

For a uniform labelled tournament, a fixed pair is homogeneous with probability
`2^(-(n-2))`. The union bound gives, for `n>=6`,

    c_n >= (1-binom(n,2)*2^(-(n-2))) * 2^binom(n,2)/n!.

Therefore any fixed-length identifier that distinguishes every ordinary
pair-free n-class needs at least

    ceil(log2(c_n))
    >= ceil(binom(n,2)-log2(n!)+log2(1-binom(n,2)*2^(-(n-2))))

bits in its worst case. The leading scale is still quadratic, despite the
lossless and useful canonical-chain compression on structured families.
Again this is an object-encoding lower bound, not a time lower bound for
computing the scalar number c_n or a_n.

## 4. Independent computation

Reproduce from this directory:

    python3 verify.py > receipt.json
    python3 -O verify.py > receipt-optimized.json
    cmp receipt.json receipt-optimized.json

The script uses no inherited implementation. All checks are active under `-O`.
Both executions completed, with byte-identical output and **2,239,411 checks**.

Coverage:

1. All 2,035 odd cycle types of total weight <=30: signed coefficient
   nonnegativity/integrality, subset bound, exact positive inverse.
2. Every odd cycle type through order six; every nonidentity odd cycle type
   at orders seven and eight: literal orientation-orbit enumeration and direct
   homogeneous-pair tests (23 permutation types total).
3. Every labelled tournament parent through order six, every one-vertex
   attachment: 2,131,018 extension instances. Direct row comparisons verify
   the entire inherited core-size transition kernel, nondecreasing core size,
   accepted alternating count, and the one-point universality corollary.
4. Independent isomorphism-orbit exhaustion through order six, not merely
   a Burnside specialization: ordinary counts `1,1,1,2,4,12,56`; pair-free
   counts `1,1,0,1,2,5,26`.
5. Independently produced ordinary/core series through30 match the inherited
   response-DP implementation coefficientwise in both modes. This last check
   is a cross-implementation comparison, not another independent census.
6. The two counting inequalities above hold for every available coefficient
   (through30, or through29 for the one-point shifted inequality).

SHA-256:

- `verify.py`: `1e31e7b92a44eb25aeb63dfb4fb00a4c0f67db397cc9d68c207939914686aab4`
- `receipt.json` and `receipt-optimized.json`:
  `77dd9234cf48fd7ad8d26a843f6faa61950057a0359b1e2af1355ee040a5ca0c`

## 5. Stopping boundary

No proof defect was found in the inherited canonical decomposition, signed
inversion, extension kernel, or finite-horizon DP statement. This independent
audit does **not** rerun the inherited degree-100/200 jobs, classify all
order-nine cores, certify metagraph incidence, prove a new complexity bound,
or show that a nonlinear scalar recurrence cannot exist.

The strongest practical next step is to state the user's goal in two tracks:
all-order exact **counts** already have closed procedures; structural
**family classification and transport** still need operation-stable object
information beyond cycle indices. The earlier rigid-six-core hostile already
proves that even a complete core cycle index is insufficient for individual
child-automorphism profiles. Reproving species inversion does not close that
missing object-level information.
