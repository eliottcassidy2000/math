# What a famous open hypothesis would have to supply

Chan–Pak is the exact range statement in [Chan and Pak, §6.5](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf). The following is a conditional reduction with an explicit *unproved physical hypothesis*, not a proof from Goldbach alone.

## A three-row Goldbach reduction that does not silently use 7

Let HH be the set of directed Hamiltonian-path counts of tournaments. The following is a structural gluing target, an independently unproved physical construction target: if one drops the physical gluing requirement and retains only bare existence, this is not an independent route.

> **Overlapping prime-gluing hypothesis (unproved).** For some explicit B, for every pair of odd primes p,q other than 7 that lie in HH, and every c in {−3,−1,1} with p+q−c ≥ B, one can construct a tournament with count p+q−c from witnesses for p,q. The rule may overlap and prepare the witnesses; it is not required to preserve their disjoint union.

**Conditional proposition.** Binary Goldbach, this physical hypothesis, and verification of admissible counts below max(B,15) imply Chan–Pak.

For a putative least missing admissible odd N above the boundary, choose c from those three offsets so that N+c−7 is divisible by 3. The even integer E=N+c satisfies E−7≥5, so E−7 is composite. A Goldbach pair for E cannot use 7, or the other leg would be E−7; it cannot use 2 for E>4. Each odd prime leg is at most E−3≤N−2; hence both have witnesses by minimality (21 is not a prime). Apply the row for c, giving E−c=N. The three offsets are important: ordinary Goldbach alone does not say that its pair can avoid 7 at an arbitrary even target.

## Why the known joins cannot serve as those rows

For two vertex-disjoint induced subtournaments (U,V) in a tournament (T), the positive odd-cycle packing formula gives
`H(T) ≥ H(U)H(V)`.
Indeed, choose any pack of pairwise disjoint directed odd cycles in (U), and independently one in (V), then unite them. The weight (2^k) for (k) cycles multiplies and the union is an injection; all other contributions are nonnegative. Thus any construction that preserves two disjoint induced prime witnesses has count at least (pq), larger than (p+q-c) for all sufficiently large targets. More generally, if the two witnesses share a separator (S), the same argument bounds the output below by (H(U-S)H(V-S)). An overlapping prime-gluing theorem would have to arrange a large count drop at that separator. Neither generic union, opaque substitution, nor the existing two-apex bridge provides it.

The same failure is visible in the known exact two-apex reverse bridge, whose absolute count is `hj+PQ−1`: both source counts (h,j) multiply, even if a local swap provides an additive *difference*. A result about a signed difference does not replace the gluing theorem.

## Other named conjectures: missing physical implications

- **Ordinary Collatz.** Convergence would supply a route on integers, not a backward construction of a tournament on that route. The existing signed arithmetic has a decreasing selector and a bypass of 7 *unconditionally*. The open arrow is its physical inverse, including bypasses, so imposing the usual Collatz conjecture on that argument adds no needed arithmetic premise.
- **Zaremba.** Bounded continued-fraction digits could supply exact scalar denominators in a positive matrix semigroup. The ribbon language currently has no authenticated functor from those digit operations carrying that scalar to a physical Hamiltonian count. Even equidistribution in every fixed congruence class fails to give exact saturation (e.g. the odd integers `n!+2n+1`).
- **GRH/Artin and ordinary prime-density conjectures.** Primes in phases can authenticate an auxiliary input or ensure one exists, but cannot by themselves turn the absolute count `hj+PQ−1`, or a known modular wheel, into a pointwise inverse. In particular the established 7-free Riesel rows only transfer along reducible rays; ordinary primes and the (7p) bottleneck survive.

I found no implication from these famous open conjectures without a new independently unproved physical hypothesis. The Goldbach reduction identifies such an obligation and an unconditional obstruction to the most immediate way of attempting it.

The [induced-gap theorem](../triangular_backbone_chanpak_20261116/CORE.md)
adds a test even for overlapping constructions. If the rule keeps a
q-path tournament induced and raises its count by only p−c, then
p−c≥2⌈log_3 q⌉. A rule retaining the full large witness for every
pair p,q is impossible already at any fixed small p and unbounded q.
A genuine prime-gluing row must change the large witness in those
cases, or force a balanced arithmetic decomposition. Ordinary
Goldbach gives no such balance.

### What a Zaremba compiler would have to accelerate

Suppose a one-letter-per-step compiler of the bounded
continued-fraction semigroup were to append at least one vertex per
digit, end in a strong tournament, and read H=2q+1 for the denominator
q of the processed continued fraction. On a string of t ones that
denominator is Fibonacci, so H=O(φ^t). Busch's strong lower bound
instead forces H≳5^(t/3), and φ<5^(1/3). Thus that literal compiler
cannot work. This does not rule out a compiler with blocked,
non-induced or decomposable operations; Zaremba's conjecture alone
does not supply one. See the
[original density-one approach of Bourgain and Kontorovich](https://arxiv.org/abs/1107.3776)
for the arithmetic, and [Busch's bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3).

### Even a one-vertex overlap cannot implement the three offsets

If prime-count induced witnesses for p,q≥3 intersect in at most
one vertex, their union lies in no tournament of count below
p+q+3. Indeed each witness has just one nontrivial strong component
by multiplicativity and primality. If these cores are disjoint,
H≥pq≥p+q+3. Otherwise their common vertex v belongs to a directed
triangle in each core. Such a triangle exists through any vertex
of a nontrivial strong tournament: otherwise the entire in-neighbor
set beats the out-neighbor set, contradicting strongness. The two
triangles make a strong five-vertex induced tournament (they share
v), whose count is at least 9 by Busch's bound and odd parity.
Among its odd-cycle packs, only the empty pack and the two
individual triangles are wholly contained in one witness,
with total weight 5. Thus there is mixed weight at least 4 in
addition to the union of packs inside the two witnesses,
which weighs p+q−1. This proves the claimed lower bound.
Of the proposed offsets only c=−3 meets the lower bound;
the c=−1 and c=+1 rows require more overlap or modification.

### A two-vertex overlap at the low offsets forces a Fermat summand

Assume the two induced witnesses intersect in at most two
vertices and the joined count is at most p+q+1. By the previous
bound, their overlap must be a pair {r,s}. The packs internal to
either witness weigh p+q−1 in total, since their only common
pack is empty. If a triangle in one witness avoids {r,s}, it is
disjoint from any triangle of the other, giving another weight
4. Thus every triangle in both witnesses meets {r,s}. If each
witness has a triangle meeting the pair in just one vertex,
choose these two triangles. If they use different vertices,
the disjoint pack again contributes 4; if they use the same
one, the two-triangle strong wedge contributes mixed weight
at least 4 by the preceding proof. Both are impossible.
Consequently in at least one witness **all triangles use both
r and s**.

Deleting either r or s from that witness leaves a transitive
tournament. They therefore have two insertion positions in
the same transitive spine. To have any triangle, the arc
between them must run backward across a nonempty interval of
m spine vertices. All directed cycles are the backward arc
closed by an increasing subset of those m vertices; the odd
cycles use odd subsets, and none can be packed together.
The witness count is 1+2^m. Primality implies m is a power of
two (an odd divisor of m would give an elementary factorization).
Thus at least one summand must be a Fermat prime.
This restriction is for the low offsets and for **induced**
witnesses with at most two common vertices, not for arbitrary
physical preparation of the gluing inputs.
