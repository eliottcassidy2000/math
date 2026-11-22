# The signed fold: the blue path spectrum and its linear-forest coefficients

**Status: PROVED.** The signed fold controls the number of fixed Hamiltonian
paths, **not** the total path count. In particular, this does not settle
Chan–Pak. It separates a surjective fixed-point channel from the still
constrained counts of mirror pairs.

Closest prior work: [Session29 coherent half-paths](../session29_framed_profiles_20261022/session29-coherent-halfpaths/REPORT.md) §§1–5 and §§11–13 already prove the controller tournament, the turn product and its mean and parity consequences, and the binary construction (with strong witnesses). The equivalent sign conventions below are included for a self-contained coefficient proof. Other links: THM-582/583 (oddness and the signed half-path count), THM-5136 (the half-size tournament for invariant-cycle parities), Session132 `converse_folding/THEORY.md` (the paired-tournament carrier), and THM-4757 (linear-forest expansion for the ordinary Hamilton count). The additions here are the **complete Walsh coefficient formula and its mod-four response**, and the path-slot operator for the known even binary family. No literature-priority claim.

## 1. The signed fold and its exact turn expansion (known input)

Let `rho` be an anti-involution of a tournament `T`: `a → b` iff
`rho(b) → rho(a)`. Its `m ≥ 1` two-element orbits are called pairs. Label
their two choices by signs, arbitrarily. There is at most one fixed vertex
`c`, present precisely when the order is odd. Write `F(T,rho)` for the
number of Hamiltonian paths fixed by reversal followed by `rho`.

A fixed path consists of a signed path through all m pairs, then the
center if present, then the reverse reflected signed path. The half-path
must end at a specified sign `b_i` of its last pair `i`: the sign pointing
to the center at odd order, or the internal winner at even order.

For distinct pairs `i,j`, the two-by-two indicator matrix for arcs from
`i` to `j` has form `[[a,d],[1-d,1-a]]`. Thus its rows are either
constant and different, or identical and nonconstant. Exactly one pair
owns the relation: one sign of that pair beats both signs of the other
pair and its opposite loses to both. Let `D` be the tournament on pairs
with `i → j` precisely when `i` owns the relation. Write `x_ij ∈ {−1,1}`
for the sign in owner `i` beating both signs of `j`. The `x_ij` for
`i → j` in `D` and the `b_i` are completely independent coordinates
for any fixed labelled `D`: prescribe these arcs, and at odd order one
may orient the internal arcs separately. The owner relation is the
converse of the convention for `Z` in THM-5136.

Take any order of the pairs. Call a step backward if its owner is its
second pair. A forward step constrains the sign of its first pair to the
owner's preferred sign, and a backward step constrains the sign of its
second pair to the opposite of its preferred sign. Each maximal run of
backward steps creates one unconstrained sign at its start and exactly
one doubled constraint at its end; all remaining signs are constrained
once. If a run ends at `j,i,k` (consecutive pairs, with `i` doubly
constrained), its compatibility demands that `−x_ij = x_ik`. If it ends
at the last pair `i` after `j`, compatibility demands `−x_ij = b_i`.
The factor from the free sign and its corresponding compatibility is
therefore respectively `1 − x_ij x_ik` or `1 − x_ij b_i`.
Consequently **the contribution of the order is the product of these
factors over its backward runs**, with empty product one. Sum over
the `m!` orders to obtain `F`.

In particular `F = H(D) (mod 2)`: each nonempty run makes the
contribution even; the orders with no backward step are exactly the
directed Hamiltonian paths of `D`. This gives a constructive lift for
the parity reduction (each directed path of D has a unique completion),
not a claim of an exact ordinary-tournament quotient: some other orders
have two or more completions. Taking all `x_ij=b_i=+1` kills every
backward-run contribution and yields equality `F=H(D)`.

## 2. Exact Walsh coefficients: alternating linear forests

Expand F in squarefree monomials in the signs `x_ij,b_i`. A monomial can
occur only if, at each marked owner, it selects either two distinct
outgoing edge variables, or one outgoing edge variable and its terminal
variable `b_i`. There can be at most one marked terminal. Put an
undirected edge in a graph G for every selected edge variable and
retain all m vertices. The coefficient is zero unless G is a linear
forest (a disjoint union of paths and isolated vertices), and if a
terminal is marked its owner has degree one in G.

In the allowed case let `p` be the number of nontrivial components,
`c` the number of all components, and `s` the number of marked owners.
The coefficient is

```text
(-1)^s 2^p c!                  if no terminal is marked;
(-1)^s 2^(p-1) (c-1)!          if a terminal is marked.
```

**Proof.** A chosen turn factor enforces adjacency to the two selected
neighbors of its marked owner (or enforces adjacency to the selected
neighbor and the end of the order). These adjacencies force G to be a
subgraph of the linear order. Conversely, turn each component of G into
a contiguous block. Every marked interior owner is automatically a
backward-to-forward turn in either orientation of its block: it owns
both incident edges. The terminal block, if present, must be last and
oriented with the marked owner at its end. Order and orient all other
blocks freely. Every compatible order contributes `(-1)^s` once to
this particular monomial and no other orders contribute. ∎

Consequences:

- For independent fair signs, `E[F | D] = m!` for **every** owner
  tournament D. The dependence on D resides entirely in nonconstant
  alternating-forest terms. In particular the coefficient on
  `x_ij b_i` is `−(m−2)!`, and on `x_ij x_ik` it is
  `−2(m−2)!` when those edges exist.
- For `m≥2`, the only **odd nonconstant** coefficients are terminal
  monomials whose G has just one nontrivial component, an alternating
  path on `2 floor(m/2)` vertices, with owners at the second, fourth,
  ..., last positions when read toward the terminal. In fact without
  a terminal there is a factor two. With a terminal there are
  `2s−1` edges, so `c−1=m−2s`, and the coefficient is odd
  exactly when `p=1` and `m−2s∈{0,1}`.
- Hence F modulo four is an affine function of which sign coordinates
  differ from the all-positive assignment. Its constant term is
  `H(D)` modulo four. Each sign contributes twice the parity of its
  incidence among the alternating paths in the preceding bullet
  (with a terminal variable included at each path's last vertex).
  Indeed flipping a collection of signs changes a monomial by
  `-2 times the number of flipped coordinates in its support (mod 4)`;
  even coefficients disappear, and the sign of an odd coefficient
  does not matter. In particular, in the absence of such alternating
  paths, the fixed count is constant modulo four even though it need
  not be constant as an integer.

## 3. The known exact binary construction at both parities

**Session29 construction (known).** For every m≥1 and every odd `f` with `1≤f≤2^m−1`,
there is an explicit self-converse tournament on **2m** vertices and
one on **2m+1** vertices with fixed Hamiltonian-path count f. Therefore
the fixed-path spectrum over all self-converse tournaments is precisely
the positive odd integers, including 7 and 21.

**Proof.** Use one pair (and a fixed center if requested) at the outset.
The number F of admissible signed half-paths is 1. Let U count signed
half-paths if the terminal sign requirement is discarded. Initially
U=2. Add a new source q in the owner tournament, whose positive sign
beats both signs of every older pair. The new pair cannot be internal:
a step into q demands its negative sign and a step out demands its
positive sign. Without a terminal condition, either placement (first
or last) contributes U, so `U↦2U`. With the terminal condition, placing
q first contributes the old F; placing q last contributes U iff its
terminal sign is negative. This sign is free: select the direction of
the internal arc in even order, or the arc to the center in odd order.
Thus `F↦F` or `F↦F+U`. After m pairs, the value is 1 plus an arbitrary
subset sum of `2,4,…,2^(m−1)`, exactly the claimed interval of odd
values. The impossibility of even F is the parity reduction and Rédei.
∎

In these witnesses the fixed-path count is also the blue-tiling count: the controller is transitive and the terminal coordinates fix the gauge, so the anti-involution centralizer is trivial (Session29 §11). In the odd construction one may orient each pair internally positive
to negative; deleting the fixed center then leaves a transitive
tournament. In the even construction the only backward arcs relative
to the outside-in linear order are precisely selected antipodal
internal arcs. These facts follow from the rule that each newest
positive vertex beats the entire older body and that body beats the
newest negative vertex.

## 4. Total count in the even binary family: a central difference

Let `P_T(q)` count assignments of the vertices of T to q labelled
path slots, allowing empty slots, together with one directed Hamiltonian
path in each nonempty slot. Thus `P_T(1)=H(T)`; this is the existing
path-cover profile. Append outer vertices `x,y` with `x → T → y`
and write b=0 for `x → y`, b=1 for `y → x`. For positive q,

```text
P_new(q) = q² P_T(q) + b q (P_T(q+1) − P_T(q−1)).
```

For b=0, independently choose the slots of x and y. When b=1,
different slots contribute `q(q−1)P_T(q)`. In the same slot, paths
with x before y require a nonempty old path between them and contribute
`q(P_T(q)−P_T(q−1))`. Paths with y before x must use `y → x`
consecutively: split their old portion into the two slots before y
and after x, giving `qP_T(q+1)`. Summing proves the formula.
The base profile with one pair is `q²`; m pairs are thus an explicit
binary word of operators applied to this polynomial. In particular a
single active bit added outside a transitive body with `2k` vertices
gives `(F,H)=(1+2^k,1+4^k)`. This is a concrete bridge, not a claim
that every prescribed F can be made equal to H.

**Boundary for Chan–Pak.** For every anti-involution, `H = F + 2R`,
where R is the number of two-element orbits of paths. Since both 7 and
21 occur for F, no obstruction based only on the nonexistence of a
particular fixed count can explain the two exceptional H-values.
The pairing term, or a constraint coupling it to F, is essential.

## Reproduction

Run the accompanying standard-library checker under Python both normally
and with `-O`. It computes fixed counts from actual cross arcs and explicit
signed permutations, independently obtains all Walsh coefficients by a
fast Walsh transform (6,674 sign assignments), and checks the slot transfer
by direct permutation and subset partition counts for 14 binary seams
at three positive slot counts. These finite controls do not establish the
all-order scope; the proofs do.
