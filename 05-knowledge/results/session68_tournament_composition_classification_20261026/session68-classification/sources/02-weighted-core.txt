# Counting tournaments through weighted pair-free cores

**Status:** PROVED structural decompositions and formal identities below;
FINITE-EXACT class reconstruction through order seven and weight-orbit
checks through total weight twelve. No closed scalar recurrence for
A000568 is claimed.

## Inheritance, objects, and the missing coordinate

[THM-4463 — homogeneous tournament contacts](../../../01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md)
proves that the graph of homogeneous pairs is a disjoint union of directed
paths, whose components are the unique maximal transitive modules.
[The preceding inflation audit](../metagraph_growth_20261020/growth.md)
uses their lengths to recover uniform transitive roots. We now retain the
whole length vector rather than only its gcd.

The closest hostile is the preceding
[order-seven local-data collision](../metagraph_growth_20261020/enumeration.md):
even ten local numerical marks do not determine metagraph line partners.
The corrected near miss here is that the **number** of vertices and even
the **order** of the automorphism group do not determine the number of
weightings. The needed sidecar for weight counts is the permutation cycle
index, not merely group size.

The concept board is: **canonical pair-free quotient; positive chain
weights; automorphism orbits; gcd scale; fixed-point-free symmetry boundary**.
The map from a tournament to a core and weights preserves its entire
isomorphism class. Projecting further to a cycle index preserves weight
counts but loses the actual arcs needed for class-specific transport.

Here “pair-free” means having no homogeneous pair. It does not mean
indecomposable: a pair-free tournament may have nontransitive modules.
It also differs from gcd-primitive, which refers only to the chain lengths.

## 1. The exact weighted-core classification

For a pair-free tournament `Q` on `q` vertices and a positive integer
vector `w=(w_1,...,w_q)`, let `Q[w]` replace vertex `i` by a transitive
chain of length `w_i`.

**Proposition.** Every nonempty tournament is uniquely of the form `Q[w]`,
up to an isomorphism of `Q` carrying its weights. For a fixed representative
`Q`, its distinct classes are precisely

```text
{positive integer vectors w}/Aut(Q).
```

**Proof.** Contract the canonical maximal transitive modules. The quotient
has no homogeneous pair, since two homogeneous quotient vertices would
give a larger transitive module. Conversely, inflating a pair-free quotient
cannot create a homogeneous pair between different inflated chains:
such a pair would imply a homogeneous pair in the quotient. Thus the
chains are exactly the canonical maximal transitive modules.
An isomorphism must permute these modules, preserve their lengths, and
preserve the quotient arcs. Inside each transitive chain it has a unique
possible order-preserving action. This proves uniqueness and the orbit
description. It also proves

```text
Aut(Q[w]) = Stab_(Aut Q)(w).                              (1)
```

Let `W_Q(n)` be the number of positive weight orbits of total weight `n`.
Writing `C_q` for the set of pair-free core classes of order `q`,

```text
A(n) = sum_(q<=n) sum_(Q in C_q) W_Q(n).                 (2)
```

This is a class bijection, not an approximate or average count.

### A constructive recurrence without orbit averaging

Fix a canonical representative of every `Aut(Q)` orbit of weight vectors.
Let `O_Q(n)` be the set of these representatives with sum `n`. Then

```text
O_Q(q) = {(1,...,1)},
O_Q(n+1) = {canon_G(w+e_i): w in O_Q(n), 1<=i<=q}.       (3)
```

The braces mean a set, with duplicates removed. Every positive vector of
sum greater than `q` has a coordinate that can be decreased, proving
surjectivity of this construction. Hence (3) is an exact constructive
recurrence on weight orbits. It is not a recurrence on the number
`|O_Q(n)|` alone, because different increments can reach the same orbit.

## 2. Rational response of a fixed core

Let `G=Aut(Q)`. If a permutation `g` has cycles of lengths
`d_1,...,d_s`, a fixed weight vector is constant on each cycle. Its positive
weight generating function is

```text
x^q / product_j(1-x^d_j).
```

Counting fixed weight vectors and dividing by `|G|` therefore gives

```text
F_Q(x) := sum_n W_Q(n)x^n
       = x^q/|G| * sum_(g in G) product_(cycles c of g)
                                      (1-x^|c|)^(-1).    (4)
```

This alternate closed form is explicitly **a local Burnside calculation
inside the core automorphism group**. The implementation uses (3), and
uses (4) as an independent check. It would be misleading to call (4)
Burnside-free.

Every fixed core consequently has a rational generating function and an
eventually quasipolynomial weight response. For example, if `L` is the lcm
of the cycle lengths occurring in `G`, its denominator divides
`(1-x^L)^q`. Multiplication by this denominator gives a polynomial of
degree at most `qL`, so

```text
W_Q(n) = sum_(j=1..q) (-1)^(j+1) binom(q,j) W_Q(n-jL),
                                                   n>qL. (5)
```

This uniform but usually nonminimal recurrence concerns one fixed core.
The total sum (2) introduces new cores at every order; it is not thereby
shown rational.

For the cyclic triangle,

```text
W_C3(n) = [binom(n-1,2)+2*[3 divides n]]/3,  n>=3.
```

The exceptional equal-weight orbit is the source of the divisibility term.

## 3. Exact finite core bank and a minimal symmetry obstruction

The independent weighted-core construction reproduces these class counts:

| Total order | Contributions by core order | Total `A(n)` | Pair-free cores |
|---:|---|---:|---:|
| 1 | `1:1` | 1 | 1 |
| 2 | `1:1` | 1 | 0 |
| 3 | `1:1, 3:1` | 2 | 1 |
| 4 | `1:1, 3:1, 4:2` | 4 | 2 |
| 5 | `1:1, 3:2, 4:4, 5:5` | 12 | 5 |
| 6 | `1:1, 3:4, 4:8, 5:17, 6:26` | 56 | 26 |
| 7 | `1:1, 3:5, 4:16, 5:47, 6:132, 7:255` | 456 | 255 |

The 255 order-seven cores occupy seven automorphism-cycle-index types,
including 224 rigid cores. This is a real compression for **weight
enumeration**: all cores with the same cycle index have the same function
(4). It is not a claim that they have identical extension or line data.

Order and group order alone fail already for two order-six pair-free
cores with adjacency row masks

```text
Q_A = (0,33,11,19,7,29),
Q_B = (10,20,33,38,13,19).
```

Bit `j` in row `i` means `i -> j`. Both have automorphism group of order
three. The two nonidentity automorphisms have respective cycle types

```text
Q_A: (3,1,1,1),
Q_B: (3,3).
```

Thus their rational responses are

```text
F_A(x) = x^6/3 * [(1-x)^(-6)
                     +2/((1-x)^3*(1-x^3))],

F_B(x) = x^6/3 * [(1-x)^(-6)+2/(1-x^3)^2].
```

Their counts at total weights six through twelve are respectively

```text
1,4,11,26,54,102,180,
1,2, 7,20,42, 84,156.
```

The first separation has a simple mechanism: at total weight seven,
exactly one vertex receives weight two. Its orbit is a vertex orbit of
the core. The two groups have four and two vertex orbits. The script
checks that no same-order, same-group-order response collision occurs
at smaller core orders.

The bank of all cores through order seven generates exactly 2,278
classes of total order eight. This is a **core-size truncation**, not
the full order-eight count. Newly appearing pair-free order-eight cores
are excluded by construction.

## 4. Primitive scale, and a number-theoretic example

The gcd of the weights is invariant under `G`. Let `P_Q(n)` count
weight orbits of total `n` with gcd one. Dividing the weight vector by
its gcd is an equivariant bijection, so

```text
W_Q(n) = sum_(d|n) P_Q(n/d),
P_Q(n) = sum_(d|n) mu(d) W_Q(n/d).                       (6)
```

In generating-function language,

```text
F_Q(x)=sum_(d>=1) P_Q(x^d),
P_Q(x)=sum_(d>=1) mu(d) F_Q(x^d).
```

These are coefficientwise formal sums. Rationality of `F_Q` does not imply
rationality of the primitive series.

For the cyclic triangle, count ordered primitive positive triples by
Möbius inversion, then divide by cyclic rotations. Nonidentity rotations
fix a primitive triple only at `(1,1,1)`. Consequently, for `n>1`,

```text
P_C3(n) = [J_2(n)-3*phi(n)]/6 + (2/3)*[n=3],             (7)
```

with `P_C3(1)=0` and `J_2(n)=sum_(d|n) mu(d)(n/d)^2`.
The fractional-looking terms add to an integer, including at `n=3`.
This gives an exact arithmetic scale law for one family of weighted
tournament cores. It does not provide a recurrence for newly appearing
pair-free cores. Direct primitive composition orbits through total sixty
independently verify (7).

## 5. Converse and merged counts also have a local response

An inflation can be self-converse only when its canonical core is
self-converse. For such a core, let `G^-` be its anti-automorphism coset.
It has `|G|` elements. The number `S_Q(n)` of weight orbits fixed by
converse satisfies the twisted orbit formula

```text
sum_n S_Q(n)x^n
 = x^q/|G| * sum_(h in G^-) product_(cycles c of h)
                                         (1-x^|c|)^(-1). (8)
```

One proof counts pairs `(weight vector, anti-automorphism fixing it)`.
Each converse-fixed `G`-orbit contributes exactly `|G|` such pairs;
other orbits contribute none. This also follows by taking the difference
between ordinary orbit counting for `G` and for `G union G^-`.

Thus the marked core bank computes `S(n)=sum_Q S_Q(n)` and
`M(n)=(A(n)+S(n))/2`. Direct weight-orbit testing and the anti-coset
formula independently reconstruct every `S(n),M(n)` through order seven.
Blue-line distributions require additional path information; (8) does
not claim to reconstruct them.

## 6. A formal recurrence with an explicit free boundary

Define cycle-index series from actual automorphism actions:

```text
T(p) = sum_[all tournament classes R] (1/|Aut R|)
                     sum_(g in Aut R) product_i p_i^m_i(g),
C(p) = sum_[pair-free core classes Q] (1/|Aut Q|)
                     sum_(g in Aut Q) product_i p_i^m_i(g).
```

Include the empty object as constant one. Variable `p_i` has weight `i`.
All identities here are **formal**: the ordinary and exponential
tournament generating series have no positive analytic radius of
convergence.

A core automorphism cycle of length `i`, with equal chain weight `l`
on its vertices, induces `l` vertex cycles of length `i` on the inflated
tournament. Summing `l>=1` replaces `p_i` by `p_i/(1-p_i)`. Therefore

```text
T(p)=C(p_i/(1-p_i)),
C(p)=T(p_i/(1+p_i)).                                    (9)
```

Distinguishing one vertex of a tournament leaves a smaller tournament
and a binary incidence choice on every vertex. An automorphism fixes
those choices precisely when they are constant on its cycles. Hence

```text
partial_1 T(p)=T(2p_1,2p_2,...).                         (10)
```

Substitute (9) into (10), then replace `p_i` by `q_i/(1+q_i)`. The chain
rule gives the exact core equation

```text
(1+q_1)^2 partial_1 C(q)
       = C(2q_i/(1-q_i)).                               (11)
```

This is an all-order structural relation. Its boundary
`C(0,q_2,q_3,...)` remains free from the viewpoint of this equation:
it records **fixed-point-free automorphisms of cores**. The derivative
does not determine those terms.

More explicitly, write `C=sum_(r>=0) q_1^r C_r(z)`, where
`z=(q_2,q_3,...)`, and let `H(z)_i=2z_i/(1-z_i)`. Then

```text
C_1(z)=C_0(H(z)),

(r+1) C_(r+1)(z)
 = sum_(j=1..r) 2^j binom(r-1,j-1) C_j(H(z))
     -2r C_r(z)-(r-1) C_(r-1)(z),             r>=1.       (12)
```

Once the fixed-point-free boundary is given, every fixed-point tail is
determined. This recovers the same boundary distinction as
[THM-514 — the A000568 unit-tail core kernel](../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md),
but now after removing canonical transitive modules. These two meanings
of “core” must not be conflated: THM-514's core is a permutation cycle
partition without `1`-parts; ours is an actual pair-free tournament.

The computed boundary through weighted degree seven is

```text
C(0,p_2,p_3,...)
 = 1+(2/3)p_3+(4/5)p_5+(10/9)p_3^2+(8/7)p_7+O_weight(8).
```

For comparison, the coefficient of `p_3^2` in `T` is `16/9`.
The missing `2/3` comes from a clone-bearing family rather than a new
pair-free core. All coefficients of (9) are independently checked
through weighted degree seven, and (10)--(11) through degree six.

A useful partially explicit boundary survives. If an automorphism has
all its cycle lengths distinct, every canonical chain has length one:
a chain orbit of length `i` and chain weight `l>1` would give `l` cycles
of the same length. Thus `C` and `T` have equal coefficients on squarefree
cycle monomials. In particular, an invariant tournament under one odd
`k`-cycle is automatically pair-free, and

```text
[p_k] C = [p_k] T = 2^((k-1)/2)/k,       k odd.
```

Repeated cycle lengths are the only place clone corrections can enter;
the `p_3^2` example is the first such correction.

## 7. Closing the construction on core objects

Equation (2) alone leaves the new core bank `C_n` unspecified. There is
an exact constructive augmentation using the same weighted state.
Take any parent `Q[w]` of order `n`. To obtain a pair-free child by adding
one vertex, its adjacency must **alternate along every parent chain**:
otherwise an old adjacent homogeneous pair survives. There are only
`2^q` starting-phase vectors. Retain those producing no homogeneous
pair involving the new vertex, and canonicalize their child tournaments.
Every pair-free child is reached by deleting any one of its vertices,
so this generates the entire next core bank.

Together, weighted-orbit increment (3) and this filtered alternating
augmentation give a closed recurrence on **finite sets of exact core
objects and weights**. They do not require Burnside's sum over `S_n`.
They still require isomorphism handling and new core construction.

Even the exact count can be written without deduplicating children if
the appropriate symmetry weights are retained. For `N=n+1`,

```text
|C_N| = (1/N) sum_[parents P of order n]
               (1/|Aut P|) sum_[admissible alternating masks e]
                                            |Aut(P+e)|. (13)
```

For a fixed child `U` and parent type `P`, let `d_U(P)` be the number of
vertices `v` for which `U-v` has type `P`. The number of extension masks
on one representative `P` producing `U` is
`d_U(P)|Aut P|/|Aut U|`: count pairs consisting of a deleted vertex and
an isomorphism of its complement with `P`. Consequently `U` contributes
`sum_P d_U(P)/N=1` to (13).

The automorphism factor is essential. The transitive two-vertex parent
has two alternating candidates, but only one is a directed triangle;
the other creates new homogeneous pairs and is rejected. Dividing its
one admissible mask by three gives `1/3`; multiplying by the triangle's
automorphism order gives the correct one core.

This is the exact boundary of the progress: the structural recurrence
is closed, whereas a recurrence involving only previous numerical class
counts is not established.

## Reproduction

```text
python 05-knowledge/results/weighted_core_transport_20261020/core_count.py
python -O 05-knowledge/results/weighted_core_transport_20261020/core_count.py
```

The script imports the prior census and chain-decomposition modules from
the committed `metagraph_growth_20261020` folder. It independently builds
weight orbits by coordinate increments, compares them with direct
composition orbits and the cycle formula, and reconstructs the complete
class sets through order seven. It tests 208,126 weight orbits through
total twelve across the core bank, primitive gcd inversion, anti-coset
counts, and formal cycle-index identities. Explicit exceptions keep every
check active under optimization. The companion `.out` is literal output.
