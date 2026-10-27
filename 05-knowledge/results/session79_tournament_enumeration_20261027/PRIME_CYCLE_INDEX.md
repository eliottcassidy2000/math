# Prime tournament cycle indices can be recovered, not enumerated

**Status:** PROVED deduction from classical decomposition and cycle-index identities;
FINITE-EXACT verification through degree 15, independent exhaustive verification through 6.
No claim of literature novelty. Independent audits are described below.

## Inheritance and exact improvement

The modular composition identity is inherited from the earlier prime-family
work and its [corrected composition review](../tournament_composition_classification_20261026/README.md).
The former supplies a prime census through order 9 and prime counts at order 10.
The latter corrects the false claim that unsigned cyclic triples preserve
automorphism order. None of those corrected-away claims is used here. A finite
prime census alone leaves higher-order prime contributions undetermined.

The improvement here is an exact **aggregate prime cycle-index recovery** from the
known tournament cycle index, avoiding a prime-class census altogether. It does
not reconstruct prime isomorphism representatives, individual automorphism groups,
or the self-converse/anti-automorphism data used by the intransitivity quotient.
It is circular if advertised as an independent way to discover A000568: its
input already encodes Burnside's answer. Its new useful output is the structural
prime inventory and its exact aggregate permutation-action sidecar.

Concept board: ordered SCCs; canonical strong modules; prime quotient action;
Burnside odd partitions; plethystic inversion. Hostile: treating a prime quotient
as an unordered set of slots, or using ordinary univariate composition/inversion.

## 1. Conventions and the classical structural input

Let T include the empty tournament, A=T minus the empty structure, S be nonempty
strongly connected tournaments (the singleton is strong), and P be prime tournaments
of order at least 3. A module is a vertex set seen uniformly by each outside vertex;
it is strong if it overlaps no module (overlap means intersection and both differences
are nonempty). Prime means no module of size 2 through n−1. In particular, prime
here excludes the singleton and the two-vertex tournament.

The classical Gallai theorem says that a non-singleton tournament's maximal proper
strong modules form a canonical partition. The quotient is either transitive,
with the blocks its ordered strong components, or prime of order at least 3.
See the cited Gallai theorem, Theorem 16 in
[Boussaïri–Chergui–Ille–Zaidi (2018)](https://arxiv.org/html/1805.04607v1).
This is a cited structural dependency, not an independent proof of Gallai here.

Consequently, as species,

    T = SEQ(S),             S = X + P ∘ A,             A = T−1.       (1)

Why no overcount: an isomorphism preserves the canonical strong-module partition,
and hence the quotient and the isomorphism types attached to its vertices. At a
prime root, two assignments give the same isotype precisely when an automorphism
of the prime quotient carries one assignment to the other. At a transitive root,
the order is intrinsic and has no nontrivial automorphism. There is no arbitrary
parenthesization of consecutive linear factors: SEQ uses the maximal SCCs.

Conversely every inflation of a prime Q of order at least 3 really has these
blocks as its canonical root blocks. To see this, any module meeting two blocks
has a support (the set of blocks it meets) that is a module of Q. Primeness forces
that support to be all of Q. If this module omits a vertex of a block B_i, that
vertex must see the module uniformly, forcing vertex i of Q to dominate or lose
to every other quotient vertex. That is impossible in a prime tournament of
order at least 3. So every proper module is confined to one block, and the blocks
are exactly the maximal proper strong modules. Also Q is strongly connected:
a nontrivial SCC order in an order-at-least-3 prime would yield a proper
nonsingleton module. Paths around the strong quotient connect every pair of
vertices of its inflation, even when they lie in the same child block.

## 2. Exactly what sidecar is sufficient

For each prime isotype Q on k vertices write

    Z_Q(p1,p2,...) = (1/|Aut Q|) Σ_g∈Aut Q ∏_d p_d^{c_d(g)},
    Z_P = Σ_[Q] Z_Q.

All cycles of a tournament automorphism have odd length: an even cycle contains
an antipodal pair whose orientation is reversed by an appropriate power.
Only odd p_d occur. The **sum Z_P** is sufficient; separate Qs, separate groups,
and a histogram of group types are unnecessary for unrefined composition counts.
Group order alone, however, is not generally a substitute for its vertex action.

Let t(x)=Σ t_n x^n and s(x)=Σ s_n x^n be ordinary isotype generating functions.
Burnside on the prime root gives exactly

    t(x) = 1/(1−s(x)),
    s(x) = x + Z_P(t(x)−1,t(x²)−1,t(x³)−1,...).                 (2)

Proof of the root term: an assignment fixed by g is constant on every d-cycle
of g. Choosing a nonempty child type of size m on this cycle contributes x^(dm),
so the cycle contributes t(x^d)−1. Multiplication handles different cycles;
averaging counts assignment orbits. Internal child automorphisms do not create
additional choices because assignments already use child isomorphism types.

The scalar prime count series is merely Z_P(x,x²,...). This specialization loses
the slot permutation action before (2) can use it. Concrete hostile: the three
five-vertex prime tournaments have automorphism orders 1,1,5, with vertex-orbit
counts 5,5,1. Replacing exactly one quotient vertex by the unique two-vertex
tournament gives 5,5,1 inflation types, respectively, not the uniform 5 each.
These are verified without literature assumptions in `verify_prime_cycle_index.py`.
This proves insufficiency of replacing each root's action by its size/count;
it is not a formal impossibility theorem for all imaginable manipulations of the
one fixed sequence of prime counts.

## 3. Recovery theorem in the completed cycle-index ring

Work over Q in power-sum variables p_i with weight i. Let

    B = Z_T,       F = B−1.

Products are ordinary multiplication of symmetric functions. Plethysm U∘V means
replacing p_d in U by V(p_d,p_(2d),p_(3d),...). It is defined for V with zero
constant term. Since F=p1+terms of weight ≥2, F has a unique plethystic inverse
F^⟨−1⟩=p1+terms of weight ≥2.

**Theorem (aggregate prime recovery).**

    Z_P ∘ F = 1−B^(−1)−p1,                                  (3)
    Z_P = p1/(1+p1) − F^⟨−1⟩.                                (4)

Proof: the cycle index of SEQ is 1/(1−p1). Apply the standard cycle-index
substitution rule to (1). This gives B=1/(1−Z_S), and
Z_S=p1+Z_P∘F, proving (3). Right-compose (3) by F^⟨−1⟩. The right side is
F/(1+F)−p1, which becomes p1/(1+p1)−F^⟨−1⟩. This proves (4).
Existence/uniqueness can also be seen directly by weighted-degree induction:
plethystic substitution of F into a homogeneous degree-n polynomial preserves
that polynomial in degree n. Therefore

    (Z_P)_n = [1−B^(−1)−p1 − (Z_P)_(<n) ∘ F]_n.              (5)

Equation (5), used by the script, requires no symbolic inverse routine.
The universal cycle-index substitution rule is classical; an explicit statement
is Theorem 2.5 of [*Enumeration of bipartite graphs and bipartite blocks* (2014)](https://emis.de/ft/31891).
The elementary fixed-assignment proof in section 2 explains the scalar specialization.

**Crucial warning:** F^⟨−1⟩ is a multivariate plethystic inverse, **not** the ordinary
compositional inverse of t(x)−1. The two operations do not commute with the
unlabelled specialization p_i=x^i.

## 4. Explicit input: no tournament-class list is necessary

For an odd partition λ=(l1,...,lk), put

    e(λ)=Σ_i (li−1)/2 + Σ_(i<j) gcd(li,lj),
    z_λ=∏_d d^(m_d) m_d!,       p_λ=∏_i p_(li).

Then the classical tournament Burnside cycle index is

    B = Σ_(λ all odd, including empty) 2^e(λ) p_λ/z_λ.        (6)

Self-contained fixed-point explanation: an even vertex cycle forbids an invariant
orientation, as above. For odd cycles, no pair-orbit contains an orientation
reversal, so each unordered-pair orbit has two independent orientations. Within
a vertex cycle of length l there are (l−1)/2 pair-orbits; between cycles of lengths
l,m there are gcd(l,m). Finally, cycle-type λ has n!/z_λ permutations. This proves
(6) directly from Burnside. The same formula appears in the repository's
THM-514, [THM-514](../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md).

Thus (5)+(6) recover Z_P to any finite requested weight using only exact rational
arithmetic indexed by odd partitions. This is partition-scale storage rather
than isomorphism-class-scale storage; no runtime complexity bound is asserted.

### Sparse support theorem: coprime-isolated cycles are forbidden

**PROVED.** Let g be an automorphism of a prime tournament, and let C be a
nontrivial vertex cycle of g. If C is not the whole vertex set, some other cycle
D of g must satisfy gcd(|C|,|D|)>1.

Proof: suppose gcd(|C|,|D|)=1 for every other cycle D. For any outside vertex
v in a cycle of length m, g^m fixes v and acts transitively on C, because m is
coprime to |C|. Thus v sees every vertex of C in the same direction. This holds
for every outside vertex, making C a proper nonsingleton module, contradiction.

Equivalently, in a nonzero monomial p_λ of Z_P, every part d>1 must either be
the sole part or have another part e with gcd(d,e)>1. This is a necessary
support condition, not sufficient for positive coefficient. Identity types
1^n and whole-set odd cycles n are not excluded. It strengthens the simpler
single-nontrivial-cycle-plus-fixed-points exclusion and rules out types such as
(3,5), (1,3,5), and (3,3,5), without examining any tournament. All recovered
coefficients through 15 pass the stronger filter. A targeted search in the
inherited Session67 packet found no statement of this filter; no literature
novelty is asserted.

## 5. Exact evidence and usable small bank

    (Z_P)_3 = (p1³+2p3)/3
    (Z_P)_4 = 0
    (Z_P)_5 = (11p1⁵+4p5)/5
    (Z_P)_6 = (43p1⁶+2p3²)/3
    (Z_P)_7 = (4057p1⁷+56p1 p3²+24p7)/21
    (Z_P)_8 = 4000p1⁸+8p1²p3²
    (Z_P)_9 = (1229279p1⁹+432p1³p3²+172p3³+12p9)/9
    (Z_P)_10 = (118583443/15)p1¹⁰ + (1400/3)p1⁴p3²
               + (368/3)p1p3³ + (44/5)p5²

Summing coefficients of each homogeneous part gives the unlabelled prime counts.
The exact outputs are:

| n | primes, excluding orders 1 and 2 |
|---|---:|
| 3 | 1 |
| 4 | 0 |
| 5 | 3 |
| 6 | 15 |
| 7 | 197 |
| 8 | 4008 |
| 9 | 136655 |
| 10 | 7906161 |
| 11 | 796876429 |
| 12 | 142980435791 |
| 13 | 46453867277461 |
| 14 | 27687212017612715 |
| 15 | 30571668704134867685 |

The entire aggregate bank through 15 has 55 nonzero monomials and is in
the output of `recover_prime_cycle_index.py 15`, including the recovered degree-10 aggregate data.
This statement does not supply individual prime classes or anti-automorphisms.

Reproduce from this directory:

```sh
python3 verify_prime_cycle_index.py
python3 recover_prime_cycle_index.py 15
python3 -O recover_prime_cycle_index.py 15
```

The first script exhausts all labelled orientations through 6, marks their entire
S_n orbits, tests every possible nonsingleton proper module, enumerates all
automorphisms, and independently checks (2) against all isotype/strong counts.
The second script checks (3) coefficientwise through 15, nonnegative coefficients,
integral prime counts, the coprime-isolated-cycle forbidden support filter, agreement of the **whole
cycle index** with that independent exhaustive census through 6, and agreement
of prime counts with the prior stored census through 10. Values at 11–15 are
exact outputs of the proved recurrence, not independent brute-force censuses.
Both scripts use explicit RuntimeError gates rather than Python `assert`, so
verification remains active under `python3 -O`. Normal and optimized runs were
repeated and their JSON outputs compared byte-for-byte, with equality. Recovery
runs the companion exhaustive checker itself and does not depend on a prior
receipt file.

## 6. A closed scalar boundary: primes fixed by a specified full cycle

For odd n define f_n=2^((n−1)/2), including f_1=1. This counts labelled tournaments
on Z/nZ fixed by the **specified** rotation v↦v+1: choose one direction for each
pair of nonzero opposite differences {d,−d}. Let g_n count prime tournaments fixed
by that same rotation, with g_1=0 under the order-at-least-3 convention.

**PROVED corollary.** On odd positive integers,

    g_n = f_n − Σ_(d|n, 1<d<n) g_d f_(n/d),          n>1,
    g = δ − f^(−1)                                   (Dirichlet inverse).

Proof: the coefficient of the lone monomial p_n in B is b_n=f_n/n. The
coefficient of p_n in Z_P is c_n=g_n/n, since z_(n)=n in the definition of a
species cycle index. For n>1, the coefficient of lone p_n in 1−1/B−p1 is b_n:
writing B=1+F, every power F^k with k≥2 contains at least two power-sum factors.
A lone p_n on the left side of (3) can arise only from an outer lone p_d and
an inner lone p_(n/d), so

    b_n = Σ_(d|n,d≥3) c_d b_(n/d).

Multiplying by n gives f_n=Σ_(d|n,d>1)g_d f_(n/d); isolate d=n for the recurrence.
At n=1, g_1=0, so the convolution identity at all odd n is g*f=f−δ. Since
f_1=1 it has the stated Dirichlet inverse solution.

Examples: g_9=16−2·2=12; g_15=128−2·4−4·2=112. For every odd prime p there
is no proper divisor term, so **every tournament fixed by a specified p-cycle
is prime**, and g_p=2^((p−1)/2). This conclusion concerns modular primeness, not
merely the arithmetic primality of the vertex count.

Independent hostile and positive controls: `verify_prime_cycle_index.py` now enumerates every
one of the 2^((n−1)/2) difference-sign assignments for each odd n=3,...,15 and
tests every proper nonsingleton vertex subset for being a module. It gives

| n | all fixed tournaments f_n | prime fixed tournaments g_n |
|---|---:|---:|
| 3 | 2 | 2 |
| 5 | 4 | 4 |
| 7 | 8 | 8 |
| 9 | 16 | 12 |
| 11 | 32 | 32 |
| 13 | 64 | 64 |
| 15 | 128 | 112 |

The recovery script checks these direct counts both against the divisor
recurrence and against n[p_n]Z_P. In particular, the composite-order hostile
shows that rotational symmetry does not itself force modular primeness.

**Scope:** these are labelled fixed points of one specified full cycle, not
unlabelled circulant isotypes, all vertex-transitive prime isotypes, or all
prime tournament isotypes. No division by n, Euler φ(n), or n! on its own
converts them into those orbit counts. Literature on unlabelled circulant
enumeration includes Klin–Liskovets–Pöschel, *Analytical Enumeration of Circulant
Graphs with Prime-Squared Number of Vertices*:
https://www.mat.univie.ac.at/~slc/wpapers/s36klp.html.
The targeted search did not establish an attribution for this precise fixed-cycle
prime recurrence, so this is presented as a proved consequence of the classical
cycle-index machinery, with no novelty claim.

## 7. Metagraph relevance and limits

A canonical modular-tree construction generates tournament isotypes, hence the
vertices of the ordinary tournament metagraph. Root action orbits are the exact
duplicate-removal mechanism. The aggregate cycle index counts these vertices
without constructing them. It does not recover edge-flip incidence, adjacency,
loops, a canonical label for each class, or the C3 quotient and its anti-action.
These require additional marked data. Therefore this is an exact enumeration
and compression bridge, not a solution of the metagraph adjacency problem.


## Independent audit

A separate implementation recovered the full plethystic inverse in the opposite
direction, checked both inverse identities, and matched all 55 prime-index
coefficients through weight 15. A separate fixed-permutation orientation census
tested every permutation cycle type through six, including even-cycle zero
controls, and filtered every possible proper module. A third direct circulant
check verified the full-cycle recurrence through 15. Normal and optimized
outputs agreed. These checks supplement, rather than replace, the all-order
proof and its explicitly cited structural dependency.
