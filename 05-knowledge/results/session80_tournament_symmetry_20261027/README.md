# A cohesive calculus of tournament symmetry and composite families

Session 80, begun 2026-10-27. **PROVED deductions**, with explicitly cited
structural dependencies and independently checked finite computations. No
historical-priority or polynomial-time claim.

## Answer in one paragraph

The recovered aggregate prime cycle index generates the full composite cycle
index, canonical-tree statistics, and depth-filtered family counts without a
composite-class census. But it cannot generate every symmetry-resolved family
statistic: actual prime tournaments with the same cycle index can give different
composite automorphism-order distributions. The exact finer formula retains
child isomorphism types and subgroup stabilizers. Adding the converse action
closes a second aggregate calculus, including realizable cyclic-triple
hypergraphs. At fixed excess above a deepest prime chain, both the counting
recurrence and the possible moving symmetry actions become bounded.

These are compatible levels of information, not competing formulas.

## 1. The exact symmetry law on a canonical tree

[Proof and controls](GROUP_LAW.md) · [Checker](check_groups.py)

For a canonical prime-root decomposition T=Q[A_i], define

    H = {g in Aut(Q): A_i is isomorphic to A_(g i) for every i}.

Then the canonical exact sequence splits after coherent child identifications:

    Aut(T) ≅ (product_i Aut(A_i)) ⋊ H.

At an ordered strong-component node the quotient action is trivial, even when
components repeat. At a prime node, repeated child types permit only the actual
quotient symmetries that preserve those types. Arbitrary modular partitions
compute only their partition stabilizer, not necessarily the whole group.

If Z_A is a child's normalized permutation cycle index and ψ_d replaces p_j by
p_(dj), then

    Z_T = (1/|H|) sum_(g in H) product_(C cycle of g) ψ_|C|(Z_A_C).

The return map around a quotient cycle explains the multiplication of vertex
cycle lengths. In the inherited twelve-vertex example, replacing just one
child by its converse changes C3 wr C3 (order 81, exponent nine) to C3³
(order 27, exponent three), despite leaving all labelled cyclic triples intact.

This gives the group, its order, its vertex action and its cycle index at
clearly separated levels. It does not reconstruct individual groups from an
aggregate index alone.

## 2. Composite families from the aggregate prime bank

[Universal formulas](COMPOSITE_FAMILIES.md) · [Exact tables](TABLES_THROUGH15.md)

Let P_k be the recovered degree-k aggregate prime cycle index, A the nonempty
tournament index, and S the strong index including the singleton. Mark each
prime k-node by u_k and each linear k-node by v_k:

    S = p_1 + sum_(k>=3) u_k (P_k ∘_w A),
    A = S + sum_(k>=2) v_k S^k.

The weighted Adams operation multiplies child statistics as well as vertex
sizes: ψ_d(u_k)=u_k^d and ψ_d(v_k)=v_k^d. A repeated child on a d-cycle occurs
d times, not once. Ignoring that rule gives fractional, incorrect node-count
distributions while still passing the unmarked count check.

The unique weighted-degree solution generates full cycle indices and the joint
prime/linear-node profile. Prime objects, composite strong objects and linear
root objects have indices P, P∘A−P and A−S respectively. At order eleven this
splits the 903,753,248 ordinary classes into:

| Root family | Classes |
|---|---:|
| Modular primes | 796,876,429 |
| Strong composites | 87,588,161 |
| Nonstrong ordered sums | 19,288,658 |

These are exact formula outputs, not a new exhaustive order-eleven atlas.
Full marked arithmetic was checked through fifteen; independent actual
canonical-module extraction checks every marked coefficient through six.

## 3. Where cycle indices stop, and the exact replacement

[Information theorem and explicit witnesses](SYMMETRY_INFORMATION.md) ·
[Standalone verifier](check_symmetry_collision.py)

Two explicit prime tournaments on 81 vertices have full regular automorphism
groups C3⁴ and Heisenberg27×C3. Both have exactly

    Z = (p_1^81 + 80 p_3^27)/81.

Inflate each slot by either TT1 or TT2. Both resulting families have the same
ordinary generating function

    ((x+x²)^81 + 80(x³+x⁶)^27)/81.

Nevertheless, at each of orders 108 and 135, they have **40 versus 13** classes
whose automorphism group has order 27. Explicit matrices, complete translation
checks, discrete rooted refinement, module closures and subgroup-coset checks
certify the example. No regular-representation existence theorem is assumed.

To compute this finer distribution, retain the slot action and subgroup marks.
For each actual subgroup K of a quotient group G, fixed child-color assignments
have polynomial

    F_K(t) = product_(O a K-orbit) sum_(child type c) t_c^|O|.

Subtract larger stabilizers to obtain E_K, the polynomial for exact stabilizer
K. The number of decoration orbits of its conjugacy type is

    |K| E_K / |N_G(K)|.

Child-type multiplicities then determine the product of internal group orders;
multiply by |K|. Actual group types require the semidirect-product action too.
This is a precise sufficient replacement, not a claim that every possible
symmetry question needs the entire molecular inventory.

## 4. Converse symmetry and ternary intransitivity

[Converse recovery](CONVERSE_RECOVERY.md) · [Cyclic-triple families](C3_FAMILIES.md)

The ordinary index alone does not contain anti-automorphisms. The second
Burnside index B_τ counts permutations taking a tournament to its converse.
Such permutations have only cycles of length 2 modulo 4, plus at most one
singleton. Their fixed orientations are explicitly counted.

With ordinary and twisted strong indices S_e,S_τ,

    S_e = 1−1/B_e,
    S_τ = B_τ(1−ψ_2 S_e)−1.

The second formula reverses the ordered SCC list. Prime substitution then uses
the converse child action on odd quotient cycles and the ordinary action on
even ones. Weighted-degree inversion recovers the complete aggregate prime
anti-index P_τ without a supplied prime census.

A specified full anti-cycle gives a small scalar boundary: for odd m, with
f_m=2^((m−1)/2), the prime anti-fixed counts satisfy

    h*f = 2^m−2f_m.

At cycle lengths 2,6,10,14,18 these counts are 0,4,24,112,472. They are fixed
points of a prescribed permutation, not self-converse isotype counts.

For realizable cyclic-triple hypergraphs, the cited prime-uniqueness theorem
identifies the prime root index as R=(P_e+P_τ)/2. The order of the SCC blocks is
forgotten, giving

    D = X + R∘(C−1),             C = SET(D).

This now generates the full aggregate C3-hypergraph index without an unknown
prime bank. It does not count arbitrary 3-uniform hypergraphs, and it does not
restore the false assertion that ordinary tournament automorphism order is a
C3 invariant. Independent full-index arithmetic reaches fifteen; actual
hypergraph-image checks reach six.

## 5. Two genuine bounded theories inside the unbounded problem

### Fixed excess: polynomial counts and finitely many moving actions

Let h be the maximum number of prime nodes on a canonical root-to-leaf path,
and write n=2h+1+r. At each fixed r, the number of exact-depth-h tournaments is
eventually a polynomial in h of degree r with leading coefficient 4^r/r!.
For positive h, the first three boundaries are

    r=0: 1,
    r=1: 4h−1,
    r=2: 8h²+11h−10 for h>=2 (the h=1 value is 10).

The stabilized transfer is multiplication by a fixed series L on
Q[x]/(x^(r+1)). It has one unipotent Jordan block, so the depth sequence obeys
an exact (r+1)-st forward-difference recurrence after its stated threshold.
This is a closed, bounded family calculation, not a bounded unrestricted atlas.

The same transfer has a rooted-count description:

    x²L = 1−B(x)²/D(x),
    D(x) = sum_(N>=1) p_N x^(N−1),

where p_N counts pointed tournament classes from the previous re-rooting work.
It is not N times the ordinary count. Thus the attachment symmetry count and
long-depth composition kernel meet in one exact identity.

There is also a structural bound, not merely an aggregate count. For h>=1,
[the support theorem](BOUNDED_SYMMETRY.md) proves

    number of vertices moved by Aut(T) <= floor(3r/2)+3.

It is sharp for every even r, at the specified sufficiently large depths.
Therefore only finitely many moved permutation actions can occur at a fixed
excess, although arbitrarily many globally fixed vertices may be present.
No optimality at odd r or complete equality classification is claimed.

### Maximum prime-node density: a complete symmetry-order recurrence

Let b count all prime nodes, not maximum depth. At n=2b+1, every internal node
is C3 and there are no linear nodes. If u marks the exponent e in |Aut|=3^e,
then the entire family satisfies

    C(z,u) = 1 + z[(C(z,u)³−C(z³,u³))/3 + u C(z³,u³)].

Identical child triples gain a quotient rotation; nonconstant triples do not.
This restricted family admits a complete symmetry-order recurrence despite
the general cycle-index obstruction above. Its unique chain has Aut=C3; the
upper bound 3^b is attained exactly for complete balanced ternary trees.

Node count and depth must not be confused. At order six, the maximal two-node
population has eight types, but only seven have prime depth two; the ordered
sum C3→C3 is the eighth, with depth one.

## Reproduction and scope

All public scripts use Python's standard library and nonremovable checks:

    python3 check_groups.py
    python3 check_symmetry_collision.py
    python3 check_marked_composites.py 9
    python3 check_depth_diagonals.py
    python3 check_density_symmetry.py
    python3 recover_converse.py 15
    python3 recover_c3.py 15
    python3 check_equivariant.py
    python3 check_support.py

Repeat with `python3 -O`. The marked order-fifteen run is optional and takes
roughly eight minutes per mode on the development host; both modes were run
and their outputs agreed. The small default deliberately supports quick replay.
The larger constructed support examples perform direct module tests and can
also take longer than the small-order checks. No saved receipt or external
class atlas is an input to the public scripts.

The coherent conclusion is not that one numerical invariant determines every
tournament symmetry. It is that canonical composition has exact compatible
calculi at the group-action, subgroup-mark, equivariant-cycle-index and
counting levels, with explicit maps and counterexamples marking each loss of
information. The unrestricted small structural carrier and full prime-class
atlas remain open in this work.
