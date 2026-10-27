# What prime cycle indices forget, and an exact sufficient replacement

**Status:** PROVED elementary witnesses/formulas; FINITE-EXACT standard-library
checks. Scope distinctions below are essential. No claimed new prime census,
no claim that the constructed inventories are the full prime family.

## 1. Two distinct losses

For a permutation group `G≤S_m`, put
`Z_G=(1/|G|) Σ_g ∏_j p_j^{c_j(g)}`. A prime-family cycle index is
`P_m=Σ_[Q prime, |Q|=m] Z_Aut(Q)`.

* **Aggregation loss:** even the orders of individual automorphism groups
  cannot generally be recovered from a sum of their cycle indices.
* **Action-index loss:** even an individual `Z_G` need not determine colored
  stabilizer distributions. It remembers cyclic-element fixed data, not
  simultaneous fixed data for all subgroups.

These are not a denial of ordinary Pólya substitution: aggregate cycle
indices are sufficient for the corresponding ordinary composition counts.
For one individual faithful action, the coefficient of `p_1^m` is `1/|G|`;
thus its group order *is* recoverable. That fact does not survive aggregation.

## 2. Actual tournament inventory counterexample, all order eight

Notation: `A→B` is the ordered sum; `1` is a singleton; `C=C3`;
`D=C3[TT2,TT2,TT2]`. Consider these two four-element sets of tournament
isomorphism classes:

```
L = { C→C→1→1, C→1→C→1, 1→C→C→1, TT8 },
R = { C→TT5, 1→C→TT4, D→TT2, 1→D→1 }.
```

Their ordered strong-component size words distinguish all eight classes.
The first three members of L have group `C3×C3`; its last member is rigid.
Every member of R has group `C3`. Both aggregate cycle indices are exactly

```
(4/3)(p1^8 + p1^5 p3 + p1^2 p3^2).
```

Consequently their class counts, all ordinary color-orbit counts, and all
ordinary substitution counts agree. Their automorphism-order inventories
are `{1:1,9:3}` and `{3:4}`; their first group-order sums are `28` and `12`,
and second sums `244` and `36`. Even adding the class count does not repair
this loss. The basic polynomial relation is

```
3 Z_(C→C) + Z_TT6 = 2 Z_(C→TT3) + 2 Z_D.
```

Padding supplies distinct actual classes instead of artificial repeated
copies. **Boundary:** these are decomposable tournaments, not prime
quotients. This disproves a general inverse claim for cycle-index inventories;
it does not establish a collision between two genuine complete prime banks.

## 3. Explicit prime witnesses at order 81

**FINITE-EXACT:** `check_symmetry_collision.py` supplies two explicit prime
tournaments with full regular automorphism groups `E81=C3^4` and
`J=Heisenberg27×C3`. They have the identical *individual* action cycle index

```
(p1^81 + 80 p3^27)/81,
```

but nonisomorphic automorphism groups and different decorated-stabilizer
histograms. This establishes the information-loss claim for actual prime quotients
(at order 81, not within a census through order 15).

Enumerate `(a,b,c,d)∈F3^4` lexicographically. For E81 use coordinate addition;
for J use

```
(a,b,c,d)(x,y,z,w)=(a+x,b+y,c+z+a y,d+w) mod 3.
```

Orient `i→j` when `i^-1 j` is in the following index set:

```
E81: 1,6,7,8,9,10,11,15,16,17,24,25,26,27,36,40,41,43,44,48,
     49,51,52,55,56,57,58,59,60,61,62,63,64,65,67,70,73,74,75,78
J:   1,5,6,8,9,10,13,14,15,16,17,19,24,28,29,30,31,32,34,35,
     42,43,46,48,49,51,53,54,57,65,67,69,70,72,73,74,75,76,77,79
```

For each, all 81 translations preserve arcs. Root-pinned equitable color
refinement (counting outneighbors in each current color class) becomes discrete
in three rounds: the numbers of cells are 2,3,15,81 and 2,3,25,81 respectively.
Every zero-fixing automorphism preserves each refinement color, so it must
fix every vertex. No backtracking is needed for this rigidity certificate. Every one
of the 3,240 pair-module closures is full, independently proving primality.
The standalone verifier reconstructs every adjacency row from the displayed constants.

E81 has 40 subgroups of index three; J has 13. Indeed every index-three
subgroup of a finite 3-group is normal and is the kernel of a nonzero
homomorphism to C3. The abelianizations have ranks four and three, so the
counts are `(3^4−1)/2=40` and `(3^3−1)/2=13`. For J, the commutator subgroup
is precisely the c-coordinate. The verifier constructs all these kernels
and checks the homomorphism identities and all their coset stabilizers.

Each such subgroup fixes six nonconstant binary colorings, all with exact
stabilizer order 27; each coloring orbit has size three. Thus the actual
prime actions have **80 versus 26 binary-decoration orbits with stabilizer
order 27**, despite identical individual cycle indices. Their total binary
color-orbit counts agree, by Burnside.

For the sharp size-graded tournament consequence, take the two colors to be
`TT1` and `TT2`, both rigid. A coloring with exact stabilizer 27 doubles
either one or two cosets of an index-three subgroup: 27 or 54 slots. For
each subgroup the three single-coset choices form one orbit, as do their
complements. Substitution therefore gives:

| total tournament order | E81-quotient classes with Aut-order 27 | J-quotient classes with Aut-order 27 |
|---:|---:|---:|
| 108 = 81+27 | 40 | 13 |
| 135 = 81+54 | 40 | 13 |

At *every* size the two families have equal ordinary counts, since their
cycle indices agree and both use `A(x)=x+x^2`. Yet their automorphism-order
histograms differ at these exact sizes. Prime quotient uniqueness makes the
substituted blocks intrinsic, so no unintended identifications occur within
either family, and the internal groups are trivial.

An equal-size variant uses rigid colors `TT4` and `C3[TT2,1,1]`, giving
80 versus 26 classes of automorphism order 27 at the single size 324.

## 4. Exact enriched inventory and transfer

For prime substitution `T=Q[T_i:i∈V(Q)]`, retain, for every prime isomorphism
class, its actual permutation group `G_Q≤Sym(V(Q))`, not merely `Z_G`.
For children retain their isomorphism types and automorphism orders (and
actual groups/actions if group type rather than just order is requested).
Choose one representative of each child type. A decoration is a map
`f:V(Q)→C` into this palette. Then

```
G_f = {g∈G_Q : f(g i)=f(i) for every i},
|Aut(T)| = |G_f| ∏_i |Aut(T_i)|.
```

For a canonical prime substitution the blocks are intrinsic. The exact
group is `(∏_i Aut(T_i)) ⋊ G_f`, where the action permutes equal-type factors
using the chosen representatives. Sum once per G-orbit of decorations.
At ordered-sum nodes the group is the direct product; do not introduce
permutations of equal ordered components.

### A table-of-marks implementation

Fix G and its slot action. For each *actual subgroup* K define the weighted
fixed-decoration polynomial

```
F_K(t) = ∏_[O∈V/K] (Σ_[c∈C] t_c^|O|).
```

For a finite palette this is literal; infinite palettes are interpreted
coefficientwise with a positive grading. If `E_K(t)` counts decorations
whose stabilizer is exactly K, then

```
F_K = Σ_[J≥K] E_J,
E_K = F_K - Σ_[J>K] E_J.
```

For one representative K of each subgroup conjugacy class the polynomial
counting decoration orbits of that stabilizer conjugacy type is

```
(|K|/|N_G(K)|) E_K(t).
```

Monomials retain the child-type multiplicities. Thus they determine total
size and the product of child automorphism orders; multiply that product
by `|K|` to obtain the resulting automorphism order, or raise to any desired
moment before summing. This is an exact finite recipe in each size degree.

Equivalently, retain the table of marks `M_[K],[J]=|(G/J)^K|` and the vector
`F_K`. Solving `F_K=Σ_[J] M_[K],[J] a_[J]` yields the weighted number of
orbits of stabilizer type J directly. **The bare abstract group's table of
marks is not by itself the required input:** one also needs the slot action
(or enough subgroup orbit-partition data to construct every F_K), plus
child-type weights. A concrete permutation group provides all of these.

For only automorphism-order distributions, full group structure may be
more information than necessary, but this is a precise sufficient inventory.
For actual automorphism-group isomorphism types, keep the semidirect-product
action as well. For correlations under augmentation, retain how decorations
and subgroups change under the specified operation; present marks alone do
not encode a future operation that has been forgotten.

The traditional molecular species representation retains a transitive
`S_m`-set `S_m/G` for each unlabelled prime Q, hence the subgroup embedding
up to conjugacy, with multiplicity. This is another sufficient structural
inventory. Replacing every molecule by its cycle index is the lossy map
exhibited here.

## 5. Inherited prospective-symmetry boundary

The inherited prospective-symmetry obstruction proves that
`C3[TT_(3k-3),1,1]` and `C3[TT_(k-1),TT_k,TT_k]`, for every `k≥2`,
are rigid tournaments of order `3k-1`, with the same actual pair-free core,
SCC-size data and individual cycle index `p1^(3k-1)`. Their one-vertex
extension masses to core order three differ by exactly `2/3`.
This is a different loss: present action data alone does not encode
prospective stabilizers after the weight vector changes. It does not
contradict prime substitution from the complete decorated tree.

## Verification and source scope

Run from this directory:

    python3 check_symmetry_collision.py
    python3 -O check_symmetry_collision.py

This standalone, independently written verifier reconstructs the nonabelian
operation via upper-triangular matrices over F3 and the abelian operation by
vector addition. It checks associativity on all triples, commutators, the 81
translations, discrete refinement, every one of the 3,240 pair-module closures,
all action cycle types, the 40 and 13 subgroup kernels, exact coset-subset
stabilizers, and the stated inflation counts. It also checks every permutation
of each order-eight witness and verifies their distinct SCC words. All guards
remain active under optimization; normal and optimized outputs are identical.
The all-order subgroup-inversion formulas are proved above, not inferred from
these finite witnesses. No external package or input data file is required.

The group-theoretic recipe uses the [canonical symmetry law](GROUP_LAW.md).
The proof of prime-cycle-index sufficiency for ordinary counts remains in the
[previous recovery](../session79_tournament_enumeration_20261027/PRIME_CYCLE_INDEX.md).
Our counterexample does not refute it: a distribution of stabilizer orders is
strictly more information than an ordinary orbit count.

### Why the tempting smaller abelian witness was rejected

The regular actions of C3^3 and the order-27 Heisenberg group have the same
cycle index, but C3^3 cannot be the full regular automorphism group of a
tournament on 27 vertices. There are only 2^13=8192 connection sets, fewer
than the 11232 elements of GL(3,3); a trivial stabilizer would require a
larger orbit than that entire set. Thus a nonidentity linear automorphism
always survives. This classical exclusion is recorded, with credit to
Godsil's 1986 correction, in [Wong et al. (2025), introduction](https://arxiv.org/html/2501.13406v1).
The order-81 witnesses above were explicitly constructed and independently
verified rather than justified by an overbroad regular-representation claim.
