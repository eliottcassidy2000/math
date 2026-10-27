# Squarefree symmetry and the sharp fifteen vertex boundary

**PROVED** using the **CITED** common strong-module tree and C3 realization
theorem stated precisely in [the inherited fibre result](../tournament_composition_classification_20261026/fibre/RESULT.md).
That source supplies the prime threshold4p and the recursion

```text
|Aut(T)| = |Gamma(Q; [T_i])| product_i |Aut(T_i)|,           (1)
```

where Gamma preserves the ordinary oriented isomorphism colours of the
canonical children. At a linear node the quotient group is trivial.
Equal labelled cyclic-triple hypergraphs have the same canonical child
sets; prime quotients agree up to converse and have the same permutation
automorphism group. These are the structural hypotheses used below.

## S1 The squarefree threshold is 5p

**Theorem.** Let T,T' have the same labelled cyclic triples on n vertices.
Suppose at least one of their automorphism orders is squarefree. For every
odd prime p,

```text
n < 5p  implies  v_p(|Aut(T)|)=v_p(|Aut(T')|).             (2)
```

This is sharp for every odd prime p, with both orders squarefree in the
sharp examples. It strengthens the unconditional4p threshold only under
the stated squarefree hypothesis.

**Proof.** Name the squarefree realization T and induct on n<5p. Every
canonical child T_i also has squarefree automorphism order, because its
order divides |Aut(T)| in (1). Induction gives equality of the p-valuations
of corresponding child orders. Linear nodes are therefore settled.

At a prime node suppose the p-parts of the two Gamma groups differ. A
Sylow p-subgroup P of the group with larger p-part cannot be contained in
the other group. It consequently has an orbit O whose child colours all agree in the
larger group but fail to agree in the other realization. P preserves child
sizes. Such an orbit is nontrivial, hence has at least p ports. Fibres of
size at most three have only one ordinary isomorphism type, so each child
in O has size at least four. Since n<5p, that common size must be four.

The only ambiguous C3 type on four vertices has realizations
`C3 -> 1` and `1 -> C3`, both of automorphism order3. Because P preserves
the child type in its realization, all children in O have this C3 type;
the same holds in T. Thus the child product in (1) for T is divisible by
`3^|O|`, contradicting its squarefree automorphism order. The quotient
p-parts agree, and (1) completes the induction. The p=2 valuations are
always zero, since tournament automorphism groups have odd order. QED.

The tiny four-vertex classification used here is checked on all64 labelled
orientations and all24 permutations. The two other C3 types each have one
ordinary realization type. This is a finite lemma, not a numerical
extrapolation in the number of vertices.

**Sharpness.** Let S4 be the unique strong four-tournament and put
`D=S4 -> 1`, `E=D*`. They are rigid, nonisomorphic five-tournaments with
identical labelled cyclic triples (choose the same four-set and isolated
vertex). Let Q_p be the carousel on Z/p, with i beating i+j for
`1<=j<=(p-1)/2`. It is prime with automorphism group C_p, as proved in the
[carousel packet](../cyclictriple_composition_20261026/README.md), carousel lemma.
Its vertex stabilizer is trivial. Formula (1) gives

```text
U_p = Q_p[D,...,D],       |Aut(U_p)|=p
V_p = Q_p[E,D,...,D],     |Aut(V_p)|=1.                    (3)
```

Changing D to E preserves every labelled cyclic triple: triples within a
block are unchanged, those meeting exactly two blocks are transitive,
and those meeting three blocks depend only on Q_p. The common size is5p.
Both orders are squarefree and their p-parts differ. Direct adjacency
backtracking verifies the expanded tournaments at p=3,5,7,11,13, on
15,25,35,55,65 vertices. QED.

The mechanism is rigid ambiguity. An oriented child colour first changes
at size4, but its internal group has order3. Repeating that child in a
nontrivial quotient orbit forces a square divisor. Size5 first permits a
changing colour whose internal group is trivial.

## S2 The first Möbius failure is fifteen

**Corollary.** The observer `mu(|Aut(T)|)` is constant on every labelled
C3 fibre through n=14 and first fails at n=15. The observer
`lambda(|Aut(T)|)` first fails at n=12.

**Proof.** If both orders are nonsquarefree their mu values are zero.
Otherwise apply (2): n<=14 is below5p for every odd prime p, so all prime
valuations agree and the orders are equal. At15, (3) for p=3 has orders
3 and1, and mu values -1,+1. The inherited4p theorem gives full order
invariance through11; at12 its81/27 pair has lambda values +1,-1. QED.

**Uniqueness at15.** Up to hypergraph isomorphism there is exactly one
realizable C3 hypergraph at15 on which mu of the automorphism order varies.
Take a pair with differing mu values and name its squarefree-order member
T. Each proper child has size<=14, so the preceding argument makes
corresponding child orders equal. Any difference must occur in Gamma at
a prime root. The unconditional4p theorem excludes p>=5 from that change
at15. Thus the two Gamma groups have different3-parts.

Choose a Sylow3 subgroup of the larger Gamma as in S1. A colour-changing
orbit has at least three equal-size children, each of size>=4. Size4
forces `3^3` in the child product of T and is impossible. The remaining
possibility at total size15 is exactly three size5 children and quotient
C3. Their common automorphism order c contributes c³ to the child product,
so squarefreeness forces c=1.

The complete labelled order5 census shows that the only C3 type with
more than one ordinary realization and a rigid realization is D/E above;
both are rigid. There are30 labelled hypergraphs of this type: choose
the isolated vertex in5 ways and the shared pair of the two triples on
the other four vertices in6 ways. The checker examines all1024 labelled
five-tournaments and all120 permutations, rather than filtering to a
preselected score sequence. This finite classification proves uniqueness
of the hypergraph `C3(C3[D,D,D])`. QED.

## S3 The complete first fibre

Write E=D*. The fixed labelled C3 hypergraph has1024 tournament lifts:
two root orientations and eight labelled D/E lifts in each of three
blocks. The four ordinary classes and their symmetries are

| Child types around C3 | Aut order | mu | lambda | Orbit size under Aut(C3 hypergraph) |
|---|---:|---:|---:|---:|
| DDD | 3 | -1 | -1 | 128 |
| DDE | 1 | 1 | 1 | 384 |
| DEE | 1 | 1 | 1 | 384 |
| EEE | 3 | -1 | -1 | 128 |

The hypergraph has `5³+3*2=131` triples. Cross-block pairs have codegree5,
whereas within-block pairs have codegree at most2. Thus its three blocks
are intrinsically recoverable. Inside a block, the isolated vertex is
intrinsic, and the two cyclic triples on the four-set share one pair.
Permutations preserving that shared pair and its complementary pair give
the internal hypergraph group `S2 x S2`, of order4. Therefore

```text
Aut(C) = (S2 x S2)^3 semidirect S3,       |Aut(C)|=384
sum_[T] 1/|Aut(T)| = 1/3+1+1+1/3 = 8/3 = 1024/384.
```

Canonical blocks and oriented child types separate the four rows; cyclic
rotation identifies all placements with the same number of E children.
All1024 labelled lifts and all384 hypergraph permutations are explicitly
generated; the four disjoint orbits exhaust the lift set.

The common Hamilton-path count is **152377413**, a **FINITE-EXACT** count
independently checked for all four representatives by forward subset DP
and backward memoized recurrence. Equality itself follows from the
inherited all-order path-cover invariance in a C3 fibre. Automorphism
orders are independently verified by adjacency backtracking and by
exhaustive permutations within colour-refined cells. The numerical count
is not needed for the symmetry threshold proof.

## Scope of the composition connection

Source: a labelled tournament with its canonical decorated module tree.
Target: the cyclic-triple hypergraph. Map: forget the permissible child
orders and prime orientation choices. Preserved: strong-module sets,
path-cover counts, and the p-parts below their applicable4p or5p boundary.
Lost: the oriented child-colour stabilizer. Necessary sidecar: the action
on those colours, with internal group orders. Cheapest decisive tests:
the81/27 twelve-vertex pair and the3/1 fifteen-vertex pair.

These are classification thresholds for an observer under module
substitution. A Hamilton-path complement line can destroy the old proper
module sets; it is not one of the moves that preserve this C3 fibre.
