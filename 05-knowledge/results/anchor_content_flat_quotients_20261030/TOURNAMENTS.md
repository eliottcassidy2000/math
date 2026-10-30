# Forget a marked base by taking the union before counting

**PROVED:** affine-flat closure and exact orbit-count formula below.
**FINITE-EXACT:** the five-vertex native experiment. A compact unrestricted
A000568 continuation grammar remains **OPEN**; this is not a new scalar
recurrence for a sequence already computable by Burnside methods.

## Inheritance first

The strongest nearby mechanisms are already stronger than a fresh generic
relabeling search:

- [Hamilton relocation](../hamilton_path_relocation_connectivity_20261020/THEOREM.txt)
  connects every ordinary class's fixed-path presentations.
- [Frobenius/LCA, E1–E3](../frobenius_lca_metagraph_20261028/THEORY.md) has the
  exact reset projector and all-order reconstruction; its trace is A000568.
- [Root switch](../session79_tournament_enumeration_20261027/ROOT_SWITCH.md)
  gives `Q^2=NQ` and `tr Q=N a_N` using pointed parents and aligned deletions.
- [Exact generation](../session87_four_vertices_enumeration_20261028/EXACT_GENERATION.md)
  already uses canonical roots and exact transporter ties, through order8.
- [Native path corrections](../endpoint_shadows_path_corrections_20261030/TOURNAMENTS.md)
  preserve endpoints and quantify complement-program spectra, but retain the
  original marked base. The present interface removes that base from the
  *family membership* predicate, not from a weighted history.

The canonical hostile is quotienting a family by the normalizer of its first
bank after legal moves have left that bank. The least-used sidecar is the
intersection of differently marked copies. The 81/27 transporter controls in
the exact-generation work remain necessary; degree and spectral summaries
are not substitute transporters.

## 1. Exact closed family carrier

For fixed labelled vertices, choose ascending edge coordinates in
`V=F_2^(choose(n,2))`. A bit1 means the smaller vertex points to the larger.
There are no ties. Retain any named ports as separately fixed labels.

Represent a tournament family by a finite union of affine flats `A x=b` over
F2, with redundant descriptions merged as sets. The following operations
are exact on this carrier:

1. **Relabeling:** an affine edge permutation, with a sign bit whenever the
   order of an edge's endpoints reverses. Literal converse is translation
   by the all-ones vector.
2. **Native complement along a legal Hamilton word:** intersect with the
   n-1 path-edge equations, then translate by the indicator of all nonpath
   edges. Taking the union over allowed words preserves the exact legal
   image, even when several words lead to the same tournament.
3. **Block substitution:** for fixed block labels, take the product of
   quotient and internal flats and copy each quotient bit to its entire
   interblock rectangle. This is an affine image. Unions over partitions
   and alternative descriptions retain shared ports and merge overlaps.
4. **Linear attachment predicates:** add their equations before taking the
   image. Arbitrary finite predicates can be unions of singleton flats,
   so closure alone makes no compactness claim.

If the original base is only a construction witness, forget it by taking
the union of the resulting sets. Future operations above inspect the current
tournament and retained ports, so membership in this union is sufficient.
This operation destroys which base, Hamilton word, and history produced a
point. It is invalid for observables that still request those data or their
multiplicities. In particular it does not license forgetting a Collatz source.

## 2. Counting without canonicalizing each tournament

For flats F_i, build the exact indicator

    1_U = sum_F c_F 1_F,

combining repeated intersections in inclusion–exclusion. Intersections are
computed by binary row reduction. For an invariant union U and a vertex
permutation pi, its fixed set is another affine flat, possibly empty. Hence

    |U/S_n| = (1/n!) sum_pi sum_F c_F |F intersect Fix(pi)|,             (1)

where each nonempty intersection has size `2^(choose(n,2)-rank)`.
Cycle-type representatives are valid because the *union* is invariant;
individual flats need not be. Saturate under S_n before using that reduction.
With named ports use their pointwise stabilizer instead, or its orbit-saturated
family. No graph canonicalization is required in this producer formula.

A specified nontrivial vertex set S is a module precisely when every outside
vertex sees all members of S the same way. These are affine bit equations.
The complement of the union of all such module flats is the prime-tournament
indicator. Multiplying that indicator into (1) gives an exact prime filter,
including nonprime boundaries rather than discarding them by a score proxy.

Two parallel distinct affine flats of equal dimension have an affine-flat
union one dimension larger over F2. This gives an exact coalescing rule.
Also, if the exact union cardinality is the whole cube's cardinality, replace
it by the cube. These rules can drastically simplify a family, but greedy
pairing is not a minimum-cover algorithm. Intersection arrangements and
permutation saturation may still grow exponentially or factorially.

## 3. Complete five-vertex test

Start with the four labelled orientations fixed by the cycle `(0 1 2 3 4)`;
take all relabelings. This gives the 24 labelled copies of the regular
five-vertex class. At each depth take **all** legal native complements.
The depth is exact, not a union of all earlier depths.

| depth | labelled tournaments | ordinary classes | prime classes | flats after coalescing |
|---:|---:|---:|---:|---:|
| 0 | 24 | 1 | 1 | 6 |
| 1 | 280 | 3 | 1 | 37 |
| 2 | 784 | 8 | 3 | 148 |
| 3 | 1024 | 12 | 3 | 1 |

Thus every five-vertex class is reached after three native moves from this
unmarked starting class. A single full-cube flat certifies the final family.
The intermediate inclusion–exclusion representation at depth2 has374 terms,
which is a visible warning that set covers and intersection complexity differ.

The original bank normalizer has order20. On the full five-vertex cube it
counts52 orbits, whereas the correct full symmetric quotient has12. Retaining
one named vertex gives48 classes. These three answers measure three different
objects, not alternate implementations of the same quotient.

[flat_compiler.py](flat_compiler.py) computes the rank formulas. Independent
literal bitmask traversal audits every native image and every union indicator
on all1024 tournaments. The prime oracle reads outside-to-inside arrow values
directly, separately from the affine module equations. Full-group bitmask
orbits audit the quotient. Controls at orders3 and4 reproduce2 and4 classes
and prime counts1 and0. The tests use no disabled Python assertions.

## 4. Interface to the incoming ramified-bank recurrence

The incoming [shared-quotient bank recurrence](../session111_scale_connections_20261030/BANK_RECURRENCE.md)
counts a specified odd subgroup's orientation-bank orbits at all depths.
Its common quotient `C_n x_(C_(n/5)) C_n` is essential: the two top coordinate
actions cannot be made independent at the next level. That is a strong compact
action-level component; it is not yet the full graph quotient.

Use it before base-forgetting. At the forgetting boundary, the present flat
intersection formula supplies the exact overlap correction and module filter.
The unresolved target is a bounded-size description of the resulting
intersection arrangement, or an exact transporter theorem that avoids that
arrangement. Neither is supplied by the experiment at n=5.
