# Session213: endpoint equality and limits of cycle-to-path conversion

**Status: PROVED structural classification, scoped identities and construction barriers; FINITE-EXACT controls. The unrestricted Chan–Pak Hamiltonian-path spectrum remains OPEN. No literature-priority claim.**

Write H(T) for directed Hamiltonian vertex words and C(T) for directed spanning cycles modulo cyclic rotation. [Session212](../session212/README.md) realizes every nonnegative cycle count. Session213 proves an all-tournament endpoint equality classification, then identifies why that binary cycle compiler cannot simply deliver every odd total-path count.

## Main results

1. **All equality cases of the endpoint half-bound.** Let e_a(T)>0 count paths ending at a. The [equality theorem](cycle-arithmetic/EQUALITY.md) proves H(T)=2e_a(T)+1 exactly when T is a transitive prefix followed by R_m, m>=3, with a penultimate in R_m. Here R_m is the transitive m-vertex order with only its extremal arc reversed. In these cases e_a=2^(m-3) and H=2^(m-2)+1. Reversing all arcs gives the starting-endpoint version. The zero-endpoint boundary is separate: transitive tournaments can have H=1=2·0+1. The endpoint half-bound was inherited; its previously finite-only equality observation now has an all-order proof.

   Consequently, deleting the sole-exit vertex from **any** positive-cycle tournament can return H=2C+1 only when C is a power of two; every power is possible. This is a universal obstruction to that endpoint-preserving conversion, not to unrelated tournaments realizing 2C+1. A power-of-two endpoint count alone does not force equality.

2. **The full binary cycle compiler has sparse total-path support.** The [explicit path-allocation proof](context-transform/PROOF.md) bounds the number of its distinct values H<=X by O(X^alpha), alpha=log_3(2)<1. Both apex orientations and all leading-zero padding are included. The [odd-cycle-collection proof](path-gadgets/PROOF.md) independently gives the same upper-bound exponent. Neither claims a sharp exponent, an asymptotic equivalence, or a result about every two-apex tournament.

3. **Finite extension menus, SCC products and finite seed banks do not repair that sparsity.** The [finite-menu corollary](context-transform/PROOF.md) permits any fixed finite collection of deterministic recipes on marked comb presentations, with any finite number of added vertices, provided each output retains its input as an induced subtournament. The [multiplicative-closure theorem](context-transform/MULTIPLICATIVE_CLOSURE.md) then allows arbitrary finite dominance products of those outputs and any fixed finite seed bank. For every beta>alpha, the resulting distinct values at most X are O_beta(X^beta), hence density zero by choosing beta<1. The endpoint exponent alpha is not asserted for products.

   These statements do not cover deletion, reorienting old arcs, unboundedly many outputs per input, new cross-block interactions, an infinite independent seed bank, or a new operation applied to a product. They do not rule out cycle-to-path constructions in general.

## Exact tools and narrower consequences

- **Quadratic response and complete family inverse.** The [closed comb formula](cycle-to-path/PROOF.md) evaluates every needed ordered-two-path-cover coefficient. General attachment quadraticity and submodularity are inherited from Session203. The bound H>=2^m+1 gives an exact all-length inverse for this family. It attains 5,072 values through 65,535; H=11 and 35 are absent at every order within the family. The first missing odd target at least 23 is 35. These are not ambient exclusions. Fixed positive C fibers also have unbounded H, blocking scalar-only bounds for extensions retaining every padded input.
- **Exact deletion error.** The [one-apex arithmetic](cycle-arithmetic/PROOF.md) gives H=2h+D(h-1), with D the residual one-apex path count. The [four-branch addendum](cycle-arithmetic/ADDENDUM.md) classifies the positive input cycle values for which its specified deletion gives 2C+1: precisely {2^k:k>=0} union {3·2^j+1:j>=0}. The extra forward-apex cases do not contradict the universal endpoint theorem: there C differs from the retained starting-endpoint count.
- **A unit cycle increment has an exponential path cost.** The [apex reversal identity](context-transform/ARC_FLIP.md) has sharp lower bound H(B+)−H(B−)>=(2−tau)2^m, although C(B+)−C(B−)=1. This is family-specific, not a general arc-reversal monotonicity law. The [scalar-growth collision](path-gadgets/PROOF.md) separately exhibits equal parent H but unequal H after the same bit insertion.

## Audit and reproduction

The [independent audit](audit/AUDIT.md) separates proof review, author replays and independent exact checks. Its [reproduction guide](audit/README.md) gives commands; [approved file hashes](audit/APPROVED_PACKETS.json) pin the reviewed bytes.

- Endpoint equality: all 33,867 labeled tournaments through order six, plus the inherited 456 and 6,880 isomorphism classes at orders seven/eight and 1,981 explicit paths outside the endpoint injection image. The higher-order class universe is inherited, not newly enumerated.
- Comb formula: 2,040 author cases; independent cover counts and 262,136 naive quadratic evaluations reproduce the complete 5,072-value support without the author's Gray traversal.
- [Arithmetic](cycle-arithmetic/README.md), [quadratic response](cycle-to-path/README.md), [path allocations](context-transform/README.md) and [odd-cycle certificates](path-gadgets/README.md) supply standalone Python verifiers, pinned sources and deterministic receipts. Equality additionally uses `equality_verify.py` and `EQUALITY_RESULT.json`.

Run each verifier normally and with `python3 -O`, directing output outside the release files and comparing it with its receipt. All-order claims follow from the written proofs, not finite extrapolation.

## Remaining target

Construct a tournament with every odd total-path count at least 23. This session does not prove eventual completeness or find a new ambient forbidden value. The structural endpoint classification closes one tempting conversion route, while the family formulas and sparsity theorem identify what a successful replacement must overcome.
