# Uniform-incidence extraction: audited limits and finite-menu size

**PROVED, independently audited in the stated scopes.** No eventual odd Hamiltonian-path spectrum threshold is proved. These results rule out a specific way to convert a prescribed binary index into an affine path count; they do not rule out individual tournament counts or the separate two-index occupancy conjectures.

## Results

- [Single-core obstruction](uniform/PROOF.md): any fixed uniform binary-incidence table fails to realize either 6u-1 or 6u+1 throughout 8193<=u<=16384, irrespective of core order.
- [Finite-menu obstruction](structural/FINITE_MENU.md): even a target-dependent selection from any bounded number q of such cores cannot cover every sufficiently large index. The cores themselves may change with digit length. For 6u±1 and m>=28, a covering menu must satisfy `(q+1)(2q-1)>m-8`. Thus q grows at least on the order of sqrt(m), or sqrt(log u). Any menu of at most three cores fails on the block 2^28+1<=u<=2^29.
- [Structural specialization](structural/PROOF.md): a fixed letter-local graft retaining the marked-edge E(u) graph cannot produce any affine count au+b on a tail. Explicit transitive-module pumps force an ordered deletion graph, after which the coefficient of the inherited nonlinear source function g is exactly the positive integer H(B)H(S).

The common physical mechanism is important: two sensitive core vertices allow three independent path slots on a common pinned subcube, while one sensitive vertex has strictly negative pair capacities. A colored-rectangle argument handles arbitrary target-dependent selection from a finite menu. The specialized pumping proof additionally identifies the residual source coordinate explicitly.

## Verification

All code uses the Python standard library and explicit exception gates, so checks remain enabled under optimization. Run each script from its own subdirectory:

    cd uniform
    python3 verify.py
    python3 -O verify.py
    python3 anchor_controls.py
    python3 -O anchor_controls.py

    cd ../structural
    python3 verify.py
    python3 -O verify.py
    python3 verify_menu.py
    python3 -O verify_menu.py

The independently written audit scripts in `audit/` import no author code. Normal and optimized receipts agree. The all-order claims follow from the proofs, not the bounded controls.

## Source and audit map

The original scratch provenance labels `session203-complement-gadget/` and `session203-structural-induction/` in frozen proofs correspond respectively to `uniform/` and `structural/` here. Their mathematical bytes are preserved after independently approved status promotion. Pre-promotion versions are retained under `originals/`.

[Independent audit](audit/CONSTRUCTOR_AUDIT.md) records exact reviewed input hashes and all scope checks. The promotion receipt records approved status-only changes. Frozen inherited statements in `dependencies/` provide the edge-sector construction, ordered root-amalgamation law and original weighted-cut formula; their mechanisms are attributed, not claimed as new.

## What remains open

This packet leaves positional incidence rules, changes to digit-to-digit arcs, nonlinear encodings, an unbounded target-dependent menu, two-source arithmetic, and unrestricted graph rebuilding outside its obstruction. In particular the twelve-form and forty-form source-pair occupancy conjectures remain open and are not contradicted. The target of proving every sufficiently large odd integer realizable, then lowering the threshold toward 23, is unchanged.
