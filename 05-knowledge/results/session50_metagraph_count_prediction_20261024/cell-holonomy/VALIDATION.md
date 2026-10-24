# Validation scope and reproduction

All proof assertions in `check.py` use explicit exceptions, not Python
assert statements, and remain enabled with optimization. The script imports
no Session50 growth code or canonicalizer. Its adjacency representation is a
tuple of outgoing bitsets. The `receipt.json` is deterministic.

## Exhausted universes

1. All labelled tournaments of orders1 through5: counts1,2,8,64,1024, total1099.
   The checker verifies internal regularity and odd sizes of stable cells,
   crossing equitability, coprime crossing uniformity, SCC restriction of
   refinement, error additivity, and zero error. All applicable one/three-cell
   groups are compared with literal automorphisms.
2. All labelled regular tournaments of orders1,3,5,7: counts1,2,24,2640.
   Pair orientations are recursively exhaustive. The sole pruning condition
   is that each endpoint's target outdegree remains attainable from its
   unoriented incident pairs. Complete label sets are partitioned into full
   permutation orbits. Every representative has all n! maps independently
   tested for arc preservation. The class counts are1,1,1,3. At order7 their
   group orders are3,7,21, and their errors14,0,0.
3. All2744=14^3 crossings on three specified directed triangles. Each pair
   crossing can be either uniform orientation or one of six permutation
   matrices or six complements. Of these,434 have one stable cell of size9
   and are explicitly OUTSIDE the theorem's one/three-cell domain. The other
   2310 have three stable triangles and every signed-kernel automorphism
   agrees with a literal degree-pruned arc-preserving backtracking search.
   Their group orders have presentation counts648,1512,144,6 at orders1,3,9,27;
   their errors are18 on648 presentations and0 on1662.
4. All216=6^3 degree-one matching configurations on three triangles, independently
   quotiented by all27 literal within-cell rotations. There are16 ordinary
   classes, with12 groups C3 and4 trivial groups; their total error is72.
   This reproduces an inherited family count and adds its error accounting.
5. Every9! permutation is literally tested on each explicit balanced and
   unbalanced nine-vertex control, finding Aut sizes3 and1. Both are native
   tournaments with three intrinsic degree-separated cells. The rigid one is
   strong and has error18.
6. All eight independent choices of generator signs and all27 rotation
   parameters are checked on every one/three-cell configuration in item3,
   verifying covariance of the signed equations under generator replacement.
7. Thirty-six ordered joins of the empty, singleton, two-chain, directed
   triangle, regular-seven hostile, and rigid-nine hostile are checked for
   stable-cell restriction and error additivity; their maximum order is18.
   Automorphism checking here is literal degree-compatible backtracking, not
   the holonomy formula or an assumed SCC product.
8. The inherited two-prime-family arithmetic gives error lower bounds160 and
   8568 for p=5,7, and the displayed order-ten lower bound arithmetic is388.
   This item evaluates proved inherited formulas; it is not a new p=5/7 census.

The values E_s,n=0 through6 and E_s,8=0 rest on the mathematical proof in
THEOREM.md, not on the small finite universes above. The order-seven input
is the complete regular-seven census. Its transfer to E8=28 rests on the
proved SCC identity. The family formula for arbitrary r and the lower bounds
on unrestricted error coefficients rest on the displayed proofs and scoped
inherited formulas, not on extrapolating finite data.

## Reproduction

From this directory:

    python3 check.py > normal.out
    python3 -O check.py > optimized.out
    cmp normal.out optimized.out
    python3 replay.py

`replay.py` verifies the SHA256 hash of every manifest-listed payload file,
copies precisely those files into a fresh temporary directory, then runs
normal and optimized Python there and compares both generated receipts with
the pinned receipt. It also compares both stdout streams with the pinned
normal/optimized outputs. Temporary directories are removed afterward.
No parent-atlas file, repository checkout or unrelated scratch code is needed.

This is a finite mathematical control suite, not a benchmark. Independent
adversarial verification is maintained in the separate audit packet.
