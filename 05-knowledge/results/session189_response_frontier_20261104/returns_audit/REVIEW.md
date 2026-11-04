# Independent review: four-core joint-entry adapter

**Verdict: PASS for the stated all-order construction, exact 48-field response, observer lattice, physical bulk update, and fixed-component spacing obstruction.** No two-branch return closure or global Hamilton-spectrum conclusion is approved or implied.

The reviewed source hashes are recorded in `audit_receipt.json`. This review did not edit the producer packet. It independently read the incoming joint-entry definitions and JE3--JE5 at commit `8f14b96ab020f7d5ceeb88e7bda413bb6f8bc9b7`, inspected the all-order proof and gap enumerator, and used a separate simple-path DFS implementation and rational reconstruction.

## 1. General-core formulas

The condition m>=1 is essential and is retained. There are two disjoint classes of Hamilton paths in D_m:

- If the path starts at w, every Z vertex must occur before the first U vertex, and the remaining path is a Hamilton path of U, contributing h.
- Otherwise w is immediately preceded by r. Removing w splits the U vertices into an ordered two-bin cover with first bin nonempty and ending at r. Each Z vertex independently belongs to the initial transitive segment of the first or second part, giving exactly 2^m choices per core cover. This contributes 2^m f.

This proves H(D_m)=h+2^m f directly. Starting at w or deleting it gives h as claimed. The cut-after-w injection from Hamilton paths to covers ending their first bin at w is surjective: a failed join would require the second bin to start at r, forcing the first bin to be {w}, but D_m-w=TT_m->U has no Hamilton path starting at r. Hence E_w(D_m,2)=H(D_m).

Every simple old v-to-w path lies in U until its final r->w edge; its complement is Z followed by the core complement. Therefore the asserted joint identity follows with the same-bin meaning of J. The inherited N^3 dial identity then gives 5h+6*2^m f and the exact restored old starting counts.

The strongness arguments are valid for every allowed nontrivial strong U. The generic formulas were separately tested on every strong labeled core of orders3 and4, every retained root, and m=1,2: **204 literal core/root/clock cases**. These tests support, rather than replace, the preceding all-order argument.

## 2. All 48 named coordinates and the all-order certificate

The marked four-core rows (2,4,9,3), h=5, f=6, and the three joint values1 agree with direct enumeration and the inherited JE4 atom. The physical construction keeps all old arcs.

The exterior-gap proof is complete. For each fixed exterior word and bin cut, allowable Z gaps are distinguishable positions. A backward exterior adjacency within a bin forces that gap to be occupied; a forbidden boundary gap is removed for an endpoint mark. Assigning the m ordered Z names to k gaps produces k^m choices, and inclusion-exclusion on the required gaps gives the stated formula. Empty exterior bins contribute one gap and are handled correctly. No unmarked exterior-only path or empty-bin case is omitted. Every relevant field is therefore an exact linear combination of b^m for 0<=b<=10.

The producer's gap enumerator preserves both endpoint constraints jointly in W2, while repeated-name cross-bin coordinates remain zero. The same-name J used for entry is not confused with the cross-bin W2 observer.

A separate DFS enumerated all simple directed paths by their support, start and end, then formed all bin products on disjoint complementary supports. It reproduces **all 48 coordinates at m=0,1,2,3**, orders8 through11. A separate exact rational inverse of the 11-by-11 exponential evaluation matrix recovers **all 48 complete base0..10 coefficient rows** from the m=0..10 DP samples and matches the independently derived symbolic table exactly. This inverse check is not used as a substitute for the combinatorial span proof above.

The formulas for F_a,F_b,F_c and the three deletion counts are correct. The base-zero boundary is real: H(X_0)=57, whereas the positive-m formula would give61. Thus the positive-mode family and every later lattice assertion concern m>=1; the m=0 data remain retained separately. All 48 rows, including zeros and all canceled bases4 through10, are preserved.

## 3. Divisor-breaking and exact sibling spacing

For m>=1, H_m is1 modulo3. The a wrapper is1 modulo3, and b,c give0 and2 in opposite orders according to clock parity. This proves access to both nonzero residues from the same physical component, not the two desired exact siblings.

The strict ordering 0<F_c<F_b<F_a<2F_c is valid. The smallest same-clock gap is F_a-F_b=64*2^m-18. The nearest possible different-clock gap is bounded below by 2F_c-F_a, and

    (2F_c-F_a)-(F_a-F_b)=38*3^m-32*2^m-27>0

for every m>=1; its first value is23 and its consecutive increments are positive. Multiplying by the smallest allowed dyadic factor2 proves the exact minimum128*2^m-36, attained at a,b and k=1. Consequently no two wrappers of this fixed component have difference2, regardless of the two clock values. The proof does not require a finite clock cutoff.

As separate controls, the audit enumerates all three marks and clocks1 through12 for each m=1 through20, finding exactly the asserted minimum every time.

## 4. Integer lattice and modular observations

For the positive-m three-mode vector, every displayed row lies in Z x 2Z x 38Z. The three supplied integer row combinations generate (1,0,0), (0,2,0), and (0,0,38), so the equality of lattices is exact, not just a rational-rank assertion or a divisibility upper bound.

The bulk update diag(1,2,3) preserves this row lattice. Equivalently the kernel of the instantaneous modular observer is forward invariant. Thus adding every future bulk observation does not distinguish additional arbitrary modular mode weights. The image size is exactly

    N*(N/gcd(N,2))*(N/gcd(N,38)).

The audit also enumerates every input triple modulo each N=1,...,12, comparing all48 observed coordinates, image sizes, and one-step kernel stability. These finite checks agree with the unbounded lattice argument. The scope remains this linear three-mode system, not nonlinear representations or arbitrary tournaments.

## 5. Physical bulk transport and arithmetic safety

The new last Z vertex beats exactly a,b,c,r and loses to w, all previous Z vertices, and A,B,C. This is precisely the specified birth incidence. The recorded canonical-to-physical permutation merely accounts for shifting the displayed newborn labels; all original physical names and arcs persist. This is a legitimate source-preserving operation, not replacement by a scalar-equivalent graph.

The C++ subset recurrences count starts and ends in the correct directions. At m<=10 the old order is at most18. Every retained count is at most the number of permutations with a chosen cut, (n+1)n!, which is below 2^64, and all bitmasks fit. No installation is required.

## 6. Reproduction and limitations

    python audit.py > normal.stdout
    python -O audit.py > optimized.stdout
    cmp normal.stdout optimized.stdout

The audit records complete independent values, recovered coefficient rows, generic-core controls, modular images, spacing controls, and source hashes. All checks are explicit exceptions and remain active under optimized Python.

No mathematical correction was required in this packet. Its scoped conclusions are approved. It correctly retains the JE3 no-second-ready-return obstruction, does not identify graphs from Hamilton counts alone, and does not claim targets for a saved frontier source or reusable3H+2/3H+4 closure.
