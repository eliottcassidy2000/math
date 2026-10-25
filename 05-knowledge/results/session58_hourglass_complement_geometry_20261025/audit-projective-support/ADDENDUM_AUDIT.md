# Independent audit: uniform-exponent addendum

Date: 2026-10-25. Verdict: PASS. The base supplement and older absolute-geometry packet remain unchanged.

For an odd prime ell and a=(ell-1)/2, the integer interval [-a,a] is a complete residue system modulo ell. It is equivariant for negation. Thus the inherited Gaussian gap bijection gives exactly the nonzero finite-field vectors modulo the subgroup {1,-1}, not the full scalar quotient. For k>=1, every projective scalar orbit is free and has ell-1 elements; partitioning into sign pairs gives exactly a gap-index classes over each projective point. A bijection therefore occurs precisely at ell=3. The explicit k>=1 restriction correctly excludes the empty-set degeneracy.

The support counts follow by first choosing r nonzero coordinates and then quotienting by the free sign action: (ell-1)^r/2 versus (ell-1)^(r-1). Full scalar multiplication preserves zeros but can change the absolute values of balanced representatives. Exact local valuations at nonzero j_i and truncated lower bounds at j_i=0 together prove gcd(e²,d)=product p_i^(2(a-|j_i|)). This proves that scalar projectivization can destroy arithmetic information, not merely increase the number of forgotten indices.

The e=25 hostile is correct: j=1 gives z=15+20i and d=600, with gcd(625,600)=25; j=2 gives z=-7+24i and d=336, with gcd(625,336)=1. Both map to the unique projective point of P^0(F_5). The endpoints 5,35 and 17,31 independently verify both gaps. Thus a naive replacement of the ternary projective bijection by another auxiliary prime is false already in one coordinate.

The auxiliary ell is correctly restricted to an odd prime and is not confused with the actual integer split primes. Uniform exponents are essential for the displayed full vector-space indexing. Coordinate-zero deletion still divides by p_i^(2a), which follows directly from the corresponding Gaussian rational factor. No integer-addition, phase-multiplication, Schur, or F1-descent compatibility is inferred.

Independent computation: independent_uniform.py uses positive Pythagorean legs, without Gaussian generation, for every ell=3,5,7 and k=1,2. It verifies the full distribution of truncated gcd values predicted from all absolute exponent vectors. This includes the largest center e=5³·13³=274625, where all 24 gaps are enumerated independently. It constructs full scalar orbits and partitions them into sign classes, confirming the fiber sizes 1,2,3 and total multiplicities. Producer normal/-O replays and independent normal/-O outputs are byte-identical. A wrong fiber-size mutation is rejected by explicit RuntimeError in both modes.

The addendum's claims are proved from the inherited Gaussian theorem; no new external theorem or primary paper is needed. Source hashes are separately pinned in ADDENDUM_REVIEWED_SOURCE_HASHES.json. This audit does not change the prior base-supplement verdict or hashes.
