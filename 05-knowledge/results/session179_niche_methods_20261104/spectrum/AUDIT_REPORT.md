# Independent audit: one-apex modular synthesis

Verdict: PASS. Audited on 2026-11-04 against REPORT.md SHA-256 `9c675656e578c2a83153431387edff4630681f253bfd5c16be68455bf906b470` and modular_apex.py SHA-256 `58735eba65cceb982aecfb1ccd1cc1a11cc6a0de9d4974d8de95696f3e713fd5`.

## Mathematical checks

The direct path partition is exhaustive and disjoint: paths ending at the new chain vertex account for old paths not ending at the apex; all others must use the new vertex followed immediately by the apex. The choice of the suffix subset gives 1+s, including the empty subset. Thus the three-register recurrence holds at the empty boundary as well as for nonempty words.

The three append gates are exact compositions of that recurrence. Over every odd modulus their linear parts are invertible, so the finite permutation semigroup is a group. Direct affine composition gives F C F^-1 C^-1 as vertical translation by 3/4 and C^3 F^-1 as vertical addition (3x+22)/8. At any fixed x, E^2 Z^-x translates vertically by 22/4. Combining with Z using 15*3-2*22=1 yields translation by 1/4. This is a unit modulo every odd modulus, including multiples of 3. C translates x by -2, completing the transitivity proof.

Starting at (2,1), reaching (2,r-2), then appending zero realizes the required residue. The first and last bits give an explicit directed Hamiltonian cycle. A simple shortest path on q^2 states has at most q^2-1 edges, each of length at most three, yielding order at most 3q^2. The q=1 case is the directed triangle.

## Independent execution

Reran the supplied checker in separate normal and optimized Python processes. Both exited zero. Their outputs matched each other and the frozen modular_apex.out byte-for-byte. Output SHA-256: `b07c27b3e3330ff3eca5337fd20629996b13f06fdc42500718b788b3ff147410`.

The checker compares recurrence, cycle sum, and independently formulated Hamiltonian-path subset DP for all 1,023 words through length nine; tests full q^2-state reachability and all target residues for odd q through 101; and verifies the affine identities on all states for odd q through 25. Runtime assertions are explicit exceptions and remain active under optimization.

This audit independently checks the all-order algebra and reruns the independently formulated DP controls in the supplied program; it does not claim a separately authored second modular-synthesis program, a literature priority result, exact spectrum completeness, or any fixed bound on the lifted count itself.
