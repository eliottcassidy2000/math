# The residual phase gives the exact departure from the repeated 27 block

25 October 2026. **PROVED elementary specialization of inherited parity-isometry and guarded-fuel mechanisms; FINITE-EXACT controls.** This records the first native branch after a periodic comparison ceases to agree. It supplies no general stopping-time or convergence theorem.

## Inheritance and object

Use shortcut Collatz T(n)=n/2 when n is even, (3n+1)/2 when n is odd. The [Session53 repeated-block theorem](../session53_g27_collatz_cooper_20261024/collatz27-repeated-block/REPORT.md) gives the 59-bit word

    w=11011111010110111011110100111011011111100111100010101000100,

with 37 odd bits and 96 ordinary operations. Write

    q=2^59=576460752303423488,
    a=3^37=450283905890997363,
    b=23q-27a=1100931843921811423,
    d=q-a=126176846412426125.

Its rational periodic comparison point is rho=b/d, with 8<rho<9 and d odd. The old whole-block fuel is K(n)=floor(v2(dn-b)/59). The general parity-isometry argument and p-adic fuel are inherited, including [the guarded pumping note](../collatz_guarded_pumping_memory_20261004.md) and the first-mismatch reasoning in [the Mahler source packet](../mahler_thue_morse_source_exclusion_20261019/THUE_MORSE_FINITE_SOURCE_EXCLUSION.md). The latter uses its own map convention; only the elementary equal-parity valuation argument is reused here, and it is reproved for the shortcut convention below.

The canonical hostile is the inherited two-block integer source that returns modulo seven without returning as an integer. The corrected near miss is to read block exhaustion as immediate disagreement: n=27 completes one block but shares two more bits. The least-used sidecar is the residual phase v2(dn-b) modulo 59.

## Exact bit departure and ordinary clock

For every positive integer n put

    ell(n)=v2(dn-b), k=floor(ell/59), r=ell mod59.

Then the actual shortcut parity sequence agrees with w repeated forever at precisely positions 0,...,ell-1 and first disagrees at position ell. In particular

    T^ell(n) mod2 = 1-w_r.                              (1)

The ordinary Collatz clock at that branch point is exactly

    C(n)=ell+37k+sum_(j<r)w_j.                          (2)

Proof. Extend T to rationals with odd reduced denominator by using the numerator's parity. This extension preserves odd denominators. The displayed rho has itinerary w repeated forever; its 59-step return follows either from the inherited affine block or a direct rational replay. If two such sources have the same parity, their next difference is their difference times 1/2 or 3/2, so its 2-adic valuation decreases by exactly one. If their difference has valuation zero, their parities differ. Since d is odd,

    v2(n-rho)=v2(dn-b)=ell.

Induction gives the exact agreement length and (1). Each common odd shortcut step contributes two ordinary steps, each common even shortcut step contributes one. Counting the common odd bits gives (2). The numerator dn-b is never zero for positive integer n because rho lies strictly between 8 and 9, so ell is always finite.

After the mismatching branch, this periodic comparison supplies no itinerary or termination control. It gives the actual next branch at a known source, not a rule that every future branch decreases a common rank.

## Sharp controls

- n=27 has ell=61, k=1, r=2. It reaches 53 at shortcut time 61 and ordinary time 100. The comparison bit w_2 is zero, whereas 53 is odd. Hence K=1 loses two lawful bits.
- n=27+2^59 has ell=59, k=1, r=0. It reaches 3^37+23, which is even, while w_0 is one. The same whole-block fuel K=1 gives a different residual phase and opposite next branch.
- The inherited two-copy source 209923714958077886127437625648218139 has ell=118, k=2, r=0, and departure endpoint 128083994787268184645564261786926256. Its residue return modulo seven remains separate from this native departure.

## Exact phase densities

Since d is odd, ell=j means a single residue class modulo 2^(j+1), namely the class agreeing with rho modulo 2^j and disagreeing at the next bit. Therefore for every j>=0,

    density(ell=j)=2^(-j-1).

For k>=0 and 0<=r<59,

    density(K=k and residual phase=r)=2^(-59k-r-1).

Summing the phases recovers the inherited density of K=k:

    (1-2^-59)2^(-59k).

Under uniform sampling from 1,...,X followed by X tending to infinity, the conditional residual-phase law is

    Pr(r | K=k)=2^(-r-1)/(1-2^-59),

independent of k. This is a density of initial cylinders, not independence of arbitrary orbit statistics. No prime-restricted density is claimed.

## Connection and next target

The exact map sends a positive source to its valuation fuel, phase and actual departure endpoint. It preserves the common-prefix length, next native branch and both clocks. Passing to K alone loses phase and the next branch; the examples above certify that loss. The phase restores precisely this consumer, not the entire post-departure orbit.

This is a useful comparison with the planar quartic/sextic proof: two first integrals there lose scalar Jacobian data, while a whole-block count here loses a residual transition. The carriers and open conjectures are different; there is no reduction between them.

The next genuine Collatz obligation is to control or join the native trajectories after these phase-labelled departures. The formula gives no uniform return to this template and no coverage theorem for a finite template bank.

## Reproduction

    python3 04-computation/shape_anchors_20261025/verify_collatz_phase.py
    python3 -O 04-computation/shape_anchors_20261025/verify_collatz_phase.py

The standard-library checker directly replays the rational periodic point, every positive integer through 4096, and every phase 0..58 at block depths 0..4 with several positive lifts. It checks first disagreement, both clocks, the next branch, exact cylinder residues and finite counts. Explicit checks survive optimized Python. The proof, not this finite universe, supplies the all-source statement.
