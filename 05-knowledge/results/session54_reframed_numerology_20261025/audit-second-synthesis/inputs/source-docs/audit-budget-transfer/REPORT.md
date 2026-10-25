# Independent audit: finite digit-budget transfer supplement

Verdict: PASS. No correction required.

Frozen target: session54-fractal-budget-transfer, manifest SHA256 5d9519992f6ff49e9e8765a7d50928634179f890492586a60fa8818aa76b841d. This is a separate audit; the producer packet and original cross-base audit were not modified.

## Proof audit

The state coordinates count exact debt, not cumulative debt. At depth one the unit least digits are1,2,4,5; the first pair has debt zero and the second pair debt one. Thus the initial vector is (2) for budget zero and (2,2,0,...) for positive budget. No depth-zero clock is defined or used.

On appending one base-nine digit, each existing state contributes three lifts to itself and, if permitted, three to the next debt state. Consequently T=3(I+N), with N the nilpotent subdiagonal shift. Overflow is discarded. The truncation imposes no boundary correction on the allowed coordinates. Since N^(b+1)=0, the binomial matrix expansion and Pascal identity give v_r(j)=2*3^(r-1)*binom(r,j), including b=0 and j>r via the zero convention. Summing gives every stated total.

For j>=1 the series identity sum_{r>=1}binom(r,j)x^(r-1)=x^(j-1)/(1-x)^(j+1) is correct, and the j=0 term must be handled separately. With x=3z this gives exactly the reported F_b. For b>=1, multiplying by (1-x)^(b+1) leaves a polynomial whose value at x=1 is2: only the j=b term survives. Therefore the pole at z=1/3 has exact order b+1 and the reduced denominator is precisely (1-3z)^(b+1) up to constant over Q. For b=0 the numerator is2 and the pole is simple. In particular F_1=(4-6z)/(1-3z)^2 is correct. The field qualifier is necessary; modulo2 the scalar sequence vanishes.

N is a single nilpotent Jordan block of size b+1. Scaling its nonzero subdiagonal entries by3 preserves a single-block similarity class for T at eigenvalue3. The scalar observation sees the top direction because the all-ones row applied to N^b v_1 equals2. Thus a lower-order scalar constant-coefficient recurrence cannot explain the sequence; equivalently the reduced generating-function denominator gives the minimal characteristic polynomial (X-3)^(b+1). The asymptotic leading term is 2*3^(r-1)*r^b/b! for fixed positive b, with the stated error. Budget zero remains exact and separate.

The polynomial factor reflects Jordan multiplicity, not a changed exponential growth rate. The report correctly does not use cylinder-count asymptotics alone as a Hausdorff lower bound. The inherited no-defect Cantor subset gives dimension at least1/2, while the countable union over finitely many defect positions gives dimension at most1/2. The inherited parity-sheet logarithmic similarity transports the result. Haar measure zero follows because polynomial growth in r is dominated by3^r. Fixed b is essential throughout; growing b changes the problem.

For b>=2 these are combinatorial languages only. No coefficient congruence modulo3^(b+1), ordinary-exponent infinitude/finiteness, or full exponent classification follows or is claimed.

## Independent verification

The new checker does not import producer functions. It enumerates unit digit strings directly through r=6, rather than scanning modular exponent orbits. For budgets0 through8, generic matrix-vector multiplication through r=30 agrees with those digit counts and all closed-form exact states. Exact rational Gaussian elimination gives full (b+1)-by-(b+1) scalar Hankel rank and verifies the resolvent generating function at four rational arguments for each budget. Weighted differences verify annihilation by (E-3)^(b+1) but not its previous power. Top nilpotent directions, budget-zero initialization, forbidden cyclic wraparound, missing shifts, and wrong initialization are explicitly checked.

check.out and check-O.out match byte-for-byte. These finite tests support, rather than replace, the all-budget proof above.

## Provenance and closure

Before copying, the target manifest was checked against the assigned hash. Its source cross-base-proof.md was verified byte-for-byte against the original frozen Session54 REPORT.md, with unchanged parent manifest SHA256 92db3b170771188fac5c9d55faec8c4c5d18e4fc8f34d93c74b78abfd0794e95. Source content hash is2326e437bd53a6c0bfc2e2b9bdabb5e6bb0dcf6d9049341be484d46a6516f357. This supplements, rather than duplicates, the original independent audit of the inherited clock and dimension proof.

frozen-target contains only the target's manifest-listed files plus its manifest. The target replay and hostile programs were actually executed from a temporary isolated copy and unrelated working directory in normal and optimized Python. Both modes recomputed producer normal/-O fixtures and rejected corrupt, missing, and unexpected nested members. See packet.out. No network, repository, or external source was used in those isolated executions.

Run this audit's replay.py for closed-file verification, independent normal/-O fixture checks, and a new isolated producer replay/hostile run. This audit is self-contained apart from the already-recorded historical provenance comparison.
