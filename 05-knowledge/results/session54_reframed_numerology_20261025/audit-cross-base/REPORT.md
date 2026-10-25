# Independent audit: Session54 cross-base fractal

Target: session54-cross-base-fractal, frozen manifest SHA256 92db3b170771188fac5c9d55faec8c4c5d18e4fc8f34d93c74b78abfd0794e95.

## Verdict and scope

PASS for the elementary consequences, conditional on the explicitly inherited Session53 coefficient identities. No correction to the producer packet is required. The whole modular Frobenius theorem is not independently recertified by this audit. The old foundry's general method is correctly acknowledged; this is a coefficient-specific adaptation, not a novelty claim. No producer or repository file was edited.

## Mathematical audit

1. The exact clock is M_r=phi(9^r)=6*9^(r-1). LTE gives v3(2^(2t)-1)=1+v3(t), and odd exponents cannot give 1 modulo3. Therefore the order is exactly M_r, not just a divisor. Exponentiation bijects exponent classes and units.
2. Each choice of j defect locations has two unit choices for the least digit, in either digit class, and three choices for every other digit. Hence the count is 2*3^(r-1)*sum_{j<=min(b,r)} binom(r,j). For b=0 and 1 this gives the claimed 2*3^(r-1) and 2(r+1)3^(r-1). The period-normalized densities are respectively 3^-r and (r+1)3^-r. These are finite-window densities only.
3. The inverse limit of Z/M_r is Z/2Z times Z3. On the sheet k=epsilon+2z, Phi(z)=2^epsilon*4^z has distance ratio 1/3 in the usual 3-adic metric. Cross-sheet distance is 1. Keeping parity is necessary: the sequence of odd exponents 3^j converges to 0 in Z3, yet 2^(3^j) stays 2 modulo3, not 1.
4. The base-nine Cantor set has Hausdorff dimension exactly 1/2. Level-r covers give the upper bound. The independent uniform digit measure gives mass at most 3^-floor(m/2) on a ball of radius 3^-m, at most sqrt(3) times radius^(1/2), giving the lower bound. The two allowed unit least digits give two scaled copies. Fixed finitely many defect positions split into finitely many prefixes followed by a full Cantor tail, so have dimension 1/2 when nonempty. Their countable union, including the no-defect set, has dimension exactly 1/2. The inverse similarity transports the result. Cylinder probabilities tend to zero, giving Haar measure zero. No ordinary-integer cardinality follows from these facts.
5. For any interval of K consecutive nonnegative exponents, take the smallest r with M_r>=K. Distinct members have distinct residues. For r>=2, M_r<9K and (2*3^(r-1))^2=(2/3)M_r<6K; for r=1, the same weak bound follows from K>=1. Multiplying by r+1 proves the budget-one estimate. Thus the claimed uniform O(sqrt K) and O(sqrt K log(2K)) bounds and zero upper Banach density are valid. They do not interchange limits. Converting k to n=2^k gives the stated logarithmic index bounds, with harmless small-X conventions.
6. The positive-index restriction on the mod4 statement is essential. k=0 corresponds to n=1 and is valid in both intersections; n=0 is not a power of two and has c(0)=1. The report and its frozen source handle this boundary correctly.
7. The full expansion 2^13=8192=(12212)_9 proves k=13 survives both languages. k=6 passes depth one but fails globally. The witnesses k=M_r move with r; they cannot establish a fixed ordinary solution. Neither finiteness, infinitude, nor completeness of the finite exponent lists is proved. The report explicitly states this stopping boundary.

## Independent tests

check_independent.py uses admissible digit enumeration and baby-step/giant-step inversion, rather than the producer's orbit scan or exponent-lifting algorithm. Through r=6, its (budget0,budget1,budget2) counts are:

- r1: (2,4,4)
- r2: (6,18,24)
- r3: (18,72,126)
- r4: (54,270,594)
- r5: (162,972,2592)
- r6: (486,3402,10692)

Actual coefficients through n=120 are independently constructed from the defining s7 binomial sum and square-root convolution, not the producer's recurrence. All three support laws pass. Normal and optimized Python outputs agree exactly. Full ordinary digit checks through k=2000 reproduce E0={0,1,13} and E1={0,1,2,7,13}, strictly as finite guards. Negative and positive sheet coordinates pass finite-precision LTE checks; odd/even precision Cantor cylinders pass the Frostman inequality; interval constants pass all 1<=K<=100000. Hostile controls cover nonunits, nonprimitive bases, parity, finite/global confusion, and distinct defect budgets.

## Provenance and closure

check_packet.py independently pins the assigned manifest hash and exact recursive membership, then copies only manifest-listed files plus manifest.json into a temporary directory. The producer's replay and hostile programs run there from an unrelated working directory under normal and optimized Python, without repository imports. The replay itself recomputes and compares both normal/-O producer fixtures, including survivor lifting through r=10 and recurrence coefficients through4096. Missing, corrupt, and unexpected nested files are rejected.

Provenance checks compare all three coefficient-source files byte-for-byte with the corresponding frozen Session53 REPORT.md and their MANIFEST.json entries. Both historical sources are compared with exact git-show bytes at published commit bbe2ae278aeabdf707fdfca3c1bf0d3208564b03; the foundry source also matches original clone commit 1b91742ce55a202b73c80d56561ba56d00b503a6 and original worktree bytes. The foundry section3 supplies the ternary count/logarithm mechanism, not the new coefficient alphabet. These provenance comparisons require the local source folders and repository; the isolated producer replay and the independent mathematical checker do not.

See check.out/check-O.out for independent test fixtures and packet.out for closure/provenance/replay outcomes. All proof judgments above are separate from finite experimental evidence.
