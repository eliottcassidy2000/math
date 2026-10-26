# Independent audit of the characteristic-three reduction packet

2026-10-26. **CONFIRMED**, with one precision clarification requested and applied by the producer: the quotient-stack fibre product is over [X/G], not over the coarse base P2. Auditor: parallel f1_extension agent, distinct from the reduction producer. This audit imports no producer Python code.

## Proof review

R1 is valid for every complementary pair H,K in F2^n and every inertia subspace I, without freeness. The dimension computation uses H+K=G and (I intersect H) intersect (I intersect K)=0. The missing-sheet dimension d is nonnegative, and r+d<=n follows from the target dimension. It includes I=0, H, K and G; none is silently excluded.

R2 keeps the correct four binary markings. The new characteristic-three triples are 136 and245. For each named complementary pair exactly one quotient-kernel intersection at each new triple has rank one; the other has rank zero. Hence d=2. This loses freeness but does not lose complementarity. Characteristic seven preserves the original incidence. No surface classification in characteristic three was imported.

The normality and conductor extension is valid in every odd characteristic. All line forms remain distinct primes. After an etale extension adjoining roots of local units, a single or double crossing is regular. At a triple the remaining local equation is w^2=a u^2+b v^2 with a,b nonzero. Its gradient vanishes only at the origin in odd characteristic, so the two-dimensional hypersurface is regular in codimension one and Cohen-Macaulay, hence normal. These local cases exhaust the arrangement. Distinct-prime valuations still prove the generic squareclass independence and the tensor image injection. The BOUNDARY B1 divisibility computation uses no characteristic-zero-specific step once 2 is invertible. Every code coordinate is active on both packets in every chart, giving the product of six chart roots as conductor; the unchanged seven-chart transition gives O_X(-7H). None of this requires a free kernel action.

R3 counts reduced rational points, not scheme lengths. For constant abelian deck group, Frobenius on the transitive geometric fibre is a translation. Root rationality is membership in I; product rationality is membership in J=(H+I) intersect (K+I). The three-event partition and its probabilities follow. Direct expansion gives E|X-Y|=2(1-2^(-d)) and E(X-Y)^2=2^(n-r)-2^(n-r-d). Summing pointwise yields 7q-11/4 in both odd-characteristic regimes, and 112q-251 versus112q-215 for the square discrepancy. These are not moments of total surface-count differences; the producer retains that distinction.

The quotient-stack statement is valid as the explicit 2-fibre-product identity [X/H] times_[X/G] [X/K] = X when G=H direct_sum K. The ramified coarse base cannot replace [X/G]. This was the only proof-scope clarification requested; no numerical or mathematical claim was retracted.

## Independent computation

`audit_f1.py` reconstructs all30 Fano designs by exact-cover recursion for Steiner triple systems on seven labels, rather than permuting the producer's label rows. The120 complementary pairs are disjoint designs. Avoiding the actual triple sets independently gives16 free kernels/32 both-free pairs away from characteristic3, and12/18 in characteristic3.

The twist census uses surviving **even character words on seven coordinates**, annihilator tests and all64 even twist words. It does not construct producer cosets or use the claimed rank formula to produce counts. It reproduces the new-triple distribution (8,2):8, (0,2):24, (0,0):32 and the original-triple distribution (8,1):8, (0,1):56.

Actual projective-point enumeration at q=3,5,7,9,25 verifies the stratum counts and pointwise moments for both norm-seven pairs, including extension fields. In this enumeration all coefficients are0,+1,-1, so the addition/negation structure of F_(p^2) suffices after using its standard normalized projective charts. All464 subspaces of F2^n for n=1,...,5 are independently checked against R3 for fixed complementary coordinate blocks; this includes nonfree and extreme inertia. Finally5,376 componentwise divisibility maximizations verify the conductor exponents for both named pairs in all seven charts.

Reproduce from the repository root:

```sh
python3 -B 05-knowledge/results/normseven_arithmetic_f1_20261026/reduction/audit_f1.py > audit_normal.json
python3 -B -O 05-knowledge/results/normseven_arithmetic_f1_20261026/reduction/audit_f1.py > audit_optimized.json
```

The committed receipts are byte-identical, and the checker contains no Python assertions. The all-field normality and conductor conclusions rest on the proof review above, not only these finite controls.
