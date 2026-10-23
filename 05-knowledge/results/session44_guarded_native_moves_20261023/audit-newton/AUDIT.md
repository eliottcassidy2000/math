# Independent audit of the ternary Newton compiler

Verdict: PASS for the corrected files pinned below. The all-ring interpolation, cyclic-wrap identities, sharp single-coordinate nilpotency and minimal linear-observer theorem are valid in their declared scopes. No native electrical predicate, source network, coupled move or all-size wiring is authenticated merely by this table encoding.

Audited producer files:
- THEOREM.md SHA256 a85633d6a3c2db83b7bab736d400772c6f899ffe731862a1c3e919ad2fbcd891.
- check.py SHA256 23697e70313b4db11f190adcb7caa5d85f9dbc6017720c5e56bbb9ff6289bf7f.

The producer lives in ../session44-ternary-newton. The two producer replays pass7548 optimization-safe checks and agree byte-for-byte. The independent integer-matrix/tensor checker passes78810 checks in both modes, also byte-identically.

## Proof review

1. The lower-triangular three-point evaluation matrix has integral inverse, not an inverse that depends on division by2 or3. Its tensor powers therefore give unique reconstruction over every commutative ring, including rings with zero divisors. The essential-axis iff statement follows from uniqueness and the fact that representatives0,1,2 form one complete cyclic axis.

2. Conjugating the literal permutation matrix by the Newton evaluation matrix produces exactly the stated cyclic-wrap matrix. Its cube is identity. The identities D^3=-3TD and the two even/odd power formulas follow over the integers, not merely modulo a selected modulus. A delta basis vector makes D^2 primitive at3. Consequently D^(2k) remains nonzero modulo3^k while D^(2k+1) vanishes. The theorem correctly limits this statement to one axis.

3. On the integral zero-sum basis (1,0,-1),(0,1,-1), the literal difference has matrix [[-1,1],[-1,-2]], determinant3. Thus it is invertible whenever3 is a unit. Constants and zero-sum functions split using1/3. The corrected wording 'splits into eigensummands' avoids falsely calling the base ring Z/8 semisimple.

4. For an arbitrary output set, literal equality of translated tables defines a subgroup of F3^r. An observation C is sufficient only if ker(C) is contained in that subgroup; quotienting by the entire subgroup is sufficient. Rank-nullity proves minimal rank, independently of Fourier inversion. 'Every B-value' is understood as every attained value in im(B); unachieved codomain points have no fibres.

5. The Newton support-index shortcut is genuinely false: f(a,b)=(a+b mod3), as a table with values0,1,2 in Z/9, has translation kernel {(0,0),(1,2),(2,1)} and minimal observation rank1. Its coefficients at(1,0) and(0,1) are both1, so their index span has rank2. Coefficient indices are not frequency characters.

## Independent exact controls

The audit does not import the producer checker. It verifies:
- integral conjugacy, inverse matrices and the cyclic cubic identity;
- exact even/odd matrix powers through k=12, including their minimum3-adic entry valuations and sharp vanishing modulo3^k;
- tensor inverse identities over the integers in dimensions1,3,9,27,81, authenticating all tensor basis vectors rather than sampled tables;
- bijectivity of the zero-sum difference at moduli2,4,5,8,10,16;
- all19683 ternary-valued tables on F3^2, with axis tests and minimal linear-observer rank recomputed directly by all nonzero row observations. The rank census is3 constants,96 rank-one tables and19584 rank-two tables;
- the full explicit nine-entry Newton support hostile, agreeing with the corrected producer receipt.

The original checker output bug was real but non-mathematical: a final one-axis loop overwrote the variable used to print the earlier two-axis hostile coefficients. The current producer freezes hostile_coefficients before that loop. The corrected printed table is independently verified here. The theorem's interpolation, essential-axis and period claims did not rely on the erroneous old output.

## Reproduction

    python3 -B independent.py > independent.normal.json
    python3 -O -B independent.py > independent.optimized.json
    cmp independent.normal.json independent.optimized.json

Run the producer in normal and optimized modes separately to reproduce producer.normal.json and producer.optimized.json. All input/output hashes are recorded in MANIFEST.json.

An optional multivariable augmentation-ideal extension was sent to the parent separately. It is not part of this audit's accepted producer statement or these source hashes.
