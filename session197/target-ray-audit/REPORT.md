# Independent review of the prescribed-ray inverse algorithm

**PASS, no correction.** The finite copy bound, sign/integrality tests, bounded parity allocations, nonempty-module condition and explicit offset-65 classification are sound. The source and independent controls agree, including separately reconstructed physical source graphs.

## Frozen source binding

Under `session197-odd-cycle-path-polynomials/`:

- `TARGET_RAY_INVERSE.md`: 68b68319b19ac826221adf1d4ef17688cfbcf928fc89f71cf92e6c65882c8288.
- `verify_target_ray_inverse.py`: bf13c8b78bf19fe5fa9f14fd4f38e534702a21c40d69d9b39de2131042afeee3.
- Normal and optimized output: 42eabf561f5391ce6445101e7fb913a5ea70d81a76e91401a920cd16b5953db4.

Dependencies read and copied for replay:

- `RATIO_LATTICE.md`: a6f9c0b09d715865d01bbd62d95f2261c59e1f4a7eea2ed05cc09ac90aff60a3.
- `verify_ratio_lattice.py`: 3b40fd4e066f3b9d9eb3b6be44645a140035fd6df1d7adf9d9c29e757814a2e8.
- `UNIT_OFFSET_RIGIDITY.md`: 3c116086bfd697eaeaf5fcbf46308f1e31cdf653eb2a651053f01bce46ce374b.
- `verify_unit_offset_rigidity.py`: 756ea2c09f32d4830a05ebf92bee079812b3b013d9526b173a04ca2f70bbab0f.
- Graph helper `verify.py`: 10c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb.

The physical coefficient theorem was independently reviewed in `session197-unit-offset-audit/REPORT.md`. The full general Smith and cancellation audit is `session197-bridges/ratio-lattice-audit/REPORT.md`, SHA0c2d6d4f2566c93105dee79abd68ea3a78e5646c8ff72a502ecc1e5d1353c8e3; it independently checks 2,250 integer systems and physical graphs. This target-ray review concentrates on the additional bounded inverse step, not just those inherited cancellation facts.

## 1. The target equations are both necessary and sufficient

For fixed nonempty leaves P,R and a finite nonempty bank of nonempty graphs, the audited coefficient theorem gives A=A0 product b_i^(u_i+v_i) and C=C0 product h_i^(u_i+v_i). Here A0>0. The cut bound B_X(2)>3H(X) for |X|>=2, together with the singleton value 2-3=-1, proves C0 is nonzero. Thus zero target offset is correctly rejected.

The sign of C must be the sign of C0; C/C0 must be a positive integer, even when both are negative. The amplitude quotient A/A0 must likewise be a positive integer. These tests precede all logarithms and integer factorizations.

The binary coefficient has exactly the previously derived zero/same-sign/opposite-sign cases. In the opposite-sign case it is equivalent to the prime-valuation equation Mz=t for z=u-v, with negative denominator valuations retained. If both leaf e-values are zero there is no lattice restriction; if exactly one is zero or both have the same nonzero sign, cancellation is impossible.

For a fixed total-copy vector w, the conditions |z_i|<=w_i and z_i congruent to w_i modulo2 are exactly equivalent to nonnegative integral u_i=(w_i+z_i)/2 and v_i=(w_i-z_i)/2. Requiring both sums to be positive enforces the physical domain. Conversely these counts build actual ordered joins, and the exact coefficient identity supplies the entire family, not merely one scalar value.

## 2. Termination and no hidden padding step

Every b_i>=2, so X=product b_i^w_i>=2^(sum w_i). The explicit bound sum w_i<=floor(log2 X) is therefore valid without any bound on the order of individual bank graphs. This also handles X=1, for which w must be zero and nonemptiness prevents a solution.

The recursive factorization search is complete: for each bank index it tries every exponent until the residual X is no longer divisible by b_i or the residual Y by h_i. All residual quantities stay positive integers. Even when h_i=1, each successful exponent step divides X by b_i>=2, so no infinite loop arises. Duplicate bank graphs, coincident response profiles and a nontrivial valuation kernel merely create additional finite allocations; they do not invalidate completeness.

Every bounded parity-compatible z is checked, rather than only the first unconstrained Smith solution. An unrestricted lattice witness can fail the bounds or nonemptiness conditions. Conversely common padding that preserves cancellation changes the prescribed A and generally C, so it cannot repair those failures. The proof and implementation explicitly forbid such a repair. The displayed source-order bound counts the two fixed singleton quotient blocks, the leaves, the root module and every bank copy correctly.

Smith normal form is optional for the mathematical finiteness proof; it provides only an early exact rejection. The separately audited general Smith dependency validates that rejection in the supplied implementation. No polynomial-time claim is made.

## 3. Worked parity and nonempty boundaries

For P=K1, R=Tr4 and bank(K1,W) with (H(W),B_W(2))=(5,32), the exact cancellation matrix forces z=(4,-1). Fixing C=-65 then forces one W copy in total, placed in V. The singleton total is 4+2t, t>=0, because it must be at least four and have the same parity as its difference four. Hence the full amplitude set at that offset is exactly 15360*4^t. The proposed intermediate amplitude30720 has the required singleton total five and cannot satisfy the parity equation. This obstruction survives unrestricted lattice membership and is therefore genuinely part of the new bounded target problem.

At P=R=Tr2 and bank(K1), A0=18 and C0=1. Amplitude72 permits one copy on each side, while amplitude36 supplies only one copy in total. No parity allocation makes both modules nonempty. This correctly excludes the tempting empty-module shortcut.

All claims are relative to the stated bank and quotient. The offset-65 examples remain composite multiples of five, and no prime-spectrum closure is asserted.

## 4. Independent computations

The new `check.py` obtains small bank responses by literal enumeration of Hamilton words on every vertex subset and then sums the two-bin split counts. It directly enumerates physical copy-count pairs, computes the three coefficients without Smith normal form, and compares their complete prescribed-target support against the source inverse routine.

Controls include five banks, including response-dependent and rank-deficient choices; all 25 ordered leaf pairs from K1, Tr2, Tr3, Tr4 and W; exact arithmetic/sign failures; both zero and opposite-sign cancellation; and nonempty allocations. The independent enumeration tests 6,175 direct copy recipes and 35 pure target pairs. There are 8,284 source inverse cross-checks, all matching. The small target range bounds the total copies independently, so a missing larger-copy representation cannot explain an apparent negative in these controls. Memoization is applied only to the unchanged pure graph-profile helper for speed; no solver, Smith or allocation logic is replaced.

A separate enumeration checks every allocation at fixed offset-65 for singleton totals zero through28, reproducing the exact parity ray. The three source graphs are independently reconstructed by explicit indexed quotient blocks, not copied from the source constructors. Their literal rows agree with the source receipts.

`count.cpp` then counts those graphs using reverse endpoint subset dynamic programming, independent of the source's forward path recurrence:

- order16: H=46015;
- order18: H=184255;
- order9: H=217.

Optimized and debug C++ counts agree. Normal and Python -O independent receipts agree; copied source normal and -O receipts match the frozen output. Original source files remain unchanged.

## Artifacts

- `check.py`: 056c6f48e985c0fafd096ccda3b32d502f689a8c2fe33d2cda4077231c52a61d.
- `independent-normal.json`: 8b85758663775bcaf9d3d2b3532edb9fea96d6dfaceb8be2a9f1080bb16adc0d.
- `count.cpp`: 7fb4c581ed9f545fcb65562161b4fa8eb4fec76b21eef3b07d3e893f64ad7dca.
- `count-normal.txt`: 0ff003b30bfa3720fc6fea3b47389ae4c0b6951b0bcb8c14f3f5df0ab8245b12.
- `source-normal.json`: 42eabf561f5391ce6445101e7fb913a5ea70d81a76e91401a920cd16b5953db4.

No correction or outstanding blocker was found.
