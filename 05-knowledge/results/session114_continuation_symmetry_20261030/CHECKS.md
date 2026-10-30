# Verification and reproduction

Run `python run_checks.py` from the packet directory. It runs twelve scripts normally and under `python -O`, requires successful exits, and compares their output bytes. Python 3 and SymPy are required; the arithmetic/interval/Collatz checkers otherwise use the standard library. The scripts use explicit exception gates, not optimization-removable assertions. `check_receipt.json` records the actual replay and output hashes.

## Universes and independent routes

- **Lonely Runner, fixed core:** 4,425 producer checks; 484 tails (1..480 and four larger controls); 1,867,031 exact candidate-clock visits. The core is also reconstructed from every wall and intervening cell. [Independent audit](lrc_audit/AUDIT.md): 1,077 gates including all 154 core walls, affine branch certificates and endpoint owners.
- **Lonely Runner, all q:** 1,364 producer checks, core windows q=3..40 at five exact radii, 261 complete maximum rows and large structured witnesses. [Independent audit](lrc_audit/GENERAL_AUDIT.md): 1,308 gates for q=3..80 using independently generated interval intersections and owner data. All-q scope comes from the Farey-basis proof, not from these finite ranges.
- **Collatz producer:** 12,192 word subgroup-membership tests, 13,936 synthetic affine tests through modulus625, 372 direct-walk mean comparisons and all 4,095 nonempty controllers of clock at most12. Of these controllers, 3,299 meet the primitive-root sufficient condition; the other796 are not covered by that asymptotic bound theorem.
- **Collatz independent:** [Proof audit](collatz_audit/AUDIT.md), 548,000 arbitrary affine membership tests, 1,016 order checks, and160 unbounded-integer prefix controls. The proof explicitly handles the infinite-layer limit and fixed exit. The tests do not establish convergence of Collatz trajectories.
- **Keller:** 3,977 producer gates. [Independent symmetry audit](audit/SYMMETRY_AUDIT.md) adds3,130 algebraic gates and proof review. Hostiles include the degree-one scalar exception, strict finite-bound alias, and mod3 degree-four two-shear map.
- **Tournament charts:** 38,707 producer gates across six finite wreath-product cases, with direct literal-permutation enumeration in the independent audit. Counts concern named parent subgroups. Neither tests nor proofs claim a universal output-cover complexity lower bound.
- **Aliquot:** 35,153 producer gates plus independent divisor-sum checks on30,000 source/multiplier pairs and42 ordered distinguishing pairs. Infinite splitting is proved only when every prime-dilation letter is allowed; the fixed-alphabet finite count is separately proved and checked.

## Audit corrections and boundaries

Prepublication review corrected a factor2 in the exponential parent example and an exponent index in the aliquot source family. The published values are K=2·3^(r−2) and n_p=2^p, respectively. The q=12,m=1,2 exceptions are explicit. No failed intermediate version is used as a theorem dependency.

The last session's old Collatz epoch-budget retraction remains in its own source and is not reused. This packet's proofs use the repaired predecessor statements and retain the limits of every inherited result.

No Lean build or CI run is claimed. Small exact computations are regression controls for ordinary written proofs; they do not certify arbitrary instances outside the declared scopes.
