# Independent collision-sieve audit

Verdict: PASS, within the explicitly bounded and conditional scope of the packet. No correction required.

Reviewed source: sibling session55-collision-sieve, certificate SHA256 32b59f10a2ce8ba1a5e5aeb9a8d83c06ba8bf36cea614b06b513d04eee6ee6e2, verifier SHA256 ea4e43bbb12a0fdc81dd5a00d248fce485c930bc0bae06528425568117b78b32. All files in its SHA256SUMS passed. This audit does not independently reprove the upstream Parker-prime classification.

- Independent enumeration of all 36 grid-pair differences yields eight rational directions, multiplicities 9,9,5,5,2,2,2,2, with eight extra factors of two. Squaring their product gives exactly the stated discriminant constant 2^16 and exponents 18,18,10,4,4.
- The degree-eight collision polynomial is necessary for square-entry integer arrays modulo each classified Parker prime. It is not a squarehood or existence certificate. Its formulation correctly preserves arbitrary offset content: dividing by gcd(u,v) loses a valid implication when a selected prime divides that gcd.
- Independent local linear-form counts reproduce both exact CRT fractions. The primitive fraction conditions on nonzero residue pairs at every selected prime; neither density is a probability model for the bounded gap banks.
- The mod-11 gap sets at normalized square centers 0 and 1 are respectively {0} and {0,3,8}. Their ordered pairs all collide. Scaling by nonzero e^2 and treating e=0 separately correctly proves redundancy for actual exact gap pairs.
- Independent integer Pythagorean-bank construction, without importing producer code, reproduces 347561 candidate pairs, all sequential and standalone rejection counts, and the exact ordered list/hash of 474 survivors. Every CSV opposite-entry root and gcd is independently checked.
- Most importantly, independent direct membership tests find zero pairs with u-v in the same bank across ALL 347561 pairs. Thus no pair reaches a modular sieve placed after the cheap difference lookup. The finite computation proves no speedup, and the packet accurately treats observed timing differences as run-specific, not guarantees. It does not claim a new integer search bound or an integer solution.
- Both the complete certificate/CSV replay and supplied controls pass under python -O. Runtime checks use explicit RuntimeError rather than removable assert statements. Supplied mutation tests check exact dictionary inequality, not full mutated-file entry-point executions; the actual verifier does use exact complete certificate and CSV equality, so this limited test design does not undermine the claimed replay guarantee.

Artifacts: independent.py, independent.json, replay-optimized.log, controls-optimized.log. No producer files or remotes were edited.
