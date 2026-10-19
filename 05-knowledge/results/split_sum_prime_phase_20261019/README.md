# Quadratic trapping and Fermat divisor balance

Research packet, 19 October 2026. No claim of global novelty.

Read [the mathematical note](PRIME_BRANCHING_PHASE_AND_FERMAT_BALANCE.md).

The main proved reduction concerns Chua's greedy **split-sum** process,
not the original Euclid–Mullin sequence. Once q is the smallest missing
prime, it is selected exactly when (−N/q)=1. Omission therefore forces an
eventual tail entirely in quadratic residues modulo q. At q≡3 mod4 this
is confinement to one Paley-tournament out-neighborhood. A first escaping
output forces q on the next step. Published divisor-support input is
explicitly credited to Booker–Pomerance.

The independent elementary theorem classifies every uniform binary
subset-product distribution on a cyclic group. At Fermat primes it gives
a sharp layer criterion for divisor multiplicities. Minimum-size balanced
sources have the wrong quadratic phase; generating sources need one more
prime factor. At q=65537 the exact smallest possible largest prime factor
of a balanced squarefree source is 327941. An explicit extremal generating
source has 16 prime factors and exactly one unordered qualifying split.

The complete kernels at q=5 for sources N=42 and N=66 are identical,
although only the latter generates q. This prevents replacing the target
phase with any statistic of that kernel. The [kernel-recovery companion](KERNEL_PHASE_RECOVERY.md) identifies the precise regimes where reflection symmetry and the period subgroup recover it, and separates remaining root-class ambiguity.

## Verification

Run:

    python3 -B verify_prime_branching.py
    python3 -B -O verify_prime_branching.py

Both pass 79,609 always-active exact checks and produce identical JSON.
`prime_branching_receipt.json` contains all candidate rejections, concrete
divisor witnesses, layer counts, and the 65537 extremal certificate.
`normal_run.txt` and `optimized_run.txt` record the two runs.

Prime completeness and the quadratic nontrapping condition remain open.
No Wilson-prime finiteness, normality assumption, or rough-sieve-to-prime
inference is used.

## Kernel recovery and independent reproduction

Run `python3 -B verify_kernel_phase.py` and repeat under `python3 -B -O`; its odd-order uniqueness guard is enabled explicitly. The [independent review](review/REVIEW.md) covers both notes and the appended odd-order theorem. Run `python3 -B review/independent_controls.py` and `python3 -B review/odd_order_controls.py`, then repeat both under `-O`. Compare stdout and generated receipts with the saved normal/optimized outputs and canonical JSON files. All programs use only the Python standard library.

The reviewed 553,351 independent checks and 28,196 supplemental odd-order checks supplement the mathematical proofs and finite certificates. No prime-completeness conclusion follows. `SOURCE_PROVENANCE.json` preserves all source hashes and documents status, cross-link, and presentation updates.
