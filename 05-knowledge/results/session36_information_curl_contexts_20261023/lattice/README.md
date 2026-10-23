# Session36: integral circulation geometry and observation depth

Frozen checkpoint2026-10-23. No repository or Git changes were made in this lane.

## Headline results

1. The integral metric alternating augmentation lattice determines a tournament up to cyclic rotation of its ordered strong components. It therefore determines H, and classifies strong tournaments. Smith type without the metric and real metric-operator spectrum without the integral root placement both lose H.
2. The entire norm2-root-pair response is exactly(c3,c4). Consequently every same-order skew conference tournament has identical two-root response. A genus3 coefficient separates the Paley8 H189/333 split, but that particular coefficient already fails at order12.
3. Joint root-frame coefficients at genus n-1 reconstruct the full simultaneous pair. Nevertheless every fixed root-shell genus cutoff fails to determine H at every sufficiently large order, even on strong tournaments with identical labelled scores. Moon interval chains supply the multiplicatively independent H-values behind this information-budget obstruction.

## Essential proof files

- `METRIC_CURL_RECOVERY.md`: fibre classification, integral isometry, matrix-pencil identity, exact two-root distribution, finite joint-theta recovery and Paley8 response.
- `FINITE_GENUS_OBSTRUCTION.md`: all-depth no-uniform-cutoff theorem, stronger fixed-score strong construction, eventual every-order corollary, and explicit order9 hostile.
- `RECOVERY_AND_ROUNDS.md`: older-source recovery, corrected Hessian boundaries, rank88 primary-source check, and five conjecture/test rounds.

## Essential reproducibility files

- `explore.py`: shared exact tournament, SCC, Hamiltonian and root-pair primitives used by both final scripts. Its exploratory main is not needed by the final tests.
- `verify_metric_curl.py`, `verify_metric_curl.json`, `verify_metric_curl.out`.
- `finite_genus_controls.py`, `finite_genus_controls.json`, `finite_genus_controls.out`.
- `VERIFICATION_RECEIPT.json`, `SHA256SUMS`.

The other exploratory files in this scratch directory are not needed for publication of these proofs and certificates.

## Run

From the current workspace root:

python3 session36-lattice/verify_metric_curl.py
python3 session36-lattice/finite_genus_controls.py

Python standard library and SymPy are required. The scripts use explicit exception checks, not removable assertions. Both ran under normal Python and python3 -O, with byte-identical JSON and stdout in each pair of runs.

## Exact test universes

- Every labelled tournament through n6:33867, complete restricted-form fibres and SCC rotations. Independent permutation Hamiltonian counts and direct root-pair formula checks through n5:1099 tournaments.
- All128 specified Paley8 conference switches, two independent genus3 counting paths.
- Sixteen conference controls at orders4,8,12,20; all root-star bases through rank5;132 exact matrix-pencil evaluations.
- Session35 Smith-type collision; explicit order12 single-coefficient hostile.
- Nearly transitive blocks a1..12;45 pairwise prime gcd checks;1024 distinct multiplicity products;232 interval-chain compositions through order13.
- Strong order9 fixed-score pair IC(2,2,4)/IC(2,3,3): same full root-pair response, H81/75.

## Audit and limits

Independent audit confirmed the fibre proof, the A-root isometry group, theta coefficient reconstruction, conference distribution, all128 Paley8 coefficients, the general two-root formula, Moon-chain product use, and eventual every-order horizon. It independently recovered the order9 hostile.

Root-shell means squared norm2. No fixed-genus assertion is made about the full theta series over unbounded norms. The sufficient genus n-1 is not claimed minimal. The eventual-order pigeonhole witnesses are existential and astronomically large; the order9 witness is explicit. No new literature-priority claim, LRC closure, Riemann-zero theorem, or rank88 tournament construction is asserted.
