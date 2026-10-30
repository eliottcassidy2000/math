# Independent interface audit

**Status: FINITE-EXACT PASS, 2026-10-30.** A separate checker reconstructs LRC survivors by endpoint/midpoint evaluation rather than the producer's interval-intersection algorithm. It checks 500 fixed-seed random source rows, affine boundary germs, and the finite nonemptiness gate. For composite symmetry it explicitly generates six permutation groups, compares them with the stated congruence kernels, and checks marking counts by Burnside fixed-point moments. At order nine it enumerates all 9! permutations to verify the claimed two-closure exactly.

From the packet directory run:

```sh
python audit/check_interfaces.py .
python -O audit/check_interfaces.py .
```

Only the Python standard library is required. All checks use explicit exceptions and remain active under optimization. Publication reruns in both modes give byte-identical `check_interfaces.out`.

The all-parameter arguments are in the lane notes; this audit is a finite independent control, not a substitute for those proofs or a major-open-problem closure.
