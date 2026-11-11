# Reproduction and claim scope

From this packet's directory, with Python3 and SymPy installed:

    python3 reproduce.py

This runs six programs in ordinary and optimized Python modes, checks
both exit statuses, and requires identical standard output. Every
load-bearing control uses explicit exceptions, not assertions removed
by `-O`. The saved run is `FINAL_CHECKS.json`.

Individual programs:

    python3 gaussian/gaussian_bridge_exact.py
    python3 planar/check.py
    python3 planar/check_extension.py
    python3 topology/permutation_recursion_controls.py
    python3 degree/check.py
    python3 degree/state_check.py

## Exact finite universes

- Gaussian: primitive two-charge examples with exponents1..8;
  explicit parameter null families; covariance ranks1..30; covariance
  hostile; two-shear coefficient controls. 6,655 counted checks.
- Rotation: rows1..8 and product contexts of grades1..3, plus the row16
  rigid-debt hostile. 334 controls.
- Rational carrier and intertwiners: row shifts through11, exact finite
  expansions of specified positive carriers and delayed-pole hostiles.
  444 controls. The rational-cone theorem is proved algebraically for
  all rows, not inferred from these examples.
- Permutations: adjacent-transposition actions through degree7; two
  degree-six wreath actions; square-dihedral and marked-survivor controls;
  the21-vertex Paley7[C3] tournament and its explicit generators. The
  45,927 group order follows from the proved semidirect structure and is
  not claimed as a full enumeration.
- Degree inflation:48 explicit composed-shear cases, maximum degree100;
  180 controls. The degree theorem is all-N by its proof.
- State dimension: exact repair ranks at rows5..80, Gaussian tomography
  through12 coordinates, and scalar-blindness hostiles. 351 controls.

## Inherited inputs and unrun work

NC2/GMC(2), the Session254 homogeneous source theorem, and the finite
survivor/normal-generation mechanism are inherited at the sources in
`SOURCE_PINS.json`. No inherited Lean build, large census replay, or
general external-proof audit was undertaken. No actual G_m row16 lift
or universal rational source realization was computed.

The all-parameter statements have proofs in the notes. Computation
checks only the specified finite controls. Open conjectures remain open
after a successful run. This is not a repository-wide test-suite result.
