# Reproduction and scope

Run from this directory:

```sh
python check_image.py > /tmp/jacobian_image.json
python -O check_image.py > /tmp/jacobian_image_optimized.json
cmp /tmp/jacobian_image.json /tmp/jacobian_image_optimized.json
cmp /tmp/jacobian_image.json receipt.json
python check_lifting.py > /tmp/jacobian_lifting.json
python -O check_lifting.py > /tmp/jacobian_lifting_optimized.json
cmp /tmp/jacobian_lifting.json /tmp/jacobian_lifting_optimized.json
cmp /tmp/jacobian_lifting.json lifting_receipt.json
```

Python 3 and SymPy (author run: 1.14.0). All checks use explicit exception
gates; no assertion-based test can disappear under optimization.
`check_image.py` takes approximately two minutes in the author environment.

The main receipt has 145 exact bracket identities, 11 coefficient-matrix
intersection tests, three multi-root residue controls and five critical-value
resultant cases. Matrix image intersections are independently constructed
from monomial Poisson derivatives, not by invoking the classification formula.
These are bounded regression controls, not all-parameter proofs.

The lifting receipt checks 15 numerical fibres (including zero), the generic
rational identity, its exact Laurent pole, a positive regular-lift control,
and the incomparable generic/special image generators. Nonexistence of a
regular lift is proved in Section 8, not inferred from finite searching.

A targeted canon search recovered substantial prior content:
- THM-3326-linear-in-z-unit-response-trichotomy-and-jet-torsion;
- THM-3348-linear-z-generic-puncture-response-and-one-root-valuation;
- THM-3412-hamiltonian-principal-part-differential-and-prufer-torsion-arms.
The first two already contain the nonzero-coefficient classification used
here. THM-3581-critical-value-multiarm-keller-compiler-and-A13-carrier also
uses deduplicated critical values in a different compiler. We therefore
make no priority claim for the critical-value operation itself. The present
packet isolates the A=0 Hamiltonian-image boundary and two explicit failures
of a specialization-to-regular-lifting inference.

No major-conjecture closure, new planar Keller chart, or higher-dimensional
Jacobian claim is established or used.
