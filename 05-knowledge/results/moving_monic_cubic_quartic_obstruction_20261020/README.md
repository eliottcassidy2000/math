# Moving entire cubic/quartic source obstruction

20 October 2026. **PROVED; independently reviewed PASS.**

For R=C[x] and global monic cubic/quartic generators h,k, the entire
source algebra R[h,k] has no generating pair with nonzero constant
(x,y)-Jacobian. The normal form is h=y³+a(x)y and
k=y⁴+b(x)y²+c(x)y; arbitrary global monic degree-three/degree-four
presentations reduce to it by a polynomial source shear of Jacobian one.

Mechanism: the immersion resultant is a polynomial unit; nongeneric
singularity regimes freeze the coefficients; algebraic degree reduction
selects one of two generic embedding classes; monic conductor division
and all-fiber cotangent invertibility force the ambient Jacobian scalar
to be constant. Three moving node periods then expose the degree-gap
contradiction, separately for each embedding class.

- [Proof and exact boundary](MOVING_CUBIC_QUARTIC_OBSTRUCTION.md)
- [Dependencies and primary-source scope](PRIMARY_SOURCES.md)
- [Exact symbolic controls](verify_moving_cubic.py)
- [Control output](results.json)
- [Independent proof and source audit](review/AUDIT.md)

Run from this directory:

    python3 verify_moving_cubic.py
    python3 -O verify_moving_cubic.py

The 64 producer checks run identically with and without optimization.
They include the exact degree-six conductor, raw resultant normalization,
all differential coefficient identities, second-generator recovery,
three singularity controls, a nonconstant immersed family with both
global embedding presentations, and a standard-class-only moving hostile.

Boundary: no generic-to-global cubic assertion is made. Combined with the
previously audited saturation dichotomy, an actual noninvertible generic
conductor-six source, under that dichotomy’s hypotheses, must have a
quadratic special fiber; the full
conductor-six sector and the planar Jacobian conjecture remain open.
