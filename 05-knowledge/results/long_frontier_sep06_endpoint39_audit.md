# Independent endpoint-39 proof and full polynomial audit

**Status: VERIFIED — analytic proof and entire exact parameter certificate
PASS.** Root independently read the candidate proof, producer source and
the exact real-root input
[THM-4436 / complete factorial-row simple negative roots](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md).
The theorem's stated hypotheses include all parameters used here.

The audited statement is precisely the one in
[the endpoint-39 proof](long_frontier_sep06_endpoint39.md): for integral
g>=20 with gcd(g,39)=1, and arbitrary nonzero complex coefficients on
support (-39,2g-39,3g-39), the first nonzero constant-term moment occurs
at g or 2g. The normalized doubled response is strictly negative at
all six first cancellation phases. No raw complex moment sign or
general all-channel closure is asserted.

## Analytic audit

The charge equation modulo g forces g to divide 39m, hence to divide
m under the stated gcd. At m=g, 2n_beta+3n_gamma=39 gives exactly
(x+j,18-3j,1+2j), j=0,...,6, with x=g-19>=1. At 2g the exact fibre
is (2x+e,36-3e,2+2e), e=-1,...,12. The lower carry is present at
every allowed parameter. The two multinomial normalizations in the
proof are exact, with a common actual phase alpha*gamma^2/beta^3.

The scalar first polynomial is the THM-4436 factorial row for
(A,B,h,r,z,x)=(2,3,6,0,1,x), so its six roots are distinct, negative
and nonzero. Dividing the lower carry by the first constant coefficient
cancels six even factors and leaves the polynomial displayed in
equation (6). Thus the actual quotient response has polynomial
coefficients of weighted degree at most12, including its inverse term.
Ordinary monic reduction does not raise weighted degree.

The multiplication matrix entry bound deg M_ij<=12+j-i implies
deg c_k<=12k by principal minors; row and column sums cancel in every
determinant term. This is a mathematical degree bound before evaluation.
The strictly positive shifted characteristic coefficients therefore
exclude every nonnegative real eigenvalue, and evaluation at each
first root makes the actual doubled response such an eigenvalue.
The conclusion is nonzero doubled moment at exactly the original root.
There is no phase retuning, missing carry or amplitude relaxation.

## Full independent reconstruction

The standalone
[audit source](../../04-computation/long_frontier_sep06_endpoint39_audit.py)
imports no producer. At each of the73 distinct integer points x=1,...,73,
it enumerates the full actual count fibres and reconstructs their exact
multinomial weights. It uses SymPy polynomial inverse and remainder
operations and the Berkowitz characteristic algorithm. This differs
from the producer's polynomial-matrix trace/Newton recursion and its
selected rational principal-minor controls.

All six reconstructed characteristic coefficients agree with all six
frozen positive-coefficient arrays at all73 points. Since both sides
have degree at most12k<=72, these equalities certify the complete
parameter polynomial identities. This is a full independent identity
certificate, not extrapolation from a short sample table. There are
45 primitive support-return parameters among these73; the other points
legitimately check the coefficient identities and no first-mass claim.

The source separately checks that all258 shifted coefficients are
strictly positive, that the carry survives, and that the gcd hostile
g=21 really has an earlier support return at mass7. All744 always-active
gates pass. Normal and optimized audit outputs agree byte-for-byte;
the producer's421-gate normal replay also matches its frozen output.

Reproduction:

```bash
python3 -B 04-computation/long_frontier_sep06_endpoint39_audit.py
python3 -B -O 04-computation/long_frontier_sep06_endpoint39_audit.py
```

The [frozen audit output](long_frontier_sep06_endpoint39_audit.out) names
the producer output hash and the audit semantic hash. The session
manifest records all source/proof/output byte hashes. The inherited
real-root theorem is used as cited proved canon; no new external
literature or formal proof-assistant claim is made.

No mathematical correction was required. The candidate may be promoted
to an analytically proved, finite-exact, independently audited result
with the scope above.
