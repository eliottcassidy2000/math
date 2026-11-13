# Independent audit of the L/K norm, Jordan and phase classifiers

**PASS.** Reviewed sibling [anchor/PROOF.md](../anchor/PROOF.md), sections B-F,
including the all-order integral argument and the precise finite-computation
dependency for reachability of the plane. Reviewed file SHA256:
`2823f4f7a1e2f4ca7b770d8023c3d36fc8cd95cc47d40b95ecfb6b870e2161fe`. The scope is the inherited physical L/K language, not arbitrary
selective tournament moves.

The F49 model checks algebraically: for c=A+Bt with t²=4t+1,
Tr(c)=2A+4B, Tr(ct)=4A+4B and Norm(c)=A²+4AB-B².
Eliminating A,B gives 6u²+3uw+w². The discriminant must be nonsquare or
zero; nonzero square roots lie in F7 and do not have the requested
extension trace and norm. The zero-norm boundary is handled correctly.
The even/odd t-power cosets have norms1/-1 because t has order16 and
Norm(t)=-1.

For F11, ell=2x+6b+d is a left eigenform of eigenvalue7; its zero
plane is exactly the span of the displayed7/9 eigenvectors. The equation
2^j=6c/a has the stated Legendre phase. Off that plane the nonzero
Jordan coefficient B=9ell produces a unique residue modulo11 for every
chosen j modulo10, so both parities occur. At77 the common exponent
condition is exactly parity since gcd(48,110)=2. The proof retains
positive-word reachability of the b=0 plane before using this CRT step;
it does not infer physical reachability merely from residue CRT.

The independent script imports no anchor implementation. It checks every
local vector against its explicit L orbit, independently performs complete
positive L/K BFS at7,11,14,77, verifies each classifier on every vector in
the relevant cube, and confirms every claimed plane point is reachable.
The counts are337,1311,434,429807, with exactly12000 incompatible local
phase pairs. It also checks all49 norm-coordinate pairs and all4095
integral words through11 births. The exact defect-one language is K*L;
the first false-ready witness is KLLLKL at state(38,15,53), as claimed.
Normal and optimized runs are identical; see [receipt](audit_anchor.json).

The integral proof is sound beyond the bounded check: s=d-b becomes
x+b>=2 after L, and later K cannot restore s=1 because d has already
increased. Thus only a pure K prefix can precede a defect-one final L.
Every birth increases full H by x+2b>=2, which correctly rules out an
integer count return even when the complete residue state returns.
No correction requested.

Reproduce with `python -B audit_anchor.py` and `python -O -B audit_anchor.py`
in this directory. The sibling audit of our separate transitive-context
orbit is in [anchor/AUDIT_LOOPS.md](../anchor/AUDIT_LOOPS.md).
