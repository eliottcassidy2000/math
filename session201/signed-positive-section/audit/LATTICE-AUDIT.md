# Addendum: cubic-order sidecar audit

**PASS.** The sidecar's algebraic statements and physical caveat are correct. No correction is needed. This audit concerns only ALGEBRAIC-LATTICE.md and verify_lattice.py; the original physical-section audit remains unchanged.

Verified frozen hashes:
- ALGEBRAIC-LATTICE.md: 9a6737cb2515dfe87a32c28bcce9405c061f8533cc219fb453d460459b379e18
- verify_lattice.py: 10ffb58b582f91842eec4092d62ffadb7fcdd774d197451f25c8a9f4001acb77
- Both frozen receipts: 8d5ed7f0c3866dc65c20d745c01ef4edc0622b9d9ebdcca27990020539f6731f

## Mathematical review

The monic cubic X^3−X^2−X−1 is irreducible over Q by the rational-root test. Thus Z[theta] is indeed a rank-three order in the associated cubic number field, without needing to identify the maximal order. The coefficient transform Psi has determinant 1, and its conjugacy with the companion multiplication matrix proves the stated J action integrally. Theta's displayed inverse proves it is a unit.

The restricted raw pairing agrees with ell(Psi(v)Psi(w)) as a polynomial identity in six independent integer coordinates. Its Gram determinant −1 establishes perfection over Z. Commutative multiplication proves self-adjointness; it does not imply positivity on the whole lattice.

The quotient-ring identity theta^6(13−7theta)=3+theta is exact. The formal z and actual ZZ source both lie in the zero-f lattice, with the specified images. The singleton and edge images also check directly; the singleton's a=−1 is not mistaken for a nonnegative physical-cone coordinate.

For q=(5,1,5,1), the unique zero-f inverse is (1,0,4,3), with a=−2 and U<2d, hence nonphysical. The physical t.O predecessor (1,0,4,1) has f=2 and the same J image. This genuinely distinguishes invertibility on the selected integral lattice from invertibility within physical sources. The all-order section proof's forced faces provide the missing guard at precisely the stages where it uses them.

No maximal-order identification, Diophantine definability theorem, physical inverse from bare lattice membership, or universal target-coverage claim is made or needed.

## Independent checks

Both candidate replay modes agree byte-for-byte with the supplied frozen receipts. The separate `independent_lattice.py` uses SymPy 1.14.0 exact quotient-ring remainders and matrix conjugacy, unlike the candidate's hand-written polynomial reduction. It checks irreducibility, characteristic polynomial, integral unit inverse, coordinate conjugacy, the boundary identity, the full six-variable pairing identity, Gram determinant, and the unique hostile lattice inverse. Normal and optimized independent outputs also agree byte-for-byte.

No candidate or repository files were edited. No GitHub writes were performed; the parent owns integration.
