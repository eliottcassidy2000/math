# Joint dyadic control of a physical rooted profile

For every power of two, actual integer gadgets can prescribe L arbitrarily and g,Q independently as odd residues, while B is zero. The construction has quadratic order in the modulus bit length and respects the minimal-width binary convention.

Its actual strong singleton-prefix outputs have independently arbitrary odd H and R residues, with both endpoint counts and root-deletion count exactly one. This is simultaneous residue coverage, not exact-target sufficiency.

- [Proof](PROOF.md)
- [Standalone verifier](verify.py) and both receipts
- [Independent physical audit](audit/REPORT.md), checker and both receipts
- [Inherited rooted grammar](dependencies/ROOTED.md)
- [Earlier modular construction](dependencies/MODULAR_GADGET_UNIVERSALITY.md)

Run each checker normally and with `-O`. All files require only Python's standard library. The L=1 transitive-edge exception is explicitly excluded from the strong-output statement. The later self-converse C2 application is a separate addendum and is not part of this snapshot.
