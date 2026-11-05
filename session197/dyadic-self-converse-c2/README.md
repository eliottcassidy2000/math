# Exact dyadic image of self-converse C2 tournaments

For every k>=2 the exact Hamilton-count residue image modulo2^k is all residues1 modulo4, together with the exceptional residue of15. Every nonregular allowed residue has infinitely many distinct physical counts, with a witness of order O(k^2). The exception is the single regular five-tournament.

This is a residue theorem, not exact integer sufficiency. The construction uses actual weak start-one sources inside a strong diagonal C2 graph.

- [Sharp residue proof](SELF_CONVERSE_C2.md)
- [Verifier](verify_self_converse_c2.py) and both receipts
- [Independent audit](audit/REPORT.md), standalone checker and both receipts
- [Joint-profile dependency](PROOF.md), with its executable helper [verify.py](verify.py)
- [Structural self-converse classification](dependencies/SELF_CONVERSE_CLASSIFICATION.md)
- [Joint-profile independent review](dependencies/JOINT_PROFILE_REVIEW.md)

Run `python3 verify_self_converse_c2.py` and its `-O` form, then run `audit/check.py` similarly. The source verifier imports only the adjacent verify.py. All dependencies are included unchanged; no external code or data is needed for either check.
