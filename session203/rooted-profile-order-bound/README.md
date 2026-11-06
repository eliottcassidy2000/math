# Rooted profile order bound

**PROVED, independently audited.** Every actual positive rooted profile (B,s), with B=H(U)+L-start_root(U), has an induced representative of order at most s+1+floor(log2 B). Delete only singleton strong-component tails; both B and s remain unchanged.

For each specified base X in the staggered criterion, any existing F(8,137) certificate can therefore be found among tournaments of order at most 27+floor(log2 X). This is a complete finite search cap, not a proof that the search succeeds or that F(8,137) holds.

- PROOF.md gives the two strong-source cover injections, weak-tail bijection, and boundary cases.
- verify.py and NORMAL.json/OPTIMIZED.json give exact producer controls.
- audit/ contains the frozen candidate and independent permutation-counter audit (PASS, no corrections).
- reviewed-originals/PROOF.md preserves the pre-promotion proof.

Reproduce with `python3 verify.py` and `python3 -O verify.py`; their outputs must match the supplied receipts. The audit has its own reproduction commands and hashes.
