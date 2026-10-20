# Sphaleron–Steinitz investigation

The main result is a typed reconstruction, not a new physical identification:

**degree-one sphaleron loop → integer winding subgroup → profinite completion Ẑ ≅ Gal(F₂bar/F₂) → supernatural order Ω.**

Read [SPHALERON_STEINITZ_BRIDGE.md](SPHALERON_STEINITZ_BRIDGE.md) for the hypotheses, explicit maps, preserved structures, obstructions, primary-source links, and relation to binary linear registers and universal UHF algebras.

## Files

- `SPHALERON_STEINITZ_BRIDGE.md`: sourced research note and self-contained abstract proofs.
- `verify_bridge.py`: reproducible symbolic and finite-field checks; requires SymPy.
- `verification.json`: checked exact output.
- `SHA256SUMS`: integrity checks for these deliverables.

From the parent directory, run `python3 session9-sphaleron/verify_bridge.py`. Its output must match `verification.json`. The optimized run with `python3 -O` has also been checked byte-for-byte.

## Verified scope

The symbolic boundary map has degree +1; reversing orientation gives −1. For every degree 1 through 8, finite-field Frobenius has operator order n, is a cyclic shift in a computed normal basis, and has the exact binary-necklace cycle counts. Primitive multiplication has order 2^n−1 and therefore misses every factor 2. Modular checks show coherent winding residues and the obstruction to a half of 1 in Ẑ. The infinite statements are proved in the note, not extrapolated from this census.

## Source-access scope

The Klinkhamer–Rupp 2003, Manton 2019, Milne, and Ghasemi–Kubiś PDFs and Stacks lemma 58.2.2 were inspected directly. Original 1983/1984 physics papers and Steinitz 1910 were checked through their official bibliographic/abstract records; the latter's historical full text was not used as a proof source. Ribes §4 definitions were checked through the search-indexed university-hosted excerpt after direct PDF retrieval failed. All substantial algebraic steps used here have explicit arguments in the note.

These files were produced outside the repository checkout. No claim of publication, independent human review, or a new electroweak solution is made.
