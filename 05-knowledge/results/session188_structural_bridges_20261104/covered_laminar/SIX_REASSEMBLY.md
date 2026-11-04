# Lossless six-name modulo-four bank

Run `python reassemble_six.py` before the validators or degree diagnostic. This restores exactly 15,859,712 bytes in `mod4_six_vectors.bin`: all 484 named coefficient residues, including all zeros, cubics and quartics, for all 32,768 physical tournament codes. The retained values are canonical residues modulo four, not full integer coefficients.

The coefficient-major layout and physical labels are specified in `MOD4_SIX.md`. No multi-byte integers or endian convention are used. Every chunk and the complete restored bank are hash- and length-checked before any write. A differing existing output is never replaced.

The compressed chunks are a lossless representation of the complete bank, not a summary or equality hash. `SIX_REASSEMBLY.json` gives exact lengths and integrity hashes. Python uses only the standard library.
