# Restore the complete six-name signature bank

Run `python reassemble.py` from any directory, then `sha256sum -c SHA256SUMS` in this directory. This restores all 131,072 original bytes of h_vectors.bin exactly. DATA.md specifies the 32,768 unsigned little-endian words, quotient-code orientation convention, and all thirty physically named signature bits; zeros are preserved. This is a complete quotient-signature bank, not a full two-million-state response bank.

Every chunk and restored output is size/hash checked. All existing outputs are checked before writes, and differing existing files are never replaced. Mathematical sources and complete literal.json physical counts remain unchanged.
