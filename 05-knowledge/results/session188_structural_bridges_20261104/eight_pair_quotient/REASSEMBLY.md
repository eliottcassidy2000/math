# Restore complete eight-quotient certificate banks

Run `python reassemble.py` from any directory, then `sha256sum -c SHA256SUMS` in this directory. The complete 2,097,152-byte seven-name signature bank and 122,880-byte eight-name residual list are restored byte for byte, with all physical labels and zero coordinates retained. DATA.md specifies the bit and little-endian word conventions. These are finite quotient-certificate data, not a complete sixteen-name response bank.

Every chunk, raw output length and hash, and existing output is validated before any writes. Differing existing files are never overwritten. The mathematical sources and original data are unchanged. The standard-library verifier compiles its C++17 controls in a temporary directory and reconstructs both banks.
