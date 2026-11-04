# Restore the complete prime-five audit records

Run `python reassemble.py` from any directory, then `sha256sum -c SHA256SUMS` in this directory. Both complete 130,076-byte original audit receipts are restored exactly, including all 264 labeled quotient maps and every one of the 96 retained profile certificates. No count, polynomial or labeled transport record is replaced by a digest.

The reassembler validates all chunks, all output hashes and sizes, and existing output contents before writing. It refuses conflicting existing files. The mathematical sources and original receipts are unchanged. The audit verifier requires the sibling prime_five_triples packet and an already installed g++; it compiles into a temporary directory without installation.
