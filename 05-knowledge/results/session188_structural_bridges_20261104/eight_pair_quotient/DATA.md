# Complete retained data and coordinate conventions

`t7.bin` has2^21=2,097,152 bytes. Byte c is the complete seven-bit vector t(R) for the seven-name tournament with edge code c; its high bit is zero. Bit i is the coordinate at physical quotient name i. Every vector is retained, including zeros.

`blind8.bin` has30,720 consecutive unsigned32-bit little-endian integers,122,880 bytes total, in increasing order. It is the complete labelled list for which every shared-cubic entry vanishes. Its four high bits are zero. Exact code values, not digests, are tested in the complete residual and orbit checks.

For either quotient order, edge codes enumerate pairs(i,j), i<j, with i outermost and j innermost; bit1 means i->j. In particular the first seven bits at order eight are the incidences0->1,...,0->7, and all later bits encode the deletion of name0 in sorted name order. The relabelling routine is used only to report the four orbits after the complete labelled enumeration. No physical response comparison relabels the attachments.

`census.json` retains the exhaustive coverage counts and initial residual codes. `check_endpoint.json` records the separate subset-DP verification of every retained residual and the complete histogram of its full start/end XOR masks. `orbits.json` gives all four representatives, complete adjacency rows, orbit and automorphism sizes, score sequences, primeness and their converse partners.

`literal.json` retains the exact modulo-eight physical eighteen-name counts for the specified prime example, its physical masks and full quotient endpoint parities. It also records570 literal small-query controls of the endpoint-inflation identity. These modular counts are not claimed to be full integer counts.

All source, both binary banks and complete JSON records are bound in the paired receipts and manifest. The driver reconstructs them deterministically; compiled executables live only in a temporary directory. The complete2^28 full response functions are not enumerated or claimed to be retained.
