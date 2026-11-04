# Lossless retained data and exact conventions

`h_vectors.bin` contains32,768 consecutive unsigned32-bit little-endian integers, one for each six-name quotient code in increasing order. Its exact size is131,072 bytes. Quotient-code bits run through unordered pairs(i,j), i<j, in lexicographic order; bit1 means i->j. The thirty signature bits run through ordered pairs(i,j), i outermost and j innermost, omitting i=j. Bit k is h_ij, and the two high bits are zero. Thus the bank retains each complete named signature, including zeros. Exact tuples, not digests, are compared when constructing the16,384 converse classes.

`literal.json` retains:

- all188 coarse battery masks, in lexicographic order;
- all188 exact integer Hamilton counts for each of four selected twelve-name states, with quotient code and six-bit pair-direction mask;
- six exact decoder counts for each of all64 pair-direction masks of quotient81, in increasing mask order;
- the six decoder masks on actual physical names;
- a two-query strong ambient control, together with its exact core and ambient counts.

Physical pair i is{2i,2i+1}; pair-flip bit0 means2i->2i+1, bit1 reverses it. A coarse mask bit expands onto both physical names. The core's ports are x=12,y=13. In the ambient control, two further exterior names14,15 are added with augmented-core->14->15->augmented-core. Therefore the ports retain the core's exterior incidences, and their direction remains x->y. All labels are physical throughout.

The literal checker uses uint64 counts, with at most16 names. Every count is bounded by16!, which fits in uint64. Modulo-eight assertions are reductions of these exact integer counts, not floating-point estimates. The separate endpoint-refinement twenty-name hostile used modular-only counts; it is not part of this data packet.

Reproduction writes the two deterministic data files and emits source/data hashes in the receipt. It compiles the C++17 checker into a temporary directory and installs nothing. The manifest binds this description, proof, source, full retained data and both normal/optimized receipts. No binary executable is needed or included.
