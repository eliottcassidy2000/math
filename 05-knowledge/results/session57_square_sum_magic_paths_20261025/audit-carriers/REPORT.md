# Session 57 independent carrier and inheritance audit

Verdict: **PASS** for the scoped mathematical statements and complete publication packets. No change to either original packet is required. All audit work is confined to this new audit directory; source files were copied without alteration. No remote writes were performed.

## Results established independently

- The fixed-total 3-uniform hypergraph formulation is exactly an injective labelled eight-line grid embedding. Extra triples are permitted: an induced embedding is not required. The additive path formulation is also exact only with named-cell increment agreement and all-different labels. Literal intersection of graph edges is neither necessary nor sufficient.
- The center identity T=3c, opposite-pair sum 2c, and displayed three-parameter array work in every commutative ring without division. The written subtraction argument is valid also in characteristics 2 and 3. Independent finite checks cover every parameter triple over Z/nZ for 2<=n<=19; the ring proof, rather than these finite checks, establishes universality.
- The line-intersection graph has 22/28 edges. This is geometry only and cannot establish arithmetic existence. The cell-line incidence graph is connected with 17 vertices and 24 edges, hence cycle rank 8. The stated O(n^(-1/2)) upper bound for the density of G_n is valid for n>1 and gives asymptotic sparsity, not an obstruction to Hamiltonicity.
- Square-shell generation and canonical path reversal, without importing the supplied verifier, reproduce the counts for n=15,...,25: 1,1,1,0,0,0,0,0,3,0,10. Counts refer to labelled Hamiltonian paths up to reversal, not graph automorphisms.
- Cain's explicit arrays over F29 and F59 are correctly transcribed. The center-zero example remains included. In F59 the opposite sum 2 is a nonsquare, so magic opposite pairs need not be edges of the value square-sum graph. The F29 shadow triangle (0,1,4) remains a genuine failure of forgetting triple witnesses.

## Independent short-gap proof and inverse

Fix e>0 and a strictly increasing positive AP triangle (x,2e²,z) in the square-sum graph. Write p²=x+2e² and r²=z+2e². The AP condition gives p²+r²=8e²; reduction modulo four forces p,r even. Put a=p/2 and b=r/2. Then a²+b²=2e² and a<e<b. Therefore d=e²-a² is a positive integral member of D_e, and direct substitution gives

    (x,2e²,z)=(2e²-4d,2e²,2e²+4d).

Positivity is equivalent to 2d<e². Conversely any such gap yields the three pair sums 4(e²-d), 4e², 4(e²+d), hence the triangle. These constructions are inverse. No reliance on an unproved source parametrization is needed.

For normalization, a,b have the same parity. Set P=(b+a)/2, Q=(b-a)/2, g=gcd(P,Q), and (A,B,C)=(P/g,Q/g,e/g). This is a primitive positive Pythagorean triple with A>B. Then e=gC and d=2g²AB. The vertex pairwise gcd is 2g² and the pair-root pairwise gcd is 2g. These follow from pairwise coprimality and opposite parity in a primitive triple. Positivity is C²>4AB, equivalently A/B>2+sqrt(3). The equality case is impossible for integral e>0. No assumption that every magic diamond can be moved into this positive cone is made.

Independent direct opposite-square enumeration, primitive Euclid generation with all gcd scales, and pair-root AP enumeration agree for every 1<=e<=1105, covering 998 gaps. Every inverse, gcd identity, and positivity condition is checked. There are 68 positive petals for e<=500. At e=1105 the result is 13 total gaps, 5 short gaps, and no additive diamond. At e=5 one long gap survives but no positive petal. This refutes capacity-only and positivity-forgetting transfers, not existence at other centers.

## Source and status discipline

The 14 inheritance snapshots match their declared SHA-256 pins. All three carrier dependency files were checked against the original paths and then copied byte-for-byte for standalone replay. Cain v2 was read directly online; source pins and permitted scope are recorded in SOURCE_REVIEW.json. Definitions, elementary identities, and displayed witnesses are adopted; no classification conjecture is promoted. Other recovered connection-board material is assessed as correctly typed inherited context, not independently reproved in its entirety by this bounded audit.

MISTAKE-668's corrected automatic-nesting boundary is retained: monotonic nesting of graph or modular objects is not by itself a new structural correspondence. Neither density nor the repeated 3-4-5 vocabulary is substituted for an operation-preserving arithmetic map. The current packets make no such promotion. The exact context is snapshotted solely as audit evidence.

## Packaging and replay

The complete inheritance packet is FILES.json itself, its listed payload, and SHA256SUMS. With that convention, its manifest closes and all hashes pass. An initial warning about omitting FILES.json was resolved by confirming this convention; no manifest edit is necessary. Unlisted build_sources.py and other recovery scratch are excluded. The carrier payload consists of the four SHA256SUMS-listed files plus SHA256SUMS; exploratory scratch is excluded. Both original verifiers replay successfully in normal and optimized Python; inheritance output matches its stored certificates.

Run from any directory:

    python3 /path/to/audit/independent.py
    python3 -O /path/to/audit/independent.py

The outputs are byte-identical and match independent.normal.json and independent.optimized.json. There are 44,519 explicit exception-based checks and six rejection controls, including two malformed shared-label paths, fake Hamilton counts, a capacity/short-gap conflation, and opposite-shell/squared-total conflations. No substantive check uses assert. The independent verifier imports neither supplied verifier and now uses only the included snapshots.

These finite tests support the scoped elementary proofs; they do not solve the integer magic-square-of-squares problem, prove Hamiltonicity asymptotics, or establish any universal compression theorem.
