# Audit addendum: distinct AP triples in the Lagrange A-part

Verdict: PASS. This audit adds to AUDIT.md without requiring another large integer biclique search.

The twelve AP equations are exhaustive: each of four three-index subsets has three possible middle indices. The stated factorization of each balance was independently checked by plain integer coefficient convolution, beginning with the four source quadratic formulas and squaring them. This route does not invoke SymPy or a polynomial factorizer. It confirms precisely three occurrences of the linear factor 2u+3v, all on index triple {0,1,3}; the other factors are the listed nine quartics and three cubics.

The supplied certificates use genuine primes 7,11,13,17. Independent integer evaluation confirms all 148 projective values, including all twelve infinity values, and agrees with every stored residue. As an additional independent route, all nonzero vectors over the respective finite fields were checked, totaling 1632 nonzero residual values. The primitive-coordinate reduction proof correctly excludes every nonzero rational residual root, without restrictions on denominators.

The only surviving rational parameter line is (u,v)=(-3t,2t). Independent coefficient evaluation gives A=(900,900,12996,900)t^4. Homogeneity supplies the full-line identity; the supplied verifier also checks this symbolically. Therefore every surviving indexed AP has three equal values, never three distinct values. The origin is likewise degenerate.

The scope statement is correct: choosing original rows and columns cannot give an AP-K3,3 certificate requiring three distinct AP labels in both parts. Relabeling by within-part permutations, gauge translation, nonzero square scaling, and exchanging parts cannot evade this. This is not an exclusion of arbitrary regroupings of nine ambient square values, fresh additive coordinates, re-embeddings, other families, or MSOS in general.

Normal and optimized Python runs agree for both the supplied verifier and independent_ap.py. Deliberately changing the exceptional-line tuple or substituting the unsuitable prime 2 produces a nonzero failure under both normal and optimized Python. Explicit exception guards remain active; no assertions are used.

The final PRIMARY_CONSTRUCTION.md attribution correction is accurate: the numerical example is taken from Guy, while the parametric family is attributed to Jean Lagrange through Guy in Choudhry's introduction and equation (1.1), https://arxiv.org/html/2508.07806v1#S1.E1. No stronger separate inventor claim is needed for the numerical example. The AP result is independently derived, not attributed to the source.

The new file ap/input_receipts.json records hashes for the AP inputs, corrected primary-construction note, and original obstruction proof/verifier/certificate/source receipt. The four original obstruction files have unchanged hashes compared with the first audit. The old input_receipts.json remains a historical receipt for the earlier version of PRIMARY_CONSTRUCTION.md; the later AP receipt records the final attribution-only revision.

Reproduce the new checks with python ap/verify_ap.py and python ap/independent_ap.py, each also with -O. The ap/ folder contains the audited small proof/certificate files and replay outputs, not the article HTML. No original source or repository file was modified by this audit.
