# Naming and source scope

The phrase Hammersley point set is used only to identify the auxiliary grid
{(A/2^k,rev_k(A)/2^k):0<=A<2^k}. No discrepancy theorem from the literature is
needed for this packet's elementary estimates.

A primary-source definition is in:

Peter Kritzer and Friedrich Pillichshammer, *An exact formula for the L2
discrepancy of the shifted Hammersley point set*, author-hosted paper:
https://www.dmg.tuwien.ac.at/nfn/KriPill_UDT06.pdf

The paper defines the binary radical inverse and shifted two-dimensional
Hammersley sets on pages1--2. At zero shift and N=2^k, the radical inverse
of A is exactly rev_k(A)/2^k; exchanging the two coordinates gives our grid. We import this identification only, not the paper's L2 estimates,
which are a different discrepancy norm.

The actual-source map and run-error derivation are proved and checked
locally in PROOF.md. The correction
from k/B to(k+1)/B for arbitrary real rectangle cutoffs arose from the explicit
k=1 hostile control; the original grid-aligned claim survives.
