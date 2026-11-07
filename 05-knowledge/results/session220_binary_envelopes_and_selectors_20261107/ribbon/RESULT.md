# Session 220: constrained ribbon induction, exact spectral probe

Status: PROVED explicit infinite family; FINITE-EXACT spectral census. Chan–Pak and HYP-9473 remain OPEN. No completeness proof or constructive odd-count induction was found.

## Main positive result

Encode a span-2 ribbon letter by its backward-arc bitmask: 0={}, 1={1}, 2={2}, 3={1,2}. Let B be either word 02032333 or 02133232. For every integer k>=1,

H(B^k) = (63*265^k + 200*2^k)/263.

In particular the first four values are 65, 16825, 4457825, 1181322025. This is an exact infinite family, not a numerical asymptotic. It establishes that even integer eigenvalues do occur in the ribbon transfer semigroup, despite their absence from short blocks. It does not supply a pure doubling construction: the displayed physical seed has a nonzero 265-mode.

Proof. In the inherited exact 14-coordinate integer representation of THM-4880, each of these block matrices Q has rank 2 and characteristic polynomial x^12(x-265)(x-2). Exact matrix multiplication verifies Q^3=267 Q^2-530 Q. Thus H(B^k) obeys the associated order-two recurrence from k>=1. The exact initial values H(B)=65 and H(B^2)=16825 give the formula. The first value was independently checked by subset DP on the actual tournament, and the matrices were checked on 100 additional independent random words in each span against both profile DP and subset DP. The quotient representation counts the empty word as zero, so the displayed formula deliberately excludes k=0.

## Exhaustive finite screen

- Span 2, every nonempty block of length <=9: 349524 words. Exactly 16 matrices have eigenvalue 2. All have length 8 and are cyclic rotations of the two displayed words. Each has the characteristic polynomial above, in particular has no eigenvalue 1.
- Span 3, every nonempty block of length <=5: 37448 words. None has eigenvalue 2.
- Separate smaller probe: among span-2 blocks of length <=7 no matrix has both eigenvalues 3 and 1. This last probe is exploratory and not part of the main verifier's assertions.

Method: exact modular nonsingularity tests for Q-2I at prime 1009 reject all but 50 / 30 blocks respectively. A second test at prime 1013 rejects all but the 16 actual span-2 survivors. All survivors then receive exact integer determinant and characteristic-polynomial checks. Modular nonsingularity is a proof of rational nonsingularity, not a probabilistic inference. Matrix products in modular screens are reduced after every multiplication and stay within signed-int64 range. The main verifier reruns every word, rather than trusting the stored candidate lists.

## Scoped constrained-state consequence

There is no block in either finite box, and no prefix u, for which

H(u B^k) = 2^k*(H(u)+c)-c for all k>=0

with c!=0 and H(u)+c!=0. Such a scalar sequence has minimal recurrence polynomial (x-1)(x-2); every matrix representation producing it must have both 1 and 2 as eigenvalues. The finite screen rules that out. Equivalently, a fixed-block 2N+c rule cannot be iterated nontrivially even on a constrained physical state family in these block-length boxes. This is stronger than testing the rule globally across all prefixes, but remains a finite box statement. Constant sequences (H(u)+c=0), pure doubling (c=0), larger blocks, variable blocks, and nonlocal edits are excluded from this claim.

THM-4880(H)'s old eigenvalue-1 simplicity test only excludes nonconstant arithmetic progressions, not affine exponential recurrences. The present probe therefore tests a different missing prerequisite; it should not be presented as a consequence of that earlier test.

## Open boundary and stopping decision

The target requested was a positive construction covering all odd counts >=23, preferably via a finite invariant cone or backward reachable sets. No such invariant was established. The exact rank-two family reveals a signed doubling mode coupled to much larger physical growth and is not a rescue of the induction. The no-go applies only to iterating a single bounded block; it does not exclude collective reachable-set tiling or state-dependent grammars. Nor does it recover the already-known one-cross-edge or root-amalgamation construction.

A cone over R would still need an integer reachable-state coverage theorem; neither a spectral eigenmode nor a real-cone invariant supplies that. This is the exact unresolved step, not a denied tool action or missing data.

## Reproduction and provenance

From this directory run: python3 verify.py

The two span*_intrep.txt files, ribbon.py, and intrep.py are verbatim read-only copies from repository commit acbad6c8636cf8cbb675b4c947fa19e33ac0df23, directory 04-computation/experiments/ribbon_spectrum_20261104/. No shared repository edits or publication were performed. The remaining code and notes are owned session artifacts. Source theorem: THM-4880-ribbon-tournaments-bounded-span-back-arcs-give-a-finite-state-transfer-for-hamiltonian-paths-strong-connectivity-is-regular-the-span-2-series-has-a-rank-14-integer-representation-irreducible-on-nonempty-words-47-needs-span-3.md. Target hypothesis: HYP-9473-ribbon-completeness-and-the-ribbon-tiling-law-span-3-ribbon-words-realise-every-odd-number-except-7-and-21-which-would-prove-chan-pak.md.
