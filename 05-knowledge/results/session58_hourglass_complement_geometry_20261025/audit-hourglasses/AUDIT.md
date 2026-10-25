# Independent audit: overlapping hourglass anchor

Status: PASS for the structural and finite-exact claims in `/tmp/session58-hourglasses/THEOREM.md`. This is not a proof of existence or nonexistence for distinct integer-square hourglasses or full magic squares, nor an endorsement of unrelated repository theorems.

## Independent checks

- Derived all eight line identities from the explicit label map `(A,B,C,D,E,F,G,H,I) -> ((1,0),(-1,-1),(0,1),(-1,1),(0,0),(1,-1),(0,-1),(1,1),(-1,0))`, where address `(i,j)` means `E+i*x+j*z`. Verified the H and V edge sets individually, rather than comparing only multisets. Shared five edges are exactly the central row/column AP tree.
- Gaussian elimination over every prime field through 59 gives horizontal/vertical/full ranks 5/5/7, including characteristics 2 and 3. Algebraically, the elimination uses unit pivots, so the theorem's arbitrary-characteristic claim is justified beyond this finite check.
- Independent exhaustive enumeration of the three actual core values E,A,C over every prime field through 59 reproduces every labelled seven-distinct-hourglass and nine-distinct-completion count in the producer certificate. In particular 13:12/0, 29:84/56, 47:276/0, 59:1044/232. No isomorphism quotient, integer lifts, or positivity inference is involved.
- Independent weighted root-cone and square-entry counts reproduce all fiber-product producer fields through 31. These are separate counts; root multiplicity is retained, and no distinct-entry condition is silently inferred from these weights.
- The independent checker has explicit exception guards and deliberately triggers/catches a failed guard. Normal and optimized execution outputs are byte-identical (`normal.json`, `optimized.json`).

## Proof audit

The seven-cell five-line definition, three-coordinate normal form, unique scalar completion, and gluing over arbitrary additive abelian groups are correct. Agreement at the five labels is equivalent to agreement at the three independent core labels. Characteristic 2 forbids injectivity; characteristic 3 does not invalidate the equations. The rational collision polynomials are correct, with their stated characteristic/order scope. The original draft incorrectly claimed two compatible injective hourglasses could have a noninjective union. That claim was withdrawn by the producer and the initial manual audit missed it. The corrected implication is valid: a private-private collision forces 2x=0 or 2z=0, which already repeats shared opposite corners. One-copy injectivity alone remains insufficient. The independent checker now exhaustively verifies the corrected two-copy implication for all offsets modulo every prime through 59.

The F13 center-zero witness has seven distinct square residues, but forced values 8 and 5 are nonsquares. It is not an integer counterexample. The integer repeated-entry witness is correctly scoped. The five triples are a maximal partial Steiner triple system: the unused pair graph is a six-cycle and has no possible additional triangle. This is not a seven-square arithmetic obstruction.

The historical 3+1+3 correction was read at its cited fragment and theorem. The relevant correction heading in MISTAKES.md is `(root audit of klein-S691)`; MISTAKE-319 is duplicated elsewhere in that large file, so the number alone is an ambiguous citation. The anchor properly refuses any transfer from cardinality or terminology alone. No grand Keller conclusion is needed or audited here.

Primary-source boundary: independently opened Cain, arXiv:1908.03236v2, Section 2, https://arxiv.org/html/1908.03236v2 on 2026-10-25. Its five equations, Theorem 2.1, Corollary 2.1, and displayed AP-only countercontrol match the anchor. No later parametrization claims are needed. No full paper is included in the audit packet.

## Reproduction

Run `python3 check.py` and `python3 -O check.py`. The checker reads the two producer JSON certificates at their original /tmp locations and compares results. It is intentionally an audit harness, not a portable replacement for the producer checks.

## Fiber-product sidecar audit

PASS for `/tmp/session58-hourglass-fiber-product/fiber_product.md`, including the C2 refinement. Forgetting added roots has no base locus: when all five shared roots vanish the quotient has nilpotent irrelevant ideal and empty Proj. Localizing at any shared root and taking degree-zero ratios proves the global projective fiber-product identity; no invalid general tensor-product/Proj rule is invoked. Monic relations give finite locally free ranks 4,4,16. Distinct linear valuations in the normalized parameter plane prove eight independent square classes in every odd characteristic, so the generic degrees and geometric integrality assertions follow. Relative non-etaleness is exactly where a new root vanishes. Characteristic 2 is correctly excluded.

The root-weight character sums, projective vertex subtraction, square-entry indicator counts, and covariance formula have correct normalization. Independent direct projective enumeration, choosing the unique representative whose first nonzero shared root is 1, gives p=3 counts (16,64,64,256), p=5 (36,88,88,320), p=7 (64,160,160,384). Its even/odd energies are respectively (256,0), (324,4), (448,64). These reproduce the corrected sidecar figures without deriving projective counts from cone division. The C2 involution exchanging d and e and the even-minus-odd energy identity are valid. Generic field-extension independence does not imply finite-field statistical independence.

The sidecar checker originally used removable assertions. With producer approval, the auditor replaced every assert by an explicit RuntimeError guard, added an intentional failure probe, and made the output path relative to the script. The corrected sidecar checker passes normal and optimized Python with byte-identical outputs. The independent `projective_check.py` also passes both modes with byte-identical outputs. External attribution of the nine-root surface is inherited from the prior session and explicitly labelled as such by the sidecar; this audit establishes the displayed coordinate statements independently, not historical novelty.

Final state: the initial two-copy-injectivity oversight is corrected, and both producer packets now pass the checks described here. Audit results apply to the recorded file hashes, not arbitrary future modifications.

## Second-wave addendum: mitres and the affine plane

PASS for MITRE_AFFINE_BRIDGE.md and check_geometry.py. Independently inspected Brouwer's primary CWI publication record and PDF, https://ir.cwi.nl/pub/6839 and https://ir.cwi.nl/pub/6839/6839D.pdf, Section 3, printed pages 6–7 (PDF indices 8–9), on 2026-10-25. The source explicitly names the five-block configuration a mitre and identifies deletion of two points and their incident lines from AG(2,3). No diagram or later existence theorem is needed. This classical attribution is correct.

Independent `mitre_check.py` does not import producer definitions. It constructs the affine Steiner lines from the unique third-point operation, enumerates every five-line subset to find all 36 mitres, and checks all 9! labelled arrays directly against the eight line predicates. It obtains 54 distinct same-center union patterns, 432 layouts, and exactly eight layouts over every pattern. This verifies the correspondence as actual labelled pattern equality, not merely equal cardinalities. The minimum Steiner host order 9 follows from the seven-point maximal-partial obstruction, the eight-point parity obstruction, and the explicit nine-point host.

The ordered additional-line defects were checked symbolically by integer coefficient vectors: AFH, BDI, BFG, CDH give +3x, −3x, −3z, +3z respectively. Hence all twelve lines force constant labels over Q (or any torsion-free additive group); characteristic 3 makes the added conditions automatic. Neither assertion supplies the extra square-membership conditions.

A second independent search normalizes every nonzero square center to E=1 by dividing all entries by that square. It finds no seven-distinct-square hourglass with both forced entries nonsquare at primes below 37, and exactly four normalized offset pairs at 37. The stated F37 witness is valid. This proves only the listed prime-field minimality, not a prime-power or integer minimum.

The standalone independent mitre checker contains an intentional exception probe and passes normal and optimized Python with byte-identical outputs. No producer source was edited in this second-wave audit.

## Root synthesis review

Reviewed Sections 1–4 of `/tmp/session58-root/first/SYNTHESIS.md`. The carrier, mitre, finite-field, fiber-product, and spectral summaries match the audited packets. Requested one precision edit in the historical comparison: identify the two hourglass-to-base ranks as four and the joint M-to-base rank as sixteen, rather than calling all present forgetful maps rank four. All are finite and hence proper. Contrasting these specified maps with the recovered historical fragment's escaping-sheet narrative is legitimate and explicitly does not transfer solution sets or endorse broader Keller assertions. Sections 5–6 are outside this auditor's assigned lane.
