# Review of the reconstructed graft theorem

**PASS.** This review concerns the reconstructed argument and fresh calculations,
not the unavailable original Session225 artifacts.

## Proof review

The general four-row closure proof was independently reviewed with particular
attention to:

- ordered path-cover bins and the absence of a factorial quotient;
- empty ordinary bins, singleton sources and singleton gadgets;
- the unique possible B-to-A transition into the old root;
- the two possible bins containing that transition in a rooted two-path cover;
- the replacement fourth row when deleting a singleton gadget's root;
- strongness under the stated reachability hypothesis;
- separation from one fixed universal O/Z/J word, rather than individual
  matching counts or adaptive constructions.

The three-bin and two-bin subtractions remove precisely the forbidden empty
runs. No defect was found in these cases or in the stated scope.

## Fresh independent arithmetic and finite enumeration

A separately written standard-library program enumerates path permutations
and vertex-to-bin assignments. It does not import verify.py or the original
source bank. It freshly recovers

    (H,Q,s,e,alpha,delta,T,q,b) = (5,32,1,1,5,7,17,1,8),

the four physical basis vectors, det(K)=-8, and the basis determinant -4.
Leibniz determinant evaluation and integer matrix multiplication confirm

    K^4 - 27K^3 - 37K^2 + 31K - 8I = 0,

and H=1,37,1113,31349,886465 from the source-rooted-edge seed.

The independent script and receipt are included as independent.py and
INDEPENDENT.txt. These are fresh reproducible checks, not recovered old logs.

## Scope of the larger fresh replay

verify.py separately generates all 75 labeled rooted sources and 75 labeled
rooted gadgets of orders 1..4 with root at label 0. Its 5625 graft-profile
checks, 225 generator checks and 1161 applicable strongness controls passed.
It also compares subset DP with direct permutation counts, and counts literal
iterates at orders 6,10,14. Normal and optimized receipts are identical.

The independent reviewer did not rerun that entire 5625-pair suite. The proof
review and the independent small enumeration are separate from that larger
fresh replay; this distinction is intentional.

## Exclusions

No literature priority or smallest-gadget assertion was certified. No later
inverse refinement, all-target construction, finite-coverage extension, or
eventual cutoff is part of this recovery. Chan-Pak completeness remains OPEN.
