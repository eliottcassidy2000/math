# Proof and certificate audit

**PASS within E1-E4's stated scopes.** Same-session proof review with
independent algorithmic paths; no separate-agent review is claimed.

## Adversarial proof audit

* **Types/division:** Delta/2 and (C-C^t)/2 are formed over integers before
  reducing modulo2. D is the actual mixed quadratic of delta; h is the
  actual shared-label cubic. Repeated names across a and c are permitted;
  repeated names in disjoint old path supports are forbidden.
* **Two proofs:** the even-order proof uses the support cutoff plus the
  surviving singleton. The all-order proof complements the direction at
  the last endpoint, retaining one odd contribution for every nonempty
  even subset of the deletion card. It uses the full endpoint-path sum
  rather than the truncated sum and independently explains the constant.
* **Failure boundary:** at n=2 that subset count is0, so the row identity
  fails and is not used. The actual3-versus5 query closes n=2 separately.
  Orders0/1 have no reversed arc; no separation is claimed there.
* **Actual consequence:** E2 is derived by Boolean finite differences,
  retaining every lower term until cancellation. At even order the final
  two coefficients have opposite parity, so exactly one final query is
  included. At twelve this gives exactly12 distinct masks, each with
  c=e_v. The conclusion is an actual count difference2 modulo4, not just
  a nonzero formal coordinate or a nonzero norm.
* **Endpoint balance:** needed only to make the twelve-query sum2 at every
  anchor. E1 and E3 do not assume it. If b_v=1 the named singleton query
  already distinguishes. Primality and strongness are likewise unnecessary
  for the new obstruction.
* **Ambient transfer:** T1 is applied to the port-extended modules, using
  identical root incidences outside. The odd factor H(Q) preserves residue2.
  Neither exact core cancellation nor cubic truncation of an internal
  ambient reversal is invoked. The quotient may be arbitrarily large/strong.
* **Names/symmetry:** no relabeling of queries is allowed. An isomorphism
  between a core and its converse does not imply fixed-query equality.
  No inference to arbitrary unrelated pairs or multiple module changes
  is made.

## Exact computation and hostile controls

`census.py` reads the SHA256-pinned complete bank of198947 classes in
orders2..9 and compiles `matrix.cpp`. Its bit-packed simple-path endpoints
compute D, the full even-support endpoint sum, and the even-order short
support sum separately. Every even-order short/full comparison passes;
the row identity passes for every n>=3, and the inherited odd-order global
cover identity passes independently. There are no inherited structural
filters besides the stated complete bank and order range. This is a
regression certificate; the all-order proof does not use census completeness.

`certificate.py` uses integer forward endpoint DP and actual ordered
two-cover convolution to recover D and every deletion current. All1099
labelled tournaments through5 pass, with literal word enumeration of
sidecars through4 and708 actual extended-graph Hamilton-word counts.
The separately compiled predecessor DP evaluates1268 actual query counts,
including both prime twelve-core controls and their internal reversals in
a strong cyclic quotient; port-extended order is at most16 and16! fits
uint64. Thirty selected counts are independently recomputed by forward DP.

The retained full integer cubics of the two prime twelve controls each
have1753 labelled coordinates. Their source polynomials are hash-pinned;
both mirrored shared-cubic families match exact deletion currents, all
mixed quadratics match D, and312 newly evaluated actual queries match the
whole retained polynomials. Remaining named cubic coordinates are retained
unchanged, including all zeros. The frozen earlier full responses are
not overwritten by the new sparse certificate.

Hostiles are mandatory and pass: C3[TT2,TT2,TT2] has D=0 but nonzero shared
cubics; Paley11 plus source has zero shared cubics but nonzero D; the
endpoint-balanced H249, H300453 and H84599 examples remain separated.
A same-module comparison in a strong quotient gives identical counts,
and the singleton core is an equality control. The n=2 exception is
tested literally.

Normal and optimized outputs are byte-identical for both engines and the
full cubic tables. Gates use exceptions, not removable Python asserts.
The validator pins inputs, sources, outputs, theorem dependencies and the
bank; optional replay writes only temporary outputs.

## Correction lineage

No prior proved obstruction is retracted. The old cover-only and
shared-cubic-only universal guesses remain false, with their original
hostiles. The new theorem proves these two failure mechanisms cannot
occur together. Previous UNKNOWN SAT runs remain UNKNOWN historical
records; they are not retrospectively described as UNSAT certificates.
The previous open status is superseded by E1-E4, not by finite evidence.
