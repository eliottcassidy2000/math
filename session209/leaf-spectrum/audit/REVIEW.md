# Independent review: joint paired-black leaf spectrum

**Verdict: APPROVED in the stated scope.** The exact joint generating
function, coefficient formulas, complete leaf H-spectrum, minimum leaf
orders, fixed-order maxima, labelled weighting, exponential moments and
Gaussian limits are valid consequences of the inherited leaf classification.
No unresolved mathematical defect remains in the reviewed packet.

## Proof review

- The Session14 rigid-word and center reduction is not independently
  reclassified here. Session15's complete deficit-two certificate proves
  the unique strong center is S4. Its current publication and subsequent
  metagraph classification continue to carry that result; targeted current
  searches found no superseding retraction. The older open-center prose is
  historical, not the governing status.
- The palindrome term correctly doubles both triangle count and order.
  Mixed half-words are not identified with their reversals. H=1 gives no
  leaf. The two power families and their minimum LEAF orders follow in
  both directions, without claiming global realizing minima.
- The maximum comparison is correct on every residue class, including
  orders4,6,8. Rare mixed leaves do maximize H at orders0 and4 modulo6.
- The labelled generating function correctly counts two ordinary classes
  over each NS quotient leaf and one over each mixed leaf. Forgetting this
  factor already fails at order4:16 and24 labelled lifts, respectively.
- The unique dominant simple pole is justified for positive weights;
  the reflection/mixed poles have a strict larger radius. Analytic
  continuation and root separation support differentiated coefficient
  asymptotics, the exact moment prefactor and the characteristic-function
  proof. Both reported mean/variance constants are correct.
- Review requested one clarification, now applied: bounded difference
  between log H and K log3 alone does not prove an O(1) variance remainder.
  The exponentially small mixed probability gives Cov(K,1_mixed)
  =O(n eta^n), which supplies the claimed remainder.

## Independent controls and replay

`independent_leaf_audit.py` imports no producer code and enumerates literal
singleton/C3 words through order30, selecting reversal representatives by
tuple comparison. All176 joint coefficient cells agree with the closed
formulas. It checks the full inherited rational numerator/denominator and
integrality of labelled weights. It also enumerates every vertex permutation
to count H and automorphisms on all15 constructed leaf classes through
order8; every count agrees with its stated product formula.

Run `python3 independent_leaf_audit.py` and
`python3 -O independent_leaf_audit.py` in this audit directory. Both outputs
are byte-identical to `independent_leaf_audit.json`.

The complete producer `verify.py` was separately replayed in an isolated
copy, normally and with `-O`. Both receipts and both stdout streams agree;
the generated receipt matches the published RESULT.json byte for byte.
That replay includes the literal path/anti-automorphism-orbit census through
order5, which decides leaf status without assuming the SCC classification.

Reviewed SHA256 values:

- PROOF.md: d7f466d8347c0ef94f381254531aa313c264caafe7198f3007460da6191fec7c
- verify.py: 7250ae757581601d47404f43f7afcc030fa50e4065afef48bf2649a84f87ba11
- RESULT.json: d14fdd4759f6ffa70cd34aa7734d776402ef7c56cf446812ae3098ef0160c446

## Retained boundaries

This is the degree-one stratum of the reflection-paired black multigraph,
not simple-support degree, native-labelled degree or the entire Hamiltonian
spectrum. Leaf vertices are not distinct bridges: two order4 leaves share
one edge. Uniform quotient leaves and uniformly labelled preimages have
different laws. None of the results reconstructs wiring, proves all bridges
pendant, or establishes an eventual Hamilton-count threshold.
