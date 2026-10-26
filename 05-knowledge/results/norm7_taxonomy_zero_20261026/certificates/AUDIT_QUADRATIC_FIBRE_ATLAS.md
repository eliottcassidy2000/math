# Independent audit of the quadratic-field resolved-fibre atlas

2026-10-26. **PASS.** No correction required for Theorem QA or its declared scope.

Reviewed `../arithmetic/QUADRATIC_FIBRE_ATLAS.md`, SHA256
`bee684647eb9476bdc90d920272106b35ab550e35315ce456a0e8870fba802fe`,
and `../arithmetic/check_quadratic_atlas.py`, SHA256
`4647c9759940aed132a9cb375310213f717ed0bf4fa17c4109091c0e07b24d28`.
The inherited `../../normseven_descent_moments_20261026/RESOLVED_FIBRES.md`
was reread; its SHA256 is
`bed3077279fc2ddb765d39147fe0bfca8d6321d005a0e1f5e134470efee72a32`.

## Proof audit

The inherited fibre is four conics over the actual quadratic node field,
regarded as Q-schemes by composition of structure morphisms. The audit
retains this carrier, rather than replacing it by a Weil restriction or
by the entire resolved surface. An F-point really forces an embedding of
the node field into F, so the lower degree bounds are valid.

The conic/quaternion map `[A:B:C] -> [delta C:A:B]` is invertible for the
stated nonzero delta and takes the equation to the pure quaternion norm
for `(-1,delta)`. The signed squarefree representative includes denominator
valuation parity. Odd-place, dyadic and real Hilbert symbols give exactly
the displayed ramification set; the dyadic formula retains the odd part
of delta even when delta is even or negative.

I independently consulted the primary textbook statements:
[Voight, Theorem5.5.3](https://link.springer.com/chapter/10.1007/978-3-030-56694-4_5),
[§12.4.12–16](https://link.springer.com/chapter/10.1007/978-3-030-56694-4_12),
and [Proposition14.6.7](https://link.springer.com/chapter/10.1007/978-3-030-56694-4_14).
Their characteristic and finite separable extension hypotheses hold here.
The last criterion uses *no ramified quaternion place splitting in the
quadratic field*, rather than merely the existence of one nonsplit place.
Together with the conic criterion it proves sufficiency as well as necessity.

The square-D branch is correctly separated. In the nonsquare branch,
degree2 is the minimum if the conic has a K-point. Otherwise every possible
degree is a multiple of2 and cannot equal2, while K(i) provides degree4.
If the conic is anisotropic then i cannot already lie in K, since
`[1:i:0]` would be a point. This validates the upper bound and excludes a
hidden degree collapse. The infinity point is permitted because the fibre
is projective; no affine norm equation is silently substituted.

The split-torus criterion is relative to K: a split torus contributes its
identity point, while a soluble conic is P1 and decomposes into G_m and
two points. For nonsquare D, the full Q-fibre has no Q-point even when
its conics split over K. Thus the proof's relative and absolute F1 claims
have the right field and no extra implication for the whole surface.

Hostile boundaries checked separately: delta3,D17 retains the dyadic
obstruction after the odd-prime test passes; delta−1,D2 retains the real
obstruction after all finite tests pass. D13 and D−7 test, respectively,
the ternary and dyadic embeddings. The formula for delta3 has no real-place
condition because its quaternion splits over R.

## Computation audit

Replaying the author's checker in normal and optimized Python gives
16,773 gates and byte-identical JSON of SHA256
`296406290b1003a43bc4feae0a6e927f63121dc967a2e3620338ddd526e9c5ec`.
The code's unit-coordinate charts exhaust all primitive local solutions.
At odd primes a primitive solution modulo p² has a unit derivative for
squarefree delta. At2 a primitive solution modulo8 has a derivative of
valuation1, so the stated strict Hensel inequality applies. Neither finite
absence nor a point count is used as an unsupported global theorem.

The separate `audit_quadratic_atlas.py` imports no producer implementation.
It searches every primitive triple modulo8,9,25,49 for every nonzero
signed squarefree delta in [−20,20], checks the local-symbol prediction,
and tests explicit high-order local roots for D17,D13,D−7. This is a second
enumeration route, not a re-proof of the cited local-global theorem.

Reproduce from this directory:

```sh
python3 -B ../arithmetic/check_quadratic_atlas.py > /tmp/atlas.replay.json
python3 -B -O ../arithmetic/check_quadratic_atlas.py > /tmp/atlas.optimized.json
cmp /tmp/atlas.replay.json /tmp/atlas.optimized.json
python3 -B audit_quadratic_atlas.py > quadratic_audit.normal.json
python3 -B -O audit_quadratic_atlas.py > quadratic_audit.optimized.json
cmp quadratic_audit.normal.json quadratic_audit.optimized.json
```
