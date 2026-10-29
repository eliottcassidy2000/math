# Exact 2-primary saturation of the two rank-17 binary-bank exceptions

**Status: FINITE-EXACT / VERIFIED.** This new packet extends, but does not edit,
the frozen parent halving packet. The prior packet's open boundaries describe
its own bounded stage; the present certificate resolves its 13 missing binary
directions for precisely these two supplied point banks.

## Result

Let L be the integer span of the 30 named source points on ICARM curve 3567
or 3568. We construct a rank-30 subgroup L' of E(Q) satisfying

- `L ⊂ L'`, with exact index `[L':L] = 8192 = 2^13`;
- `L'/L ≅ (Z/2Z)^13`;
- **L' is 2-saturated in E(Q): if 2Q belongs to L', then Q belongs to L'.**

Consequently L' is exactly the **2-primary saturation of L**, not merely an
index-two improvement. We do not claim it is the whole Mordell–Weil group, do
not prove rank at most 30, and do not certify saturation at odd primes.

| Curve | Replacements | Binary ranks encountered | Final maximum prime | Rootless good prime | Maximum final coordinate bits |
|---|---:|---|---:|---:|---:|
| 3567 | 13 | 17,18,...,30 | 223 | 23 | 603 |
| 3568 | 13 | 17,18,...,30 | 251 | 3 | 630 |

All replacements were tested as exact rational group identities, not inferred
from local character vanishing. The completed-square half coordinates are
stored with their signs. Each step replaces one supported basis vector P_i
by Q satisfying `2Q = sum(P_j : j in support)`. The old P_i is recovered as
`2Q − sum(P_j : j in support except i)`. Each replacement therefore enlarges
the free subgroup by exactly index two.

The smallest-support candidate at each subsequent step lifted successfully;
there was no exhaustive rational arithmetic over the binary kernel. Supports
had sizes:

- 3567: `6,6,6,6,7,7,8,8,8,8,9,11,11`;
- 3568: `5,6,6,6,6,6,8,8,8,8,9,9,10`.

The two processes terminated at full binary rank, strictly within the stated
20-replacement and 4096-coordinate-bit limits. Every quartic factorization
finished within its 60-second ceiling. These are computational budget facts,
not a claim that this strategy always terminates on arbitrary point banks.

## Why the binary certificate proves 2-saturation

For each retained odd prime p of good reduction and rational cubic root r,
use the two-torsion local character on E(F_p):

- `O` has character `+1` (binary 0);
- `(r,0)` on the completed-square model has character `g'(r)` modulo squares;
- other affine points have character `x−r` modulo squares.

These are homomorphisms and annihilate doubles. Their combined 30-column
matrix on the final basis has rank 30 over F₂. The independently checked
rootless good prime excludes rational two-torsion. It follows that the basis
is independent: any integer relation would have all coefficients even, and
absence of rational two-torsion allows division of the relation by two.
Repeating this contradicts a nonzero finite coefficient.

If `2Q = Σ a_i P'_i`, every character vanishes on the right, so every a_i is
even. Thus `2(Q−Σ(a_i/2)P'_i)=O`. The independently checked rootless good prime
excludes rational two-torsion, hence Q is in L'. This proves 2-saturation
without assuming that E(Q) has rank 30.

## Exact transport and the meaning of the original deficit

The certificate retains every half, support, replaced index, preceding local
bank, and updated rational basis row. With rows as the convention, the final
basis is `P' = M P`. We independently verify:

- `det(M)=2^-13`;
- every entry of M has denominator 1 or 2;
- the recorded inverse N is integral and `N M = I`.

Thus `2L' ⊂ L ⊂ L'` and the quotient has order `2^13`, proving the displayed
elementary abelian quotient rather than merely guessing it from the number
of steps. The inverse transports explicitly recover all original points.

There is now a precise retrospective interpretation of the original binary
bank: its rank-17, 13-dimensional kernel on L/2L equals the kernel of
`L/2L → E(Q)/2E(Q)`. Indeed, because L' is 2-saturated, the latter kernel is
`2L'/2L`, of dimension 13; it is contained in the finite character kernel,
which has that same dimension. **All thirteen missing directions represent
actual global joint divisibility for these inputs.** This conclusion required
the exact lifts and terminal full-rank certificate. It did not follow from
the original rank-17 bank alone.

The smallest examples illustrate why observations of individual points miss
joint information: no original P_i is a rational double, while selected sums
of five or six of them are. The obstruction becomes visible only in joint
parity coordinates; this is an exact group/code interface, not an asserted
identification with a named graphic, Golay, or Hamming code.

## Reproduction and audit

```
python verify_extension.py
python -O verify_extension.py
```

The verifier is **standard-library only** and uses native generalized
Weierstrass group arithmetic, not the producer's completed-square group law.
It checks every step, all retained local rows, source hash, final matrix,
inverse transport, determinant, good/rootless primes, and six hostile
mutations (half, row, and matrix for each curve). Local characters are checked
with explicit finite-field square sets rather than Euler exponentiation.
The source bank is `../halving_certificate.json`; no full raw database is
used. Final curves and source order remain exactly those of the frozen packet.

The `O` and two-torsion special cases are exercised by the actual data; they
are not silently dropped. A rational point with denominator divisible by p
in x reduces to O at good reduction; otherwise the affine character is used,
with the derivative value at the two-torsion point. Both normal and `-O`
runs have byte-identical output. Exact counts are emitted by the verifier.

`extend_saturation.py` is the bounded SymPy-1.14 producer. It preserves each
step and all stopping metadata in `extension_3567.json` and
`extension_3568.json`; its nonessential progress transcript is `progress.log`.
`inverse_transport_*.json` provides explicit integral recovery maps.

This research made use of the Elliptic Curve Rank Leaderboard, maintained by
the Institute for Computer-Aided Reasoning in Mathematics (ICARM) under NSF
Grant DMS 2425401. Source records: https://elliptic-rank.icarm.cloud/curve/3567
and https://elliptic-rank.icarm.cloud/curve/3568 .
