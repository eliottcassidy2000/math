# The binary-character gap is not individual doubling

**Update, 2026-10-29:** [Phase2](phase2/README.md) now proves the exact
2-primary saturation, index8192 and quotient(C2)^13 for both named banks.
This first-stage note retains its original bounded scope; its open
2-primary-saturation questions are resolved by that later certificate.

**Status: FINITE-EXACT / VERIFIED.** Bounded source-sensitive follow-up to the
Session97 rank certificates. No rank upper bound, no 2-saturation claim, and no
claim that a finite prime bank detects every rational square class.

## Results

For the 30 supplied points on each of ICARM curves **3567 and 3568**:

1. **None of the 60 individual points is a rational double.** Every individual
   point has a nonzero verified local F₂ character. Independently, each exact
   rational halving quartic is irreducible of degree four over Q.
2. The retained 17-row character bank has rank 17, hence a 13-dimensional
   binary kernel. Exhaustive enumeration of its 8192 vectors gives:

   | Curve | Minimum nonzero kernel weight | Number attaining it |
   |---|---:|---:|
   | 3567 | 6 | 4 |
   | 3568 | 5 | 1 |

   Consequently **no signed sum of between 1 and 5 distinct supplied points on
   3567**, or **between 1 and 4 on 3568**, can be a rational double. This follows
   because signs disappear under a homomorphism to F₂. More generally these
   bounds concern the support of odd coefficients, not arbitrary integer
   coefficient size. A combination with all coefficients even is trivially a
   double and is not covered by a nonzero-support assertion.
3. Minimum-support binary-bank kernel vectors are, in **zero-based source-list
   indices**:

   - 3567: `{1,14,16,27,28,29}`, `{6,7,18,24,27,28}`,
     `{4,8,9,13,26,28}`, `{3,6,8,12,16,28}`.
   - 3568: `{19,20,21,23,28}`.

   The finite-bank kernel alone gives **locally invisible candidates only**.
   In the optional bounded continuation, exactly the first listed support for
   each curve was tested with all signs positive. **Both are genuine rational
   doubles**, with the exact native-model halves stored in
   `joint_candidate_3567.json` and `joint_candidate_3568.json`. Independent
   native-model arithmetic verifies each identity. The other three minimum
   supports on 3567 remain untested. There was no 8191-combination
   elliptic-curve arithmetic hunt.

The prior independently verified F₃ certificates establish 30 independent
points for both curves. These new results do not weaken that conclusion. They
rule out a tempting explanation of the rank-17 binary bank (individual doubled
points) and prove **genuinely joint divisibility** instead. Given the already
certified independence of each 30-point list, the displayed half is outside
its integer span: otherwise doubling would contradict an odd coefficient in
an independent relation. Thus each list generates a subgroup with a genuine
order-two coset in its saturation. Adjoining the displayed half enlarges that
subgroup by exactly index two; replacing any one supported point by the half
gives an independent 30-generator bank for the enlarged subgroup. This is an
exact consequence of the rational identity and independence, **not** an
inference from rank 17. It does not establish full 2-saturation or prove that
the entire dimension-13 finite-bank kernel lifts globally.

## Exact proof mechanism

For the source model

`y² + a1 xy + a3 y = x³ + a2 x² + a4 x + a6`,

put `Y=y+(a1 x+a3)/2` and `g(x)=x³+A x²+B x+C`, with
`A=a2+a1²/4`, `B=a4+a1 a3/2`, and `C=a6+a3²/4`.
The tangent slope at `(X,Y)` is `g'(X)/(2Y)`, so for a finite non-2-torsion
half of a point with x-coordinate u,

`g'(X)² - 4g(X)(A+2X+u)=0`.

Expansion gives the monic quartic

`X⁴−4uX³+(−2B−4Au)X²+(−8C−4Bu)X+(B²−4AC−4Cu)`.

A rational half necessarily gives a rational root. All 60 quartics have no
linear factor over Q. The factors are stored exactly and the verifier checks
both reconstruction and irreducibility. The no-halves conclusion is also
certified without relying on polynomial factorization: at a stored odd good
prime p and cubic root r, the nontrivial quadratic character of `x(P)−r` is a
homomorphism `E(F_p) → F₂`, hence annihilates doubles. We verify the good
reduction, root, all nonzero evaluation values, and characters by an explicitly
constructed finite-field square set. None of these rows uses a special value
at r or the point at infinity. The character extension is the usual one
attached to the rational two-torsion point `(r,0)` on the completed-square
model. Reduction of a rational half is legitimate at good reduction, whether
or not that half's affine coordinates have p-integral denominators.

The quartic alone forgets the sign of Y. A positive control on
`Y²=X³−2`, with `Q=(3,5)` and
`2Q=(129/100,−383/1000)`, verifies the exact doubling formula and explicitly
rejects the opposite sign as a half of the same target. No individual target point has a rational candidate needing this final sign
test. The two joint targets do: their quartics factor as degrees 1+3, and both
rational-root candidates have square rational g-values. Their signs are
retained and verified directly by `check_joint_halves.py`.

## Source preservation and reproducibility

`halving_certificate.json` embeds precisely the two relevant source curves,
point lists in their original order, and the 17 retained character rows per
curve. It records SHA256 digests of:

- `../ranks/phase2/highrank_inputs.json`;
- `../ranks/phase2/all_highrank_character_certificates.json`;
- `../ranks/phase2/odd_primary_certificates.json`.

No full raw leaderboard database was read. Embedded source indices and model
coordinates are unchanged. `binary_kernel_certificate.json` stores full
weight enumerators, a kernel basis, and every minimum-support vector.

```
python check_halving.py
python -O check_halving.py
python check_kernel.py
python -O check_kernel.py
python check_joint_halves.py
python -O check_joint_halves.py
```

The halving checker needs **SymPy 1.14.0** and reports 66,628 checks, rejecting
four hostiles (point, character, factor, source index). The independent kernel
checker uses only the standard library and reports 921 checks, rejecting
three hostiles (basis, minimum weight, support). Normal and `-O` outputs match
byte for byte. The independent standard-library native group-law checker
reports 203 checks and rejects six hostiles (half, index, sign on each curve).
It confirms the joint sums and halves without using SymPy or the producer's
completed-square addition formula. The kernel checker uses a different elimination direction and
mask-to-sum enumeration from the producer. The local-character check is
independent of SymPy and uses square sets rather than Euler exponentiation.

`produce_halving.py` and `inspect_kernel.py` are bounded regenerators, not
additional mathematical dependencies of the checkers. The former needs the
adjacent compact source packet; the latter needs the embedded certificate.
`factorization.log` is a nonessential 60-line audit transcript.

This research made use of the Elliptic Curve Rank Leaderboard, maintained by
the Institute for Computer-Aided Reasoning in Mathematics (ICARM) under NSF
Grant DMS 2425401. Source records:
https://elliptic-rank.icarm.cloud/curve/3567 and
https://elliptic-rank.icarm.cloud/curve/3568 .

## Exact limit / stopping condition

The bounded individual-halving investigation, finite-bank kernel census, and
two optional minimum-support joint-halving tests are complete. Full
2-saturation, lifting of all 13 kernel directions, the three untested minimum
supports on 3567, and completeness of the finite character bank remain
**OPEN in this packet**. The two codes cannot be identical in the retained
point order or permutation-equivalent: their minimum weights differ.
