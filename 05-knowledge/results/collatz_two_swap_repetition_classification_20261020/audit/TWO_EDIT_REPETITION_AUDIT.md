# Independent audit of the two-edit repetition classification

20 October 2026. **Verdict: PASS. No mathematical gap or required correction
found in Theorem R2 or its stated mechanical corollary.**

The elementary theorem is valid for every displayed repetition `v^d` with a
nonconstant binary base, including nonminimal and unbalanced bases. The
mechanical consequence additionally uses the already audited positive and
negative primitive one-swap theorems. This audit does not claim a primitive
two-edit exclusion, an arbitrary-source two-edit exclusion, or a result about
unrestricted Collatz trajectories.

## 1. Actual operations, cuts, and prime-power content

The rotation identity

    2 B(rho(w)) = 3^(w_0) B(w) + w_0 D

is exact. Since `H | D` and `gcd(H,6)=1`, it preserves the entire gcd with
`H`, not merely its prime support. The source after rotation remains a
`d`-fold repetition of a rotated base. Therefore the favorable cuts used in
both cases of the proof are legitimate without a balance or minimal-period
assumption.

There are at least four positions. Thus two different edges either overlap in
exactly one position or are disjoint; the two-position-circle ambiguity is
absent from the full source. When `q=4`, antipodal edges are still two disjoint
pairs. The base itself may have length two, where both cyclic edge labels
describe the same actual base transposition; this causes no problem for the
statement that the output is a square of a one-edited base.

In the overlap case, requiring the second move to be actual is essential.
It leaves the length-three chains stated in the proof and their reversals,
with local carry differences `±3` or `±5`. The outside prefix contributes
no residual difference because the local weight is unchanged; an outside
suffix multiplies the difference by its power of three. Position contributes
a power of two before recutting. Hence the exact contents are one and
`gcd(H,5)`, respectively. The lower bound `H>=7` suffices; there is no
incorrect claim that five is always invertible modulo `H`.

For disjoint edits, both are already legal in the source, and their
transpositions commute. More importantly, the earlier pair has the same
weight after it is edited, while the later edit does not touch the other
pair. This directly verifies the second perturbation's suffix exponent.
With the proof's interval convention, the exact difference is

    3^s (epsilon_0 3^m + epsilon_ell 2^ell).

The interval defining `m` includes the entire second unequal pair, so `m>=1`.
Consequently the binomial is never zero, including when its signs differ.

## 2. Bounds, equality, and nonminimal bases

An interval of length at most `(d-1)q0` can be extended inside the periodic
sequence to one of exactly that length; the latter has `(d-1)p0` ones.
This is a periodicity argument, not a balanced-word estimate.

For `d>=3`, the shorter-edge distance satisfies
`ell <= d q0/2 < (d-1)q0`. Thus the power of two is strictly below its endpoint
term, while the power of three is at most its endpoint term. Together with
nonvanishing this proves `0<|U|<H`, regardless of either sign.

For `d=2`, all inequalities in

    |U| <= 2^ell + 3^m <= 2^q0 + 3^p0 = H

are equalities precisely when `ell=q0`, `m=p0`, and the signs agree.
Periodicity forces the last two conditions when the edges are `q0` apart.
The edited word is literally `(v')^2` at the original cut, including edits
across copy boundaries. Conversely, two distinct non-antipodal edges cannot
retain all of `H`. The undo case is exactly equality of the two edge labels.

The theorem is representation-relative as stated. A nonminimal base may
give a smaller displayed `H` than the multiplier of its least-period
representation. Retaining the displayed `H` in a doubled edit does not
assert that the complete clock divides the carry. The identity
`B((v')^2)/D = B(v')/D0` preserves the remaining integrality question.
There is therefore no contradiction between different repetition
representations of the same source.

## 3. Mechanical corollary and small boundaries

The primitive balanced-base freeness used here has a self-contained
elementary proof. For `gcd(p0,q0)=1`, increasing the numerator of the
mechanical intercept from zero to one changes exactly one internal adjacent
`01` to `10`. These two words are rotations because rotation by
`p0^(-1) mod q0` produces the intercept-one word. A common clock divisor of
their carries must divide their signed `2,3`-smooth difference, which is
coprime to the clock. Thus the base carry has content one.

The positive unit-clock factorization is valid: reduction modulo three
makes `q0` even, and the coprime odd factors
`2^(q0/2)-1`, `2^(q0/2)+1` are powers of three differing by two.
They must be one and three. For the negative clock, `0<p0<q0` forces
`q0>=3`; reduction modulo eight makes `p0` even, and the two neighboring
even factors have gcd two and must be two and four. The only nonconstant
unit-clock shapes are therefore `(p0,q0)=(1,2),(2,3)`.

For the doubled-edit branch, the positive and negative one-swap dependency
statements cover primitive bases of the required signs. The positive
theorem's formal `q>=4` boundary leaves only the directly handled `p0=1`
cases and the small unit clocks; no small shape is omitted. At the two
surviving shapes every actual base edit is a rotation, so the full output
also lies on the original necklace. Exact necklace distance two is therefore
excluded for every nonprimitive mechanical shape, not merely the finite
ones tested below.

The proof correctly distinguishes the elementary repetition theorem from
the mechanical corollary's inherited analytic dependency. Only the
elementary freeness assertion of THM-4716
(`THM-4716-balanced-collatz-cycles-are-three-smooth-ladders-and-are-free.md`)
is needed, and it is proved directly here; unrelated claims in that file
are not dependencies.

## 4. Independent implementation and exact receipts

The new `audit_boundary_transport.py` imports neither existing verifier.
It represents words as integer masks, generates each actual ordered
two-swap path directly, and computes carry by a right-to-left positional
sum. It adds cases outside both prior universes:

- all nonconstant bases of lengths 2 through 7, with every `d=2,...,7`;
- all 510 nonconstant nine-bit bases with `d=2`;
- every primitive mechanical shape with length at most 200, for the
  intercept, rotation, smooth-difference, freeness, and unit-clock checks.

Normal and optimized runs are byte-identical and report:

- 3,013,615 explicit runtime checks;
- 1,950 sources and 437,288 actual ordered paths;
- 27,584 nonminimal-base paths;
- 28,252 paths crossing the original cut and all 16 length-four paths;
- 220,538 distinct-output rotation/cocycle checks;
- 12,618 overlap paths of each local factor type, three and five;
- 31,488 full-`H` retentions, including 6,144 mirrored edits;
- 12,231 primitive mechanical shapes.

The pre-existing parent checker was preserved unchanged. Its normal and
optimized reruns both reproduce the original receipt: 240 bases, 720 sources,
81,592 paths, 6,912 undo retentions, and 1,536 mirrored retentions.

Both hostile controls were independently regenerated. The `0011` square
with ordered edges `(3,2)` gives `H=25`, `B=680`, `D=175`, and exact
`H`-content five. The `11101111000` square with edges `(3,14)` gives
`H=4235`, `B=10007305`, `D=-588665`, and fixed point `-17`.
The latter has exact necklace distance two from its nonmechanical source.
The independent program additionally follows all 22 shortcut steps, checks
each prescribed parity, and returns to `-17`.

These computations are finite checks of the proof and its boundaries, not a
replacement for its all-length argument.

## 5. Reproduction and frozen inputs

From this audit directory:

    python3 audit_boundary_transport.py > boundary-receipt.replay.json
    python3 -O audit_boundary_transport.py > boundary-receipt.optimized.replay.json
    cmp boundary-receipt.json boundary-receipt.replay.json
    cmp boundary-receipt.json boundary-receipt.optimized.replay.json
    python3 check_periodic_two_edits.py > receipt.reviewer-replay.json
    python3 -O check_periodic_two_edits.py > receipt.reviewer-optimized-replay.json
    cmp receipt.json receipt.reviewer-replay.json
    cmp receipt.json receipt.reviewer-optimized-replay.json

Audited production verifier SHA-256:

    a03cbad0d25b72ee2f0132958d63801846ac553371f54cf917b3c91e3a28e10b

Both production receipt files have SHA-256:

    00b7b45a78c2612efcefdce1574df8d2eed871d609c6d7689d710f5de889d0f0

The theorem text as reviewed, including the finalized production counts,
has SHA-256:

    fd7419b3cfa199a161631736365e68cb9c5ee8ad2ed24e3288c12c14e2c7b1bf

Independent boundary-audit verifier SHA-256:

    d91849913eabeeeb8846e4a3fc75595a68e8b3b3f0b7d3778898dc90613e47bd

Both boundary receipt files have SHA-256:

    0b853d64a8eaa5fb93c5110e8fec397dbe750e4a6045923feecfddb4d0f1354a

The preserved parent checker has SHA-256:

    88ed2811b67d113db5acb4b9c8af3c0a582c80b646cb511c1e97469f1eebe22a

Its preserved normal and optimized receipts have SHA-256:

    a4cd732c97aa7e86eaa629bcd3e6becdc532ffe409e420db8c4933e04ccb1a0d

No publication or modification of a previously frozen packet was performed
in this audit.
