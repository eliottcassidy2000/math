# Independent audit: complete paired-charge dyadic compiler

23 October 2026. **PASS within the stated electrical and conditional-realization scope.**
The producer packet is `session45-dyadic-constraints`. This audit neither
imports its mathematical checker nor imports the original filtration checker.
It reads the pinned original edge export and the saved literal lift vectors,
then independently proves that those vectors have the required completeness.
All output files are separate from the publisher checkout; no Git operation
was performed and no publisher file was modified.

## Main finding

The original critical 2-primary group is independently certified as

    K_(2) = (Z/2)^3358 + (Z/4)^4 + Z/256.

Consequently its doubling image is **entirely dyadic**:

    Q_1 = 2 K_(2) = (Z/2)^4 + Z/128,  |Q_1| = 2048.

The complete owner-marked coordinate map may be taken to be

    (Y0, Y1, Y2, Y3 mod 2; X mod 128).

All 3,528 original coordinates and all 1,695 occupied fibres match the
producer table. The fibre-size distribution is exactly
`1:642, 2:545, 3:310, 4:139, 5:44, 6:15`.

## Independent filtration, not a trusted factor assertion

1. Rebuild the full root-reduced binary Laplacian directly from all 438,216
   original edge rows, grounding vertex 3527. Connectivity is checked. A
   high-pivot bit elimination, distinct from the old low-pivot routine,
   obtains binary rank 164 and a 3,363-dimensional literal nullspace.
2. For its binary basis B, compute the full half form `B^T L B / 2 mod 2`.
   Every actual conductance w is split as `2 floor(w/2) + (w mod 2)`.
   The even contribution is assembled from rank-one cut vectors over F2;
   the odd contribution is accumulated over the integers and only then
   divided by two. This is a complete bilinear form, not a diagonal test.
   Its independently computed rank is 3358 and radical dimension is five.
3. Multiply each saved Y vector by the literal original L. All five are
   harmonic mod 4. Their parity vectors are independent and are explicitly
   located in the independently computed half-form radical. Hence they
   span the entire next stage; their provenance is not taken on faith.
4. Direct integer quarter energies produce rows `[8,12,14,31,8]`, rank four.
   The nonzero parity class `X = Y0+Y4 mod 2` lies in the one-dimensional
   radical. Literal multiplication verifies `LX = 0 mod 256`, and

       X^T L X = 6,819,178,752,
       (X^T L X)/256 = 26,637,417, which is odd.

   Therefore that unique surviving primitive mode persists through levels
   three to seven and is killed by the level-eight form. The valuation
   multiplicities are independently `[3358,4,0,0,0,0,0,1]`.
5. The saved Z is independently checked harmonic mod 8, and every coordinate
   verifies `X = Z+2Y0+2Y3 mod 4`.

The source filtration theorem was also proof-checked. For symmetric
nonsingular A, let `V_j` be the parities of vectors x with `Ax=0 mod 2^j`.
The divided form `x^T A y / 2^j mod 2` is independent of the chosen lifts:
a change by 2z contributes an even integer because Ay is divisible by 2^j.
If `Ax=2^j b`, lifting the same parity to the next level is equivalent to
solving `Az=-2^(j-1)b mod 2^j`. The finite perfect pairing says this is
possible exactly when `w^T b` is even for every harmonic w modulo 2^j,
which is exactly the radical condition. Smith reduction then identifies
`dim V_j` with the number of invariant factors divisible by `2^j`.
Thus the independent ranks and literal lifts above certify the complete
filtration and do not merely verify a selected bank.

## Harmonic completeness and the exact tower

For each k from two through eight, the audit enumerates the actual full
3,528-coordinate vectors generated modulo `2^(k-1)` by

    2^(k-2)Y0, ..., 2^(k-2)Y3, X.

The resulting orders are respectively

    32, 64, 128, 256, 512, 1024, 2048.

These equal the full reduced harmonic-image orders obtained independently
from the verified filtration. Indeed, the reduction kernel of grounded
`H_(2^k)` consists of `2^(k-1)` times the binary kernel and has order
`2^3363`. Subtracting its exponent from the full harmonic order gives
`k+3`. No extra characters are silently omitted.

Therefore, for an integral symmetric zero-row-sum F, preservation by `L+2F`
of the complete original harmonic bank is automatic at k=1, while at k>=2
it is exactly

    FY_i = 0 mod 2 (i=0,1,2,3),
    FX = 0 mod 2^min(k-1,7).

For k>8 the general row-image theorem supplies saturation; this is not an
extrapolation from finite tests. At k=8 the first four generators have order
two and the final one has order 128, with no nontrivial relation among them.
This validates the displayed Q_1 coordinates as complete character-dual
coordinates, not merely the invariant-factor orders.

## General theorem and shared physical rows

The all-size statement is correct for a nonnegative integer b. A zero-sum
integer row d annihilates `ker(L mod 2^k)` exactly when d belongs to the row
image modulo `2^k`, equivalently `[d]` belongs to `2^k K`. The intersection
over k is precisely the odd-primary subgroup. Applying this to each row of
`2^b F` proves the criterion in `Q_b = 2^b K_(2)` in both directions. It
proves inclusion of original kernels, not equality.

For a rectangle sum, the full physical-row contribution is

    sum_t a_t (u_t[j] ell_b(v_t) + v_t[j] ell_b(u_t)).

Both terms are necessary when source and target roles overlap. Disjoint
supports justify the gcd reduction by Bezout; no such separated rule is
used for overlapping roles. Multi-source debts may cancel only after
aggregation at the same physical coordinate.

The binary-cut identity follows from zero row sums at each vertex. Hence
zero current modulo two is exactly even crossing degree at every row.
For higher bits, after t successful integer divisions, the accumulated
lower-bit current is `2^t c_t`. The next numerator `F x_t+c_t` must be even;
only then may it be divided. This proves the signed carry compiler, including
negative entries, without cancellation of a nonunit in a residue ring.

## Exact tests and compiler checks

The independent checker passes **530,454 checks** in both ordinary and
optimized Python, with byte-identical output. In addition to the original
matrix and tower checks, these include:

- 4,320 exhaustive finite row-image/full-harmonic comparisons for b=0..3,
  moduli 2,4,8, mixed odd/dyadic graph groups, signed amplitudes, and
  overlapping source/target supports.
- 5,184 physical-row binary-cut tests.
- 78,336 exact signed-carry tests.
- 840 deterministic multi-rectangle tower tests on original marked labels.
- The old 1238/2246 by 881/2762 full-bank rejection; the electrically safe
  913/2701 by 398/528 rectangle; the shared-target cancellation between
  source pairs 0/1 and 55/46, with each component separately rejected; and
  the higher-carry obstruction first appearing at bit one.

The separate CLI audit passes **17 checks** in both modes, with identical
output. It verifies the shared circuit and individual failures, the
higher-bit hostile, a doubled-overlap diagonal example, unchanged-resource
reservations, conflicting shared keys/pins, and a consistent key ledger.
It also verifies that every result explicitly disclaims native
resource authentication. Each CLI case was itself run in both modes.

Both producer-checker replays pass **60,205 checks**, agree byte-for-byte
with one another, and equal the original producer receipt. These replays
are supplementary; the independent result above is not inferred from them.

## Scope and non-blocking cautions

No mathematical correction is required for publication. The all-size
native wiring statement is valid **conditionally on the supplied,
authenticated finite menu and its full endpoint/complement records**.
It does not prove the existence of sufficiently many native options.

The menu CLI only checks electrical currents plus consistency of supplied
resource keys and shared values. In particular it does not authenticate
paths, owners, converse orbits, pin semantics, or the correspondence between
an amplitude and the number of actual supplied copies. Those obligations
must remain in the separate native certificate. Algebraic cancellation is
not by itself an injective native move.

Odd-primary preservation, equality of kernels or critical groups, a global
native bijection, and unconditional all-order native supply are not proved
or claimed. Retaining these limitations is necessary when this packet is
combined with the independent four-source native work.

## Reproduction

From a directory containing the producer packet and this audit:

    python3 -B session45-audit-dyadic/audit.py --repo /path/to/math --packet session45-dyadic-constraints --out session45-audit-dyadic/normal.json
    python3 -O -B session45-audit-dyadic/audit.py --repo /path/to/math --packet session45-dyadic-constraints --out session45-audit-dyadic/optimized.json
    cmp session45-audit-dyadic/normal.json session45-audit-dyadic/optimized.json
    python3 -B session45-audit-dyadic/compiler_cases.py --packet session45-dyadic-constraints --out session45-audit-dyadic/compiler-cases.normal.json
    python3 -O -B session45-audit-dyadic/compiler_cases.py --packet session45-dyadic-constraints --out session45-audit-dyadic/compiler-cases.optimized.json
    cmp session45-audit-dyadic/compiler-cases.normal.json session45-audit-dyadic/compiler-cases.optimized.json

`normal.json` and `optimized.json` contain every pinned original source
hash and every producer-packet file hash seen during the audit. The local
manifest and SHA256SUMS freeze this audit's deliverables.
