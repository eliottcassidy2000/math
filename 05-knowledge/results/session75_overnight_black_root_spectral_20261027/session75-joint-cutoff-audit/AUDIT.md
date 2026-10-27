# Independent audit: native joint depth/clock cutoffs

**Decision: PASS, in the finite scopes below. No mathematical correction to the
producer packet is required.** This is an independent-agent review of the
frozen producer manifest
`09f5e8372a8133109cec03447284c49881440b36521eb22d2ef2ea7c67cf0ef7`.
The producer directory contains 129 manifest-listed files plus the manifest
itself. Every listed byte hash, exact directory membership, and the four
repository source extractions at commit
`96e35642110104cc547516a63f37eebcf930f707` were checked.

## Numerical acceptance and exact coverage

No additional native scales were tested. The universe is precisely the 19
scales in the two frozen preregistrations, not all integers between 4096 and
262144.

| Check | Independent coverage |
|---|---:|
| Uncached source reconstruction, all source fields | 1,000,908 records, all 19 scales |
| Global depth/clock breakpoint Cartesian completeness | 245,826 recorded rectangles, all 19 scales |
| Fresh pair-birth recomputation of every rectangle field | 42,425 rectangles at N=4096, 18275, 65536 |
| Fresh direct binning | 43 cutoffs, including each scale's maximizer and larger-field controls |
| N4096 boundary clipping | All 16 removed records, all 951 boundary pairs, all ordered port aggregates |
| N4096 single clock arrival | Every same-terminal interaction, including zero beta |
| Generic scalar controls | The pre-existing powers 4 through 4096 only; NON-NATIVE |

The fresh pair audits comprise:

| N | Rectangles | Same-terminal pairs | Nonzero pairs | Maximum D+ |
|---:|---:|---:|---:|---:|
| 4096 | 5,917 | 91,710 | 76,697 | 495 |
| 18275 | 13,348 | 1,705,728 | 1,411,131 | 4,132 |
| 65536 | 23,160 | 21,818,296 | 17,979,937 | 4,790 |

At these three scales the independent surface verifies **every** saved row's
arrived count, arrived norm, terminal energy, and debt, as well as the complete
maximizer list. The signed level-vector surface is hash-bound in the receipt.
At the other 16 scales this audit does **not** independently recompute every
numerical rectangle: it verifies all source records, the complete rectangle
universe, the producer algorithm, every saved summary maximum against the saved
surface, and fresh direct maximum/control ledgers. The producer's 19-scale
normal and optimized receipts are byte-equal and hash-bound; that is distinct
from independently executing every calculation.

The accepted producer finite result is that all 245,826 listed rectangles
satisfy `D+ <= S_N/2`; the largest listed ratio is exactly `4132/44381` at
N18275. Its maximizer and its entire N18275 surface were independently
recomputed. No finite result establishes an all-N theorem or a convergence
claim.

## Independent methods and algorithm review

`independent_audit.py` does not import any producer implementation. Source
charges are rebuilt by multiplying Eisenstein integers according to
`omega^(B^j(s)) (omega-omega^2)` or `omega^(B^j(s)) omega`, rather than copying
the producer's six coordinate rows. Every source uses a fresh uncached route;
primitive divisibility, each inverse braid, the exact affine clock identity,
cycle/repetition failure, first-entry terminal, final exterior port, integer
clock, depth, charge, sign and `(s-1)/2` address are retained and compared.
All four pinned source files also match `git show` at the specified commit.
The 627-state frontier's closure and descending ranks were checked.

The full source denominator is recomputed before clipping and includes every
moving label and zero-charge label. The source JSON's canonical hash is also
recomputed from the independent integer records. Hence a field comparison
cannot silently identify a boolean with an integer, omit a zero record, or
replace a nonarrived label by a smaller denominator.

The producer's sweep processes each global depth and every global integer
clock, never only selected endpoints. For every depth it starts a fresh set
of bins and adds each eligible source exactly once at its clock. Its update
is exactly `Q(x+v)-Q(x)=Q(v)+B(x,v)` and skips only zero-charge arithmetic, not
the label count. Terminals remain part of every bin key. Whole-level energies
are completed before taking positive parts. The final modulus is the first
power of two at least `2N`, so all distinct source addresses are separated.
The final energy must equal the arrived source norm.

For any integer cutoff, the selected set is constant between successive
global resource values. Below either minimum it is empty; beyond either
maximum it is constant. Therefore crossing all distinct depth values with
all distinct clock values, plus the zero-debt empty case, is a complete
classification of integer cutoff rectangles for that finite source family.
No route is silently omitted: an unresolved route fails the run instead of
reducing the denominator.

The independent surface does not reuse the producer's energy recurrence.
For each named same-terminal pair it uses the lowest set bit of `m xor n`
(equal to the lowest set bit of `m-n`) and puts the signed polarization at
the pair's coordinatewise resource maximum. Flat arrays are prefixed in two
one-dimensional passes. Single-source count/norm surfaces are independently
prefixed too. This checks every field, not just the peak or telescoping sum.

`supplementary_checks.py` additionally replays the producer's **N4096** full
sweep, full literal source reconstruction, detailed route/B-chain hash, and
all-rectangle pair routine. This deliberately limited producer replay is
reported separately from the independent implementation; the entire
19-scale `run_all.sh` was not rerun by this auditor.

## Proof review: pair allocation and clipping repair

For a fixed terminal-labelled refinement tree, expand each quadratic bin
energy. The same diagonal `Q(v_s)` appears at both adjacent levels, so it
cancels. A named pair contributes to the coarse bin but not the next fine
bin precisely when that level is its first differing low bit. Different
terminals never share a vector sum. Both pair members are selected by a
rectangle precisely when it contains their coordinatewise maximum depth
and clock. This proves the stated signed pair-birth identity without a
positivity assumption or a change of source law.

Write `F` for a selected family, `K` for its retained subfamily, and keep
removed labels as zero leaves in the **same fixed partition**. Let
`delta_k = Delta_k(F)-Delta_k(K)`. For every pair of real numbers `x,delta`,

`(x-delta)_+ - x_+ <= (-delta)_+`.

For `delta>=0` the left side is nonpositive. For `delta<0`, the positive-part
map can increase by at most `-delta`. Substituting `x=Delta_k(F)` and summing
proves `D+(K)-D+(F) <= sum_k (-delta_k)_+`. The pair expansion makes each
`delta_k` exactly the signed sum of pairs touching the removed family at that
level. Pairs with two removed labels occur once. No source diagonal is owed.
The fixed-tree/zero-leaf convention is necessary: arbitrary deletion followed
by compressing or renumbering levels is not the comparison being proved.

The authentic N4096 receipt is correct:

- Full original `S_N=9949` is unchanged.
- At fixed clock64, depth67 gives debt439; depth36 gives debt495.
- Exactly 16 removed labels have norm22, so the gain56 refutes payment by
  removed norm alone for this native clipping operation.
- Their exact boundary delta is
  `[-34,-19,-1,-11,13,-6,-8,2,2,0,-3,0,0]`.
- All 951 boundary pairs, including 808 nonzero pairs, match by source IDs,
  terminal, both ordered entry ports, beta, low-bit owner, and whether both
  endpoints were removed. The negative-delta invoice is82 and slack is26.

This is not an R3 counterexample, not a violation of the joint half-source
target, and not a calculation of paid forcing.

For the **NON-NATIVE** scalar balanced tree with `n=2^h>=4`, one charge is
`-1` and all other charges are `+1`. In bins of size `a`, the exceptional bin
has sum `a-2`, changing to `a-1` on zeroing that leaf. Thus
`E_F(a)=(n-4)a+4` and `E_K(a)=(n-2)a+1`. All adjacent drops are nonnegative.
Their debt difference is `2(n-1)`, while each boundary delta is `-a` and the
invoice is `n+n/2+...+2=2(n-1)`. The removed norm is1. This proves sharpness
and unboundedness of the invoice/self-norm ratio only in the stated generic
class. It supplies no native all-order counterexample.

## Hostiles and integer safety

At N4096 `(D,T)=(36,64)`, independent direct binning gives the correct debt495.
Premature per-parent positive parts give11227 and terminal erasure gives56715;
both are false manufactured violations against S9949. Losing the clock sign
gives209. Using odd `s` instead of `(s-1)/2` changes the full energy path but
accidentally leaves debt495: debt-only comparison is insufficient. All source
metadata mutations, a zero-label deletion, and a moving-source deletion are
rejected by comparison to independently rebuilt records.

Source13759 is the sole clock64-to66 arrival at depth36. Its three nonzero
interactions have `(source,q,beta)` equal to `(9,1,-1)`, `(2831,8,+1)`, and
`(10065,1,-1)`. The signed whole-level changes are `-2` at q1 and `+1` at q8;
the debt changes495 to494. All zero-beta interactions are checked as well.

Python arbitrary-precision integers are used throughout both independent
programs and the reviewed producer algorithms; explicit exceptions rather
than `assert` enforce all gates. At N262144, independent source1042431 has
depth116, clock194, terminal921, and exterior entry port9827. Its route reaches
28167930305, exceeding unsigned32-bit range. `LARGE_INTEGER_WITNESS.json`
retains all116 exact edges. Fixed-width32-bit arithmetic is not safe here.

## Preregistration provenance

Both registrations are hash-bound into the relevant receipts. Their bytes
are unchanged in the earlier `core17` and `core17_verified` snapshots. Local
mtime evidence places the original registration at07:20:08Z, before the
N4096 checkpoint log at07:21:08Z, and the phase2 registration at07:26:58Z,
before the first N11001 receipt at07:27:38Z. The source-pinned inherited README
already identifies N18275 and the N11002 control, consistent with the stated
selection rationale.

These are corroborating local-file observations, **not a trusted timestamp or
cryptographic proof of write-before-execution chronology**. The JSON flags
alone cannot prove chronology, and file mtimes can be altered. The finite
identities and completeness checks do not rely on a stronger provenance
claim. `ORDERING_EVIDENCE.json` records the precise evidence and its limits.

## Reproduction and deliverable binding

Run `sh run_audit.sh PRODUCER_PACKET REPOSITORY_CHECKOUT`. This runs each
independent/replay checker under normal Python and `python -O`, requiring
byte-equal receipts. It reads but does not modify the producer packet or the
repository. Only standard-library Python and local `git show` are needed.
`REVIEWED_PRODUCER_BYTES.json` binds every reviewed producer member;
`MANIFEST.json` binds this audit's scripts, report and completed receipts.

The audit made no shared repository edits, remote actions, numerical-range
expansions, source substitutions, or global-convergence claims. Both open
native targets—depth-only R3 and its stronger joint-cutoff extension—remain
**OPEN**.
