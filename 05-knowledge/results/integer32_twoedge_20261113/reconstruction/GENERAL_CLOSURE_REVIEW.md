# Independent review of projective-to-affine observer closure

**Verdict: PASS for G1, G2, and the exact checker below.** This is an
independent review; no producer file was edited or promoted. The proof's
RESERVED status remains a root integration decision. Finite replays have
only their declared universes.

| Reviewed packet path | SHA256 |
|---|---|
| `census/GENERAL_CLOSURE.md` | `a7f664adbadd16efc75f9d4ade404a9358971b1fd8d2a6eea56c6347f0aa31a1` |
| `census/check_general_closure.py` | `d0b02864e5c99bf9856ecfb4a89b2c770fdf9307b4a87a58aac74052277c58d2` |
| `census/general.normal.json` | `8efa40c340ec2d12a1c15cc389b29ecb26c0be5fb48fc848a20622b64d0755bf` |
| `census/general.optimized.json` | `8efa40c340ec2d12a1c15cc389b29ecb26c0be5fb48fc848a20622b64d0755bf` |

## G1: exact first partition and two-round recovery

Linearity makes the initial projective label of Bx depend only on the
initial label of x; invertibility ensures a nonzero line stays a nonzero
line. Hence B cannot split an initial cell in round one. On a line different
from F_q v, equality of the observations of cx+v and ex+v forces
`cx+v=lambda(ex+v)`. Linear independence of x and v forces lambda=1 and c=e.
Thus every point outside the distinguished line is already a singleton.
On that line, zero is initially separate and translation singles out −v;
all remaining q−2 points have the same first signature. These statements
prove the exact partition, not just its number of blocks.

When q>=4 the surviving q−2 points are a genuine nonsingleton cell. The
line-moving condition puts their distinct B images outside that line, where
they already have distinct singleton labels. One further response step
separates all of them. The first round is strict because d>=2 supplies a
nonzero scalar orbit outside the line, and the second is strict because
q−2>=2. The claimed exact two-round depth follows.

The q=2 boundary is correctly separated: every nonzero scalar orbit already
has size one, so the initial partition is discrete. At q=3 the first residual
has size one and the first round is sufficient. The formula q^d−q+3 must not
count an empty residual at q=2; the text explicitly avoids doing so.

The bounded-word interpretation is correct by induction for the declared
simultaneous refinement retaining the current label. It supplies the coarsest
stable observation refinement for those deterministic operations. It does
not supply a universal two-round theorem for other action systems.

## G2: preserving the line may stall or take arbitrarily many rounds

For B equal to the identity, all words are translations by jv, with j in
the prime subfield. Within the distinguished line, only the points −jv
can eventually encounter zero. At depth k their individual zero-occurrence
positions separate precisely the listed points; all remaining points have
the same nonzero observation at every queried time.

For q=p prime, the last residual becomes a singleton after p−2 strict rounds
(the p=2 initial-discrete case is separate). This includes the displayed
three- and five-round examples at p=5 and p=7 and proves unbounded repair
depth as p grows. For q>p, after p−1 rounds all p prime-subfield points are
singled, and the remaining q−p points are a stable common cell. Its size is
at least two because q>=p². Together with the q^d−q already separated points
outside the line, the stable count is q^d−q+p+1. These conclusions concern
the identity-B family, not every line-preserving B. The text correctly states
that the line-moving hypothesis is sufficient and does not claim necessity.

The relation to the earlier 64-state example is typed correctly. Its **joint**
norm observer has F4*-scalar fibres; F64 is a three-dimensional F4-space,
and multiplication by b outside F4 is an invertible F4-linear map moving
the line through 1. This specializes G1 to q=4,d=3. It does not identify
field norm quotient fibres in general with arbitrary projective quotients,
or supply a tournament realization of these operations.

## Code and replay

The checker uses an independent Cartesian coordinate model, with explicit
F4 arithmetic modulo X²+X+1. It enumerates every 2-by-2 matrix and tests
invertibility by its complete image. The counts 6,48,180,480 agree with the
GL2 orders for q=2,3,4,5. Every nonzero marked vector is considered, giving
12,288,2160,9600 line-moving cases. The proof's exact first residual is checked
as a vertex set, not only a block count. Canonical labels make the partition
fixed-point test valid, and the old label in every refinement prevents
merging. Every correctness gate is an explicit exception and survives -O.

The reviewer independently replayed normal and optimized runs; both pass
**38,705 gates**, agree byte for byte, and reproduce the source receipts.
An additional direct implementation used cross products to recognize
projective directions over prime fields, rather than the producer's
first-nonzero-coordinate normalization. With v=(1,0), the observed sequences
were:

| q | Coordinate-swap B | Identity B |
|---:|---|---|
| 3 | 5→9 | 5→9 |
| 5 | 7→23→25 | 7→23→24→25 |
| 7 | 9→45→49 | 9→45→46→47→48→49 |
| 11 | 13→113→121 | 13→113→114→115→116→117→118→119→120→121 |

This small control corroborates both the line-moving theorem and the
unbounded-depth hostile. It is not a further complete GL2 census at q=7 or
q=11. The all-q,d conclusions rest on the reviewed elementary proof.
