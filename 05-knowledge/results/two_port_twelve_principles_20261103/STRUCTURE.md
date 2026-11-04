# The twelve-vertex structural split and its resolved pair branch

**PROVED structural implications and T5 identity; FINITE-EXACT order-twelve
reduction using the explicitly inherited complete banks.** The original
necessary-congruence split below is retained. The incoming exact twin sign
closes its homogeneous-pair branch for hidden ambient module converses;
the prime case remains **OPEN**. This is not an exhaustive enumeration of
twelve-vertex tournaments.

Write R(T) for the necessary conditions S=E modulo two and C=C^t modulo four.
By inherited B3-B4, any module whose converse is invisible in some ambient
tournament must satisfy R. A combined survivor would still need to pass
the full response, including [T1-T4](PROOF.md).

## Structure that follows before any order-twelve search

1. **R(T) forces strong connectivity when |T|>=2.** A Hamilton path in a
   nonstrong tournament starts in the first strongly connected component and
   ends in the last. Those components are disjoint. Since sum_u S_u=H is odd,
   some vertex in the first has odd S and zero E, contradicting R.
2. **Every module inherits C=C^t modulo four.** The uniform-run derivative
   formula in B4, independent of its later collision assumption, is
   C_uv(T)-C_vu(T)=H(Q)(C_uv(M)-C_vu(M)) modulo four for u,v in M.
   H(Q) is odd. In particular no module can have odd size at least three,
   by the inherited odd-module current theorem A1.
3. **There are no modules of sizes four or eight.** The unfiltered complete
   C-skew census below has no symmetric-mod-four core at either size. At
   size six it has only C3[TT2,TT2,TT2], which contains homogeneous pairs.
   These are finite bank consequences, not extrapolations to other even sizes.

Dependencies: [B3-B4](../two_port_seven_core_20261103/algebra/PROOF.md),
[A1 odd-module current](../two_port_beyond_nine_20261103/algebra/PROOF.md),
and the complete-bank receipt [cover_scout.json](cover_scout.json).

## T5. Adding the two exterior vertices of a directed triangle swaps endpoint parity

For any nonempty M, form T=C3[M,{a},{b}] with M->a->b->M.
Let E2_u and B2_u count ordered two-path covers of M with one fixed slot
ending or starting at u respectively. Then exactly, for u in M,

    E_u(T)=E_u(M)+E2_u(M),
    S_u(T)=S_u(M)+B2_u(M),                         (1)
    S_a(T)=E_a(T)=S_b(T)=E_b(T)=H(M).

**Proof.** A Hamilton word ending in M has either the skeleton a,b,M or
M,a,b,M. Their contributions are the two terms of the first identity.
The first and last M-runs in the latter skeleton occupy distinct fixed
slots. Starting in M is the reversed skeleton argument. If an endpoint
is a or b, the remaining M vertices occupy one run; the unique compatible
exterior order gives H(M). There are no other legal skeletons. QED.

Inherited squarefree P2 gives E2_u=B2_u=E_u+S_u modulo two. Therefore

    E_u(T)=S_u(M),  S_u(T)=E_u(M) (mod 2).         (2)

Thus T is endpoint-parity balanced iff M is. If T also has symmetric C
modulo four, M inherits that property by the internal derivative identity.
Hence R(T) implies R(M). This is a necessary descent, not an iff for all
of R and not an equality of exact response functions.

**Twelve-vertex consequence.** A strong twelve-vertex R(T) cannot have a
ten-vertex module: its three-vertex quotient would be strong, hence C3.
T5 would give R on the ten-vertex core, excluded by the complete
[9,733,056-class sieve](../two_port_seven_core_20261103/residual/RESULT.md).
Every proper nontrivial module therefore has size 2 or 6. In the latter
case it contains a homogeneous pair, which remains a module of T.

Consequently **any twelve-vertex hidden-module candidate is strongly
connected and either prime (no proper nontrivial module) or contains a
homogeneous transitive pair.** This is a necessary structural split, not a
proof that either branch is empty using these congruences alone.
Self-converse candidates are already excluded by inherited P4/B3; that
condition is retained.

## Incoming exact closure of the homogeneous-pair branch

[C2 in the independently audited eight-core packet](../two_port_eight_core_20261103/algebra/PROOF.md#c2-an-exact-ordered-cover-sign-recovers-the-orientation-of-twins)
supplies the missing exact obstruction. If the actual ambient tournament
has named twins a->b, contract them to z in Q. For every exterior w,

    C_bw-C_aw+C_wa-C_wb = C_zw(Q)+C_wz(Q) > 0.

For the opposite orientation b->a of the same named twin pair the contrast
is strictly negative, even if the remaining quotient changes. Full two-port
response equality determines ordered C and thus cannot identify these two
orientations. A homogeneous pair in a module remains a homogeneous pair
in every uniform ambient substitution; reversing that module reverses the
pair. This excludes every hidden module converse containing such a pair.
At ambient order two the exact endpoint identities separate the orientations.

Consequently [C4 of that packet](../two_port_eight_core_20261103/algebra/PROOF.md#c4-a-genuine-twelve-vertex-converse-collision-must-have-a-prime-core)
closes the pair branch of the split above: a twelve-vertex module whose
reversal is invisible inside any ambient tournament must be strong, prime,
and nonselfconverse. The prime residual remains OPEN. The [independent
review](../two_port_eight_core_20261103/recovery/ALGEBRA_RESPONSE_REVIEW.md)
records the exact identity and finite dependencies.

This conclusion does not upgrade the necessary congruences R(T) into a
complete response test. Nor does the cubic battery of T4 apply directly
to arbitrary ambient internal-module reversals: T4 compares a whole core
with its whole converse. Ambient equality gives intrinsic response equality
modulo four by T1; exact ambient separation here comes from the named twin
contrast itself.

## Complete small cover-skew census

|Order|Inherited complete classes|C=C transpose modulo4|
|---:|---:|---:|
|2|1|1|
|4|4|0|
|6|56|1|
|8|6880|0|

No endpoint-parity or strongness filter is applied. The sole order-six
class has code22693 in the inherited pair-bit convention and H45; literal
permutations and every two-path cut independently verify its cover matrix.
The separate structure checker explicitly finds its isomorphism to
C3[TT2,TT2,TT2]. The bank is SHA256 pinned; its isomorphism completeness is
inherited from the authenticated representative-bank packet, not reproved
by this script. Normal and optimized receipts are separate replays of the
same census implementation.

## Named order-twelve controls and a stopping witness

The [declared-family scan](twelve.json) checks all 2048 named extensions
of each of transitive11, cyclic11 (steps1..5), and Paley11, plus every
TT2 inflation of the 56 order-six quotient classes: 6200 presentations.
There is no input filter or claim of disjoint isomorphism classes. All fail
endpoint balance. The new shared-label cubic test excludes 6198; its two
survivors are Paley11 with a universal source or a universal sink. These
are exact counterexamples to using that cubic projection alone. Neither
passes the combined residual, and neither is a full-response collision.

T5 gives a positive control for the endpoint filter at twelve: apply the
cyclic extension twice to the known parity-balanced eight-vertex hostile.
Its rows are

    (1284,1297,1322,1411,1325,1355,1311,1399,1536,1279,2048,1023).

It has S=E modulo two but C_02-C_20=28638=2 modulo four and nonzero shared
cubics. This control exercises the endpoint-survivor branch; an empty
survivor list in the four families is not the only code path tested.
The C++ parity engine is checked against independent integer endpoint and
cover DP and the unsimplified cubic formula in [audit.json](audit.json).

Reproduce `scout.py`, `structure.py`, `run_twelve.py`, and `audit.py`, each
with `--output <receipt>`; use `python -O` for the companion receipts.
These finite tests motivate the structural split; only the proofs and
complete-bank scopes support its stated consequences.
