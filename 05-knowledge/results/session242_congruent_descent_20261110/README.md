# Session242: congruent arithmetic and physical Chan–Pak descent

**COMPLETE, 10 November 2026. Independently audited. Chan–Pak remains OPEN.**

The original numbers have an exact common source:

    41²−31² = 49²−41² = 720 = 5·12².

They give the area-5 rational right triangle (3/2,20/3,41/6), the primitive
triangle (9,40,41), and the earlier smallest positive AP square-sum clique
(482,3362,6242). This session explains the distinct scales 6,12,24, strengthens
the denominator statement, and tests which descent operations actually retain
the property one needs. It does not identify elliptic curves with tournament
path counts.

## Audited results

1. **Twelve is forced; one elliptic halving is all this point permits.** Every
   common denominator of a rational square AP of difference 5 is divisible by
   12. At denominator 12, the only positive increasing numerators are 31,41,49.
   The natural triangle point P=(25/4,75/8) and direct AP point
   Q=(1681/144,62279/1728) satisfy Q=−2P, with Q∈2E₅(ℚ) but Q∉4E₅(ℚ).
   The two maps have different positivity chambers. Classical correspondences
   are credited; no literature firstness is claimed.
   [THM-5106](../../../01-canon/theorems/THM-5106-gap-five-minimal-denominator-and-exact-elliptic-halving-depth.md),
   [full arithmetic proof](arithmetic/ARITHMETIC.md), [audit](audit/REVIEW.md).

2. **A new obstruction on physically valid graft inverses.** For every proper
   rooted graft, output N and gadget pair (H,Q) satisfy
   N/gcd(H,Q)≠7. Therefore a target 7p with odd prime p≠7 can use gadget content only 1 or 7,
   not the other prime p. This is stronger than merely saying the gcd divides
   N. It relies on an actual normalized gadget, not an invented integer pair.
   Both allowed branches occur in strong examples 791=7·113 and 1673=7·239.
   Removing properness breaks the theorem: a sink root gives 49/7=7.
   [THM-5107](../../../01-canon/theorems/THM-5107-physical-graft-content-excludes-cofactor-seven.md),
   [proof and sharp witnesses](exception/PROOF.md), [audit](graft_audit/REPORT.md).

3. **A complete constructive descent on a declared infinite sector.** Every
   {3,5,7,31,41}-smooth positive odd target except 7 and 21 is produced by an
   explicit ten-graph bank. The descent preserves the predicate “this target
   has no tournament witness”: every smaller witness lifts by a universal
   ordinal join. Literal graph suppliers for 217 and 287 are certified.
   This closure is an application of inherited atom reduction, not newly
   discovered global coverage. The additional exact cost theorem determines
   the minimum order within this particular bank; for fixed other exponents,
   adding three powers of 7 eventually costs exactly eleven vertices.
   [THM-5108](../../../01-canon/theorems/THM-5108-five-prime-physical-descent-and-exact-join-bank-cost.md),
   [descent and cost proof](descent/PROOF.md), [audit](audit/REVIEW.md).

## What infinite descent would still have to do

The inherited prime/7p reduction already isolates the infinite obligations.
Our smooth-sector proof supplies a working model: a strictly smaller admissible
missing target and a constructor lifting *any* smaller witness. Generic graft
inversion does not yet have that universal lift, since equal Hamilton counts
can hide different rooted responses. Its finite, decreasing search still
needs a theorem guaranteeing a physical predecessor in every required fiber.

For prime targets the gadget must be primitive. For 7p targets the sharper
permitted contents are 1 and 7. Those are necessary filters, not existence.
The congruent-number lesson is exact: Berggren descent changes class 5 to 21
and loses positivity; elliptic halving stops at a nondivisible point; integer
gcd division of a tournament interface can leave the physical category.
A proof must retain its actual predicate and lift, not only a decreasing size.

[Research synthesis and next obligations](SYNTHESIS.md) keeps the successful
sector, the sharper obstruction and the remaining global gap separate.

## Reproduction and inherited sources

From this directory, run:

    python3 -I -B verify_all.py > replay.json
    cmp replay.json verification.json
    sha256sum -c SHA256SUMS

Ten canonical normal/optimized replays and four deliberately failing controls
pass. Independent algorithms include chord addition, divisor enumeration,
forward endpoint DP, recursive start-vertex counting and full permutation
counts of both ten-vertex graft outputs. No numerical cutoff is promoted to an
unbounded theorem. [Source pins](source_pins.json) identify inherited blobs at
baseline 68e5a06fbb520a861eb2e094efbc40a82cc52b1a. No old repository-wide
computation or full imported dependency chain is claimed freshly reverified.

Early checkpoints preserved the inheritance and arithmetic proof before this
final synthesis. The finite universe, hostiles, equalities and scope limits
are given in each lane. Session242 has ended.
