# Independent audit of the finite O/Z tail-bank reduction

Status: PASS for the corrected elementary all-order reduction and its restricted necessary density condition.

Reviewed proof SHA256: `6a20f66b4a85c48c70f1436b4c7b6556c44bffbcaa93381a935dae02527f4030`.

## Exact identities

The binary orientation is correct: chronological O/Z letters are encoded low bit first, and terminal O letters are high zero bits. Removing them preserves h. For m>0, at least one Z occurs and raises h, so the nonconstant-tail assertion holds for physical nontrivial unit-start seeds.

The formula for R before each step follows by iterating R'=2R-db. Summing b_i(R_(i-1)-e_(i-1)) gives exactly the pair sum defining Q, with the separate e0 term only at the first digit. The adjacent-pair coefficient is correctly doubled. Removing a low zero bit or adding a low one proves the displayed recurrence, with no reversal of binary orientation.

Q is nonnegative and even; its upper bound Q(m)<=m-1 for m>=1 follows by induction, separating m=1 from the positive-parent odd case. The corrected stronger lower bound explicitly requires m>=1. The weaker bound includes m=0 separately, so the inverse-search bound is valid for all possible preimages. R0>d for physical edge-grammar sources is an inherited positive-start inequality. The singleton is not an allowed seed here and would violate that strict gap.

## Covering and normalization scope

The arithmetic equivalence is exact for the O/Z tails of the specified finite bank. The revised text correctly distinguishes it from the unrestricted eventual-spectrum conjecture: a covering bank is sufficient, not known necessary. No physicality claim is made for arbitrary nonnegative integer states.

Chronological J followed by Z and O followed by J both produce (h+2R,h+2R,h,h). Thus JZ=OJ is a response identity, without an asserted graph isomorphism. The revised directional statement about the rule JZ -> OJ and the unreduced subword JO is accurate. No finite-prefix-state bound follows merely from this identity.

For a fixed finite bank, the inverse bound gives at most X/(R_s-d_s)+O(1) distinct outputs per seed. Summing and comparing to X/2+O(1) proves the strict reciprocal-width obstruction. Overlaps cannot invalidate this upper bound. The equality boundary is not decided, and the theorem is not a general finite-bank obstruction.

## Independent controls and exact finite scope

The separate `verify.py` in this audit directory generated142 physical seed responses and tested every encoded tail0<=m<=511:72,704 cases. It checks the direct pair-sum formula independently of the Q recurrence, parity, upper/lower bounds, inverse bound, exact response orbit, and the JZ/OJ identity. `RESULT.json` records PASS.

I reran the author's verifier and reproduced its receipt. I also compiled `image.cpp` afresh and reproduced the quoted finite root-image totals:249,077 odd outputs at most999,999, and24,904,080 at most99,999,999. `ROOT_IMAGE_RECOMPILED.out` records the replay. The finite root-image computation uses the proved inverse horizon and is properly labeled exploratory rather than an asymptotic density or noncofiniteness proof.

All requested corrections are incorporated. No mathematical blocker remains. A status-only promotion to independently audited in precisely the stated scope is approved.
