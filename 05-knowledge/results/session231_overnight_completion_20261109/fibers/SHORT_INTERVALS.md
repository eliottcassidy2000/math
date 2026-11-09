# Twelve-bit movable-core intervals

Status: PROVED, independently audited. This is an efficiency improvement of `UNIT_INTERVALS.md`; its earlier 1188-bit construction remains valid. Unrestricted Chan–Pak completeness remains OPEN.

## Explicit all-order theorem

Let

    C=[4,1,18,7,75,31,175,63], cut124,
    X=101000111000, Y=111011100110.

The core state (h,q,A,G,T,S,E) is (45,1440,878,577,332,25,14). For k≥0 and 0≤j≤k take the deletion tournament

    L_(k,j)=s1 ⇒ s0 ⇒ X^j ⇒ C ⇒ Y^(k−j) ⇒ s1.

Bits are singleton root-cut labels; all arcs between listed blocks point forward. Add v,r with L→v→r and r beating exactly the cut vertices. Every resulting tournament is strong: r reaches the leading selected singleton and hence all of L, while every L vertex reaches r through v.

These k+1 graphs, each on 12k+13 vertices, have Hamilton-path counts

    D_k−2j,  j=0,...,k,
    D_k=(245051264·4096^k+1602510k+1118023)/24843.       (1)

Thus every k yields k+1 CONSECUTIVE ODD counts in one physical full-interface fiber. The displayed division is exact for all k, either by the path-count proof below or the integer recurrence giving (1). At k=0, D_0=9909.

Set J_k=(4096^k−1)/4095. Their common full interface is

    h=45, q=11520·4096^k,
    A=7516·4096^k+343260 J_k+45,
    G=2353+9532800 J_k+5760·4096^k,
    S=E=45.

Only the cut-coordinate T varies, with T_(k,j)=T_(k,0)+2j. The deletion h is FIXED 45 and gcd(h,q)=45. Saturation does not produce a growing divisor because the construction does not repeat nontrivial cores.

## Finite identity that proves every k

Exact literal core counting and the ordinal composition law give

    state(s0⇒C⇒Y⇒s1)=(45,23592960,7535881,16564034,5290654,0,45),
    state(s0⇒X⇒C⇒s1)=(45,23592960,7535881,16564034,5290656,0,45).

So moving the core once preserves h,q,A,G,S,E and raises T by 2. In L_(k,j), the context immediately left of each move ends 0 (X ends 0 and the left flank is0); the right context starts 1 (Y starts 1 and the right flank is1). The general boundary formula ΔT_total=ΔT+last_bit(left)ΔS+first_bit(right)ΔE therefore reproduces this SAME flanked difference at every position. Induction gives T_(k,j)=T_(k,0)+2j. Since closure count is2A+h−T, the count decreases by 2 at each move.

This proof uses the one finite word identity plus context algebra, not a finite-k extrapolation. The carry-loop derivation is useful discovery provenance but not a dependency of the all-order theorem.

## Exact baseline derivation

The selected wordY has state(1,4096,3814,1655,1540,1,0). The state of s1⇒s0⇒C has

    h45,q5760,A3758,G2353,T1486,S45,E14.

After k copies ofY,

    A_k=3758·4096^k+45·3814 J_k,
    G_k=2353+5760·1655 J_k,
    T_k+E_k=1500+1655[3758J_k+45·3814(J_k−k)/4095]+45·1540k.

For k=0 this uses T+E=1486+14; for k≥1 the firstY absorbs the14 boundary contribution and all later E values are0. Appending the final selected singleton and taking star closure gives

    D_k=3A_k+135−(T_k+E_k),

which simplifies to (1). Each intermediate expression is an integer. J_k−k is divisible by4095 because4096≡1 modulo4095.

## Context-stable growing windows

General principle: if physical C and singleton wordsX,Y satisfy a flanked full-interface identity with T difference±2, and X's last bit andY's first bit match the two flanks, then the same moving-core construction produces arbitrary consecutive-odd widths. An extra leading selected singleton ensures strongness when needed.

Moreover, for each fixed k, arbitrary COMMON transitive prefix/suffix cut contexts preserve the exact gap2 and length k+1, by the pure-T context law. Strongness is retained whenever the chosen context is legal (for instance a selected first prefix vertex). This supplies a family of movable baseline windows, not just one isolated sequence.

The intervals in (1) alone do not overlap between successive k: their centers grow on scale4096^k whereas their widths are only2k. Adding arbitrary singleton contexts preserves width but does not itself prove target occupancy. The remaining completeness obligation is an explicit baseline-selection theorem showing some such context places a window over each target odd integer outside 7,21. No such theorem is asserted here.

## Exact controls

- `SHORT12_UNIT_LOOP_WITNESS.json`: physical core, cut, X,Y, and carry gain−2.
- `verify_short12_unit.py`: literal permutation core count, every carry edge, exact flanked identity, and all positions through k=8.
- `verify_short12_independent.py`: separate six-state path-cover transfer derived from literal category counts, checking every position through k=15 against the CLOSED baseline formula. It does not use the seven-coordinate composition routine for its scalar counts.
- Normal and optimized Python outputs agree for both paths.

No minimum certificate-length claim is made. The twelve-bit loop is a direct directed closed walk; unlike the earlier 1188-bit witness it needs no Bézout concatenation.
