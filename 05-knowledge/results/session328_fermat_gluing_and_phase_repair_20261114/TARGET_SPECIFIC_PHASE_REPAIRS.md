# A target-specific phase repair costs less than a universal companion

Let the slope-99 chart read a source pair (H,Q), with digit c the odd representative of N/H modulo198 when H is coprime to99. S327's complementary-pair theorem asks for two sources whose eligibility classes partition the entire first-phase ring. A **target-specific repair** asks less: at this particular N, the first source forces an impossible effective clock, and a different source with the same exact available Q can hit N. Distinguishing the two notions shows why large group orbits are not always needed.

## Two actual prime repairs

| N | Old source (H,Q) | Old forced digit and effective Q | New source (H,Q) | New digit, phase |
|---:|---|---|---|---|
| 3169 | (5,32) | 119,26 | (1,32) | 1,0 |
| 25667 | (5,256) | 25,258 | (19,256) | 17,0 |

An effective clock 26 or258 equals 2 mod4 and is greater than two. It cannot be the two-slot count of any tournament; nor can it be such a count times a power of two, since ordinal singleton padding would make it a physical count. Hence the old H has a genuine phase failure, not just a bad witness choice. The actual available source clock is unchanged when switching to the new H, which is the exact one needed.

The source at (5,32) is the strong order-four core (2,4,9,3); at (5,256) add three singleton SCCs. The H=1 source is TT5; the H=19 source has rows (2,20,9,3,45,15). Digits1 and17 are phase-zero one-return stars with left masks respectively (2,48,3,7,13,29,63) and (8,9,3,100,15,23,55), apex outcuts33 and9, and an empty right block (the published variant has a removable terminal singleton). The two exact outputs are

    3169 = 99·32+1,
    25667 = 99·256+17·19.

Both are prime by trial division. The replay independently builds the output tournaments, takes their unique nonsingleton strong SCC, and counts H directly.

## Sharpness and the absolute obstruction

In the complete low-clock catalogue through Q=256, the only shared exact clocks with different unit H are powers of two: Q=32,64,128 have H values{1,5}, and Q=256 first has {1,5,19}. Thus Q=32 is the earliest possible budget for any such repair with unit H; Q=256 is earliest if both choices must have H>1. The finite completeness follows from the [small atlas proof](../session326_exact_source_profiles_20261114/PRIMORIAL_SWITCH.md); the replay recomputes the pairs from scratch. This is distinct from demanding a complementary pair that works *at every residue*: see the separate exact Q=4096 theorem in this packet.

At a prescribed N and chosen exact budget Q, the inverse can be made simpler by factoring the residual R=N−99Q. A source count H must divide R with c=R/H an odd bank digit, and the raw pair (H,Q/2^r_c) must be physical. Replacing H without replacing Q is exactly an intersection with this divisor menu. If two candidate H values are congruent modulo198 and units, they force the same c for that N; they cannot both hit N with the same Q unless they are equal. Complementary phase partners can choose which H is eligible, but exact budget reuse also needs the residual to land on the chosen source's actual cH. One cannot infer budget-free arbitrary-target entry from full residue coverage alone.
