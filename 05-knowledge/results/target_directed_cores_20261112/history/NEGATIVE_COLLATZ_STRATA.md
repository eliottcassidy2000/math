# A fixed-valuation conjugacy, with exact physical boundaries

2026-11-13. **PROVED:** the arithmetic conjugacy and branch conditions
below. **FINITE-EXACT:** the attached arithmetic and literal graph controls.
The two forbidden counts and the ready module families are inherited from
the explicitly linked sources. No Collatz cycle classification, universal
graph lift, or complete tournament spectrum is claimed.

## 1. The retained coordinate is the valuation of H+1

For every positive odd integer H, set

    a=v_2(H+1)≥1, u=(H+1)/2^a odd, m=(u+1)/2.

For fixed a the map

    Φ_a(m)=2^a(2m−1)−1                                  (1)

is a bijection from positive integers onto the stratum
S_a={H>0 odd:v_2(H+1)=a}. Use the **ordinary, unshortened** minus map

    C_−(m)=3m−1 if m odd, and m/2 if m even.

Then Φ_a conjugates C_− to the following deterministic map on S_a:

    H ↦ 3H+2                  if u=1 mod4,
    H ↦ (H−1−2^a)/2           if u=3 mod4.              (2)

Proof: m odd is equivalent to u=1 mod4, and
Φ_a(3m−1)=3Φ_a(m)+2. If m is even, then u=3 mod4 and
Φ_a(m/2)=(Φ_a(m)−1−2^a)/2. Both outputs are positive and
their H+1 valuations are a: their normalized odd parts are 3u and
m−1, respectively. This proves the conjugacy, including its guards.

Negation m↦−m conjugates C_− to ordinary 3n+1 on negative integers.
The old [inverse-tree transport, Theorem6](../../collatz_procgen_20260924_inverse_tree_mod192.md)
already proves the signed transport for the shortcut convention. Here
the odd halving is deliberately kept as a separate edge. Likewise
[THM-4472](../../../../01-canon/theorems/THM-4472-four-vertex-reading-of-3n-plus-minus-1.md)
keeps sheet change distinct from time reversal and full converse.

This encoding differs from the older moving-table modulus m=4x±1 in
[mutation arithmetic A1](../../mutation_collatz_chanpak_20261112/arithmetic/PROOF.md).
Even a=1 here gives H=4m−3. The affine offset and guard cannot be
discarded merely because an odd branch has coefficient3.

## 2. Consecutive ternary target steps are not consecutive Collatz odd steps

The map T(H)=3H+2 preserves a on **every** odd H. In the m-coordinate
it is the unconditional affine map m↦3m−1, which reverses parity.
Consequently, after any legal odd branch of (2), the next legal branch
is necessarily its halving branch. Two consecutive uses of T can never
both be consecutive odd steps of this Collatz system.

The independent integer coordinate

    r=v_3(H+1)

counts how many inverse T steps are possible before reaching an odd
H_0 with3 not dividing H_0+1:

    H=3^r(H_0+1)−1.                                   (3)

This is an exact target-address clock, not a count of physical tournament
returns or a Collatz stopping-time theorem. In particular parity guards
must still be checked on an arithmetic T ladder. For the inherited462
ladder, H=462·2^k−5 has a=2 and even m when k≥3; in that range even its
first proposed T step is not the deterministic odd branch of (2).

## 3. The failure to lift is already literal at23 and45

At H=23, a=3 and m=2, so the legal halving step gives7.
At H=45, a=1 and m=12, so it gives21. Both source counts are attained:

* H23 has the six-vertex outgoing rows (48,1,3,7,14,30), independently
  counted by the attached checker.
* H45 is C3→C3→T4, where T4 has rows (2,4,9,3) and H5.

[THM-1370, the all-order exclusions](../../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md)
excludes both target counts. Thus there is no tournament operation
defined on **every** attained source count which realizes the guarded
halving branch. Changing root or source at those same counts cannot
repair this particular obstruction. This is stronger than a failed
small mutation menu and weaker than excluding useful restricted lifts.

The m-cycle1→2→1 embeds in every stratum as

    2^a−1 → 3·2^a−1 → 2^a−1.

Its a=3 copy contains the forbidden count7. This single explicit cycle
is enough for the boundary; no exhaustion of minus-map cycles is used.
The old [powers-of-two plus-sheet lift](../../mutation_collatz_chanpak_20261112/returns/COLLATZ_LIFT.md)
is an actual, differently guarded birth/deletion construction. It does
not supply the halving operation in (2).

## 4. The405 ray records a second physical rearm failure

The inherited [module return, §§1–4](../../cycle_opening_six_return_20261112/compiler/MODULE_RETURN.md)
supplies indegree-one unit-ready sources for every k≥1 at both

    H_0(k)=45·2^k+7,
    H_1(k)=135·2^k+23=3H_0(k)+2.

Two more births attain H_2(k)=405·2^k+71. At k=17, however, the current
[complete inverse](../TARGET_INVERSE.md) excludes **every** indegree-one
unit-ready presentation of H_2(17)=53084231. Its count remains attained.
Moreover

    H_2(17)+1=9(45·2^17+8),  v_3(H_2(17)+1)=2.

Thus the failure occurs precisely at the second requested rearm of this
ternary target ladder, after two existing ready source charts. This is
a physical-interface failure at a specified arithmetic depth; valuation
alone does not imply that all depth-two targets fail.

Here a=3 throughout. The m-coordinates of H_0,H_1,H_2 at k17 are
368641,1105922,3317765, with parities odd,even,odd. Only the first ternary
step is the legal odd Collatz edge. From H_1 the legal arithmetic halving
instead gives135·2^16+7=8847367. This explicitly prevents conflating the
two-step tournament rearm question with a negative-Collatz orbit.

The already constructed [selective-port bridge](../AT_TARGET.md) reaches
the next ready target3H_2(17)+2 despite the missing ready interface at
H_2(17). It uses a different physical operation, so neither the present
conjugacy nor the unit-ready obstruction rules it out.

## Map and reproduction

Source: positive ordinary3m−1 dynamics. Target: odd integer labels in a
fixed H+1 valuation stratum. Map: Φ_a. Preserved: exact directed edges,
their parity guards, positive sign, cycles and full orbit words. Lost if
only H is retained: the tournament, marked root, defect, and legal graph
operation. Required sidecar for a graph lift: an actual construction at
each guarded edge. Cheapest decisive tests:23→7 and45→21.

Run `python3 -B negative_strata.py` and
`python3 -O -B negative_strata.py`. The script checks every a=1,...,12 and
m=1,...,4096, the literal H23/H45 sources, the displayed cycle and both
ray clocks. These finite controls accompany the all-input algebra above.
