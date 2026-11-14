# Session 326: exact physical sources, clock gates, and least resistance

**Chan–Pak remains open.** This session constructs actual (H,Q₂) source pairs and proves an exact inverse in named grammars. It also separates absolute chart failure from missing residues or from a genuine spectrum hole.

## Constructing the pair

- [Arbitrary-Q inverse](ARBITRARY_Q_INVERSE.md): read either end of a strong one-apex cut and enforce the actual two-slot carries. Given H and Q₂, at most two chain lengths can work; the guaranteed branch bound at length n is O(6^(n/4)). A divisor-lattice recurrence finds the minimum-order source in the ordinal closure, keeping both coordinates. The easy h=3 factor can lead from the attainable (9,336) into the impossible (3,28), so local least-cost selection is not a valid inverse.
- [Growing physical source alphabet](PHYSICAL_SOURCE_ALPHABET.md): symmetric changes of a fixed binary word cancel the two-slot response after arbitrarily much wrapping, while retaining distinct H and a nontrivial third-slot difference. The literal pairs produce an infinite usable alphabet and a complete factor-free subpolynomial entry algorithm in its ordinal closure. That same capacity bound proves this alphabet alone is insufficient.
- [Matched-Q source switches](SAME_Q_SWITCHES.md): two different H values at one shared Q can pick different digits and meet at the exact same target. The literal pairs produce, among other values, prime common targets 25,439 and 7,509,013. Their ordinal products make growing rows of selectable H values at a fixed clock.
- [Source sieve](SOURCE_GATES.md): the real digit-specific phases and parity leave only 1050 of 5760 unit classes of prospective H, for any fixed target coprime to 99. The survivor still needs an exact physical source. At a fixed H=23 and order six, Q=272 is locally smaller but only Q=296 hits the prescribed 30,017.

## 30,000, the primorial, and the near-100 comparison

[The 30,030 boundary](PRIMORIAL_SWITCH.md) proves an exact obstruction to the entire relaxed 99 chart at 30,013, 30,029, 30,031 and 30,047, not just at one selected source. A different exact star cut supplies strong realizers of all four. The intermediate 30,017 really comes from source (23,296), digit 31. Q=30,000 is physically possible, for instance at (1827,30000), but forces H to be divisible by three for any tournament. Q=30,030 and Q=198 cannot be two-slot counts at all. In contrast 100 and 200 are clocks of one exact H=11 core, before/after a neutral singleton shift.

[The 99/100 near rails](NINETY_NINE_NEAR_RAILS.md) revisit the repo's older 30,000-range Paley and least-resistance work. Two strong families of exactly the same order realize 99·2^m+23 and 100·2^m+23. The cheaper H choice becomes the more expensive after either selective continuation. The exact 99/198 doubling is a retiming, but the step from 99 to 100 changes the needed structure: all-transitive sides cannot hit the +23 ray in the slope-99 class, whereas one cyclic SCC suffices.

## What order and arithmetic can, and cannot, tell us

[The base-four carry cone](SLOT_CARRY_INVERSE.md) exactly solves the nonnegative polynomial relaxation of the two moments at fixed size. Zero slack determines the full response, and at first slack all higher responses are totally ordered. [A physical next-slot collision](NEXT_SLOT_COLLISION.md) shows why this is not general: two strong order-eight graphs share (345,5640) yet Q₃ differs by 720, seen by either selective continuation.

[Joint modular realization](JOINT_RESIDUE_SOURCE.md) prescribes both coordinates simultaneously modulo any modulus prime to six, even by a strong tournament and with zero divisors permitted. This is not an absolute source theorem: the very same constructive grammar has a subpolynomial absolute capacity. The next unresolved proof obligation is to place an unbounded and sufficiently rich family of physical profiles on each target's exact affine line. The natural least-resistance rule applies to compatible full profiles, not isolated numbers.

Run the Python replay files in this directory; the phase table is pinned here and the S325 replay validates the literal templates. All finite witness statements are scoped to their stated catalogue or controls.
