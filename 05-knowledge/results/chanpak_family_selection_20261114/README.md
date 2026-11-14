# Target-first Chan–Pak selection: switch the quotient or break the seam

**Status: exact new witness families and exact rejections within the
specified families; universal coverage is not proved.**

The [previous factorization](../chanpak_absolute_reflection_20261114/README.md)
asks for physical factors in a singleton cyclic join. Instead of
trying to make a fixed family solve every target, this session derives
different exact selectors, each with a literal witness on acceptance
and a proof of what its rejection means.

1. **[All three uniform cyclic blocks](CYCLIC_RUNS.md).** Compress a
   Hamiltonian path to its word of module runs, then count actual
   ordered path partitions for each label. The nonnegative formula
   allows a third path-cover channel when the singleton becomes
   a pair, reaches 45, 137 and 159 outside all singleton joins,
   and still proves a structural failure at 63. A split seam
   reaches the prime 97 outside the entire uniform cyclic class.
2. **[Prime quotient factorization](PRIME_QUOTIENT.md).** One
   five-vertex quotient has two different exact inverses. With
   three transitive widths it asks for an actual factorization
   of N+2^c−2, and realizes 63 via 65=5·13. With two transitive
   widths it asks for 4N−11 and realizes prime 97 via 377=13·29.
   The least strong quotient order at 63 is exactly five, even
   allowing arbitrary tournament modules in smaller quotients.
   A different six-vertex quotient gives 56·2^n−15, including
   both 97 and 881. The same word method gives a necessary
   recursion-depth scale for any fixed strong-quotient bank with
   transitive leaves.
3. **[Change exactly one arc](ARC_READER_SELECTOR.md).** A broken
   cyclic seam gives a uniform two-path repair, and a different
   endpoint reduces the absolute target to its 2-adic width and
   a strictly increasing odd mantissa. Literal examples include
   primes 59 and 83 and 119=7·17. Both endpoint formulae are
   replayed over every width through 11 and every reversal site.
4. **[Two sources the old reader cannot distinguish](THREE_SLOT_FALLBACK.md).**
   At fixed (H,F₂)=(345,5640), two literal sources have different F₃.
   The third-slot reader gives 202425 and 205305, and padding
   realizes the exact separation 2880·3^s. This explains why
   extending just the old two-coordinate factor bank loses
   real choices and cannot emulate the third-slot test.

The run-coefficient calculation is finite by construction. The
quotient words enforce which runs really meet. The contents of each
run are actual path covers; no residue or formal signed support
is substituted for an absolute realization.

The remaining issue is target selection among an expanding bank.
The new equations are sufficient and exactly invertible on each
listed chart, not a proof of a successful chart at every prescribed
odd target. In particular a fixed menu of transitive-width charts
is too sparse. A viable coverage proof must establish a genuinely
growing substitution or provide a different exact construction
for its remaining targets.

Run the companion JavaScript or Python replays in this folder.
For example **node replay_prime_quotient.js 63** returns the
transitive widths of the first prime quotient chart at 63.

**Portfolio executable:** run `node select_family.js TARGET` for an
exact inverse on the included prime-quotient, three-slot transitive
and single-reversed-seam charts. It returns a named construction and
literal outgoing masks (decimal strings for arbitrary precision), or
reports that the tested charts missed. Direct path DP checks thirteen
accepted examples: 1, 5, 11, 13, 45, 59, 63, 83, 97, 119, 137, 159,
and 881. The executable is a partial selector, not a universal search
or a certificate of impossibility for a rejected scalar target.
