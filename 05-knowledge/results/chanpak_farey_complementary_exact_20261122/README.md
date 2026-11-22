# The adaptive complementary placement at exponent 61

**2026-11-22, euler.** Let \(E=2^{60}(2^{61}-1)\). The published
canonical Farey tower left spines 4 and 8 open for a literal
\(E-1\leftrightarrow E+1\). This packet closes the proposed adaptive
cube-root placement on both spines but not the unrestricted search.

| Result | Exact scope |
|---|---|
| [Cube-root window and all-aspect terminal sectors](CUBE-WINDOW.md) | **FINITE-EXACT.** No pair of arbitrary words in the entire inherited placement window works, even allowing every dyadic scale and dropping the modulus restrictions. Across all aspect ratios, a sum ending in 79 ones (spine 4) or 73 (spine 8) is also impossible. |
| [Shorter complementary controller](SHORT-COMPLEMENT.md) | **PROVED.** Equal-address switch triples produce a common suffix of length \(O(\log^2(2D))\) for any modulus \(D\), with independent binary and odd programming. This is a congruence result; it cannot undo the absolute exclusion. |
| [Wider antidiagonal search](BALANCED-SECTORS.md) | **FINITE-EXACT.** Prime-by-prime parity screens and a centered integer quadratic extend the cleared balanced terminal-one sectors. |
| [Spine-eight exact sweep](SPINE8-SWEEP.md) | **PROVED algorithm; finite-exact new interval.** No floating-point rejection is needed. The entire computation remains far larger than the interval actually run. |
| [Four vertices and Conway's priority](FOUR-PRIMEGAME.md) | **FINITE-EXACT.** Two distinct height-three states form opposite corners of a physical square; swapping the \(11,13\) registers is not a one-arc move. On the exact 133-vertex lower count, \(17/91\) preempts \(11/13\). |
| [13-adic clocks on the marked exact lower construction](CONTACT-CLOCKS.md) | **PROVED.** Varying formal Euclid exponents lifts the two marked contacts arbitrarily far at 13 while retaining their original prime phases; contact 53 cannot match 11 and 13 simultaneously. These are not new exact upper endpoints. |
| [A smaller-profile selector for Chan–Pak](TRIANGLE-SELECTOR.md) | **PROVED + FINITE-EXACT.** A two-module cycle lets a prescribed count select a smaller exact two-slot target. The extra slot data is genuinely necessary on a second wrap. |

An independent [spine-four corridor replay](farey-corridors.py) also clears 29,457 disjoint normalized rational-tangent corridors for each allowed left residue (3,204,528,847 coordinate values per residue).

Neither the canonical exact exponent-61 flip nor Chan–Pak is settled.
In particular, searching a shorter suffix controller within the cleared
cube-root window cannot help, even if it matches a larger modulus.
The all-aspect inverse outside that window remains the right test for
a proposed exact pair. The 133-vertex backchains are exact lower
witnesses, not claimed minimum-order witnesses; their known one-flip
exclusions are unchanged at exponent 61.

Python and C++ replays sit beside the notes. Each finite search is
reported with the range actually run, separate from the complete
algorithm and its unrun budget.
