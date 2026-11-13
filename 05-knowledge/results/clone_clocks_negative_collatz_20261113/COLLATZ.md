# 3H+2, the Collatz maps, and the two anchors

klein, 2026-11-13. **PROVED (elementary) / VERIFIED.** No Collatz statement is claimed; every identity below is an
exact algebraic fact checked in `verify.py` section 1. The shadowing lemma and the 3x−1 ⇔ negative-integer
equivalence are classical (Terras 1976, Everett 1977, Lagarias 1985); the Størmer/Gersonides reading of the balanced
cycles is inherited from [HYP-9339](../../hypotheses/HYP-9339-the-two-three-axis-numbers-fixed-points-and-relations.md).

## 1. C2 is a conjugate of the Collatz step

Write `T(x)=(3x+1)/2` for odd `x`, `T(x)=x/2` for even `x` (on all of Z), and `D(x)=2x+1` (append a binary 1).
Hamiltonian path counts are odd (Rédei), and on odd `H`

    C2(H) := 3H+2 = T(D(H)) = D(T(H)).

So the C2 child of a unit-ready source is "Collatz after appending a 1 bit", and equally "append a 1 bit after
Collatz". In the coordinate `r=(H+1)/2` (used in [mrb_seven_cutoffs](../mrb_seven_cutoffs_20261113/README.md), where
the forbidden counts 7, 21 become 4, 11):

| map | on H | on r |
|---|---|---|
| C2, digit σ (two births at a σ-ready root) | 3H+2σ | 3r+(σ−1) |
| D (one dyadic clone, on the −1 anchor) | 2H+1 | 2r |
| T on odd H | (3H+1)/2 | 3r/2 |

The three C2 digits σ=0,1,2 are therefore exactly `r↦3r−1, 3r, 3r+1`. The digit-2 branch is the odd step of
the classical `3x+1` map; the digit-0 branch is the odd step of the `3x−1` map, which is `3x+1` on the negative
integers under `x↦−x`. Pure tripling, the digit-1 branch, sits between them. Every `r≥2` has exactly one C2
parent `(r−σ+1)/3`, and the affine offset `σ−1` is nonzero exactly when `3∤r`. Hence the sequence
`2,4,5,7,8,10,11,13,14,16,…` (non-multiples of three) is precisely the set of `r` whose last C2 step is a Collatz
step (`3r±1`); the multiples of three are reached by pure tripling. The same residue split governs Collatz
itself: an odd number has a `3x+1` (or `3x−1`) predecessor iff it is prime to 3.

The two forbidden counts sit at the very top of this tree: `r(7)=4=3·1+1` is the `3x+1` child of the
transitive root `r(1)=1` (the trivial Collatz step `1↦4`), and `r(21)=11=3·4−1` is the `3x−1` child of the
missing node 4. This locates the exceptions; it explains nothing about why they are exceptions.

## 2. Dyadic families shadow Collatz cycles, with Størmer exchange rates

A dyadic family `c·2^k+d` tends 2-adically to its anchor `d` as `k→∞`. By the classical parity-vector lemma,
`T^i(c·2^k+d)=T^i(d)+3^{a_i}c·2^{k−i}` for `i≤k`, where `a_i` counts odd steps among the first `i` steps of
`d`. If `d` lies on a `T`-cycle of length `L` with `A` odd steps, then

    T^{Lt}(c·2^{Lt}+d) = c·3^{At}+d     (all integers c, t ≥ 0).

| anchor `d` | cycle | `(L,A)` | exchange | Størmer pair |
|---|---|---|---|---|
| −1 | `{−1}` | (1,1) | `2→3` | 3−2=1 |
| +1 | `{1,2}` | (2,1) | `4→3` | 4−3=1 |
| −5 | `{−5,−7,−10}` | (3,2) | `8→9` | 9−8=1 (Catalan) |
| −17 | 11-cycle | (11,7) | `2^11→3^7` | unbalanced (3^7−2^11=139) |

The inherited constructor chain [FAMILIES F3](../multicore_metagraph_20261113/FAMILIES.md)

    462·2^k−5 → 1386·2^k−13 → 4158·2^k−37

moves its anchors along the C2 orbit `d↦3d+2`, i.e. `(d+1)↦3(d+1)`: `−5→−13→−37`. The first anchor lies on the
balanced Catalan cycle, so `T^{3t}(462·8^t−5)=462·9^t−5`: Collatz turns that tournament ray's dyadic clock into
a ternary one at the exchange rate `8→9`. The third anchor `−37` lies on the unbalanced `−17` cycle, and `−13` is
preperiodic into the `−5` cycle (`−13→−19→−28→−14→−7`). These memberships are verified facts; nothing
is claimed to depend on them. The source chain itself is reproduced in the funnel compiler (`verify.py` §2).

## 3. The anchors are fixed points, and their sign selects the Collatz variant

Every exact chart considered below has the form `H = ε + α·2^k·3^j` or `H+1 = 3^j F(k)`. In the coordinate
`y=H−ε` the three physical clocks are multiplications: dyadic growth `y↦2y`, ternary growth `y↦3y`, and the
clock transfer `(j,k)↦(j+1,k−1)` is `y↦3y/2`.

* On the **−1 anchor** (`ε=−1`): `y=H+1`, ternary growth is `H↦3H+2` (C2, digit 1), dyadic growth is `D`, and
  `y↦3y/2` is the `3x+1` odd step `H↦(3H+1)/2`.
* On the **+1 anchor** (`ε=+1`): `y=H−1`, ternary growth is `H↦3H−2`, dyadic growth is `H↦2H−1`, and
  `y↦3y/2` is the `3x−1` odd step `H↦(3H−1)/2`. Reflecting `x=−H`, these are `x↦3x+2`, `x↦2x+1` and the
  `3x+1` step **on the negative integers**: `T(−H)=−(3H−1)/2`.

The classical odd-branch criterion makes the residues explicit. For positive odd `n` and `s≥2`:
`n≡−1 (mod 2^s)` iff the first `s` applications of `T` use its odd branch, and `n≡+1 (mod 2^s)` iff the
first `s` applications of the analogous minus map use its odd branch. These are also exactly `s` strict
ascents for the plus map, and for the minus map when `n>1`; its fixed point `n=1` is not a strict ascent.
Thus `4x±1` and `8x±1` specify two and three consecutive odd-branch applications, respectively, with the
sign choosing the variant (`verify.py` checks `s≤5`, odd `n<4000`, and strict growth separately).
Indeed along such a run the iterates are `3^j(n+1)/2^j−1` or `3^j(n−1)/2^j+1`; the first `s` inputs
are odd exactly under the stated congruence. Positive even inputs decrease, so these parity guards
are necessary for strict growth as well. This corrects the earlier `s−1` ascent wording and fixed-point omission.
The doubly multiplicative
two-clock grids of [CLOCKS](CLOCKS.md) (`σ_2=0, σ=−1`) have amplitude divisible by 8 and anchor +1, so they live
in `8x+1`: the `3x−1` (negative-integer) ascent class. The C2 line with a free dyadic parameter would need `8x−1`.
It is absent, although the digit-3 family Ψ (anchor −3) shows the −anchors are not excluded wholesale. Section 4 of CLOCKS realizes the negative-integer Collatz step by **one arc reversal**.
