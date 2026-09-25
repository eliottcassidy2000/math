# Root decoding, ordered carries, and the Mahler comparison

**Status: PROVED elementary identities and obstructions; FINITE-EXACT
controls.** This note does not prove universal Collatz convergence or
Mahler's Z-number conjecture. The contribution is an explicit decoder
interface and a scoped obstruction to identifying the two problems by
their parity words. No novelty or canon promotion is claimed.

## Inheritance and live board

Closest proved mechanism: the affine block/cylinder identity in
[THM-4469, Mahler bridge adjacent block pairs](../../01-canon/theorems/THM-4469-mahler-bridge-adjacent-block-pairs.md),
and the independent real-tail/ordinary-integer conditions of
[THM-2228, Mahler carry tail and integral stabilization](../../01-canon/theorems/THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization.md).
The more exact sidecar implementation is
[THM-4072, Mahler safe-terminal fibre product](../../01-canon/theorems/THM-4072-mahler-safe-terminal-fibre-product-and-finite-state-obstruction.md):
finite local terminal-prefix tests cannot replace ordinary stabilization,
and the conversion needs an unbounded integer sidecar.

Canonical hostile: a valid infinite 2-adic parity word need not be an
ordinary positive integer; the safe word `(100)^infinity` has Mahler
state `-9/19`. Corrected near miss: Catalan's unit-gap classification is
not a cycle classification; the ordered carry cancels the non-unit gap
`-139` in the known negative Collatz seven-cycle. See
[Pillai/convergents/cycle gates](collatz_mod6_20260921_pillai_convergents_cycle_gates.md)
and the 2026-09-21 Catalan entry in `01-canon/MISTAKES.md`.
Least-used relevant sidecar: the ordered carry, rather than only the
number of odd and halving steps.

Anchor: decode an arithmetic route to root 4. Niche: the two halving
systems' opposite adjacent-exchange laws. Wildcard: their common parity
space sends the entire known Collatz root basin outside Mahler's positive
ordinary domain.

| Live concept | Preserved predicate / operation | Missing or destroyed information |
|---|---|---|
| Odd chain, doubles seam, higher powers of 2 | dyadic valuation and forced initial zero bits | later odd/even order |
| Root decoder | exact equality `3^a n+R=4*2^L` | clock `(L,a)` alone loses R |
| Tournament halving | proposed size reduction by 2 | must supply a parity/residue-compatible decoder |
| Catalan / near powers | size of `2^L-3^a` | ordered carry and its divisibility |
| Mahler | same binary itinerary space | ordinary positivity and real suffix safety |
| Generalized Mahler bridge | special block carries and confinement interval | does not transfer the original 3/2 interval automatically |

## 1. An exact finite decoder (PROVED)

Use distinct names for the maps:

```text
C(n) = n/2                    (n even),
       (3n+1)/2               (n odd);       shortcut Collatz

H(n) = (3n+(n mod 2))/2 = ceil(3n/2);       Mahler ceiling map.
```

For a length-L parity word `w=(e_0,...,e_(L-1))`, first bit first, put

```text
a = sum e_j,
R(w) = sum_j e_j 2^j 3^(sum_(i>j)e_i),
M(w) = sum_j e_j 2^j 3^(L-1-j).
```

Then on the appropriate parity cylinders,

```text
2^L C^L(n) = 3^a n + R(w),
2^L H^L(n) = 3^L n + M(w).                      (1)
```

The two recurrences are `R(we)=3^e R(w)+e*2^L` and
`M(we)=3M(w)+e*2^L`. In either case the cylinder is exactly the single
residue class that makes the right-hand side divisible by `2^L`.
One direct proof is induction: the first bit is forced modulo 2; the
odd numerator coefficient is a unit modulo every power of 2; each next
bit splits one cylinder into its two lifts modulo the next power.

Therefore the following verifier is exact:

> Given n and a binary word w, a shortcut-Collatz certificate to 4 holds
> if and only if `3^a n+R(w)=4*2^L`.

Indeed the equality implies divisibility, hence membership in w's
cylinder, and then (1) yields the endpoint. Conversely an actual route
satisfies (1). This certifies a supplied route; it does not construct one
for every n.

**Convention boundary.** Under C, `1 <-> 2` is a cycle that does not visit
4. Thus universal certificates to 4 are asserted only for `n>=3`, with
1 and 2 handled as the terminal cycle. For the ordinary map with odd
step `n -> 3n+1`, the cycle is `1 -> 4 -> 2 -> 1` instead. Mixing these
step conventions would invalidate a literal all-n decoder assertion.

**Minimal clock-loss witness.** Words `01` and `10` have identical
`(L,a)=(2,1)` but carries 2 and 1. Their proposed root starts are

```text
(4*4-2)/3 = 14/3,     (4*4-1)/3 = 5.
```

The second is the valid route `5 -> 8 -> 4`; the first is not an
integer route. The slope `3^a/2^L`, number of halvings, and block
cardinality do not distinguish these cases.

## 2. What the n+2 and doubling modes preserve (PROVED)

If `n=2^m q` with q odd, its C word begins with exactly m zeros.
Thus the odd chain is the `m=0` row, the doubles seam is `m=1`, and the
higher sea is `m>=2`; pure halving peels off those forced zero bits.
This is a genuine decomposition, but it does not specify the parity
word after q reaches an odd step.

The exact local translation rule is

```text
C(n+2)=C(n)+1          for even n,
C(n+2)=C(n)+3          for odd n.                    (2)
```

So the n+2 move preserves the current bit and reverses the next bit.
At a deeper dyadic scale, if w is the first L-bit word of n, then

```text
word_L(n+k*2^L) = w,
C^L(n+k*2^L) = C^L(n) + k*3^a                       (3)
```

whenever the arguments lie in the chosen integer domain. Equation (3)
is the precise residue information a structural halving operation must
respect. For H the same shift changes the output by `k*3^L` instead.
The valuation row and the branch count alone lose this distinction.

## 3. The adjacent-exchange law has opposite sign (PROVED)

Let u have length j; let v have length h and b ones. Directly applying
the recurrences in (1) gives

```text
R(u01v)-R(u10v) =  2^j 3^b,
M(u01v)-M(u10v) = -2^j 3^h.                         (4)
```

For the first identity, the shared prefix contribution is multiplied
by 3 in either two-bit ordering; the new contributions are `2^(j+1)`
and `2^j`. Each later odd bit multiplies their difference by 3. For
the second identity the shared prefix gets multiplied by 9; the new
contributions are `2^(j+1)` and `3*2^j`, and every later bit multiplies
their difference by 3.

This explains the change, rather than merely finding unequal values:
**Collatz's even branch only halves; Mahler's even branch triples and
then halves.** Moving an odd bit right increases the Collatz carry but
decreases the Mahler carry. One cannot identify their carry order by
keeping the same bit word.

An immediate Collatz corollary, by bubbling ones to the two extremes,
is

```text
3^a-2^a <= R(w) <= 2^(L-a)(3^a-2^a).               (5)
```

The extremes occur at `1^a 0^(L-a)` and `0^(L-a)1^a`, respectively.
At fixed `(L,a)`, all carries are distinct: equality of two carries
would give the same parity cylinder in (1), hence the same word.
Thus replacing all words by their fixed clock loses exactly
`binomial(L,a)` distinguishable possibilities, not a cosmetic ordering.

## 4. Root certificates become explicit Mahler hostiles (PROVED)

Both maps are conjugate on the 2-adic integers to the one-sided binary
shift by their parity vectors. The shared-word change of coordinates
is consequently exact as a 2-adic map. It does not preserve the
positive ordinary integers.

On odd-denominator rational inputs below, H denotes the extension
`H(x)=(3x+e(x))/2`, where `e(x)` is its residue modulo 2. This is not
the literal real ceiling function on rational inputs; the ceiling
description applies on ordinary integers.

Here is the complete root calculation. The C itinerary of 1 is
`(10)^infinity`. Under H a period `10` has multiplier 9, denominator 4,
and carry 3, so its fixed starting state is

```text
x = (9x+3)/4,          x=-3/5.
```

Root 4 has itinerary `00(10)^infinity`, and two H-even steps multiply
by `9/4`. Its shared-word state is therefore

```text
Psi(4) = (4/9)(-3/5) = -4/15.                       (6)
```

Suppose w is a finite C route from n to 4, with length L. The unique
H-state having the same full itinerary is then

```text
Psi(n) = (2^L*(-4/15)-M(w))/3^L < 0.                (7)
```

It is an odd-denominator negative rational, hence a valid 2-adic
integer but not an ordinary positive integer. Formula (7) proves the
domain loss for every certified member of the root basin, without
assuming all integers lie in that basin.

There is a second, independent failure. Mahler's suffix criterion is

```text
Y_j(e) = sum_(k>=0) e_(j+k)(2/3)^(k+1) < 1
         for every j.                              (8)
```

Every C root-basin word eventually has the suffix `(10)^infinity`.
Its full tail is `6/5`; already the five-bit prefix `10101` contributes
`266/243>1`. This exact hostile value is inherited from THM-4072,
section 6, centered/oriented trap; its application here is that every
root certificate supplies such a rejecting suffix. Reusing a Collatz
root word therefore fails both Mahler gates: positivity and real-tail
safety.

**Connection card.** Source: a finite C root certificate. Target: the
H parity-state space. Map: preserve the entire parity itinerary, with
formula (7). Preserved predicate: orbit shift and each parity bit.
Destroyed information: ordinary positive domain; Mahler's safe suffix
constraint is violated. Needed sidecar for a different transfer: a new
block coding plus its exact real interval and ordinary-realization
condition. Cheapest decisive test: root 4 itself, (6), together with
the `10101` suffix.

This does not contradict THM-4469. Its exact transfer uses restricted
equal-length/equal-weight block alphabets, adjacent carries, the
ratio `3^a/2^L`, and a specially aligned confinement interval. The
smallest recorded instance is `(L,a)=(10,7)` and ratio `2187/1024`,
with carries 4726,4727 and interval `[4726,4727]/1163. It is still
[OPEN HYP-9134](../hypotheses/HYP-9134-mahler-bridge-smallest-instance.md).
Mahler's original ratio 3/2 and interval `[0,1/2)` are not automatically
the same instance.

## 5. Catalan's role is the clock, not the decoder (CITED + PROVED)

Mihailescu's theorem settles the consecutive-perfect-powers equation:
with bases and exponents greater than 1, the unique positive solution
to `x^a-y^b=1` is `3^2-2^3=1`. Primary publication:
[Mihailescu, 2004, DOI 10.1515/crll.2004.048](https://doi.org/10.1515/crll.2004.048).
The publisher record was checked 2026-09-25; the statement and the
cycle-use boundary were already audited in the inherited Pillai note.

For positive exponents K,L this gives the familiar small unit gaps
between `2^K` and `3^L`, including the exponent-one boundary. It says
nothing about the general ordered carry. A C-cycle word satisfies

```text
(2^L-3^a)n = R(w),                                  (9)
```

whereas a route to 4 satisfies the different endpoint equality in
section 1. A small power gap is neither the route equality nor the
divisibility gate in (9). For the inherited signed odd seven-cycle,
the exact data are

```text
17 -> 25 -> 37 -> 55 -> 41 -> 61 -> 91 -> 17,
halving exponents (1,1,1,2,1,1,4),
2^11-3^7=-139,       carry B=2363=17*139.
```

Here the map on odd numbers is `(3n-1)/2^v`; the carry cancels a
non-unit gap. This is a hostile control against using Catalan as an
exhaustion of possible Collatz cycles. Near powers can locate clocks
with weak expansion/contraction. The decoder still needs carry order,
sign, and integrality.

For primary Mahler context, see
[Mahler's original 1968 paper, EMS reprint](https://ems.press/books/dms/252/4994)
and the authors' current
[Andrieu--Eliahou--Vivion v2](https://arxiv.org/abs/2510.11723v2), revised
2026-04-07. The latter still treats Z-number existence as open and
offers a conditional normality route. Neither is a claimed solution.
The precise carry and suffix formulas used above are proved in this
note or the cited repository canon, not imported from a search snippet.

## 6. Reproduction, controls, and stopping boundary

Run:

```text
python 04-computation/experiments/decoder_mahler_catalan_20260925.py
python -O 04-computation/experiments/decoder_mahler_catalan_20260925.py
```

The [frozen output](decoder_mahler_catalan_20260925.out) agrees in normal
and optimized execution. Universe: all 8190 nonempty binary words of
length at most 12, 24570 integer lifts, and all 20481 occurrences of
the `01 -> 10` exchange in that universe. Both carry formulas are
checked by a closed sum independent of their recurrences; cylinder
classes are independently checked by direct integer iteration.
All 4094 starts `3..4096` reach 4 within the stated 10000-step control
budget (maximum 148), and their negative rational H-images reproduce
the same parity words through the root route and 18 further steps.
This bounded root census is a control, not new convergence evidence.

Positive controls: exact route `5 -> 8 -> 4`, all dyadic lifts of the
sampled cylinders, and the known signed seven-cycle. Hostiles: clock
`(2,1)` loses integrality, the root H-image is `-4/15`, and `10101`
violates the Mahler tail bound. No random sampling or floating-point
arithmetic is used.

The board after these checks: the halving layer has an exact arithmetic
interface, but the tournament construction must still produce it.
Catalan constrains one clock coordinate. Mahler provides a productive
language of carries, real tails, and integer realization, while the
same-word transfer provably leaves the desired domain. The next
decisive experiment is an orientation rule whose contraction supplies
a certified `(L,a,R)` and whose residual labels satisfy (3). A mere
decrease in vertex count or a finite core cannot replace that predicate.
