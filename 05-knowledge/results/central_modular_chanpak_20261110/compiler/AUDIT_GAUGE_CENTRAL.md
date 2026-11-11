# Independent audit of the affine gauge and source cylinders

**VERIFIED — algebraic audit, 2026-11-10.** Audited
[gauge G1–G3](../gauge/PROOF.md) and [central C1–C3](../central/PROOF.md)
against their stated domains. No new computational universe or classification
claim is introduced. No mathematical correction is required.

## Affine gauge and the tournament bridge

For the ordinary clock, the even branch gives
`a*x/2+b = (a*x+b+1)/2`, hence `b=1`. The odd branch then gives
`3*a*x+a+1 = 3*a*x+3+delta`, hence `a=delta+2`.
Both equalities hold conversely. Positivity is consistent for the three
specified values `delta=0,2,4`.

The forbidden labels occur inside the claimed positive-integer domains:

| Gauge | Positive input with forbidden label | Allowed even-source boundary |
|---|---|---|
| `2*x+1` | `x=3 -> 7`, `x=10 -> 21` | `x=6 -> 3`: `13 -> 7`; `x=20 -> 10`: `41 -> 21` |
| `4*x+1` | `x=5 -> 21` | `x=10 -> 5`: `41 -> 21` |
| `6*x+1` | `x=1 -> 7` | `x=2 -> 1`: `13 -> 7` |

The all-order exclusion used here is inherited from **THM-1370**,
[h-spectrum-omits-7-21-all-n](../../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).
Its stated attained interval also certifies existence at 13 and 41. Thus
universal assignments `H(T_x)=a*x+1` fail, and merely deleting the forbidden
source labels does not produce an invariant domain under both branch laws.
The proof does not exclude a smaller guarded invariant domain, non-affine
counts, a larger observer, or an arithmetic action on constructor addresses.

In particular, `4*x+1` is exactly the plus moving-modulus conjugacy, but the
fact that a label is a modulus supplies no tournament having that many paths.
The repair in G3 correctly retains this distinction: the carry lift computes
the new count `Phi(A,B')` of an actual reconstructed graph. It neither claims
`Phi(A,B)=4*B+1` nor closure of the old joint-port observer.

## Central source-cylinder induction and clock

For either sign, the ordinary odd branch maps every positive odd integer to
a positive even integer. Hence `OO` is impossible and an `E` after `O` is
forced. If the current word ends in `O`, appending this forced `E` changes
`(b,delta)` from `(b,1)` to `(b+1,0)`: the exponent `b+delta` stays fixed.

Otherwise the existing source class has modulus `2^b`. Its two lifts modulo
`2^(b+1)` differ by `2^b`; the endpoint formula makes their current outputs
differ by `3^a`, which is odd. Exactly one lift has either required next
parity. Appending `O` or `E` increases `b+delta` by one. This proves both
necessity and sufficiency of the claimed single source class, including the
empty word. Positivity presents no boundary exception for either sign, since
`3*x-1 >= 2` on positive odd inputs. Every such class has infinitely many
positive representatives.

The carry recurrences are the direct compositions: an ordinary `O` sends
`B` to `3*B+2^b`, whereas `E` changes only the denominator exponent. The
endpoint-parity bit is indispensable: the one-letter word `O` already needs
modulus 2 despite having no halving. Order is also indispensable, as `OEE`
and `EOE` have different carries and source classes.

The text correctly distinguishes three operations:

- ordinary odd step `x -> 3*x+epsilon`;
- one-halved odd step `2*n+1 -> 3*n+(3+epsilon)/2`;
- full acceleration through every subsequent factor of 2, which is not used.

The no-`OO` language and exponent `b+delta` belong to the first clock.
The multiplier-three carry channels in compiler C4 use the second. Neither
clock statement is substituted for the other. C2's even-source division by
6 is a decoder of the displayed central pair, and is also correctly
distinguished from reversing a directed edge or taking an odd shortcut.

## Audit conclusion

The arithmetic maps and cylinders are exact within their specified universes.
The remaining tournament obligation is still a native chart switch carrying
named ports and readiness, not an unproved scalar interpretation of the
Collatz labels. The files state that boundary explicitly.
