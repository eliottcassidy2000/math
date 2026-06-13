# HYP-2214: Pi/e Algebraic Shadows Are Lonely

**Status:** OPEN / proved conditional fallout.

## Claim

Let

```text
S = e + pi,    P = e*pi,    D = e - pi,    L = log(pi),
p_k = e^k + pi^k.
```

The Vieta carrier `T^2 - S*T + P` does more than prove that `S` and `P`
cannot both be algebraic.

1. If `L=log(pi)` is algebraic, then `S`, `P`, and `D` are all
   transcendental.
2. If any one of `S`, `P`, or `D` is algebraic, then `log(pi)` is
   transcendental.
3. If `S` is algebraic, then every `p_k` for `k>=2` is transcendental.
4. If `P` is algebraic, then every `p_k` for `k>=1` is transcendental.

Thus an algebraic `e+pi` or `e*pi` would be a lonely algebraic shadow: it
would force the transverse symmetric shadows out of `Qbar`.

## Proof Notes

If `L` were algebraic, then `pi=e^L`.  Since `pi` is neither `1` nor `e`, the
exponents `0`, `1`, and `L` are distinct.  Lindemann-Weierstrass makes `1`,
`e`, and `e^L=pi` linearly independent over `Qbar`; hence `e+pi` and `e-pi`
cannot be algebraic.
Also `e*pi=e^(1+L)` is transcendental by Hermite-Lindemann.  This proves the
first claim and its contrapositive.

For power sums, Newton's recurrence gives

```text
p_0 = 2,  p_1 = S,  p_k = S*p_{k-1} - P*p_{k-2}.
```

So `p_k=f_k(S,P)`.  If `S` is algebraic and `p_k` is algebraic for some
`k>=2`, then `P` is algebraic over `Qbar` because `f_k(S,X)` is nonconstant in
`X`; then `T^2-S*T+P` makes `e` and `pi` algebraic.  The same argument with
`S` and `P` exchanged proves the `P`-algebraic fallout.

## Repo Meaning

This extends HYP-2211/HYP-2212 from "two shadows reconstruct the hidden pair"
to "one descended shadow makes transverse shadows forbidden."  The useful
analogy is the incoming H=21 window closure: the scalar alone is too thin; the
result becomes sharp only after retaining side conditions (`strong`, `c3<=10`
there; log-commensurability, discriminant sheet, and transverse symmetric
polynomial fallout here).

## Next Tests

- Build a generic transverse-shadow helper for any two-root carrier.
- Apply the same lonely-shadow logic to LRC reset-period/relation-lattice
  ledgers.
- Search tournament fibers with fixed scalar `H` but varying transverse packets
  (`c3`, SCC, beta, OCF packet) as a finite analogue.

**See also:** `04-computation/pi_e_lonely_shadow_s638.py`,
`05-knowledge/results/pi_e_lonely_shadow_s638.out`,
`07-reflections/pi-e-lonely-shadow-fallout-s638.md`, HYP-2211, HYP-2212,
HYP-2213, HYP-2200.
