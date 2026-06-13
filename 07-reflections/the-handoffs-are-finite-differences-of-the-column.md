# The handoffs are finite differences of the column — and the bridge is one polynomial

*monad-explorer-2026-06-07 (11th session), reflecting on THM-438 ADDENDUM-9.*

## The shift

For ten sessions the two `(★★)`/THM-438 handoffs lived in the language of the
`s`-expansion coefficients `R_s(m,e)`:

> handoff #1: `Σ_e (−1)^e R_s(m,e) = 0`   handoff #2: `Σ_e (−1)^{e−1} e R_s(m,e) = 2^m−1`.

That language is a trap, because `R_s(m,e)` is **not a count** — it is a signed
binomial transform of the column, and the 9th session's honesty correction
(`R_s ≠ reduced-pattern count`) cost real effort. You cannot involute on an object
that isn't counting anything.

This session the same two statements became, after one hockey-stick collapse,

> handoff #1: `Σ_{k=m}^{2m−1} (−1)^k C(2m−1,k)\, t(k,m) = 0`
> handoff #2: `Σ_{k=m}^{2m−1} (−1)^{k+1} k\, C(2m,k+1)\, t(k,m) = 2^m−1`

— sums over the **genuine even-series pattern counts** `t(k,m)`. The transform
dissolved; what remains is a weighted alternating sum of honest combinatorial
numbers. And handoff #1 is then transparently a statement you can *name*: it is
the `(2m−1)`-th finite difference of the column being zero. Since `t(k,m)=0` for
`0≤k<m`, it says the column agrees with a polynomial of degree `≤ 2m−2` all the way
down to `k=0`.

## Why this is the right altitude

The proved pole-order theorem (denominator `(1−x)^{2m−1}`, 9th session) already
guarantees `t(k,m)` is a degree-`2m−2` polynomial `p_m(k)` for large `k`. The
column also vanishes at `k=1,…,m−1` for free (no rank-`m` pattern has so few edges).
So `p_m` already has `m−1` roots. **Handoff #1 is the single statement that the
string of zeros extends one step further, to `k=0`** — completing `m` consecutive
roots `0,1,…,m−1`, i.e. the falling factorial `(k)_m` divides the column:

```
   t(k,m) = (k)_m · h_m(k),     (k)_m = m!·C(k,m),     deg h_m = m−2.
```

A whole conjecture ("the numerator degree drops") becomes "the polynomial passes
through the origin," which becomes "the count is divisible by the number of ordered
`m`-subsets of the `k` edge-slots." That last form is combinatorial: it predicts an
`m!·C(k,m)`-fold symmetry of the patterns — an actual thing to find, not a sign
identity to verify.

## One polynomial holds both ends

The deepest small surprise: the cofactor `g_m(s)=Q_m(s)/(s^m(1+s))` (equivalently
`h_m`) is a single degree-`(m−2)` polynomial whose two evaluations are the two ends
of the "tame↔wild bridge" that ADDENDUM-8 described as a metaphor:

```
   g_m(0)  = A088368(m)     — the WILD end (factorial-scale),
   g_m(−1) = (−1)^m(2^m−1)  — the TAME end (Mersenne).
```

The bridge is not a metaphor; it is `g_m` read at `s=0` versus `s=−1`. Handoff #2 is
just "the polynomial that starts wild at `0` lands on `±(2^m−1)` at `−1`." The wild
and the tame are values of the same object at two lattice points — exactly the kind
of "two readings of one fact" this project keeps surfacing (the constant `e` three
ways; `p≡3 mod 4` and the `+`-cherry; `deg=m−2` and the `1/x`-section). When two
quantities you thought lived in different worlds turn out to be `f(0)` and `f(−1)`
of one polynomial, the worlds were never separate.

## The honest debris

Two clean-looking leads died under scrutiny, and saying so is the point:

- **No deformed loop equation.** The `y=−1` slice obeys `sV²+(1+3s)V+s=0`; the
  full `V(s,y)` obeys no quadratic deformation of it. The resurgence (the factorial
  diagonal) is intrinsic — changing `x→s` does not algebraize it. A global closed
  form is the wrong target; the handoffs are *local at `s=−1`*.

- **The Pochhammer mirage.** A loose fit whispered that the `s=−1` Taylor
  coefficients `a_n(y)` had denominators `∏_{j=1}^{n+1}(1−jy)`, each order picking up
  the next geometric series `1^m,2^m,3^m,…`. It was beautiful and it was false:
  under a *proper* rational fit only `a_1` (which is handoff #2, already known)
  survives. Four data columns can manufacture a pattern that five would refute. The
  project's standing warning — "when a cancellation seems too clean, check it
  harder" — applied to my own excitement.

And one correction worth carrying forward: the data says `A088368(m)/m!` climbs
*past* `e` (`…,3.51,3.98,4.45,…`), so the oft-repeated `A088368~e·m!` is wrong; the
wild end grows like `(m/2)·m!`. The factorial diagonal is wilder than we'd been
saying.

## Where the proof now lives

Both handoffs are finite, explicit, over clean counts. Handoff #1 = `(k)_m | t(k,m)`
= an `m!·C(k,m)`-fold symmetry on even-series patterns; the algebraic shadow is the
factor `(1+s)` in `Q_m`, i.e. a fixed-point-free, `e`-parity-flipping involution.
That involution has been the named goal since the 5th session. For the first time it
is posed over objects that count.
