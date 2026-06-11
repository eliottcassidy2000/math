# HYP-2440 - d16+ is excluded from the order-5 fixed projection

**Status:** CLAIMED / exact finite proof. Codex 2026-06-11.

In the `d16+` Type II `[16,8]` code, after pairing coordinates as

```text
(0,1), (2,3), ..., (14,15),
```

the weight-4 codewords are exactly the unions of two coordinate pairs. Hence
there are `C(8,2)=28` tetrads.

This tetrad system covers every coordinate pair:

- two coordinates in the same `d16+` pair occur together in `7` tetrads;
- two coordinates in different `d16+` pairs occur together in exactly `1`
  tetrad.

Therefore any choice of two fixed marks in a `d16+` projected fixed code lies
inside at least one tetrad. That tetrad lifts under a `5-(14,2)` automorphism
to an original fixed codeword of weight `12`, contradicting extremality of a
`[72,36,16]` code.

## Consequence

The skew-doubling / `d+` branch from THM-480 cannot be the fixed projection of
an order-5 automorphism of an extremal `[72,36,16]` code. The order-5 branch, if
it exists, must use the decomposable `e8+e8` fixed projection, even though
`d16+` is the more natural object in the tournament-doubling tower.

This is a small but real fork: doubling-forgets naturally produces `d+`, while
the order-5 extremal fixed gate demands split `e8`.

## Links

- Sharpens HYP-2409/THM-480 and HYP-2439.
- Computation: `04-computation/order5_fixed_projection_72_codex.py`.
