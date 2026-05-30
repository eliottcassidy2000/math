# HYP-1804: Pinned triangle variance is the additive-energy axis

**Status:** OPEN proof target; exhaustively verified for circulant orientation cubes through `p=23`.

## Statement

For a circulant tournament on `Z_p` with connection set `S`, let

- `E(S)=|{(a,b,c,d) in S^4 : a+b=c+d mod p}|` be additive energy.
- `J_3(0,v)` be the number of directed 3-cycles passing through both `0` and `v`.
- `Var(J_3)` be the variance of `J_3(0,v)` over `v != 0`.

Then, conjecturally/proof-target:

```text
Var(J_3) = E(S)/(p-1) - (p^2 - 2p + 5)/16
```

Equivalently, pinned triangle localization, additive energy, and spectral IPR
are all the same one-dimensional axis on the circulant orientation cube.

## Evidence

`04-computation/fejer_kernel_wild_session.py` fits the affine law over every
valid circulant tournament for

```text
p = 7, 11, 13, 17, 19, 23
```

with slope `1/(p-1)`, intercept `-(p^2 - 2p + 5)/16`, and roundoff-scale
maximum error.

The same run verifies the Fejer/autocorrelation identity for interval sets up
to `p=31` and records that `corr(IPR,J3var)=1.0000` on the `p=7,11,13`
Hamiltonian-path cubes.

## Why It Matters

This sharpens the S64 proof gap. The first bridge

```text
Fejer / additive-energy concentration -> triangle spatial localization
```

is probably not the hard part. The hard part is the next bridge:

```text
triangle localization -> higher simple-cycle localization -> Omega packing
```

or, more concretely, proving that localization survives the repeated-vertex
corrections in pinned `J_5,J_7,...` profiles strongly enough to control the
OCF hard-core partition function.

## Next Steps

- Prove the affine law by expressing `J_3(0,v)` as an oriented wedge transform
  of the connection-set indicator.
- Build the analogous decomposition for `J_5`: pinned walk profile minus
  repeated-vertex correction diagrams.
