# Session143 bounded follow-up: Frobenius-cube magic over F5

**PROVED from cited Frobenius facts; FINITE-EXACT census independently audited.**
This is one precise candidate meaning of the unresolved lambda-cubed cue,
**not a user-confirmed interpretation**. It is not a magic square of squares.

For an elliptic curve over F5, put `t = 6 - #E(F5)`. Then

    Tr(Frob_5^3) = t^3 - 15t,
    #E(F125) = 126 - t^3 + 15t.

This explains the same exceptional cubic in Session143's complete real cubic
magic classification. All nine traces `t=-4,...,4` occur among the 20
nonsingular short Weierstrass models over F5. Direct extension-field counts give

    108  148  122
    140  126  112
    130  104  144

Every row, column and main diagonal sums to **378**. Only 144 is a square.
There are eight arrangements of these nine distinct numbers, related by the
symmetries of a square. This is a **new magic rearrangement** of the extension
counts, not pointwise preservation of a base-field magic layout. No fixed
placement of all nine traces can make both layouts magic: a trace-magic grid
has center zero, while an outer line has nonzero `a,b,c` with `a+b+c=0`.
Thus its power-cube trace sum is `a^3+b^3+c^3=3abc != 0`; the linear correction
`-15(a+b+c)` vanishes and cannot restore that magic line. The corresponding
F125-count line sums to `378-3abc != 378`, contrary to extension-count magic.

## What the connection does and does not preserve

- It maps the scalar labels `n=#E(F5)` in `2,...,10` to point counts over F125.
  The user's cycle is relabeled, not realized by a curve-level successor map.
  This does not derive its `+4 mod 9` ordering or canonical curve representatives.
  Distinct base-field point counts also rule out F5-defined isogenies between
  adjacent cycle states; no Hecke/isogeny meaning is assigned to the arrows.
  [Cited isogeny invariant](https://public.websites.umich.edu/asnowden/teaching/2013/679/L03.html).
- The Frobenius **power trace** is not the symmetric-cube trace or the Hecke
  prime-cubed coefficient: those latter expressions are `t^3-10t` at `q=5`.
- For a formal common-determinant family with all nine specified traces,
  the inherited cubic classification forces `q=5` for power-cube magic;
  symmetric-cube magic would require `q=15/2`, not a finite-field size.
- Both elliptic curves and five-vertex tournaments have 12 classes and nine
  scalar invariants here, but their fiber partitions disagree. Quadratic
  twisting has six class orbits; tournament converse has ten. These exclude
  the proposed invariant- and involution-preserving identifications.
- No squarehood transport, octonion equivalence, general MSOS solution,
  all-extension search, new conjecture, or global originality claim is made.

## Recovered lambda provenance, still not an identification

A prior repository annotation uses
`lambda=(1+3 phi)^(-1/3)`, satisfying `5 lambda^6+5 lambda^3-1=0`,
where `phi=(1+sqrt(5))/2`. Its underlying quadratic element has norm `-5`,
not the `+5` determinant of the Frobenius operator here. Its inverse
`r=1/lambda` satisfies `r^6-5r^3-5=0`. The universal two-dimensional identity
`tr(M^3)=tr(M)^3-3 det(M) tr(M)` explains why the determinant sign matters;
it does not identify the old lambda with a Frobenius eigenvalue. The exact older source
is pinned in [SOURCE_PINS.json](SOURCE_PINS.json). Neither expression is assumed
to be the user's intended lambda; the two constructions are kept distinct.

## Evidence and replay

- [Proof, connection contract and cited inputs](producer/PROOF.md)
- [Independent audit](independent/PROOF.md)
- [Verification scope](VERIFICATION.md)
- [Exact source pins](SOURCE_PINS.json)
- [Replay output](REPLAY.json) and [payload manifest](MANIFEST.json)

Run from any working directory with Python 3.9 or newer, without third-party
packages, Git, network, or the original repository:

```sh
python3 /path/to/this/packet/replay.py
python3 -O /path/to/this/packet/replay.py
python3 /path/to/this/packet/verify_manifest.py
```

The two replay outputs must equal `REPLAY.json` byte for byte. Each replay
runs both independent checkers and their positive/mutation controls in normal
and optimized modes. Checks use explicit exceptions, not removable asserts.

The [original Session143 packet](../session143_magic_square_cyclic_bridges_20261103/README.md)
is preserved byte for byte at its 183-file published subtree. This addendum
extends the result, not its original frozen evidence or unresolved referent.
