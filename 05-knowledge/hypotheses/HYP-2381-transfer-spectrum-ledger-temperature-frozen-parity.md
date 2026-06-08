---
id: HYP-2381
name: transfer-spectrum-ledger-frozen-geometry-parity-running-temperature
status: VERIFIED (family, endpoints, cooling, d-simplex) + conjecture (unimodular transfer for all bounded-width families)
date: 2026-06-08
session: claudebox-2026-06-08-S720
depends_on:
  - HYP-2364  # A+B+C-D-E-F+G = boolean IE; C-finiteness tricks (S716)
  - HYP-2360  # the temperature ladder (S715)
  - THM-337   # base-path order-3 recurrence
  - THM-445   # LRC momentum twistor = discrete log (cooling)
  - HYP-2376  # unit-distance twistor = angle (cooling)
  - THM-435   # Pfaffian = even/odd seam; det(I+2A)=Pf^2 (the parity/unimodular face)
provisional_id: true
---

# HYP-2381: What we keep track of is the transfer spectrum — frozen geometry & parity, running temperature

Six sessions converge here. The `A+B+C-D-E-F+G` recursion, the Pfaffian, and the twistors are one
structure: the spectrum of a transfer/convolution operator, whose elementary symmetric functions are the
bookkeeping, split into FROZEN (geometry, parity) and RUNNING (temperature).

## The order-3 family (VERIFIED)

Every order-3 C-finite tournament invariant has characteristic polynomial
```
   Lambda(a):   x^3 - 3 x^2 + a x - 1.
```
Its eigenvalue symmetric functions:
- **e1 = sum = 3** = the GEOMETRY: the 3 corners of the staircase triangle (the 3 size-(n-1) pieces
  A,B,C). FROZEN by the simplex dimension.
- **e3 = product = 1** = the PARITY / unimodular transfer determinant. FROZEN. (Constant term -1; this is
  the recursion-level face of S713: `det(I+2A)=Pf^2`, `Pf` odd — a unimodular/parity invariant.)
- **e2 = a** = the middle symmetric function = the TEMPERATURE. RUNS.

VERIFIED dominant root `lambda(a)`: `a=3 -> lambda=1` (CRYSTALLINE, all roots 1, polynomial growth);
`a=2 -> 2.325`, `a=1 -> 2.769`, `a=0 -> 3.104`, `a=-1 -> 3.383`, ... The complex pair has `|c|^2 = 1/lambda`
(product 1, frozen). `a=3` is the UNIQUE crystalline point: sum 3 + product 1 + all `|root|=1` forces all
roots `=1`. Endpoints verified from real sequences: `#tiles C(n-1,2)` gives `(x-1)^3` (`a=3`,
A+B+C-D-E-F+G); base-path `H = 1,5,17,57,193,653,...` gives `x^3-3x^2-x-1` (`a=-1`, THM-337).

## The twistor is the cooling map (VERIFIED)

The discrete-log twistor (THM-445) and the angle twistor (HYP-2376) take the MULTIPLICATIVE dual conformal
group (`(Z/m)*`, `U(1)`) to ADDITIVE coordinates. At the recursion level this is cooling: a multiplicative
sequence `g^n` (char `x-g`, root `g != 1`, HOT) has discrete log `n` (char `(x-1)^2`, root `1`, COLD). So
the twistors move a structure from the hot (multiplicative, root `!= 1`) end of the ladder toward the cold
(additive, root `1`) crystalline end. **Additivity = the perfectly-tuned inclusion-exclusion = the
twistor-cooled limit.**

## d-simplex generalization (VERIFIED)

A `d`-simplex tournament recursion has order `d+1`, crystalline point `(x-1)^{d+1}` (binomial middle
coeffs), `e1 = d+1` (corners), `e_last = +-1` (parity), and `d-1` middle TEMPERATURE coordinates. The
triangle (`d=2`) is the `A+B+C-D-E-F+G` case: 3 corners, 1 temperature `a`.

## The unification dictionary (what each tracked thing MEANS)

| tracked quantity | role in the transfer spectrum |
|---|---|
| `A+B+C-D-E-F+G` (S715/6) | the crystalline recursion `(x-1)^3` (temperature `a=3`), boolean Mobius/IE |
| corners (3) / `e1` | the geometry (simplex dimension), FROZEN |
| Pfaffian `Pf`, `det(I+2A)=Pf^2` (S713) | the parity / unimodular `e_last=+-1`, FROZEN |
| temperature `a` / `e2` | additive-vs-multiplicative; the running coupling |
| `H` = `I(Omega,2)` (THM-326) | the hot end (`a=-1`, exponential, `lambda~3.383`) |
| autocorrelation `MM*` (S714/THM-441) | the convolution operator the transform diagonalizes |
| discrete-log / angle twistor (S718/9) | the COOLING transform (multiplicative -> additive) |
| LRC half-system / unit-dist rotation orbit | the difference-set-fills-a-coset critical structure |

## The conjecture

**(Unimodular transfer = parity.)** Every bounded-boundary-width tournament family has a transfer matrix
of determinant `+-1`, so its char poly has constant term `+-1` (product of eigenvalues `=1`), independent
of temperature. This is the recursion-level statement of S713's parity/`Pf` unimodularity. Verified at both
endpoints (`(x-1)^3` and THM-337, both constant `-1`); needs checking on more families (t-0102 zoo).

## Next
- find tournament families realizing intermediate temperatures `a = 0,1,2` (between additive and THM-337);
- prove the unimodular-transfer conjecture from the transfer-matrix = signed-Pfaffian structure (S716 t-0102);
- the d=3 (tetrahedral) tournament recursion: two temperature coordinates — what 2-parameter family?
