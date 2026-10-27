# Session73 changing-modulus carrier packet

Read [REPORT.md](REPORT.md) for proofs and scope. The packet combines:

- exact least guarded-affine source precision, including arbitrary-modulus
  Collatz cylinders and unchanged integer clocks;
- a sharp same-source weighted Gram interval for paid pointwise forcing;
- exact CRT marginal rank and a two-debt decomposition separating source
  aliases from joint-incidence loss;
- a pinned refutation of incoming HYP-9419(4)'s polynomial-degree exclusion.

The uniform finite examples are separate illustrative source models, not a
replacement for the original native ROOT weights. No native ROOT positivity,
Collatz convergence, or all-word tournament linear exception bound is claimed.
The general theory retains supplied masses and actual guard support.

## Reproduce

From this directory (Python standard library only):

```sh
python3 -B check.py > /tmp/session73-moduli-normal.json
python3 -B -O check.py > /tmp/session73-moduli-optimized.json
cmp /tmp/session73-moduli-normal.json /tmp/session73-moduli-optimized.json
cmp /tmp/session73-moduli-normal.json receipt.json
python3 -B seal.py --verify
sha256sum -c SHA256SUMS
```

All mathematical validity gates raise explicit exceptions; none use `assert`.
The normal and optimized committed receipts are byte-identical. Exact arithmetic
uses integers and `fractions.Fraction`; primality controls use integer trial
division. There is no atlas filter, external package, random sample inference,
or floating-point rank. The 400 affine-word controls are reproducible finite
samples using the fixed seed731042, not an exhaustive word theorem.

## Exact finite controls

- 400 affine guarded words, including149 zero-slope words, with all22,678
  shifts in two full periods; cancelled-endpoint guard hostile.
- Empty word plus all84 valuation words of lengths1–3 and exponents1–4;
  q=1–18, giving1,530 word/modulus pairs,14,535 endpoint image tests,
  5,634 guarded split compositions and77,976 literal odd transitions.
- q,m=1–12:144 full residue frames and6,084 compatible/incompatible cell pairs;
  independent exact incidence elimination, explicit uniform residuals.
- 36 nonuniform q,m<=6 frames,103 named targets,1,545 exact extremal witnesses,
  and103 independently grounded Laplacian solves;101 one-cell thresholds.
- All unordered modulus banks of sizes1–3 drawn from1–8, plus four declared
  controls, each at N=L and N=2L:168 banks/336 fixed-N frames. Nested,
  repeated, noncoprime, and the joint repair bank{2,3,6} are included.
  Twenty-five nonuniform frames test the weighted orthogonal alias/incidence
  split and the single named joint-cell repair. Empty banks are not defined.
- A256-source synthetic nonuniform prior tests unchanged weighted pushforward
  through the native word(1,2) and joint endpoint moduli6,10, at clock5.
- Every odd prime p<=199 and exponent d=1–60:2,700 monomial criteria;
  3,540 odd-degree opposite-pair controls at moduli3–120;60 polynomial
  degree-alias controls; explicit x^4 mod7 and nonunit-square mod15 hostiles.

## Provenance and validation status

`SOURCE_PINS.json` gives20 exact commit/blob/SHA256 source pins; snapshots
retain native bytes. Main inheritance is pinned to92e7f12f3b95. The bounded
incoming correction is pinned to64a7985e705. `MANIFEST.json` and `SHA256SUMS`
seal this packet without recursively hashing themselves.

The producing lane supplied self-contained proofs and separate exact
formulations; the independent session auditor reviewed the draft's P1–P3 and
W1–W3 and requested the odd-map naming clarification, which was applied.
The parent session owns the definitive independent audit/replay record and
publication. This packet alone does not claim a separately sealed final audit.
